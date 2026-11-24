"""F1's run: the tree's harmonizer on the published datasets, its answers and its outputs written for the scorer.

Usage: run_f.py --organism <o> --clone <checkout of the tree> --graph <graph file> --out <families>/f/runs/<tree>/<o>
                [--set f1_v1]

Runs in an environment with idtrack's `harmonize` extra (anndata 0.13 or later; the package's own environment has none,
docs/families/f.md "Running it"), with the tree's checkout first on the path (as the common adapter puts it). Nothing is
judged here: what the harmonizer did is written down, and `score_f.py` judges it against the truth.

For each configuration (`gene`: `final_database="ensembl_gene"`, the Ensembl genes the truth names; `symbol`: the
organism's default final database), one `HarmonizeFeatures` over the set's datasets with `from_release_by_dataset` the
datasets' source releases, `to_release` 115, the loaded graph given as `api`:

- `<config>/answers.jsonl.gz`: one row per feature and dataset in the runner's answer form (`got`: the gene of a
  converted feature, the genes of an ambiguous one, nothing otherwise; `status`: the harmonizer's `feature_table`
  status; `seconds`: the dataset's conversion time divided by its features), with the conversion's own status,
  match, reason and `own` beside;
- `<config>/outputs/<merge>/<dataset>.tsv.gz` and `.meta.json` (`merge` "sum" and "flag", `keep_unconverted` True):
  every output feature (its name, its input features, status, gene, final-database name, the digests of its values in
  `X` and in the layer `unspliced`, its input `var` columns) and the object's `obs` (as text, with dtypes), `uns`,
  the digests of `obsm`, `obsp` and `varm`, the layers' names;
- `<config>/unify/<features>_<merge>.*` (`features` "union" and "intersect"): the merged object's features, per
  dataset the block of each feature (its input features, status, digests), its `obs` (`obs.tsv.gz`), `uns`, the
  digests of every dataset's block of `obsm`, `obsp` and the layers;
- `<config>/shuffle/union_<merge>.*`: the same union of the datasets in another order, each with its features
  permuted (seeded; the permuted files in `shuffle_datasets/`, one for both configurations), harmonized again from a
  copy of the result cache
  (whose stored feature names no longer match: every dataset converted again);
- `<config>/cache.json`: the result cache: its headers, the constructor again on it (seconds, the conversions equal),
  with a stored content digest changed (the file written again with the graph's digest, the conversions equal), and
  without `api` (whether the graph was loaded, seconds, equal).

`load.json` (the graph's load record), `input_check.json` (each dataset as anndata reads it: its features in order and
the digests of X and of the layer, for the scorer to compare with what the generator wrote), `run.json` (commit,
Python, anndata, seconds, peak memory). The result caches are removed when the run ends (their headers are kept in
`cache.json`); the permuted datasets are kept (inputs of the shuffle relation).
"""
import argparse
import gzip
import json
import os
import random
import resource
import shutil
import sys
import time

TARGET = 115
CONFIGS = {"gene": "ensembl_gene", "symbol": None}
MERGES = ("sum", "flag")
FEATURES = ("union", "intersect")
LAYER = "unspliced"
INPUT_VAR = ("feature_type", "highly_variable", "mean_input")  # the generator's var columns
BLOCK_HEADER = ["name", "dataset", "query_ids", "status", "x_digest", "u_digest"] + [f"var:{c}" for c in INPUT_VAR]


def write_json(path, d):
    with open(path + ".tmp", "w") as fo:
        json.dump(d, fo, indent=1, default=str)
    os.replace(path + ".tmp", path)


def write_tsv_gz(path, header, rows):
    with gzip.open(path + ".tmp", "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
    os.replace(path + ".tmp", path)


def text(v):
    from idtrack_reproducibility.families.f import h5ad
    return h5ad.as_text(v)


def obs_table(obs):
    """(header, rows, dtypes) of an obs frame as text."""
    cols = list(obs.columns)
    rows = [[str(i)] + [text(v) for v in vals] for i, vals in zip(obs.index, zip(*[obs[c].tolist() for c in cols]))] \
        if cols else [[str(i)] for i in obs.index]
    return ["obs_name"] + cols, rows, {c: str(obs[c].dtype) for c in cols}


def meta_of(adata, extra=None):
    from idtrack_reproducibility.families.f import h5ad
    _h, rows, dtypes = obs_table(adata.obs)
    return dict(n_obs=int(adata.n_obs), n_vars=int(adata.n_vars), x_dtype=str(adata.X.dtype), x_kind=type(adata.X).__name__,
                obs_columns=list(adata.obs.columns), obs_dtypes=dtypes, obs_names=[str(x) for x in adata.obs_names],
                obs_values={c: [r[1 + i] for r in rows] for i, c in enumerate(adata.obs.columns)},
                uns=h5ad.plain(dict(adata.uns)),
                obsm={k: h5ad.digest(v) for k, v in adata.obsm.items()},
                obsp={k: h5ad.digest(v) for k, v in adata.obsp.items()},
                varm={k: h5ad.digest(v) for k, v in adata.varm.items()},
                layers=sorted(str(k) for k in adata.layers.keys() if k is not None),
                var_columns=list(adata.var.columns), **(extra or {}))


def describe_dataset_output(adata, final_db):
    """Rows of one harmonized dataset's output features (module docstring)."""
    from idtrack_reproducibility.families.f import h5ad
    xd = h5ad.column_digests(adata.X)
    ud = h5ad.column_digests(adata.layers[LAYER]) if LAYER in adata.layers else [""] * adata.n_vars
    v = adata.var
    inputs = [c for c in v.columns if c not in (final_db, f"{final_db}_alternatives", "Query ID", "matching_ensembl_gene",
                                                  "idtrack_status", "idtrack_reason")]
    header = ["name", "query_ids", "status", "gene", "final_id", "reason", "x_digest", "u_digest"] + [f"var:{c}" for c in inputs]
    rows = []
    for j, name in enumerate(adata.var_names):
        rows.append([name, text(v["Query ID"].iloc[j]), text(v["idtrack_status"].iloc[j]),
                     text(v["matching_ensembl_gene"].iloc[j]), text(v[final_db].iloc[j]), text(v["idtrack_reason"].iloc[j]),
                     xd[j], ud[j]] + [text(v[c].iloc[j]) for c in inputs])
    return header, rows


def describe_merged(m, handles, n_cells):
    """(feature rows, block rows, meta) of a merged object: per feature its merged var columns, per dataset and
    feature present its block (input features, status, digests of X and the layer); blocks of absent features must
    be zeros (counted)."""
    import numpy as np

    from idtrack_reproducibility.families.f import h5ad
    v = m.var
    handle_col = m.obs["handle_anndata"].astype(str).to_numpy()
    rows_of = {h: np.where(handle_col == h)[0] for h in handles}
    feat_cols = [c for c in v.columns if not any(c.endswith(f"_{h}") for h in handles)]
    feat_rows = [[name] + [text(v[c].iloc[j]) for c in feat_cols] for j, name in enumerate(m.var_names)]
    blocks, zero_bad, obsm, obsp, layer_blocks = [], 0, {}, {}, {}
    X = m.X.tocsr() if hasattr(m.X, "tocsr") else m.X
    U = m.layers[LAYER] if LAYER in m.layers else None
    import scipy.sparse as sp
    off_block = 0
    for h in handles:
        r = rows_of[h]
        xd = h5ad.column_digests(X[r])
        ud = h5ad.column_digests(U[r]) if U is not None else [""] * m.n_vars
        z = h5ad.zero_digest(len(r))
        q = v.get(f"Query ID_{h}")
        st = v.get(f"idtrack_status_{h}")
        ins = [v[f"{c}_{h}"] if f"{c}_{h}" in v.columns else v.get(c) for c in INPUT_VAR]
        for j, name in enumerate(m.var_names):
            qq = text(q.iloc[j]) if q is not None else ""
            if qq:
                blocks.append([name, h, qq, text(st.iloc[j]) if st is not None else "", xd[j], ud[j]] +
                              [text(c.iloc[j]) if c is not None else "" for c in ins])
            elif xd[j] != z or (U is not None and ud[j] != z):
                zero_bad += 1
        obsm[h] = {k: h5ad.digest(np.asarray(a)[r]) for k, a in m.obsm.items()}
        obsp[h] = {k: h5ad.digest(sp.csr_matrix(a)[r][:, r]) for k, a in m.obsp.items()}
    for k, a in m.obsp.items():
        a = sp.csr_matrix(a)
        inside = sum(sp.csr_matrix(a[rows_of[h]][:, rows_of[h]]).nnz for h in handles)
        off_block += a.nnz - inside
    meta = meta_of(m, dict(feature_columns=feat_cols, nonzero_blocks_of_absent_features=zero_bad, obsm_blocks=obsm,
                           obsp_blocks=obsp, obsp_off_block_nnz=off_block, handles=handles))
    meta.pop("obsm")
    meta.pop("obsp")
    obs_h, obs_rows, _ = obs_table(m.obs)
    meta.pop("obs_values")
    return feat_cols, feat_rows, blocks, meta, (obs_h, obs_rows)


def peak_gb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", default="f1_v1")
    ap.add_argument("--configs", default="gene,symbol")
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
    t_start = time.time()
    clone = os.path.abspath(a.clone)
    sys.path.insert(0, clone)
    import anndata as ad
    import numpy as np

    import idtrack
    import idtrack.harmonize.features as HF
    from idtrack.harmonize import HarmonizeFeatures
    from idtrack_reproducibility.families import adapter as F0
    from idtrack_reproducibility.families import common as C
    from idtrack_reproducibility.families import lines as L
    from idtrack_reproducibility.families.f import h5ad
    log = lambda *x: print(time.strftime("%H:%M:%S"), *x, file=sys.stderr, flush=True)  # noqa: E731
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    man = L.manifest("f", a.organism, a.set)
    fdir = C.family_dir("f")
    datasets = {n: os.path.join(fdir, d["file"]) for n, d in sorted(man["datasets"].items())}
    releases = {n: d["release"] for n, d in man["datasets"].items()}
    handles = list(datasets)
    # the graph, linked into the run folder (nothing written beside the original)
    link = F0.link_graph(os.path.abspath(a.graph), out)
    repo = os.path.dirname(link)
    before = set(os.listdir(repo))
    t0 = time.time()
    api = idtrack.API.load(link)
    rec = F0.graph_manifest(api.graph)
    load = dict(generation=C.generation(clone), idtrack_file=idtrack.__file__, python=sys.version.split()[0],
                pythonhashseed=os.environ.get("PYTHONHASHSEED"), load_seconds=round(time.time() - t0, 1),
                graph=os.path.abspath(a.graph), graph_real=os.path.realpath(a.graph),
                files_written_in_repo=sorted(set(os.listdir(repo)) - before), content_digest=rec.get("content_digest"),
                package_commit=rec.get("package_commit"), slurm_job=os.environ.get("SLURM_JOB_ID"))
    write_json(os.path.join(out, "load.json"), load)
    log("graph loaded", load["load_seconds"], "s", rec.get("content_digest"))
    # each dataset as anndata reads it
    check = {}
    for n, p in datasets.items():
        x = ad.read_h5ad(p)
        check[n] = dict(sha256=L.sha256(p), features=[str(f) for f in x.var_names], n_obs=int(x.n_obs),
                        x=h5ad.digest(x.X), layer=h5ad.digest(x.layers[LAYER]), obs_columns=list(x.obs.columns),
                        obs_dtypes={c: str(x.obs[c].dtype) for c in x.obs.columns},
                        uns_keys=sorted(x.uns.keys()), obsm=sorted(x.obsm.keys()))
    write_json(os.path.join(out, "input_check.json"), check)
    # the permuted datasets of the shuffle relation (one order of the datasets, one permutation of each dataset's
    # features, seeded; written compressed, kept beside the run)
    sdir = os.path.join(out, "shuffle_datasets")
    os.makedirs(sdir, exist_ok=True)
    rng = random.Random(f"20261124:{a.organism}:shuffle")
    order = list(handles)
    rng.shuffle(order)
    shuffled, perms = {}, {}
    for d in handles:
        x = ad.read_h5ad(datasets[d])
        perm = list(range(x.n_vars))
        rng.shuffle(perm)
        sp_ = os.path.join(sdir, f"{d}.h5ad")
        x[:, perm].copy().write_h5ad(sp_ + ".tmp.h5ad", compression="gzip")
        os.replace(sp_ + ".tmp.h5ad", sp_)
        shuffled[d] = sp_
        perms[d] = L.sha256(sp_)
    info = dict(organism=a.organism, set=a.set, clone=clone, commit=C.git_commit(clone)[0], python=sys.version.split()[0],
                anndata=__import__("importlib.metadata").metadata.version("anndata"), numpy=np.__version__,
                started=time.strftime("%Y-%m-%d %H:%M:%S"), configs={})
    n_cells = man["n_cells"]
    for cfg in a.configs.split(","):
        fd = CONFIGS[cfg]
        base = os.path.join(out, cfg)
        if os.path.isdir(base):
            shutil.rmtree(base)  # a run is whole: an earlier attempt's outputs of this configuration go
        os.makedirs(base)
        cache = os.path.join(base, "cache")

        def make(dsets, folder, api_=api, order=None):
            names = order or list(dsets)
            return HarmonizeFeatures(f"f1_{cfg}", {k: dsets[k] for k in names}, folder, repo, to_release=TARGET,
                                     final_database=fd, organism=a.organism, snapshot_release=TARGET,
                                     from_release_by_dataset={k: releases[k] for k in names}, api=api_, verbosity=0)
        t = time.time()
        h = make(datasets, cache)
        sec = time.time() - t
        final_db = h.final_database
        log(cfg, "constructed", round(sec, 1), "s; final database", final_db)
        # the answers
        n_all = sum(len(h.feature_names[d]) for d in handles)
        with gzip.open(os.path.join(base, "answers.jsonl.gz.tmp"), "wt") as fo:
            for d in handles:
                tab = h.feature_table(d)
                for i, (name, m) in enumerate(zip(h.feature_names[d], h.matchings[d])):
                    row = tab.iloc[i]
                    st = row["status"]
                    best = sorted({c.gene for c in m.view("best")})
                    got = [row["ensembl_gene"]] if st in ("converted", "converted_no_final_id") else \
                        (best if st == "ambiguous" else [])
                    outcome, et, emro, err = "answer", None, None, None
                    if st == "error":
                        et = m.error.type if m.error is not None else None
                        outcome = "timeout" if et == "SearchBudgetExceededError" else "error"
                        emro, err = [et], (m.error.message if m.error is not None else str(row["reason"]))[:500]
                    a_ = dict(qid=f"f1|{d}|{name}", view="best", outcome=outcome, got=got,
                              apart=sorted(str(x) for x in (m.alternative_locus_genes or [])), status=st, match=m.match,
                              tie=bool(m.tie), from_release=m.from_release, fallback=bool(m.from_release_fallback),
                              final_database=final_db, flags=json.dumps(dict(no_corresponding=m.status in ("not_found", "unknown_version"),
                                                                             no_conversion=m.status == "no_conversion",
                                                                             no_target=bool(m.no_target))),
                              error_type=et, error_mro=emro, error=err, traceback=None, seconds=round(sec / n_all, 6),
                              calls=1, path=None, conv_status=m.status, reason=text(row["reason"]), graph_id=m.graph_id,
                              own=bool(row["own"]), final_id=text(row["final_id"]), dataset=d)
                    fo.write(json.dumps(a_, default=str) + "\n")
        os.replace(os.path.join(base, "answers.jsonl.gz.tmp"), os.path.join(base, "answers.jsonl.gz"))
        cinfo = dict(final_database=final_db, construct_seconds=round(sec, 2), n_features=n_all, outputs={}, unify={})
        # the outputs of each dataset
        for merge in MERGES:
            od = os.path.join(base, "outputs", merge)
            os.makedirs(od)
            for d in handles:
                with __import__("warnings").catch_warnings(record=True) as w:
                    __import__("warnings").simplefilter("always")
                    o, t_in, t_out = h.feature_harmonizer(d, merge=merge)
                hd, rows = describe_dataset_output(o, final_db)
                write_tsv_gz(os.path.join(od, f"{d}.tsv.gz"), hd, rows)
                write_json(os.path.join(od, f"{d}.meta.json"), meta_of(o, dict(n_in=t_in, n_out=t_out,
                                                                                 warnings=[str(x.message) for x in w])))
            log(cfg, merge, "outputs written")
        # the merged objects
        ud = os.path.join(base, "unify")
        os.makedirs(ud)
        for feats in FEATURES:
            for merge in MERGES:
                t = time.time()
                m = h.unify_multiple_anndatas(features=feats, merge=merge)
                us = time.time() - t
                cols, frows, blocks, meta, (oh, orows) = describe_merged(m, handles, n_cells)
                stem = os.path.join(ud, f"{feats}_{merge}")
                write_tsv_gz(stem + ".features.tsv.gz", ["name"] + cols, frows)
                write_tsv_gz(stem + ".blocks.tsv.gz", BLOCK_HEADER, blocks)
                write_tsv_gz(stem + ".obs.tsv.gz", oh, orows)
                meta["seconds"] = round(us, 2)
                write_json(stem + ".meta.json", meta)
                cinfo["unify"][f"{feats}_{merge}"] = dict(seconds=round(us, 2), n_vars=int(m.n_vars), n_obs=int(m.n_obs))
                del m
                log(cfg, feats, merge, "merged", round(us, 1), "s")
        # the result cache
        cache_info = {}
        headers = {d: json.load(open(h.result_cache_path(d)))["header"] for d in handles}
        cache_info["headers"] = headers
        cache_info["graph_content_digest"] = rec.get("content_digest")
        first = {d: [x.to_dict() for x in h.matchings[d]] for d in handles}
        t = time.time()
        h2 = make(datasets, cache)
        cache_info["again"] = dict(seconds=round(time.time() - t, 2), construct_seconds=round(sec, 2),
                                   equal=all([x.to_dict() for x in h2.matchings[d]] == first[d] for d in handles),
                                   files_rewritten=[d for d in handles if os.path.getmtime(h2.result_cache_path(d)) > t])
        del h2
        tam = os.path.join(base, "cache_tampered")
        shutil.copytree(cache, tam)
        victim = handles[0]
        p = h.result_cache_path(victim).replace(cache, tam)
        stored = json.load(open(p))
        stored["header"]["graph_content_digest"] = "0" * 64
        json.dump(stored, open(p, "w"))
        t = time.time()
        h3 = make(datasets, tam)
        back = json.load(open(p))["header"]
        cache_info["changed_digest"] = dict(dataset=victim, seconds=round(time.time() - t, 2),
                                            rewritten_with=back.get("graph_content_digest"),
                                            equal=all([x.to_dict() for x in h3.matchings[d]] == first[d] for d in handles))
        del h3
        calls = []
        orig_build = HF.API.build
        HF.API.build = lambda *x, **k: (calls.append(x), orig_build(*x, **k))[1]
        try:
            t = time.time()
            h4 = make(datasets, cache, api_=None)
            cache_info["without_api"] = dict(seconds=round(time.time() - t, 2), graph_loaded=len(calls),
                                             equal=all([x.to_dict() for x in h4.matchings[d]] == first[d] for d in handles))
            del h4
        except Exception as e:  # noqa: BLE001  a failure here is a finding of the cache, recorded
            cache_info["without_api"] = dict(error=f"{type(e).__name__}: {e}"[:500], graph_loaded=len(calls))
        finally:
            HF.API.build = orig_build
        # the shuffle relation: datasets in another order, features permuted, from a copy of the result cache
        sd = os.path.join(base, "shuffle")
        os.makedirs(sd)
        scache = os.path.join(sd, "cache")
        shutil.copytree(cache, scache)
        t = time.time()
        hs = make(shuffled, scache, order=order)
        sh_sec = time.time() - t
        sheaders = {d: json.load(open(hs.result_cache_path(d)))["header"] for d in handles}
        same_answers = {}
        for d in handles:
            mine = {n: x.to_dict() for n, x in zip(h.feature_names[d], h.matchings[d])}
            theirs = {n: x.to_dict() for n, x in zip(hs.feature_names[d], hs.matchings[d])}
            same_answers[d] = sorted(n for n in mine if mine[n] != theirs.get(n))[:50]
        for merge in MERGES:
            m = hs.unify_multiple_anndatas(features="union", merge=merge)
            cols, frows, blocks, meta, (oh, orows) = describe_merged(m, order, n_cells)
            stem = os.path.join(sd, f"union_{merge}")
            write_tsv_gz(stem + ".features.tsv.gz", ["name"] + cols, frows)
            write_tsv_gz(stem + ".blocks.tsv.gz", BLOCK_HEADER, blocks)
            write_tsv_gz(stem + ".obs.tsv.gz", oh, orows)
            write_json(stem + ".meta.json", meta)
            del m
        cache_info["shuffle"] = dict(order=order, seconds=round(sh_sec, 2), datasets_sha256=perms,
                                     recomputed=[d for d in handles if sheaders[d]["feature_names_sha256"] != headers[d]["feature_names_sha256"]],
                                     answers_differ=same_answers)
        write_json(os.path.join(base, "cache.json"), cache_info)
        for folder in (cache, tam, scache):
            shutil.rmtree(folder)  # rebuildable; the headers are in cache.json
        del h, hs
        cinfo["peak_gb"] = peak_gb()
        info["configs"][cfg] = cinfo
        log(cfg, "done; peak", cinfo["peak_gb"], "GB")
    info.update(finished=time.strftime("%Y-%m-%d %H:%M:%S"), seconds=round(time.time() - t_start, 1), peak_gb=peak_gb())
    write_json(os.path.join(out, "run_f.json"), info)
    print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
