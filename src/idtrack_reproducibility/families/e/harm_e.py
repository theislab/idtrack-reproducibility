"""E5 and E8 for the harmonizer (`idtrack.harmonize.HarmonizeFeatures`, trees of the convert generation, 71690ab on).

Usage: harm_e.py --organism <o> --clone <clone> --graph <graph> --out <run>/<organism>
       (idtrack_py313: anndata 0.13.4, the package's anndata>=0.13; idtrack_dev_env's anndata 0.10.8 cannot read .h5ad)

Datasets are built from the truth (GTF 115 and BioMart 115 of the store, never idtrack): three datasets of 40 cells,
each with 300 current display names of one gene, 200 versioned and 50 unversioned gene ids (overlapping between the
datasets), strings that are nothing, the missing-value spellings of rule `missing_values` (`nan`, `NA`, `NaT`,
`unknown`), and make-unique copies `X-1` / `X.1` of names `X` of the dataset (rule `made_unique_copies`); counts are seeded
integers. Checks (`<out>/harm/report.tsv`, details in `<out>/harm/details.jsonl`):

- arguments (E5): every documented refusal of the constructor and of `unify_multiple_anndatas` raises the
  documented class (`InvalidArgumentError`, `ArgumentTypeError`); valid arguments raise nothing;
- order (E8): `unify_multiple_anndatas(features="union")`, `merge` "sum" and "flag", with the datasets in another
  order and with each dataset's features shuffled: the same features, each with the same values in every cell
  and the same `converted_id`, `matching_ensembl_gene` and `intersection`;
- missing values (rule `missing_values`) and make-unique copies (rule `made_unique_copies`), both judges: with `merge="sum"`, a
  missing-value feature whose counts end in a gene's feature is wrong under the decided judge (71690ab's docs read it
  as a name: open there); a copy `X-1` whose counts end in another gene's feature than X's is wrong under both (safe:
  X's gene or unconverted is right).
"""
import argparse
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.clone)
    from idtrack_reproducibility.families import common as C
    from idtrack_reproducibility.harness import harness_lib as H
    import anndata as ad
    import numpy as np
    import pandas as pd
    import idtrack
    from idtrack.harmonize import HarmonizeFeatures
    from idtrack_reproducibility.families import adapter as F0
    out = os.path.join(os.path.abspath(a.out), "harm")
    os.makedirs(out, exist_ok=True)
    link = F0.link_graph(os.path.abspath(a.graph), os.path.abspath(a.out))
    repo = os.path.dirname(link)
    api = idtrack.API.load(link)
    rng = random.Random(f"20261118:{a.organism}:harm")
    g115 = H.gtf_genes(a.organism, 115, "primary") or {}
    bm = H.biomart_genes(a.organism, 115) or []
    by_name = {}
    for r in bm:
        if r["name"]:
            by_name.setdefault(r["name"], set()).add(r["gene_id"])
    names = sorted(n for n, g in by_name.items() if len(g) == 1 and next(iter(g)) in g115 and "-" not in n and "." not in n
                   and "_" not in n)
    ids = sorted(g115)
    ver = {g: f"{g}.{r['gene_version']}" for g, r in g115.items()}
    pool_n = rng.sample(names, min(900, len(names)))
    pool_i = rng.sample(ids, min(600, len(ids)))
    junk = ["NOT_A_REAL_GENE", "zz_nothing_zz"]
    missing = ["nan", "NA", "NaT", "unknown"]
    datasets = {}
    meta = {}
    for k, h in enumerate(["d1", "d2", "d3"]):
        nn = pool_n[k * 200: k * 200 + 300]
        ii = pool_i[k * 150: k * 150 + 200]
        uu = [x for x in pool_i[k * 150 + 200: k * 150 + 250]]
        mu = [f"{x}-1" for x in nn[:5]] + [f"{x}.1" for x in nn[5:10]]
        feats = nn + [ver[x] for x in ii] + uu + junk + missing + mu
        feats = list(dict.fromkeys(feats))
        X = np.array([[rng.randrange(0, 20) for _ in feats] for _ in range(40)], dtype=np.float32)
        obs = pd.DataFrame({"donor": [f"{h}_{i % 3}" for i in range(40)]}, index=[f"cell{i}" for i in range(40)])
        adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=feats))
        adata.uns["study"] = h
        p = os.path.join(out, f"{h}.h5ad")
        adata.write_h5ad(p)
        datasets[h] = p
        meta[h] = dict(names=nn, mu={f: f.rsplit("-", 1)[0] if f.endswith("-1") else f.rsplit(".", 1)[0] for f in mu})
        # the same dataset with its features shuffled
        perm = list(range(len(feats)))
        rng.shuffle(perm)
        sh = ad.AnnData(X=X[:, perm], obs=obs, var=pd.DataFrame(index=[feats[i] for i in perm]))
        sh.uns["study"] = h
        sh.write_h5ad(os.path.join(out, f"{h}_shuffled.h5ad"))
    rows, details = [], open(os.path.join(out, "details.jsonl"), "w")

    def note(check, ok, detail=None, judge="both"):
        rows.append((check, judge, ok))
        if not ok or detail is not None:
            details.write(json.dumps(dict(check=check, judge=judge, ok=ok, detail=detail), default=str) + "\n")

    def H_(dsets, cache, **kw):
        return HarmonizeFeatures("e", dsets, os.path.join(out, cache), repo, to_release=115, organism=a.organism,
                                 snapshot_release=115, api=api, verbosity=0, **kw)

    # ---- E5: the documented refusals ----
    def expect(check, cls, fn):
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            note(check, type(e).__name__ == cls or cls in [k.__name__ for k in type(e).__mro__],
                 None if cls in [k.__name__ for k in type(e).__mro__] else f"raised {type(e).__name__}: {str(e)[:300]}")
            return
        note(check, False, f"no error; the docs promise {cls}")
    expect("args_empty_dict", "InvalidArgumentError", lambda: H_({}, "c0"))
    expect("args_missing_file", "InvalidArgumentError", lambda: H_({"x": os.path.join(out, "nope.h5ad")}, "c0"))
    expect("args_verbosity", "InvalidArgumentError",
           lambda: HarmonizeFeatures("e", datasets, os.path.join(out, "c0"), repo, to_release=115, organism=a.organism,
                                     snapshot_release=115, api=api, verbosity=5))
    expect("args_to_release_str", "ArgumentTypeError",
           lambda: HarmonizeFeatures("e", datasets, os.path.join(out, "c0"), repo, to_release="115", organism=a.organism,
                                     snapshot_release=115, api=api, verbosity=0))
    expect("args_to_release_after_snapshot", "InvalidArgumentError",
           lambda: HarmonizeFeatures("e", datasets, os.path.join(out, "c0"), repo, to_release=116, organism=a.organism,
                                     snapshot_release=115, api=api, verbosity=0))
    expect("args_from_release_unknown_dataset", "InvalidArgumentError", lambda: H_(datasets, "c0", from_release_by_dataset={"zz": 100}))
    expect("args_from_database_unknown_dataset", "InvalidArgumentError",
           lambda: H_(datasets, "c0", from_database_by_dataset={"zz": "UniProtKB Gene Name"}))
    rep = ad.AnnData(X=np.ones((2, 2), dtype=np.float32), var=pd.DataFrame(index=["A", "B"]))
    rep.var_names = ["TP53", "TP53"]
    rp = os.path.join(out, "repeated.h5ad")
    rep.write_h5ad(rp)
    expect("args_repeated_features", "InvalidArgumentError", lambda: H_({"r": rp}, "c0"))
    h = H_(datasets, "c1")
    expect("unify_features", "InvalidArgumentError", lambda: h.unify_multiple_anndatas(features="x"))
    expect("unify_obs_join", "InvalidArgumentError", lambda: h.unify_multiple_anndatas(obs_join="x"))
    expect("unify_merge", "InvalidArgumentError", lambda: h.unify_multiple_anndatas(merge="x"))
    expect("unify_obs_keep", "InvalidArgumentError", lambda: h.unify_multiple_anndatas(obs_columns_to_keep=["no_such_column"]))
    expect("unify_handle_key", "InvalidArgumentError", lambda: h.unify_multiple_anndatas(handle_anndata_key="donor"))
    expect("plan_merge", "InvalidArgumentError", lambda: h.output_plan("d1", merge="x"))

    # ---- E8: order of datasets and of features ----
    def table(m):
        X = m.X.toarray() if hasattr(m.X, "toarray") else np.asarray(m.X)
        cols = {}
        for j, f in enumerate(m.var_names):
            cols[f] = tuple(float(v) for v in X[:, j])
        cells = list(m.obs_names)
        order = np.argsort(cells)
        per = {f: tuple(v[i] for i in order) for f, v in cols.items()}
        v = m.var
        keep = [c for c in ("converted_id", "matching_ensembl_gene", "intersection") if c in v.columns]
        return per, {f: tuple(str(v.at[f, c]) for c in keep) for f in m.var_names}
    for merge in ("sum", "flag"):
        ref = table(h.unify_multiple_anndatas(merge=merge))
        h2 = H_({k: datasets[k] for k in ("d3", "d1", "d2")}, "c2")
        other = table(h2.unify_multiple_anndatas(merge=merge))
        note(f"order_datasets_{merge}_features", set(ref[0]) == set(other[0]), sorted(set(ref[0]) ^ set(other[0]))[:20] or None)
        diff = [f for f in ref[0] if f in other[0] and (ref[0][f] != other[0][f] or ref[1][f] != other[1][f])]
        note(f"order_datasets_{merge}_values", not diff, diff[:20] or None)
        h3 = H_({k: os.path.join(out, f"{k}_shuffled.h5ad") for k in ("d1", "d2", "d3")}, "c3")
        other = table(h3.unify_multiple_anndatas(merge=merge))
        note(f"order_features_{merge}_features", set(ref[0]) == set(other[0]), sorted(set(ref[0]) ^ set(other[0]))[:20] or None)
        diff = [f for f in ref[0] if f in other[0] and (ref[0][f] != other[0][f] or ref[1][f] != other[1][f])]
        note(f"order_features_{merge}_values", not diff, diff[:20] or None)

    # ---- rules `missing_values` and `made_unique_copies` in the harmonizer ----
    for hd in ("d1", "d2", "d3"):
        m = h.feature_harmonizer(hd, merge="sum")[0]
        v = m.var
        q = v["Query ID"].astype(str)
        for tok in missing:
            hit = [f for f in m.var_names if tok in str(q[f]).split(";")]
            st = [str(v.at[f, "idtrack_status"]) for f in hit]
            into_gene = any(s in ("converted", "converted_no_final_id", "merged", "merged_no_final_id") for s in st)
            rows.append((f"missing_{tok}", "now", None))
            note(f"missing_{tok}", not into_gene, dict(dataset=hd, features=hit, status=st), judge="decided")
        names_genes = {x: ver[next(iter(by_name[x]))] for x in meta[hd]["names"]}
        for f, x in meta[hd]["mu"].items():
            hit = [g for g in m.var_names if f in str(q[g]).split(";")]
            genes = {str(v.at[g, "matching_ensembl_gene"]) for g in hit if str(v.at[g, "idtrack_status"]).startswith(("converted", "merged"))}
            ok = genes <= {names_genes.get(x)}
            note("makeunique_copy", ok, None if ok else dict(dataset=hd, feature=f, base=x, base_gene=names_genes.get(x), went_to=sorted(genes)))
    details.close()
    import collections
    agg = collections.defaultdict(lambda: [0, 0, 0])
    for check, judge, ok in rows:
        k = (check.split("_")[0] if check.startswith(("args", "unify", "plan")) else check, judge)
        agg[k][0] += 1
        if ok is True:
            agg[k][1] += 1
        elif ok is False:
            agg[k][2] += 1
    with open(os.path.join(out, "report.tsv"), "w") as fo:
        fo.write("check\tjudge\torganism\tn\tholds\tviolations\n")
        for (c, j), (n, ok, bad) in sorted(agg.items()):
            fo.write(f"{c}\t{j}\t{a.organism}\t{n}\t{ok}\t{bad}\n")
    json.dump(dict(organism=a.organism, clone=a.clone, commit=C.git_commit(a.clone)[0], finished=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(out, "info.json"), "w"))
    print(open(os.path.join(out, "report.tsv")).read())


if __name__ == "__main__":
    main()
