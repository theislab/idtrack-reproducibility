"""F1's judge: every feature's assignment by the common scorer, then the matrices, the metadata, the order of features
and datasets and the result cache by F's own checks (docs/families/f.md, "Scoring").

Usage: score_f.py --organism <o> --run <families>/f/runs/<tree>/<o> [--set f1_v1] [--configs gene,symbol]

1. The run's answers (`<config>/answers.jsonl.gz`, written by run_f in the harmonizer's environment, which has no
   pyarrow) become the runner's answer files (`<config>/answers/shard_0_of_1.0.parquet`, `load_0_of_1.json`), and the
   common scorer judges them (`families/scorer.py --family f --sets <the configuration's sets> --scores-to
   <run>/scores_<config>`): one row of the record per feature and dataset, `report.tsv` per stratum.
   Configuration `gene` is judged by the set (and its revision sets); `symbol` also by the `symbol` judge's lines
   (`truth_f.py`, judge `symbol`).
2. F's checks, written as summary rows (`f_checks.tsv`: one row per check, configuration, merge and scope with the
   items checked, those that hold and those that do not, and an example; every violation in `f_violations.jsonl.gz`):
   - `input_*`: each dataset as anndata read it in the run is what the generator wrote (features, X, the layer);
   - `plan_*` (no truth): each harmonized dataset is what the documented rules make of the harmonizer's own
     assignment of its features (`truth_f.plan`): the output features in order, their inputs and status, their
     values in X and in the layer (the sums of their inputs), the input `var` columns, and `obs`, `obsm`, `obsp`, `uns`
     kept as they are, the layers, `varm` subset or left out where features were summed;
   - `truth_*`: on the decided features (`truth_f.decided`), the outputs the truth implies: each expected output
     feature is in the harmonized dataset with the same inputs, status and values, judged where every input on its
     gene (expected or as the harmonizer assigned it) and every input of the output holding it is decided, else apart
     (`truth_feature_apart`);
   - `unify_*` (no truth): the merged objects of `unify_multiple_anndatas` (union, intersection; sum, flag): their
     features the union or intersection of the datasets' outputs (with the documented renaming), each dataset's block
     of a feature its output's values and absent features zero, `intersection`, `obs` (every input column with each
     cell's value, missing where a dataset lacks it, the documented dtypes, the dataset column, the cell names),
     `obsm`, `obsp` and the layers block by block, `uns` the entries equal in every dataset;
   - `unify_truth_*`: every output feature that matched the truth in a dataset (`truth_feature`) is in the union, and
     in the intersection where it matched in every dataset;
   - `shuffle_*` (E8's relation): the union of the datasets in another order with permuted features: the same features,
     each with the same inputs and values in every dataset's block, the same merged columns and the same answer per
     feature;
   - `cache_*` (rule `cache_by_digest`): the result cache's headers name the graph's content digest and the rules'
     version and nothing of the package's version; the constructor again on it gives the same conversions, faster; a
     stored digest changed is computed again and written with the graph's; without `api` the graph is not loaded;
     changed feature names are converted again.
3. `f_counts.tsv`, count conservation by class: per configuration, merge and stratum, the decided input features and
   their counts (the sums of their X columns) `kept` (in the output feature the truth implies), `moved` (in another),
   `dropped` (in none), and the undecided ones apart; the expected is all kept.
4. `f_assign.tsv`: the assignment rate per configuration, form and stratum, with Wilson's and the cluster bootstrap's
   intervals (unit: the source gene's lineage or the name).
"""
import argparse
import collections
import glob
import gzip
import json
import os
import subprocess
import sys

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as R
from idtrack_reproducibility.families.f import h5ad
from idtrack_reproducibility.families.f import truth_f as T

CONFIGS = ("gene", "symbol")
MERGES = ("sum", "flag")
LAYER = "unspliced"
INPUT_VAR = ("feature_type", "highly_variable", "mean_input")


def sets_for(org, config, set_name="f1_v1"):
    """The sets a configuration is judged by, in publication order: the set and its revision sets; `symbol` also the
    lines of the `symbol` judge."""
    out = []
    for s in L.set_names("f", org):
        m = L.manifest("f", org, s)
        if s == set_name or set_name in (m.get("revises") or []):
            if m.get("judge") == "symbol" and config != "symbol":
                continue
            out.append(s)
    return out


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, x.rstrip("\n").split("\t"))) for x in fh]


def to_answer_files(run_cfg):
    """`answers.jsonl.gz` -> the runner's answer Parquet and load record (step 1)."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    from idtrack_reproducibility.families import runner
    sch = runner.answer_schema()
    rows = []
    with gzip.open(os.path.join(run_cfg, "answers.jsonl.gz"), "rt") as fh:
        for x in fh:
            r = json.loads(x)
            r["h"] = C.crc(r["qid"])
            rows.append({f.name: r.get(f.name) for f in sch})
    d = os.path.join(run_cfg, "answers")
    os.makedirs(d, exist_ok=True)
    for old in glob.glob(os.path.join(d, "*.parquet")):
        os.remove(old)
    pq.write_table(pa.Table.from_pylist(rows, schema=sch), os.path.join(d, "shard_0_of_1.0.parquet"), compression="zstd")
    load = json.load(open(os.path.join(os.path.dirname(run_cfg), "load.json")))
    load.update(shard="0/1")
    json.dump(load, open(os.path.join(d, "load_0_of_1.json"), "w"), indent=1)
    return len(rows)


class Checks:
    def __init__(self):
        self.rows = collections.OrderedDict()
        self.violations = []

    def add(self, check, config, merge, scope, ok, detail=None, n=1):
        k = (check, config, merge, scope)
        r = self.rows.setdefault(k, [0, 0, 0, ""])
        r[0] += n
        if ok:
            r[1] += n
        else:
            r[2] += n
            if not r[3]:
                r[3] = json.dumps(detail, default=str)[:300]
            self.violations.append(dict(check=check, config=config, merge=merge, scope=scope, detail=detail))

    def write(self, out):
        with open(os.path.join(out, "f_checks.tsv.tmp"), "w") as fo:
            fo.write("check\tconfig\tmerge\tscope\tn\tholds\tviolations\texample\n")
            for (c, cfg, m, s), (n, ok, bad, ex) in self.rows.items():
                fo.write(f"{c}\t{cfg}\t{m}\t{s}\t{n}\t{ok}\t{bad}\t{ex.replace(chr(9), ' ')}\n")
        os.replace(os.path.join(out, "f_checks.tsv.tmp"), os.path.join(out, "f_checks.tsv"))
        with gzip.open(os.path.join(out, "f_violations.jsonl.gz.tmp"), "wt") as fo:
            for v in self.violations:
                fo.write(json.dumps(v, default=str) + "\n")
        os.replace(os.path.join(out, "f_violations.jsonl.gz.tmp"), os.path.join(out, "f_violations.jsonl.gz"))


class Inputs:
    """The datasets as the generator wrote them (h5py, `h5ad.read`), with the digests of sums of their columns."""

    def __init__(self, path):
        import scipy.sparse as sp
        d = h5ad.read(path)
        self.features = list(d["var"]["index"])
        self.pos = {f: i for i, f in enumerate(self.features)}
        self.X = sp.csc_matrix(d["X"])
        self.U = sp.csc_matrix(d["layers"][LAYER])
        self.n_cells = self.X.shape[0]
        self.obs = d["obs"]
        self.obs_values = h5ad.frame_values(d["obs"])
        self.var_values = h5ad.frame_values(d["var"])
        self.obsm = {k: h5ad.digest(v) for k, v in d["obsm"].items()}
        self.obsm_raw = d["obsm"]
        self.obsp = {k: h5ad.digest(v) for k, v in d["obsp"].items()}
        self.obsp_raw = d["obsp"]
        self.varm = d["varm"]
        self.uns = h5ad.plain(d["uns"])
        self.totals = {f: float(self.X[:, i].sum()) for f, i in self.pos.items()}

    def sum_digest(self, members, layer=False):
        import hashlib

        import numpy as np
        M = self.U if layer else self.X
        v = np.zeros(self.n_cells, dtype=np.float64)
        for f in members:
            j = self.pos[f]
            s, e = M.indptr[j], M.indptr[j + 1]
            v[M.indices[s:e]] += M.data[s:e]
        return hashlib.sha256(v.tobytes()).hexdigest()[:16]

    def var_expected(self, members):
        """The input `var` values of an output feature: kept where all its inputs agree, missing otherwise."""
        out = {}
        for c in INPUT_VAR:
            vals = {self.var_values[c][self.pos[f]] for f in members}
            out[c] = vals.pop() if len(vals) == 1 else ""
        return out


def status_class(st):
    st = str(st)
    if st.startswith("merged"):
        return "merged"
    if st == "unconverted_merge_partner":
        return "merge_partner"
    if st.startswith("unconverted"):
        return "unconverted"
    return "converted"


def kinds_ok(inp, merged, partial):
    """The documented dtype of a merged obs column (unify_multiple_anndatas): text becomes categorical, integers and
    booleans a column some dataset lacks become nullable, numbers stay numbers."""
    i, m = str(inp), str(merged)
    if i in ("category",):
        return m == "category"
    if i in ("object", "str", "string") or i.startswith("string"):
        return m == "category"
    if i.startswith("int"):
        return m.startswith("Int") if partial else m.startswith("int")
    if i == "bool":
        return m == "boolean" if partial else m == "bool"
    if i.startswith("float"):
        return m.startswith("float")
    return i == m


def score_config(a, cfg, ck, inputs, lines_by_q, man, counts_rows):
    run = a.run
    base = os.path.join(run, cfg)
    answers = {}
    with gzip.open(os.path.join(base, "answers.jsonl.gz"), "rt") as fh:
        for x in fh:
            r = json.loads(x)
            answers[r["qid"]] = r
    handles = sorted(man["datasets"])
    expected_names = {m: {} for m in MERGES}  # merge -> dataset -> [(name, status)] of the decided expected outputs
    for merge in MERGES:
        for d in handles:
            inp = inputs[d]
            outs = read_tsv_gz(os.path.join(base, "outputs", merge, f"{d}.tsv.gz"))
            meta = json.load(open(os.path.join(base, "outputs", merge, f"{d}.meta.json")))
            # ---- plan_*: the documented rules applied to the harmonizer's own assignment ----
            own_in = []
            for f in inp.features:
                ans = answers[f"f1|{d}|{f}"]
                st = ans["status"]
                gene = ans["got"][0] if st in T.CONVERTED and ans["got"] else None
                own_in.append((f, gene, bool(ans.get("own")), st))
            exp = T.plan(own_in, merge)
            ck.add("plan_names_in_order", cfg, merge, d, [o["name"] for o in exp] == [o["name"] for o in outs],
                   dict(expected=[o["name"] for o in exp][:5], got=[o["name"] for o in outs][:5]))
            got_by = {o["name"]: o for o in outs}
            for o in exp:
                g = got_by.get(o["name"])
                if g is None:
                    ck.add("plan_feature", cfg, merge, d, False, dict(missing=o["name"]))
                    continue
                same = (g["query_ids"] == ";".join(o["members"]) and g["status"] == o["status"] and
                        (g["gene"] or None) == o["gene"])
                ck.add("plan_feature", cfg, merge, d, same, dict(name=o["name"], expected=o, got=g))
                ck.add("plan_values_x", cfg, merge, d, g["x_digest"] == inp.sum_digest(o["members"]), dict(name=o["name"]))
                ck.add("plan_values_layer", cfg, merge, d, g["u_digest"] == inp.sum_digest(o["members"], True), dict(name=o["name"]))
                ve = inp.var_expected(o["members"])
                ck.add("plan_var_inputs", cfg, merge, d, all(g.get(f"var:{c}", "") == ve[c] for c in INPUT_VAR),
                       dict(name=o["name"], expected=ve, got={c: g.get(f"var:{c}") for c in INPUT_VAR}))
            members_out = collections.Counter(m for g in outs for m in (g["query_ids"].split(";") if g["query_ids"] else []))
            ck.add("plan_every_input_once", cfg, merge, d, all(members_out.get(f) == 1 for f in inp.features),
                   dict(not_once=[f for f in inp.features if members_out.get(f) != 1][:10]))
            # metadata of the dataset (kept as they are)
            ck.add("plan_obs", cfg, merge, d, meta["obs_names"] == list(inp.obs["index"]) and
                   all(meta["obs_values"].get(c) == v for c, v in inp.obs_values.items()) and
                   set(meta["obs_columns"]) == set(inp.obs_values), dict(columns=meta["obs_columns"]))
            ck.add("plan_obs_dtypes", cfg, merge, d, meta["obs_dtypes"] == a.input_check[d].get("obs_dtypes", meta["obs_dtypes"]),
                   dict(got=meta["obs_dtypes"], input=a.input_check[d].get("obs_dtypes")))
            ck.add("plan_obsm", cfg, merge, d, meta["obsm"] == inp.obsm, dict(got=meta["obsm"], input=inp.obsm))
            ck.add("plan_obsp", cfg, merge, d, meta["obsp"] == inp.obsp, dict(got=meta["obsp"], input=inp.obsp))
            ck.add("plan_uns", cfg, merge, d, meta["uns"] == inp.uns, dict(got=meta["uns"], input=inp.uns))
            ck.add("plan_layers", cfg, merge, d, meta["layers"] == [LAYER], dict(got=meta["layers"]))
            summed = any(len(o["members"]) > 1 for o in exp)
            if summed:
                ck.add("plan_varm_left_out_when_summed", cfg, merge, d, meta["varm"] == {}, dict(got=meta["varm"]))
            else:
                import numpy as np
                rows = [inp.pos[o["members"][0]] for o in exp]
                want = {k: h5ad.digest(np.asarray(v)[rows]) for k, v in inp.varm.items()}
                ck.add("plan_varm_subset", cfg, merge, d, meta["varm"] == want, dict(got=meta["varm"], expected=want))
            # ---- truth_*: the outputs the truth implies on the decided features ----
            dec, tin = {}, []
            for f in inp.features:
                q = lines_by_q[f"f1|{d}|{f}"]
                kind, g = T.decided(q)
                dec[f] = (kind, g)
                own = T.own_by_truth(f, g, cfg) if kind == "gene" else False
                tin.append((f, g if kind == "gene" else None, own, "converted" if kind == "gene" else "kept"))
            texp = T.plan(tin, merge)
            undecided = {f for f, (k, _g) in dec.items() if k is None}
            dest_exp = {}
            for o in texp:
                for m_ in o["members"]:
                    dest_exp[m_] = o
            dest_got = {}
            for g in outs:
                for m_ in (g["query_ids"].split(";") if g["query_ids"] else []):
                    dest_got[m_] = g
            # the inputs the harmonizer put on each gene: an expected output is judged only where every input on its
            # gene (as expected and as assigned) and every input of the output holding it is decided
            on_gene = collections.defaultdict(set)
            for f, gene, _own, st in own_in:
                if gene is not None:
                    on_gene[gene].add(f)
            names, apart = [], set()
            for o in texp:
                if any(m_ in undecided for m_ in o["members"]) or o["status"] == "undetermined":
                    apart |= set(o["members"])
                    continue
                g = got_by.get(o["name"])
                around = set(on_gene.get(o["gene"], set())) if o["gene"] else set()
                for m_ in o["members"]:
                    h = dest_got.get(m_)
                    around |= set(h["query_ids"].split(";")) if h and h["query_ids"] else set()
                    if h and h["gene"]:
                        around |= on_gene.get(h["gene"], set())
                if around & undecided:
                    ck.add("truth_feature_apart", cfg, merge, d, True)
                    apart |= set(o["members"])
                    continue
                ok = g is not None and set(g["query_ids"].split(";")) == set(o["members"]) and \
                    status_class(g["status"]) == status_class(o["status"]) and g["x_digest"] == inp.sum_digest(o["members"])
                ck.add("truth_feature", cfg, merge, d, ok, dict(name=o["name"], expected=o, got=g))
                if ok:
                    names.append((o["name"], o["status"]))
            expected_names[merge][d] = names
            # count conservation by class
            for f in inp.features:
                q = lines_by_q[f"f1|{d}|{f}"]
                k = (cfg, merge, q["stratum"])
                c = counts_rows.setdefault(k, collections.Counter())
                tot = inp.totals[f]
                o = dest_exp.get(f)
                if f in undecided or f in apart or o is None or o["status"] == "undetermined":
                    c["apart"] += 1
                    c["apart_counts"] += tot
                    continue
                g = dest_got.get(f)
                what = "dropped" if g is None else ("kept" if g["name"] == o["name"] else "moved")
                c["features"] += 1
                c["counts"] += tot
                c[what] += 1
                c[f"{what}_counts"] += tot
    # ---- unify_* and unify_truth_* ----
    for feats in ("union", "intersect"):
        for merge in MERGES:
            stem = os.path.join(base, "unify", f"{feats}_{merge}")
            meta = json.load(open(stem + ".meta.json"))
            frows = read_tsv_gz(stem + ".features.tsv.gz")
            blocks = read_tsv_gz(stem + ".blocks.tsv.gz")
            per = {d: [(o["name"], o["status"]) for o in read_tsv_gz(os.path.join(base, "outputs", merge, f"{d}.tsv.gz"))]
                   for d in handles}
            renamed, names = T.merged_names(per, feats)
            got_names = [r["name"] for r in frows]
            ck.add("unify_features", cfg, merge, feats, set(got_names) == names and len(got_names) == len(set(got_names)),
                   dict(missing=sorted(names - set(got_names))[:10], extra=sorted(set(got_names) - names)[:10]))
            outs = {d: {o["name"]: o for o in read_tsv_gz(os.path.join(base, "outputs", merge, f"{d}.tsv.gz"))} for d in handles}
            seen = set()
            for b in blocks:
                d = b["dataset"]
                orig = renamed[d].get(b["name"])
                o = outs[d].get(orig) if orig is not None else None
                seen.add((b["name"], d))
                ok = o is not None and o["query_ids"] == b["query_ids"] and o["x_digest"] == b["x_digest"] and \
                    o["u_digest"] == b["u_digest"] and all(o.get(f"var:{c}") == b.get(f"var:{c}") for c in INPUT_VAR)
                ck.add("unify_block", cfg, merge, feats, ok, dict(name=b["name"], dataset=d, block=b, output=o))
            missing = [(n, d) for d in handles for n in renamed[d] if n in names and (n, d) not in seen]
            ck.add("unify_every_block", cfg, merge, feats, not missing, dict(missing=missing[:10]))
            ck.add("unify_absent_blocks_zero", cfg, merge, feats, meta["nonzero_blocks_of_absent_features"] == 0,
                   dict(nonzero=meta["nonzero_blocks_of_absent_features"]))
            if feats == "union":
                col = {r["name"]: r.get("intersection") for r in frows}
                want = {n: "1" if all(n in renamed[d] for d in handles) else "0" for n in names}
                bad = [n for n in names if col.get(n) != want[n]]
                ck.add("unify_intersection_column", cfg, merge, feats, not bad, dict(bad=bad[:10]))
            # obs
            obs = read_tsv_gz(stem + ".obs.tsv.gz")
            by_cell = {r["obs_name"]: r for r in obs}
            all_cols = []
            for d in handles:
                all_cols += [c for c in inputs[d].obs_values if c not in all_cols]
            ok_names = [r["obs_name"] for r in obs] == [f"{d}_{c}" for d in handles for c in inputs[d].obs["index"]]
            ck.add("unify_obs_names", cfg, merge, feats, ok_names, dict(first=[r["obs_name"] for r in obs[:3]]))
            bad = []
            for d in handles:
                inp = inputs[d]
                for i, c in enumerate(inp.obs["index"]):
                    r = by_cell.get(f"{d}_{c}")
                    if r is None or r.get("handle_anndata") != d:
                        bad.append((d, c, "cell or dataset column"))
                        continue
                    for col in all_cols:
                        want = inp.obs_values[col][i] if col in inp.obs_values else ""
                        if r.get(col) != want:
                            bad.append((d, c, col, want, r.get(col)))
            ck.add("unify_obs_values", cfg, merge, feats, not bad, dict(bad=bad[:5]), n=1)
            dt = meta["obs_dtypes"]
            dbad = []
            for col in all_cols:
                partial = not all(col in inputs[d].obs_values for d in handles)
                src = next(a.input_check[d]["obs_dtypes"][col] for d in handles if col in inputs[d].obs_values
                           and "obs_dtypes" in a.input_check[d])
                if not kinds_ok(src, dt.get(col), partial):
                    dbad.append((col, src, dt.get(col), partial))
            ck.add("unify_obs_dtypes", cfg, merge, feats, not dbad and set(dt) == set(all_cols) | {"handle_anndata"},
                   dict(bad=dbad, columns=sorted(dt)))
            ck.add("unify_obsm_blocks", cfg, merge, feats, all(meta["obsm_blocks"][d] == inputs[d].obsm for d in handles),
                   dict(got=meta["obsm_blocks"][handles[0]], input=inputs[handles[0]].obsm))
            ck.add("unify_obsp_blocks", cfg, merge, feats, all(meta["obsp_blocks"][d] == inputs[d].obsp for d in handles)
                   and meta["obsp_off_block_nnz"] == 0, dict(off_block=meta["obsp_off_block_nnz"]))
            ck.add("unify_layers", cfg, merge, feats, meta["layers"] == [LAYER], dict(got=meta["layers"]))
            want_uns = T.uns_kept([inputs[d].uns for d in handles])
            ck.add("unify_uns", cfg, merge, feats, meta["uns"] == want_uns, dict(got=meta["uns"], expected=want_uns))
            # the merge of the decided features: every output that matched the truth in a dataset is in the union, and
            # in the intersection where it matched in every dataset
            # (under the name the merge gives it from every dataset's outputs: `renamed`, the documented renaming)
            merged_of = {d: {n: k for k, n in renamed[d].items()} for d in handles}
            matched = {d: {merged_of[d][n] for n, _ in expected_names[merge][d] if n in merged_of[d]} for d in handles}
            want = set().union(*matched.values()) if feats == "union" else set.intersection(*matched.values())
            gotn = set(got_names)
            for n in sorted(want):
                ck.add(f"unify_truth_{feats}", cfg, merge, feats, n in gotn, dict(missing=n))
    # ---- shuffle_* (E8's relation) ----
    cache = json.load(open(os.path.join(base, "cache.json")))
    for merge in MERGES:
        a_stem = os.path.join(base, "unify", f"union_{merge}")
        b_stem = os.path.join(base, "shuffle", f"union_{merge}")
        fa = {r["name"]: r for r in read_tsv_gz(a_stem + ".features.tsv.gz")}
        fb = {r["name"]: r for r in read_tsv_gz(b_stem + ".features.tsv.gz")}
        ck.add("shuffle_features", cfg, merge, "union", set(fa) == set(fb), dict(only_one=sorted(set(fa) ^ set(fb))[:10]))
        keys = [k for k in ("converted_id", "matching_ensembl_gene", "intersection") if any(k in r for r in fa.values())]
        bad = [n for n in fa if n in fb and any(fa[n].get(k) != fb[n].get(k) for k in keys)]
        ck.add("shuffle_merged_columns", cfg, merge, "union", not bad, dict(bad=bad[:10], columns=keys))
        ba = {(r["name"], r["dataset"]): r for r in read_tsv_gz(a_stem + ".blocks.tsv.gz")}
        bb = {(r["name"], r["dataset"]): r for r in read_tsv_gz(b_stem + ".blocks.tsv.gz")}
        bad = [k for k in ba if k not in bb or set(ba[k]["query_ids"].split(";")) != set(bb[k]["query_ids"].split(";"))
               or ba[k]["x_digest"] != bb[k]["x_digest"] or ba[k]["u_digest"] != bb[k]["u_digest"]
               or ba[k]["status"] != bb[k]["status"]]
        ck.add("shuffle_blocks", cfg, merge, "union", not bad and len(ba) == len(bb), dict(bad=bad[:10], n=(len(ba), len(bb))))
        oa = {r["obs_name"]: r for r in read_tsv_gz(a_stem + ".obs.tsv.gz")}
        ob = {r["obs_name"]: r for r in read_tsv_gz(b_stem + ".obs.tsv.gz")}
        ck.add("shuffle_obs", cfg, merge, "union", oa == ob, dict(n=(len(oa), len(ob))))
    diff = {d: v for d, v in cache["shuffle"]["answers_differ"].items() if v}
    ck.add("shuffle_answers", cfg, "-", "union", not diff, dict(differ=diff))
    # ---- cache_* (rule `cache_by_digest`) ----
    allowed = {"format", "rules_version", "organism", "final_database", "target_ensembl_release", "source_release",
               "source_database", "graph_snapshot", "graph_content_digest", "feature_names_sha256"}
    for d, hd in cache["headers"].items():
        ck.add("cache_header_digest", cfg, "-", d, hd.get("graph_content_digest") == cache["graph_content_digest"], hd)
        ck.add("cache_header_rules_version", cfg, "-", d, hd.get("rules_version") is not None, hd)
        ck.add("cache_header_no_package_version", cfg, "-", d, set(hd) <= allowed, sorted(set(hd) - allowed))
    ag = cache["again"]
    ck.add("cache_again_equal", cfg, "-", "all", ag["equal"] and not ag["files_rewritten"], ag)
    ck.add("cache_again_faster", cfg, "-", "all", ag["seconds"] < ag["construct_seconds"], ag)
    cd = cache["changed_digest"]
    ck.add("cache_changed_digest_recomputed", cfg, "-", cd["dataset"], cd["rewritten_with"] == cache["graph_content_digest"]
           and cd["equal"], cd)
    wa = cache["without_api"]
    ck.add("cache_without_api_no_graph_load", cfg, "-", "all", wa.get("graph_loaded") == 0 and wa.get("equal") is True, wa)
    sh = cache["shuffle"]
    ck.add("cache_changed_features_recomputed", cfg, "-", "all", sorted(sh["recomputed"]) == handles,
           dict(recomputed=len(sh["recomputed"]), of=len(handles)))


def assign_table(a, out):
    """`f_assign.tsv`: per configuration, form and stratum, the assignment rate with Wilson's and the cluster
    bootstrap's intervals (unit: the line's unit)."""
    import numpy as np
    rows = []
    for cfg in a.configs:
        folder = os.path.join(a.run, f"scores_{cfg}")
        df = R.read_outcomes(folder, columns=["qid", "kind", "stratum", "unit", "outcome"])
        for key, g in [((k, s), x) for (k, s), x in df.groupby(["kind", "stratum"])] + \
                [(("*", s), x) for s, x in df.groupby("stratum")] + [(("*", "*"), df)]:
            c = g["outcome"].value_counts().to_dict()
            right = c.get("right", 0) + c.get("expected_error", 0)
            den = right + c.get("wrong", 0) + c.get("error", 0)
            lo, hi = R.wilson(right, den)
            u = g.assign(_r=g["outcome"].isin(["right", "expected_error"]).astype(int),
                         _d=g["outcome"].isin(["right", "expected_error", "wrong", "error"]).astype(int))
            per = u.groupby("unit")[["_r", "_d"]].sum()
            blo, bhi = R.cluster_bootstrap(per["_r"].to_numpy(), per["_d"].to_numpy(), f"f|{cfg}|{key[0]}|{key[1]}")
            wrong_units = g.loc[g["outcome"] == "wrong", "unit"].nunique()
            rows.append([cfg, key[0], key[1], len(g), right, c.get("wrong", 0), c.get("undetermined", 0), c.get("timeout", 0),
                         c.get("error", 0), f"{right / den:.4f}" if den else "nan", f"{lo:.4f}", f"{hi:.4f}", f"{blo:.4f}",
                         f"{bhi:.4f}", g["unit"].nunique(), wrong_units])
        _ = np
    with open(os.path.join(out, "f_assign.tsv.tmp"), "w") as fo:
        fo.write("config\tform\tstratum\tn\tright\twrong\tundetermined\ttimeout\terror\trate\twilson_low\twilson_high\t"
                 "boot_low\tboot_high\tunits\tunits_wrong\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    os.replace(os.path.join(out, "f_assign.tsv.tmp"), os.path.join(out, "f_assign.tsv"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", default="f1_v1")
    ap.add_argument("--configs", default="gene,symbol")
    ap.add_argument("--checks-only", action="store_true", help="F's checks only (the common scorings are made)")
    a = ap.parse_args()
    a.run = os.path.abspath(a.run)
    a.configs = a.configs.split(",")
    log = lambda *x: print(*x, file=sys.stderr, flush=True)  # noqa: E731
    man = L.manifest("f", a.organism, a.set)
    if not a.checks_only:
        for cfg in a.configs:
            n = to_answer_files(os.path.join(a.run, cfg))
            sets = sets_for(a.organism, cfg, a.set)
            cmd = [sys.executable, "-m", "idtrack_reproducibility.families.scorer", "--family", "f", "--organism", a.organism,
                   "--run", os.path.join(a.run, cfg), "--sets", ",".join(sets), "--scores-to", os.path.join(a.run, f"scores_{cfg}")]
            log(cfg, n, "answers; scoring:", " ".join(cmd[1:]))
            subprocess.run(cmd, check=True)
    a.input_check = json.load(open(os.path.join(a.run, "input_check.json")))
    ck = Checks()
    inputs = {}
    fdir = C.family_dir("f")
    for d, info in sorted(man["datasets"].items()):
        p = os.path.join(fdir, info["file"])
        inputs[d] = Inputs(p)
        ic = a.input_check[d]
        ck.add("input_file", "-", "-", d, ic["sha256"] == info["sha256"] == L.sha256(p), dict(run=ic["sha256"], manifest=info["sha256"]))
        ck.add("input_features", "-", "-", d, ic["features"] == inputs[d].features)
        import scipy.sparse as sp
        ck.add("input_x", "-", "-", d, ic["x"] == h5ad.digest(sp.csr_matrix(inputs[d].X)))
        ck.add("input_layer", "-", "-", d, ic["layer"] == h5ad.digest(sp.csr_matrix(inputs[d].U)))
    counts_rows = {}
    for cfg in a.configs:
        lines, _, _ = L.load("f", a.organism, sets_for(a.organism, cfg, a.set))
        lines_by_q = {q["qid"]: q for q in lines}
        log(cfg, "checks")
        score_config(a, cfg, ck, inputs, lines_by_q, man, counts_rows)
    ck.write(a.run)
    with open(os.path.join(a.run, "f_counts.tsv.tmp"), "w") as fo:
        cols = ["features", "counts", "kept", "moved", "dropped", "kept_counts", "moved_counts", "dropped_counts", "apart",
                "apart_counts"]
        fo.write("config\tmerge\tstratum\t" + "\t".join(cols) + "\texpected\n")
        for (cfg, merge, st), c in sorted(counts_rows.items()):
            fo.write(f"{cfg}\t{merge}\t{st}\t" + "\t".join(str(int(c.get(k, 0))) for k in cols) +
                     f"\t{'holds' if c.get('moved', 0) == c.get('dropped', 0) == 0 else 'differs'}\n")
    os.replace(os.path.join(a.run, "f_counts.tsv.tmp"), os.path.join(a.run, "f_counts.tsv"))
    assign_table(a, a.run)
    print(open(os.path.join(a.run, "f_checks.tsv")).read())


if __name__ == "__main__":
    main()
