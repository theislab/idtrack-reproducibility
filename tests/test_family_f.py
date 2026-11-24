"""Family F (harmonization): the generator's draw, the datasets' files, the truth and its expected outputs, and every
check of the scorer, on invented datasets and an invented run written by a harmonizer that follows the documented
rules (no data root, no idtrack, no anndata)."""
import collections
import gzip
import json
import os
import random

import numpy as np
import pytest
import scipy.sparse as sp

from idtrack_reproducibility.families.f import gen_f as G
from idtrack_reproducibility.families.f import h5ad
from idtrack_reproducibility.families.f import score_f as S
from idtrack_reproducibility.families.f import truth_f as T

GA, GB, GC = "ENSSSCG00000000001.3", "ENSSSCG00000000002.1", "ENSSSCG00000000003.2"


# ---------------- the generator ----------------

def test_releases_cover_every_era_and_spread_by_releases_per_slot():
    era = lambda r: "a" if r < 60 else ("b" if r < 80 else "c")  # noqa: E731
    rels = list(range(48, 116))
    got = G.choose_releases(rels, era)
    assert len(got) == 5 and 115 not in got
    assert {era(r) for r in got} == {"a", "b", "c"}
    assert sum(era(r) == "c" for r in got) >= sum(era(r) == "a" for r in got)  # the largest era gets more slots
    assert G.choose_releases(rels, era) == got  # deterministic
    assert G.choose_releases([110, 112], era) == [110, 112]


def test_a_line_whose_broad_walk_reaches_more_ids_is_resurrected_and_retired_none_is_retired():
    assert G.class_of_a1({"expect": [GA], "broad": [GA, GB]}, "unchanged") == "resurrected"
    assert G.class_of_a1({"expect": [GA], "broad": [GA]}, "merged") == "merged"
    assert G.class_of_a1({"expect": [], "broad": None}, "retired_none") == "retired"


def _item(feature, expect, family="a", unit=None, hist="unchanged", **kw):
    line = dict(qid=f"{family}1|{feature}|90|115", set=f"{family}1_v1", expect=expect, rule="set", cmp="versioned", **kw)
    return (feature, line, family, unit or feature, hist)


def _pools():
    ids = [_item(f"ENSSSCG{i:011d}", [f"ENSSSCG{i:011d}.1"]) for i in range(100, 160)]
    names = [_item(f"N{i}", [f"ENSSSCG{i:011d}.1"], family="b", hist=None) for i in range(100, 160)]
    merged = [_item("ENSSSCG00000000900", [GA], hist="merged"), _item("ENSSSCG00000000901", [GB], hist="merged")]
    partner = _item("ENSSSCG00000000910", [GA])
    by_gene = collections.defaultdict(list)
    for it in ids + merged + [partner]:
        by_gene[it[1]["expect"][0]].append(it)
    return {("id", "unchanged"): ids, ("id", "merged"): merged, ("name", "name_unique"): names,
            ("absent", "absent"): [_item("ENSSSCG99999999999", [], family="harness")],
            ("by_gene", "id"): by_gene, ("copy_ok", None): {f"N{i}.1" for i in range(100, 160)} | {f"N{i}-1" for i in range(100, 130)},
            ("pairs", "id_and_name"): [(ids[i], names[i]) for i in range(40, 50)],
            ("pairs", "versioned_and_unversioned"): []}


def test_a_dataset_meets_its_quotas_with_distinct_features_and_copies_carry_their_names_truth():
    quotas = {("id", "unchanged"): 10, ("id", "merged"): 2, ("name", "name_unique"): 10, ("copy", "made_unique"): 4,
              ("absent", "absent"): 1, ("dup", "id_and_name"): 4, "partners": 2, "fill": ("name", "name_unique")}
    lines, drawn = G.dataset_lines("mix", 90, _pools(), [], random.Random(1), n=40, quotas=quotas)
    feats = [q["id"] for q in lines]
    assert len(lines) == 40 and len(set(feats)) == 40
    st = collections.Counter(q["stratum"] for q in lines)
    assert st["f1_merged"] == 2 and st["f1_made_unique"] == 4 and st["f1_absent"] == 1 and st["f1_duplicate"] == 4
    by = {q["id"]: q for q in lines}
    for q in lines:
        assert q["qid"] == f"f1|r090_mix|{q['id']}" and q["from_release"] == 90 and q["to_release"] == 115
        assert q["strategy"] == "best" and q["family"] == "f" and q["source"].count(":") == 2
        if q["stratum"] == "f1_made_unique":
            base = by[q["copy_of"]]
            assert q["id"] in (base["id"] + ".1", base["id"] + "-1") and q["expect"] == base["expect"]
            assert q["unit"] == base["unit"] and q["group"] == base["group"] == f"copy:{base['id']}"
    # the merged id has its partner (another id of the release expecting the same gene) in the dataset
    m = next(q for q in lines if q["id"] == "ENSSSCG00000000900")
    assert m["group"] == f"partner:{GA}" and "ENSSSCG00000000910" in by
    dups = [q for q in lines if q["stratum"] == "f1_duplicate"]
    assert len({q["group"] for q in dups}) == 2 and all(len([x for x in dups if x["group"] == q["group"]]) == 2 for q in dups)


def test_a_copy_is_drawn_only_where_the_truth_store_holds_no_string_of_its_form():
    p = _pools()
    p[("copy_ok", None)] = set()
    lines, _ = G.dataset_lines("symbols", 90, p, [], random.Random(2), n=20,
                               quotas={("name", "name_unique"): 10, ("copy", "made_unique"): 5, "fill": ("name", "name_unique")})
    assert not [q for q in lines if q["stratum"] == "f1_made_unique"]


def test_quotas_that_exceed_the_dataset_are_refused():
    with pytest.raises(SystemExit):
        G.dataset_lines("ids", 90, _pools(), [], random.Random(3), n=5, quotas={("id", "unchanged"): 10, "fill": ("id", "unchanged")})


def test_counts_are_seeded_negative_binomial_draws():
    X1, U1, m1 = G.counts(7, 300, 40)
    X2, U2, m2 = G.counts(7, 300, 40)
    assert (X1 != X2).nnz == 0 and (U1 != U2).nnz == 0 and np.allclose(m1, m2)
    assert X1.dtype == np.float32 and U1.dtype == np.int32 and X1.shape == (300, 40)
    assert (X1.data == np.round(X1.data)).all() and X1.data.min() > 0
    # the mean of a column is near the feature's mean times the mean size factor (log-normal, sigma 0.3)
    col = np.asarray(X1.mean(axis=0)).ravel()
    big = m1 > 5
    assert np.all(np.abs(col[big] / (m1[big] * np.exp(G.SIGMA_CELL ** 2 / 2)) - 1) < 0.5)
    X3, _, _ = G.counts(8, 300, 40)
    assert (X1 != X3).nnz > 0


# ---------------- the files ----------------

def _write(tmp_path, name, feats, seed=5, form="mix"):
    X, U, _ = G.counts(seed, 12, len(feats))
    obs, var, obsm, obsp, varm, uns = G.metadata("sus_scrofa", name, form, 90, seed, feats, X)
    p = str(tmp_path / f"{name}.h5ad")
    h5ad.write(p, X, obs, var, layers={"unspliced": U}, obsm=obsm, obsp=obsp, varm=varm, uns=uns)
    return p


def test_a_written_dataset_reads_back_whole(tmp_path):
    feats = ["TP53", "ENSSSCG00000000001", "X.1"]
    p = _write(tmp_path, "r090_mix", feats)
    d = h5ad.read(p)
    assert d["var"]["index"] == feats and d["obs"]["index"][0] == "c000" and d["X"].shape == (12, 3)
    assert set(d["obs"]["columns"]) == {"cell_type", "donor", "n_counts", "batch_int", "is_doublet", "stage"}
    assert d["uns"]["organism"] == "sus_scrofa" and d["uns"]["params"]["form"] == "mix" and d["uns"]["source_release"] == 90
    assert set(d["obsm"]) == {"X_pca", "X_umap"} and set(d["layers"]) == {"unspliced"} and d["varm"]["PCs"].shape == (3, 5)
    vals = h5ad.frame_values(d["obs"])
    assert vals["cell_type"][0] in ("acinar", "ductal", "alpha", "beta", "delta") and vals["batch_int"][:3] == ["0", "1", "2"]
    g = d["X"]
    assert g.indices.dtype == g.indptr.dtype  # one index dtype, as scipy requires


def test_digests_ignore_the_index_dtype_and_column_digests_are_sums_of_inputs():
    a = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0]], dtype=np.float32))
    b = sp.csr_matrix((a.data, a.indices.astype(np.int64), a.indptr.astype(np.int64)), shape=a.shape)
    assert h5ad.digest(a) == h5ad.digest(b) != h5ad.digest(a.astype(np.float64))
    assert h5ad.column_digests(a) == h5ad.column_digests(a.toarray())
    summed = sp.csr_matrix(np.array([[3], [3]], dtype=np.float32))
    assert h5ad.column_digests(summed)[0] != h5ad.column_digests(a)[0]
    assert h5ad.as_text(float("nan")) == "" and h5ad.as_text(np.int64(3)) == "3" and h5ad.as_text(True) == "True"


# ---------------- the truth ----------------

def test_a_feature_is_decided_only_where_its_truth_names_one_destination():
    assert T.decided({"expect": [GA]}) == ("gene", GA)
    assert T.decided({"expect": []}) == ("nothing", None)
    assert T.decided({"expect": [GA, GB]})[0] is None
    assert T.decided({"expect": [GA], "broad": [GA, GB]})[0] is None
    assert T.decided({"expect": [], "locus": {"seqname": "1"}})[0] is None
    assert T.decided({"expect": [GA], "open": "x"})[0] is None


def test_the_plan_follows_the_documented_rules_of_sum_and_flag():
    ins = [("ENSSSCG00000000001", GA, True, "converted"), ("GA_NAME", GA, False, "converted"),
           ("ENSSSCG00000000002", GB, True, "converted"), ("NOTHING", None, False, "not_in_graph"),
           ("ENSSSCG00000000009", GC, False, "converted_no_final_id"), ("ENSSSCG00000000008", GC, False, "converted_no_final_id")]
    s = {o["name"]: o for o in T.plan(ins, "sum")}
    assert s[GA]["members"] == ["ENSSSCG00000000001", "GA_NAME"] and s[GA]["status"] == "merged"
    assert s[GB]["status"] == "converted" and s["NOTHING"]["status"] == "unconverted_not_in_graph"
    assert s[GC]["status"] == "merged_no_final_id"
    f = T.plan(ins, "flag")
    by = {o["name"]: o for o in f}
    assert by[GA]["members"] == ["ENSSSCG00000000001"] and by["GA_NAME"]["status"] == "unconverted_merge_partner"
    # none of a group names its gene: every input keeps its name, the gene has no feature
    assert GC not in by and by["ENSSSCG00000000009"]["status"] == "unconverted_merge_partner"
    assert [o["name"] for o in f][:2] == [GA, "GA_NAME"]  # in the order of each output's first input
    # an own flag the truth cannot read leaves the group undetermined under flag, never under sum
    und = T.plan([("A", GA, None, "converted"), ("B", GA, False, "converted")], "flag")
    assert {o["status"] for o in und} == {"undetermined"}
    assert T.plan([("A", GA, None, "converted"), ("B", GA, False, "converted")], "sum")[0]["status"] == "merged"


def test_an_unconverted_name_equal_to_a_converted_feature_is_renamed():
    out = T.plan([("X", GA, True, "converted"), (GA, None, False, "no_conversion")], "sum")
    assert [o["name"] for o in out] == [GA, f"unconverted:{GA}"]
    renamed, names = T.merged_names({"d1": [(GB, "converted"), ("Z", "unconverted_not_in_graph")],
                                     "d2": [("Z", "converted"), (GB, "converted")]}, "union")
    assert names == {GB, "Z", "unconverted:Z"} and renamed["d1"]["unconverted:Z"] == "Z"
    _, inter = T.merged_names({"d1": [(GB, "converted"), ("Z", "x")], "d2": [(GB, "converted")]}, "intersect")
    assert inter == {GB}


def test_own_is_read_for_ids_and_for_names_only_in_the_gene_configuration():
    assert T.own_by_truth("ENSSSCG00000000001", GA, "symbol") is True
    assert T.own_by_truth("ENSSSCG00000000001.2", GA, "gene") is True
    assert T.own_by_truth("ENSSSCG00000000005", GA, "gene") is False
    assert T.own_by_truth("TP53", GA, "gene") is False and T.own_by_truth("TP53", GA, "symbol") is None


def test_uns_keeps_the_entries_equal_in_every_dataset():
    u = T.uns_kept([{"o": "pig", "d": "a", "p": {"s": 1, "f": "x"}}, {"o": "pig", "d": "b", "p": {"s": 2, "f": "x"}}])
    assert u == {"o": "pig", "p": {"f": "x"}}
    assert T.uns_kept([{"d": "a"}, {"d": "b"}]) == {}


def test_the_symbol_judge_opens_exact_multi_gene_lines_within_their_genes():
    lines = [dict(qid="q1", expect=[GA, GB], rule="set_exact"), dict(qid="q2", expect=[GA, GB], rule="set"),
             dict(qid="q3", expect=[GA], rule="set_exact")]
    out = T.symbol_judge_lines(lines)
    assert [q["qid"] for q in out] == ["q1"] and out[0]["open_within"] == [GA, GB] and "symbol_tie" in out[0]["revised_by"]


def test_the_symbol_configuration_reads_the_symbol_judges_sets(monkeypatch):
    from idtrack_reproducibility.families import lines as L
    mans = {"f1_v1": {}, "f1_symbol_v1": {"revises": ["f1_v1"], "judge": "symbol"}, "f1_r1": {"revises": ["f1_v1"]}}
    monkeypatch.setattr(L, "set_names", lambda family, organism: list(mans))
    monkeypatch.setattr(L, "manifest", lambda family, organism, name: mans[name])
    assert S.sets_for("sus_scrofa", "gene") == ["f1_v1", "f1_r1"]
    assert S.sets_for("sus_scrofa", "symbol") == ["f1_v1", "f1_symbol_v1", "f1_r1"]


# ---------------- the checks, on an invented run ----------------

FEATS = {"r090_ids": ["ENSSSCG00000000001", "ENSSSCG00000000011", "ENSSSCG00000000002", "ENSSSCG99999999999"],
         "r090_mix": ["ENSSSCG00000000001", "GA_NAME", "GA_NAME.1", "NOTHING", "GC_NAME"]}
TRUTH = {"ENSSSCG00000000001": [GA], "ENSSSCG00000000011": [GA], "ENSSSCG00000000002": [GB], "ENSSSCG99999999999": [],
         "GA_NAME": [GA], "GA_NAME.1": [GA], "NOTHING": [], "GC_NAME": [GC]}


def _write_tsv_gz(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")


def _json(path, d):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(d, open(path, "w"))


def _assignment(f):
    """The invented harmonizer's assignment: the truth's gene (GC_NAME: converted to GB, a wrong assignment)."""
    if f == "GC_NAME":
        return GB
    e = TRUTH[f]
    return e[0] if len(e) == 1 else None


def _harmonize(inp, d, merge, cfg):
    """A dataset harmonized by the documented rules from the invented assignment: (rows, meta)."""
    ins = [(f, _assignment(f), T.own_by_truth(f, _assignment(f), "gene") if _assignment(f) else False,
            "converted" if _assignment(f) else "not_in_graph") for f in inp.features]
    plan = T.plan(ins, merge)
    rows = []
    for o in plan:
        ve = inp.var_expected(o["members"])
        rows.append([o["name"], ";".join(o["members"]), o["status"], o["gene"] or "", "", "", inp.sum_digest(o["members"]),
                     inp.sum_digest(o["members"], True)] + [ve[c] for c in S.INPUT_VAR])
    summed = any(len(o["members"]) > 1 for o in plan)
    varm = {} if summed else {k: h5ad.digest(np.asarray(v)[[inp.pos[o["members"][0]] for o in plan]]) for k, v in inp.varm.items()}
    meta = dict(obs_names=list(inp.obs["index"]), obs_values=inp.obs_values, obs_columns=list(inp.obs_values),
                obs_dtypes={"x": "y"}, obsm=inp.obsm, obsp=inp.obsp, uns=inp.uns, layers=["unspliced"], varm=varm)
    return rows, meta, plan


def _merged(inputs, outs, feats, handles):
    """The merged object's features, blocks, obs and meta from the datasets' outputs, as the docs say."""
    per = {d: [(r[0], r[2]) for r in outs[d]] for d in handles}
    renamed, names = T.merged_names(per, feats)
    names = sorted(names)
    frows = [[n, "1" if all(n in renamed[d] for d in handles) else "0"] for n in names]
    blocks = []
    for d in handles:
        by = {r[0]: r for r in outs[d]}
        for n in names:
            if n in renamed[d]:
                r = by[renamed[d][n]]
                blocks.append([n, d, r[1], r[2], r[6], r[7]] + r[8:])
    all_cols = []
    for d in handles:
        all_cols += [c for c in inputs[d].obs_values if c not in all_cols]
    obs = [[f"{d}_{c}", d] + [inputs[d].obs_values[col][i] if col in inputs[d].obs_values else "" for col in all_cols]
           for d in handles for i, c in enumerate(inputs[d].obs["index"])]
    meta = dict(nonzero_blocks_of_absent_features=0, obs_dtypes={**{c: "x" for c in all_cols}, "handle_anndata": "category"},
                obsm_blocks={d: inputs[d].obsm for d in handles}, obsp_blocks={d: inputs[d].obsp for d in handles},
                obsp_off_block_nnz=0, layers=["unspliced"], uns=T.uns_kept([inputs[d].uns for d in handles]))
    return frows, blocks, (["obs_name", "handle_anndata"] + all_cols, obs), meta


@pytest.fixture
def run(tmp_path):
    """Two invented datasets, their lines, and a run of the invented harmonizer in run_f's format (configuration gene)."""
    handles = sorted(FEATS)
    inputs, lines, man = {}, {}, {"datasets": {}}
    for k, d in enumerate(handles):
        p = _write(tmp_path, d, FEATS[d], seed=11 + k, form=d.split("_")[1])
        inputs[d] = S.Inputs(p)
        man["datasets"][d] = {"file": p}
        for f in FEATS[d]:
            lines[f"f1|{d}|{f}"] = dict(qid=f"f1|{d}|{f}", expect=TRUTH[f], stratum="f1_x" if f != "GA_NAME.1" else "f1_made_unique",
                                        rule="set")
    base = tmp_path / "run" / "gene"
    outs = {m: {} for m in S.MERGES}
    with gzip.open(os.path.join(os.makedirs(base, exist_ok=True) or base, "answers.jsonl.gz"), "wt") as fo:
        for d in handles:
            for f in FEATS[d]:
                g = _assignment(f)
                fo.write(json.dumps(dict(qid=f"f1|{d}|{f}", got=[g] if g else [], status="converted" if g else "not_in_graph",
                                         own=bool(g and T.own_by_truth(f, g, "gene")))) + "\n")
    for merge in S.MERGES:
        for d in handles:
            rows, meta, _ = _harmonize(inputs[d], d, merge, "gene")
            outs[merge][d] = rows
            _write_tsv_gz(str(base / "outputs" / merge / f"{d}.tsv.gz"),
                          ["name", "query_ids", "status", "gene", "final_id", "reason", "x_digest", "u_digest"] +
                          [f"var:{c}" for c in S.INPUT_VAR], rows)
            _json(str(base / "outputs" / merge / f"{d}.meta.json"), meta)
        for feats in ("union", "intersect"):
            for where, stem in (("unify", f"{feats}_{merge}"), ("shuffle", f"union_{merge}")):
                if where == "shuffle" and feats != "union":
                    continue
                frows, blocks, (oh, orows), meta = _merged(inputs, outs[merge], feats, handles)
                s = str(base / where / stem)
                _write_tsv_gz(s + ".features.tsv.gz", ["name", "intersection"], frows)
                _write_tsv_gz(s + ".blocks.tsv.gz", ["name", "dataset", "query_ids", "status", "x_digest", "u_digest"] +
                              [f"var:{c}" for c in S.INPUT_VAR], blocks)
                _write_tsv_gz(s + ".obs.tsv.gz", oh, orows)
                _json(s + ".meta.json", meta)
    digest = "d" * 64
    hdr = lambda d: dict(format=4, rules_version=4, organism="sus_scrofa", final_database="ensembl_gene",  # noqa: E731
                         target_ensembl_release=115, source_release=90, source_database=None, graph_snapshot=115,
                         graph_content_digest=digest, feature_names_sha256=d)
    _json(str(base / "cache.json"), dict(headers={d: hdr(d) for d in handles}, graph_content_digest=digest,
                                         again=dict(seconds=1, construct_seconds=5, equal=True, files_rewritten=[]),
                                         changed_digest=dict(dataset=handles[0], rewritten_with=digest, equal=True),
                                         without_api=dict(graph_loaded=0, equal=True),
                                         shuffle=dict(recomputed=handles, answers_differ={d: [] for d in handles})))

    class A:
        pass
    a = A()
    a.run = str(tmp_path / "run")
    kinds = {"cell_type": "category", "donor": "str", "n_counts": "float64", "batch_int": "int64", "is_doublet": "bool",
             "stage": "str", "plate": "int64"}
    a.input_check = {d: {"obs_dtypes": {c: kinds[c] for c in inputs[d].obs_values}} for d in handles}
    for d in handles:
        for m in S.MERGES:
            p = base / "outputs" / m / f"{d}.meta.json"
            meta = json.load(open(p))
            meta["obs_dtypes"] = {c: a.input_check[d]["obs_dtypes"][c] for c in inputs[d].obs_values}
            _json(str(p), meta)
    for where in ("unify", "shuffle"):
        for p in (base / where).glob("*.meta.json"):
            meta = json.load(open(p))
            meta["obs_dtypes"] = {"cell_type": "category", "donor": "category", "n_counts": "float64", "batch_int": "int64",
                                  "is_doublet": "bool", "stage": "category", "plate": "Int64", "handle_anndata": "category"}
            _json(str(p), meta)
    return a, inputs, lines, man, base


def _checks(run_):
    a, inputs, lines, man, _base = run_
    ck, counts = S.Checks(), {}
    S.score_config(a, "gene", ck, inputs, lines, man, counts)
    bad = {k[0]: v for k, v in ck.rows.items() if v[2]}
    return ck, bad, counts


def test_a_harmonizer_that_follows_its_docs_passes_every_check_but_the_truths(run):
    ck, bad, counts = _checks(run)
    # the invented harmonizer converts GC_NAME to GB: the truth's checks see it, the mechanical ones do not
    assert "truth_feature" in bad and set(bad) <= {"truth_feature", "unify_truth_union", "unify_truth_intersect"}
    names = {k[0] for k in ck.rows}
    for c in ("plan_feature", "plan_values_x", "plan_values_layer", "plan_obs", "plan_obsm", "plan_obsp", "plan_uns",
              "plan_varm_left_out_when_summed", "unify_features", "unify_block", "unify_obs_values", "unify_uns",
              "shuffle_blocks", "shuffle_answers", "cache_header_digest", "cache_without_api_no_graph_load"):
        assert c in names, c
    # count conservation: GC_NAME's counts moved; the others kept; nothing dropped
    c = counts[("gene", "sum", "f1_x")]
    assert c["moved"] == 1 and c["dropped"] == 0 and c["kept"] == c["features"] - 1


def _edit_rows(path, fn):
    rows = S.read_tsv_gz(path)
    head = list(rows[0])
    rows = fn(rows)
    _write_tsv_gz(path, head, [[r[h] for h in head] for r in rows])


@pytest.mark.parametrize("what,check", [
    ("x", "plan_values_x"), ("members", "plan_feature"), ("obs", "plan_obs"), ("uns", "plan_uns"), ("obsm", "plan_obsm"),
    ("lost", "plan_every_input_once")])
def test_each_dataset_check_sees_its_violation(run, what, check):
    a, inputs, lines, man, base = run
    d = "r090_ids"
    tsv, meta_p = str(base / "outputs" / "sum" / f"{d}.tsv.gz"), str(base / "outputs" / "sum" / f"{d}.meta.json")
    meta = json.load(open(meta_p))
    if what == "x":
        _edit_rows(tsv, lambda rows: [dict(r, x_digest="0" * 16) if i == 0 else r for i, r in enumerate(rows)])
    elif what == "members":
        _edit_rows(tsv, lambda rows: [dict(r, query_ids=r["query_ids"].split(";")[0]) for r in rows])
    elif what == "lost":
        _edit_rows(tsv, lambda rows: rows[:-1])
    elif what == "obs":  # issue #85: obs lost between input and output
        meta["obs_values"] = {}
        meta["obs_columns"] = []
    elif what == "uns":
        meta["uns"] = {}
    elif what == "obsm":
        meta["obsm"] = {}
    _json(meta_p, meta)
    _, bad, _ = _checks(run)
    assert check in bad


@pytest.mark.parametrize("what,check", [
    ("block", "unify_block"), ("absent", "unify_absent_blocks_zero"), ("obs", "unify_obs_values"), ("uns", "unify_uns"),
    ("feature", "unify_features"), ("obsp", "unify_obsp_blocks")])
def test_each_merged_object_check_sees_its_violation(run, what, check):
    a, inputs, lines, man, base = run
    stem = str(base / "unify" / "union_sum")
    meta = json.load(open(stem + ".meta.json"))
    if what == "block":
        _edit_rows(stem + ".blocks.tsv.gz", lambda rows: [dict(r, x_digest="f" * 16) if i == 0 else r for i, r in enumerate(rows)])
    elif what == "absent":
        meta["nonzero_blocks_of_absent_features"] = 3
    elif what == "obs":
        _edit_rows(stem + ".obs.tsv.gz", lambda rows: [dict(r, donor="") for r in rows])
    elif what == "uns":
        meta["uns"] = {}
    elif what == "feature":
        _edit_rows(stem + ".features.tsv.gz", lambda rows: rows[1:])
    elif what == "obsp":
        meta["obsp_off_block_nnz"] = 1
    _json(stem + ".meta.json", meta)
    _, bad, _ = _checks(run)
    assert check in bad


def test_the_shuffle_relation_sees_a_feature_whose_values_or_answer_change(run):
    a, inputs, lines, man, base = run
    stem = str(base / "shuffle" / "union_flag")
    _edit_rows(stem + ".blocks.tsv.gz", lambda rows: [dict(r, x_digest="e" * 16) if i == 1 else r for i, r in enumerate(rows)])
    c = json.load(open(base / "cache.json"))
    c["shuffle"]["answers_differ"]["r090_mix"] = ["GA_NAME"]
    _json(str(base / "cache.json"), c)
    _, bad, _ = _checks(run)
    assert "shuffle_blocks" in bad and "shuffle_answers" in bad
    # a permutation changes the order of a merged feature's inputs, never the relation
    _edit_rows(stem + ".blocks.tsv.gz", lambda rows: [dict(r, x_digest=r["x_digest"]) for r in rows])


def test_the_cache_checks_see_a_package_version_a_stale_digest_and_a_graph_loaded(run):
    a, inputs, lines, man, base = run
    c = json.load(open(base / "cache.json"))
    c["headers"]["r090_ids"]["idtrack_version"] = "0.1.1"
    c["headers"]["r090_mix"]["graph_content_digest"] = "x"
    c["without_api"]["graph_loaded"] = 1
    c["again"]["seconds"] = 9
    _json(str(base / "cache.json"), c)
    _, bad, _ = _checks(run)
    assert {"cache_header_no_package_version", "cache_header_digest", "cache_without_api_no_graph_load",
            "cache_again_faster"} <= set(bad)


def test_merged_obs_dtypes_follow_the_documented_rules():
    assert S.kinds_ok("str", "category", False) and not S.kinds_ok("str", "str", False)
    assert S.kinds_ok("int64", "Int64", True) and S.kinds_ok("int64", "int64", False) and not S.kinds_ok("int64", "float64", True)
    assert S.kinds_ok("bool", "boolean", True) and S.kinds_ok("category", "category", False)


def test_the_plugin_registers_fs_kinds_and_its_causes():
    from idtrack_reproducibility.families import causes as K
    from idtrack_reproducibility.families.f import causes_f  # noqa: F401  registers on import
    assert "copy_of_another_feature" in K.CAUSE_CLASS and "batch_differs_from_single" in K.CAUSE_CLASS
    names = [t.__name__ for t in K.TESTS]
    assert names.index("t_copy_of_another_feature") < names.index("t_rule_decided_after_tree")

    class Ctx:
        r = {"stratum": "f1_made_unique", "reason": "x", "query": json.dumps({"id": "X.1"})}
        same_answer = False
    assert causes_f.t_copy_of_another_feature(Ctx()) and causes_f.t_batch_differs_from_single(Ctx())
    Ctx.same_answer = True
    assert causes_f.t_copy_of_another_feature(Ctx()) is None


def test_the_best_view_of_a_shared_name_is_a_cause_only_within_the_carriers_and_the_wider_set():
    from idtrack_reproducibility.families import causes as K
    from idtrack_reproducibility.families.f import causes_f

    def ctx(got, expected, broad=""):
        r = dict(qid="q", kind="f1_symbols", stratum="f1_name_shared", variant="best,from", unit="X", reason="",
                 query=json.dumps({"id": "X", "from_release": 90, "to_release": 115}), got=got, expected=expected, broad=broad)
        return K.Ctx(r, None, None)
    assert causes_f.t_best_keeps_part_of_shared_name(ctx(GA, f"{GA},{GB}"))
    assert causes_f.t_best_keeps_part_of_shared_name(ctx(GC, f"{GA},{GB}", f"{GA},{GB},{GC}"))
    assert causes_f.t_best_keeps_part_of_shared_name(ctx(GC, f"{GA},{GB}")) is None  # another gene: not this cause
    assert causes_f.t_best_keeps_part_of_shared_name(ctx(GA, GA)) is None


def test_a_string_of_the_store_takes_every_form_that_folds_like_it():
    idx = G.fold_index({"Ighv5-9-1", "Ighv5-9.1", "TP53.1"})
    assert idx["ighv5-9-1"] == {"Ighv5-9-1", "Ighv5-9.1"} and idx["tp53-1"] == {"TP53.1"}
