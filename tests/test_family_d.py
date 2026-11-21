"""Family D (real-world inputs): the reduction of the matrices' features, the place truth and its chain mapping, the
release of a list, the mixed lists, D's revisions and its limitations of open lines, on invented models and lists."""
import datetime
import importlib
import itertools
import types

import numpy as np
import pytest

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families.d import gen_d2 as G2
from idtrack_reproducibility.families.d import gen_d4 as G4
from idtrack_reproducibility.families.d import limits_d as LD
from idtrack_reproducibility.families.d import revise_d1 as RD1
from idtrack_reproducibility.families.d import revise_d4 as RD4


@pytest.fixture
def sources(tmp_path, monkeypatch):
    """A configuration naming the atlas project's inputs (read when D1's and D3's generators are imported)."""
    cfg = tmp_path / "idtrack-repro.toml"
    cfg.write_text('[data]\nroot = "data"\n[sources]\nrosetta = "atlas"\nrosetta_m02 = "atlas_m02"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(cfg))
    CFG.reset()
    yield lambda name: importlib.import_module(f"idtrack_reproducibility.families.d.{name}")
    CFG.reset()


def models(G1, genes):
    """Gene models of one release: gene -> merged exons on chromosome 1 (strand `-` for names ending in `_minus`)."""
    M = object.__new__(G1.Models)
    M.genes = {g: dict(seq="1", strand="-" if g.endswith("_minus") else "+", start=ex[0][0], end=ex[-1][1], exons=ex,
                       bases=sum(e - s + 1 for s, e in ex), multi=False, version="1", biotype="protein_coding", name=g)
               for g, ex in genes.items()}
    M.names, M.lower, M._ix, M._sx = {}, {}, None, None
    return M


def place(G1, target, exons):
    """The place class of a model of release 100 against the target's models (one assembly: no chain)."""
    loc = object.__new__(G1.Locus)
    loc.asm, loc.T, loc.target = {100: "GRCh38", 115: "GRCh38"}, 115, models(G1, target)
    return loc._of(dict(seq="1", strand="+", exons=exons, bases=sum(e - s + 1 for s, e in exons), multi=False), 100)


def test_features_are_reduced_as_the_atlas_project_reduces_them(sources):
    G1 = sources("gen_d1")
    g = dict(type="ensembl_gene_versioned", collapse_suffix=r"_PAR_Y$", other_features=["EGFP"])
    assert G1.reduce("EGFP", g) is None
    assert G1.reduce("ENSG00000141510.17", g) == ("ENSG00000141510", "id", "ENSG00000141510.17")
    assert G1.reduce("ENSG00000002586.20_PAR_Y", g) == ("ENSG00000002586", "id", "ENSG00000002586.20")
    assert G1.reduce("TP53", dict(type="symbol_else_id")) == ("TP53", "symbol", None)
    assert G1.reduce("ENSG00000141510", dict(type="symbol_else_id")) == ("ENSG00000141510", "id", None)
    assert G1.reduce("ENSG00000141510", dict(type="symbol")) == ("ENSG00000141510", "symbol", None)


def test_the_place_of_a_model_names_one_gene_a_dominant_gene_or_none(sources):
    G1 = sources("gen_d1")
    model = [(1000, 1999)]  # 1,000 bases
    x = place(G1, {"A": [(1000, 1999)], "B_minus": [(1000, 1999)]}, model)
    assert (x["cls"], x["expect"]) == ("one", ["A"])  # a gene on the other strand shares nothing
    x = place(G1, {"A": [(1000, 1399)]}, model)
    assert (x["cls"], x["expect"], x["open_within"]) == ("below_share", [], ["A"])
    assert place(G1, {"A_minus": [(1000, 1999)]}, model)["cls"] == "none"
    x = place(G1, {"A": [(1000, 1899), (3000, 9999)], "S": [(1950, 2100)]}, model)  # 900 and 50 shared bases
    assert (x["cls"], x["expect"], x["dominant_by"]) == ("dominant", ["A"], "bases")
    x = place(G1, {"A": [(1000, 1899)], "S": [(1950, 2100)]}, model)  # A also holds half of its own bases
    assert (x["cls"], x["expect"], x["dominant_by"]) == ("dominant", ["A"], "bases,mutual")
    # a readthrough over the gene: both share the whole model, only the gene holds half of its own bases
    x = place(G1, {"GENE": [(1000, 1999)], "RT": [(1000, 1999), (5000, 8999)]}, model)
    assert (x["cls"], x["expect"], x["dominant_by"]) == ("dominant", ["GENE"], "mutual")
    # a gene grown over the model and a small new gene inside it: the small one is no answer (revision d1_r2)
    x = place(G1, {"BIG": [(500, 4999)], "SMALL": [(1000, 1599)]}, model)
    assert (x["cls"], x["expect"]) == ("several", []) and set(x["open_within"]) == {"BIG", "SMALL"}


def test_models_of_another_assembly_are_mapped_exon_by_exon_through_the_chain(sources):
    G1 = sources("gen_d1")

    class Chain:
        def __init__(self, blocks):
            blocks = sorted(blocks)
            self.index = {"1": (blocks, list(itertools.accumulate((b[1] for b in blocks), max)))}

    fwd = Chain([(0, 1000, "1", "+", 100000, 5000)])  # 0-based target block -> query start
    assert G1.map_exons(fwd, "1", "+", [(101, 200), (301, 400)], 200) == (("1", "+", [(5101, 5200), (5301, 5400)], 1.0), "")
    rev = Chain([(0, 1000, "2", "-", 10000, 0)])
    assert G1.map_exons(rev, "1", "+", [(101, 200)], 100) == (("2", "-", [(9801, 9900)], 1.0), "")
    short = Chain([(0, 250, "1", "+", 100000, 0)])
    res, why = G1.map_exons(short, "1", "+", [(101, 200), (301, 1300)], 1100)
    assert res is None and "aligned" in why  # under half of the model's bases
    assert G1.map_exons(fwd, "9", "+", [(1, 10)], 10)[0] is None


def test_a_symbol_of_several_genes_a_made_unique_name_and_an_unknown_name(sources):
    G1 = sources("gen_d1")
    per = {"G1": dict(cls="one", expect=["T1"], open_within=None, hits={"T1": 900}),
           "G2": dict(cls="dominant", expect=["T2"], open_within=None, hits={"T2": 600, "T3": 10}),
           "G3": dict(cls="several", expect=[], open_within=["T4", "T5"], hits={"T4": 300, "T5": 300})}
    names = {"MALAT1": ["G1", "G2"], "MIX": ["G1", "G3"], "SOLO": ["G1"]}
    loc = types.SimpleNamespace(source=lambda R: types.SimpleNamespace(names=names, lower={k.lower(): v for k, v in names.items()}),
                                of=lambda g, R: per[g])
    cls, E, within, genes, match, _ = G1.symbol_truth(loc, "MALAT1", 77)
    assert (cls, E, within, genes, match) == ("name_several", ["T1", "T2"], None, ["G1", "G2"], "exact")
    assert G1.symbol_truth(loc, "Solo", 77)[:5] == ("one", ["T1"], None, ["G1"], "case")
    cls, E, within, genes, match, _ = G1.symbol_truth(loc, "MIX", 77)
    assert (cls, E, within) == ("name_mixed", [], ["T1", "T4", "T5"])
    cls, E, within, genes, match, _ = G1.symbol_truth(loc, "SOLO.1", 77)
    assert (cls, E, within, genes) == ("made_unique", [], ["T1"], ["G1"])
    assert G1.symbol_truth(loc, "NOPE", 77)[:3] == ("name_not_at_source", [], None)


def presence(rels, kinds):
    """D2's presence matrices of invented releases: kind -> {release: items}."""
    P = object.__new__(G2.Presence)
    P.rels, P.np, P.index, P.matrix, P.no_xref = rels, np, {}, {}, []
    for kind, sets in kinds.items():
        keys = sorted(set().union(*sets.values()))
        ix = {x: i for i, x in enumerate(keys)}
        m = np.zeros((len(keys) + 1, len(rels)), dtype=np.int32)
        for j, r in enumerate(rels):
            m[[ix[x] for x in sets[r]], j] = 1
        P.index[kind], P.matrix[kind] = ix, m
    return P


def test_the_release_of_a_list_is_its_tie_and_its_band():
    P = presence([1, 2, 3], {"ids": {1: {"A", "B", "C"}, 2: {"A", "B", "C", "D"}, 3: {"A", "D", "E"}},
                             "symbols": {1: {"X"}, 2: {"Y"}, 3: {"Y"}}})
    assert P.truth(["A", "B", "C"], "ids")[0] == [1, 2]  # a tie is the whole set
    assert P.truth(["A", "B", "C", "D"], "ids")[0] == [2]
    E, B, shares = P.truth(["A", "B", "NOT_A_GENE_1"], "ids")
    assert E == [1, 2] and shares[1] == pytest.approx(2 / 3) and shares[3] == pytest.approx(1 / 3)
    assert B == [1, 2, 3]  # 3 lacks one item of the best: no significant difference
    assert P.truth(["A", "B", "B", "C"], "ids")[2][1] == pytest.approx(1.0)  # a distinct item counts once
    assert P.truth(["A", "X"], "mixed")[0] == [1]  # an item counts as an id or as a name of the release
    q = G2.line("sus_scrofa", "d2m|x|0|all", "d2_matrix_ids_stated", ["A", "B", "C", "D"], P, "ids", unit="x", stated=3)
    assert q["expect"] == ["2"] and q["open"] and q["stated_in_expected"] is False and q["call"] == "infer_source"
    q = G2.line("sus_scrofa", "d2m|x|0|all", "d2_matrix_ids_stated", ["A", "B", "C"], P, "ids", unit="x", stated=2)
    assert q["expect"] == ["1", "2"] and "open" not in q and q["rule"] == "within_band" and q["band"][:2] == ["1", "2"]


def test_the_band_is_the_documented_sign_test():
    """Rule `inference_per_release`: a release is out of the band when the items found at the best and not there significantly outnumber
    those found there and not at the best (one-sided sign test at 0.05). Fails if a lone missing item drops a release, or
    if many missing items keep it."""
    items = [f"I{k}" for k in range(20)]
    P = presence([1, 2, 3], {"ids": {1: set(items), 2: set(items[:-4]), 3: set(items[:-5])}})
    T, B, _ = P.truth(items, "ids")
    assert T == [1] and B == [1, 2]  # 4 items against none: p = 1/16, kept; 5: p = 1/32, out
    assert not G2.significantly_fewer(4, 0) and G2.significantly_fewer(5, 0) and not G2.significantly_fewer(3, 3)
    assert G2.significantly_fewer(2000, 1800) and not G2.significantly_fewer(1810, 1800)  # the normal approximation


def test_within_band_accepts_the_tie_and_the_band_and_nothing_else():
    """Rule `within_band`: an answer holds every release of the tie and none outside the band. Fails if the tie alone or the
    band is wrong, or if a release outside the band or a missing tied release passes."""
    from idtrack_reproducibility.families import scorer as S
    line = dict(qid="d2s|x", rule="within_band", expect=["10", "11"], band=["9", "10", "11", "12"], strategy="best", cmp="versioned")

    def ans(got):
        return dict(outcome="answer", got=got)
    assert S.score_line(line, "best", ans(["10", "11"]), None)[0] == "right"
    assert S.score_line(line, "best", ans(["9", "10", "11", "12"]), None)[0] == "right"
    assert S.score_line(line, "best", ans(["10"]), None) == ("wrong", "missing from the tie: 11")
    assert S.score_line(line, "best", ans(["10", "11", "13"]), None) == ("wrong", "outside the band: 13")


def test_a_blank_is_an_error_of_its_own_item(sources):
    """Revision d4_r3: the items of a list holding a blank expect their own truth again; the blank keeps the error.
    Fails if an item keeps the whole list's error, or the blank loses it."""
    R4 = sources("revise_d4")
    first = {"d4|1|0": dict(qid="d4|1|0", id="TP53", expect=["ENSG00000141510"]), "d4|1|1": dict(qid="d4|1|1", id=" ", expect=[], expect_error=R4.ERR),
             "d4|1|2": dict(qid="d4|1|2", id="nan", expect=[], expect_error="MissingValueError")}
    latest = [dict(first["d4|1|0"], expect=[], expect_alone=["ENSG00000141510"], expect_error=R4.ERR, revised_by=["blank_refuses_batch"]),
              dict(first["d4|1|1"], revised_by=["blank_refuses_batch"]),
              dict(first["d4|1|2"], expect_alone=[], expect_error=R4.ERR, revised_by=["blank_refuses_batch"])]
    out = {q["qid"]: q for q in R4.revise_blank_items("homo_sapiens", latest, first)}
    assert out["d4|1|0"]["expect"] == ["ENSG00000141510"] and "expect_error" not in out["d4|1|0"]
    assert out["d4|1|1"]["expect_error"] == R4.ERR
    assert out["d4|1|2"]["expect_error"] == "MissingValueError" and out["d4|1|2"]["open"] == R4.OPEN_52


def test_a_marker_list_is_read_at_the_latest_release_published_before_its_date(sources):
    G3 = sources("gen_d3")
    dates = {98: datetime.date(2019, 9, 1), 99: datetime.date(2020, 1, 1), 100: datetime.date(2020, 4, 1)}
    assert G3.release_at(datetime.date(2020, 3, 27), dates, [98, 99, 100]) == 99
    assert G3.release_at(datetime.date(2020, 4, 1), dates, [98, 99, 100]) == 100
    assert G3.release_at(datetime.date(2010, 1, 1), dates, [98, 99, 100]) == 98  # before the window: its first release
    N = object.__new__(G3.Names)
    N.exact, N.lower = {99: {"Cd4": {"G1"}}}, {99: {"cd4": {"G1"}}}
    assert N.genes("Cd4", 99) == (["G1"], "exact") and N.genes("CD4", 99) == (["G1"], "case")
    assert N.genes("CD45", 99) == ([], "none")


def test_mixed_list_items_keep_their_harness_truth_and_their_place_in_the_list():
    q = dict(qid="h|1", kind="name_unique_to_gene", id="TP53", expect=["ENSG00000141510"], rule="set", cmp="base", broad=None,
             support={"ENSG00000141510.17": 3})
    it = G4.item(q, "symbol")
    assert it["expect"] == q["expect"] and it["support"] == q["support"] and it["src_qid"] == "h|1" and "broad" not in it
    P = {"symbol": [q], "nothing": [dict(qid="h|2", kind="nonexistent", id="NOT_A_GENE", expect=[], rule="set", cmp="versioned")]}
    lists = list(G4.make_lists("sus_scrofa", P, 40))
    assert lists == list(G4.make_lists("sus_scrofa", P, 40))  # seeded per list
    items = [x for _, xs in lists for x in xs]
    missing = [x for x in items if x["category"] == "missing"]
    assert missing and all(x["id"] in G4.MISSING and x["expect_error"] == G4.EXPECT_MISSING for x in missing)
    assert any(x.get("duplicate") for x in items)
    k, xs = lists[0]
    lines = list(G4.lines_of("sus_scrofa", k, xs))
    assert [x["pos"] for x in lines] == list(range(len(xs))) and all(x["batch"] == [y["id"] for y in xs] for x in lines)
    assert lines[0]["qid"] == f"d4|{k}|0" and lines[0]["call"] == "convert_many" and lines[0]["kind"] == f"d4_{xs[0]['category']}"


def test_revisions_of_d1_and_d4():
    assert RD1.name_several_set(dict(truth_class="name_several", rule="set_exact"))["rule"] == "set"
    assert RD1.name_several_set(dict(truth_class="one", rule="set")) is None
    lines = [dict(qid="1", id="nan", expect_error="E"), dict(qid="2", id=" NA ", expect_error="E"),
             dict(qid="3", id="", expect_error="E"), dict(qid="4", id="TP53", expect=["G"]),
             dict(qid="5", id="null", expect_error="E", open="already open")]
    out = RD4.revise_52(lines)
    assert [x["qid"] for x in out] == ["1", "2"]  # a blank refuses the whole batch: another revision
    assert all(x["open"] == RD4.OPEN_52 and x["expect_error"] == "E" and len(x["revised_by"]) == 1 for x in out)


def test_the_open_reasons_of_the_truth_map_to_their_limitations(sources):
    G1 = sources("gen_d1")
    want = {"several": "L-place-several", "below_share": "L-place-shrunk", "none": "L-place-empty",
            "unmapped": "L-place-unmapped", "patch_at_source": "L-patch-source", "not_at_source": "L-id-not-at-source",
            "name_not_at_source": "L-name-not-at-source", "made_unique": "L-made-unique", "name_mixed": "L-name-several-place"}
    for cls, lim in want.items():
        assert LD.classify(G1.OPEN_TEXT[cls])[0] == lim, cls
    assert LD.classify(RD4.OPEN_52)[0] == "L-missing-values"
    assert LD.classify("the stated release 92 fits the list worse (0.9) than 89..89 (1.0) by Ensembl's GTFs")[0] == "L-stated-release"
    assert LD.classify("a reason no row names")[0] == "L-undetermined-other"
