"""The judge's rules on invented lines and answers: the harness's function and the common scorer's own rules."""
import math

from idtrack_reproducibility.families import scorer as SC
from idtrack_reproducibility.harness import score as S
from idtrack_reproducibility.harness import harness_lib as H


class Places:
    """Gene places of an invented assembly: {release: {gene: row}}."""

    def __init__(self, genes, assembly="asm1"):
        self.genes, self.asm = genes, assembly

    def get(self, rel, gid):
        return self.genes.get(rel, {}).get(gid)

    def core(self, rel, gid):
        return None

    def assembly(self, rel):
        return self.asm

    def projected(self, rel, g):
        return g


def line(**kw):
    rec = dict(qid="q", kind="k", id="X", from_release=1, to_release=3, final_database=None, strategy="best",
               rule="set", cmp="versioned", expect=[], broad=None, locus=None)
    rec.update(kw)
    return rec


def answer(got, **kw):
    return dict(got=list(got), **kw)


NO = Places({})


def test_set_best_accepts_a_part_of_the_expected_set():
    rec = line(expect=["G1.1", "G2.1"])
    assert S.score_one(rec, answer(["G1.1"]), NO)[0] == "right"
    assert S.score_one(rec, answer(["G1.1", "G3.1"]), NO)[0] == "wrong"
    assert S.score_one(rec, answer([]), NO) == ("wrong", "no answer")


def test_all_and_set_exact_need_the_whole_set():
    rec = line(expect=["G1.1", "G2.1"], strategy="all")
    assert S.score_one(rec, answer(["G1.1"]), NO)[0] == "wrong"
    assert S.score_one(rec, answer(["G2.1", "G1.1"]), NO)[0] == "right"
    rec = line(expect=["G1.1", "G2.1"], rule="set_exact")
    assert S.score_one(rec, answer(["G1.1"]), NO)[0] == "wrong"


def test_nothing_expected_and_nothing_answered_is_right():
    assert S.score_one(line(expect=[]), answer([]), NO)[0] == "right"


def test_broad_makes_an_answer_undetermined_not_right():
    rec = line(expect=["G1.1"], broad=["G1.1", "G2.1"])
    assert S.score_one(rec, answer(["G2.1"]), NO)[0] == "undetermined"
    assert S.score_one(rec, answer(["G3.1"]), NO)[0] == "wrong"


def test_base_comparison_drops_versions():
    rec = line(expect=["ENSSSCG00000000001"], cmp="base")
    assert S.score_one(rec, answer(["ENSSSCG00000000001.4"]), NO)[0] == "right"
    rec = line(expect=["ENSSSCG00000000001.3"], cmp="versioned")
    assert S.score_one(rec, answer(["ENSSSCG00000000001.4"]), NO)[0] == "wrong"


def test_locus_rule_by_place():
    loc = dict(assembly="asm1", seqname="1", start=100, end=200, strand="+")
    places = Places({3: {"NEAR": dict(seqname="1", start=150, end=400, strand="+"),
                         "OTHERSTRAND": dict(seqname="1", start=150, end=400, strand="-"),
                         "FAR": dict(seqname="2", start=150, end=400, strand="+")}})
    rec = line(expect=[], locus=loc)
    assert S.score_one(rec, answer(["NEAR"]), places)[0] == "undetermined"
    assert S.score_one(rec, answer(["OTHERSTRAND"]), places)[0] == "wrong"
    assert S.score_one(rec, answer(["FAR", "NEAR"]), places)[0] == "wrong"


def test_locus_on_another_assembly_without_a_mapped_place_is_open():
    loc = dict(assembly="asm0", seqname="1", start=100, end=200, strand="+")
    places = Places({3: {"FAR": dict(seqname="2", start=1, end=9, strand="+")}})
    out, why = S.score_one(line(expect=[], locus=loc), answer(["FAR"]), places)
    assert out == "undetermined" and "not comparable" in why
    mapped = dict(assembly="asm1", seqname="1", start=100, end=200, strand="+", via="asm0 -> asm1", aligned=0.9)
    assert S.score_one(line(expect=[], locus=loc, locus_mapped=mapped), answer(["FAR"]), places)[0] == "wrong"


def test_copies_must_be_reported_apart():
    rec = line(expect=["G1.1"], apart=["COPY.1"])
    assert S.score_one(rec, answer(["G1.1"], apart=["COPY.1"]), NO)[0] == "right"
    assert S.score_one(rec, answer(["G1.1"], apart=[]), NO)[0] == "wrong"
    assert S.score_one(rec, answer(["G1.1"]), NO) == ("right", S.APART_UNCHECKED)


def test_support_rule_leaves_the_better_supported_part_open():
    rec = line(expect=["G1.1", "G2.1"], rule="set_exact", support={"G1.1": 5, "G2.1": 1})
    assert S.score_one(rec, answer(["G1.1"]), NO)[0] == "undetermined"
    assert S.score_one(rec, answer(["G2.1"]), NO)[0] == "wrong"  # the less supported gene alone
    tied = line(expect=["G1.1", "G2.1"], rule="set_exact", support={"G1.1": 2, "G2.1": 2})
    assert S.score_one(tied, answer(["G1.1"]), NO)[0] == "wrong"


def test_two_spellings_of_one_value_count_once():
    rec = line(expect=["IFN-DELTA-8"], spellings=[["IFN-DELTA-8", "IFN-delta-8"]], cmp="exact")
    assert S.score_one(rec, answer(["IFN-delta-8"]), NO)[0] == "right"


def test_timeouts_and_errors_are_their_own_outcomes():
    assert S.score_one(line(expect=["G"]), {"timeout": 30}, NO)[0] == "timeout"
    assert S.score_one(line(expect=["G"]), {"error": "KeyError: x"}, NO)[0] == "error"


def ans(outcome="answer", got=(), **kw):
    d = dict(outcome=outcome, got=list(got), apart=None, flags=None, error_type=None, error_mro=None, error=None)
    d.update(kw)
    return d


def test_expect_error_accepts_the_promised_class_or_a_subclass_only():
    rec = line(expect_error="UnknownReleaseError|d4c2dda:ValueError")
    sub = ans("error", error_type="Sub", error_mro=["Sub", "UnknownReleaseError", "InputError"])
    assert SC.score_line(rec, "best", sub, NO)[0] == "expected_error"
    other = ans("error", error_type="KeyError", error_mro=["KeyError", "LookupError"])
    assert SC.score_line(rec, "best", other, NO)[0] == "wrong"
    assert SC.score_line(rec, "best", ans(got=["G"]), NO)[0] == "wrong"


def test_safe_rule_nothing_or_the_gene_is_right():
    rec = line(rule="safe", expect=["G1.1"])
    assert SC.score_line(rec, "best", ans(got=[]), NO) == ("right", "nothing")
    assert SC.score_line(rec, "best", ans(got=["G1.1"]), NO) == ("right", "recovered")
    assert SC.score_line(rec, "best", ans(got=["G9.1"]), NO)[0] == "wrong"


def test_open_lines_are_undetermined_but_errors_stay_errors():
    rec = line(expect=["G1.1"], open="the walk and REST differ")
    assert SC.score_line(rec, "best", ans(got=["G9.1"]), NO)[0] == "undetermined"
    assert SC.score_line(rec, "best", ans("error", error="x"), NO)[0] == "error"


def test_relations():
    rel = line(rule="relation", relation=dict(op="equal", **{"with": "q2"}))
    assert SC.relation_outcome(rel, "best", ans(got=["A"]), ans(got=["A"]))[0] == "right"
    assert SC.relation_outcome(rel, "best", ans(got=["A"]), ans(got=["B"]))[0] == "wrong"
    sub = line(rule="relation", relation=dict(op="subset", **{"with": "q2"}))
    assert SC.relation_outcome(sub, "best", ans(got=["A"]), ans(got=["A", "B"]))[0] == "right"


def test_the_fast_path_key_matches_the_rule():
    assert SC.norm_key(["ENSG00000000001.2"], "base") == "ENSG00000000001"
    assert SC.norm_key(["b", "a", "a"], "versioned") == SC.SEP.join(["a", "b"])


def test_wilson_matches_between_harness_and_report():
    from idtrack_reproducibility.families import report as R
    assert R.wilson(37, 40) == H.wilson(37, 40)
    assert all(math.isnan(x) for x in R.wilson(0, 0))


def test_a_synonym_read_at_a_release_whose_synonyms_the_truth_did_not_read_is_undetermined():
    rec = dict(qid="q", kind="k", rule="safe", cmp="versioned", expect=["G1.1"], names_read_at=[98, 115])
    a = ans(got=["G9.1"], match="synonym:exact", from_release=97)
    assert SC.score_line(rec, "best", a, NO)[0] == "undetermined"
    assert SC.score_line(rec, "best", dict(a, from_release=98), NO)[0] == "wrong"  # the truth read that release
    assert SC.score_line(rec, "best", dict(a, match="exact"), NO)[0] == "wrong"     # a name: read at every release
    assert SC.score_line(dict(rec, names_read_at=None), "best", a, NO)[0] == "wrong"
    # through the separator or case step another string is read: judged as it is (the fault IN-3, not the truth's reach)
    assert SC.score_line(rec, "best", dict(a, match="synonym:separator"), NO)[0] == "wrong"
    assert SC.score_line(rec, "best", dict(a, match="synonym:case"), NO)[0] == "wrong"
