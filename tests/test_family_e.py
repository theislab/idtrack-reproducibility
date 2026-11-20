"""Family E (robustness): the encoding of non-string inputs and the folding of names, on invented values."""
import math

from idtrack_reproducibility.families.e import adapter_e as AE
from idtrack_reproducibility.families.e import gen_e as GE
from idtrack_reproducibility.families.e import truth_e as TE


def test_non_string_values_round_trip_through_their_encoding():
    for v in (float("nan"), float("inf"), float("-inf")):
        back = AE.decode(GE.enc(v))
        assert (math.isnan(back) and math.isnan(v)) or back == v
    assert GE.enc("TP53") == "TP53" and AE.decode("TP53") == "TP53"
    assert AE.decode({"$": "bytes_hex", "v": "414243"}) == b"ABC"


def test_names_are_folded_as_the_case_and_separator_steps_fold():
    assert TE.fold("Hla-A") == TE.fold("HLA_A") == TE.fold("hla.a") == "hla-a"
    assert TE.fold("TP53") != TE.fold("TP5-3")


def test_expected_error_classes_name_both_api_generations():
    for e in (GE.E_RELEASE, GE.E_DB, GE.E_INVALID, GE.E_UNKNOWN_ID):
        assert "|" in e and e.split("|")[1].startswith("d4c2dda:")


def test_the_decided_judge_takes_every_later_revision_of_its_sets(monkeypatch):
    from idtrack_reproducibility.families import lines as L
    from idtrack_reproducibility.families.e import score_e as SE
    monkeypatch.setattr(L, "set_names", lambda family, organism: ["e1", "e1_dec", "e5", "e5_r1", "e1_dec2", "e8_identity",
                                                                  "e8_identity_dec", "e2", "e2_dec2", "e2_r1"])
    assert SE.sets_for("x", "now") == ["e1", "e5", "e5_r1", "e8_identity", "e2", "e2_r1"]
    # the decided revisions come after every set of the judge `now`: `e2_dec2` (the decided line, revised alike)
    # stays the truth of its qids although `e2_r1` (the same revision for the judge `now`) was published after it
    assert SE.sets_for("x", "decided") == ["e1", "e5", "e5_r1", "e8_identity", "e2", "e2_r1", "e1_dec", "e1_dec2",
                                           "e8_identity_dec", "e2_dec2"]
    assert SE.sets_for("x", "decided", only=["e1"]) == ["e1", "e1_dec", "e1_dec2"]


def test_a_string_carried_in_two_letter_cases_is_the_graphs_merged_spelling(monkeypatch):
    from idtrack_reproducibility.families.e import causes_e as CE
    monkeypatch.setattr(CE, "case_carriers", lambda org, rel: {"gc": {("ENSMUSG00000000001", "GC"), ("ENSMUSG00000000002", "gc")}})
    line = dict(id="GC", to_release=115, expect=["ENSMUSG00000000001.1"], organism="mus_musculus")
    assert "ENSMUSG00000000002:gc" in CE.spellings_merged(line, ["ENSMUSG00000000001.1", "ENSMUSG00000000002.3"], {"from_release": 115})
    # an extra gene that carries the string in no case is not explained so
    assert CE.spellings_merged(line, ["ENSMUSG00000000001.1", "ENSMUSG00000000009.1"], {"from_release": 115}) is None
    # the lower-case form answered by the case step where its exact carrier was expected
    low = dict(id="gc", to_release=115, expect=["ENSMUSG00000000002.3"], organism="mus_musculus")
    assert CE.spellings_merged(low, ["ENSMUSG00000000001.1"], {"from_release": 115})
