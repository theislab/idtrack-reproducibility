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


# ---------------- E8: a batch against single calls, one release deciding a reading (batch_e) ----------------

def test_the_copy_form_is_another_item_with_a_numbered_suffix_never_an_ensembl_version():
    from idtrack_reproducibility.families.e import batch_e as BE
    assert BE.copy_form("SRF.1", {"SRF"}) == "SRF" and BE.copy_form("OCLN_2", {"OCLN"}) == "OCLN"
    assert BE.copy_form("SRF.1", {"MB"}) is None
    assert BE.copy_form("ENSG00000141510.5", {"ENSG00000141510"}) is None


def _display(truth_store, R, rows):
    from tests.conftest import write_gz
    write_gz(str(truth_store / "core_display" / f"{R}.tsv.gz"),
             ["gene_id", "db_name", "db_display_name", "acc", "label", "info_type"],
             [(g, "HGNC", "HGNC", "1", lab, "DIRECT") for g, lab in rows])


def test_batch_lines_relate_each_item_to_its_single_call_and_copies_to_their_item(truth_store, monkeypatch):
    from idtrack_reproducibility.families.e import batch_e as BE
    from tests.conftest import GENE_COLS, write_gz
    write_gz(str(truth_store / "gtf" / "genes_3.tsv.gz"), GENE_COLS, [
        ("ENSSSCG00000000002", 2, "GB", "protein_coding", "1", 1000, 1250, "+", "primary"),
        ("ENSSSCG00000000004", 1, "GB-1", "protein_coding", "1", 100, 180, "+", "primary"),
        ("ENSSSCG00000000006", 1, "GF", "protein_coding", "3", 5, 50, "+", "primary")])
    _display(truth_store, 1, [("ENSSSCG00000000001.1", "GA")])
    _display(truth_store, 3, [("ENSSSCG00000000002.2", "GB"), ("ENSSSCG00000000004.1", "GB-1"), ("ENSSSCG00000000006.1", "GF")])
    monkeypatch.setattr(BE, "SAMPLE", 1)
    monkeypatch.setattr(BE, "MIXED", 1)
    monkeypatch.setattr(BE, "MIXED_N", 20)
    b = BE.Batch("sus_scrofa")
    lines = b.lines()
    by = {q["qid"]: q for q in lines}
    assert len(by) == len(lines)  # every qid once
    assert b.lists["release|3"] == ["GB", "GB-1", "GF"]
    rel = [q for q in lines if q["rule"] == "relation"]
    assert rel and all(q["relation"]["with"] in by and by[q["relation"]["with"]].get("open") for q in rel)
    # Ensembl's own name GB-1 beside GB: asked in both modes, related to its own single call (it is no copy)
    for mode, fr in (("atR", 3), ("default", None)):
        q = by[f"e8b|batch|{mode}|3|GB-1"]
        assert q["stratum"] == f"e8_batch_release_{mode}_suffix" and q["batch_ref"] == {"set": "e8_batch", "key": "release|3"}
        assert b.lists["release|3"][q["pos"]] == "GB-1"
        single = by[q["relation"]["with"]]
        assert single["id"] == "GB-1" and single["from_release"] == fr and single["to_release"] == fr
    # in a mixed list a made-unique copy no table holds is related to its item's single call
    for q in rel:
        if q["stratum"] == "e8_batch_mixed_copy":
            item = by[q["relation"]["with"]]["id"]
            assert q["id"] != item and BE.copy_form(q["id"], {item}) == item
    # the lines of one list and arguments come together (the adapter converts a batch once)
    keys = [(q["batch_ref"]["key"], q["from_release"]) for q in rel]
    assert keys == sorted(keys, key=lambda k: (k[0], "-" if k[1] is None else str(k[1])))


def test_one_release_pairs_strings_whose_reading_changes_between_releases(truth_store):
    from idtrack_reproducibility.families.e import batch_e as BE
    from tests.conftest import write_gz
    _display(truth_store, 1, [("ENSSSCG00000000001.1", "GA"), ("ENSSSCG00000000003.1", "GC")])
    _display(truth_store, 2, [("ENSSSCG00000000001.1", "GA")])
    _display(truth_store, 3, [("ENSSSCG00000000002.2", "GB"), ("ENSSSCG00000000004.1", "GD")])
    syn = ["synonym", "gene_id", "db_name", "label", "object_type"]
    write_gz(str(truth_store / "core_synonyms" / "3.tsv.gz"), syn, [("GA", "ENSSSCG00000000004", "HGNC", "GD", "Gene")])
    write_gz(str(truth_store / "core_synonyms" / "1.tsv.gz"), syn, [("GB", "ENSSSCG00000000001", "HGNC", "GA", "Gene")])
    lines, extra = BE.one_release("sus_scrofa")
    by = {q["qid"]: q for q in lines}
    # GA: a name at 1 and 2, no name and a synonym at the last release
    to = [q for q in lines if q["stratum"] == "e8_onerelease_name_then_synonym"]
    assert len(to) == 1 and to[0]["id"] == "GA" and to[0]["from_release"] is None and to[0]["to_release"] in (1, 2)
    ref = by[to[0]["relation"]["with"]]
    assert ref["from_release"] == ref["to_release"] == to[0]["to_release"] and ref.get("open")
    # GB: a synonym at 1 (no name there), a name at the last release
    back = [q for q in lines if q["stratum"] == "e8_onerelease_synonym_then_name"]
    assert [q["id"] for q in back] == ["GB"] and back[0]["to_release"] == 1
    assert extra["candidates"] == {"name_then_synonym": 1, "synonym_then_name": 1}


def test_a_listed_batch_is_converted_once_per_arguments():
    from idtrack_reproducibility.families.e import adapter_e as AE
    ad = AE.AdapterE.__new__(AE.AdapterE)
    ad.timeout, ad.organism = 30, "sus_scrofa"
    ad._lists_cache = {"e8_batch": {"release|3": ["GB", "GB-1"], "release|2": ["GA"]}}
    calls = []

    def fn(limit):
        calls.append(limit)
        return ["result"], None, 1.0
    line = dict(batch_ref={"set": "e8_batch", "key": "release|3"})
    assert ad._kept_batch(line, {"from_release": 3}, fn) == (["result"], None, 1.0, False)
    assert ad._kept_batch(line, {"from_release": 3}, fn)[3] is True
    assert len(calls) == 1 and calls[0] == 30 + 30 + AE.BATCH_ITEM_SECONDS * 2
    ad._kept_batch(line, {}, fn)  # other arguments: another conversion
    ad._kept_batch(dict(batch_ref={"set": "e8_batch", "key": "release|2"}), {}, fn)
    assert len(calls) == 3
    ad._kept_batch(line, {"from_release": 3}, fn)  # only the two latest are kept
    assert len(calls) == 4
    assert ad._kept_batch(dict(batch=["x"]), {}, lambda lim: (lim, None, 0.0)) == (None, None, 0.0, False)


def test_wrong_relations_are_classed_by_the_answer():
    import re
    from idtrack_reproducibility.families import limits_outcomes as LO
    m = [dict(limitation="L-batch-copies", strata="^e8_batch_", when="copy"),
         dict(limitation="L-batch-other", strata="^e8_batch_", when="any"),
         dict(limitation="L-display-name-not-held", strata="^b9_", when="read_elsewhere")]
    for r in m:
        r["rx"] = re.compile(r["strata"])
    assert LO.classify(m, "e8_batch_release_atR", dict(match="copy:exact"))["limitation"] == "L-batch-copies"
    assert LO.classify(m, "e8_batch_mixed", dict(match="exact"))["limitation"] == "L-batch-other"
    assert LO.classify(m, "b9_early_hgnc", dict(fallback=True, status="converted"))["limitation"] == "L-display-name-not-held"
    assert LO.classify(m, "b9_early_hgnc", dict(status="not_found"))["limitation"] == "L-display-name-not-held"
    assert LO.classify(m, "b9_early_hgnc", dict(fallback=False, status="converted")) is None
