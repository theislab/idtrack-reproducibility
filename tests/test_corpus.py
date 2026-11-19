"""The regression corpus: its canonical answer form, the comparator's classes of difference and the steps that read
answer files, on invented answers (no data root, no idtrack)."""
import gzip
import json
import math

import numpy as np
import pytest

from idtrack_reproducibility.corpus import compare as CMP
from idtrack_reproducibility.corpus import pipeline as P
from idtrack_reproducibility.corpus import run_harmonize as RH
from idtrack_reproducibility.corpus.run_corpus import canon, error_type, norm_message, old_message

# ----------------------------------------------------------------------------------------------- the canonical form


def test_canon_sorts_dict_keys_and_set_items_and_keeps_list_order():
    x = canon({"b": [3, 1, 2], "a": {"z", "x", "y"}})
    assert list(x) == ["a", "b"]
    assert x["a"] == {"__set__": ["x", "y", "z"]}
    assert x["b"] == [3, 1, 2]  # a list is the engine's ranking: its order is part of the answer


def test_canon_writes_non_string_keys_as_sorted_pairs():
    assert canon({2: "b", 1: "a"}) == {"__dict__": [[1, "a"], [2, "b"]]}
    # one string key among numbers is still a dict with non-string keys (pairs sorted by their JSON text)
    assert canon({1: "a", "k": 2}) == {"__dict__": [["k", 2], [1, "a"]]}


def test_canon_writes_nan_and_infinities_as_strings_and_numpy_as_python():
    assert canon([math.nan, math.inf, -math.inf, 1.5]) == ["NaN", "Infinity", "-Infinity", 1.5]
    v = canon({"n": np.int64(3), "f": np.float32(0.5), "arr": np.array([1, 2]), "b": np.bool_(True)})
    assert v == {"arr": [1, 2], "b": True, "f": 0.5, "n": 3}
    assert type(v["n"]) is int and type(v["b"]) is bool


def test_canon_tuples_become_lists_and_unknown_objects_their_repr():
    class Thing:
        def __repr__(self):
            return "Thing()"

    assert canon((1, (2, 3))) == [1, [2, 3]]
    assert canon(Thing()) == {"__repr__": "Thing()"}
    assert canon(len) == "<callable len>"
    # bool is not turned into int, None stays None
    assert canon([True, None, 0]) == [True, None, 0]


def test_canon_is_stable_under_insertion_order():
    a = {"x": {3, 1}, "y": {"q": 1, "p": 2}}
    b = {"y": {"p": 2, "q": 1}, "x": {1, 3}}
    assert json.dumps(canon(a)) == json.dumps(canon(b))


def test_norm_message_sorts_set_literals_leaves_dicts_and_drops_addresses():
    assert norm_message("unknown {'b', 'a'} at 0x7f00ab") == "unknown {'a', 'b'}"
    assert norm_message("got {'k': 1, 'a': 2}") == "got {'k': 1, 'a': 2}"  # a dict literal is left as it is
    assert norm_message("no braces") == "no braces"


def test_typed_errors_are_written_as_the_class_d4c2dda_raised():
    class InputError(ValueError):
        pass

    class UnknownVersionError(Exception):
        pass

    assert error_type(InputError("x")) == "ValueError"  # nearest built-in base
    assert error_type(UnknownVersionError("x")) == "UnknownVersionError"  # a class name d4c2dda had itself
    assert error_type(KeyError("x")) == "KeyError"


def test_renamed_arguments_are_named_as_d4c2dda_named_them():
    assert old_message("`allow_bridges` must be a bool", "convert_identifier") == "`go_external` must be a bool"
    assert old_message("`identifier` is empty", "track_convert") == "`from_id` is empty"
    assert old_message("`identifier` is empty", "convert_identifier") == "`identifier` is empty"  # only the engine's


# -------------------------------------------------------------------------------------------- the comparator classes


def ok(entry, result):
    return {"qid": "q", "entry": entry, "ok": True, "result": result}


def err(entry, typ, msg):
    return {"qid": "q", "entry": entry, "ok": False, "error": {"type": typ, "message": msg}}


def test_the_same_answer_is_no_difference():
    a = ok("convert_identifier", {"target_id": ["G1"], "match": "exact"})
    assert CMP.classify(a, dict(a)) is None


def test_timeouts_are_their_own_class_never_an_answer_difference():
    a, b = ok("convert_identifier", {"target_id": ["G1"]}), err("convert_identifier", "CallTimeout", "exceeded 60 s")
    b["timeout"] = 60
    assert CMP.classify(a, b) == "timeout_b"
    assert CMP.classify(b, a) == "timeout_a"
    assert CMP.classify(b, dict(b)) == "timeout_both"


def test_error_classes():
    e = "convert_identifier"
    assert CMP.classify(ok(e, {"target_id": []}), err(e, "ValueError", "bad")) == "ok_vs_error"
    assert CMP.classify(err(e, "ValueError", "bad"), err(e, "TypeError", "bad")) == "error_type"
    assert CMP.classify(err(e, "ValueError", "bad"), err(e, "ValueError", "worse")) == "error_message"
    # the items of a printed set in another order: the same normalised message
    assert CMP.classify(err(e, "ValueError", "{'a', 'b'}"), err(e, "ValueError", "{'b', 'a'}")) == "error_message_order_only"
    # the same characters in another order (a printed dict): order only, not a new message
    assert CMP.classify(err(e, "ValueError", "{'a': 1, 'b': 2}"),
                        err(e, "ValueError", "{'b': 2, 'a': 1}")) == "error_message_order_only"
    assert CMP.classify(err(e, "ValueError", "same"), err(e, "ValueError", "same")) is None


def test_order_only_is_told_apart_from_a_change():
    e = "infer_identifier_source"
    assert CMP.classify(ok(e, [[1, "x"], [2, "y"]]), ok(e, [[2, "y"], [1, "x"]])) == "order_only"
    ci = "convert_identifier"
    assert CMP.classify(ok(ci, {"target_id": ["A", "B"]}), ok(ci, {"target_id": ["B", "A"]})) == "order_only"


@pytest.mark.parametrize("a,b,kind", [
    (["A", "B"], ["A"], "target_id:subset"),
    (["A"], ["A", "B"], "target_id:superset"),
    (["A"], ["C"], "target_id:changed"),
    (["A", "B"], ["A", "C"], "target_id:changed"),
])
def test_targets_of_convert_identifier(a, b, kind):
    ci = "convert_identifier"
    assert CMP.classify(ok(ci, {"target_id": a, "tie": False}), ok(ci, {"target_id": b, "tie": False})) == kind


def test_the_first_differing_field_is_read_answer_first():
    ci = "convert_identifier"
    a = {"target_id": ["A"], "no_conversion": False, "last_node": [["G", "A"]], "match": "exact"}
    b = {"target_id": ["A"], "no_conversion": False, "last_node": [["H", "A"]], "match": "case"}
    # `match` comes first alphabetically, but how the answer was reached (last_node) is read before the rest
    assert CMP.classify(ok(ci, a), ok(ci, b)) == "value:last_node"
    b2 = dict(a, no_conversion=True, match="case")
    assert CMP.classify(ok(ci, a), ok(ci, b2)) == "value:no_conversion"
    # other entries: alphabetical
    assert CMP.classify(ok("list", {"z": 1, "b": 1}), ok("list", {"z": 2, "b": 2})) == "value:b"
    assert CMP.classify(ok("list", [1]), ok("list", [2])) == "value"


def test_track_convert_classes():
    tc = "track_convert"
    g1 = {"score": 1.0, "the_path": [1], "final_conversion": {"final_elements": {"T1": {}}, "confidence": 0}}
    assert CMP.classify(ok(tc, {"G1": g1}), ok(tc, {"G2": g1})) == "candidates_changed"
    other_targets = dict(g1, final_conversion={"final_elements": {"T2": {}}, "confidence": 0})
    assert CMP.classify(ok(tc, {"G1": g1}), ok(tc, {"G1": other_targets})) == "candidate:final_elements_changed"
    details = dict(g1, final_conversion={"final_elements": {"T1": {}}, "confidence": 1})
    assert CMP.classify(ok(tc, {"G1": g1}), ok(tc, {"G1": details})) == "candidate:final_conversion_details"
    both = dict(g1, score=2.0, the_path=[2])
    assert CMP.classify(ok(tc, {"G1": g1}), ok(tc, {"G1": both})) == "candidate:score+the_path"


def test_compare_records_groups_by_entry_and_class_with_strata():
    A = {"q1": dict(ok("convert_identifier", {"target_id": ["A"]}), qid="q1"),
         "q2": dict(ok("convert_identifier", {"target_id": ["A"]}), qid="q2"),
         "q3": dict(ok("list", [1]), qid="q3")}
    B = {"q1": dict(ok("convert_identifier", {"target_id": ["B"]}), qid="q1"),
         "q2": dict(ok("convert_identifier", {"target_id": ["A"]}), qid="q2"),
         "q4": dict(ok("list", [1]), qid="q4")}
    corpus = {"q1": {"stratum": "gene_split", "args": {"identifier": "X"}}}
    rep = CMP.compare_records(A, B, corpus)
    assert rep["calls"] == 4 and rep["differing"] == 3
    kinds = {(g["entry"], g["kind"]): g for g in rep["groups"]}
    assert set(kinds) == {("convert_identifier", "target_id:changed"), ("list", "missing_in_b"), ("list", "missing_in_a")}
    assert kinds[("convert_identifier", "target_id:changed")]["strata"] == {"gene_split": 1}
    assert kinds[("convert_identifier", "target_id:changed")]["examples"][0]["args"] == {"identifier": "X"}
    assert CMP.text(rep).splitlines()[0] == "DIFFERENT: 3 of 4 calls"
    assert CMP.text(CMP.compare_records(A, A)).splitlines()[0] == "IDENTICAL (3 calls)"


# ------------------------------------------------------------------------------------------ answer files and the steps


def write_run(path, records, meta=None):
    with gzip.open(path, "wt") as fo:
        fo.write(json.dumps({"meta": meta or {"tree": "x"}}) + "\n")
        for r in records:
            fo.write(json.dumps(r) + "\n")
        fo.write(json.dumps({"meta_end": {"end": "t"}}) + "\n")


def test_load_ignores_dt_and_merges_the_meta_lines(tmp_path):
    p = str(tmp_path / "a.jsonl.gz")
    write_run(p, [{"qid": "q1", "entry": "list", "ok": True, "result": [1], "dt": 0.5}])
    answers, meta = CMP.load(p)
    assert answers == {"q1": {"qid": "q1", "entry": "list", "ok": True, "result": [1]}}
    assert meta["meta"] == {"tree": "x"} and meta["meta_end"] == {"end": "t"}


def test_slow_pass_merge_replaces_timeouts_unless_the_rerun_timed_out_again(tmp_path):
    base, rerun, out = (str(tmp_path / f) for f in ("base.jsonl.gz", "slow.jsonl.gz", "merged.jsonl.gz"))
    t1 = {"qid": "q1", "entry": "e", "ok": False, "timeout": 60, "error": {"type": "CallTimeout", "message": "x"}}
    t2 = dict(t1, qid="q2")
    write_run(base, [t1, t2, {"qid": "q3", "entry": "e", "ok": True, "result": 1}])
    assert P.timed_out(base) == ["q1", "q2"]
    write_run(rerun, [{"qid": "q1", "entry": "e", "ok": True, "result": 7},
                      dict(t2, timeout=60)], meta={"tree": "slow"})  # q2 hit a cap not larger than the base's again
    assert P.merge(base, rerun, out) == 1
    merged, meta = CMP.load(out)
    assert merged["q1"]["result"] == 7 and "timeout" in merged["q2"] and merged["q3"]["result"] == 1
    assert meta["meta"] == {"tree": "x"} and meta["meta_rerun"]["meta"] == {"tree": "slow"}
    assert P.timed_out(out) == ["q2"]


def test_gzip_parts_join_into_one_stream(tmp_path):
    parts = []
    for i in range(3):
        p = str(tmp_path / f"p{i}.jsonl.gz")
        write_run(p, [{"qid": f"q{i}", "entry": "harmonize", "ok": True, "result": i}])
        parts.append(p)
    out = str(tmp_path / "all.jsonl.gz")
    P.join_gzip(parts, out)
    answers, _ = CMP.load(out)
    assert sorted(answers) == ["q0", "q1", "q2"]


def test_harmonize_datasets_are_drawn_deterministically_and_leave_out_slow_ids(tmp_path):
    corpus = tmp_path / "c.jsonl"
    lines = []
    for s, n in (("gene_graph_at_release", 500), ("external_per_database", 500), ("nonexistent", 5), ("malformed", 5)):
        lines += [{"qid": f"{s}{i}", "entry": "convert_identifier", "stratum": s, "args": {"identifier": f"{s}_{i}"}}
                  for i in range(n)]
    lines.append({"qid": "nb", "entry": "convert_identifier", "stratum": "malformed", "args": {"identifier": None}})
    corpus.write_text("".join(json.dumps(x) + "\n" for x in lines))
    a = RH.dataset_ids(str(corpus), 7)
    b = RH.dataset_ids(str(corpus), 7)
    assert a == b and set(a) == {"D1", "D2", "D3"}
    assert all(len(v) == len(set(v)) for v in a.values())
    assert set(a["D1"]) <= {f"gene_graph_at_release_{i}" for i in range(500)}
    assert len(a["D3"]) == 400 + 400 + 5 + 5
    c = RH.dataset_ids(str(corpus), 7, exclude={"gene_graph_at_release_0"})
    assert "gene_graph_at_release_0" not in c["D1"] + c["D3"]


def test_harmonize_configurations_follow_the_corpus_summary():
    summary = {"releases": [56, 70, 90, 100, 110, 115], "databases_main": ["EntrezGene", "VGNC Symbol"]}
    cfgs = RH.configurations("sus_scrofa", 115, summary)
    assert cfgs == [("symbol_latest", 115, "VGNC Symbol"), ("ensembl_older", 90, "ensembl_gene"),
                    ("ncbi_latest", 115, "EntrezGene")]
    summary["databases_main"].append("NCBI gene")
    assert RH.configurations("sus_scrofa", 115, summary)[2] == ("ncbi_latest", 115, "NCBI gene")
