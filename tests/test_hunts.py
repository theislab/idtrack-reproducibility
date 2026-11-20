"""Tests of the weakness hunts' reproducers (`idtrack_reproducibility.hunts.reproduce`) and their findings tables.

No graph and no data root: the reproducers' logic is tested on invented answers, the table and the findings tables
for their form and their agreement.
"""
import csv
import os
import re

import pytest

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.hunts import reproduce as R

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLUMNS = ["id", "severity", "organisms", "how_many", "reproducer", "expected", "got", "cause", "disposition",
           "measured_by"]


def rule_names():
    text = open(os.path.join(ROOT, "docs", "rules.md"), encoding="utf-8").read()
    return set(re.findall(r"^### ([a-z0-9_]+)$", text, flags=re.M))


def corrections():
    """The judge's corrections named by rule `judge_corrections` (A*'s revision `copies_apart_early` among them)."""
    text = open(os.path.join(ROOT, "docs", "rules.md"), encoding="utf-8").read()
    part = text.split("### judge_corrections", 1)[1]
    return set(re.findall(r"`([a-z0-9_]+)`", part))


def findings_table(angle):
    with open(os.path.join(ROOT, "hunts", angle, "findings.tsv"), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))  # quotes are part of the text


# ---------------- compare: one call's answer against the recorded one ----------------

def test_same_answer_is_as_recorded():
    rec = dict(status="converted", match="exact", got=["ENSG00000000001.2"])
    obs = dict(outcome="answer", status="converted", match="exact", got=["ENSG00000000001.2"], apart=[])
    assert R.compare(rec, obs) == []


def test_other_gene_is_a_change():
    rec = dict(status="converted", got=["ENSG00000000001.2"])
    obs = dict(outcome="answer", status="converted", got=["ENSG00000000009.1"], apart=[])
    diff = R.compare(rec, obs)
    assert len(diff) == 1 and diff[0].startswith("got:")


def test_ids_compared_as_sets_not_in_order():
    rec = dict(got=["ENSG00000000002.1", "ENSG00000000001.1"])
    assert R.compare(rec, dict(got=["ENSG00000000001.1", "ENSG00000000002.1"])) == []
    assert R.compare(rec, dict(got=["ENSG00000000001.1"])) != []  # one of two left out is a change


def test_got_base_ignores_versions_only():
    rec = dict(got_base=["ENSG00000000001"])
    assert R.compare(rec, dict(got=["ENSG00000000001.7"])) == []
    assert R.compare(rec, dict(got=["ENSG00000000001.7", "ENSG00000000002.1"])) != []
    assert R.compare(rec, dict(got=[])) != []


def test_got_has_is_a_subset_and_n_got_a_count():
    assert R.compare(dict(got_has=["A.1"]), dict(got=["A.1", "B.1"])) == []
    assert R.compare(dict(got_has=["A.1", "C.1"]), dict(got=["A.1", "B.1"])) != []
    assert R.compare(dict(n_got=3), dict(got=["a", "b", "c"])) == []
    assert R.compare(dict(n_got=3), dict(got=["a", "b"])) != []


def test_status_match_and_copies_apart():
    obs = dict(status="tie", match="synonym:exact", got=["a", "b"], apart=["c"])
    assert R.compare(dict(status="tie", match="synonym:exact", apart_any=True), obs) == []
    assert R.compare(dict(status="converted"), obs) != []
    assert R.compare(dict(match="exact"), obs) != []
    assert R.compare(dict(apart_any=False), obs) != []


def test_an_error_where_an_answer_was_recorded_is_a_change():
    obs = dict(outcome="error", error_type="ValueError", error="ValueError: x", got=[])
    assert R.compare(dict(got=[]), obs) != []  # the same ids (none), but the call raised
    assert R.compare(dict(outcome="error", error_type="ValueError"), obs) == []
    fixed = dict(outcome="error", error_type="InvalidIdentifierError", error="...", got=[])
    assert R.compare(dict(outcome="error", error_type="ValueError"), fixed) != []  # the fix: an input error


def test_a_recorded_error_that_now_answers_is_a_change():
    assert R.compare(dict(outcome="error", error_type="ValueError"), dict(outcome="answer", got=[], status="not_found")) != []


def test_truth_fact():
    assert R.compare(dict(present=False), dict(present=False)) == []
    assert R.compare(dict(present=False), dict(present=True)) != []
    missing = dict(outcome="error", error_type="MissingTruth", error="no file")
    assert R.compare(dict(present=True), missing) != []


def test_unknown_recorded_field_is_refused():
    with pytest.raises(ValueError):
        R.compare(dict(gene="x"), dict(got=[]))


# ---------------- evaluate and the exit code, on invented reproducers ----------------

def invented(recorded_a, recorded_b):
    return [dict(angle="inputs", id="IN-0", title="invented", calls=[
        R.conv("homo_sapiens", "X", recorded_a, "nothing", "a test"),
        R.conv("mus_musculus", "Y", recorded_b, "nothing", "a test", from_release=100)])]


def test_finding_holds_when_every_call_answers_as_recorded():
    sel = invented(dict(status="converted", got=["G1.1"]), dict(status="not_found", got=[]))
    answers = {"X": dict(outcome="answer", status="converted", got=["G1.1"]),
               "Y": dict(outcome="answer", status="not_found", got=[])}
    res = R.evaluate(sel, lambda f, i, c: answers[c["id"]])
    assert res[0]["holds"] and R.exit_code(res) == 0


def test_a_reproducer_whose_answer_changed_exits_1():
    sel = invented(dict(status="converted", got=["G1.1"]), dict(status="not_found", got=[]))
    answers = {"X": dict(outcome="answer", status="converted", got=["G1.1"]),
               "Y": dict(outcome="answer", status="converted", got=["G2.1"])}  # the second call is fixed or changed
    res = R.evaluate(sel, lambda f, i, c: answers[c["id"]])
    assert not res[0]["holds"]
    assert [bool(r["differences"]) for r in res[0]["calls"]] == [False, True]
    assert R.exit_code(res) == 1


def test_one_changed_finding_among_several_exits_1():
    a = invented(dict(got=["G1.1"]), dict(got=[]))[0]
    b = dict(a, id="IN-00")
    same = lambda f, i, c: dict(got=["G1.1"] if c["id"] == "X" else [])  # noqa: E731
    assert R.exit_code(R.evaluate([a, b], same)) == 0
    changed = lambda f, i, c: dict(got=["G9.1"]) if f["id"] == "IN-00" and c["id"] == "X" else same(f, i, c)  # noqa: E731
    res = R.evaluate([a, b], changed)
    assert [r["holds"] for r in res] == [True, False] and R.exit_code(res) == 1


def test_report_marks_changes(capsys):
    sel = invented(dict(got=["G1.1"]), dict(got=[]))
    res = R.evaluate(sel, lambda f, i, c: dict(outcome="answer", got=["G3.1"]))
    import io
    buf = io.StringIO()
    R.report(res, "atree", buf)
    text = buf.getvalue()
    assert "CHANGED" in text and "inputs IN-0 (invented): CHANGED on atree" in text


# ---------------- selection and description ----------------

def test_select_by_angle_id_and_organism():
    assert [f["id"] for f in R.select("inputs", "IN-1")] == ["IN-1"]
    mouse = R.select("inputs", "IN-1", "mus_musculus")[0]["calls"]
    assert mouse and all(c["organism"] == "mus_musculus" for c in mouse)
    assert {f["angle"] for f in R.select("all", "all")} == set(R.ANGLES)
    with pytest.raises(SystemExit):
        R.select("inputs", "IN-999")
    with pytest.raises(SystemExit):
        R.select("names", "all")


def test_describe_writes_ascii_python():
    c = R.conv("homo_sapiens", "﻿TP53", {}, "x", "y", to_release=115)
    assert R.describe(c) == 'convert("\\ufeffTP53", to_release=115)'
    long = R.conv("homo_sapiens", "ENSG00000141510." + "9" * 5000, {}, "x", "y")
    assert R.describe(long) == 'convert("ENSG00000141510." + "9" * 5000)'
    allv = R.conv("homo_sapiens", "A", {}, "x", "y", view="all", from_release=75)
    assert R.describe(allv) == 'convert("A", from_release=75, strategy="all")'


def test_list_needs_no_graph(capsys):
    assert R.main(["--list"]) == 0
    assert "IN-1" in capsys.readouterr().out


def test_submit_one_job_per_organism_with_every_tree(tmp_path, monkeypatch):
    toml = tmp_path / "idtrack-repro.toml"
    toml.write_text('[data]\nroot = "data"\n[trees.d4c2dda]\nsource = "t1"\n[trees.71690ab]\nsource = "t2"\n'
                    '[jobs]\nrunner = "local"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(toml))
    CFG.reset()
    from idtrack_reproducibility import jobs
    seen = []
    monkeypatch.setattr(jobs, "submit", lambda what, cmd, log, **kw: seen.append((what, cmd, kw)) or f"j{len(seen)}")
    try:
        assert R.main(["inputs", "IN-1", "--submit", "--tree", "d4c2dda,71690ab", "--out", str(tmp_path / "o")]) == 0
    finally:
        CFG.reset()
    orgs = sorted({c["organism"] for c in R.select("inputs", "IN-1")[0]["calls"]})
    assert [w for w, _, _ in seen] == [f"hunts_inputs_in_1_{o}" for o in orgs]
    for (what, cmd, kw), org in zip(seen, orgs):
        assert "--tree d4c2dda" in cmd and "--tree 71690ab" in cmd and f"--organism {org}" in cmd
        assert kw["mem_gb"] == CFG.MEMORY_GB[org]


# ---------------- the table and the findings tables agree ----------------

def test_table_is_well_formed():
    seen = set()
    for f in R.FINDINGS:
        assert f["angle"] in R.ANGLES and (f["angle"], f["id"]) not in seen and f["calls"] and f["title"]
        seen.add((f["angle"], f["id"]))
        for c in f["calls"]:
            assert c["kind"] in ("convert", "infer", "truth")
            assert c["organism"] in CFG.ORGANISMS
            assert c["expected"] and c["source"] and c["recorded"]
            assert set(c["recorded"]) <= set(R.FIELDS)
            if c["kind"] == "convert":
                assert c["view"] in ("best", "all")
                assert set(c["kwargs"]) <= {"from_release", "to_release", "final_database", "from_database"}


@pytest.mark.parametrize("angle", R.ANGLES)
def test_findings_table_form(angle):
    rows = findings_table(angle)
    assert rows and list(rows[0]) == COLUMNS
    rules, corr = rule_names(), corrections()
    for r in rows:
        assert all(r[k].strip() for k in COLUMNS), r["id"]
        kind, _, what = r["disposition"].partition(":")
        if kind == "rule":
            assert set(what.split(",")) <= rules, r["id"]
        elif kind == "truth_error":
            assert set(what.split(",")) <= corr, r["id"]
        elif kind == "limitation":
            assert re.match(r"^L-[a-z0-9-]+$", what), r["id"]
        else:
            assert r["disposition"] == "not_a_defect", r["id"]


@pytest.mark.parametrize("angle", R.ANGLES)
def test_every_finding_with_a_reproducer_is_in_the_table(angle):
    table = {f["id"] for f in R.FINDINGS if f["angle"] == angle}
    rows = findings_table(angle)
    for r in rows:
        runnable = r["reproducer"].startswith(f"reproduce {angle} {r['id']}")
        assert runnable == (r["id"] in table), r["id"]
    assert table <= {r["id"] for r in rows}
