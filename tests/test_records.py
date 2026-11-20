"""The record (docs/records.md): one row per line, view and subject, with every column, from scorings of either record
version."""
import gzip
import json
import os

import pyarrow as pa
import pyarrow.parquet as pq

from idtrack_reproducibility import records as R
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import scorer as SC


def scored_folder(root, version, outcome_rows, answers):
    run = root / f"run_v{version}"
    (run / "answers").mkdir(parents=True)
    (run / "outcomes").mkdir()
    pq.write_table(pa.Table.from_pylist(answers), run / "answers" / "shard_0_of_1.0.parquet")
    pq.write_table(pa.Table.from_pylist(outcome_rows, schema=SC.out_schema()), run / "outcomes" / "part_0_of_1.parquet")
    json.dump(dict(record_version=version, family="t", organism="sus_scrofa", run=str(run), commit="d4c2dda25b61",
                   graph_digests=["abc"], sets={"t_v1": "x"}, scored_at="2026-11-19T20:00:00",
                   load=[dict(python="3.11.12")]), open(run / "run.json", "w"))
    with gzip.open(run / "causes.tsv.gz", "wt") as fo:
        fo.write("organism\tqid\tview\tkind\tstratum\tunit\tcause\tclass\tdetail\n")
        fo.write("sus_scrofa\tq2\tbest\tk\ts\tu\texternal_bridge\tb\t-\n")
    return run


def row(qid, outcome, **kw):
    base = dict(qid=qid, view="best", family="t", stratum="s", kind="k", variant="best,default", unit="u", set="t_v1",
                outcome=outcome, reason="", n_got=1, got=None, seconds=0.1, calls=1, rule="published", query=None,
                expected=None, broad=None, status=None, match=None, error_type=None)
    base.update(kw)
    return base


def test_every_row_holds_the_query_the_truth_and_the_answer(data_root):
    L.publish("t", "sus_scrofa", "t_v1", [dict(qid="q1", kind="k", id="TP53", from_release=None, to_release=115,
                                                final_database=None, strategy="best", rule="set", cmp="versioned",
                                                expect=["G1.1"], broad=None)], "a test set", __file__)
    ans = [dict(qid="q1", view="best", got=["G1.1"], status="converted", match="exact", error_type=None)]
    v0 = scored_folder(data_root, 0, [row("q1", "right"), row("q2", "wrong", got="G9.1", expected="G2.1",
                                                              query='{"id":"X"}')], ans)
    df = R.load(folders=[str(v0)], complete=True)
    assert list(df.columns) == R.COLUMNS and len(df) == 2
    r1 = df[df.qid == "q1"].iloc[0]
    assert r1.got == "G1.1" and r1.status == "converted" and r1.expected == "G1.1" and '"id":"TP53"' in r1.query
    assert r1.subject_kind == "tree" and r1.subject == "d4c2dda" and r1.graph_digest == "abc" and r1.python == "3.11.12"
    r2 = df[df.qid == "q2"].iloc[0]
    assert r2.cause == "external_bridge" and r2.cause_class == "b" and r2.outcome == "wrong"
    assert set(df.outcome) <= set(R.OUTCOMES)


def test_another_tool_is_a_subject_of_its_own(data_root):
    run = scored_folder(data_root, 1, [row("q1", "not_supported", got="", expected="G1.1", query="{}")], [])
    info = json.load(open(run / "run.json"))
    info["subject"] = dict(kind="tool", name="gprofiler", version="e113", asked_at="2026-11-20", endpoint="https://x")
    json.dump(info, open(run / "run.json", "w"))
    df = R.load(folders=[str(run)])
    r = df.iloc[0]
    assert (r.subject_kind, r.subject, r.subject_version, r.asked_at, r.outcome) == \
        ("tool", "gprofiler", "e113", "2026-11-20", "not_supported")
