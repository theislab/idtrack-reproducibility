"""The record (docs/records.md): one row per line, view and subject, with every column; the corpus's calls, H's timed calls
and another tool's raw responses in the same terms."""
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
    run = scored_folder(data_root, 1, [row("q1", "right", got="G1.1", expected="G1.1", query='{"id":"TP53"}',
                                           status="converted"),
                                       row("q2", "wrong", got="G9.1", expected="G2.1", query='{"id":"X"}')], [])
    df = R.load(folders=[str(run)])
    assert list(df.columns) == R.COLUMNS and len(df) == 2
    r1 = df[df.qid == "q1"].iloc[0]
    assert r1.got == "G1.1" and r1.status == "converted" and r1.expected == "G1.1" and '"id":"TP53"' in r1.query
    assert r1.subject_kind == "tree" and r1.subject == "d4c2dda" and r1.graph_digest == "abc" and r1.python == "3.11.12"
    r2 = df[df.qid == "q2"].iloc[0]
    assert r2.cause == "external_bridge" and r2.cause_class == "b" and r2.outcome == "wrong"
    assert set(df.outcome) <= set(R.OUTCOMES)


def test_a_scoring_of_an_earlier_record_version_is_refused(data_root):
    import pytest
    v0 = scored_folder(data_root, 0, [row("q1", "right")], [])
    with pytest.raises(SystemExit, match="rescore"):
        R.load(folders=[str(v0)])


def test_the_corpus_and_h_calls_are_rows_of_the_record(data_root):
    corpus = data_root / "corpus"
    (corpus / "corpus").mkdir(parents=True)
    (corpus / "corpus" / "sus_scrofa.jsonl").write_text(json.dumps(dict(qid="c1", entry="convert_identifier",
                                                                         stratum="gene_latest", args={"identifier": "X"})) + "\n")
    run = corpus / "runs" / "d4c2dda" / "sus_scrofa"
    run.mkdir(parents=True)
    with gzip.open(run / "run1.jsonl.gz", "wt") as fo:
        fo.write(json.dumps({"meta": {"start": "2026-11-17 15:50:12"}}) + "\n")
        fo.write(json.dumps(dict(qid="c1", entry="convert_identifier", ok=True, dt=0.01,
                                 result={"target_id": ["G1.1"], "match": "exact"})) + "\n")
        fo.write(json.dumps({"meta_end": {}}) + "\n")
    df = R.load_corpus(str(run / "run1.jsonl.gz"))
    assert list(df.columns) == R.COLUMNS and len(df) == 1
    r = df.iloc[0]
    assert (r.part, r.subject, r.outcome, r.got, r.stratum, r.seconds) == ("corpus", "d4c2dda", "unjudged", "G1.1",
                                                                            "gene_latest", 0.01)
    from idtrack_reproducibility.families.h import calls as HC
    from idtrack_reproducibility.families.h import common as HH
    rep = dict(tree="d4c2dda", commit="d4c2dda25b61", organism="sus_scrofa", python="3.11.12", date="2026-11-19")
    hdir = data_root / "performance" / "runs" / "d4c2dda" / "sus_scrofa"
    out = HH.dump_calls(str(hdir / "calls_rep1.json"), [
        HC.call_row(rep, "A", "h|A|X", None, {"id": "X"}, 0.002, {"no_conversion": False, "tie": False}),
        HC.call_row(rep, "F", "c1", "gene_latest", {"identifier": "X"}, 0.003, None, "ValueError")])
    assert out.endswith("calls_rep1.calls.parquet")
    h = R.load_h(str(hdir))
    import pandas as pd
    assert list(h.columns) == R.COLUMNS and list(h.kind) == ["A", "F"] and h.status.iloc[0] == "converted"
    assert pd.isna(h.status.iloc[1]) and pd.isna(h.error_type.iloc[0]) and h.error_type.iloc[1] == "ValueError"
    assert set(h.outcome) == {"unjudged"} <= set(R.OUTCOMES)


def test_another_tools_raw_responses_have_an_index(data_root):
    folder = data_root / "families" / "g" / "raw" / "tool" / "1.0" / "sus_scrofa"
    folder.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([dict(request_id="r1", tool="tool", tool_version="1.0", organism="sus_scrofa",
                                              endpoint="https://x", method="POST", url="https://x/convert", params="{}",
                                              body='{"ids":["TP53"]}', qids=["q1"], asked_at="2026-11-20", http_status=200,
                                              seconds=0.5, bytes=10, sha256="0" * 64, file="r1.json", error=None)],
                                        schema=R.raw_schema()), folder / "index.parquet")
    df = R.load_raw(str(folder))
    assert list(df.columns) == R.RAW_COLUMNS and df.iloc[0].qids.tolist() == ["q1"]


def test_another_tool_is_a_subject_of_its_own(data_root):
    run = scored_folder(data_root, 1, [row("q1", "not_supported", got="", expected="G1.1", query="{}")], [])
    info = json.load(open(run / "run.json"))
    info["subject"] = dict(kind="tool", name="gprofiler", version="e113", asked_at="2026-11-20", endpoint="https://x")
    json.dump(info, open(run / "run.json", "w"))
    df = R.load(folders=[str(run)])
    r = df.iloc[0]
    assert (r.subject_kind, r.subject, r.subject_version, r.asked_at, r.outcome) == \
        ("tool", "gprofiler", "e113", "2026-11-20", "not_supported")


def test_the_columns_come_from_one_table_that_the_document_shows():
    assert R.COLUMNS == [c for c, *_ in R.SCHEMA] and len(set(R.COLUMNS)) == len(R.COLUMNS)
    assert [f.name for f in R.schema()] == R.COLUMNS
    doc = open(os.path.join(os.path.dirname(__file__), "..", "docs", "records.md"), encoding="utf-8").read()
    assert R.schema_table() in doc  # docs/records.md, "Columns": regenerate with `idtrack-repro records --schema-table`


def exported(data_root, tool=False):
    a = scored_folder(data_root, 1, [row("q1", "right", got="G1.1", expected="G1.1", query='{"id":"TP53"}', n_got=1,
                                         seconds=0.25, calls=2, status="converted"),
                                     row("q2", "wrong", got="G9.1", expected="G2.1", n_got=None, seconds=None)], [])
    folders = [str(a)]
    if tool:  # family G's rows: another tool as the subject, written by the same scorer
        g = data_root / "g"
        g.mkdir()
        b = scored_folder(g, 1, [row("q1", "not_supported", got="", expected="G1.1", n_got=0)], [])
        info = json.load(open(b / "run.json"))
        info["subject"] = dict(kind="tool", name="g:Profiler", version="e113", asked_at="2026-11-20", endpoint="https://x")
        json.dump(info, open(b / "run.json", "w"))
        folders.append(str(b))
    out = data_root / "export"
    return out, R.export(str(out), folders=folders), folders


def test_the_export_is_partitioned_typed_indexed_and_read_back(data_root):
    out, index, folders = exported(data_root, tool=True)
    files = sorted(r["file"] for r in index)
    assert files == ["part=t/organism=sus_scrofa/subject=d4c2dda/records.parquet",
                     "part=t/organism=sus_scrofa/subject=g_Profiler/records.parquet"]
    for r in index:
        t = pq.read_table(out / r["file"])
        assert t.schema.equals(R.schema(), check_metadata=False) and t.num_rows == r["rows"]
        import hashlib
        assert hashlib.sha256(open(out / r["file"], "rb").read()).hexdigest() == r["sha256"]
    doc = json.load(open(out / "schema.json"))
    assert doc["record_version"] == R.RECORD_VERSION
    assert [(c["name"], c["type"]) for c in doc["columns"]] == [(c, t) for c, t, *_ in R.SCHEMA]
    assert {s["part"] for s in doc["kept_as_summaries"]} == {"a0", "e7", "e9"}
    tsv = (out / "index.tsv").read_text().splitlines()
    assert tsv[0].split("\t") == R.INDEX_COLS and len(tsv) == 3
    back = R.load(export=str(out))
    direct = R.load(folders=folders)
    assert list(back.columns) == R.COLUMNS and len(back) == len(direct) == 3
    assert str(back.n_got.dtype) == "Int32" and back.n_got.isna().sum() == 1
    tree = back[back.subject == "d4c2dda"].sort_values("qid")
    assert tree.seconds.tolist()[0] == 0.25 and tree.calls.tolist()[0] == 2 and tree.cause.isna().tolist() == [True, False] and \
        tree.cause.tolist()[1] == "external_bridge"
    tool = back[back.subject_kind == "tool"].iloc[0]  # passes through unchanged
    assert (tool.subject, tool.subject_version, tool.asked_at, tool.endpoint, tool.outcome) == \
        ("g:Profiler", "e113", "2026-11-20", "https://x", "not_supported")
    only = R.load(export=str(out), subjects=["g:Profiler"])
    assert len(only) == 1 and R.load(export=str(out), parts=["x"]).empty


def test_the_export_command_writes_the_folder(data_root):
    a = scored_folder(data_root, 1, [row("q1", "right")], [])
    out = data_root / "ex"
    assert R.main(["--folders", str(a), "--out", str(out)]) == 0
    assert (out / "schema.json").exists() and (out / "index.tsv").exists()
    assert len(R.load(export=str(out), organisms=["sus_scrofa"])) == 1
