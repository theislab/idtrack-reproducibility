"""Intervals, the report of a scored run, and the known-limitations builder, on invented outcomes."""
import csv
import gzip
import json
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from idtrack_reproducibility.families import limitations as LIM
from idtrack_reproducibility.families import report as R
from idtrack_reproducibility.families import scorer as SC


def test_cluster_bootstrap_is_deterministic_and_bounds_the_rate():
    r = np.array([1, 1, 0, 1, 1, 1, 1, 0, 1, 1])
    d = np.ones(10, dtype=int)
    a = R.cluster_bootstrap(r, d, "fam|org|s|*")
    assert a == R.cluster_bootstrap(r, d, "fam|org|s|*")
    assert a[0] <= 0.8 <= a[1]
    assert R.cluster_bootstrap(np.array([3, 4]), np.array([3, 4]), "k") == (1.0, 1.0)
    assert all(np.isnan(R.cluster_bootstrap(np.array([]), np.array([]), "k")))


def test_clustered_interval_is_wider_than_wilson_when_errors_cluster():
    """50 units with one right line each and one unit with 50 wrong lines: a bootstrap over lines (or Wilson's
    interval, which treats lines as independent) is far narrower than the one over units."""
    r = np.array([1] * 50 + [0])
    d = np.array([1] * 50 + [50])
    lo, hi = R.cluster_bootstrap(r, d, "k")
    wlo, whi = R.wilson(50, 100)
    assert (hi - lo) > 2 * (whi - wlo)


def write_outcomes(folder, rows):
    os.makedirs(os.path.join(folder, "outcomes"), exist_ok=True)
    full = [dict({c: None for c in SC.OUT_COLS}, **r) for r in rows]
    pq.write_table(pa.Table.from_pylist(full, schema=SC.out_schema()), os.path.join(folder, "outcomes", "part_0_of_1.parquet"))
    json.dump(dict(family="x", organism="sus_scrofa", parts=1, missing={}), open(os.path.join(folder, "run.json"), "w"))


def test_report_counts_rates_units_and_missing_text(tmp_path):
    rows = []
    for i in range(8):
        rows.append(dict(qid=f"q{i}", view="best", family="x", stratum="s1", kind="k", variant="best,from", unit=f"u{i % 4}",
                         outcome="right", reason="", n_got=1, seconds=0.001 * (i + 1), calls=1, rule="published"))
    rows.append(dict(qid="q8", view="best", family="x", stratum="s1", kind="k", variant="best,from", unit="u9",
                     outcome="wrong", reason="not expected: G", n_got=1, got="G", seconds=0.5, calls=1, rule="published",
                     query="{}", expected="H", broad=None))
    rows.append(dict(qid="q9", view="best", family="x", stratum="s1", kind="k", variant="best,from", unit="u9",
                     outcome="undetermined", reason="within", n_got=1, got="G", seconds=0.5, calls=1, rule="published",
                     query="{}", expected="H", broad=None))
    write_outcomes(str(tmp_path), rows)
    R.write(str(tmp_path), "x", "sus_scrofa")
    rep = {(r["stratum"], r["view"]): r for r in csv.DictReader(open(tmp_path / "report.tsv"), delimiter="\t")}
    star = rep[("s1", "*")]
    assert (star["n"], star["right"], star["wrong"], star["undetermined"]) == ("10", "8", "1", "1")
    assert star["rate"] == f"{8 / 9:.4f}"
    assert star["units"] == "5" and star["units_wrong"] == "1"
    wrong = gzip.open(tmp_path / "wrong.tsv.gz", "rt").read().split("\n")
    fields = wrong[1].split("\t")
    assert fields[7] == ""  # `broad` missing: an empty field under pandas 2 and 3 alike, never `nan`
    scores = list(csv.DictReader(open(tmp_path / "scores.tsv"), delimiter="\t"))
    assert scores[0]["variant"] == "*" and scores[0]["n"] == "10"


def test_limitations_rows_are_measured_never_typed(tmp_path, data_root):
    run = data_root / "families" / "x" / "runs" / "t" / "sus_scrofa"
    os.makedirs(run)
    json.dump(dict(organism="sus_scrofa", commit="d4c2dda25b61cf0b", graph_digests=["abc", "None"]),
              open(run / "run.json", "w"))
    with open(run / "report.tsv", "w") as fo:
        fo.write("stratum\tview\tn\nretired_none\t*\t200\nretired_none\tbest\t100\ncase_variant\t*\t50\n")
    with gzip.open(run / "causes.tsv.gz", "wt") as fo:
        fo.write("organism\tqid\tview\tkind\tstratum\tunit\tcause\tclass\tdetail\n")
        for i in range(4):
            fo.write(f"sus_scrofa\tq{i}\tbest\tretired_none\tretired_none\tu{i % 2}\texternal_bridge\tb\t\n")
        fo.write("sus_scrofa\tq9\tbest\tcase_variant\tcase_variant\tu9\tnot_a_known_cause\t?\t\n")
    rows = LIM.build("x", [str(run)], LIM.load_map(None))
    by = {r["cause"]: r for r in rows}
    b = by["external_bridge"]
    assert (b["answers"], b["units"], b["of_answers"], b["share"]) == (4, 2, 200, "0.02000")
    assert b["limitation"] == "L-bridges" and b["graph_digest"] == "abc"
    assert b["run"] == "families/x/runs/t/sus_scrofa"  # a place in the data root, never a machine path
    assert by["not_a_known_cause"]["title"] == "UNMAPPED"
    out = tmp_path / "lim.tsv"
    LIM.write(rows, str(out))
    head = open(out).readline().rstrip("\n").split("\t")
    assert head == LIM.HEADER and "rule" in head and "decision" not in head
