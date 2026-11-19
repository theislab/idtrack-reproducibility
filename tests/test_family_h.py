"""Family H (performance): the latency convention, the interleaved order, the api calls of both generations, and the
result tables with their provenance, on invented measurements (no data, no idtrack)."""
import csv
import hashlib
import json
import os

import pytest

from idtrack_reproducibility.families.h import calls as HC
from idtrack_reproducibility.families.h import common as H
from idtrack_reproducibility.families.h import run as HR
from idtrack_reproducibility.families.h import tables as T


def test_percentiles_are_the_element_at_n_times_q_of_the_sorted_times():
    times = [i / 1000 for i in range(100, 0, -1)]  # 0.1 s down to 1 ms, unsorted
    lat = H.latency(times)
    assert lat["n"] == 100 and lat["p50_ms"] == 51.0 and lat["p90_ms"] == 91.0 and lat["p99_ms"] == 100.0
    assert lat["max_ms"] == 100.0 and lat["total_s"] == 5.05 and lat["mean_ms"] == 50.5
    assert H.pct([5], 0.99) == 5


def test_trees_are_interleaved_with_the_order_reversed_every_other_repetition():
    assert HR.order(["d4c2dda", "71690ab"], 3) == [("d4c2dda", 1), ("71690ab", 1), ("71690ab", 2), ("d4c2dda", 2),
                                                   ("d4c2dda", 3), ("71690ab", 3)]
    assert HR.order(["a"], 2) == [("a", 1), ("a", 2)]


class FakeTypedApi:
    def __init__(self):
        self.calls = []

    def convert(self, identifier, **kw):
        self.calls.append(("convert", identifier, kw))
        return "r"

    def convert_many(self, identifiers, **kw):
        self.calls.append(("convert_many", list(identifiers), kw))
        return ["r"]


def typed_tree():
    t = object.__new__(H.Tree)
    t.name, t.gen, t.cfg = "71690ab", "convert", {"commit": "c"}
    return t


def test_the_typed_api_is_asked_without_strategy_and_with_its_own_argument_names():
    t, api = typed_tree(), FakeTypedApi()
    t.convert(api, "X", strategy="all", go_external=False, to_release=100, explain=True)
    t.convert_many(api, ["X"], strategy="best", final_database="HGNC Symbol")
    assert api.calls == [("convert", "X", {"allow_bridges": False, "to_release": 100, "explain": True}),
                         ("convert_many", ["X"], {"final_database": "HGNC Symbol"})]


def test_typed_outcomes_are_named_as_d4c2dda_names_them():
    class Cand:
        def __init__(self, target):
            self.target = target

    class Result:
        status, no_target, tie = "converted", False, True

        def view(self, strategy):
            assert strategy == "best"
            return [Cand("A"), Cand("B")]

    out = HC.outcome_of(typed_tree(), Result())
    assert out == {"no_corresponding": False, "no_conversion": False, "no_target": False, "ambiguous": True, "tie": True}
    Result.status = "unknown_version"
    assert HC.outcome_of(typed_tree(), Result())["no_corresponding"] is True
    assert HC.outcome_of(typed_tree(), {"error": "ValueError: x"}) == {}


def test_median_rows_are_the_lower_median_of_the_repetitions_per_key():
    rows = [{"rep": 1, "m": "a", "x": 3, "s": "t"}, {"rep": 2, "m": "a", "x": 1, "s": "t"}, {"rep": 3, "m": "a", "x": 2, "s": "u"},
            {"rep": 1, "m": "b", "x": 1.5}, {"rep": 2, "m": "b", "x": 2.5}, {"rep": 1, "m": "c", "x": 9}]
    med = {r["m"]: r for r in T.median_row(rows, ["m"])}
    assert med["a"]["x"] == 2 and med["a"]["rep"] == "median"
    assert med["b"]["x"] == 1.5  # the lower median: a measured value, never a mean of two
    assert med["a"]["s"] == ""  # not a number: no median
    assert "c" not in med  # one repetition: no median row
    assert T.median_row([{"rep": 1, "f": True}, {"rep": 2, "f": False}], [])[0]["f"] == ""  # flags are no numbers


def env():
    return {"date": "2026-11-19 10:00:00", "memory_unit": H.MEMORY_UNIT, "cpu_model": "cpu", "python": "3.13",
            "commit": "c0ffee", "content_digest": "dig"}


def load_measurement(rep):
    return {**env(), "import_idtrack_s": 0.03, "import_api_s": 1.0 + rep, "load_s": 2.0 * rep,
            "rss_after_load": {"rss_GB": 0.38, "shared_GB": 0.02}, "peak_after_load_GB": 0.4, "nodes": 10, "edges": 20,
            "graph_file_MB": 85.2, "files_written_beside_graph": [], "gc_tracked_objects": 5, "gc_collect_s": 0.1,
            "modules_after_load": ["numpy"]}


def calls_measurement(rep):
    lat = {"n": 3, "total_s": 0.3, "p50_ms": 0.2, "p90_ms": 0.3, "p99_ms": 1.0, "p999_ms": 1.0, "max_ms": 1.0, "mean_ms": 0.1}
    gc_ = {"s": 0.01, "n": [1, 0, 0], "s_by_gen": [0.01, 0, 0]}
    return {**env(), "n_ids": 3, "A_wall_s": 0.3 * rep, "A_latency": lat, "A_gc": gc_, "A_errors": {},
            "B_wall_s": 0.1, "B_gc": gc_, "C_symbol_wall_s": 0.2, "C_gc": gc_, "D_both_s": 0.4, "E_threads": 4,
            "E_threads_wall_s": 0.5, "E_serial_again_s": 0.1, "F_latency": dict(lat, n=19, total_s=1.9), "F_gc": gc_,
            "F_errors": {"ValueError": 2, "KeyError": 1}, "F_by_stratum": {"gene_split": lat},
            "S_scaling": [{"size": 1, "wall_s": 0.001, "ms_per_id": 1.0, "rss_GB": 0.4}]}


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_tables_are_generated_from_the_measurements_with_one_provenance_row_per_source(data_root, tmp_path):
    runs = os.path.join(str(data_root), "performance", "runs", "d4c2dda", "sus_scrofa")
    for rep in (1, 2, 3):
        write(os.path.join(runs, f"load_rep{rep}.json"), load_measurement(rep))
        write(os.path.join(runs, f"calls_rep{rep}.json"), calls_measurement(rep))
    write(os.path.join(runs, "build_warm_rep1.json"), {
        **env(), "kind": "warm", "profiled": False, "build_s": 130.5, "peak_GB": 1.2, "bytes_received": 0,
        "repo_bytes": 640_000_000, "graph_file_bytes": 85_229_259, "same_content_as_configured": True})
    write(os.path.join(runs, "notes.json"), {"not": "a measurement"})
    results = tmp_path / "repo" / "results"
    stale = results / "h" / "d4c2dda" / "sus_scrofa" / "server.tsv"
    os.makedirs(stale.parent)
    stale.write_text("old\n")
    written = T.generate(str(results), ["d4c2dda", "71690ab"])
    assert sorted(written) == sorted(f"results/h/d4c2dda/sus_scrofa/{t}.tsv" for t in
                                     ("load", "calls", "calls_by_stratum", "batch_scaling", "build"))
    assert not stale.exists()  # a table whose measurement is gone is removed
    load = read_tsv(results / "h" / "d4c2dda" / "sus_scrofa" / "load.tsv")
    assert [r["rep"] for r in load] == ["1", "2", "3", "median"]
    assert load[-1]["load_s"] == "4.0" and load[-1]["import_api_s"] == "3.0"
    assert "host" not in load[0]  # the processor, never the machine
    calls = read_tsv(results / "h" / "d4c2dda" / "sus_scrofa" / "calls.tsv")
    corpus = [r for r in calls if r["measurement"] == "corpus_calls" and r["rep"] == "1"][0]
    assert corpus["n"] == "19" and corpus["errors"] == "3" and corpus["wall_s"] == "1.9"
    assert {r["measurement"] for r in calls} == {"one_call_cold", "batch_warm", "batch_symbols", "best_and_all_explained",
                                                 "threads_4", "serial_again", "corpus_calls"}
    build = read_tsv(results / "h" / "d4c2dda" / "sus_scrofa" / "build.tsv")
    assert build == [{"rep": "1", "kind": "warm", "profiled": "false", "cpu_model": "cpu", "build_s": "130.5",
                      "peak_gb": "1.2", "bytes_received_mb": "0.0", "repo_mb": "640.0", "graph_file_mb": "85.2",
                      "same_content_as_configured": "true"}]
    prov = read_tsv(str(results / "h" / "d4c2dda" / "sus_scrofa" / "load.tsv") + ".provenance")
    assert [r["source"] for r in prov] == [f"performance/runs/d4c2dda/sus_scrofa/load_rep{k}.json" for k in (1, 2, 3)]
    with open(os.path.join(runs, "load_rep2.json"), "rb") as fh:
        assert prov[1]["source_sha256"] == hashlib.sha256(fh.read()).hexdigest()
    with open(results / "h" / "d4c2dda" / "sus_scrofa" / "load.tsv", "rb") as fh:
        assert prov[0]["sha256"] == hashlib.sha256(fh.read()).hexdigest()
    assert {r["commit"] for r in prov} == {"c0ffee"} and {r["graph_digest"] for r in prov} == {"dig"}
    assert prov[0]["generated_by"] == T.COMMAND and prov[0]["scored_at"] == "2026-11-19 10:00:00"


def test_a_measurement_without_its_memory_unit_is_refused(data_root, tmp_path):
    runs = os.path.join(str(data_root), "performance", "runs", "d4c2dda", "sus_scrofa")
    m = load_measurement(1)
    del m["memory_unit"]  # e.g. ru_maxrss read as kB: 2.4 % off; never mixed with GB of 10^9 bytes
    write(os.path.join(runs, "load_rep1.json"), m)
    with pytest.raises(SystemExit):
        T.generate(str(tmp_path / "results"), ["d4c2dda"])


def test_measurement_file_names():
    assert T.FILE_RE.match("build_warm_rep2.json").group("measure") == "build_warm"
    assert T.FILE_RE.match("load_rep10.json").group("rep") == "10"
    assert T.FILE_RE.match("load.json") is None and T.FILE_RE.match("load_rep1.json.tmp") is None


def test_server_rows_take_the_extremes_over_the_workers():
    mem = lambda pss, priv: {"Rss": 1.0, "Pss": pss, "Private_Clean": 0.0, "Private_Dirty": priv}  # noqa: E731
    d = {**env(), "freeze": True, "workers": 2, "per_worker": 10, "import_s": 0.1, "start_s": 12.0,
         "load": {o: {"s": 1.0, "rss_after_GB": 2.4} for o in H.ORGANISMS}, "parent_end": {"Pss": 1.1},
         "pss_total_GB": 3.3, "private_per_worker_GB": [0.8, 0.9], "workers_memory_together": [mem(1.1, 0.8), mem(1.2, 0.9)],
         "workers_report": [{"wall_s": 10.0, "p50_ms": 0.3, "p99_ms": 40.0, "max_ms": 900.0, "errors": 1},
                            {"wall_s": 12.0, "p50_ms": 0.4, "p99_ms": 50.0, "max_ms": 500.0, "errors": 2}]}
    row = T.server_rows([(1, d)])[0]
    assert (row["worker_private_gb_min"], row["worker_private_gb_max"]) == (0.8, 0.9)
    assert (row["worker_pss_gb_min"], row["worker_pss_gb_max"]) == (1.1, 1.2)
    assert row["worker_max_ms"] == 900.0 and row["errors"] == 3 and row["load_s_homo_sapiens"] == 1.0
