"""H: the result tables, generated from the stored measurements (docs/families/h.md, "Results").

    python -m idtrack_reproducibility.families.h.tables [--results <repository>/results]

Reads every measurement `<performance>/runs/<tree>/<organism>/<measure>_rep<k>.json` of a configured tree and writes
`results/h/<tree>/<organism>/<table>.tsv`; the server's tables go to `results/h/<tree>/all/`. Every table has one
row per repetition and, per quantity, a row `median` (the lower median of the repetitions, so every median is a
measured value and nothing is rounded anew). The tables:

| table | from | rows |
| --- | --- | --- |
| `load.tsv` | `load` | import and load time, resident and peak memory after the load, nodes, edges, the graph file, the collector; for the d4c2dda generation also the store's memory by part |
| `calls.tsv` | `calls` | per measurement (`one_call_cold`, `batch_warm`, `batch_symbols`, `best_and_all_explained`, `threads_<n>`, `serial_again`, `corpus_calls`): wall time, latency percentiles, time in the collector |
| `calls_by_stratum.tsv` | `calls` | the corpus calls' latency per stratum of the corpus |
| `batch_scaling.tsv` | `calls` | time per identifier and resident memory by batch size |
| `harmonize.tsv` | `harmonize` | features in and out, constructor, unify and constructor-from-cache times with the collector's share, peak memory, a digest of the outputs |
| `cycles.tsv` | `cycles` | per kind of call: conversions that leave reference cycles, unreachable objects, memory they held |
| `build.tsv` | `build_<kind>` | build time, peak memory, bytes received, sizes, whether the graph is the configured one |
| `server.tsv` | `server` | start time, the memory of the parent and the workers, the workers' latency |

Memory is in GB of 10^9 bytes (a measurement that does not state this unit is refused). The tables name the processor
of each repetition, never the machine.

Provenance, as `idtrack_reproducibility.results` reads it for tables a command writes straight into `results/`:
beside each table, `<table>.tsv.provenance` holds one row per measurement the table rests on (the command, the
measurement's path in the data root with its sha256, when it was measured, the tree's commit, the graph's content
digest, the time of generation); `idtrack-repro results` (registry entry `[[generated]]` of family `h`) folds them into
`results/provenance.tsv`. Every table and provenance file this command wrote before is removed first, so a
measurement removed from the data removes its table.
"""

import argparse
import csv
import glob
import hashlib
import json
import os
import re
import statistics
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families.h import common as H

COMMAND = "python -m idtrack_reproducibility.families.h.tables"
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
FILE_RE = re.compile(r"^(?P<measure>[a-z]+(?:_[a-z]+)?)_rep(?P<rep>\d+)\.json$")
PROVENANCE_COLUMNS = ["file", "sha256", "generated_by", "source", "source_sha256", "scored_at", "commit", "graph_digest",
                      "generated_at"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def median_row(rows, key_cols):
    """One `median` row per value of the key columns: the lower median of every numeric column over the repetitions."""
    out = []
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r.get(k) for k in key_cols), []).append(r)
    for key, rs in groups.items():
        if len(rs) < 2:
            continue
        m = {"rep": "median", **dict(zip(key_cols, key))}
        for col in rs[0]:
            if col in m or col == "rep":
                continue
            vals = [r.get(col) for r in rs if isinstance(r.get(col), (int, float)) and not isinstance(r.get(col), bool)]
            m[col] = statistics.median_low(vals) if len(vals) == len(rs) else ""
        out.append(m)
    return out


def gb(statm):
    return None if statm is None else statm["rss_GB"]


def load_rows(ms):
    rows = []
    for rep, d in ms:
        rows.append({
            "rep": rep, "cpu_model": d.get("cpu_model"), "python": d["python"],
            "import_idtrack_s": d["import_idtrack_s"], "import_api_s": d["import_api_s"], "load_s": d["load_s"],
            "rss_after_load_gb": gb(d["rss_after_load"]), "peak_after_load_gb": d["peak_after_load_GB"],
            "nodes": d["nodes"], "edges": d["edges"], "graph_file_mb": d["graph_file_MB"],
            "files_written_beside_graph": len(d["files_written_beside_graph"]),
            "gc_tracked_objects": d["gc_tracked_objects"], "gc_collect_s": d["gc_collect_s"],
            "modules_after_load": ",".join(d["modules_after_load"]),
            "store_arrays_mb": d.get("store_arrays_total_MB"), "names_list_mb": d.get("names_list_MB"),
            "index_dict_mb": d.get("index_dict_MB"), "index_int_objects_mb": d.get("index_int_objects_MB"),
            "ntype_list_mb": d.get("ntype_list_MB"), "connections_mb": d.get("connections_MB"),
            "versions_of_first_call_s": d.get("versions_of_first_call_s"),
            "versions_of_second_call_s": d.get("versions_of_second_call_s"),
        })
    return rows + median_row(rows, [])


def calls_rows(ms):
    rows = []
    for rep, d in ms:
        def add(name, wall=None, lat=None, gc_=None, errors=None):
            lat = lat or {}
            rows.append({"rep": rep, "measurement": name, "n": lat.get("n", d["n_ids"]), "wall_s": wall,
                         "total_s": lat.get("total_s"), "p50_ms": lat.get("p50_ms"), "p90_ms": lat.get("p90_ms"),
                         "p99_ms": lat.get("p99_ms"), "p999_ms": lat.get("p999_ms"), "max_ms": lat.get("max_ms"),
                         "mean_ms": lat.get("mean_ms"), "gc_s": (gc_ or {}).get("s"),
                         "errors": None if errors is None else sum(errors.values())})
        add("one_call_cold", d["A_wall_s"], d["A_latency"], d["A_gc"], d["A_errors"])
        add("batch_warm", d["B_wall_s"], None, d["B_gc"])
        add("batch_symbols", d["C_symbol_wall_s"], None, d["C_gc"])
        add("best_and_all_explained", d["D_both_s"])
        add(f"threads_{d['E_threads']}", d["E_threads_wall_s"])
        add("serial_again", d["E_serial_again_s"])
        add("corpus_calls", d["F_latency"]["total_s"], d["F_latency"], d["F_gc"], d["F_errors"])
    return rows + median_row(rows, ["measurement"])


def stratum_rows(ms):
    rows = []
    for rep, d in ms:
        for s, lat in d.get("F_by_stratum", {}).items():
            rows.append({"rep": rep, "stratum": s, "n": lat["n"], "total_s": lat["total_s"], "p50_ms": lat["p50_ms"],
                         "p90_ms": lat["p90_ms"], "p99_ms": lat["p99_ms"], "max_ms": lat["max_ms"],
                         "mean_ms": lat["mean_ms"]})
    return rows + median_row(rows, ["stratum"])


def scaling_rows(ms):
    rows = [{"rep": rep, "size": s["size"], "wall_s": s["wall_s"], "ms_per_id": s["ms_per_id"], "rss_gb": s["rss_GB"]}
            for rep, d in ms for s in d.get("S_scaling", [])]
    return rows + median_row(rows, ["size"])


def harmonize_rows(ms):
    rows = []
    for rep, d in ms:
        dg = d["digests"]
        rows.append({
            "rep": rep, "cpu_model": d.get("cpu_model"), "freeze": d["freeze"],
            "features_in_d115": d["features_in"]["d115"], "features_in_dold": d["features_in"]["dold"],
            "features_in_s115": d["features_in"]["s115"], "features_out": d["features_out"],
            "constructor_s": d["constructor_s"], "constructor_gc_s": d["constructor_gc"]["s"],
            "unify_s": d["unify_s"], "unify_gc_s": d["unify_gc"]["s"],
            "from_cache_s": d["constructor_from_cache_s"], "from_cache_gc_s": d["from_cache_gc"]["s"],
            "rss_before_run_gb": gb(d["rss_before_run"]), "peak_gb": d["peak_GB"],
            "result_cache_mb": round(sum(d["result_cache_MB"].values()), 2),
            "output_digest": hashlib.sha256(json.dumps({k: dg[k] for k in sorted(dg) if k != "cache_files"},
                                                       sort_keys=True).encode()).hexdigest()[:16],
        })
    return rows + median_row(rows, ["freeze"])


def cycles_rows(ms):
    rows = []
    for rep, d in ms:
        for kind in ("default", "symbol", "all_explain"):
            k = d[kind]
            rows.append({"rep": rep, "kind": kind, "n": d["n"], "conversions_leaving_cycles": k["conversions_leaving_cycles"],
                         "unreachable_objects_total": k["unreachable_objects_total"],
                         "per_conversion_median": k["per_conversion_median"], "per_conversion_max": k["per_conversion_max"],
                         "freed_total_mb": k["bytes_freed_by_collection_total_MB"], "freed_max_mb": k["bytes_freed_max_MB"]})
    return rows + median_row(rows, ["kind"])


def build_rows(ms):
    rows = []
    for rep, d in ms:
        rows.append({"rep": rep, "kind": d["kind"], "profiled": d["profiled"], "cpu_model": d.get("cpu_model"), "build_s": d["build_s"],
                     "peak_gb": d["peak_GB"], "bytes_received_mb": round(d["bytes_received"] / 1e6, 1),
                     "repo_mb": round(d["repo_bytes"] / 1e6, 1),
                     "graph_file_mb": None if d["graph_file_bytes"] is None else round(d["graph_file_bytes"] / 1e6, 1),
                     "same_content_as_configured": d["same_content_as_configured"]})
    return rows + median_row(rows, ["kind", "profiled"])


def server_rows(ms):
    rows = []
    for rep, d in ms:
        w, mem = d["workers_report"], d["workers_memory_together"]
        row = {"rep": rep, "cpu_model": d.get("cpu_model"), "freeze": d["freeze"], "workers": d["workers"], "per_worker": d["per_worker"],
               "import_s": d["import_s"], "start_s": d["start_s"]}
        for org in H.ORGANISMS:
            row[f"load_s_{org}"] = d["load"][org]["s"]
        row.update({
            "rss_after_load_gb": d["load"][H.ORGANISMS[-1]]["rss_after_GB"],
            "parent_pss_gb": d["parent_end"]["Pss"], "pss_total_gb": d["pss_total_GB"],
            "worker_private_gb_min": min(d["private_per_worker_GB"]), "worker_private_gb_max": max(d["private_per_worker_GB"]),
            "worker_pss_gb_min": min(m["Pss"] for m in mem), "worker_pss_gb_max": max(m["Pss"] for m in mem),
            "worker_wall_s_min": min(x["wall_s"] for x in w), "worker_wall_s_max": max(x["wall_s"] for x in w),
            "worker_p50_ms_min": min(x["p50_ms"] for x in w), "worker_p50_ms_max": max(x["p50_ms"] for x in w),
            "worker_p99_ms_min": min(x["p99_ms"] for x in w), "worker_p99_ms_max": max(x["p99_ms"] for x in w),
            "worker_max_ms": max(x["max_ms"] for x in w), "errors": sum(x["errors"] for x in w),
        })
        rows.append(row)
    return rows + median_row(rows, ["freeze"])


# table name -> (measure prefix, row builder)
TABLES = {
    "load": ("load", load_rows),
    "calls": ("calls", calls_rows),
    "calls_by_stratum": ("calls", stratum_rows),
    "batch_scaling": ("calls", scaling_rows),
    "harmonize": ("harmonize", harmonize_rows),
    "cycles": ("cycles", cycles_rows),
    "build": ("build_", build_rows),
    "server": ("server", server_rows),
}


def measurements(folder):
    """measure -> [(rep, data, path)], sorted by rep."""
    out = {}
    for path in sorted(glob.glob(os.path.join(folder, "*_rep*.json"))):
        m = FILE_RE.match(os.path.basename(path))
        if not m:
            continue
        with open(path) as fh:
            d = json.load(fh)
        if d.get("memory_unit") != H.MEMORY_UNIT:
            raise SystemExit(f"{path}: no memory_unit {H.MEMORY_UNIT!r}; not a measurement of this package")
        out.setdefault(m.group("measure"), []).append((int(m.group("rep")), d, path))
    for v in out.values():
        v.sort(key=lambda x: x[0])
    return out


def write_tsv(path, rows):
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow(["" if r.get(c) is None else str(r.get(c)).lower() if isinstance(r.get(c), bool) else r.get(c)
                        for c in cols])


def generate(results_dir, trees=None):
    """Write every table of every configured tree (or of `trees`) and its provenance; returns the tables written."""
    cfg = CFG.get()
    runs = os.path.join(H.root(), "runs")
    hdir = os.path.join(results_dir, "h")
    for old in glob.glob(os.path.join(hdir, "*", "*", "*.tsv")) + glob.glob(os.path.join(hdir, "*", "*", "*.tsv.provenance")):
        os.remove(old)
    written = []
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    for tree in sorted(trees or cfg.trees):
        if not os.path.isdir(os.path.join(runs, tree)):
            continue
        for organism in sorted(os.listdir(os.path.join(runs, tree))):
            folder = os.path.join(runs, tree, organism)
            if not os.path.isdir(folder):
                continue
            ms = measurements(folder)
            for table, (prefix, builder) in TABLES.items():
                picked = [(rep, d, p) for name, v in sorted(ms.items()) if name == prefix or
                          (prefix.endswith("_") and name.startswith(prefix)) for rep, d, p in v]
                if not picked:
                    continue
                rows = builder([(rep, d) for rep, d, _ in picked])
                if not rows:
                    continue
                out = os.path.join(hdir, tree, organism, f"{table}.tsv")
                write_tsv(out, rows)
                rel = os.path.relpath(out, os.path.dirname(os.path.abspath(results_dir)))
                prov = []
                for _, d, p in picked:
                    digests = sorted(({d.get("content_digest")} | {x.get("content_digest") for x in (d.get("load") or {}).values()
                                                                   if isinstance(x, dict)}) - {None})
                    prov.append({"file": rel, "sha256": sha256(out), "generated_by": COMMAND,
                                 "source": os.path.relpath(p, cfg.root), "source_sha256": sha256(p),
                                 "scored_at": d.get("date", ""), "commit": d.get("commit", ""),
                                 "graph_digest": ";".join(digests), "generated_at": stamp})
                with open(out + ".provenance", "w", newline="") as fh:
                    w = csv.DictWriter(fh, PROVENANCE_COLUMNS, delimiter="\t", lineterminator="\n")
                    w.writeheader()
                    w.writerows(prov)
                written.append(rel)
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.tables", description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", default=os.path.join(REPO_ROOT, "results"), help="the repository's results folder")
    ap.add_argument("--trees", default=None, help="comma-separated (default: every configured tree)")
    a = ap.parse_args(argv)
    for rel in generate(a.results, a.trees.split(",") if a.trees else None):
        print(rel)


if __name__ == "__main__":
    main()
