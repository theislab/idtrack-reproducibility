"""H, measure `calls`: per-call latency and batch throughput of one organism's conversions through a tree's api.

Usage (in the tree's Python; `run.py` starts it): python -m idtrack_reproducibility.families.h.calls --tree <tree>
    --organism <organism> --out <file.json> [--n 10000] [--threads 4] [--batch-sizes 1,10,100,1000,10000]

The identifiers: `--n` distinct identifier strings of the organism's corpus `convert_identifier` calls (every
stratum, seeded; the corpora hold about 7,500 to 8,100, so 10,000 is all of them). In one process, after the load:

- A: one call per identifier with the defaults (the graph's last release, `ensembl_gene`, the best answer), the
  per-node caches cold: wall time, latency percentiles, the 25 slowest, time in the cyclic collector, outcomes.
- B: the same identifiers as one batch (warm caches); its answers must equal A's (`B_differs_from_A`).
- C: the same batch to the organism's symbol database.
- D: the best and all answers with explained paths. The d4c2dda generation needs two batches (`best`, then `all`
  with `explain`); the typed api ranks every candidate in one call, so one batch with `explain` gives both views.
- E: the batch of B on `--threads` threads (the identifiers dealt round-robin), against B serially again; the answers
  must equal B's.
- F: every `convert_identifier` call of the corpus with its own arguments (release, database, strategy, explain):
  latency percentiles, overall and per stratum of the corpus, the slowest calls, resident memory every 2,000 calls.
- S: batch scaling: batches of each size of `--batch-sizes` (the first identifiers of the sample, warm caches): time
  per identifier and resident memory after each.

The bounded per-node caches of the graph are reported at the end (`cache_info`, where the tree exposes them).
"""

import argparse
import gc
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor

from idtrack_reproducibility.families.h import common as H


def outcome_of(tree, r):
    """The outcome flags of one answer, named as the d4c2dda generation names them."""
    if isinstance(r, dict) and "error" in r and len(r) == 1:
        return {}
    if tree.gen == "d4c2dda":
        return {k: bool(r.get(k)) for k in ("no_corresponding", "no_conversion", "no_target", "ambiguous", "tie")}
    best = r.view("best")
    return {"no_corresponding": r.status in ("not_found", "unknown_version"), "no_conversion": r.status == "no_conversion",
            "no_target": bool(r.no_target), "ambiguous": len({c.target for c in best}) > 1, "tie": bool(r.tie)}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.calls", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--batch-sizes", default="1,10,100,1000,10000")
    a = ap.parse_args(argv)
    org = a.organism
    report = {"measure": "calls", "tree": a.tree, "organism": org, **H.environment()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    link = tree.link(org, H.run_dir(a.tree, org))
    t = H.Timer()
    api = tree.load(link)
    report["load_s"] = t.lap()
    report["rss_after_load"] = H.statm()
    report["content_digest"] = tree.manifest(api).get("content_digest")
    ids = H.distinct_ids(org, a.n)
    report["n_ids"] = len(ids)
    gct = H.GCTimer()

    # A: one call per identifier, cold caches
    times, results_a, errors, outcomes = [], [], {}, {}
    gct.reset()
    t.lap()
    for i in ids:
        t0 = time.perf_counter()
        try:
            r = tree.convert(api, i)
        except Exception as error:  # noqa: BLE001 - recorded
            r = {"error": f"{type(error).__name__}: {error}"}
            errors[type(error).__name__] = errors.get(type(error).__name__, 0) + 1
        times.append(time.perf_counter() - t0)
        results_a.append(r)
    report["A_wall_s"] = t.lap()
    report["A_latency"] = H.latency(times)
    report["A_gc"] = gct.take()
    report["A_errors"] = errors
    report["A_slowest"] = sorted(([round(dt * 1e3, 1), i] for dt, i in zip(times, ids)), reverse=True)[:25]
    for r in results_a:
        for k, v in outcome_of(tree, r).items():
            outcomes[k] = outcomes.get(k, 0) + int(v)
    report["A_outcomes"] = outcomes
    report["rss_after_A"] = H.statm()

    # B: one batch, warm caches
    gct.reset()
    t.lap()
    results_b = tree.convert_many(api, ids)
    report["B_wall_s"] = t.lap()
    report["B_gc"] = gct.take()
    raised = {k for k, x in enumerate(results_a) if isinstance(x, dict) and set(x) == {"error"}}
    plain_a = [json.dumps(tree.plain(x), sort_keys=True, default=str) for x in results_a]
    plain_b = [json.dumps(tree.plain(x), sort_keys=True, default=str) for x in results_b]
    report["B_differs_from_A"] = sum(plain_a[k] != plain_b[k] for k in range(len(ids)) if k not in raised)

    # C: to the symbol database
    gct.reset()
    t.lap()
    tree.convert_many(api, ids, final_database=H.SYMBOL_DB[org])
    report["C_symbol_wall_s"] = t.lap()
    report["C_gc"] = gct.take()

    # D: the best and all answers with explained paths
    t.lap()
    if tree.gen == "d4c2dda":
        tree.convert_many(api, ids, strategy="best")
        report["D_best_s"] = t.lap()
        tree.convert_many(api, ids, strategy="all", explain=True)
        report["D_all_explain_s"] = t.lap()
        report["D_both_s"] = round(report["D_best_s"] + report["D_all_explain_s"], 4)
    else:
        tree.convert_many(api, ids, explain=True)
        report["D_both_s"] = t.lap()

    # E: threads
    chunks = [ids[k :: a.threads] for k in range(a.threads)]
    t.lap()
    with ThreadPoolExecutor(a.threads) as pool:
        parts = list(pool.map(lambda chunk: tree.convert_many(api, chunk), chunks))
    report["E_threads"] = a.threads
    report["E_threads_wall_s"] = t.lap()
    by_id = {}
    for chunk, part in zip(chunks, parts):
        by_id.update({i: json.dumps(tree.plain(r), sort_keys=True, default=str) for i, r in zip(chunk, part)})
    report["E_equals_B"] = [by_id[i] for i in ids] == plain_b
    t.lap()
    tree.convert_many(api, ids)
    report["E_serial_again_s"] = t.lap()

    # F: every corpus call with its own arguments
    calls = [c for c in H.corpus_calls(org) if isinstance(c["args"].get("identifier"), str)]
    ftimes, slow, rss_trace, ferr = [], [], [], {}
    gct.reset()
    for k, c in enumerate(calls):
        args = dict(c["args"])
        identifier = args.pop("identifier")
        t0 = time.perf_counter()
        try:
            tree.convert(api, identifier, **args)
        except Exception as error:  # noqa: BLE001 - recorded
            ferr[type(error).__name__] = ferr.get(type(error).__name__, 0) + 1
        dt = time.perf_counter() - t0
        ftimes.append(dt)
        slow.append((round(dt * 1e3, 1), c["qid"], c["stratum"]))
        if k % 2000 == 0:
            rss_trace.append([k, H.rss_gb()])
    report["F_latency"] = H.latency(ftimes)
    by_stratum = {}
    for (_ms, _qid, stratum), dt in zip(slow, ftimes):
        by_stratum.setdefault(stratum, []).append(dt)
    report["F_by_stratum"] = {s: H.latency(v) for s, v in sorted(by_stratum.items())}
    report["F_gc"] = gct.take()
    report["F_errors"] = ferr
    report["F_slowest"] = sorted(slow, key=lambda x: -x[0])[:30]
    report["F_rss_trace"] = rss_trace

    # S: batch scaling (warm caches)
    scaling = []
    for size in (int(x) for x in a.batch_sizes.split(",") if x):
        batch = ids[:size]
        if len(batch) < size:
            continue
        gc.collect()
        t.lap()
        tree.convert_many(api, batch)
        wall = t.lap()
        scaling.append({"size": size, "wall_s": wall, "ms_per_id": round(wall / size * 1e3, 4), "rss_GB": H.rss_gb()})
    report["S_scaling"] = scaling

    report["cache_info"] = tree.cache_info(api)
    report["rss_end"] = H.statm()
    report["peak_GB"] = H.peak_gb()
    H.dump(a.out, report)
    print(json.dumps({k: v for k, v in report.items() if not k.endswith("slowest")}, default=str)[:4000])


if __name__ == "__main__":
    main()
