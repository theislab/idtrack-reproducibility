"""H, measure `server`: what a small server of the five organisms needs: start time, memory, and the memory that forked
worker processes share.

Usage (in the tree's Python; `run.py` starts it with the organism `all`): python -m
    idtrack_reproducibility.families.h.server --tree <tree> --out <file.json> [--workers 4] [--per-worker 4000]
    [--no-freeze]

1. Load the five graphs in one process, each timed, resident memory after each.
2. `gc.collect()`, then `gc.freeze()` (unless `--no-freeze`), then fork the workers. Each worker converts a seeded
   sample of the corpus's `convert_identifier` calls of all five organisms with their own arguments (as requests
   would come) and reports its wall time, latency percentiles and its memory from `/proc/self/smaps_rollup` (Rss,
   Pss, private and shared pages) before and after. When every worker is done and still alive, the memory of all of
   them and of the parent is read at one moment, so that shared pages are divided among all: the Pss summed over the
   parent and the workers is the memory the server really uses.
"""

import argparse
import gc
import json
import os
import random
import signal
import time

from idtrack_reproducibility.families.h import common as H


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.server", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--per-worker", type=int, default=4000)
    ap.add_argument("--no-freeze", action="store_true")
    a = ap.parse_args(argv)
    report = {"measure": "server", "tree": a.tree, "organism": "all", "workers": a.workers, "per_worker": a.per_worker,
              "freeze": not a.no_freeze, **H.environment()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    links = {org: tree.link(org, H.run_dir(a.tree, org)) for org in H.ORGANISMS}
    t0 = time.perf_counter()
    import idtrack  # noqa: F401

    report["import_s"] = round(time.perf_counter() - t0, 3)
    apis, report["load"] = {}, {}
    for org in H.ORGANISMS:
        t = time.perf_counter()
        apis[org] = tree.load(links[org])
        report["load"][org] = {"s": round(time.perf_counter() - t, 3), "rss_after_GB": H.rss_gb(),
                               "content_digest": tree.manifest(apis[org]).get("content_digest")}
    report["start_s"] = round(time.perf_counter() - t0, 3)
    report["parent_after_load"] = H.smaps_rollup()
    calls = []
    for org in H.ORGANISMS:
        calls += [(org, c["args"]) for c in H.corpus_calls(org) if isinstance(c["args"].get("identifier"), str)]
    rng = random.Random(H.SEED_IDS)
    gc.collect()
    if not a.no_freeze:
        gc.freeze()
    report["gc_tracked_parent"] = len(gc.get_objects())
    children = {}
    for _w in range(a.workers):
        sample = rng.sample(calls, a.per_worker)
        read_fd, write_fd = os.pipe()
        pid = os.fork()
        if pid == 0:  # worker
            os.close(read_fd)
            before = H.smaps_rollup()
            times, errors = [], 0
            start = time.perf_counter()
            for org, args in sample:
                args = dict(args)
                identifier = args.pop("identifier")
                t = time.perf_counter()
                try:
                    tree.convert(apis[org], identifier, **args)
                except Exception:  # noqa: BLE001 - counted
                    errors += 1
                times.append(time.perf_counter() - t)
            wall = time.perf_counter() - start
            out = {"wall_s": round(wall, 2), "errors": errors, **H.latency(times), "before": before,
                   "after": H.smaps_rollup()}
            with os.fdopen(write_fd, "w") as fh:
                json.dump(out, fh)
            time.sleep(900)  # stay alive until the parent has read every worker's memory at one moment
            os._exit(0)
        os.close(write_fd)
        children[pid] = read_fd
    report["workers_report"] = []
    for _pid, fd in children.items():
        with os.fdopen(fd) as fh:
            report["workers_report"].append(json.loads(fh.read()))
    report["parent_end"] = H.smaps_rollup()
    report["workers_memory_together"] = [H.smaps_rollup(pid) for pid in children]
    for pid in children:
        os.kill(pid, signal.SIGKILL)
        os.waitpid(pid, 0)
    report["pss_total_GB"] = round(report["parent_end"]["Pss"] + sum(m["Pss"] for m in report["workers_memory_together"]), 3)
    report["private_per_worker_GB"] = [round(m["Private_Clean"] + m["Private_Dirty"], 3)
                                       for m in report["workers_memory_together"]]
    report["rss_per_worker_GB"] = [m["Rss"] for m in report["workers_memory_together"]]
    H.dump(a.out, report)
    print(json.dumps(report, default=str)[:4000])


if __name__ == "__main__":
    main()
