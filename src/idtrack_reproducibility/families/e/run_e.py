"""Family E's runner: F0's runner (`../runner.py`: sharded, resumable, Parquet answers) with family E's adapter.

Usage (as `../runner.py`; the tree's Python is chosen by `common.python_for`):
  run_e.py --organism <o> --clone <clone> --graph <graph file> --out <run folder>/<organism> [--sets a,b]
           [--shard i/n] [--limit N] [--timeout 30]
  run_e.py ... --submit <n>      n shard jobs, then exit (score with score_e.py once they are done)

Lines whose `generations` leave out the tree's api generation are not asked (counted in the log and in
`<out>/answers/skipped_<i>_of_<n>.json`).
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import runner as R

FAMILY = "e"


def run(a):
    i, n = (int(x) for x in a.shard.split("/"))
    mine, shas, excluded = L.load(FAMILY, a.organism, a.sets, a.strata, None, shard=(i, n), limit=a.limit)
    gen = C.generation(a.clone)
    from idtrack_reproducibility.families.e import adapter_e
    skipped = [q["qid"] for q in mine if not adapter_e.asked_on(q, gen)]
    mine = [q for q in mine if adapter_e.asked_on(q, gen)]
    os.makedirs(os.path.join(a.out, "answers"), exist_ok=True)
    json.dump(dict(generation=gen, n=len(skipped), qids=skipped[:2000]),
              open(os.path.join(a.out, "answers", f"skipped_{i}_of_{n}.json"), "w"))
    done = R.answered(a.out, (i, n))
    todo = [q for q in mine if q["qid"] not in done]
    sys.stderr.write(f"shard {i}/{n}: {len(mine)} lines ({len(skipped)} not asked on {gen}), {len(todo)} to ask\n")
    if not todo:
        return
    ad = adapter_e.AdapterE(a.clone, a.graph, a.organism, a.release, a.out, timeout=a.timeout)
    info = dict(ad.info, shard=a.shard, sets=shas, host=os.uname().nodename, slurm_job=os.environ.get("SLURM_JOB_ID"),
                started=time.strftime("%Y-%m-%d %H:%M:%S"), n_lines=len(mine), n_todo=len(todo), adapter="e/adapter_e.py")
    a.family = FAMILY
    todo, info["reused"] = R.reuse(a, todo, shas, info)
    info["n_asked"] = len(todo)
    json.dump(info, open(os.path.join(a.out, "answers", f"load_{i}_of_{n}.json"), "w"), indent=1)
    sys.stderr.write(f"shard {i}/{n}: reused {sum(info['reused'].values())} ({info['reused']}), {len(todo)} to ask\n")
    buf, t0, tf = [], time.time(), time.time()
    for k, q in enumerate(todo):
        buf.extend(ad.ask(q))
        if len(buf) >= a.flush or time.time() - tf > 300:
            R.write_part(a.out, i, n, buf)
            buf, tf = [], time.time()
        if k % 2000 == 0:
            sys.stderr.write(f"{k}/{len(todo)} {time.time() - t0:.0f}s\n")
    if buf:
        R.write_part(a.out, i, n, buf)
    sys.stderr.write(f"shard {i}/{n} done in {time.time() - t0:.0f}s\n")


def submit(a):
    logs = os.path.join(a.out, "logs")
    os.makedirs(logs, exist_ok=True)
    args = (f"--organism {a.organism} --clone {a.clone} --graph {a.graph} --out {a.out} --timeout {a.timeout} "
            f"--release {a.release} --flush {a.flush}" + (f" --sets {','.join(a.sets)}" if a.sets else "")
            + (f" --limit {a.limit}" if a.limit else ""))
    from idtrack_reproducibility import jobs
    mem = a.mem or (jobs.memory_for(a.organism) if C.generation(a.clone) == "d4c2dda" else 16)
    ids = [jobs.submit(f"e_{a.tag or a.organism}_{i}",
                       jobs.module(C.python_for(a.clone), "idtrack_reproducibility.families.e.run_e",
                                   f"{args} --shard {i}/{a.submit}"),
                       logs, cpus=2, mem_gb=mem, time=a.time) for i in range(a.submit)]
    print(" ".join(ids))
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sets", default=None)
    ap.add_argument("--strata", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--release", type=int, default=None)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--flush", type=int, default=20000)
    ap.add_argument("--submit", type=int, default=0)
    ap.add_argument("--mem", type=int, default=None)
    ap.add_argument("--time", default="04:00:00")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0" and not a.submit:
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
    a.out, a.clone, a.graph = os.path.abspath(a.out), os.path.abspath(a.clone), os.path.abspath(a.graph)
    a.sets = a.sets.split(",") if a.sets else None
    a.strata = a.strata.split(",") if a.strata else None
    a.release = a.release or C.release_of_graph(a.graph)
    os.makedirs(a.out, exist_ok=True)
    if a.submit:
        return submit(a)
    run(a)


if __name__ == "__main__":
    main()
