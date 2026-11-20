"""Family D's runner for the lines that are not one `convert` call (D2 `infer_source`, D4 `convert_many`): F0's runner
(`../runner.py`: sharded, resumable, Parquet answers) with family E's adapter (`../e/adapter_e.py`, which asks every
call of the API on both generations; imported, never copied). D1's lines are plain `convert` lines: F0's runner.

Usage:
  run_d.py --organism <o> --clone <clone> --graph <graph file> --out <run folder>/<organism> --sets d2_v1 [--shard i/n]
  run_d.py ... --submit <n> [--tag <tag>]     n shard jobs (the configured runner), named <name_prefix>d_<tag>_s<i>of<n>
Score with F0's scorer (`../scorer.py --family d ...`). A resubmission asks only what has no answer.
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import runner as R

FAMILY = "d"


def run(a):
    from idtrack_reproducibility.families.e import adapter_e
    i, n = (int(x) for x in a.shard.split("/"))
    mine, shas, _ = L.load(FAMILY, a.organism, a.sets, a.strata, None, shard=(i, n), limit=a.limit)
    os.makedirs(os.path.join(a.out, "answers"), exist_ok=True)
    done = R.answered(a.out, (i, n))
    todo = [q for q in mine if q["qid"] not in done]
    sys.stderr.write(f"shard {i}/{n}: {len(mine)} lines, {len(todo)} to ask\n")
    if not todo:
        return
    ad = adapter_e.AdapterE(a.clone, a.graph, a.organism, a.release, a.out, timeout=a.timeout)
    info = dict(ad.info, shard=a.shard, sets=shas, host=os.uname().nodename, slurm_job=os.environ.get("SLURM_JOB_ID"),
                started=time.strftime("%Y-%m-%d %H:%M:%S"), n_lines=len(mine), n_todo=len(todo), adapter="e/adapter_e.py")
    json.dump(info, open(os.path.join(a.out, "answers", f"load_{i}_of_{n}.json"), "w"), indent=1)
    buf, t0, tf = [], time.time(), time.time()
    for k, q in enumerate(todo):
        buf.extend(ad.ask(q))
        if len(buf) >= a.flush or time.time() - tf > 300:
            R.write_part(a.out, i, n, buf)
            buf, tf = [], time.time()
        if k % 200 == 0:
            sys.stderr.write(f"{k}/{len(todo)} {time.time() - t0:.0f}s\n")
    if buf:
        R.write_part(a.out, i, n, buf)
    sys.stderr.write(f"shard {i}/{n} done in {time.time() - t0:.0f}s\n")


def submit(a):
    logs = os.path.join(a.out, "logs")
    os.makedirs(logs, exist_ok=True)
    args = (f"--organism {a.organism} --clone {a.clone} --graph {a.graph} --out {a.out} --timeout {a.timeout} "
            f"--release {a.release} --flush {a.flush}" + (f" --sets {','.join(a.sets)}" if a.sets else "")
            + (f" --strata {','.join(a.strata)}" if a.strata else "") + (f" --limit {a.limit}" if a.limit else ""))
    from idtrack_reproducibility import jobs
    mem = a.mem or (jobs.memory_for(a.organism) if C.generation(a.clone) == "d4c2dda" else 16)
    ids = [jobs.submit(f"d_{a.tag}_s{i}of{a.submit}", jobs.module(C.python_for(a.clone), "idtrack_reproducibility.families.d.run_d",
                                                                 f"{args} --shard {i}/{a.submit}"),
                       logs, cpus=2, mem_gb=mem, time=a.time) for i in range(a.submit)]
    json.dump(dict(cmd=args, tag=a.tag, n=a.submit, jobs=ids, mem=mem, time=a.time), open(os.path.join(logs, f"submit_{a.tag}.json"), "w"))
    print(" ".join(ids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sets", default=None)
    ap.add_argument("--strata", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--timeout", type=float, default=600.0)
    ap.add_argument("--release", type=int, default=None)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--flush", type=int, default=2000)
    ap.add_argument("--submit", type=int, default=0)
    ap.add_argument("--mem", type=int, default=None)
    ap.add_argument("--time", default="08:00:00")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0" and not a.submit:
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)
    a.out, a.clone, a.graph = os.path.abspath(a.out), os.path.abspath(a.clone), os.path.abspath(a.graph)
    a.sets = a.sets.split(",") if a.sets else None
    a.strata = a.strata.split(",") if a.strata else None
    a.release = a.release or C.release_of_graph(a.graph)
    a.tag = a.tag or f"{'_'.join(a.sets or ['all'])}_{a.organism}"
    os.makedirs(a.out, exist_ok=True)
    if a.submit:
        return submit(a)
    run(a)


if __name__ == "__main__":
    main()
