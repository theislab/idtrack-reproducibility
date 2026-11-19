"""H: run measurements of trees, interleaved and repeated, as jobs (docs/families/h.md, "Reproducing from scratch").

    python -m idtrack_reproducibility.families.h.run --organisms <o,...|all> --measures <m,...> [--trees d4c2dda,71690ab]
        [--reps 3] [--submit] [--cpus 4] [--mem 16] [--time 12:00:00] [--build cold|warm] [--warm-from <folder>]

Measures: `load`, `calls`, `harmonize`, `cycles` and `build` per organism; `server` for the organism `all` (the five
graphs in one process). One job per organism (with `--submit`, through the configured runner, `jobs.submit`, named
`h_<organism>_<measures>`, a build named with its kind; without it, here) runs, for repetition k = 1..`--reps`, every tree in turn, the order
reversed every other repetition (two trees, three repetitions: A B, B A, A B), each measure in a fresh process of the
tree's Python (`[trees.<tree>] python`), so that the trees are measured on the same node at nearly the same time.
A measurement writes `<performance>/runs/<tree>/<organism>/<measure>_rep<k>.json` (`build_<kind>_rep<k>.json`,
`server_rep<k>.json`); one that exists is not run again, so a preempted job is simply submitted again.
"""

import argparse
import os
import subprocess
import sys
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import jobs
from idtrack_reproducibility.corpus import _common as CC
from idtrack_reproducibility.families.h import common as H

PER_ORGANISM = ("load", "calls", "harmonize", "cycles", "build")


def order(trees, reps):
    """(tree, rep) in the interleaved order: the order of the trees reversed every other repetition."""
    out = []
    for k in range(1, reps + 1):
        out += [(t, k) for t in (trees if k % 2 else list(reversed(trees)))]
    return out


def measure_name(measure, a):
    return f"build_{a.build}" if measure == "build" else measure


def run_one(tree, organism, measure, rep, a, log):
    name = measure_name(measure, a)
    out = H.out_file(tree, organism, name, rep)
    if os.path.exists(out):
        return "exists"
    args = ["--tree", tree, "--out", out]
    if measure != "server":
        args += ["--organism", organism]
    if measure == "build":
        args += ["--cold"] if a.build == "cold" else ["--warm-from", a.warm_from or os.path.dirname(
            os.path.realpath(CC.graph_file(tree, organism)))]
    t = CC.tree(tree)
    env = dict(os.environ, **jobs.job_env({"PYTHONHASHSEED": "0"}))
    env["PYTHONPATH"] = CFG.SRC + (os.pathsep + os.environ["PYTHONPATH"] if os.environ.get("PYTHONPATH") else "")
    cmd = [t["python"], "-m", f"idtrack_reproducibility.families.h.{measure}"] + args
    with open(log, "a") as fo:
        fo.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(cmd)}\n")
        fo.flush()
        r = subprocess.run(cmd, stdout=fo, stderr=subprocess.STDOUT, env=env)
    return "ok" if r.returncode == 0 else f"failed ({r.returncode})"


def run_organism(organism, a):
    logs = os.path.join(H.root(), "runs", "logs")
    os.makedirs(logs, exist_ok=True)
    log = os.path.join(logs, f"h_{organism}_{'_'.join(measure_name(m, a) for m in a.measures)}.log")
    status = []
    for tree, rep in order(a.trees, a.reps):
        for m in a.measures:
            s = run_one(tree, organism, m, rep, a, log)
            status.append((tree, rep, m, s))
            print(f"{time.strftime('%H:%M:%S')} {organism} {tree} rep{rep} {m}: {s}", flush=True)
    return status


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="h.run", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organisms", required=True, help="comma-separated, or `all` for the server")
    ap.add_argument("--measures", required=True)
    ap.add_argument("--trees", default="d4c2dda,71690ab")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--build", choices=("cold", "warm"), default="warm")
    ap.add_argument("--warm-from", default=None, help="the folder of a build's tables (default: the graph's folder)")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--mem", type=int, default=16)
    ap.add_argument("--time", default="12:00:00")
    a = ap.parse_args(argv)
    a.organisms = a.organisms.split(",")
    a.measures = a.measures.split(",")
    a.trees = a.trees.split(",")
    for m in a.measures:
        if m not in PER_ORGANISM + ("server",):
            raise SystemExit(f"unknown measure {m}")
        if (m == "server") != (a.organisms == ["all"]):
            raise SystemExit("the measure `server` runs alone, with --organisms all")
    if a.submit:
        logs = os.path.join(H.root(), "runs", "logs")
        rest = [x for x in argv if x != "--submit"]
        i = rest.index("--organisms")
        for org in a.organisms:
            args = rest[:i] + ["--organisms", org] + rest[i + 2:]
            cmd = jobs.module(CFG.get().python, "idtrack_reproducibility.families.h.run", " ".join(args))
            what = "_".join(measure_name(m, a) for m in a.measures)
            print(org, jobs.submit(f"h_{org}_{what}", cmd, logs, cpus=a.cpus, mem_gb=a.mem, time=a.time))
        return 0
    failed = [s for org in a.organisms for s in run_organism(org, a) if s[3].startswith("failed")]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
