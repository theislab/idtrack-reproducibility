"""`idtrack-repro modest`: the modest set, a measurement of an idtrack tree in a few hours on one job per organism.

Per organism it holds (docs/modest.md):

- the accuracy harness's published sets, every line (the judge's fixed queries; through the common runner, family
  `harness`);
- every reproducer of the weakness hunts (`idtrack_reproducibility.hunts.reproduce`, the organism's calls);
- per family (A*, A, B, C, D, E) a seeded stratified sample: up to `--per-stratum` lines (default 500) of every stratum
  of the family's own sets, the lines with the smallest hash of (seed, family, organism, qid), with the family's hunt
  lines. It is a published selection of qids (`<families>/<family>/panels/<organism>/modest_v<k>.parquet`, built by
  the panel machinery, `families/panel.py`); the truth of a line is the latest line of its qid in the family's sets.

Usage:
  idtrack-repro modest size [--per-stratum 500]          the size, from the sets' manifests (reads no set)
  idtrack-repro modest build [--per-stratum 500] [--families a,e] [--organisms ...] [--dry]
  idtrack-repro modest run --tree <name> --out <folder> [--organisms ...] [--families ...] [--shards d=8]
                                                           one job per organism (the configured runner)
  idtrack-repro modest organism --tree <name> --organism <o> --out <folder> [--shards d=8] [--score-only]
                                                           the job's body: every part in turn (--score-only: the
                                                           stored answers scored again under today's judge)
  idtrack-repro modest compare <folder a> <folder b> [--out <tsv>]
                                                           every part of two trees' runs: rates, regressions, moves

A run's folder: `<out>/<organism>/<part>/` with the answers and the scored tables of each part (`harness`, the
families, `hunts.json`), and `<out>/<organism>/modest.json` (the parts, their exit codes, times, the answers taken
from other runs of the same tree (`reused`: the runner's `reuse`), a part killed by a signal). A part may be asked in
shards at once (`--shards d=8`, `[jobs.modest_shards]`): the runner's own sharding, each shard a process of the job.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time

from . import config as CFG

FAMILIES = ["astar", "a", "b", "c", "d", "e"]
PER_STRATUM = 500
SELECTION = "modest"
# The runner of each family's lines (A*, A, B, C use the common runner; D and E have their own, with the same options).
RUNNER = {"d": "idtrack_reproducibility.families.d.run_d", "e": "idtrack_reproducibility.families.e.run_e"}


def selection(family, organism):
    """The latest modest selection of a family and organism, or None."""
    from .families import lines as L
    names = [n for n in L.panel_names(family, organism) if n.startswith(f"{SELECTION}_v")]
    return names[-1] if names else None


def stratum_sizes(family, organism):
    """{stratum: lines} of the family's sets (revision sets left out: they repeat qids), from their manifests."""
    from .families import lines as L
    d = L.query_dir(family, organism)
    sizes = {}
    for p in glob.glob(os.path.join(d, "*.manifest.json")):
        m = json.load(open(p))
        if m.get("revises"):
            continue
        for s, n in (m.get("counts") or {}).items():
            sizes[s] = sizes.get(s, 0) + n
    return sizes


def harness_lines(organism):
    from .harness import harness_lib as H
    d = H.query_dir(organism)
    seen = 0
    for s in H.set_names(d) if os.path.isdir(d) else []:
        m = json.load(open(os.path.join(d, f"{s}.manifest.json")))
        if not m.get("revises"):
            seen += m.get("n") or sum((m.get("counts") or m.get("kinds") or {}).values())
    return seen


def size(per_stratum=PER_STRATUM, families=FAMILIES, organisms=CFG.ORGANISMS):
    rows = []
    for o in organisms:
        rows.append(("harness", o, harness_lines(o), "all"))
        for f in families:
            sz = stratum_sizes(f, o)
            if sz:
                rows.append((f, o, sum(min(per_stratum, n) for n in sz.values()), f"{len(sz)} strata"))
    return rows


def cmd_size(a):
    rows = size(a.per_stratum, a.families, a.organisms)
    tot = {}
    for f, o, n, what in rows:
        print(f"{f}\t{o}\t{n}\t{what}")
        tot[f] = tot.get(f, 0) + n
    print("total\t" + "\t".join(f"{f} {n:,}" for f, n in tot.items()) + f"\tall {sum(tot.values()):,} lines")


def cmd_build(a):
    from .families import panel as P
    for f in a.families:
        for o in a.organisms:
            if not stratum_sizes(f, o):
                continue
            P.build(f, o, folders=[], per_stratum=a.per_stratum, dry=a.dry, kind=SELECTION)


def run(cmd, log):
    t0 = time.time()
    with open(log, "a") as fo:
        fo.write(f"$ {cmd}\n")
        fo.flush()
        rc = subprocess.run(cmd, shell=True, stdout=fo, stderr=subprocess.STDOUT).returncode
    return rc, round(time.time() - t0, 1)


def killed_by(rc):
    """The signal that killed a command (`subprocess` gives -N; a shell that ran it gives 128 + N), else None."""
    if rc < 0:
        return -rc
    if rc > 128:
        return rc - 128
    return None


def asked(cmd, n):
    """A runner's command for `n` shards: one command (`--shard 0/1`), or a list of the n shards' commands, run at once."""
    return f"{cmd} --shard 0/1" if n <= 1 else [f"{cmd} --shard {i}/{n}" for i in range(n)]


def shards_of(spec=None):
    """{part: shards} of the configuration (`[jobs.modest_shards]`), updated by `spec` (`d=8,e=2`)."""
    out = dict(CFG.get().modest_shards)
    for x in (spec or "").split(","):
        if x.strip():
            k, v = x.split("=")
            out[k.strip()] = int(v)
    return out


def organism_parts(tree, organism, out, families=FAMILIES, shards=None):
    """The commands of one organism's modest run, in order: (part, [commands]); a command that is a list is the
    shards of one part's runner, run at once (`shards`: {part: n})."""
    cfg = CFG.get()
    shards = shards or {}
    t = cfg.tree(tree)
    graph = cfg.graph(tree, organism)
    py_tree, py = t["python"], cfg.python
    base = os.path.join(out, organism)
    common = f"--organism {organism} --clone {t['source']} --graph {graph}"
    parts = [("harness", [asked(f"{py_tree} -m idtrack_reproducibility.families.runner --family harness {common} "
                                f"--out {base}/harness", shards.get("harness", 1)),
                          f"{py} -m idtrack_reproducibility.families.scorer --family harness --organism {organism} "
                          f"--run {base}/harness"])]
    for f in families:
        sel = selection(f, organism)
        if not sel:
            continue
        runner = RUNNER.get(f, "idtrack_reproducibility.families.runner")
        fam = "" if f in RUNNER else f"--family {f} "
        cmds = [asked(f"{py_tree} -m {runner} {fam}{common} --out {base}/{f} --sets {sel}", shards.get(f, 1))]
        if f == "e":  # family E is judged twice: its sets without the `_dec` revision sets, and with them
            from .families.e import score_e as SE
            for judge in ("now", "decided"):
                sets = ",".join([sel] + SE.sets_for(organism, judge))
                cmds.append(f"{py} -m idtrack_reproducibility.families.scorer --family e --organism {organism} "
                            f"--run {base}/e --sets {sets} --scores-to {base}/e/scores_{judge}")
        else:
            cmds.append(f"{py} -m idtrack_reproducibility.families.scorer --family {f} --organism {organism} "
                        f"--run {base}/{f} --sets {sel}")
        parts.append((f, cmds))
    parts.append(("hunts", [f"{py} -m idtrack_reproducibility.hunts.reproduce all all --tree {tree} "
                            f"--organism {organism} --json {base}/hunts.json"]))
    return parts


def run_shards(cmds, log):
    """The shards of one part's runner at once (each to its own log beside `log`); the worst exit (a signal first)."""
    t0 = time.time()
    procs = []
    with open(log, "a") as fo:
        for k, c in enumerate(cmds):
            fo.write(f"$ {c}   (shard log {log}.s{k})\n")
    for k, c in enumerate(cmds):
        procs.append(subprocess.Popen(c, shell=True, stdout=open(f"{log}.s{k}", "a"), stderr=subprocess.STDOUT))
    rcs = [p.wait() for p in procs]
    killed = [rc for rc in rcs if killed_by(rc)]
    return (killed[0] if killed else max(rcs, key=abs)), round(time.time() - t0, 1)


def reused_of(folder):
    """{source run folder: lines} the part's runner took from other runs of the tree (its load records' `reused`)."""
    out = {}
    for p in glob.glob(os.path.join(folder, "answers", "load_*.json")):
        for k, v in (json.load(open(p)).get("reused") or {}).items():
            out[k] = out.get(k, 0) + v
    return out


def cmd_organism(a):
    base = os.path.join(a.out, a.organism)
    os.makedirs(base, exist_ok=True)
    log = os.path.join(base, "modest.log")
    done, killed = [], None
    for part, cmds in organism_parts(a.tree, a.organism, a.out, a.families, shards_of(getattr(a, "shards", None))):
        if getattr(a, "score_only", False):  # the stored answers scored again (the judge of today); nothing asked
            cmds = [c for c in cmds if isinstance(c, str) and "families.scorer" in c]
        for c in cmds:
            rc, secs = run_shards(c, log) if isinstance(c, list) else run(c, log)
            first = c[0] if isinstance(c, list) else c
            done.append(dict(part=part, command=first.split(" -m ", 1)[-1].split(" ", 1)[0], exit=rc, seconds=secs,
                             **({"shards": len(c)} if isinstance(c, list) else {})))
            if " --shard " in first:
                done[-1]["reused"] = reused_of(os.path.join(base, part))
            sig = killed_by(rc)
            if sig:  # killed (e.g. by the scheduler for memory): the run stops, never goes on beside a lost part
                killed = dict(part=part, command=done[-1]["command"], signal=sig)
                break
            if rc and part != "hunts":  # a reproducer exits 1 when a finding changed: a result, not a failure
                break
        if killed:
            break
    tr = CFG.get().tree(a.tree)
    rec = os.path.join(base, "modest.json")
    if not getattr(a, "score_only", False) and os.path.exists(rec):  # an earlier run of this organism keeps its record
        k = len(glob.glob(os.path.join(base, "modest.earlier_*.json"))) + 1
        os.replace(rec, os.path.join(base, f"modest.earlier_{k}.json"))
    json.dump(dict(tree=a.tree, commit=tr["commit"], organism=a.organism, parts=done, killed=killed,
                   finished=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(base, "modest_scored_again.json" if getattr(a, "score_only", False) else "modest.json"), "w"),
              indent=1)  # a scoring again keeps the run's own record
    if killed:
        sys.stderr.write(f"{a.organism}: part {killed['part']} ({killed['command']}) was killed by signal {killed['signal']} "
                         f"(SIGKILL: most often the job's memory; see `modest run` and [jobs.modest_memory_gb]); "
                         f"the run stopped, the parts after it were not run\n")
        return 2
    bad = [d for d in done if d["exit"] and d["part"] != "hunts"]
    print(f"{a.organism}: {len(done)} steps, {len(bad)} failed")
    return 1 if bad else 0


def cmd_run(a):
    from . import jobs
    cfg = CFG.get()
    ids = []
    shards = shards_of(a.shards)
    n = max([1] + [v for k, v in shards.items() if k == "harness" or k in a.families])
    for o in a.organisms:
        # the parts run one after another: the largest part's measured peak; a part's n shards at once each hold the graph
        mem = a.mem or max(cfg.modest_memory_gb[o], n * cfg.memory_gb[o] if n > 1 else 0)
        fam = f" --families {','.join(a.families)}" if a.families != FAMILIES else ""
        sh = f" --shards {a.shards}" if a.shards else ""
        ids.append(jobs.submit(f"modest_{a.tree}_{o}", jobs.module(cfg.python, "idtrack_reproducibility.modest",
                                                                   f"organism --tree {a.tree} --organism {o} "
                                                                   f"--out {os.path.abspath(a.out)}{fam}{sh}"),
                               os.path.join(os.path.abspath(a.out), "logs"), cpus=max(2, n), mem_gb=mem, time=a.time))
    print(" ".join(ids))


def cmd_compare(a):
    rows, regress = [], 0
    for org_dir in sorted(glob.glob(os.path.join(a.a, "*", "modest.json"))):
        o = os.path.basename(os.path.dirname(org_dir))
        for part_dir in sorted(glob.glob(os.path.join(a.a, o, "*"))):
            part = os.path.basename(part_dir)
            for sub in (["scores_now", "scores_decided"] if part == "e" else [""]):
                pa, pb = os.path.join(part_dir, sub), os.path.join(a.b, o, part, sub)
                if not os.path.exists(os.path.join(pa, "report.tsv")) or not os.path.exists(os.path.join(pb, "report.tsv")):
                    continue
                cmd = [sys.executable, "-m", "idtrack_reproducibility.families.compare", pa, pb,
                       "--moved", os.path.join(pb, "moved_vs_" + os.path.basename(os.path.normpath(a.a)) + ".tsv")]
                r = subprocess.run(cmd, capture_output=True, text=True)
                regress += r.returncode == 1
                rows.append((o, part + (f"/{sub}" if sub else ""), "regression" if r.returncode == 1 else
                             ("ok" if r.returncode == 0 else f"error {r.returncode}")))
    for r in rows:
        print("\t".join(r))
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("organism\tpart\tverdict\n" + "".join("\t".join(r) + "\n" for r in rows))
    return 1 if regress else 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro modest", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="what", required=True)

    def lists(p):
        p.add_argument("--families", default=",".join(FAMILIES))
        p.add_argument("--organisms", default=",".join(CFG.ORGANISMS))
    s = sub.add_parser("size")
    s.add_argument("--per-stratum", type=int, default=PER_STRATUM)
    lists(s)
    b = sub.add_parser("build")
    b.add_argument("--per-stratum", type=int, default=PER_STRATUM)
    b.add_argument("--dry", action="store_true")
    lists(b)
    r = sub.add_parser("run")
    r.add_argument("--tree", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--time", default="08:00:00")
    r.add_argument("--mem", type=int, default=None, help="GB per job (default: the configuration's [jobs.modest_memory_gb])")
    r.add_argument("--shards", default=None, help="parts asked in shards at once, `d=8,e=2` (beside [jobs.modest_shards])")
    lists(r)
    o = sub.add_parser("organism")
    o.add_argument("--tree", required=True)
    o.add_argument("--organism", required=True)
    o.add_argument("--out", required=True)
    o.add_argument("--families", default=",".join(FAMILIES))
    o.add_argument("--shards", default=None)
    o.add_argument("--score-only", action="store_true", help="score the stored answers again (today's judge); ask nothing, "
                   "run no reproducer")
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    for k in ("families", "organisms"):
        if hasattr(a, k):
            setattr(a, k, getattr(a, k).split(","))
    return {"size": cmd_size, "build": cmd_build, "run": cmd_run, "organism": cmd_organism,
            "compare": cmd_compare}[a.what](a)


if __name__ == "__main__":
    sys.exit(main())
