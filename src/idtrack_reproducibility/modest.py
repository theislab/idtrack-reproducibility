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
  idtrack-repro modest run --tree <name> --out <folder> [--organisms ...] [--families ...]
                                                           one job per organism (the configured runner)
  idtrack-repro modest organism --tree <name> --organism <o> --out <folder>   (the job's body: every part in turn)
  idtrack-repro modest compare <folder a> <folder b> [--out <tsv>]
                                                           every part of two trees' runs: rates, regressions, moves

A run's folder: `<out>/<organism>/<part>/` with the answers and the scored tables of each part (`harness`, the
families, `hunts.json`), and `<out>/<organism>/modest.json` (the parts, their exit codes, times).
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


def organism_parts(tree, organism, out, families=FAMILIES):
    """The commands of one organism's modest run, in order: (part, [commands])."""
    cfg = CFG.get()
    t = cfg.tree(tree)
    graph = cfg.graph(tree, organism)
    py_tree, py = t["python"], cfg.python
    base = os.path.join(out, organism)
    common = f"--organism {organism} --clone {t['source']} --graph {graph}"
    parts = [("harness", [f"{py_tree} -m idtrack_reproducibility.families.runner --family harness {common} "
                          f"--out {base}/harness --shard 0/1",
                          f"{py} -m idtrack_reproducibility.families.scorer --family harness --organism {organism} "
                          f"--run {base}/harness"])]
    for f in families:
        sel = selection(f, organism)
        if not sel:
            continue
        runner = RUNNER.get(f, "idtrack_reproducibility.families.runner")
        fam = "" if f in RUNNER else f"--family {f} "
        cmds = [f"{py_tree} -m {runner} {fam}{common} --out {base}/{f} --sets {sel} --shard 0/1"]
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


def cmd_organism(a):
    base = os.path.join(a.out, a.organism)
    os.makedirs(base, exist_ok=True)
    log = os.path.join(base, "modest.log")
    done = []
    for part, cmds in organism_parts(a.tree, a.organism, a.out, a.families):
        for c in cmds:
            rc, secs = run(c, log)
            done.append(dict(part=part, command=c.split(" -m ", 1)[-1].split(" ", 1)[0], exit=rc, seconds=secs))
            if rc and part != "hunts":  # a reproducer exits 1 when a finding changed: a result, not a failure
                break
    tr = CFG.get().tree(a.tree)
    json.dump(dict(tree=a.tree, commit=tr["commit"], organism=a.organism, parts=done,
                   finished=time.strftime("%Y-%m-%d %H:%M:%S")), open(os.path.join(base, "modest.json"), "w"), indent=1)
    bad = [d for d in done if d["exit"] and d["part"] != "hunts"]
    print(f"{a.organism}: {len(done)} steps, {len(bad)} failed")
    return 1 if bad else 0


def cmd_run(a):
    from . import jobs
    cfg = CFG.get()
    ids = []
    for o in a.organisms:
        mem = jobs.memory_for(o) + 16
        fam = f" --families {','.join(a.families)}" if a.families != FAMILIES else ""
        ids.append(jobs.submit(f"modest_{a.tree}_{o}", jobs.module(cfg.python, "idtrack_reproducibility.modest",
                                                                   f"organism --tree {a.tree} --organism {o} "
                                                                   f"--out {os.path.abspath(a.out)}{fam}"),
                               os.path.join(os.path.abspath(a.out), "logs"), cpus=2, mem_gb=mem, time=a.time))
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
    lists(r)
    o = sub.add_parser("organism")
    o.add_argument("--tree", required=True)
    o.add_argument("--organism", required=True)
    o.add_argument("--out", required=True)
    o.add_argument("--families", default=",".join(FAMILIES))
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
