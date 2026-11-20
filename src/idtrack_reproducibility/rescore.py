"""Score the stored scorings again under the judge of today (docs/records.md, "Scoring again").

Usage: python -m idtrack_reproducibility.rescore [--parts harness,f0,astar,a,b,c,d,e] [--organisms <o>,...]
                                                  [--trees <t>,...] [--workers N] [--dry-run]

Every scored folder that `results.toml` registers (the common scorer's, and the harness's own scorer's: part
`harness`, `harness.score --merge`) is scored again from its stored answers (nothing is asked again), written as the
record of today's version (`scorer.RECORD_VERSION`), with:

- the sets it read, and every set published after it that revises them (a published set is never rewritten: a
  change of the truth is a revision set, and the latest line of a qid is its truth); a decided judge's revision set
  (`<set>_dec<k>`, or a set whose manifest names the judge `decided`) only in a decided scoring (a folder
  `scores_decided`), after the others (`score_e.sets_for`); bounded to the qids of the sets it read (`--within`: a
  revision set may also revise lines of sets the run did not ask, e.g. the harness's travel lines beside its kinds);
- the strata, kinds and limit it was scored with; the partitions of the stored scoring, `--workers` at
  once (each a process of the scorer), then the merge (`run.json`, the reports). A partition stored under the same
  scoring key (`scorer.scoring_key`: the sets with their checksums, the selection, the record version and the code
  of the judge) is not scored again, nor a scoring whose `run.json` has that key: a run that stops is run again and
  goes on. Family E's summary (`e_summary.tsv`) is written again once both of its judges are scored.

The causes, the families' own tables and the results are made from the new scorings by their own commands, as each
family's document says ("Reproducing from scratch").
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time
import types

from . import config as CFG

PARTS = ("harness", "f0", "astar", "a", "b", "c", "d", "e")


def log(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def is_decided_set(family, organism, name):
    from .families import lines as L
    from .families.e import score_e as SE
    return bool(SE.DEC.match(name)) or L.manifest(family, organism, name).get("judge") == "decided"


def current_sets(family, organism, original, decided):
    """(sets in reading order, within): the sets a scoring that read `original` reads today (module docstring)."""
    from .families import lines as L
    names = L.set_names(family, organism)
    keep, orig = [], set(original)
    for n in names:
        dec = is_decided_set(family, organism, n)
        if dec and not decided and n not in orig:
            continue
        if n in orig or set(L.manifest(family, organism, n).get("revises") or []) & set(keep):
            keep.append(n)
    if decided:
        keep = [n for n in keep if not is_decided_set(family, organism, n)] + \
               [n for n in keep if is_decided_set(family, organism, n)]
    added = [n for n in keep if n not in orig]
    return keep, (sorted(orig) if added else None)


def answers_of(folder):
    """The run folder whose `answers/` a scored folder scored (itself, or its parent for `scores_<judge>/`)."""
    for d in (folder, os.path.dirname(folder)):
        if os.path.isdir(os.path.join(d, "answers")):
            return d
    raise SystemExit(f"{folder}: no answers beside it")


def scorings(parts=None, organisms=None, trees=None):
    """[(part, tree, organism, scored folder)] of the registry's runs scored by the common scorer."""
    from .records import registered_folders
    out = []
    for part, tree, org, folder in registered_folders(parts, organisms, trees, with_rows=False):
        if part in PARTS and os.path.exists(os.path.join(folder, "run.json")):
            out.append((part, tree, org, folder))
    return out


def plan(folder):
    """The scorer's arguments for scoring `folder` again (a namespace as scorer.main makes it)."""
    from .families import scorer as SC
    info = json.load(open(os.path.join(folder, "run.json")))
    family, org = info["family"], info["organism"]
    decided = os.path.basename(folder) == "scores_decided"
    sets, within = current_sets(family, org, list((info.get("sets") or {}).keys()), decided)
    run = answers_of(folder)
    a = types.SimpleNamespace(family=family, organism=org, run=run, scores_to=folder, sets=sets, within=within,
                              strata=info.get("strata"), kinds=info.get("kinds"), limit=info.get("limit"),
                              judge=info.get("judge") or "revised", part=None, merge=False)
    # the partitions of the stored scoring (reading the row counts of thousands of answer files again takes long)
    a.parts = 1 if a.limit else (info.get("parts") or SC.n_parts(run))
    a.key = SC.scoring_key(a)
    return a


def plan_harness(folder, tree):
    """The harness's own scorer's command for scoring `folder` again (`harness.score --merge`) and its key."""
    from .harness import score as HS
    from .families import lines as L
    info = json.load(open(os.path.join(folder, "run.json")))
    org = info["organism"]
    # the answers of `<harness>/runs/<tree>/<organism>/<part>` are `<harness>/answers/<tree>/<organism>/<part>` of this
    # data root (a run.json written before the data root moved names their former place)
    answers = os.path.join(CFG.get().harness, "answers", os.path.relpath(folder, os.path.join(CFG.get().harness, "runs")))
    sets, within = current_sets("harness", org, list((info.get("sets") or {}).keys()), False)
    man = {n: L.manifest("harness", org, n)["sha256"] for n in sets}
    a = types.SimpleNamespace(organism=org, judge=info.get("judge") or "revised", kinds=None, limit=None, within=within)
    t = CFG.get().tree(tree)
    cmd = ["-m", "idtrack_reproducibility.harness.score", "--clone", t["source"], "--graph", CFG.get().graph(tree, org),
           "--organism", org, "--out", answers, "--merge", "--scores-to", folder, "--sets", ",".join(sets),
           "--judge", a.judge, "--release", str(info.get("release") or 115)] + (["--within", ",".join(within)] if within else [])
    return cmd, HS.harness_key(a, man), sets, info


def args_of(a):
    out = ["--family", a.family, "--organism", a.organism, "--run", a.run, "--scores-to", a.scores_to,
           "--judge", a.judge, "--parts", str(a.parts)]
    for k in ("sets", "within", "strata", "kinds"):
        if getattr(a, k):
            out += [f"--{k}", ",".join(getattr(a, k))]
    if a.limit:
        out += ["--limit", str(a.limit)]
    return out


def done_key(path):
    try:
        return json.load(open(path)).get("key")
    except (OSError, ValueError):
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.rescore", description=__doc__.split("\n\n")[0])
    ap.add_argument("--parts", default=None, help=f"the parts (default: {','.join(PARTS)})")
    ap.add_argument("--organisms", default=None)
    ap.add_argument("--trees", default=None, help="the registry's tree names (e.g. d4c2dda,d4c2dda_s110)")
    ap.add_argument("--workers", type=int, default=1, help="partitions scored at once")
    ap.add_argument("--dry-run", action="store_true", help="print what would be scored, score nothing")
    a = ap.parse_args(argv)
    sp = lambda x: x.split(",") if x else None  # noqa: E731
    todo = scorings(sp(a.parts), sp(a.organisms), sp(a.trees))
    py = sys.executable
    tasks, merges, summaries = [], {}, {}  # summaries: family E's run folder -> organism
    for part, tree, org, folder in todo:
        rel = CFG.in_root(folder)
        if part == "harness":
            cmd, key, sets, info = plan_harness(folder, tree)
            if done_key(os.path.join(folder, "run.json")) == key:
                log("current", rel)
                continue
            log("to score", rel, f"sets added {[x for x in sets if x not in info.get('sets', {})]}")
            merges[rel] = [py] + cmd
            continue
        p = plan(folder)
        if done_key(os.path.join(folder, "run.json")) == p.key:
            log("current", rel)
            continue
        added = [s for s in (p.sets or []) if s not in json.load(open(os.path.join(folder, "run.json"))).get("sets", {})]
        log("to score", rel, f"parts {p.parts}", f"sets added {added}" if added else "", f"within {p.within}" if p.within else "")
        for i in range(p.parts):
            if done_key(os.path.join(folder, "outcomes", f"part_{i}_of_{p.parts}.json")) != p.key:
                tasks.append((rel, [py, "-m", "idtrack_reproducibility.families.scorer"] + args_of(p) + ["--part", str(i)]))
        merges[rel] = [py, "-m", "idtrack_reproducibility.families.scorer"] + args_of(p) + ["--merge"]
        if part == "e":
            summaries[os.path.dirname(folder)] = org
    log(f"{len(tasks)} partitions and {len(merges)} merges to run")
    if a.dry_run:
        for rel, cmd in tasks[:3] + [(k, v) for k, v in list(merges.items())[:2]]:
            log("would run", " ".join(cmd))
        return 0
    pending = {rel: sum(1 for r, _ in tasks if r == rel) for rel in merges}
    queue = list(tasks)
    running, failed = [], []
    logs = os.path.join(CFG.get().root, "cache", "rescore_logs")
    os.makedirs(logs, exist_ok=True)

    def start(rel, cmd, what):
        lp = os.path.join(logs, rel.replace("/", "__") + f"__{what}.log")
        running.append((rel, what, subprocess.Popen(cmd, stdout=open(lp, "w"), stderr=subprocess.STDOUT), lp))

    ready = [rel for rel, n in pending.items() if n == 0]
    while queue or running or ready:
        while ready and len(running) < a.workers:
            rel = ready.pop(0)
            start(rel, merges[rel], "merge")
        while queue and len(running) < a.workers:
            rel, cmd = queue.pop(0)
            start(rel, cmd, "part_" + cmd[-1])
        time.sleep(5)
        for item in list(running):
            rel, what, proc, lp = item
            if proc.poll() is None:
                continue
            running.remove(item)
            if proc.returncode:
                failed.append((rel, what, lp))
                log("FAILED", rel, what, lp)
                continue
            log("done", rel, what)
            if what != "merge":
                pending[rel] -= 1
                if pending[rel] == 0 and not any(r == rel for r, _, _ in failed):
                    ready.append(rel)
    if failed:
        log(f"{len(failed)} failed:", failed)
        return 1
    from .families.e import score_e as SE
    for base, org in summaries.items():  # family E's summary over the judges of the run
        judges = sorted(os.path.basename(d)[len("scores_"):] for d in glob.glob(os.path.join(base, "scores_*"))
                        if os.path.exists(os.path.join(d, "run.json")))
        SE.summary(base, org, judges)
        log("summary", CFG.in_root(base), judges)
    log("all scored")
    return 0


if __name__ == "__main__":
    sys.exit(main())
