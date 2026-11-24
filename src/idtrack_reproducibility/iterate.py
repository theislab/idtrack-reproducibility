"""`idtrack-repro iterate <tree>`: one command per iteration of idtrack, a resumable plan of jobs (docs/iterate.md).

Usage:
  idtrack-repro iterate <tree> [--previous <tree>] [--reference <tree>] [--organisms <o,...|all>] [--tables <folder>]
                               [--snapshot 115] [--retry <step,...|failed>] [--dry]
  idtrack-repro iterate <tree> --report          the report again from what the steps left (no job)
  idtrack-repro iterate <tree> --complete [--dry]
                                                 the record of a complete run: the families' run folders of the tree

An iteration measures one tree (`[trees.<tree>]`) against the iteration before it (`[iterations] order`, or
`--previous`) and against the reference (`[iterations] reference`). Its steps, per organism where they ask a graph:

1. `graph_build_<o>`: `graph-build` of the tree's graph; done where the configured graph exists and the build record
   beside it (`build.json`) names that file and records its content digest;
2. `graph_judge_<o>`: `graph-judge --against <previous>` (after 1); done where `judge.json` judges the configured graph
   against the previous tree and is newer than the build record;
3. `modest_<o>`: the modest set's run of the organism (`modest organism`, after 1), each family's part on its panel
   where `<families>/<family>/panels/<o>/panel_v<k>` exists, else on its modest selection; the hunts' reproducers are its
   last part (`hunts.json`); done where `<modest>/<tree>/<o>/modest.json` is the configured commit's, no part killed
   and no part failed;
4. `compare_previous`, `compare_reference`: `modest compare` of the previous tree's and the reference's runs with the
   tree's (after every 3); done where `compare_vs_<other>.tsv` is newer than both runs' records;
5. `corpus_<o>`: the regression corpus's reference answers of the tree (`corpus reference`, after 1), compared with the
   latest earlier tree of the order that has references, and with the reference tree's (`corpus compare`); done where
   the tree's reference file and the comparison exist;
6. `report`: `iterate_report.md` and `iterate_report.tsv` in `<iterations>/<tree>/` (after every other step, whatever
   their end), written again by every call of the command too.

Steps 1 to 5 are chained `afterok` (a step does not run after a failed predecessor), the report `afterany`. The state
is `<iterations>/<tree>/iterate.json`: every step with its command, its dependencies, its attempts (job id, when, the
scheduler's state, elapsed time, peak memory) and its status. A call computes every step's status from its outputs
and its jobs (squeue, sacct): `done`, `active` (pending or running), `lost` (preempted, node failure: the scheduler's),
`failed` (failed, out of memory, out of time, or completed without its output), `cancelled`, `blocked` (pending on a
predecessor that can no longer succeed), `todo`; it submits the steps that are `todo`, `lost` or `blocked` (and those
named by `--retry`) under their names, cancelling a blocked job of this plan first (it can never run), and nothing
else: no job resubmits by itself. Before each submission the queue is read against `[jobs.limits]` (`max_jobs`,
`max_memory_gb`, `names`); a step that would exceed them is `held` and waits for the next call. `--dry` prints the plan
and writes and submits nothing.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import gzip
import json
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Callable

from . import config as CFG
from . import jobs

STATE = "iterate.json"
REPORT_MD, REPORT_TSV = "iterate_report.md", "iterate_report.tsv"
REPORT_COLS = ["section", "organism", "part", "against", "key", "value", "cause", "detail"]
# The scheduler's states, by what they mean for a step (sacct's and squeue's names).
ACTIVE = {"PENDING", "RUNNING", "REQUEUED", "CONFIGURING", "COMPLETING", "SUSPENDED", "RESIZING", "SIGNALING",
          "STAGE_OUT", "REQUEUE_FED", "REQUEUE_HOLD", "SPECIAL_EXIT"}
LOST = {"PREEMPTED", "NODE_FAIL", "BOOT_FAIL", "DEADLINE", "LAUNCH_FAILED", "RECONFIG_FAIL", "REVOKED"}
FAILED = {"FAILED", "OUT_OF_MEMORY", "TIMEOUT"}
# The outcomes a moved answer is read from and to (the regressions of `compare`).
WAS_RIGHT, NOW_BAD = ("right", "expected_error"), ("wrong", "error")
EXAMPLES = 100  # moved answers listed per comparison in the Markdown report (the TSV lists every one)
UNKNOWN_AFTER = 600  # seconds after its submission a job the scheduler does not know is taken as lost


@dataclass
class Step:
    """One step of the plan: a job's command, its resources, the steps it runs after, and how its output is judged."""
    name: str
    kind: str
    command: str
    organism: str | None = None
    after: list = field(default_factory=list)
    cpus: int = 2
    mem_gb: int = 4
    time: str = "04:00:00"
    dependency: str = "afterok"
    done: Callable[[], str | None] = lambda: None  # the reason the step's output is current, or None
    note: str = ""
    tree: str = ""


# ---------------------------------------------------------------------------------------------------------------------
# Places


def state_dir(tree):
    return os.path.join(CFG.get().iterations, tree)


def modest_dir(tree):
    return os.path.join(CFG.get().modest, tree)


def q(x):
    return shlex.quote(str(x))


def module(name, args):
    return jobs.module(CFG.get().python, name, args)


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def verdict_command(command, out):
    """A command whose exit 1 is a verdict written to `out` (a fault found, a regression), not a failure: the job exits 0
    when the command exits 0, or exits 1 having written `out` during the job; any other exit is a failure."""
    return (f'start=$(( $(date +%s) - 1 )); ( {command} ); rc=$?; [ "$rc" -eq 0 ] || {{ [ "$rc" -eq 1 ] && '
            f'[ -n "$(find {q(out)} -maxdepth 0 -newermt "@$start" 2>/dev/null)" ]; }}')


# ---------------------------------------------------------------------------------------------------------------------
# The steps' outputs


def graph_built(tree, organism):
    """Why the tree's configured graph of the organism is built, or None."""
    from .graphs import common as G
    graph = CFG.get().tree(tree)["graphs"].get(organism)
    record = load_json(os.path.join(G.out_dir(tree, organism), "build.json"))
    if not graph or not os.path.exists(graph) or not record or not record.get("content_digest"):
        return None
    if os.path.realpath(os.path.join(CFG.get().graphs, record.get("graph_file", ""))) != os.path.realpath(graph):
        return None
    return f"built, content digest {record['content_digest'][:12]}"


def graph_judged(tree, organism, against):
    from .graphs import common as G
    out = G.out_dir(tree, organism)
    judge, built = load_json(os.path.join(out, "judge.json")), mtime(os.path.join(out, "build.json"))
    graph = CFG.get().tree(tree)["graphs"].get(organism)
    if not judge or not graph or os.path.realpath(judge.get("graph", "")) != os.path.realpath(graph):
        return None
    if (judge.get("diff") or {}).get("against") != against or (built and mtime(os.path.join(out, "judge.json")) < built):
        return None
    return "judged: " + ("holds" if judge.get("holds") else f"faults {judge.get('faults')}")


def modest_run(tree, organism):
    """Why the tree's modest run of the organism is complete, or None (with the reason it is not, as a second value)."""
    rec = load_json(os.path.join(modest_dir(tree), organism, "modest.json"))
    if not rec:
        return None, "no run"
    commit = CFG.get().tree(tree)["commit"]
    if commit and not (rec.get("commit") or "").startswith(commit) and not commit.startswith(rec.get("commit") or "-"):
        return None, f"the stored run is of {str(rec.get('commit'))[:7]}, the configuration names {commit[:7]}"
    if rec.get("killed"):
        return None, f"part {rec['killed'].get('part')} was killed"
    bad = [p["part"] for p in rec.get("parts", []) if p.get("exit") and p.get("part") != "hunts"]
    if bad:
        return None, f"parts failed: {sorted(set(bad))}"
    return f"run of {str(rec.get('commit'))[:7]}, {len(rec.get('parts', []))} steps", None


def records_of(tree, organism):
    """The modest run's own records of an organism (the run's, its scorings again, its parts asked again)."""
    return glob.glob(os.path.join(modest_dir(tree), organism, "modest*.json"))


def compared(tree, other, organisms):
    out = os.path.join(modest_dir(tree), f"compare_vs_{other}.tsv")
    t = mtime(out)
    if t is None:
        return None
    newest = [mtime(p) for o in organisms for x in (tree, other) for p in records_of(x, o)]
    return f"compared ({os.path.basename(out)})" if all(m is None or m <= t for m in newest) else None


def corpus_reference_of(tree, organism):
    from .corpus import _common as C
    return C.reference_file(organism, tree)


def corpus_done(tree, organism, against_reference):
    from .corpus import _common as C
    if not os.path.exists(corpus_reference_of(tree, organism)):
        return None
    if against_reference and not os.path.exists(os.path.join(C.run_dir(tree, organism),
                                                             f"compare_reference_vs_{against_reference}.json")):
        return None
    return "reference answers made" + (f", compared with {against_reference}" if against_reference else "")


def panel_selections(organism, families):
    """{family: its latest panel} where a family's panel of the organism is published (docs/methods.md, "Measurement
    design"): the part asks it instead of the modest selection."""
    from .families import lines as L
    out = {}
    for f in families:
        names = [n for n in L.panel_names(f, organism) if n.startswith("panel_v")]
        if names:
            out[f] = names[-1]
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The plan


def organisms_of(tree, organisms=None):
    if organisms and organisms != ["all"]:
        return organisms
    configured = CFG.get().tree(tree)["graphs"]
    return [o for o in CFG.ORGANISMS if o in configured] or list(CFG.ORGANISMS)


def plan(tree, previous=None, reference=None, organisms=None, tables=None, snapshot=115, families=None):
    """The steps of the tree's iteration, in the order they are submitted."""
    from .corpus import _common as C
    from .graphs import common as G
    from . import modest as M
    cfg = CFG.get()
    t = cfg.tree(tree)
    organisms = organisms_of(tree, organisms)
    families = families or M.FAMILIES
    steps = []
    for o in organisms:
        if o not in t["graphs"]:
            raise SystemExit(f"tree {tree}: no graph of {o} configured; add [trees.{tree}.graphs] {o} = "
                             f"\"<graphs>/{tree}/{o}/graph_{o}_assembly<key>_min<first>_max{snapshot}.npz\" (the file "
                             f"graph-build writes there), or leave {o} out (--organisms)")
        build = f"graph_build_{o}"
        args = f"--tree {q(tree)} --organisms {o} --snapshot {snapshot}"
        if tables:
            args += f" --tables {q(os.path.abspath(tables))}"
        note = ""
        graph_dir = os.path.realpath(G.out_dir(tree, o))
        if os.path.exists(t["graphs"][o]) and not os.path.exists(os.path.join(G.out_dir(tree, o), "build.json")):
            note = ("the configured graph exists without a build record: graph-build keeps that file and records the "
                    "digest of the graph it builds beside it")
        elif not os.path.realpath(t["graphs"][o]).startswith(graph_dir + os.sep):
            note = f"the configured graph lies outside {G.out_dir(tree, o)}, where graph-build writes"
        steps.append(Step(build, "graph_build", module("idtrack_reproducibility.graphs.build", args), o,
                          mem_gb=G.MEMORY_GB.get(o, 16), time=G.TIME.get(o, "06:00:00"),
                          done=lambda o=o: graph_built(tree, o), note=note))
        against = previous if previous and o in cfg.tree(previous)["graphs"] else None
        judge = f"--tree {q(tree)} --organisms {o}" + (f" --against {q(against)}" if against else "")
        steps.append(Step(f"graph_judge_{o}", "graph_judge",
                          verdict_command(module("idtrack_reproducibility.graphs.judge", judge),
                                          os.path.join(G.out_dir(tree, o), "judge.json")), o, after=[build],
                          mem_gb=G.MEMORY_GB.get(o, 16), time=G.TIME.get(o, "06:00:00"),
                          done=lambda o=o, against=against: graph_judged(tree, o, against),
                          note="" if against or not previous else f"{previous} has no graph of {o}: no diff"))
        sel = panel_selections(o, families)
        args = f"organism --tree {q(tree)} --organism {o} --out {q(modest_dir(tree))}"
        if families != M.FAMILIES:
            args += f" --families {','.join(families)}"
        if sel:
            args += " --selections " + ",".join(f"{f}={s}" for f, s in sorted(sel.items()))
        steps.append(Step(f"modest_{o}", "modest", module("idtrack_reproducibility.modest", args), o, after=[build],
                          mem_gb=cfg.modest_memory_gb.get(o, 8), time="08:00:00",
                          done=lambda o=o: modest_run(tree, o)[0],
                          note="; ".join(x for x in (modest_run(tree, o)[1], ("panels: " + ", ".join(
                              f"{f} {s}" for f, s in sorted(sel.items()))) if sel else "") if x)))
    runs = [f"modest_{o}" for o in organisms]
    for name, other in (("compare_previous", previous), ("compare_reference", reference)):
        if not other or other == tree or (name == "compare_reference" and other == previous):
            continue
        if not os.path.isdir(modest_dir(other)):
            steps.append(Step(name, "compare", "", after=runs, done=lambda: None,
                              note=f"no modest run of {other} ({modest_dir(other)}): not compared"))
            continue
        out = os.path.join(modest_dir(tree), f"compare_vs_{other}.tsv")
        cmd = module("idtrack_reproducibility.modest", f"compare {q(modest_dir(other))} {q(modest_dir(tree))} --out {q(out)}")
        steps.append(Step(name, "compare", verdict_command(cmd, out), after=runs, mem_gb=8, time="04:00:00",
                          done=lambda other=other: compared(tree, other, organisms), note=f"against {other}"))
    order = cfg.iteration["order"]
    earlier = order[:order.index(tree)] if tree in order else ([previous] if previous else [])
    for o in organisms:
        if not os.path.exists(C.corpus_file(o)):
            continue
        with_ref = [x for x in earlier if os.path.exists(corpus_reference_of(x, o))]
        nearest = with_ref[-1] if with_ref else None
        also = reference if reference and reference != nearest and reference != tree and \
            os.path.exists(corpus_reference_of(reference, o)) else None
        cmds = []
        if not os.path.exists(corpus_reference_of(tree, o)):
            cmds.append(module("idtrack_reproducibility.corpus", f"reference --tree {q(tree)} --organism {o}"
                               + (f" --compare-with {q(nearest)}" if nearest else "")))
        if also:
            cmds.append(module("idtrack_reproducibility.corpus",
                               f"compare --tree {q(tree)} --organism {o} --tag reference --against {q(also)}"))
        steps.append(Step(f"corpus_{o}", "corpus", " && ".join(cmds), o, after=[f"graph_build_{o}"],
                          mem_gb=jobs.memory_for(o), time="04:00:00",
                          done=lambda o=o, also=also: corpus_done(tree, o, also),
                          note=f"compared with {', '.join(x for x in (nearest, also) if x) or 'no earlier reference'}"))
    rest = [s.name for s in steps]
    report = f"{q(tree)} --report" + (f" --previous {q(previous)}" if previous else "") + \
        (f" --reference {q(reference)}" if reference else "") + f" --organisms {','.join(organisms)}" + \
        (f" --families {','.join(families)}" if families != M.FAMILIES else "") + \
        (f" --tables {q(os.path.abspath(tables))}" if tables else "") + f" --snapshot {snapshot}"
    steps.append(Step("report", "report", module("idtrack_reproducibility.iterate", report), after=rest, mem_gb=4,
                      time="01:00:00", dependency="afterany"))
    return steps


# ---------------------------------------------------------------------------------------------------------------------
# The scheduler


def user():
    return os.environ.get("USER") or os.environ.get("LOGNAME") or ""


def memory_gb(text):
    """squeue's memory (`%m`: `8G`, `500M`, `1T`, a number of megabytes) in GB."""
    text = (text or "").strip()
    if not text:
        return 0.0
    units = {"K": 1 / 1024 ** 2, "M": 1 / 1024, "G": 1, "T": 1024}
    if text[-1].upper() in units:
        return float(text[:-1] or 0) * units[text[-1].upper()]
    return float(text) / 1024


def queue():
    """[(name, memory GB)] of the user's jobs in the queue (running and pending); none with the local runner."""
    if CFG.get().runner != "slurm":
        return []
    r = subprocess.run(["squeue", "-u", user(), "-h", "-o", "%j|%m"], capture_output=True, text=True, check=True)
    out = []
    for line in r.stdout.splitlines():
        if "|" in line:
            name, mem = line.rsplit("|", 1)
            out.append((name, memory_gb(mem)))
    return out


def room(mem_gb, planned=()):
    """Why a job asking `mem_gb` would exceed `[jobs.limits]` now (the queue, and `planned` jobs not yet in it), or None."""
    lim = CFG.get().limits
    if lim["max_jobs"] is None and lim["max_memory_gb"] is None:
        return None
    held = [x for x in list(queue()) + list(planned) if x[0].startswith(lim["names"])]
    n, gb = len(held), sum(m for _, m in held)
    if lim["max_jobs"] is not None and n + 1 > lim["max_jobs"]:
        return f"{n} jobs in the queue, at most {lim['max_jobs']} ([jobs.limits] max_jobs)"
    if lim["max_memory_gb"] is not None and gb + mem_gb > lim["max_memory_gb"]:
        return f"{gb:.0f} GB asked in the queue, {mem_gb} more would exceed {lim['max_memory_gb']} ([jobs.limits] max_memory_gb)"
    return None


def job_states(ids):
    """{job id: {state, reason, elapsed_s, max_rss_gb, req_mem}} from squeue (jobs in the queue) and sacct (the others)."""
    ids = [str(i) for i in ids if i and not str(i).startswith("local:")]
    out = {}
    if not ids or CFG.get().runner != "slurm":
        return out
    r = subprocess.run(["squeue", "-h", "-j", ",".join(ids), "-o", "%i|%T|%r"], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        jid, state, reason = (line.split("|") + ["", ""])[:3]
        out[jid.strip()] = dict(state=state.strip(), reason=reason.strip())
    r = subprocess.run(["sacct", "-n", "-P", "-j", ",".join(ids), "-o", "JobID,State,ElapsedRaw,MaxRSS,ReqMem"],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        f = line.split("|")
        if len(f) < 5:
            continue
        base = f[0].split(".")[0]
        rec = out.setdefault(base, {})
        if f[0] == base:
            if "reason" not in rec and f[1]:  # squeue's state, where the job is still in the queue, is the current one
                rec["state"] = f[1].split()[0]
            rec["elapsed_s"] = int(f[2]) if f[2].isdigit() else None
            rec["req_mem"] = f[4]
        if f[3]:
            rec["max_rss_gb"] = round(max(rec.get("max_rss_gb") or 0, memory_gb(f[3])), 3)
    return out


def cancel(jid):
    if CFG.get().runner == "slurm" and jid and not str(jid).startswith("local:"):
        subprocess.run(["scancel", str(jid)], check=False)


def submit(step, after_ids, log_dir):
    """Submit one step (the configured runner); returns (job id, attempt record)."""
    t0 = time.time()
    try:
        jid = jobs.submit(f"iterate_{step.tree}_{step.name}", step.command, log_dir, cpus=step.cpus, mem_gb=step.mem_gb,
                          time=step.time, after=after_ids, dependency=step.dependency)
    except SystemExit as e:  # the local runner raises when the command fails
        return None, dict(job=None, submitted=now(), state="FAILED", seconds=round(time.time() - t0, 1), error=str(e))
    rec = dict(job=jid, submitted=now(), state="SUBMITTED")
    if str(jid).startswith("local:"):  # the local runner ran it here and then
        rec.update(state="COMPLETED", seconds=round(time.time() - t0, 1))
    return jid, rec


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def age(when):
    """Seconds since an ISO time (`now()`), or a large number when it cannot be read."""
    try:
        return (datetime.datetime.now() - datetime.datetime.fromisoformat(when)).total_seconds()
    except (TypeError, ValueError):
        return float("inf")


# ---------------------------------------------------------------------------------------------------------------------
# Running the plan


def load_state(tree):
    return load_json(os.path.join(state_dir(tree), STATE))


def save_state(tree, state):
    path = os.path.join(state_dir(tree), STATE)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(state, fh, indent=1, sort_keys=False, default=str)
        fh.write("\n")
    os.replace(path + ".tmp", path)


def status_of(step, rec, live):
    """(status, why) of a step from its outputs, its latest attempt and the scheduler (`live`: job_states)."""
    attempts = (rec or {}).get("attempts") or []
    last = attempts[-1] if attempts else None
    if last and last.get("job") and not str(last["job"]).startswith("local:"):
        info = live.get(str(last["job"]), {})
        state = info.get("state") or last.get("state")
        last.update({k: v for k, v in info.items() if v is not None})
        if state in ACTIVE:
            if "DependencyNeverSatisfied" in (info.get("reason") or ""):
                return "blocked", f"job {last['job']} waits on a predecessor that did not succeed"
            return "active", f"job {last['job']} {state.lower()}"
    else:
        state = last.get("state") if last else None
    reason = step.done()
    if reason:
        return ("done" if last else "skipped"), reason
    if not step.command:
        return "skipped", step.note or "nothing to run"
    if not last:
        return "todo", step.note
    if state in LOST:
        return "lost", f"job {last.get('job')} {state.lower()} (the scheduler's)"
    if state == "COMPLETED" and step.kind == "report":
        return "done", f"written by job {last.get('job')}"
    if state in FAILED or state == "COMPLETED":
        return "failed", (f"job {last.get('job')} {state.lower()}" if state != "COMPLETED"
                          else f"job {last.get('job')} completed without its output") + \
            (f": {last['error']}" if last.get("error") else "")
    if state and state.startswith("CANCELLED"):
        return "cancelled", f"job {last.get('job')} cancelled"
    if state == "SUBMITTED" and age(last.get("submitted")) < UNKNOWN_AFTER:  # not yet in the scheduler's records
        return "active", f"job {last.get('job')} submitted {last.get('submitted')}, not yet in the scheduler's records"
    return "lost", f"job {last.get('job')} unknown to the scheduler"


def run(tree, previous=None, reference=None, organisms=None, tables=None, snapshot=115, retry=(), dry=False,
        families=None, steps=None, out=print, observe=False):
    """Compute every step's status, submit what is due, save the state and write the report; returns the state.
    `observe`: the statuses and the report only, nothing submitted (the report's job)."""
    cfg = CFG.get()
    t = cfg.tree(tree)
    steps = steps if steps is not None else plan(tree, previous, reference, organisms, tables, snapshot, families)
    for s in steps:
        s.tree = tree
    old = load_state(tree) or {}
    state = dict(tree=tree, commit=t["commit"], previous=previous, reference=reference,
                 organisms=list(dict.fromkeys(s.organism for s in steps if s.organism)),
                 created=old.get("created") or now(), updated=now(), steps={}, **({"complete": old["complete"]}
                                                                                    if old.get("complete") else {}))
    recs = old.get("steps") or {}
    live = job_states([a.get("job") for r in recs.values() for a in (r.get("attempts") or [])[-1:]])
    status = {}
    for s in steps:
        rec = dict(recs.get(s.name) or {})
        rec["attempts"] = [dict(a) for a in rec.get("attempts") or []]
        status[s.name] = status_of(s, rec, live)
        recs[s.name] = rec
    retry = set(retry or ())
    log_dir = os.path.join(state_dir(tree), "logs")
    planned, submitted_any, ran = [], False, set()  # ran: the steps (re)run by this call (or that would be)
    for s in steps:
        st, why = status[s.name]
        if observe:
            continue
        due = st in ("todo", "lost", "blocked") or (st in ("failed", "cancelled") and
                                                     (s.name in retry or "failed" in retry))
        rerun = [x for x in s.after if x in ran]
        if st in ("done", "skipped", "cancelled") and rerun and s.dependency == "afterok" and s.command:
            due, why = True, f"again: {', '.join(rerun)} runs again"  # its output would be stale after them (or the
            # scheduler cancelled it when its predecessor ended badly)
        if s.kind == "report":  # a job after the others where they run as jobs; else written below, by this call
            due = st != "active" and bool(rerun) and cfg.runner == "slurm"
        if not due:
            continue
        waits = [x for x in s.after if status[x][0] in ("failed", "cancelled", "held", "waits", "blocked", "lost", "todo")
                 and s.dependency == "afterok"]
        if waits:
            if st == "blocked":
                last = recs[s.name]["attempts"][-1]
                if not dry:
                    cancel(last.get("job"))
                last["state"] = "CANCELLED (blocked)"
            status[s.name] = ("waits", f"on {', '.join(waits)}")
            continue
        reason = room(s.mem_gb, planned)
        if reason:
            status[s.name] = ("held", reason)
            continue
        after_ids = [recs[x]["attempts"][-1]["job"] for x in s.after
                     if status[x][0] in ("active", "submitted") and recs[x].get("attempts")]
        ran.add(s.name)
        if dry:
            planned.append((cfg.name_prefix + f"iterate_{tree}_{s.name}", s.mem_gb))
            status[s.name] = ("would submit", why if st != "todo" else (s.note or ""))
            continue
        if st == "blocked":
            cancel(recs[s.name]["attempts"][-1].get("job"))
            recs[s.name]["attempts"][-1]["state"] = "CANCELLED (blocked)"
        jid, attempt = submit(s, after_ids, log_dir)
        recs[s.name]["attempts"].append(attempt)
        submitted_any = True
        if jid is None:
            status[s.name] = ("failed", attempt.get("error", ""))
        elif str(jid).startswith("local:"):
            reason = s.done() if s.command else "ran"
            status[s.name] = ("done", reason) if reason or s.kind == "report" else ("failed", "completed without its output")
        else:
            status[s.name] = ("submitted", f"job {jid}" + (f" after {','.join(map(str, after_ids))}" if after_ids else ""))
    if status.get("report", ("",))[0] not in ("submitted", "would submit") and not dry:
        status["report"] = ("written", "by this call")
    for s in steps:
        st, why = status[s.name]
        state["steps"][s.name] = dict(kind=s.kind, organism=s.organism, command=s.command, after=s.after,
                                      dependency=s.dependency, cpus=s.cpus, mem_gb=s.mem_gb, time=s.time,
                                      status=st, why=why, note=s.note, attempts=recs[s.name].get("attempts") or [])
    print_plan(state, out)
    if dry:
        out("(--dry: nothing written, nothing submitted)")
        return state
    save_state(tree, state)
    write_report(tree, state)  # now, and again by the report's job when the submitted steps end
    out(f"report: {os.path.join(state_dir(tree), REPORT_MD)}" + (" (written again when the jobs end)"
                                                                  if submitted_any and cfg.runner == "slurm" else ""))
    return state


def print_plan(state, out=print):
    out(f"iteration {state['tree']} ({str(state.get('commit'))[:7]}): previous {state.get('previous')}, "
        f"reference {state.get('reference')}")
    w = max(len(n) for n in state["steps"]) if state["steps"] else 4
    for name, s in state["steps"].items():
        job = (s["attempts"][-1].get("job") if s["attempts"] else None) or "-"
        out(f"  {name.ljust(w)}  {s['status']:<13} {str(job):<11} {s['mem_gb']:>3} GB  {s['why'] or ''}")
    if any(s["status"] == "held" for s in state["steps"].values()):
        out("steps held by [jobs.limits]: call the command again when the queue has room")


# ---------------------------------------------------------------------------------------------------------------------
# The report


def read_tsv(path):
    try:
        with open(path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    except OSError:
        return []


def graph_digest_of(tree, organism):
    """The content digest of a tree's graph: its build record, else the registry's `[graphs.<tree>]`, else None."""
    from .graphs import common as G
    rec = load_json(os.path.join(G.out_dir(tree, organism), "build.json")) or {}
    if rec.get("content_digest"):
        return rec["content_digest"]
    try:
        from .check import ROOT
        from .results import registry
        return (registry(ROOT).get("graphs", {}).get(tree) or {}).get(organism)
    except (OSError, ValueError):
        return None


def graph_rows(tree, previous, organisms):
    from .graphs import common as G
    rows = []
    for o in organisms:
        b = load_json(os.path.join(G.out_dir(tree, o), "build.json")) or {}
        j = load_json(os.path.join(G.out_dir(tree, o), "judge.json")) or {}
        digest, before = b.get("content_digest"), graph_digest_of(previous, o) if previous else None
        h, d = j.get("history") or {}, j.get("diff") or {}
        links = {k: sum(v.get(k, 0) for v in (d.get("links") or {}).values()) for k in ("only_old", "only_new",
                                                                                        "other_releases")}
        nodes = d.get("nodes") or {}
        rows.append(dict(organism=o, digest=digest or "not built", previous_digest=before or "unknown",
                         equal="yes" if digest and before and digest == before else
                         ("unknown" if not digest or not before else "no"),
                         verify=(j.get("verify") or {}).get("status", "not judged"),
                         history=(f"{h.get('held', 0) + h.get('placed_with_the_version_of_the_release', 0)} of "
                                  f"{h.get('events')} events" if h.get("status") == "checked" else h.get("status", "-")),
                         faults=",".join(j.get("faults") or []) or ("none" if j else "-"),
                         nodes=f"+{sum((nodes.get('added') or {}).values())} -{sum((nodes.get('removed') or {}).values())}"
                         if d else "-",
                         links=f"only old {links['only_old']}, only new {links['only_new']}, other releases "
                               f"{links['other_releases']}, moved {sum((d.get('links_moved_between_databases') or {}).values())}"
                         if d else "-",
                         history_edges=f"+{(d.get('history') or {}).get('n_added', 0)} -{(d.get('history') or {}).get('n_removed', 0)}"
                         if d else "-",
                         build=f"{b.get('build_s')} s, {b.get('peak_GB')} GB" if b else "-"))
    return rows


def causes_of(folder):
    out = {}
    for name in ("causes.tsv.gz", "causes_und.tsv.gz"):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            with gzip.open(p, "rt", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                    out[(r["qid"], r["view"])] = f"{r['cause']} ({r['class']})"
    return out


def part_folders(tree, organism):
    """[(part, scored folder)] of a tree's modest run of an organism (family E's two scorings apart)."""
    base = os.path.join(modest_dir(tree), organism)
    out = []
    for d in sorted(glob.glob(os.path.join(base, "*", ""))):
        part = os.path.basename(os.path.normpath(d))
        if part == "logs":
            continue
        subs = sorted(glob.glob(os.path.join(d, "scores_*", "")))
        out += [(f"{part}/{os.path.basename(os.path.normpath(s))}", os.path.normpath(s)) for s in subs] or \
            [(part, os.path.normpath(d))]
    return out


def moves(tree, other, organisms):
    """Per (organism, part): the counts of moved answers against `other`, and every answer from right (or the error
    the docs promise) to wrong or an error, with its cause or `unexplained`."""
    counts, bad = {}, []
    for o in organisms:
        for part, folder in part_folders(tree, o):
            rows = read_tsv(os.path.join(folder, f"moved_vs_{other}.tsv"))
            if not rows and not os.path.exists(os.path.join(folder, f"moved_vs_{other}.tsv")):
                continue
            c = dict(wrong=0, error=0, undetermined=0, gains=0)
            cz = None
            for r in rows:
                a, b = r.get("outcome_a"), r.get("outcome_b")
                if a in WAS_RIGHT and b in ("wrong", "error", "undetermined"):
                    c[b] += 1
                if b in WAS_RIGHT and a not in WAS_RIGHT:
                    c["gains"] += 1
                if a in WAS_RIGHT and b in NOW_BAD:
                    cz = causes_of(folder) if cz is None else cz
                    bad.append(dict(organism=o, part=part, against=other, qid=r.get("qid"), view=r.get("view"),
                                    move=f"{a} -> {b}", stratum=r.get("stratum"),
                                    cause=cz.get((r.get("qid"), r.get("view")), "unexplained"),
                                    detail=f"got {r.get('got_a') or '-'} -> {r.get('got_b') or '-'}; {r.get('reason_b') or ''}"))
            counts[(o, part)] = c
    return counts, bad


def hunts_rows(tree, organisms):
    rows = []
    for o in organisms:
        d = load_json(os.path.join(modest_dir(tree), o, "hunts.json"))
        if d:
            res = d.get("results") or []
            changed = [f.get("id") for f in res if not f.get("holds")]
            rows.append(dict(organism=o, findings=len(res), holds=len(res) - len(changed), changed=",".join(changed)))
    return rows


def corpus_rows(tree, organisms):
    from .corpus import _common as C
    rows = []
    for o in organisms:
        d = C.run_dir(tree, o)
        ref = load_json(os.path.join(d, "reference.json"))
        if not ref:
            continue
        lines = [f"determinism {ref.get('determinism')}"]
        lines += [f"against {k[len('compare_'):]}: {v}" for k, v in ref.items() if k.startswith("compare_")]
        for p in sorted(glob.glob(os.path.join(d, "compare_reference_vs_*.txt"))):
            with open(p) as fh:
                first = fh.readline().strip()
            lines.append(f"against {os.path.basename(p)[len('compare_reference_vs_'):-4]}: {first}")
        rows.append(dict(organism=o, lines=lines))
    return rows


def job_rows(state):
    live = job_states([a.get("job") for s in state["steps"].values() for a in s["attempts"]])
    rows = []
    for name, s in state["steps"].items():
        for a in s["attempts"]:
            info = live.get(str(a.get("job")), {})
            rows.append(dict(step=name, job=a.get("job") or "-", submitted=a.get("submitted"),
                             state=info.get("state") or a.get("state"),
                             elapsed=info.get("elapsed_s") if info.get("elapsed_s") is not None else a.get("seconds"),
                             max_rss_gb=info.get("max_rss_gb"), asked_gb=s["mem_gb"]))
    return rows


def write_report(tree, state):
    """`iterate_report.md` and `iterate_report.tsv` in the iteration's folder, from what the steps left."""
    previous, reference, organisms = state.get("previous"), state.get("reference"), state.get("organisms") or []
    against = [x for x in dict.fromkeys([previous, reference]) if x and x != tree]
    tsv, md = [], []
    md += [f"# Iteration {tree}", "",
           f"Tree `{tree}` at `{str(state.get('commit'))[:7]}`; previous iteration `{previous or '-'}`, reference "
           f"`{reference or '-'}`; organisms {', '.join(organisms)}. Written {now()} from the steps' outputs "
           f"(`{STATE}` beside it holds the plan and its jobs).", "", "## Steps", "",
           "| step | status | why |", "| --- | --- | --- |"]
    for name, s in state["steps"].items():
        md.append(f"| {name} | {s['status']} | {(s.get('why') or '').replace('|', '/')} |")
        tsv.append(dict(section="step", organism=s.get("organism") or "", key=name, value=s["status"],
                        detail=s.get("why") or ""))
    md += ["", f"## Graphs (against `{previous or '-'}`)", "",
           "| organism | content digest | equal to previous | verify | history | faults | nodes | links | history edges | build |",
           "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in graph_rows(tree, previous, organisms):
        md.append(f"| {r['organism']} | `{r['digest'][:12]}` | {r['equal']} | {r['verify']} | {r['history']} | "
                  f"{r['faults']} | {r['nodes']} | {r['links']} | {r['history_edges']} | {r['build']} |")
        for k in ("digest", "previous_digest", "equal", "verify", "history", "faults", "nodes", "links", "history_edges",
                  "build"):
            tsv.append(dict(section="graph", organism=r["organism"], against=previous or "", key=k, value=r[k]))
    md += ["", "## Parts", "", "Per part of the modest run (or the family's panel): the verdict of `modest compare` and the "
           "answers that moved from right (or the error the docs promise) to wrong, to an error or to undetermined, "
           "and the gains (to right).", ""]
    allbad = []
    for other in against:
        verdicts = {(r["organism"], r["part"]): r["verdict"]
                    for r in read_tsv(os.path.join(modest_dir(tree), f"compare_vs_{other}.tsv"))}
        counts, bad = moves(tree, other, organisms)
        allbad += bad
        md += [f"Against `{other}`:", "", "| organism | part | verdict | wrong · error · undetermined | gains |",
               "| --- | --- | --- | ---: | ---: |"]
        for key in sorted(set(verdicts) | set(counts), key=lambda k: (organisms.index(k[0]) if k[0] in organisms else 99, k[1])):
            c = counts.get(key)
            md.append(f"| {key[0]} | {key[1]} | {verdicts.get(key, '-')} | "
                      + (f"{c['wrong']} · {c['error']} · {c['undetermined']} | {c['gains']} |" if c else "- | - |"))
            tsv.append(dict(section="part", organism=key[0], part=key[1], against=other, key="verdict",
                            value=verdicts.get(key, "")))
            for k, v in (c or {}).items():
                tsv.append(dict(section="part", organism=key[0], part=key[1], against=other,
                                key=k if k == "gains" else f"right_to_{k}", value=v))
        if not verdicts and not counts:
            md.append("| - | - | not compared | - | - |")
        md.append("")
    md += ["## Right to wrong or error", ""]
    if not allbad:
        md += ["None.", ""]
    for other in against:
        bad = [b for b in allbad if b["against"] == other]
        if not bad:
            continue
        by = {}
        for b in bad:
            by[b["cause"]] = by.get(b["cause"], 0) + 1
        md += [f"Against `{other}`: {len(bad)} answers; by cause: "
               + "; ".join(f"{c} {n}" for c, n in sorted(by.items(), key=lambda kv: -kv[1])) + ".", "",
               "| organism | part | qid | view | move | cause |", "| --- | --- | --- | --- | --- | --- |"]
        for b in bad[:EXAMPLES]:
            md.append(f"| {b['organism']} | {b['part']} | `{b['qid']}` | {b['view']} | {b['move']} | {b['cause']} |")
        if len(bad) > EXAMPLES:
            md.append(f"| ... | | {len(bad) - EXAMPLES} more in `{REPORT_TSV}` | | | |")
        md.append("")
    for b in allbad:
        tsv.append(dict(section="moved", organism=b["organism"], part=b["part"], against=b["against"],
                        key=f"{b['qid']}|{b['view']}", value=b["move"], cause=b["cause"], detail=b["detail"]))
    md += ["## Hunts", "", "| organism | findings | hold as recorded | changed |", "| --- | ---: | ---: | --- |"]
    for r in hunts_rows(tree, organisms):
        md.append(f"| {r['organism']} | {r['findings']} | {r['holds']} | {r['changed'] or '-'} |")
        tsv.append(dict(section="hunts", organism=r["organism"], key="holds", value=f"{r['holds']} of {r['findings']}",
                        detail=r["changed"]))
    md += ["", "## Regression corpus", ""]
    rows = corpus_rows(tree, organisms)
    md += [f"- {r['organism']}: " + "; ".join(r["lines"]) for r in rows] or ["No reference answers of the tree yet."]
    for r in rows:
        for line in r["lines"]:
            tsv.append(dict(section="corpus", organism=r["organism"], key=line.split(":")[0], value=line))
    md += ["", "## Jobs", "", "| step | job | state | elapsed (s) | peak (GB) | asked (GB) |",
           "| --- | --- | --- | ---: | ---: | ---: |"]
    for r in job_rows(state):
        md.append(f"| {r['step']} | {r['job']} | {r['state']} | {r['elapsed'] if r['elapsed'] is not None else '-'} | "
                  f"{r['max_rss_gb'] if r['max_rss_gb'] is not None else '-'} | {r['asked_gb']} |")
        tsv.append(dict(section="job", key=r["step"], value=r["state"], part=str(r["job"]),
                        detail=f"elapsed_s={r['elapsed']}; max_rss_gb={r['max_rss_gb']}; asked_gb={r['asked_gb']}"))
    if state.get("complete"):
        md += ["", "## The complete run", ""] + complete_md(state["complete"])
        tsv += complete_tsv(state["complete"])
    folder = state_dir(tree)
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, REPORT_MD + ".tmp"), "w") as fh:
        fh.write("\n".join(md).rstrip() + "\n")
    os.replace(os.path.join(folder, REPORT_MD + ".tmp"), os.path.join(folder, REPORT_MD))
    with open(os.path.join(folder, REPORT_TSV + ".tmp"), "w") as fh:
        fh.write("\t".join(REPORT_COLS) + "\n")
        for r in tsv:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in REPORT_COLS) + "\n")
    os.replace(os.path.join(folder, REPORT_TSV + ".tmp"), os.path.join(folder, REPORT_TSV))


# ---------------------------------------------------------------------------------------------------------------------
# The complete run (`--complete`): what the families' run folders of the tree hold


def run_names(runs, tree):
    """The run folders of a tree in a family's `runs/`: `<tree>`, its variants `<tree>_<what>`, F0's `check_<tree>`."""
    if not os.path.isdir(runs):
        return []
    return sorted(n for n in os.listdir(runs) if n in (tree, f"check_{tree}") or n.startswith(f"{tree}_"))


def scored_in(folder):
    """The scored folders of a run folder: itself, or its scorings (`scores_now`, `scores_decided`)."""
    if os.path.exists(os.path.join(folder, "run.json")):
        return [folder]
    return sorted(os.path.dirname(p) for p in glob.glob(os.path.join(folder, "*", "run.json")))


def report_strata(folder):
    """{stratum: (n, missing)} of a scored folder's report.tsv (both views)."""
    out = {}
    for r in read_tsv(os.path.join(folder, "report.tsv")):
        if r.get("view") == "*":
            out[r["stratum"]] = (int(r.get("n") or 0), int(r.get("missing") or 0))
    return out


def published_sets(family, organism):
    """The family's published sets of the organism (revision sets left out: they revise lines of these)."""
    from .families import lines as L
    d = L.query_dir(family, organism)
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*.manifest.json"))):
        m = load_json(p) or {}
        if not m.get("revises"):
            out.append(os.path.basename(p)[: -len(".manifest.json")])
    return out


def complete(tree, reference=None, organisms=None):
    """The record of a complete (exhaustive) run of the tree: per family and organism its scored folders and their
    `run.json`, the published sets no folder asked, and every stratum with lines missing, with why (asked in another
    folder of the tree; missing on the reference as well) or `gap`."""
    from .records import RECORD_VERSION
    cfg = CFG.get()
    organisms = organisms or list(CFG.ORGANISMS)
    commit = cfg.trees.get(tree, {}).get("commit", "")
    out = dict(checked_at=now(), tree=tree, commit=commit, reference=reference, families={})
    fams = sorted(f for f in os.listdir(cfg.families) if os.path.isdir(os.path.join(cfg.families, f, "runs"))) \
        if os.path.isdir(cfg.families) else []
    for fam in fams:
        runs = os.path.join(cfg.families, fam, "runs")
        names = run_names(runs, tree)
        if not names:
            continue
        for o in organisms:
            folders = [s for n in names for s in scored_in(os.path.join(runs, n, o))]
            if not folders:
                continue
            entry = dict(folders=[], sets_not_asked=[], missing={}, problems=[])
            asked, strata = set(), {}
            for f in folders:
                info = load_json(os.path.join(f, "run.json")) or {}
                asked |= set(info.get("sets") or {})
                strata[f] = report_strata(f)
                rel = CFG.in_root(f)
                entry["folders"].append(dict(folder=rel, commit=info.get("commit"), record_version=info.get("record_version"),
                                             n_rows=info.get("n_rows"), n_missing=info.get("n_missing"),
                                             n_excluded=info.get("n_excluded"), scored_at=info.get("scored_at"),
                                             graph_digests=info.get("graph_digests"), sets=len(info.get("sets") or {})))
                if commit and info.get("commit") and not str(info["commit"]).startswith(commit[:7]):
                    entry["problems"].append(f"{rel}: of commit {str(info['commit'])[:7]}, the tree is {commit[:7]}")
                if int(info.get("record_version") or 0) != RECORD_VERSION:
                    entry["problems"].append(f"{rel}: record version {info.get('record_version')}")
                for stratum, n in (info.get("missing") or {}).items():
                    entry["missing"].setdefault(stratum, []).append(dict(folder=rel, n=n))
            entry["sets_not_asked"] = [s for s in published_sets(fam, o) if s not in asked]
            ref_missing, ref_asked = {}, set()  # strata no folder of the reference asked, and those one did
            for n in run_names(runs, reference) if reference else []:
                for f in scored_in(os.path.join(runs, n, o)):
                    for stratum, k in ((load_json(os.path.join(f, "run.json")) or {}).get("missing") or {}).items():
                        ref_missing[stratum] = ref_missing.get(stratum, 0) + k
                    ref_asked |= {x for x, (k, m) in report_strata(f).items() if k and not m}
            ref_missing = {x: k for x, k in ref_missing.items() if x not in ref_asked}
            for stratum, where in entry["missing"].items():
                elsewhere = [CFG.in_root(f) for f, s in strata.items() if s.get(stratum, (0, 1))[0] and
                             not s[stratum][1]]
                if elsewhere:
                    why = f"asked in {', '.join(elsewhere)}"
                elif stratum in ref_missing:
                    why = f"missing on the reference {reference} as well ({ref_missing[stratum]})"
                else:
                    why = "gap"
                entry["missing"][stratum] = dict(where=where, why=why)
            entry["complete"] = not entry["sets_not_asked"] and not entry["problems"] and \
                all(m["why"] != "gap" for m in entry["missing"].values())
            out["families"].setdefault(fam, {})[o] = entry
    a0 = os.path.join(cfg.families, "a", "a0")
    for n in run_names(a0, tree):
        for o in organisms:
            p = os.path.join(a0, n, o, "summary.json")
            if os.path.exists(p):
                out.setdefault("summaries", []).append(dict(part="a0", folder=CFG.in_root(os.path.dirname(p))))
    return out


def complete_md(c):
    md = [f"Checked {c['checked_at']} (`--complete`): the families' run folders of `{c['tree']}`.", "",
          "| family | organism | complete | scored folders | sets not asked | strata with missing lines |",
          "| --- | --- | --- | --- | --- | --- |"]
    for fam, orgs in c["families"].items():
        for o, e in orgs.items():
            miss = "; ".join(f"{s}: {m['why']}" for s, m in e["missing"].items()) or "-"
            md.append(f"| {fam} | {o} | {'yes' if e['complete'] else 'no'} | "
                      + ", ".join(f"`{f['folder']}`" for f in e["folders"]) + f" | {', '.join(e['sets_not_asked']) or '-'} | "
                      + miss.replace("|", "/") + (" | problems: " + "; ".join(e["problems"]) if e["problems"] else "") + " |")
    for s in c.get("summaries", []):
        md.append(f"| {s['part']} | | kept as a summary | `{s['folder']}` | | |")
    if not c["families"]:
        md.append("| - | - | no run folders of the tree | | | |")
    return md


def complete_tsv(c):
    rows = []
    for fam, orgs in c["families"].items():
        for o, e in orgs.items():
            rows.append(dict(section="complete", organism=o, part=fam, key="complete", value="yes" if e["complete"] else "no",
                             detail="; ".join(e["problems"])))
            for f in e["folders"]:
                rows.append(dict(section="complete", organism=o, part=fam, key="folder", value=f["folder"],
                                 detail=json.dumps({k: v for k, v in f.items() if k != "folder"}, sort_keys=True)))
            for s in e["sets_not_asked"]:
                rows.append(dict(section="complete", organism=o, part=fam, key="set_not_asked", value=s))
            for s, m in e["missing"].items():
                rows.append(dict(section="complete", organism=o, part=fam, key="missing", value=s, cause=m["why"],
                                 detail=json.dumps(m["where"])))
    return rows


def cmd_complete(tree, reference, organisms, dry, out=print):
    c = complete(tree, reference, organisms)
    for line in complete_md(c):
        out(line)
    if dry:
        out("(--dry: nothing written)")
        return c
    state = load_state(tree) or dict(tree=tree, commit=CFG.get().trees.get(tree, {}).get("commit"), steps={},
                                     created=now(), organisms=organisms or list(CFG.ORGANISMS), reference=reference,
                                     previous=CFG.get().previous(tree))
    state["complete"], state["updated"] = c, now()
    save_state(tree, state)
    write_report(tree, state)
    return c


# ---------------------------------------------------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro iterate", description=__doc__.split("\n\n")[0])
    ap.add_argument("tree")
    ap.add_argument("--previous", default=None, help="the iteration before (default: the tree before it in "
                    "[iterations] order)")
    ap.add_argument("--reference", default=None, help="the reference tree (default: [iterations] reference)")
    ap.add_argument("--organisms", default="all", help="comma-separated, or `all` (the tree's configured graphs)")
    ap.add_argument("--families", default=None, help="the modest run's families (default: all)")
    ap.add_argument("--tables", default=None, help="graph-build's stored tables, with `{organism}` (default: its own)")
    ap.add_argument("--snapshot", type=int, default=115)
    ap.add_argument("--retry", default="", help="steps submitted again although they failed or were cancelled "
                    "(comma-separated names, or `failed` for every one)")
    ap.add_argument("--dry", action="store_true", help="print the plan; write and submit nothing")
    ap.add_argument("--report", action="store_true", help="write the report again from the steps' outputs (no job)")
    ap.add_argument("--complete", action="store_true", help="record a complete run: the families' run folders of the "
                    "tree (submits nothing)")
    a = ap.parse_args(argv)
    cfg = CFG.get()
    cfg.tree(a.tree)
    previous = a.previous or cfg.previous(a.tree)
    if previous:
        cfg.tree(previous)
    reference = a.reference or cfg.iteration["reference"]
    organisms = organisms_of(a.tree, a.organisms.split(","))
    if a.complete:
        cmd_complete(a.tree, reference, organisms, a.dry)
        return 0
    if a.report:  # every step's status from its outputs and jobs, then the report; nothing submitted
        if not load_state(a.tree):
            raise SystemExit(f"no plan of {a.tree} yet ({os.path.join(state_dir(a.tree), STATE)})")
        run(a.tree, previous, reference, organisms, a.tables, a.snapshot, observe=True,
            families=a.families.split(",") if a.families else None)
        return 0
    if previous is None and a.tree not in cfg.iteration["order"]:
        raise SystemExit(f"tree {a.tree} is not in [iterations] order {cfg.iteration['order']}: add it, or give "
                         f"--previous <tree>")
    state = run(a.tree, previous, reference, organisms, a.tables, a.snapshot,
                retry=[x for x in a.retry.split(",") if x], dry=a.dry,
                families=a.families.split(",") if a.families else None)
    bad = [n for n, s in state["steps"].items() if s["status"] in ("failed", "cancelled")]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
