"""One command per iteration (docs/iterate.md): the plan's order and dependencies, skipping what is done, resuming after
a lost step, the limits, the report from fixture outputs, `--dry`, and the record of a complete run. Invented
configurations, stub steps, the `local` runner or a stand-in scheduler; no data root, no idtrack."""
import gzip
import json
import os
import subprocess
import sys

import pytest

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import iterate as IT
from idtrack_reproducibility import jobs

ORG = "sus_scrofa"
GRAPH = "graph_sus_scrofa_assembly111_min56_max115.npz"


def quiet(*_):
    pass


@pytest.fixture
def conf(data_root, tmp_path, monkeypatch):
    """A configuration of three trees (t0 the reference, t1, t2) of one organism, the order t0 -> t1 -> t2, the local
    runner; returns the data root."""
    lines = ["[iterations]", 'order = ["t0", "t1", "t2"]', "", "[jobs]", 'runner = "local"', f'python = "{sys.executable}"', ""]
    for t in ("t0", "t1", "t2"):
        lines += [f"[trees.{t}]", f'commit = "{t}aaaaaaa"', f'source = "trees/{t}"', f'python = "{sys.executable}"',
                  f"[trees.{t}.graphs]", f'{ORG} = "{data_root}/graphs/{t}/{ORG}/{GRAPH}"', ""]
    p = tmp_path / "idtrack-repro.toml"
    p.write_text("\n".join(lines))
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(p))
    CFG.reset()
    return data_root


def built(root, tree, digest="d1"):
    d = root / "graphs" / tree / ORG
    d.mkdir(parents=True, exist_ok=True)
    (d / GRAPH).write_bytes(b"x")
    (d / "build.json").write_text(json.dumps(dict(content_digest=digest, graph_file=f"{tree}/{ORG}/{GRAPH}", build_s=1,
                                                  peak_GB=0.1)))


def stub(name, root, after=(), kind="stub", dependency="afterok", fail=False, mem=2):
    """A step whose job writes `<root>/<name>.out` (or fails), done when that file exists."""
    out = root / f"{name}.out"
    code = "import sys; sys.exit(3)" if fail else f"open({str(out)!r}, 'w').write('x')"
    return IT.Step(name, kind, f"{sys.executable} -c {IT.q(code)}", ORG, after=list(after), mem_gb=mem,
                   dependency=dependency, done=lambda: "written" if out.exists() else None)


class Scheduler:
    """A stand-in for SLURM: jobs.submit records each job and puts it in the queue; states are set by the test."""

    def __init__(self, monkeypatch, foreign=()):
        self.queue, self.submitted, self.cancelled, self.states, self.n = list(foreign), [], [], {}, 100
        CFG.get().runner = "slurm"
        monkeypatch.setattr(jobs, "submit", self.submit)
        monkeypatch.setattr(IT, "queue", lambda: list(self.queue))
        monkeypatch.setattr(IT, "job_states", lambda ids: {i: dict(self.states[i]) for i in map(str, ids) if i in self.states})
        monkeypatch.setattr(IT, "cancel", self.cancelled.append)

    def submit(self, what, command, log_dir, cpus=2, mem_gb=8, time="04:00:00", after=None, env=None, dependency="afterany"):
        self.n += 1
        jid = str(self.n)
        self.submitted.append(dict(name=what, after=list(after or []), dependency=dependency, mem=mem_gb, id=jid))
        self.queue.append((CFG.get().name_prefix + what, mem_gb))
        self.states[jid] = dict(state="PENDING", reason="Dependency" if after else "Priority")
        return jid


def test_the_plan_builds_judges_runs_compares_and_reports_in_order(conf):
    (conf / "modest" / "t1").mkdir(parents=True)
    (conf / "modest" / "t0").mkdir(parents=True)
    (conf / "corpus" / "corpus").mkdir(parents=True)
    (conf / "corpus" / "corpus" / f"{ORG}.jsonl").write_text("")
    (conf / "corpus" / "reference").mkdir()
    (conf / "corpus" / "reference" / f"{ORG}.t0.jsonl.gz").write_bytes(b"")
    steps = IT.plan("t2", previous="t1", reference="t0", organisms=[ORG])
    names = [s.name for s in steps]
    assert names == [f"graph_build_{ORG}", f"graph_judge_{ORG}", f"modest_{ORG}", "compare_previous", "compare_reference",
                     f"corpus_{ORG}", "report"]
    by = {s.name: s for s in steps}
    assert by[f"graph_judge_{ORG}"].after == [f"graph_build_{ORG}"] and "--against t1" in by[f"graph_judge_{ORG}"].command
    assert by[f"modest_{ORG}"].after == [f"graph_build_{ORG}"] and "modest organism --tree t2" in by[f"modest_{ORG}"].command
    assert by["compare_previous"].after == [f"modest_{ORG}"] and "compare_vs_t1.tsv" in by["compare_previous"].command
    assert "compare_vs_t0.tsv" in by["compare_reference"].command
    # the corpus: compared with the nearest earlier tree that has reference answers (t1 has none: t0)
    assert "reference --tree t2" in by[f"corpus_{ORG}"].command and "--compare-with t0" in by[f"corpus_{ORG}"].command
    assert by["report"].after == names[:-1] and by["report"].dependency == "afterany"
    assert all(s.dependency == "afterok" for s in steps[:-1])
    assert by[f"graph_build_{ORG}"].mem_gb == 3  # graphs.common.MEMORY_GB of pig


def test_a_panel_replaces_the_modest_selection_of_its_family(conf):
    d = conf / "families" / "a" / "panels" / ORG
    d.mkdir(parents=True)
    for n in ("modest_v1", "panel_v1", "panel_v2"):
        (d / f"{n}.manifest.json").write_text("{}")
    m = {s.name: s for s in IT.plan("t2", previous="t1", reference="t0", organisms=[ORG])}[f"modest_{ORG}"]
    assert "--selections a=panel_v2" in m.command and "panel_v2" in m.note


def test_a_done_step_is_skipped(conf):
    built(conf, "t2")
    by = {s.name: s for s in IT.plan("t2", previous="t1", reference="t0", organisms=[ORG])}
    assert IT.status_of(by[f"graph_build_{ORG}"], {}, {})[0] == "skipped"
    assert IT.status_of(by[f"modest_{ORG}"], {}, {})[0] == "todo"
    # a build record of another file does not count
    rec = json.loads((conf / "graphs" / "t2" / ORG / "build.json").read_text())
    rec["graph_file"] = f"t2/{ORG}/other.npz"
    (conf / "graphs" / "t2" / ORG / "build.json").write_text(json.dumps(rec))
    assert IT.status_of(by[f"graph_build_{ORG}"], {}, {})[0] == "todo"


def test_local_runner_runs_the_steps_in_order_and_a_second_call_does_nothing(conf, tmp_path):
    steps = [stub("a", tmp_path), stub("b", tmp_path, ["a"]), stub("report", tmp_path, ["a", "b"], kind="report",
                                                                     dependency="afterany")]
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert state["steps"]["a"]["status"] == "done" and state["steps"]["b"]["status"] == "done"
    assert (conf / "iterations" / "t2" / "iterate.json").exists()
    assert (conf / "iterations" / "t2" / "iterate_report.md").exists()
    saved = json.loads((conf / "iterations" / "t2" / "iterate.json").read_text())
    assert [x["state"] for x in saved["steps"]["a"]["attempts"]] == ["COMPLETED"]
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert len(state["steps"]["a"]["attempts"]) == 1 and state["steps"]["b"]["status"] == "done"


def test_a_failed_step_stops_its_dependents_and_is_run_again_only_when_asked(conf, tmp_path):
    steps = [stub("a", tmp_path, fail=True), stub("b", tmp_path, ["a"])]
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert state["steps"]["a"]["status"] == "failed" and state["steps"]["b"]["status"] == "waits"
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert len(state["steps"]["a"]["attempts"]) == 1  # nothing resubmits a failure by itself
    steps[0] = stub("a", tmp_path)  # the fault fixed
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, retry=["a"], out=quiet)
    assert state["steps"]["a"]["status"] == "done" and state["steps"]["b"]["status"] == "done"
    assert len(state["steps"]["a"]["attempts"]) == 2


def test_a_lost_step_is_submitted_again_under_its_name_with_its_blocked_dependents(conf, tmp_path, monkeypatch):
    sched = Scheduler(monkeypatch)
    steps = [stub("a", tmp_path), stub("b", tmp_path, ["a"]), stub("c", tmp_path)]
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert [s["name"] for s in sched.submitted] == ["iterate_t2_a", "iterate_t2_b", "iterate_t2_c"]
    assert sched.submitted[1]["after"] == ["101"] and sched.submitted[1]["dependency"] == "afterok"
    assert state["steps"]["b"]["status"] == "submitted"
    # a is preempted, b waits for ever on it, c completes
    sched.states = {"101": dict(state="PREEMPTED"), "102": dict(state="PENDING", reason="DependencyNeverSatisfied"),
                    "103": dict(state="COMPLETED", elapsed_s=5, max_rss_gb=0.2)}
    sched.queue = [x for x in sched.queue if not x[0].endswith("_c")]
    (tmp_path / "c.out").write_text("x")
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert [s["name"] for s in sched.submitted[3:]] == ["iterate_t2_a", "iterate_t2_b"]
    assert sched.cancelled == ["102"] and sched.submitted[4]["after"] == ["104"]
    assert state["steps"]["c"]["status"] == "done" and len(state["steps"]["c"]["attempts"]) == 1
    assert [x["job"] for x in state["steps"]["a"]["attempts"]] == ["101", "104"]
    assert state["steps"]["a"]["attempts"][0]["state"] == "PREEMPTED"


def test_the_limits_hold_jobs_and_their_dependents(conf, tmp_path, monkeypatch):
    CFG.get().limits.update(max_jobs=3, max_memory_gb=100, names="idtrack_")
    sched = Scheduler(monkeypatch, foreign=[("idtrack_other", 10), ("someone_else", 500)])
    steps = [stub("a", tmp_path), stub("b", tmp_path), stub("c", tmp_path), stub("d", tmp_path, ["c"])]
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    # one job counted already (the other's name does not match), so two more fit
    assert [s["name"] for s in sched.submitted] == ["iterate_t2_a", "iterate_t2_b"]
    assert state["steps"]["c"]["status"] == "held" and "max_jobs" in state["steps"]["c"]["why"]
    assert state["steps"]["d"]["status"] == "waits"
    CFG.get().limits.update(max_jobs=None, max_memory_gb=15)
    sched.queue = [("idtrack_other", 10)]
    steps = [stub("e", tmp_path, mem=4), stub("f", tmp_path, mem=4)]
    state = IT.run("t1", "t0", "t0", [ORG], steps=steps, out=quiet)
    assert state["steps"]["e"]["status"] == "submitted" and state["steps"]["f"]["status"] == "held"
    assert "max_memory_gb" in state["steps"]["f"]["why"]


def test_dry_writes_and_submits_nothing(conf, tmp_path, monkeypatch, capsys):
    sched = Scheduler(monkeypatch)
    assert IT.main(["t2", "--dry"]) == 0
    printed = capsys.readouterr().out
    assert f"graph_build_{ORG}" in printed and "would submit" in printed and "nothing written" in printed
    assert sched.submitted == [] and not (conf / "iterations").exists()
    steps = [stub("a", tmp_path)]
    IT.run("t2", "t1", "t0", [ORG], steps=steps, dry=True, out=quiet)
    assert sched.submitted == [] and not (conf / "iterations").exists() and not (tmp_path / "a.out").exists()


def test_previous_comes_from_the_order(conf):
    cfg = CFG.get()
    assert cfg.previous("t2") == "t1" and cfg.previous("t0") is None and cfg.iteration["reference"] == "t0"


def test_a_verdict_is_no_failure_but_a_crash_is(tmp_path):
    out = tmp_path / "verdict.tsv"
    run = lambda c: subprocess.run(IT.verdict_command(c, out), shell=True).returncode  # noqa: E731
    assert run("true") == 0
    assert run(f"echo x > {out}; exit 1") == 0      # exit 1 having written its verdict: a result
    os.utime(out, (1_000_000_000, 1_000_000_000))    # the verdict of an earlier job
    assert run("exit 1") != 0                         # exit 1 with the verdict from before this job: a failure
    assert run(f"echo x > {out}; exit 2") != 0


def test_squeue_memory():
    assert IT.memory_gb("8G") == 8 and IT.memory_gb("512M") == 0.5 and IT.memory_gb("1T") == 1024
    assert IT.memory_gb("2048") == 2 and IT.memory_gb("") == 0


def write_gz(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fo:
        fo.write(text)


MOVED_HEAD = "qid\tview\tkind\tstratum\tvariant\tmove\toutcome_a\trule_a\toutcome_b\trule_b\tgot_a\tgot_b\treason_a\treason_b\n"


def test_the_report_from_the_steps_outputs(conf):
    built(conf, "t1", "same")
    built(conf, "t2", "same")
    g = conf / "graphs" / "t2" / ORG
    (g / "judge.json").write_text(json.dumps(dict(
        graph=str(g / GRAPH), holds=True, faults=[], verify=dict(status="holds"),
        history=dict(status="checked", events=10, held=9, placed_with_the_version_of_the_release=1),
        diff=dict(against="t1", nodes=dict(added={"external": 2}, removed={}), links={"RefSeq": {"only_new": 3}},
                  links_moved_between_databases={}, history=dict(n_added=0, n_removed=1)))))
    m = conf / "modest" / "t2"
    part = m / ORG / "c"
    part.mkdir(parents=True)
    (m / "compare_vs_t1.tsv").write_text(f"organism\tpart\tverdict\n{ORG}\tc\tregression\n")
    (part / "moved_vs_t1.tsv").write_text(
        MOVED_HEAD
        + "q1\tbest\tk\tc1\tbest,from\tright -> wrong\tright\tpublished\twrong\tpublished\tG1\tG2\t\tnot the gene\n"
        + "q2\tall\tk\tc1\tall,from\texpected_error -> error\texpected_error\tpublished\terror\tpublished\t\t\t\tValueError\n"
        + "q3\tbest\tk\tc1\tbest,from\tright -> undetermined\tright\tpublished\tundetermined\tpublished\tG1\tG3\t\t\n"
        + "q4\tbest\tk\tc1\tbest,from\twrong -> right\twrong\tpublished\tright\tpublished\tG9\tG1\t\t\n")
    write_gz(str(part / "causes.tsv.gz"), "organism\tqid\tview\tkind\tstratum\tunit\tcause\tclass\tdetail\n"
             f"{ORG}\tq1\tbest\tk\tc1\tu\texternal_bridge\tb\t-\n")
    (m / ORG / "hunts.json").write_text(json.dumps(dict(results=[dict(id="IN-1", holds=True), dict(id="IN-2", holds=False)])))
    state = dict(tree="t2", commit="t2aaaaaaa", previous="t1", reference="t0", organisms=[ORG],
                 steps={"report": dict(kind="report", organism=None, status="written", why="", mem_gb=4,
                                       attempts=[dict(job="local:1", state="COMPLETED", seconds=2.0)])})
    IT.write_report("t2", state)
    md = (conf / "iterations" / "t2" / "iterate_report.md").read_text()
    assert "| c | regression | 1 · 1 · 1 | 1 |" in md
    assert "`q1` | best | right -> wrong | external_bridge (b)" in md and "`q2` | all | expected_error -> error | unexplained" in md
    assert "| yes | holds | 10 of 10 events" in md and "| IN-2 |" in md
    rows = [r.split("\t") for r in (conf / "iterations" / "t2" / "iterate_report.tsv").read_text().splitlines()]
    assert rows[0] == IT.REPORT_COLS
    moved = [r for r in rows if r[0] == "moved"]
    assert [(r[4], r[5], r[6]) for r in moved] == [("q1|best", "right -> wrong", "external_bridge (b)"),
                                                    ("q2|all", "expected_error -> error", "unexplained")]
    assert ["graph", ORG, "", "t1", "equal", "yes", "", ""] in rows


def test_complete_lists_the_folders_and_finds_the_gaps(conf):
    fam = conf / "families" / "x"
    q = fam / "queries" / ORG
    q.mkdir(parents=True)
    for s, rev in (("x1_v1", None), ("x2_v1", None), ("x3_v1", None), ("x1_r1", "x1_v1")):
        (q / f"{s}.manifest.json").write_text(json.dumps(dict(revises=rev)))

    def scored(run, sets, missing, strata):
        d = fam / "runs" / run / ORG
        d.mkdir(parents=True)
        (d / "run.json").write_text(json.dumps(dict(record_version=1, commit=run[:2] + "aaaaaaa", sets=sets, n_rows=10,
                                                    n_missing=sum(missing.values()), missing=missing)))
        (d / "report.tsv").write_text("stratum\tview\tn\tmissing\n" + "".join(f"{s}\t*\t{n}\t{k}\n" for s, (n, k) in strata.items()))

    scored("t2", {"x1_v1": "h", "x1_r1": "h"}, {"tx": 4, "odd": 2, "lost": 1},
           {"x1": (8, 0), "tx": (4, 4), "odd": (2, 2), "lost": (1, 1)})
    scored("t2_tx", {"x2_v1": "h"}, {}, {"tx": (4, 0)})
    scored("t0", {"x1_v1": "h"}, {"odd": 2}, {"odd": (2, 2)})
    scored("t1", {"x1_v1": "h"}, {}, {})          # another tree's run: not t2's
    c = IT.complete("t2", reference="t0", organisms=[ORG])
    e = c["families"]["x"][ORG]
    assert [f["folder"] for f in e["folders"]] == ["families/x/runs/t2/sus_scrofa", "families/x/runs/t2_tx/sus_scrofa"]
    assert e["sets_not_asked"] == ["x3_v1"]
    assert e["missing"]["tx"]["why"].startswith("asked in families/x/runs/t2_tx")
    assert e["missing"]["odd"]["why"].startswith("missing on the reference t0")
    assert e["missing"]["lost"]["why"] == "gap" and e["complete"] is False
    IT.cmd_complete("t2", "t0", [ORG], dry=False, out=quiet)
    saved = json.loads((conf / "iterations" / "t2" / "iterate.json").read_text())
    assert saved["complete"]["families"]["x"][ORG]["sets_not_asked"] == ["x3_v1"]
    md = (conf / "iterations" / "t2" / "iterate_report.md").read_text()
    assert "## The complete run" in md and "lost: gap" in md


def test_a_done_step_runs_again_after_a_predecessor_that_runs_again(conf, tmp_path):
    (tmp_path / "b.out").write_text("x")  # b's output is there, a's is not: b would be stale after a
    steps = [stub("a", tmp_path), stub("b", tmp_path, ["a"]), stub("c", tmp_path)]
    (tmp_path / "c.out").write_text("x")
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    assert len(state["steps"]["a"]["attempts"]) == 1 and len(state["steps"]["b"]["attempts"]) == 1
    assert state["steps"]["c"]["status"] == "skipped" and not state["steps"]["c"]["attempts"]


def test_the_reports_job_observes_and_submits_nothing(conf, tmp_path, monkeypatch):
    sched = Scheduler(monkeypatch)
    steps = [stub("a", tmp_path), stub("b", tmp_path, ["a"])]
    IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet)
    n = len(sched.submitted)
    sched.states["101"] = dict(state="COMPLETED")
    (tmp_path / "a.out").write_text("x")
    state = IT.run("t2", "t1", "t0", [ORG], steps=steps, out=quiet, observe=True)
    assert len(sched.submitted) == n and state["steps"]["a"]["status"] == "done"
    assert state["steps"]["b"]["status"] == "active"
    assert (conf / "iterations" / "t2" / "iterate_report.md").exists()
