"""The modest set: its families' samples (a fixed, stratified selection of qids) and the commands of a run."""
import os

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import modest as M
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import panel as P

ORG = "sus_scrofa"


def publish(sizes):
    lines = [dict(qid=f"{st}|{i}", kind=st, stratum=st, id=f"ENSSSCG{i:011d}", from_release=None, to_release=115,
                  final_database=None, strategy="best", rule="set", cmp="versioned", expect=[], broad=None, locus=None)
             for st, n in sizes.items() for i in range(n)]
    L.publish("x", ORG, "x_v1", lines, "invented", __file__, fmt="jsonl.gz")


def test_a_stratified_sample_of_fixed_size_per_stratum(data_root):
    publish({"small": 3, "mid": 10, "big": 400})
    assert M.stratum_sizes("x", ORG) == {"small": 3, "mid": 10, "big": 400}
    man = P.build("x", ORG, folders=[], per_stratum=5, kind="modest", log=lambda *a: None)
    assert man["name"] == "modest_v1" and man["counts"] == {"big": 5, "mid": 5, "small": 3}
    assert set(man["by_reason"]) == {"sample"}  # never the failure lines of a panel
    first = L.panel_qids("x", ORG, "modest_v1")
    assert P.build("x", ORG, folders=[], per_stratum=5, kind="modest", log=lambda *a: None) is None  # nothing new
    assert M.selection("x", ORG) == "modest_v1"
    lines, _, _ = L.load("x", ORG, sets=["modest_v1"])
    assert {q["qid"] for q in lines} == first and len(lines) == 13


def test_the_sample_does_not_depend_on_the_order_of_reading(tmp_path, monkeypatch):
    picks = []
    for order in ("forward", "backward"):
        root = tmp_path / order
        monkeypatch.setenv("IDTRACK_REPRO_DATA_ROOT", str(root))
        CFG.reset()
        sizes = {"big": 300}
        lines = [dict(qid=f"big|{i}", kind="big", stratum="big", id=f"ENSSSCG{i:011d}", from_release=None,
                      to_release=115, final_database=None, strategy="best", rule="set", cmp="versioned", expect=[])
                 for i in range(sizes["big"])]
        L.publish("x", ORG, "x_v1", lines if order == "forward" else lines[::-1], "invented", __file__, fmt="jsonl.gz")
        P.build("x", ORG, folders=[], per_stratum=7, kind="modest", log=lambda *a: None)
        picks.append(L.panel_qids("x", ORG, "modest_v1"))
    CFG.reset()
    assert picks[0] == picks[1] and len(picks[0]) == 7


def test_a_run_asks_the_harness_the_families_and_the_hunts(tmp_path, monkeypatch, data_root):
    publish({"s": 20})
    P.build("x", ORG, folders=[], per_stratum=5, kind="modest", log=lambda *a: None)
    cfg = tmp_path / "idtrack-repro.toml"
    cfg.write_text(f'[data]\nroot = "{data_root}"\n[trees.t1]\ncommit = "abc"\nsource = "tree"\npython = "py-tree"\n'
                   f'[trees.t1.graphs]\n{ORG} = "g.npz"\n[jobs]\npython = "py-pkg"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(cfg))
    monkeypatch.delenv("IDTRACK_REPRO_DATA_ROOT")
    CFG.reset()
    parts = M.organism_parts("t1", ORG, "/out", families=["x"])
    assert [p for p, _ in parts] == ["harness", "x", "hunts"]
    harness, fam = parts[0][1], parts[1][1]
    assert harness[0].startswith("py-tree -m idtrack_reproducibility.families.runner --family harness")
    assert "--sets modest_v1" in fam[0] and fam[1].startswith("py-pkg -m idtrack_reproducibility.families.scorer --family x")
    assert "--tree t1 --organism sus_scrofa" in parts[2][1][0]
    CFG.reset()


def test_a_part_killed_by_a_signal_stops_the_organisms_run(data_root, monkeypatch, tmp_path):
    import json
    import sys
    import types
    py = sys.executable
    after = tmp_path / "after"
    monkeypatch.setattr(M, "organism_parts", lambda tree, organism, out, families, shards=None: [
        ("harness", [f"{py} -c 'import os, signal; os.kill(os.getpid(), signal.SIGKILL)'", f"{py} -c 'print(1)'"]),
        ("a", [f"{py} -c 'open(r\"{after}\", \"w\")'"])])
    monkeypatch.setattr(M.CFG, "get", lambda: types.SimpleNamespace(tree=lambda name: {"commit": "abc"}, modest_shards={}))
    a = types.SimpleNamespace(tree="t", organism=ORG, out=str(tmp_path / "out"), families=[])
    assert M.cmd_organism(a) == 2
    rec = json.load(open(tmp_path / "out" / ORG / "modest.json"))
    assert rec["killed"]["part"] == "harness" and rec["killed"]["signal"] == 9
    assert len(rec["parts"]) == 1 and not after.exists()  # nothing ran after the killed part
    assert M.killed_by(-9) == 9 and M.killed_by(137) == 9 and M.killed_by(1) is None and M.killed_by(0) is None


def test_a_modest_job_asks_the_measured_peak_of_its_largest_part(data_root, monkeypatch, tmp_path):
    import types
    from idtrack_reproducibility import jobs
    asked = {}
    monkeypatch.setattr(jobs, "submit", lambda name, cmd, logs, cpus, mem_gb, time: asked.setdefault(name, mem_gb) and name)
    monkeypatch.setattr(jobs, "module", lambda py, mod, args: "")
    a = types.SimpleNamespace(tree="t", out=str(tmp_path), organisms=["homo_sapiens", "sus_scrofa"], families=M.FAMILIES,
                              mem=None, time="01:00:00", shards=None)
    M.cmd_run(a)
    assert asked == {"modest_t_homo_sapiens": CFG.MODEST_MEMORY_GB["homo_sapiens"],
                     "modest_t_sus_scrofa": CFG.MODEST_MEMORY_GB["sus_scrofa"]}
    # above the scoring of the harness's sets, the largest part (human 6.14 GB measured)
    assert CFG.get().modest_memory_gb["homo_sapiens"] >= 6.14 * 1.25


def test_answers_of_the_same_tree_are_taken_not_asked_again(data_root, monkeypatch):
    import json
    import types
    import pyarrow as pa
    import pyarrow.parquet as pq
    from idtrack_reproducibility.families import common as C
    from idtrack_reproducibility.families import runner as R
    src = data_root / "families" / "x" / "runs" / "ref" / ORG
    (src / "answers").mkdir(parents=True)
    rows = [dict(qid=q, h=0, view=v, outcome="answer", got=[f"G{q}.1"], seconds=0.1, calls=1)
            for q in ("a", "b") for v in ("best", "all")]
    pq.write_table(pa.Table.from_pylist(rows, schema=R.answer_schema()), src / "answers" / "shard_0_of_1.0.parquet")
    json.dump(dict(commit="c0ffee", graph_digests=["dig"], sets={"x_v1": "sha1"}), open(src / "run.json", "w"))
    monkeypatch.setattr(C, "git_commit", lambda clone: ("c0ffee", False))
    out = data_root / "modest" / "ref" / ORG / "x"
    a = types.SimpleNamespace(clone="tree", out=str(out), organism=ORG, family="x", shard="0/1")
    todo = [dict(qid="a", set="x_v1", strategy="both"), dict(qid="b", set="x_v2", strategy="both"),
            dict(qid="c", set="x_v1", strategy="both")]
    (out / "answers").mkdir(parents=True)
    left, taken = R.reuse(a, todo, {"x_v1": "sha1", "x_v2": "sha2"}, {"content_digest": "dig"})
    # a: taken (same commit, digest, set and sha); b: its set is not one the run read; c: no answer there
    assert [q["qid"] for q in left] == ["b", "c"] and taken == {f"families/x/runs/ref/{ORG}": 1}
    assert sorted((r["qid"], r["view"]) for r in pq.read_table(out / "answers").to_pylist()) == [("a", "all"), ("a", "best")]
    assert R.reuse(a, todo, {"x_v1": "sha1", "x_v2": "sha2"}, {"content_digest": "other"})[1] == {}  # another graph
    monkeypatch.setattr(C, "git_commit", lambda clone: ("beef", False))
    assert R.reuse(a, todo, {"x_v1": "sha1", "x_v2": "sha2"}, {"content_digest": "dig"})[1] == {}  # another commit


def test_a_part_is_asked_in_shards_at_once():
    assert M.asked("run", 1) == "run --shard 0/1"
    assert M.asked("run", 3) == ["run --shard 0/3", "run --shard 1/3", "run --shard 2/3"]
