"""The configuration, the job runner's local mode and the checks' patterns."""
import os
import re

import pytest

from idtrack_reproducibility import check as CHK
from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import jobs


def test_defaults_without_a_file(data_root):
    c = CFG.get()
    assert c.truth == os.path.join(str(data_root), "truth")
    assert c.families == os.path.join(str(data_root), "families")
    assert c.runner == "local"


def test_file_paths_are_read_from_the_files_folder(tmp_path, monkeypatch):
    for k in list(os.environ):
        if k.startswith("IDTRACK_REPRO_"):
            monkeypatch.delenv(k)
    cfg = tmp_path / "conf" / "idtrack-repro.toml"
    cfg.parent.mkdir()
    cfg.write_text('[data]\nroot = "../store"\ntruth = "/elsewhere/truth"\n[trees.t1]\ncommit = "abc"\nsource = "tree"\n'
                   '[trees.t1.graphs]\nsus_scrofa = "g/pig.npz"\n[jobs]\nrunner = "slurm"\n[jobs.slurm]\nqos = "q1"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(cfg))
    monkeypatch.setenv("IDTRACK_REPRO_SLURM_QOS", "q2")
    CFG.reset()
    c = CFG.get()
    assert c.root == str(tmp_path / "store")
    assert c.harness == str(tmp_path / "store" / "harness")
    assert c.truth == "/elsewhere/truth"
    assert c.tree("t1")["source"] == str(tmp_path / "conf" / "tree")
    assert c.graph("t1", "sus_scrofa") == str(tmp_path / "conf" / "g" / "pig.npz")
    assert c.slurm["qos"] == "q2"  # the variable wins over the file
    with pytest.raises(SystemExit):
        c.tree("t2")
    CFG.reset()


def test_local_runner_runs_the_command_and_keeps_its_log(data_root, tmp_path):
    jid = jobs.submit("hello", "echo hello-$PYTHONHASHSEED", str(tmp_path / "logs"))
    assert jid.startswith("local:")
    logs = os.listdir(tmp_path / "logs")
    assert len(logs) == 1 and open(tmp_path / "logs" / logs[0]).read().strip() == "hello-0"
    with pytest.raises(SystemExit):
        jobs.submit("fails", "exit 3", str(tmp_path / "logs"))


# the examples are written in pieces, so that this file does not hold the text it tests
@pytest.mark.parametrize("text,flagged", [
    ("data in /" + "home/someone/work", "machine path"),
    ("see /" + "lustre/groups/x", "machine path"),
    ("written by an AI " + "ag" + "ent working on it", "production mention"),
    ("as deci" + "sion 71 says", "private decision number"),
    ("finding accuracy" + "_build-4", "private review finding"),
    ("built in S" + "2b", "private plan stage"),
    ("rule `support_ranks_shared_name`; hunt finding IN-3; the data root's truth/", None),
])
def test_forbidden_text(text, flagged):
    hits = [n for n, p in CHK.FORBIDDEN if re.search(p, text)]
    assert hits == ([flagged] if flagged else [])


def test_anchors_follow_github(tmp_path):
    p = tmp_path / "x.md"
    p.write_text("# The accuracy harness\n## The corrections of the judge\n## `score.py` and A*: runs\n## Same\n## Same\n"
                 "```\n# not a heading\n```\n")
    a = CHK.anchors_of(str(p))
    assert {"the-accuracy-harness", "the-corrections-of-the-judge", "scorepy-and-a-runs", "same", "same-1"} <= a
    assert "not-a-heading" not in a


def test_no_tracked_file_shows_a_machine_path_the_account_or_a_work_order_name():
    """Code, docs, configurations, results and manifests (their compressed ones too) name no path of a machine, no
    account and no agent or work order: the cluster's paths live only in the data root's own configuration."""
    from idtrack_reproducibility import check as K
    problems = []
    K.check_text(K.tracked_files(), problems)
    bad = [p for p in problems if p[2].split(":")[0] in ("machine path", "cluster account or work folder",
                                                           "agent or work-order name", "production mention")]
    assert not bad, bad[:5]
    for name, text in (("machine path", "/" + "home/someone/x"), ("machine path", "/" + "lustre/groups/x"),
                       ("cluster account or work folder", "idtrack_upgrade_" + "work"),
                       ("agent or work-order name", "cleanup_" + "repro")):
        assert any(n == name and __import__("re").search(p, text) for n, p in K.FORBIDDEN), (name, text)
