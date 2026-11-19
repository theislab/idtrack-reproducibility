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
