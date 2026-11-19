"""The panel command of the common machinery on a small invented family (`zz`) in an invented data root.

The family: a base set (strata s1 of 300 lines, s2 of 40), a revision set (3 lines with a new `expect`, 1 excluded),
a hunt set `hunt_x` (2 lines); a scored folder of one tree (4 wrong, 2 undetermined, 1 error, the rest right).
Checks: the panel holds every line not right, the hunt lines and per stratum the sample (50 of s1, all 40 of s2);
`lines.load(sets=[panel])` gives the latest lines of those qids (the revision applies, the excluded line is left
out), with a shard; the sample is the same when built again; nothing new -> no new version; a second tree with a
new wrong line -> panel_v2, a superset; a panel is never rewritten.
"""
import json
import os

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import panel as P

F, O = "zz", "sus_scrofa"
ok = []


def check(name, cond):
    ok.append((name, bool(cond)))
    print("ok " if cond else "BAD", name)


def line(i, st, expect="G1"):
    return dict(qid=f"zz|{st}|{i}", kind="k", stratum=st, id=f"X{i}", from_release=None, to_release=115, final_database=None,
                strategy="best", rule="set", cmp="versioned", expect=[expect], broad=None, unit=f"u{i % 37}")


def scored(tree, wrong, und, err, base):
    import pyarrow as pa
    import pyarrow.parquet as pq
    d = os.path.join(C.BENCH, F, "runs", tree, O)
    os.makedirs(os.path.join(d, "outcomes"), exist_ok=True)
    rows = []
    for q in [x["qid"] for x in base]:
        o = "wrong" if q in wrong else "undetermined" if q in und else "error" if q in err else "right"
        rows.append(dict(qid=q, view="best", stratum=q.split("|")[1], outcome=o))
    pq.write_table(pa.Table.from_pylist(rows), os.path.join(d, "outcomes", "part_0_of_1.parquet"))
    json.dump(dict(organism=O, family=F, parts=1, commit="test", sets={"v1": "x"}), open(os.path.join(d, "run.json"), "w"))
    return d


def test_panel(data_root):
    ok.clear()
    base = [line(i, "s1") for i in range(300)] + [line(i, "s2") for i in range(40)]
    L.publish(F, O, "v1", base, "invented base", __file__)
    rev = [dict(line(i, "s1", "G2"), revised_by=["t"]) for i in (1, 2, 3)] + [dict(line(4, "s1"), excluded="test")]
    L.publish(F, O, "v1r1", rev, "invented revision", __file__, revises=["v1"])
    L.publish(F, O, "hunt_x", [line(i, "hx") for i in range(2)], "invented hunt", __file__)


    W1 = {"zz|s1|10", "zz|s1|11", "zz|s2|5", "zz|s1|2"}
    U1, E1 = {"zz|s1|20", "zz|s2|6"}, {"zz|s1|30"}
    d1 = scored("t1", W1, U1, E1, base)
    m1 = P.build(F, O, per_stratum=50)
    q1 = L.panel_qids(F, O, "panel_v1")
    check("every line not right is in the panel", (W1 | U1 | E1) <= q1)
    check("the hunt lines are in the panel", {"zz|hx|0", "zz|hx|1"} <= q1)
    check("all 40 lines of the small stratum s2", {f"zz|s2|{i}" for i in range(40)} <= q1)
    s1 = {q for q in q1 if q.startswith("zz|s1|")}
    check("about the sample of s1: 50 sampled plus the failures", 50 <= len(s1) <= 50 + len(W1 | U1 | E1))
    check("manifest counts by reason", m1["by_reason"]["hunt"] == 2 and m1["by_reason"]["failure"] == 7)
    got, shas, excl = L.load(F, O, sets=["panel_v1"])
    by = {x["qid"]: x for x in got}
    check("load gives only the panel's qids", set(by) <= q1)
    check("the revision applies to a panel line", by["zz|s1|2"]["expect"] == ["G2"])
    check("an excluded line is left out", "zz|s1|4" not in by)
    check("the panel's sha256 is recorded", "panel_v1" in shas)
    sh = L.load(F, O, sets=["panel_v1"], shard=(0, 2))[0]
    check("a shard of the panel", 0 < len(sh) < len(got) and all(C.shard_of(x["qid"], 2) == 0 for x in sh))
    only_v1 = {x["qid"]: x for x in L.load(F, O, sets=["panel_v1", "v1"])[0]}
    check("the panel over named sets only (the revision not named: the base line)", only_v1["zz|s1|2"]["expect"] == ["G1"])
    check("load_qids reads through a panel name", L.load_qids(F, O, {"zz|s1|2"}, sets=["panel_v1"])["zz|s1|2"]["expect"] == ["G2"])
    sample_again, _, _ = P.scan_sets(F, O, 50, C.SEED, set())
    check("the sample is the same when drawn again", set(sample_again) <= q1)
    check("nothing new: no new version", P.build(F, O, per_stratum=50) is None)
    scored("t2", W1 | {"zz|s1|200"}, U1, E1, base)
    m2 = P.build(F, O, per_stratum=50)
    q2 = L.panel_qids(F, O, "panel_v2")
    check("a new wrong line on a second tree: panel_v2, a superset", m2 and m2["name"] == "panel_v2" and q1 <= q2 and "zz|s1|200" in q2)
    try:
        import pyarrow.parquet as pq  # noqa: F401
        os.chmod(os.path.join(L.panel_dir(F, O), "panel_v1.parquet"), 0o400)
        open(os.path.join(L.panel_dir(F, O), "panel_v1.parquet"), "a").close()
        check("a published panel is read-only", False)
    except PermissionError:
        check("a published panel is read-only", True)
    bad = [n for n, c in ok if not c]
    assert not bad, bad
