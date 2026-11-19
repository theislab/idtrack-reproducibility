"""D1's measure of rule `bridge_same_place` (bridges where Ensembl's history has no successor): among the D1 lines whose place names
a gene (`expect`) and whose Ensembl history from the source release has no successor (`history` empty), how often
idtrack's answer is the place's gene (a bridge that lands right), nothing, or another gene; by stratum and view, by
lines and by counts. Usage: bridges_d1.py <scored folder> --organism <o> [--sets d1_v1,d1_r1]"""
import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP

ap = argparse.ArgumentParser()
ap.add_argument("scored")
ap.add_argument("--organism", required=True)
ap.add_argument("--sets", default="d1_v1,d1_r1")
a = ap.parse_args()
lines, _, _ = L.load("d", a.organism, a.sets.split(","))
want = {q["qid"]: q for q in lines if q["qid"].startswith("d1|") and q.get("expect") and q.get("history") == [] and q["kind"] == "d1_id"}
df = RP.read_outcomes(a.scored)
df = df[df["qid"].isin(want.keys())]
c = collections.defaultdict(collections.Counter)
for qid, view, o, n_got in zip(df["qid"], df["view"], df["outcome"], df["n_got"]):
    q = want[qid]
    k = "right (bridge to the place's gene)" if o == "right" else ("nothing" if not n_got else "another gene") if o == "wrong" else o
    for key in ((q["stratum"], view), ("all strata", view)):
        c[key][k] += 1
        c[key][k + " counts"] += int(q.get("counts") or 0)
out = os.path.join(a.scored, "d1", "bridges_no_successor.tsv")
keys = ["right (bridge to the place's gene)", "nothing", "another gene"]
with open(out, "w") as fo:
    fo.write("\t".join(["stratum", "view"] + keys + [k + " counts" for k in keys]) + "\n")
    for (st, v), x in sorted(c.items()):
        fo.write("\t".join(map(str, [st, v] + [x[k] for k in keys] + [x[k + " counts"] for k in keys])) + "\n")
print(open(out).read())
