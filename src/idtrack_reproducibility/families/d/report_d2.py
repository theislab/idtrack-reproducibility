"""D2's own report (beside F0's report.tsv): source inference by list size and noise, how the answers miss, and the
matrices line by line.

Usage: report_d2.py <scored folder> --organism <o> [--sets d2_v1]

Writes `<scored>/d2/`:
  synthetic.tsv  size, noise: lines, right, wrong, rate; of the wrong: the answer holds the expected releases and more
                 (`superset`: idtrack ties releases the GTFs tell apart), a part of them (`subset`), other releases
                 (`other`), with the drawn release among the answer (`drawn_in_answer`)
  matrices.tsv   every matrix line: dataset, source, features, stated release, rosetta's release, the expected releases,
                 idtrack's answer, outcome
The answer of a right line is its expected set (F0 stores the answer of non-right rows only).
"""
import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP


def rng(xs):
    xs = sorted(int(x) for x in xs)
    if not xs:
        return "-"
    return f"{xs[0]}..{xs[-1]} ({len(xs)})" if len(xs) > 2 else ",".join(map(str, xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="d2_v1")
    a = ap.parse_args()
    lines, _, _ = L.load("d", a.organism, a.sets.split(","))
    by = {q["qid"]: q for q in lines if q["qid"].startswith("d2")}
    df = RP.read_outcomes(a.scored)
    df = df[df["qid"].isin(by.keys())]
    out = os.path.join(a.scored, "d2")
    os.makedirs(out, exist_ok=True)
    syn = collections.defaultdict(collections.Counter)
    mats = []
    for qid, o, got in zip(df["qid"], df["outcome"], df["got"]):
        q = by[qid]
        E = set(q["expect"])
        g = E if o == "right" else set((got or "").split(",")) - {""}
        if qid.startswith("d2s|"):
            for key in ((q["size"], q["noise"]), (q["size"], "*"), ("*", q["noise"]), ("*", "*")):
                c = syn[key]
                c["lines"] += 1
                c[o] += 1
                if o == "wrong":
                    c["superset" if g > E else "subset" if g and g < E else "other"] += 1
                    c["drawn_in_answer"] += str(q["drawn"]) in g
        else:
            mats.append((q["dataset"], q["source"], q["features"], q.get("stated"), q.get("rosetta_source_release"), rng(E), rng(g), o,
                         q["n_items"], q.get("best_share")))
    with open(os.path.join(out, "synthetic.tsv"), "w") as fo:
        keys = ["lines", "right", "wrong", "undetermined", "error", "timeout", "superset", "subset", "other", "drawn_in_answer"]
        fo.write("\t".join(["size", "noise"] + keys + ["rate"]) + "\n")
        for (s, n), c in sorted(syn.items(), key=lambda kv: (str(kv[0][0]).zfill(5), kv[0][1])):
            den = c["right"] + c["wrong"] + c["error"]
            fo.write("\t".join(map(str, [s, n] + [c[k] for k in keys] + [f"{c['right'] / den:.4f}" if den else ""])) + "\n")
    with open(os.path.join(out, "matrices.tsv"), "w") as fo:
        fo.write("dataset\tsource\tfeatures\tstated\trosetta_release\texpected\tidtrack\toutcome\titems\tbest_share\n")
        for r in sorted(mats, key=lambda r: (r[0], r[2])):
            fo.write("\t".join(map(str, r)) + "\n")
    print(open(os.path.join(out, "synthetic.tsv")).read())
    print(open(os.path.join(out, "matrices.tsv")).read())


if __name__ == "__main__":
    main()
