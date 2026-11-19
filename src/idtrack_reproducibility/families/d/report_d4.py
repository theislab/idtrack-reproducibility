"""D4's own report (beside F0's report.tsv): every item of a mixed list against the same line asked alone.

Usage: report_d4.py <scored folder> --organism <o> --single <scored folder of the harness's sets on the same tree>
       [--sets d4_v1]

The single calls are F0's check run of the harness's sets on the same tree and graph (`<families>/f0/runs/check_<tree>/<o>`).
Writes `<scored>/d4/batch_vs_single.tsv`: per stratum and view, items, their outcome in the list and alone: both right,
right in the list only, right alone only (an answer the list changes: a finding), both wrong, and the rest; and
`<scored>/d4/list_only.tsv` (every item whose outcome in the list differs from its single call, and whether the two
answers are the same: an outcome can differ on one answer where the single call's run was scored under an earlier
judge, `outcome_differs_same_answer`).
"""
import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--single", required=True)
    ap.add_argument("--sets", default="d4_v1")
    a = ap.parse_args()
    lines, _, _ = L.load("d", a.organism, a.sets.split(","))
    by = {q["qid"]: q for q in lines if q["qid"].startswith("d4|")}
    df = RP.read_outcomes(a.scored)
    df = df[df["qid"].isin(by.keys())]
    want = {q.get("src_qid") for q in by.values() if q.get("src_qid")}
    single = RP.read_outcomes(a.single, columns=["qid", "view", "outcome", "got"])
    single = single[single["qid"].isin(want)]
    alone = {(q, v): (o, g) for q, v, o, g in zip(single["qid"], single["view"], single["outcome"], single["got"])}
    c = collections.defaultdict(collections.Counter)
    diffs = []
    for qid, view, o, got in zip(df["qid"], df["view"], df["outcome"], df["got"]):
        q = by[qid]
        s = alone.get((q.get("src_qid"), view)) if q.get("src_qid") else None
        so = s[0] if s else "no_single"
        if o in ("right", "expected_error") and so == "right":
            k = "both_right"
        elif o in ("right", "expected_error") and so == "no_single":
            k = "right_no_single"
        elif o in ("right", "expected_error"):
            k = "right_in_list_only"
        elif so == "right":
            k = "right_alone_only"
        elif o == so:
            k = f"both_{o}"
        else:
            k = f"list_{o}_alone_{so}"
        for key in ((q["stratum"], view), ("*", view)):
            c[key][k] += 1
            c[key]["items"] += 1
        if k in ("right_in_list_only", "right_alone_only") or k.startswith("list_"):
            ga, gb = got or "", (s[1] if s else "") or ""
            same = "same_answer" if (ga and gb and ga == gb) else ("different_answer" if (ga and gb) else "unknown")
            if s and so != "no_single":
                c[(q["stratum"], view)][f"outcome_differs_{same}"] += 1
                c[("*", view)][f"outcome_differs_{same}"] += 1
            diffs.append((q["stratum"], view, qid, q.get("src_qid"), q["id"], o, ga, so, gb, same, q.get("list_size"), q.get("composition")))
    out = os.path.join(a.scored, "d4")
    os.makedirs(out, exist_ok=True)
    keys = sorted({k for x in c.values() for k in x if k != "items"})
    with open(os.path.join(out, "batch_vs_single.tsv"), "w") as fo:
        fo.write("\t".join(["stratum", "view", "items"] + keys) + "\n")
        for (st, v), x in sorted(c.items()):
            fo.write("\t".join(map(str, [st, v, x["items"]] + [x[k] for k in keys])) + "\n")
    with open(os.path.join(out, "list_only.tsv"), "w") as fo:
        fo.write("stratum\tview\tqid\tsrc_qid\tid\tin_list\tgot_in_list\talone\tgot_alone\tanswers\tlist_size\tcomposition\n")
        for r in diffs:
            fo.write("\t".join(map(str, r)) + "\n")
    print(open(os.path.join(out, "batch_vs_single.tsv")).read())
    print(len(diffs), "items whose outcome in the list differs from the single call")


if __name__ == "__main__":
    main()
