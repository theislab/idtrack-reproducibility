"""Compare two scored runs of the accuracy harness.

Usage: compare_scores.py <run a> <run b> [--out <file.tsv>] [--moved <file.tsv>] [--variants]

<run> is an output folder of score.py (scores.tsv, outcomes.tsv.gz, run.json). Prints, per kind (variant `*`), the
rates of a and b and the change, then the regressions: queries right in a and wrong in b. Also listed, apart:
right in a and error or timeout in b, and right in a and undetermined in b. Only queries scored in both runs are
compared; a different query set (sha256) or time limit is reported. Exit status 1 when there is a regression
(right -> wrong or right -> error), 0 otherwise. With --out, every changed query is written to that file.

Every answer that moved between right, wrong and undetermined (and to or from error and timeout) is counted by its
move and by the judge of its line in b (`rule`: `published`, or `revised:<corrections of the corrected judge>` as score.py
writes it; runs scored before the column existed read as `published`); --moved writes one line per moved answer
(qid, kind, variant, move, outcome and rule in a and in b, answers, reasons), so that every changed answer of a
stage can be labelled with the rule that judged it.
"""
import argparse
import csv
import gzip
import json
import os
import sys


def read_scores(d):
    with open(os.path.join(d, "scores.tsv")) as fh:
        return {(r["kind"], r["variant"]): r for r in csv.DictReader(fh, delimiter="\t")}


def read_outcomes(d):
    with gzip.open(os.path.join(d, "outcomes.tsv.gz"), "rt") as fh:
        return {r["qid"]: r for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--out", default=None)
    ap.add_argument("--variants", action="store_true", help="also show every variant row")
    ap.add_argument("--moved", default=None, help="write every answer whose outcome moved, with the judging rule")
    a = ap.parse_args()
    ra, rb = (json.load(open(os.path.join(x, "run.json"))) for x in (a.a, a.b))
    for k in ("organism", "timeout"):
        if ra.get(k) != rb.get(k):
            print(f"NOTE: {k} differs: {ra.get(k)} vs {rb.get(k)}")
    for s in sorted(set(ra["sets"]) | set(rb["sets"])):
        if ra["sets"].get(s) != rb["sets"].get(s):
            print(f"NOTE: query set {s} differs or is missing in one run")
    print(f"a: {ra.get('commit', '')[:10]} {ra['graph']}\nb: {rb.get('commit', '')[:10]} {rb['graph']}\n")
    sa, sb = read_scores(a.a), read_scores(a.b)
    print(f"{'kind':34s} {'variant':14s} {'n':>7s} {'rate a':>8s} {'rate b':>8s} {'change':>8s}  right/wrong/undet/timeout/error a -> b")
    for key in sorted(set(sa) | set(sb)):
        if key[1] != "*" and not a.variants:
            continue
        x, y = sa.get(key), sb.get(key)

        def f(r):
            return float(r["rate"]) if r and r["rate"] != "nan" else float("nan")

        def c(r):
            return "/".join(r[k] for k in ("right", "wrong", "undetermined", "timeout", "error")) if r else "-"

        d = f(y) - f(x)
        print(f"{key[0]:34s} {key[1]:14s} {(y or x)['n']:>7s} {f(x):8.4f} {f(y):8.4f} {d:+8.4f}  {c(x)} -> {c(y)}")
    oa, ob = read_outcomes(a.a), read_outcomes(a.b)
    common = set(oa) & set(ob)
    groups = {"right -> wrong (regression)": [], "right -> error (regression)": [], "right -> timeout": [],
              "right -> undetermined": [], "wrong -> right": [], "other change": []}
    for q in sorted(common):
        x, y = oa[q]["outcome"], ob[q]["outcome"]
        if x == y and oa[q]["got"] == ob[q]["got"]:
            continue
        if x == "right" and y == "wrong":
            groups["right -> wrong (regression)"].append(q)
        elif x == "right" and y == "error":
            groups["right -> error (regression)"].append(q)
        elif x == "right" and y == "timeout":
            groups["right -> timeout"].append(q)
        elif x == "right" and y == "undetermined":
            groups["right -> undetermined"].append(q)
        elif y == "right" and x != "right":
            groups["wrong -> right" if x == "wrong" else "other change"].append(q)
        else:
            groups["other change"].append(q)
    print(f"\nqueries compared: {len(common)} (only in a: {len(set(oa) - common)}, only in b: {len(set(ob) - common)})")
    for g, qs in groups.items():
        print(f"{g}: {len(qs)}")
    for g in ("right -> wrong (regression)", "right -> error (regression)", "right -> timeout"):
        for q in groups[g][:200]:
            print(f"  {g.split(' (')[0]}\t{q}\tgot a: {oa[q]['got']}\tgot b: {ob[q]['got']}\t{ob[q]['reason']}")
        if len(groups[g]) > 200:
            print(f"  ... {len(groups[g]) - 200} more (see --out)")
    moved = []
    for q in sorted(common):
        x, y = oa[q]["outcome"], ob[q]["outcome"]
        if x != y:
            moved.append((f"{x} -> {y}", q))
    tally = {}
    for move, q in moved:
        key = (move, ob[q].get("rule") or "published")
        tally[key] = tally.get(key, 0) + 1
    print(f"\nanswers whose outcome moved: {len(moved)} (by move and by the judge of the line in b)")
    for (move, rule), n in sorted(tally.items(), key=lambda kv: (kv[0][0], -kv[1])):
        print(f"  {n:7d}  {move:28s} {rule}")
    if a.moved:
        with open(a.moved, "w") as fo:
            fo.write("qid\tkind\tvariant\tmove\toutcome_a\trule_a\toutcome_b\trule_b\tgot_a\tgot_b\treason_a\treason_b\n")
            for move, q in moved:
                x, y = oa[q], ob[q]
                fo.write("\t".join([q, y["kind"], y["variant"], move, x["outcome"], x.get("rule") or "published", y["outcome"],
                                    y.get("rule") or "published", x["got"], y["got"], x["reason"], y["reason"]]) + "\n")
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("change\tqid\toutcome_a\toutcome_b\tgot_a\tgot_b\treason_b\n")
            for g, qs in groups.items():
                for q in qs:
                    fo.write("\t".join([g, q, oa[q]["outcome"], ob[q]["outcome"], oa[q]["got"], ob[q]["got"], ob[q]["reason"]]) + "\n")
    sys.exit(1 if groups["right -> wrong (regression)"] or groups["right -> error (regression)"] else 0)


if __name__ == "__main__":
    main()
