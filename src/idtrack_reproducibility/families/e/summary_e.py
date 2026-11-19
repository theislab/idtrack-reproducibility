"""Family E's table over trees, judges and organisms, from every run's `e_summary.tsv`, E7's, E6's, E8's and E9's reports.

Usage: summary_e.py [--trees d4c2dda,71690ab] --out <families>/e/e_rates.tsv

Rows: tree, judge, organism, group (a stratum without the kind of identifier its forms were made from: `e2_zw_name`
and `e2_zw_tx` are `e2_zw`), n, right, expected_error, wrong, undetermined, timeout, error, rate, wilson_low,
wilson_high, answered, undetermined_answered. Rate = (right + expected_error) / (right + expected_error + wrong + error).
"""
import argparse
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import report as RP

KIND = re.compile(r"_(name|syn|gv|gu|tx|tl|ncbi|refseq|uniprot)(_unversioned|_versioned)?$")


def group(st):
    if st.startswith("e2_") and not st.startswith("e2_case"):
        return KIND.sub(lambda m: m.group(2) or "", st)
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", default="d4c2dda,71690ab")
    ap.add_argument("--out", default=os.path.join(C.BENCH, "e", "e_rates.tsv"))
    a = ap.parse_args()
    cols = ["n", "right", "expected_error", "wrong", "undetermined", "timeout", "error", "answered", "undetermined_answered"]
    agg = collections.defaultdict(lambda: collections.Counter())
    for t in a.trees.split(","):
        for o in C.ORGANISMS:
            p = os.path.join(C.BENCH, "e", "runs", t, o, "e_summary.tsv")
            if not os.path.exists(p):
                continue
            for r in csv.DictReader(open(p), delimiter="\t"):
                for k in (o, "all"):
                    c = agg[(t, r["judge"], k, group(r["stratum"]))]
                    for x in cols:
                        c[x] += int(r[x])
    with open(a.out, "w") as fo:
        fo.write("\t".join(["tree", "judge", "organism", "group"] + cols[:7] + ["rate", "wilson_low", "wilson_high"] + cols[7:]) + "\n")
        for (t, j, o, g), c in sorted(agg.items()):
            ok, d = c["right"] + c["expected_error"], c["right"] + c["expected_error"] + c["wrong"] + c["error"]
            lo, hi = RP.wilson(ok, d) if d else (float("nan"), float("nan"))
            fo.write("\t".join([t, j, o, g] + [str(c[x]) for x in cols[:7]] +
                               ([f"{ok / d:.5f}", f"{lo:.5f}", f"{hi:.5f}"] if d else ["", "", ""]) + [str(c[x]) for x in cols[7:]]) + "\n")
    print(a.out)


if __name__ == "__main__":
    main()
