"""D4's wrong answers classified by their single call (no graph needed): an item answered in a list as alone (D4's
`report_d4.py` found no answer that differs) has the cause F0's classifier gave the same harness line on the same tree
(`<families>/f0/runs/check_<tree>/<o>/causes.tsv.gz`, same view). Rows without a classified single call (the `all` view of
a harness line asked only with `best`) are left to `causes_job.sh` (explain + classify).

Usage: causes_d4_from_single.py <scored D4 folder> --organism <o> --single <check run folder> [--sets d4_v1,d4_r1,d4_r2]
Writes `<scored>/d4/causes_from_single.tsv` (qid, view, src_qid, cause, class, detail) and prints the summary.
"""
import argparse
import collections
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L

csv.field_size_limit(1 << 30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--single", required=True)
    ap.add_argument("--sets", default="d4_v1,d4_r1,d4_r2")
    a = ap.parse_args()
    by = {q["qid"]: q for q in L.load("d", a.organism, a.sets.split(","))[0] if q["qid"].startswith("d4|")}
    single = {}
    with gzip.open(os.path.join(a.single, "causes.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            single[(r["qid"], r["view"])] = (r["cause"], r["class"], r["detail"])
    out, summ = [], collections.Counter()
    with gzip.open(os.path.join(a.scored, "wrong.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            q = by.get(r["qid"])
            view = r["variant"].split(",")[0]
            src = (q or {}).get("src_qid")
            c = single.get((src, view)) if src else None
            if c is None:
                c = ("not_classified_by_single", "?", "no classified single call in this view")
            out.append((r["qid"], view, src, *c))
            summ[(r["stratum"], c[0], c[1])] += 1
    os.makedirs(os.path.join(a.scored, "d4"), exist_ok=True)
    with open(os.path.join(a.scored, "d4", "causes_from_single.tsv"), "w") as fo:
        fo.write("qid\tview\tsrc_qid\tcause\tclass\tdetail\n")
        for row in out:
            fo.write("\t".join(str(x).replace("\t", " ") for x in row) + "\n")
    for k, v in sorted(summ.items()):
        print("\t".join(map(str, k)), v)


if __name__ == "__main__":
    main()
