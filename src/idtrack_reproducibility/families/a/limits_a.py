"""Family A's known limitations (docs/methods.md, "Known limitations"): the common rows for the wrong answers (limitations.py, from
the causes) plus rows for the answers the truth cannot decide, grouped by the judge's reason (undetermined.tsv.gz of
each scored run). Every number comes from the runs; nothing typed by hand.

Usage: limits_a.py --out <tsv> <scored run folder> ...   (writes F0's columns; the undetermined rows appended)
"""
import argparse
import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import limitations as LIM

# (limitation id, group, title, what, user action, fixed by), matched on the judge's reason in this order
UND = [
    (r"no assembly chain from", ("L-no-chain", "history", "no successor and no assembly chain: the answer cannot be judged",
     "an id Ensembl's history ends, asked across an assembly change for which Ensembl publishes no chain file (pig "
     "Sscrofa10.2 -> 11.1, rhesus, cynomolgus): idtrack's answer (a bridge, or nothing) cannot be checked by place",
     "check the answer's place in the explained path; prefer `allow_bridges=False` when only Ensembl's history counts",
     "")),
    (r"places not comparable", ("L-chain-gap", "history", "no successor and the place does not map: the answer cannot be judged",
     "the old place does not align through the chain file (no aligned block, sequence not in the chain, too little "
     "aligned)", "as L-no-chain", "")),
    (r"every returned gene overlaps the query's place", ("L-bridge-at-place", "history",
     "no successor in Ensembl's history; idtrack answers a gene at the same place",
     "Ensembl's history ends the id; idtrack bridges (rule `bridge_same_place`) to a gene that overlaps the old place: plausible, "
     "not Ensembl's statement", "read the answer as a candidate at the same locus, not as Ensembl's successor",
     "decided: rule `bridge_same_place`")),
    (r"wider set", ("L-branch", "history", "ids that branched off while the old id lived on",
     "the answer lies within the broad walk (branches Ensembl recorded while the id lived on) but outside the strict one",
     "use `all` to see every branch", "decided: rule `main_line_first`")),
    (r"^open: ", ("L-open-truth", "history", "the truth leaves the line open",
     "the gene is absent from the GTF of the release (core-only genes) or its place cannot be projected", "", "")),
]


def classify(reason):
    for pat, row in UND:
        if re.search(pat, reason or ""):
            return row
    return ("L-undetermined-other", "history", "other undetermined answers", "", "", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    rows = [{k: str(v) for k, v in r.items()} for r in
            LIM.build("a", a.runs, LIM.load_map(os.path.join(HERE, "limitations_map_a.tsv")))]
    cols = list(LIM.HEADER)
    for run in a.runs:
        info = json.load(open(os.path.join(run, "run.json")))
        org = info.get("organism")
        totals = collections.Counter()
        for r in csv.DictReader(open(os.path.join(run, "report.tsv")), delimiter="\t"):
            if r["view"] == "*":
                totals[r["stratum"]] = int(r["n"])
        agg = collections.defaultdict(lambda: [0, set(), set()])
        with gzip.open(os.path.join(run, "undetermined.tsv.gz"), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                lim = classify(r["reason"])
                x = agg[lim]
                x[0] += 1
                x[1].add(r["unit"])
                x[2].add(r["stratum"])
        for lim, (n, units, strata) in sorted(agg.items(), key=lambda kv: -kv[1][0]):
            of = sum(totals[s] for s in strata)
            row = dict(limitation=lim[0], group=lim[1], title=lim[2], what=lim[3], organism=org, answers=n,
                       units=len(units), of_answers=of, share=f"{n / of:.5f}" if of else "", family="a",
                       strata=",".join(sorted(strata)), run=LIM.where(run), commit=(info.get("commit") or "")[:12],
                       graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                       **{"class": "e", "cause": "undetermined", "rule": "", "user_action": lim[4],
                          "fixed_by": lim[5], "measured_by": f"a: undetermined.tsv.gz of {LIM.where(run)}",
                          "since": (info.get("commit") or "")[:12]})
            if cols is None:
                cols = list(row)
            rows.append({c: row.get(c, "") for c in cols})
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(len(rows), "rows ->", a.out)


if __name__ == "__main__":
    main()
