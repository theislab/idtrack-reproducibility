"""The check of F0 (docs/methods.md, F0): the harness's sets asked through the adapter and scored by the scorer give
the same outcomes and answers as the stored runs of score.py.

Usage: check_harness.py <scored run folder>/<organism> <stored score.py run folder> [<stored ...>] [--out diff.tsv]

For every qid scored in both: the outcome (right, wrong, undetermined, timeout, error), the reason, and the answer
(the ids; ours from the answer files, score.py's from outcomes.tsv.gz) must be equal. Prints the counts, the
outcome table of each side and every difference (the first 50; all with --out). Exit status 1 on any difference.
"""
import argparse
import csv
import glob
import gzip
import os
import sys

from idtrack_reproducibility.families import report

csv.field_size_limit(1 << 30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ours")
    ap.add_argument("stored", nargs="+")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    import pyarrow as pa
    import pyarrow.parquet as pq
    st = {}
    for d in a.stored:
        with gzip.open(os.path.join(d, "outcomes.tsv.gz"), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                st[r["qid"]] = r
    df = report.read_outcomes(a.ours, columns=["qid", "view", "kind", "outcome", "reason"])
    ans = pa.concat_tables([pq.read_table(p, columns=["qid", "view", "got"]) for p in
                            glob.glob(os.path.join(a.ours, "answers", "*.parquet"))]).to_pandas()
    got = {(q, v): ",".join(g) if g is not None else "" for q, v, g in zip(ans["qid"], ans["view"], ans["got"])}
    diffs, n, only_ours = [], 0, 0
    for q, v, k, o, why in zip(df["qid"], df["view"], df["kind"], df["outcome"], df["reason"]):
        s = st.get(q)
        if s is None:
            only_ours += 1
            continue
        n += 1
        g = got.get((q, v), "")
        if s["outcome"] != o or s["got"] != g or (s["reason"] or "") != (why or ""):
            diffs.append((q, k, s["outcome"], o, s["got"], g, s["reason"], why))
    ours = set(df["qid"])
    only_stored = sum(1 for q in st if q not in ours)
    print(f"compared {n} lines; only in ours {only_ours}; only in the stored runs {only_stored}; differences {len(diffs)}")
    from collections import Counter
    print("stored:", dict(Counter(st[q]["outcome"] for q in ours if q in st)))
    print("ours:  ", dict(Counter(o for q, o in zip(df["qid"], df["outcome"]) if q in st)))
    for d in diffs[:50]:
        print("DIFF", *d, sep="\t")
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("qid\tkind\toutcome_stored\toutcome_ours\tgot_stored\tgot_ours\treason_stored\treason_ours\n")
            for d in diffs:
                fo.write("\t".join(str(x) for x in d) + "\n")
    sys.exit(1 if diffs else 0)


if __name__ == "__main__":
    main()
