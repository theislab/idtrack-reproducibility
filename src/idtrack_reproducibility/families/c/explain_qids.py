"""The wrong answers of a family C run to explain (qid<TAB>view), for explain.py --qids: every stratum but c6 and
the wrong-version strata (their causes need no path: the database, or the version rule), at most N per stratum and
unit (default 5).

Usage: explain_qids.py <scored folder> [N]   writes <scored>/explain_qids.tsv
"""
import collections
import csv
import gzip
import os
import sys

csv.field_size_limit(1 << 30)
run = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 5
seen = collections.Counter()
out = 0
with gzip.open(os.path.join(run, "wrong.tsv.gz"), "rt") as fh, open(os.path.join(run, "explain_qids.tsv"), "w") as fo:
    for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
        st = r["stratum"]
        if st.startswith("c6") or "_ver_" in st:
            continue
        k = (st, r["unit"])
        if seen[k] >= n:
            continue
        seen[k] += 1
        fo.write(f"{r['qid']}\t{r['variant'].split(',')[0]}\n")
        out += 1
print(run, out, "to explain")
