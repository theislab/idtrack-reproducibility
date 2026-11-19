"""Compare two scored runs of a family (F0): the harness's compare_scores.py over Parquet outcomes.

Usage: compare.py <scored a> <scored b> [--out changes.tsv] [--moved moved.tsv] [--views]

Each argument is a scored folder (outcomes/, report.tsv, run.json; the answers are read from the run folder that
run.json names). Prints per stratum the rates of a and b with their clustered intervals and the counts, then the
moves of every (qid, view) scored in both: regressions (right -> wrong, right -> error), right -> timeout, right ->
undetermined, wrong -> right, other changes (an outcome or an answer that changed), and every moved answer counted
by its move and the judge of its line in b. --moved writes one line per answer whose outcome moved; --out every
changed answer (also those whose outcome stayed and whose ids changed). Exit status 1 on a regression.

Streaming: the outcomes and answers of both runs are read in hash buckets of qids (about 15 million rows each), so
runs of family A's size compare in the memory of one bucket; only the changed answers are kept.
"""
import argparse
import glob
import json
import os
import sys
from collections import Counter

from idtrack_reproducibility.families import report


def bucket(qids, k):
    """A stable bucket (0..k-1) per qid (pandas' fixed-key hash: the same in every process)."""
    import pandas as pd
    return pd.util.hash_array(qids.to_numpy(dtype=object)) % k


def n_rows(folder):
    import pyarrow.parquet as pq
    return sum(pq.ParquetFile(p).metadata.num_rows for p in report.outcome_parts(folder))


def outcomes_bucket(folder, cols, b, k):
    import pandas as pd
    import pyarrow.parquet as pq
    out = []
    for p in report.outcome_parts(folder):
        df = pq.read_table(p, columns=cols).to_pandas()
        out.append(df[bucket(df["qid"], k) == b] if k > 1 else df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=cols)


def answers_bucket(folder, b, k):
    """The answers (ids joined by `,`) of the run folder that run.json names, one bucket; the last answer wins."""
    import pandas as pd
    import pyarrow.parquet as pq
    info = json.load(open(os.path.join(folder, "run.json")))
    out = []
    for p in sorted(glob.glob(os.path.join(info.get("run", folder), "answers", "*.parquet"))):
        df = pq.read_table(p, columns=["qid", "view", "got"]).to_pandas()
        out.append(df[bucket(df["qid"], k) == b] if k > 1 else df)
    if not out:
        return pd.DataFrame(columns=["qid", "view", "got"])
    df = pd.concat(out, ignore_index=True).drop_duplicates(["qid", "view"], keep="last")
    df["got"] = [",".join(g) if g is not None else "" for g in df["got"]]
    return df


def read_report(folder):
    import csv
    p = os.path.join(folder, "report.tsv")
    with open(p) as fh:
        return {(r["stratum"], r["view"]): r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--out", default=None)
    ap.add_argument("--moved", default=None)
    ap.add_argument("--views", action="store_true", help="also every view row of the report")
    a = ap.parse_args()
    ia = json.load(open(os.path.join(a.a, "run.json")))
    ib = json.load(open(os.path.join(a.b, "run.json")))
    for s in sorted(set(ia["sets"]) | set(ib["sets"])):
        if ia["sets"].get(s) != ib["sets"].get(s):
            print(f"NOTE: query set {s} differs or is missing in one run")
    print(f"a: {ia.get('commit', '')[:10]} {ia.get('run')}\nb: {ib.get('commit', '')[:10]} {ib.get('run')}\n")
    ra, rb = read_report(a.a), read_report(a.b)
    print(f"{'stratum':34s} {'view':5s} {'n':>8s} {'rate a':>8s} {'boot a':>15s} {'rate b':>8s} {'boot b':>15s} {'change':>8s}  "
          "right/wrong/undet/timeout/error a -> b")
    for key in sorted(set(ra) | set(rb)):
        if key[1] != "*" and not a.views:
            continue
        x, y = ra.get(key), rb.get(key)

        def f(r):
            return float(r["rate"]) if r and r["rate"] != "nan" else float("nan")

        def iv(r):
            return f"[{r['boot_low']},{r['boot_high']}]" if r else "-"

        def c(r):
            return "/".join(str(int(r["right"]) + int(r["expected_error"])) if k == "right" else r[k]
                            for k in ("right", "wrong", "undetermined", "timeout", "error")) if r else "-"

        print(f"{key[0]:34s} {key[1]:5s} {(y or x)['n']:>8s} {f(x):8.4f} {iv(x):>15s} {f(y):8.4f} {iv(y):>15s} "
              f"{f(y) - f(x):+8.4f}  {c(x)} -> {c(y)}")
    cols = ["qid", "view", "kind", "stratum", "variant", "outcome", "reason", "rule"]
    order = ["right -> wrong (regression)", "right -> error (regression)", "right -> timeout", "right -> undetermined",
             "wrong -> right", "other change"]
    k = max(1, -(-max(n_rows(a.a), n_rows(a.b)) // 15_000_000))  # buckets of at most about 15 million rows
    changed, n_common, only_a, only_b = [], 0, 0, 0
    for b in range(k):
        oa, ob = outcomes_bucket(a.a, cols, b, k), outcomes_bucket(a.b, cols, b, k)
        m = oa.merge(ob, on=["qid", "view"], how="outer", suffixes=("_a", "_b"), indicator=True)
        only_a += int((m["_merge"] == "left_only").sum())
        only_b += int((m["_merge"] == "right_only").sum())
        m = m[m["_merge"] == "both"].drop(columns="_merge")
        n_common += len(m)
        m = m.merge(answers_bucket(a.a, b, k).rename(columns={"got": "got_a"}), on=["qid", "view"], how="left")
        m = m.merge(answers_bucket(a.b, b, k).rename(columns={"got": "got_b"}), on=["qid", "view"], how="left")
        m[["got_a", "got_b"]] = m[["got_a", "got_b"]].fillna("")
        m = m[(m["outcome_a"] != m["outcome_b"]) | (m["got_a"] != m["got_b"])]
        if len(m):
            changed.append(m)
    import pandas as pd
    ch = pd.concat(changed, ignore_index=True).sort_values(["qid", "view"]) if changed else pd.DataFrame(
        columns=[c + s for c in cols[2:] for s in ("_a", "_b")] + ["qid", "view", "got_a", "got_b"])
    x, y = ch["outcome_a"], ch["outcome_b"]
    group = pd.Series("other change", index=ch.index)
    group[(x == "right") & (y == "wrong")] = order[0]
    group[(x == "right") & (y == "error")] = order[1]
    group[(x == "right") & (y == "timeout")] = order[2]
    group[(x == "right") & (y == "undetermined")] = order[3]
    group[(x == "wrong") & (y == "right")] = order[4]
    ch["group"] = group
    ch["move"] = x + " -> " + y
    moved = ch[x != y]
    print(f"\nanswers compared: {n_common} (only in a: {only_a}, only in b: {only_b})")
    for g in order:
        print(f"{g}: {int((ch['group'] == g).sum())}")
    for g in order[:3]:
        sub = ch[ch["group"] == g]
        for r in sub.head(200).itertuples(index=False):
            print(f"  {g.split(' (')[0]}\t{r.qid}\t{r.view}\tgot a: {r.got_a}\tgot b: {r.got_b}\t{r.reason_b}")
        if len(sub) > 200:
            print(f"  ... {len(sub) - 200} more (see --out)")
    tally = Counter(zip(moved["move"], moved["rule_b"].fillna("published").replace("", "published")))
    print(f"\nanswers whose outcome moved: {len(moved)} (by move and by the judge of the line in b)")
    for (mv, rule), n in sorted(tally.items(), key=lambda kv: (kv[0][0], -kv[1])):
        print(f"  {n:8d}  {mv:28s} {rule}")
    if a.moved:
        with open(a.moved, "w") as fo:
            fo.write("qid\tview\tkind\tstratum\tvariant\tmove\toutcome_a\trule_a\toutcome_b\trule_b\tgot_a\tgot_b\treason_a\treason_b\n")
            for r in moved.itertuples(index=False):
                fo.write("\t".join(str(v) for v in [r.qid, r.view, r.kind_b, r.stratum_b, r.variant_b, r.move, r.outcome_a, r.rule_a,
                                                    r.outcome_b, r.rule_b, r.got_a, r.got_b, r.reason_a, r.reason_b]) + "\n")
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("change\tqid\tview\toutcome_a\toutcome_b\tgot_a\tgot_b\treason_b\n")
            for g in order:
                for r in ch[ch["group"] == g].itertuples(index=False):
                    fo.write("\t".join(str(v) for v in [g, r.qid, r.view, r.outcome_a, r.outcome_b, r.got_a, r.got_b, r.reason_b]) + "\n")
    sys.exit(1 if (ch["group"].isin(order[:2])).any() else 0)


if __name__ == "__main__":
    main()
