"""D2's own report (beside F0's report.tsv): source inference by list size and noise, how the answers miss, and the
matrices line by line.

Usage: report_d2.py <scored folder> --organism <o> [--sets d2_v1]

Usage (the band's price): report_d2.py <scored folder> --organism <o> --price [--sets <sets>]

Writes `<scored>/d2/`:
  synthetic.tsv  size, noise: lines, right, wrong, rate; of the wrong: the answer holds the expected releases and more
                 (`superset`: idtrack ties releases the GTFs tell apart), a part of them (`subset`), other releases
                 (`other`), with the drawn release among the answer (`drawn_in_answer`)
  matrices.tsv   every matrix line: dataset, source, features, stated release, rosetta's release, the expected releases,
                 idtrack's answer, outcome
  price.tsv      (`--price`; the judge accepts any answer between the exact tie and the band, so what the
                 band costs and gains is measured beside the rate) per stratum and per group (`d2_synth_ids`, every
                 synthetic id list; `*`): lines judged, right, wrong, rate; the number of releases per answer (mean,
                 median), of the tie and of the band (mean); the share of answers that hold the drawn release, whose first
                 (best) origin is the drawn release, and that are the exact tie
  lines.tsv.gz   (`--price`) every D2 line: qid, stratum, drawn, tie, band, answer, first origin, outcome
The answer of a right line is its expected set (F0 stores the answer of non-right rows only); `--price` reads the
answers themselves. The first origin is the first release the call returned (the answer's `flags.order`, recorded since
2026-11-20); for an answer stored without it, the documented order of `d4c2dda` and `71690ab` (the tied releases,
ascending: the first is the lowest), else unknown (empty).
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


ASCENDING = ("d4c2dda", "71690ab")  # trees whose infer_source returns the tied releases in ascending order


def price(a):
    """`--price` (module docstring)."""
    import gzip
    import json
    import statistics
    import pandas as pd
    from idtrack_reproducibility.families import common as C
    sets = a.sets.split(",") if a.sets else [x for x in L.set_names("d", a.organism) if x.startswith("d2")]
    lines = {q["qid"]: q for q in L.load("d", a.organism, sets)[0] if q["qid"].startswith("d2")}
    info = json.load(open(os.path.join(a.scored, "run.json")))
    ascending = any(str(info.get("commit") or "").startswith(t) for t in ASCENDING)
    ans = pd.concat([pd.read_parquet(f, columns=["qid", "view", "outcome", "got", "flags"]) for f in C.answer_files(a.scored)])
    ans = ans[ans.qid.isin(lines.keys())].drop_duplicates(["qid", "view"], keep="last")
    outc = RP.read_outcomes(a.scored, columns=["qid", "outcome"])
    outc = dict(zip(outc["qid"], outc["outcome"]))
    out = os.path.join(a.scored, "d2")
    os.makedirs(out, exist_ok=True)
    groups = collections.defaultdict(list)
    with gzip.open(os.path.join(out, "lines.tsv.gz"), "wt") as fo:
        fo.write("qid\tstratum\tdrawn\ttie\tband\tanswer\tfirst\toutcome\n")
        for qid, o, got, flags in zip(ans["qid"], ans["outcome"], ans["got"], ans["flags"]):
            q = lines[qid]
            g = sorted({str(x) for x in (got if got is not None else [])}, key=int) if o == "answer" else []
            order = (json.loads(flags) if isinstance(flags, str) and flags else {}).get("order")
            first = str(order[0]) if order else (g[0] if g and ascending else "")
            oc = outc.get(qid, "")
            drawn = "" if q.get("drawn") is None else str(q["drawn"])
            fo.write("\t".join([qid, q["stratum"], drawn, ",".join(q["expect"]), ",".join(q.get("band") or []), ",".join(g),
                                 first, oc]) + "\n")
            if q.get("open") or oc not in ("right", "wrong", "error"):
                continue
            rec = dict(outcome=oc, width=len(g), tie=len(q["expect"]), band=len(q.get("band") or []),
                       drawn=(drawn in g) if drawn else None, first=(first == drawn) if drawn and first else None,
                       is_tie=set(g) == set(q["expect"]))
            st = q["stratum"]
            for k in (st, "d2_synth_ids" if st.startswith("d2_synth_") and "names" not in st and "versioned" not in st else None, "*"):
                if k:
                    groups[k].append(rec)
    cols = ["stratum", "lines", "right", "wrong", "rate", "width_mean", "width_median", "tie_mean", "band_mean",
            "holds_drawn", "first_is_drawn", "is_tie"]
    with open(os.path.join(out, "price.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for k in sorted(groups):
            rs = groups[k]
            n = len(rs)
            right = sum(r["outcome"] == "right" for r in rs)
            dr = [r["drawn"] for r in rs if r["drawn"] is not None]
            fi = [r["first"] for r in rs if r["first"] is not None]
            row = [k, n, right, n - right, f"{right / n:.4f}", f"{sum(r['width'] for r in rs) / n:.2f}",
                   f"{statistics.median(r['width'] for r in rs):g}", f"{sum(r['tie'] for r in rs) / n:.2f}",
                   f"{sum(r['band'] for r in rs) / n:.2f}", f"{sum(dr) / len(dr):.4f}" if dr else "",
                   f"{sum(fi) / len(fi):.4f}" if fi else "", f"{sum(r['is_tie'] for r in rs) / n:.4f}"]
            fo.write("\t".join(map(str, row)) + "\n")
    print(open(os.path.join(out, "price.tsv")).read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="d2_v1")
    ap.add_argument("--price", action="store_true", help="the band's price: price.tsv and lines.tsv.gz")
    a = ap.parse_args()
    if a.price:
        if a.sets == "d2_v1":
            a.sets = None  # every published D2 set
        return price(a)
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
        elif qid.startswith("d2m|"):
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
