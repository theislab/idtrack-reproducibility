"""D3's own report (beside F0's report.tsv): idtrack's answers to the marker lists' symbols read against the authority
(HGNC, MGI: principle 2's second truth), and how often the authority and Ensembl's names differ.

Usage: report_d3.py <scored folder> --organism <o> [--sets d3_v1]

Writes `<scored>/d3/authority.tsv`: per stratum and view, lines; with an authority reading (`authority` non-empty);
`ensembl_vs_authority_differ` (lines with an Ensembl truth whose expected genes differ from the authority's);
idtrack's answer equal to the authority's gene(s) (`agree`), another gene (`other`), nothing (`nothing`), on all lines
and apart on the lines Ensembl's names leave open (`open_*`: the string is no Ensembl name at the release read: aliases,
old symbols, papers' spellings), by the authority's reading (`approved`, `previous`, `alias`, `withdrawn`, `synonym`).
A right answer's genes are the expected ones (F0 stores the answer of non-right rows only; a `best` answer that is a
part of several expected genes is counted as agreeing where the authority's gene is among them).
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
    ap.add_argument("--sets", default="d3_v1")
    a = ap.parse_args()
    lines, _, _ = L.load("d", a.organism, a.sets.split(","))
    by = {q["qid"]: q for q in lines if q["qid"].startswith("d3|")}
    df = RP.read_outcomes(a.scored)
    df = df[df["qid"].isin(by.keys())]
    c = collections.defaultdict(collections.Counter)
    for qid, view, o, got, n_got in zip(df["qid"], df["view"], df["outcome"], df["got"], df["n_got"]):
        q = by[qid]
        E, A = set(q.get("expect") or []), set(q.get("authority") or [])
        if o == "right":
            g = E if n_got == len(E) else None
        else:
            g = {x.split(".")[0] for x in (got.split(",") if isinstance(got, str) and got else [])}
        k = (q["stratum"], view)
        x = c[k]
        x["lines"] += 1
        if not A:
            continue
        x["with_authority"] += 1
        if E and E != A:
            x["ensembl_vs_authority_differ"] += 1
        if o in ("error", "timeout"):
            v = o
        elif g is None:
            v = "agree" if A <= E else "other"
        elif not g:
            v = "nothing"
        else:
            v = "agree" if g == A else "other"
        x[v] += 1
        if q.get("open"):
            x["open_" + v] += 1
            x[f"open_{v}_{q.get('authority_how')}"] += 1
    keys = ["lines", "with_authority", "ensembl_vs_authority_differ", "agree", "other", "nothing", "open_agree", "open_other", "open_nothing"]
    hows = sorted({k for x in c.values() for k in x if k.count("_") >= 2 and k.startswith("open_") and k.split("_")[1] in ("agree", "other", "nothing")
                   and k not in keys})
    os.makedirs(os.path.join(a.scored, "d3"), exist_ok=True)
    out = os.path.join(a.scored, "d3", "authority.tsv")
    with open(out, "w") as fo:
        fo.write("\t".join(["stratum", "view"] + keys + hows) + "\n")
        for (st, v), x in sorted(c.items()):
            fo.write("\t".join(map(str, [st, v] + [x[k] for k in keys] + [x[k] for k in hows])) + "\n")
    print(open(out).read())


if __name__ == "__main__":
    main()
