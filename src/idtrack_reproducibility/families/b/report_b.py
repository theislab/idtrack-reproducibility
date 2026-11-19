"""Family B's report of one run (tree, organism): both judges side by side per stratum, the authority measure, the
name sources the truth read, and the causes.

Usage: report_b.py --organism <o> --run <run folder>/<organism>

Writes into the run folder:
  report_b.tsv      per judge and stratum (and the groups B1..B8): n, right, wrong, undetermined, timeout, error, rate
                    with Wilson and clustered intervals (F0's report.tsv of each judge, view `*`), units behind the
                    wrong answers
  authority_b.tsv   per stratum of the lines that carry an authority's gene (B4, B5's aliases): lines, the authority
                    names a gene at 115, Ensembl's truth equals the authority's gene (agree), differs; idtrack's `best`
                    answer equals the authority's gene, holds it among others, misses it, answers nothing
  name_src_b.tsv    per stratum and the tables the truth read at the line's release (`name_src`): rows and outcomes
"""
import argparse
import collections
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP

GROUP = {"b1": "B1 every name at every release", "b2": "B2 moved names", "b3": "B3 names of several genes",
         "b4": "B4 previous symbols", "b5": "B5 synonyms and aliases", "b6": "B6 forms", "b7": "B7 spreadsheet and typography",
         "b8": "B8 names to databases"}


def group_of(stratum):
    if stratum.startswith("b1_") and stratum.endswith("_shared"):
        return "b3"
    return stratum[:2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    a = ap.parse_args()
    run, org = os.path.abspath(a.run), a.organism
    rows = []
    for judge in ("now", "decided"):
        p = os.path.join(run, f"scores_{judge}", "report.tsv")
        if not os.path.exists(p):
            continue
        rep = [r for r in csv.DictReader(open(p), delimiter="\t") if r["view"] == "*"]
        agg = collections.defaultdict(lambda: collections.Counter())
        for r in rep:
            rows.append([judge, r["stratum"], group_of(r["stratum"]), r["n"], r["right"], r["expected_error"], r["wrong"],
                         r["undetermined"], r["timeout"], r["error"], r["rate"], r["wilson_low"], r["wilson_high"], r["boot_low"],
                         r["boot_high"], r["units"], r["units_wrong"]])
            for k in ("n", "right", "expected_error", "wrong", "undetermined", "timeout", "error"):
                agg[group_of(r["stratum"])][k] += int(r[k])
        for g, c in sorted(agg.items()):
            ok = c["right"] + c["expected_error"]
            den = ok + c["wrong"] + c["error"]
            lo, hi = RP.wilson(ok, den)
            rows.append([judge, f"group {g}: {GROUP.get(g, g)}", g, c["n"], c["right"], c["expected_error"], c["wrong"],
                         c["undetermined"], c["timeout"], c["error"], f"{ok / den:.5f}" if den else "nan", f"{lo:.5f}", f"{hi:.5f}",
                         "", "", "", ""])
    with open(os.path.join(run, "report_b.tsv"), "w") as fo:
        fo.write("judge\tstratum\tgroup\tn\tright\texpected_error\twrong\tundetermined\ttimeout\terror\trate\twilson_low\t"
                 "wilson_high\tboot_low\tboot_high\tunits\tunits_wrong\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    # the authority measure (lines with an authority's gene) and the name sources read, from the lines and answers
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    info = json.load(open(os.path.join(run, "scores_now", "run.json")))
    sets = list(info["sets"])
    auth = {}
    src = {}
    d = L.query_dir("b", org)
    for s in sets:
        p = L._set_file(d, s)
        if p.endswith(".parquet"):
            t = pq.read_table(p, columns=["qid", "stratum", "expect", "extra"])
            for q, st, e, x in zip(*(t.column(c).to_pylist() for c in ("qid", "stratum", "expect", "extra"))):
                xx = json.loads(x) if x else {}
                src[q] = (st, xx.get("name_src", ""))
                if "authority" in xx:
                    auth[q] = (st, set(e or []), set(xx["authority"]))
        else:
            for ln in L.iter_set(p):
                src[ln["qid"]] = (ln.get("stratum"), ln.get("name_src", ""))
                if "authority" in ln:
                    auth[ln["qid"]] = (ln.get("stratum"), set(ln.get("expect") or []), set(ln["authority"]))
    ans = {}
    if auth:
        want = pa.array(sorted(auth))
        for p in glob.glob(os.path.join(run, "answers", "*.parquet")):
            t = pq.read_table(p, columns=["qid", "view", "outcome", "got"])
            t = t.filter(pc.and_(pc.is_in(t["qid"], value_set=want), pc.equal(t["view"], "best")))
            for q, o, g in zip(*(t.column(c).to_pylist() for c in ("qid", "outcome", "got"))):
                ans[q] = (o, {C.base_id(x) for x in (g or [])})
    acc = collections.defaultdict(collections.Counter)
    for q, (st, E, A) in auth.items():
        c = acc[st]
        c["lines"] += 1
        Ab, Eb = {C.base_id(x) for x in A}, {C.base_id(x) for x in E}
        if not Ab:
            c["authority_no_gene_at_115"] += 1
            continue
        c["authority_gene"] += 1
        c["ensembl_agrees" if Eb == Ab else ("ensembl_overlaps" if Eb & Ab else "ensembl_differs")] += 1
        o, g = ans.get(q, (None, set()))
        if o != "answer":
            c["idtrack_no_answer_" + str(o)] += 1
        elif not g:
            c["idtrack_nothing"] += 1
        elif g == Ab:
            c["idtrack_equals_authority"] += 1
        elif Ab & g:
            c["idtrack_holds_authority"] += 1
        else:
            c["idtrack_misses_authority"] += 1
    keys = ["lines", "authority_no_gene_at_115", "authority_gene", "ensembl_agrees", "ensembl_overlaps", "ensembl_differs",
            "idtrack_equals_authority", "idtrack_holds_authority", "idtrack_misses_authority", "idtrack_nothing"]
    with open(os.path.join(run, "authority_b.tsv"), "w") as fo:
        fo.write("stratum\t" + "\t".join(keys) + "\n")
        for st, c in sorted(acc.items()):
            fo.write(st + "\t" + "\t".join(str(c[k]) for k in keys) + "\n")
    # outcomes by the tables the truth read at the line's release (judge now)
    out = collections.defaultdict(collections.Counter)
    for p in RP.outcome_parts(os.path.join(run, "scores_now")):
        t = pq.read_table(p, columns=["qid", "outcome"])
        for q, o in zip(t["qid"].to_pylist(), t["outcome"].to_pylist()):
            st, ns = src.get(q, ("?", ""))
            out[(st, ns)][o] += 1
    with open(os.path.join(run, "name_src_b.tsv"), "w") as fo:
        fo.write("stratum\tname_src\trows\tright\twrong\tundetermined\terror\ttimeout\n")
        for (st, ns), c in sorted(out.items()):
            fo.write(f"{st}\t{ns}\t{sum(c.values())}\t{c['right']}\t{c['wrong']}\t{c['undetermined']}\t{c['error']}\t{c['timeout']}\n")
    print(open(os.path.join(run, "report_b.tsv")).read()[:3000])


if __name__ == "__main__":
    main()
