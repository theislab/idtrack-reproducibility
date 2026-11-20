"""Family A's known limitations (docs/methods.md, "Known limitations"): the common rows for the wrong answers (limitations.py, from
the causes) plus rows for the answers the truth cannot decide, grouped by the judge's reason (undetermined.tsv.gz of
each scored run). Every number comes from the runs; nothing typed by hand.

Usage: limits_a.py --out <tsv> <scored run folder> ...   (writes F0's columns; the undetermined rows appended)
       limits_a.py --a0 --out <tsv> <A0 folder> ...        (A0, graph fidelity: the classes of its summary.tsv that are
                                                             neither right nor by design, mapped by `A0_MAP`)
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


# A0's classes that are neither right (`ok`) nor by design (a0_fidelity.BY_DESIGN): (limitation, group, title, what, user
# action, fixed by); a class without a row is `L-a0-unexplained` (read one by one in the summary's examples).
A0_MAP = {
    "history:between_ids:no_edge": ("L-event-versions", "history", "events whose version contradicts the tables",
                                    "Ensembl's event between two ids names an old version the id did not have at that "
                                    "release; the builder leaves the event out (A0: no edge)", "none",
                                    "decided: rule `event_version`"),
    "transcripts:tl:graph_extra": ("L-id-shared-across-types", "history",
                                   "a translation id equal to a transcript id: one node holds both",
                                   "Ensembl gave a transcript's protein the stable id of another gene's transcript (human 107: "
                                   "ENST00000515292); the graph has one node per name, typed as a translation, so that "
                                   "transcript's links of its own releases are held by a translation node",
                                   "ask such an id with its own release and check the answer's match", "none decided (an open question)"),
    "transcripts:tx:not_linked_at_release": ("L-id-shared-across-types", "history",
                                             "a translation id equal to a transcript id: one node holds both",
                                             "the transcript whose id Ensembl later gave to another gene's protein is not linked "
                                             "to its gene at its own releases (the node is the translation's)",
                                             "ask such an id with its own release and check the answer's match",
                                             "none decided (an open question)"),
}
A0_OTHER = ("L-a0-unexplained", "graph", "graph mismatches without a cause, read one by one",
            "a class of A0 that is neither right nor by design (its examples: a0_summary.tsv)", "-", "-")


def a0_rows(folders):
    """The known-limitation rows of A0 runs: per organism and class that is neither `ok` nor by design."""
    from idtrack_reproducibility.families.a import a0_fidelity as A0
    rows = []
    for d in folders:
        info = json.load(open(os.path.join(d, "summary.json")))
        org, commit = info["organism"], (info.get("commit") or "")[:12]
        digest = (info.get("load") or {}).get("content_digest") or ""
        summ = list(csv.DictReader(open(os.path.join(d, "summary.tsv")), delimiter="\t", quoting=csv.QUOTE_NONE))
        totals = collections.Counter()
        for r in summ:
            totals[r["check"]] += int(r["count"])
        for r in summ:
            cls = r["class"]
            if any(p.startswith("ok") for p in cls.replace("|", ":").split(":")) or cls.startswith("edge_ok") or \
                    A0.design_of(A0.BY_DESIGN.get(r["check"], {}), cls):
                continue
            lim = A0_MAP.get(f"{r['check']}:{cls}", A0_OTHER)
            n, of = int(r["count"]), totals[r["check"]]
            rows.append(dict(limitation=lim[0], group=lim[1], title=lim[2], what=lim[3] if lim is not A0_OTHER else
                             f"{lim[3]}: {r['check']} {cls}", organism=org, answers=n, units=r["ids"], of_answers=of,
                             share=f"{n / of:.5f}" if of else "", family="a", strata=f"a0_{r['check']}", run=LIM.where(d),
                             commit=commit, graph_digest=digest, **{"class": "?" if lim is A0_OTHER else ("c" if lim[0] == "L-event-versions" else "a")},
                             cause=f"a0:{r['check']}:{cls}", rule="", user_action=lim[4], fixed_by=lim[5],
                             measured_by=f"a: A0 summary.tsv of {LIM.where(d)}", since=commit))
    return rows


def classify(reason):
    for pat, row in UND:
        if re.search(pat, reason or ""):
            return row
    return ("L-undetermined-other", "history", "other undetermined answers", "", "", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--a0", action="store_true", help="the folders are A0 runs (summary.tsv, summary.json)")
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    if a.a0:
        LIM.write([{k: str(v) for k, v in r.items()} for r in a0_rows(a.runs)], a.out)
        print("A0 rows ->", a.out)
        return
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
