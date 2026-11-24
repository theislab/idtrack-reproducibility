"""Family F's known limitations (docs/methods.md, "Known limitations"): F0's rows for the wrong assignments (limitations.py,
from the causes; F's map `limitations_map_f.tsv`), rows for the assignments the truth leaves open, by the judge's
reason (undetermined.tsv.gz of each scoring), and rows for the checks of the matrices, metadata, order and cache that
do not hold (`f_checks.tsv` of the run: each such check a class of outputs the harmonizer does not give as its
documentation says; the `truth_*` checks are left out, as an output that differs from the truth's on decided features
is a wrong assignment, which the causes count), and a row of the cost of reading the result cache (`cache.json`: the
constructor again over its valid cache against the first conversion; idtrack's documented cost). Every number comes from the runs; nothing typed by hand.

Usage: limits_f.py --out <tsv> <run folder or scored folder> ...   (a run folder: its scorings `scores_<config>`)
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import limitations as LIM  # noqa: E402
from idtrack_reproducibility.families.f import causes_f  # noqa: E402,F401  F's causes and their classes

# (pattern on the judge's reason) -> (limitation, group, title, what, user action, fixed by)
UND = [
    (r"^open: under the default final database", ("L-harmonize-tie-final-db", "harmonization",
     "a feature of several genes under the default final database",
     "the harmonizer judges a tie among the genes with an id in the final database (its documentation): which of a "
     "name's genes have one is not read by the truth, so the gene it keeps cannot be judged",
     "pass `final_database=\"ensembl_gene\"` to keep every gene of the tie, or read `idtrack_status`", "none (by design)")),
    (r"wider set", ("L-branch", "history", "ids that branched off while the old id lived on",
     "the answer lies within the broad walk (branches Ensembl recorded while the id lived on) but outside the strict one",
     "use `all` to see every branch", "decided: rule `main_line_first`")),
    (r"no assembly chain from|places not comparable", ("L-no-chain", "history",
     "no successor and no assembly chain: the answer cannot be judged",
     "an id Ensembl's history ends, asked across an assembly change without a chain file: the gene answered cannot be "
     "checked by place", "check the answer's place in the explained path", "")),
    (r"every returned gene overlaps the query's place", ("L-bridge-at-place", "history",
     "no successor in Ensembl's history; idtrack answers a gene at the same place",
     "Ensembl's history ends the id; idtrack bridges (rule `bridge_same_place`) to a gene that overlaps the old place",
     "read the answer as a candidate at the same locus, not as Ensembl's successor", "decided: rule `bridge_same_place`")),
    (r"support_ranks_shared_name", ("L-shared-names", "names", "a name of several genes: `best` keeps the better-supported",
     "", "use the `all` view (the harmonizer reads `best`)", "none (rule `support_ranks_shared_name`: by design)")),
    (r"^open: the form is a string of the truth store", ("T-made-unique-taken", "truth",
     "a made-unique form that is a name of its own (not idtrack's)",
     "a form `X.1`/`X-1` drawn as a copy of X is a string of the truth store (revision `copies_taken`): its copy's truth "
     "does not hold", "-", "the truth: the set's revision `copies_taken`")),
    (r"no place on the reference sequence", ("L-unplaced-copy", "history",
     "a copy on a haplotype or patch without a place on the reference", "", "-", "")),
    (r"^open: ", ("L-open-truth", "history", "the truth leaves the line open", "", "", "")),
]

CHECK_TITLES = {
    "plan": ("L-harmonize-plan", "a harmonized dataset that differs from its documented rules"),
    "unify": ("L-harmonize-merge", "merged objects that differ from their documented rules"),
    "shuffle": ("L-harmonize-order", "a merge that changes with the order of datasets or features"),
    "cache": ("L-harmonize-cache", "a result cache that does not behave as documented"),
    "input": ("L-harmonize-input", "a dataset read otherwise than written"),
}


def classify(reason):
    for pat, row in UND:
        if re.search(pat, reason or ""):
            return row
    return ("L-undetermined-other", "?", "other undetermined answers", re.sub(r"^\w+ \d+: ", "", reason or "")[:120], "", "")


def _row(lim, org, n, units, of, strata, run, info, cls, cause, measured_by):
    commit = (info.get("commit") or "")[:12]
    return dict(limitation=lim[0], group=lim[1], title=lim[2], what=lim[3], organism=org, answers=n, units=units,
                of_answers=of, share=f"{n / of:.5f}" if of else "", family="f", strata=strata, run=run, commit=commit,
                graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                **{"class": cls}, cause=cause, rule="", user_action=lim[4], fixed_by=lim[5], measured_by=measured_by,
                since=commit)


def rows_for(scored):
    """Rows of one scored folder (`<run>/scores_<config>`): causes, undetermined by reason, and the run's checks."""
    info = json.load(open(os.path.join(scored, "run.json")))
    org = info.get("organism")
    cfg = os.path.basename(scored)[len("scores_"):]
    out = []
    if os.path.exists(os.path.join(scored, "causes.tsv.gz")):
        out += [{k: str(v) for k, v in r.items()} for r in
                LIM.build("f", [scored], LIM.load_map(os.path.join(HERE, "limitations_map_f.tsv")))]
    totals = collections.Counter()
    for r in csv.DictReader(open(os.path.join(scored, "report.tsv")), delimiter="\t"):
        if r["view"] == "*":
            totals[r["stratum"]] = int(r["n"])
    agg = collections.defaultdict(lambda: [0, set(), set()])
    p = os.path.join(scored, "undetermined.tsv.gz")
    if os.path.exists(p):
        with gzip.open(p, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                x = agg[classify(r["reason"])]
                x[0] += 1
                x[1].add(r["unit"])
                x[2].add(r["stratum"])
    for lim, (n, units, strata) in sorted(agg.items(), key=lambda kv: -kv[1][0]):
        out.append(_row(lim, org, n, len(units), sum(totals[s] for s in strata), ",".join(sorted(strata)), LIM.where(scored),
                        info, "e", "undetermined", f"f: undetermined.tsv.gz of {LIM.where(scored)}"))
    cache = os.path.join(os.path.dirname(scored), cfg, "cache.json")
    if os.path.exists(cache):  # the cost of reading the result cache, as measured (idtrack's documented cost)
        c = json.load(open(cache))
        again, first = c["again"]["seconds"], c["again"]["construct_seconds"]
        lim = ("L-harmonize-cache-reuse", "performance", "reading the harmonizer's result cache costs nearly a conversion",
               f"constructing the harmonizer again over its valid result cache took {again} s against {first} s for "
               f"converting every feature (share: the cache read's time over the conversion's), the cost idtrack's "
               f"documentation states (usage, \"Known Costs of This Version\")",
               "keep the harmonizer, or write the merged AnnData to a file, rather than constructing it again",
               "none (the release's documented cost)")
        n = len(c["headers"])
        row = _row(lim, org, n, "", n, "cache", LIM.where(scored), info, "-", "cache_reuse",
                   f"f: cache.json of {LIM.where(os.path.join(os.path.dirname(scored), cfg))}")
        row["share"] = f"{again / first:.3f}" if first else ""
        out.append(row)
    checks = os.path.join(os.path.dirname(scored), "f_checks.tsv")
    if os.path.exists(checks):
        agg = collections.defaultdict(lambda: [0, 0, set()])
        for r in csv.DictReader(open(checks), delimiter="\t", quoting=csv.QUOTE_NONE):
            if r["config"] not in (cfg, "-"):
                continue
            k = r["check"].split("_")[0]
            agg[k][0] += int(r["violations"])
            agg[k][1] += int(r["n"])
            if int(r["violations"]):
                agg[k][2].add(r["check"])
        for k, (bad, n, which) in sorted(agg.items()):
            if not bad or k not in CHECK_TITLES:
                continue
            lid, title = CHECK_TITLES[k]
            lim = (lid, "harmonization", title, "the checks that do not hold: " + ",".join(sorted(which)), "", "")
            out.append(_row(lim, org, bad, "", n, ",".join(sorted(which)), LIM.where(scored), info, "a", f"check_{k}",
                            f"f: f_checks.tsv of {LIM.where(os.path.dirname(scored))}"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("scored", nargs="+")
    a = ap.parse_args()
    rows = []
    for s in a.scored:
        scored = [s] if os.path.exists(os.path.join(s, "report.tsv")) else \
            sorted(os.path.dirname(p) for p in glob.glob(os.path.join(s, "scores_*", "report.tsv")))
        for x in scored:
            rows += rows_for(x)
    with open(a.out, "w") as fo:
        fo.write("\t".join(LIM.HEADER) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in LIM.HEADER) + "\n")
    print(len(rows), "rows ->", a.out)


if __name__ == "__main__":
    main()
