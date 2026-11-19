"""Family D's known limitations (docs/methods.md, "Known limitations"): F0's rows for the wrong answers (limitations.py, from the
causes; D's map `limitations_map_d.tsv`) plus rows for the answers the truth leaves open, by D1's truth class (the
`open:` reason of undetermined.tsv.gz) and D2/D3's open reasons. Every number comes from the runs; nothing typed by hand.

Usage: limits_d.py --out <tsv> <scored run folder> ...
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
from idtrack_reproducibility.families import limitations as LIM  # noqa: E402

# (pattern on the judge's reason) -> (limitation, group, title, what, user action, fixed by)
UND = [
    (r"several genes of 115", ("L-place-several", "history", "a feature's model holds several genes of the target release",
     "the gene model of the source release overlaps several genes of the target release with none dominant: no one gene "
     "can be told by place", "keep the feature apart, or sum it into the gene the history names (explain)", "none (Ensembl's data)")),
    (r"under half of the model", ("L-place-shrunk", "history", "the one gene at the place holds under half of the old model",
     "the gene of the target release at the old model's place holds less than half of its exon bases (the model shrank "
     "or moved)", "as L-place-several", "none (Ensembl's data)")),
    (r"holds no gene of 115", ("L-place-empty", "history", "no gene of the target release at the old model's place",
     "Ensembl removed the gene and no gene of the target release shares an exon base with its old model", "drop the feature",
     "none (Ensembl's data)")),
    (r"no place on the assembly", ("L-place-unmapped", "history", "the old model does not map to the target assembly",
     "the chain file aligns under half of the model's bases, or its sequence is not in the chain", "-", "none")),
    (r"only on a patch", ("L-patch-source", "history", "an id of the source release only on a patch or alternative locus",
     "GENCODE-based matrices list genes of alternative loci; their place is projected to the reference sequence", "-", "none")),
    (r"no gene of the source release's GTF", ("L-id-not-at-source", "history", "an id the source release does not hold", "", "-", "none")),
    (r"name of no gene in the source release's GTF", ("L-name-not-at-source", "names",
     "a symbol that no gene carries as its Ensembl name at the source release",
     "names of another annotation (RefSeq-based tables of 2013 to 2016, a laboratory's own names): Ensembl's GTF cannot say "
     "which gene they are", "map such a table by its own annotation", "none")),
    (r"made unique", ("L-made-unique", "names", "names made unique by the matrix's maker (X.1, X-1)",
     "the second row of a name its reference gave to several genes: which gene cannot be told", "-", "decided: rule `made_unique_copies`")),
    (r"name of several genes at the source release", ("L-name-several-place", "names", "a name of several genes, not all of them with one gene at the target",
     "", "-", "none")),
    (r"stated release .* fits the list worse", ("L-stated-release", "source inference", "a stated release that fits the list worse than another",
     "the authors' statement and Ensembl's gene lists disagree: the list's ids are genes of another release", "check the stated release", "none")),
    (r"no gene carries the symbol as its Ensembl name", ("L-marker-not-a-name", "names", "a marker string that is no Ensembl name at the list's release",
     "marker lists write aliases, old symbols and protein names (CD45, Ki67): Ensembl's names cannot judge them; the authority's reading is reported beside",
     "convert the list's symbols with HGNC/MGI first, or read idtrack's match kind", "none")),
    (r"history ends for the genes the symbol names", ("L-marker-retired", "names", "a marker whose gene Ensembl retired", "", "-", "rule `bridge_same_place`")),
    (r"within the wider set", ("L-branch", "history", "ids that branched off while the old id lived on", "", "use `all`", "decided: rule `main_line_first`")),
    # the harness's reasons, on D4's items (the harness's lines asked in lists)
    (r"every returned gene overlaps the query's place", ("L-bridge-at-place", "history",
     "no successor in Ensembl's history; idtrack answers a gene at the same place",
     "Ensembl's history ends the id; idtrack bridges (rule `bridge_same_place`) to a gene that overlaps the old place: plausible, "
     "not Ensembl's statement", "read the answer as a candidate at the same locus", "decided: rule `bridge_same_place`")),
    (r"places not comparable|no assembly chain from", ("L-no-chain", "history", "no successor and the place does not map: the answer cannot be judged",
     "the old place has no mapping to the target assembly through Ensembl's chain files", "check the answer's place in the explained path", "")),
    (r"no place on the reference sequence", ("L-unplaced-copy", "history", "a copy on a haplotype or patch without a place on the reference", "", "-", "")),
    (r"a missing-value string \(nan, NA, None, null\) raises MissingValueError", ("L-missing-values", "inputs", "missing-value strings are read as identifiers before rule `missing_values`",
     "`nan`, `NA`, `None`, `null` in a list are looked up as identifiers (d4c2dda: nothing; 71690ab: `nan` -> SCN11A, `NA` -> XK "
     "in pig, by synonyms) until rule `missing_values` raises MissingValueError for them", "remove missing values before converting",
     "decided: rule `missing_values`")),
    (r"write no versions in their GTF", ("T-versionless-releases", "truth", "versions at releases whose GTF writes none (not idtrack's)",
     "pig and human GTFs before 77 (mouse 76) write no versions; idtrack's graph holds them (from the core tables), so whether "
     "such a release belongs to a versioned list's tie cannot be told from the GTFs", "-", "the truth: the core tables' versions")),
    (r"carry the name in more name sources than every expected gene left out", ("L-shared-names", "names", "a name of several genes: `best` keeps the better-supported", "", "use `all`",
     "none (rule `support_ranks_shared_name`: by design)")),
]


def classify(reason):
    for pat, row in UND:
        if re.search(pat, reason or ""):
            return row
    # the judge's reason, without the earlier notes' number of the rule it cites
    return ("L-undetermined-other", "?", "other undetermined answers", re.sub(r"^\w+ \d+: ", "", reason or "")[:120], "", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    with_causes = [r for r in a.runs if os.path.exists(os.path.join(r, "causes.tsv.gz"))]
    rows = [{k: str(v) for k, v in r.items()} for r in
            LIM.build("d", with_causes, LIM.load_map(os.path.join(HERE, "limitations_map_d.tsv")))] if with_causes else []
    cols = list(LIM.HEADER)
    for run in a.runs:
        info = json.load(open(os.path.join(run, "run.json")))
        org = info.get("organism")
        totals = collections.Counter()
        for r in csv.DictReader(open(os.path.join(run, "report.tsv")), delimiter="\t"):
            if r["view"] == "*":
                totals[r["stratum"]] = int(r["n"])
        agg = collections.defaultdict(lambda: [0, set(), set()])
        p = os.path.join(run, "undetermined.tsv.gz")
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                lim = classify(r["reason"])
                x = agg[lim]
                x[0] += 1
                x[1].add(r["unit"])
                x[2].add(r["stratum"])
        for lim, (n, units, strata) in sorted(agg.items(), key=lambda kv: -kv[1][0]):
            of = sum(totals[s] for s in strata)
            row = dict(limitation=lim[0], group=lim[1], title=lim[2], what=lim[3], organism=org, answers=n, units=len(units),
                       of_answers=of, share=f"{n / of:.5f}" if of else "", family="d", strata=",".join(sorted(strata)), run=LIM.where(run),
                       commit=(info.get("commit") or "")[:12],
                       graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                       **{"class": "e", "cause": "undetermined", "rule": "", "user_action": lim[4], "fixed_by": lim[5],
                          "measured_by": f"d: undetermined.tsv.gz of {LIM.where(run)}", "since": (info.get("commit") or "")[:12]})
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
