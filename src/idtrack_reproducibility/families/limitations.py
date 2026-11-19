"""Rows of the known-limitations table (docs/methods.md, "Known limitations") from a family's causes and reports (F0).

Usage: limitations.py --family <f> --scored <scored folder> [<scored folder> ...] [--map <family map .tsv>] --out <tsv>

A limitation is a class of answers that a family finds wrong or undeterminable and that no stage fixes: one row
per (limitation, organism) with its measured size, never typed by hand. Its text comes from a map, cause ->
(`limitation` id, `group`, `title`, `what`, `user_action`, `fixed_by`): F0's map below for the causes of the residual
reading, a family's own map (`--map`, the same columns, its rows replacing F0's) for its causes. A cause whose
`fixed_by` names a stage or decision is written as such (`until`), so that the next tree re-measures it; a cause
without a map row is written with `title` `UNMAPPED` (to be read and mapped, never left so).

Columns: limitation, group, title, what, organism, answers, units, of_answers (the lines of those strata scored on
the tree), share, family, strata, run, commit, graph_digest, class, cause, decision, user_action, fixed_by,
measured_by, since. `run.json` and `report.tsv` of each scored folder give the run, commit, digest and the totals.
"""
import argparse
import collections
import csv
import json
import os
import sys

from idtrack_reproducibility.families import causes as K

COLS = ["limitation", "group", "title", "what", "user_action", "fixed_by"]

# The causes of the reading of every wrong answer of the reference tree (the rules: docs/rules.md).
F0_MAP = {
    "entrezgene_label_rows": ("L-ncbi-labels", "external identifiers", "NCBI gene labels read as NCBI gene ids",
                              "up to release 87 Ensembl's EntrezGene holds NCBI labels beside the numbers; a gene converted to "
                              "EntrezGene at those releases also returns its label", "ask for `NCBI gene (formerly Entrezgene)` "
                              "or filter numeric ids", "decided: rule `one_kind_per_database`"),
    "mirbase_not_in_graph": ("L-mirbase", "names", "miRBase names not in the graph",
                             "miRNA genes whose display name is a miRBase name are not found by that name",
                             "use the Ensembl id or the HGNC/MGI symbol of the miRNA gene", "decided: rule `mirbase_genes`"),
    "external_bridge": ("L-bridges", "history", "no successor in Ensembl's history: a bridge through a shared id",
                        "an id Ensembl's history ends is carried through an external id to a gene that does not overlap it",
                        "`allow_bridges=False`; check the answer's place (explain)", "decided: rule `bridge_same_place`"),
    "support_ranks_shared_name": ("L-shared-names", "names", "a name of several genes: `best` keeps the better-supported",
                                  "`best` returns the gene more databases name so; the others are in `all`",
                                  "use the `all` view for names of several genes", "none (rule `support_ranks_shared_name`: by design)"),
    "rfam_family_label": ("L-rfam-labels", "names", "RFAM family labels name many genes",
                          "an RFAM family label attaches the string to genes whose name it is not",
                          "use the gene's own symbol", "decided: rule `rfam_display_label`"),
    "synonym_spelling_merged": ("L-spellings", "names", "two genes with one string in two letter cases",
                                "the exact step cannot tell `Tc` from `tc`", "type the exact spelling (decided: rule `spelling_per_edge`)",
                                "decided: rule `spelling_per_edge`"),
    "one_spelling_per_db_release": ("L-spellings", "external identifiers", "two spellings of one value",
                                    "one database writes one value in two spellings at one release; one is returned",
                                    "compare values case-insensitively", "decided: rule `spelling_per_edge`"),
    "event_version_contradicts_tables": ("L-event-versions", "history", "events whose version contradicts the tables",
                                         "Ensembl's event names a version the id did not have; the event was left out",
                                         "none", "decided: rule `event_version`"),
    "path_through_resurrection": ("L-resurrection", "history", "ids retired and brought back",
                                  "an id retired and given back later is one gene throughout (Ensembl's events end it)",
                                  "check the explained path", "none (rule `resurrection`: by design)"),
    "copies_apart_before_82": ("L-copies-before-82", "history", "copies on patches before release 82",
                               "before 82 the GTF does not mark patch and haplotype copies; they are reported apart but "
                               "cannot be checked", "read `alternative_locus_genes`", "none (Ensembl's data)"),
    "name_in_other_source": ("L-other-name-sources", "names", "names from other name sources",
                             "the string names the gene in another source (UniProt gene name, NCBI label, WikiGene)",
                             "pass `from_database`", "none (Ensembl's data; undeterminable)"),
    "variant_lacks_broad": ("L-other-name-sources", "names", "case and separator variants of names from other sources",
                            "as L-other-name-sources for spelling variants", "type the exact spelling", "none"),
    "gtf_lacks_core_gene": ("L-truth-gtf-short", "truth", "genes of a release missing from its GTF",
                            "human release 48 (NCBI36): the GTF lists 32,354 genes, the core gene table 32,561; idtrack "
                            "answers genes of the core that the truth's gene list lacks", "none needed (idtrack follows "
                            "the core)", "none (the truth: presence from the core table where the GTF is short)"),
    "decision14_reverses_kind": ("L-newer-versions", "history", "a version newer than the graph is answered",
                                 "a version above every version of a living id is answered with its gene, flagged "
                                 "`version_newer_than_graph` (rule `version_newer_than_graph`)", "read `match`", "none (rule `version_newer_than_graph`: by design)"),
}


def load_map(path):
    m = dict(F0_MAP)
    if path:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                m[r["cause"]] = tuple(r[c] for c in COLS)
    return m


HEADER = ["limitation", "group", "title", "what", "organism", "answers", "units", "of_answers", "share", "family",
          "strata", "run", "commit", "graph_digest", "class", "cause", "rule", "user_action", "fixed_by", "measured_by",
          "since"]


def where(path):
    """A folder as the table names it: its path in the data root (never a machine path)."""
    from idtrack_reproducibility import config as CFG
    root = os.path.realpath(CFG.get().root)
    full = os.path.abspath(path)
    for base in (os.path.abspath(CFG.get().root), root):
        if full.startswith(base + os.sep):
            return os.path.relpath(full, base)
    return os.path.basename(full.rstrip(os.sep))


def build(family, scored, lim_map, und=False):
    """The rows of one family's limitations from its scored folders (each with run.json, report.tsv and causes)."""
    import gzip
    rows = []
    for d in scored:
        info = json.load(open(os.path.join(d, "run.json")))
        org = info["organism"]
        rep = {}
        with open(os.path.join(d, "report.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["view"] == "*":
                    rep[r["stratum"]] = int(r["n"])
        suffix = "_und" if und else ""
        agg = collections.defaultdict(lambda: [0, set(), set()])
        with gzip.open(os.path.join(d, f"causes{suffix}.tsv.gz"), "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                lim = lim_map.get(r["cause"])
                key = (lim[0] if lim else f"UNMAPPED-{r['cause']}", r["cause"])
                x = agg[key]
                x[0] += 1
                x[1].add(r["unit"])
                x[2].add(r["stratum"])
        commit = (info.get("commit") or "")[:12]
        for (lid, cause), (n, units, strata) in sorted(agg.items()):
            lim = lim_map.get(cause) or (lid, "", "UNMAPPED", "", "", "")
            of = sum(rep.get(s, 0) for s in strata)
            cls, rule, _text = K.CAUSE_CLASS.get(cause, ("?", "-", ""))
            rows.append(dict(limitation=lim[0], group=lim[1], title=lim[2], what=lim[3], organism=org, answers=n,
                             units=len(units), of_answers=of, share=f"{n / of:.5f}" if of else "", family=family,
                             strata=",".join(sorted(strata)), run=where(d), commit=commit,
                             graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                             **{"class": cls}, cause=cause, rule=rule, user_action=lim[4], fixed_by=lim[5],
                             measured_by=f"{family}: causes{suffix}.tsv.gz of {where(d)}", since=commit))
    return rows


def write(rows, out):
    with open(out + ".tmp", "w") as fo:
        fo.write("\t".join(HEADER) + "\n")
        for r in rows:
            fo.write("\t".join(str(r[h]).replace("\t", " ") for h in HEADER) + "\n")
    os.replace(out + ".tmp", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--scored", nargs="+", required=True)
    ap.add_argument("--map", default=None)
    ap.add_argument("--und", action="store_true", help="from causes_und (the undetermined answers)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = build(a.family, a.scored, load_map(a.map), a.und)
    write(rows, a.out)
    print(f"{len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
