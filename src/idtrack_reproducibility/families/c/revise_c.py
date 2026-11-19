"""Revision sets of family C: corrections of the judge where its truth is wrong, published under the qids of the lines
they change (a published set is never rewritten; `lines.publish(revises=...)`), with the evidence in the docstring.

Usage: revise_c.py <organism> --correction <name> --set <new set name> [--dry-run]

Corrections:
- `rfam_label_accession`: an RFAM label line (`c4_rfam_label_*`, `c7_rfam_label`) whose label is written as an RFAM
  accession (`RF00003`): Ensembl writes the entry's accession as its label at some releases (core xrefs: `RFAM` rows
  with `display_label` = `dbprimary_acc`), so the string is the accession of the entry, whose question the
  accession's lines ask (`c4_rfam_acc_*`, `c7_rfam`); the label line expected only the genes of the label rows and
  is left out (`excluded`, with the reason). Evidence: the core xref extracts (`core_xref/<release>.tsv.gz`).
- `other_database_broad`: a line without `from_database` whose string is also a value (accession or label) of another
  Ensembl database, or the GTF's gene name, of other genes at the release read (the core xref extracts and the GTF of
  that release): without `from_database` the string names those genes too, so an answer within them is not refuted;
  they go to `broad` (their walks to the target for a line read at an older release), and to `expect` of a `safe`
  line (the string names them). The corrected judge's `names_broad` for family C.
"""
import argparse
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.c import truth_c as TC

FAMILY = "c"
ME = os.path.abspath(__file__)


def log(*a):
    import time
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rfam_label_accession(org, lines, ctx):
    """RFAM label lines whose label is the entry's accession: left out (the accession's lines ask the question)."""
    for rec in lines:
        if "rfam_label_accession" in (rec.get("revised_by") or []):
            continue
        k = rec["kind"]
        if (k.startswith("c4_rfam_label") or k == "c7_rfam_label") and re.match(r"^RF\d{5}$", str(rec["id"])):
            new = dict(rec)
            new["excluded"] = ("the label is written as the RFAM accession (Ensembl's label = accession at some "
                               "releases): the accession's lines ask this question")
            new["revised_by"] = list(rec.get("revised_by") or []) + ["rfam_label_accession"]
            yield new


class Carriers:
    """value -> genes (base) of every database of the core xref extract at a release and the GTF's gene names."""

    def __init__(self, org):
        self.org, self.cache = org, {}

    def at(self, rel):
        if rel not in self.cache:
            if len(self.cache) > 3:
                self.cache.clear()
            from idtrack_reproducibility.harness import harness_lib as H
            m = collections.defaultdict(set)
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(self.org, rel):
                for v in (acc, lab):
                    if v:
                        m[v].add(gid)
            for g, r in (H.gtf_genes(self.org, rel) or {}).items():
                if r.get("gene_name"):
                    m[r["gene_name"]].add(g)
            self.cache[rel] = m
        return self.cache[rel]


def other_database_broad(org, lines, ctx):
    """A string that other databases or the GTF give to other genes at the release read: those genes (walked to the
    target) in `broad` (in `expect` of a `safe` line)."""
    from idtrack_reproducibility.families.a import gen_a as GA
    G = ctx.setdefault("G", GA.GenA(org, 115))
    tr = G.tr
    car = ctx.setdefault("car", Carriers(org))
    # lines sorted by the release read, so that one release's carriers are read once
    lines = sorted(lines, key=lambda r: (r.get("from_release") or r.get("read_at") or r["to_release"]))
    for rec in lines:
        if rec.get("from_database") or rec.get("final_database") or rec.get("excluded") or rec.get("open") or \
                "other_database_broad" in (rec.get("revised_by") or []):
            continue
        q = str(rec["id"])
        R = rec.get("from_release") or rec.get("read_at") or rec["to_release"]
        t = rec["to_release"]
        exp_b = {C.base_id(x) for x in rec.get("expect") or []}
        gs = car.at(R).get(q, set())
        if not gs and "." in q:
            gs = car.at(R).get(q.rsplit(".", 1)[0], set())  # an accession with its version: its unversioned value
        others = {g for g in gs if g in (tr.genes.get(R) or {})}
        if R == t:
            extra = others - exp_b
            if not extra:
                continue
            ext = set(tr.vset(extra, t))
        else:
            ext = set()
            for g in sorted(others):
                E, B = G.walks(g, R, [t])
                ext |= set(tr.vset(E.get(t, set()) | B.get(t, set()), t))
            ext -= set(rec.get("expect") or [])
            if not ext:
                continue
        new = dict(rec)
        if rec.get("rule") == "safe":
            new["expect"] = sorted(set(rec.get("expect") or []) | ext)
        else:
            new["broad"] = sorted(set(rec.get("broad") or []) | set(rec.get("expect") or []) | ext)
        new["other_carriers"] = sorted(ext)[:50]
        new["revised_by"] = list(rec.get("revised_by") or []) + ["other_database_broad"]
        yield new


CORRECTIONS = {"rfam_label_accession": rfam_label_accession, "other_database_broad": other_database_broad}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--correction", required=True, choices=sorted(CORRECTIONS))
    ap.add_argument("--set", required=True)
    ap.add_argument("--sets", default=None, help="the sets read (default: every published set)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org = a.organism
    sets = a.sets.split(",") if a.sets else None
    names = L.resolve_sets(FAMILY, org, sets, "revised")
    # the latest line of a qid: a revision set's line replaces the base set's (sets in the order of their seq)
    latest, base_sets = {}, []
    for n in L.resolve_sets(FAMILY, org, None, "revised"):
        if L.manifest(FAMILY, org, n).get("revises"):
            for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), n)):
                latest[r["qid"]] = r
        elif n in names:
            base_sets.append(n)
    names = base_sets
    ctx, out = {}, []
    d = L.query_dir(FAMILY, org)
    for n in names:
        p = L._set_file(d, n)
        recs = (latest.get(r["qid"], r) for r in L.iter_set(p)
                if not r["stratum"].startswith("c6") and "_ver_" not in r["stratum"])
        k = len(out)
        out.extend(CORRECTIONS[a.correction](org, recs, ctx))
        log(n, "revised", len(out) - k)
    by = collections.Counter(r["stratum"] for r in out)
    log("revised", len(out), dict(by))
    if not out:
        log("nothing to revise")
        return
    doc = CORRECTIONS[a.correction].__doc__ or ""
    man = L.publish(FAMILY, org, a.set, out, f"revision `{a.correction}` (revise_c.py docstring)", ME,
                    revises=L.resolve_sets(FAMILY, org, None, "revised"), judge="c-" + a.correction, sources=["core_xref/<release>.tsv.gz", "gtf/genes_<release>.tsv.gz"],
                    manifest_extra=dict(correction=a.correction, revision=doc.strip() or a.correction), dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k != "counts"})


if __name__ == "__main__":
    main()
