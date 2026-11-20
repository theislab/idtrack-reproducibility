"""Revision sets of family C: corrections of the judge where its truth is wrong, published under the qids of the lines
they change (a published set is never rewritten; `lines.publish(revises=...)`), with the evidence in the docstring.

Usage: revise_c.py <organism> --correction <name> --set <new set name> [--dry-run] [--redo <published revision set>]

Corrections:
- `rfam_label_accession`: an RFAM label line (`c4_rfam_label_*`, `c7_rfam_label`) whose label is written as an RFAM
  accession (`RF00003`): Ensembl writes the entry's accession as its label at some releases (core xrefs: `RFAM` rows
  with `display_label` = `dbprimary_acc`), so the string is also the entry's accession and names every gene of the
  entry; the label line expected only the genes of the label rows: the entry's other genes go to `broad` (not
  refuted). Evidence: the core xref extracts (`core_xref/<release>.tsv.gz`).
- `other_database_broad`: a line without `from_database` whose string is also a value (accession or label) of another
  Ensembl database, or the GTF's gene name, of other genes at the release read, or, for a line asked without a
  release and without a version, at the target where the string is then read (rule `start_release`; an accession
  with a version is read where it had it, rule `accession_versions`; the core xref extracts and the GTF of that
  release): without `from_database` the string names those genes too, so an answer within them is not refuted;
  they go to `broad` (their walks to the target for a line read at an older release), and to `expect` of a `safe`
  line (the string names them); the copies the line expects apart are not others. The corrected judge's
  `names_broad` for family C.
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
    """RFAM label lines whose label is written as an RFAM accession (`RF00003`): Ensembl writes an entry's accession as
    its label at some releases, so the string is also the entry's accession and names every gene of the entry at the
    release read and at the target, where the string is read when no release is given (core xrefs of database `RFAM`
    with that accession). Those genes (walked to the target from an older release) are not refuted by the truth: they
    go to `broad`; the label rows' genes stay `expect`."""
    from idtrack_reproducibility.families.a import gen_a as GA
    car = ctx["car"] if "car" in ctx else ctx.setdefault("car", Carriers(org))
    for rec in sorted(lines, key=lambda r: (r.get("from_release") or r.get("read_at") or r["to_release"])):
        k = rec["kind"]
        if "rfam_label_accession" in (rec.get("revised_by") or []) or rec.get("excluded") or \
                not (k.startswith("c4_rfam_label") or k == "c7_rfam_label") or not re.match(r"^RF\d{5}$", str(rec["id"])):
            continue
        G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))  # built once, on the first such line
        R = rec.get("from_release") or rec.get("read_at") or rec["to_release"]
        t = rec["to_release"]
        genes = {g for g in car.entry(R).get(str(rec["id"]), set()) if g in (G.tr.genes.get(R) or {})}
        at_t = {g for g in car.entry(t).get(str(rec["id"]), set()) if g in (G.tr.genes.get(t) or {})}
        ext = set(G.tr.vset(at_t, t))  # the entry's genes at the target, where the string is read without a release
        if R == t:
            ext |= set(G.tr.vset(genes, t))
        else:
            for g in sorted(genes):
                E, B = G.walks(g, R, [t])
                ext |= set(G.tr.vset(E.get(t, set()) | B.get(t, set()), t))
        ext -= set(rec.get("expect") or []) | set(rec.get("broad") or []) | set(rec.get("apart") or [])
        if not ext:
            continue
        new = dict(rec)
        new["broad"] = sorted(set(rec.get("broad") or []) | set(rec.get("expect") or []) | ext)
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

    def entry(self, rel):
        """accession -> genes (base) of database `RFAM` at a release (core xref extract)."""
        key = ("rfam", rel)
        if key not in self.cache:
            m = collections.defaultdict(set)
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(self.org, rel):
                if db == "RFAM" and acc:
                    m[acc].add(gid)
            self.cache[key] = m
        return self.cache[key]


def other_database_broad(org, lines, ctx):
    """A string that other databases or the GTF give to other genes at the release read: those genes (walked to the
    target) in `broad` (in `expect` of a `safe` line)."""
    from idtrack_reproducibility.families.a import gen_a as GA
    G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))
    tr = G.tr
    car = ctx["car"] if "car" in ctx else ctx.setdefault("car", Carriers(org))
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
        at_t = set()
        # without a release an unversioned string is read at the target where it is held (rule `start_release`); a
        # version of an accession is read where it had that version (rule `accession_versions`), not here
        if rec.get("from_release") is None and R != t and not re.search(r"\.\d+$", q):
            gt = car.at(t).get(q, set()) or (car.at(t).get(q.rsplit(".", 1)[0], set()) if "." in q else set())
            at_t = set(tr.vset({g for g in gt if g in (tr.genes.get(t) or {})}, t))
        if R == t:
            extra = others - exp_b
            ext = set(tr.vset(extra, t))
        else:
            ext = set(at_t)
            for g in sorted(others):
                E, B = G.walks(g, R, [t])
                ext |= set(tr.vset(E.get(t, set()) | B.get(t, set()), t))
        # the line's own genes and the copies it expects reported apart (rule `judge_copies_apart`: a copy on a patch
        # or haplotype that carries the value is the line's, answered apart, never inside the answer) are no others
        ext -= set(rec.get("expect") or []) | set(rec.get("apart") or [])
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
    ap.add_argument("--redo", default=None, help="a published revision set of this correction, made again (the correction "
                    "fixed since): every line of it is corrected anew from its line before the set, and the lines that "
                    "differ from the published ones are published as --set")
    a = ap.parse_args()
    org = a.organism
    if a.redo:
        return redo(org, a)
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


def redo(org, a):
    """`--redo`: a published revision set of the correction made again from the lines before it (main's docstring)."""
    names = L.resolve_sets(FAMILY, org, None, "revised")
    before = names[:names.index(a.redo)]
    published = {r["qid"]: r for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), a.redo))}
    prior = L.load_qids(FAMILY, org, set(published), sets=before)
    for r in prior.values():
        r.pop("set", None)
    again = {r["qid"]: r for r in CORRECTIONS[a.correction](org, list(prior.values()), {})}
    keys = ("expect", "broad", "apart", "rule", "open", "revised_by", "other_carriers")
    out = []
    for q, old in published.items():
        new = again.get(q, prior[q])  # where the correction no longer applies, the line before it is the truth again
        if any((new.get(k) or None) != (old.get(k) or None) for k in keys):
            out.append(new)
    log("redo", a.redo, "lines", len(published), "differ", len(out),
        dict(collections.Counter((r["stratum"], q in again) for q, r in ((r["qid"], r) for r in out))))
    if not out:
        return
    doc = CORRECTIONS[a.correction].__doc__ or ""
    man = L.publish(FAMILY, org, a.set, out, f"revision `{a.correction}` made again ({a.redo}; revise_c.py docstring)", ME,
                    revises=names, judge="c-" + a.correction, sources=["core_xref/<release>.tsv.gz", "gtf/genes_<release>.tsv.gz"],
                    manifest_extra=dict(correction=a.correction, redo=a.redo, revision=doc.strip() or a.correction),
                    dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k != "counts"})


if __name__ == "__main__":
    main()
