"""Family B's set B9: every display name at its own release, by the database Ensembl filed it under, from Ensembl's
core databases only.

Usage: gen_b9.py --organism <o> [--dry] [--releases r,r,...]

**B9** (`b9_v1`): for every release R of the window whose core database is read (`core_display/<R>.tsv.gz`,
`fetch_b display`: each gene's display cross-reference, `gene.display_xref_id` -> `xref` -> `external_db`), every
label that is the display name of a gene of R's GTF (`truth_b.clean`: no Ensembl-id-like string, no number), asked
at R (`from_release=R`, `to_release=R`, kind `b9_display_at_release`, `strategy: both`). The question: which gene is
this name at the release it is a name of.

- **expect**: the genes of R's GTF on the primary assembly whose display cross-reference at R is the name (where only
  copies on patches carry it, those copies), compared as stable ids without versions (`cmp: base`: the question is
  which gene, not which version; the broad set holds genes the GTF lacks, whose versions the GTF cannot give);
- **apart**: copies on patches (where the GTF of R marks them, human and mouse from 82; rule `judge_copies_apart`), and
  before 82, where the GTFs mark none, the carriers on haplotype and patch sequences by Ensembl's assembly exceptions
  (the corrected judge's `copies_apart_early`), and the Y copy of a pseudoautosomal gene whose X gene carries the name
  (rule `par_shared_names`, decided and built in the trees measured): the accuracy harness's revisions themselves
  (`make_queries.REVISIONS`), applied to each release's lines; the stratum follows the genes expected after them;
- a display name that is an LRG record's (the carrier an LRG gene, `LRG_<n>_g<k>`) expects the Ensembl gene the
  record describes (rule `lrg_answers_gene`, decided): the gene Ensembl's core database of the release links to the
  LRG gene by its cross-reference of database `ENS_LRG_gene` (`core_xref/<R>.tsv.gz`), where it does;
- **rule `set`**, for a name of one gene and for a name several genes carry at once alike: under `best` a non-empty
  part of the carriers is right (each carries the name at R; which of them the support ranks first, rule
  `support_ranks_shared_name`, the display records cannot decide), under `all` exactly the carriers;
- **broad** (the corrected judge's `names_broad`): the other genes that carry the string at R: as a display name off
  R's GTF (genes of the core database the GTF lacks: haplotypes before 82, cynomolgus 112), as a label of a name
  database of R's core cross-references (`truth_b.NAME_DB`: HGNC, MGI, VGNC, UniProt gene name, NCBI gene, WikiGene,
  Vega, clone-based names), or as a synonym of a cross-reference of R (`core_synonyms/<R>.tsv.gz`); an answer within
  them is undetermined, outside them wrong.

Strata (the modest set samples each): `b9_<era>_<db>` for a name of one gene, `<db>` the database of its display
cross-reference (lower case, `[a-z0-9_]`), and `b9_shared_<era>` for a name several genes carry at once (the
databases of their cross-references in the line's `db`). The **era** of a release, per organism, from the display
records alone: `late` from the first release B such that at every release from B on every database filing a display
name also files display names at the last release of the window (115); `early` before B: the releases at which
Ensembl filed display names under databases it later renamed or dropped (`HGNC_curated_gene`, `Uniprot/SPTREMBL`,
`Uniprot_genename`, `Clone_based_vega_gene`, ...). The boundary is in the manifest (`era_from`).

qid `b9|<name>|<R>`; unit the name. A published set is never rewritten (F0's `lines.publish`).
"""
import argparse
import collections
import gzip
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.b import truth_b as TB
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M

FAMILY = "b"
SET = "b9_v1"
T = TB.T
RE_LRG_GENE = re.compile(r"^LRG_\d+_g\d+$")
LRG_NOTE = ("rule `lrg_answers_gene`: the display name is an LRG record's; the answer is the Ensembl gene the record "
            "describes, the gene that Ensembl's core database of the release links to the LRG gene by its cross-reference "
            "of database `ENS_LRG_gene` (LRG display in Ensembl), whose accession is the LRG gene id (core_xref/<R>.tsv.gz)")
# the accuracy harness's revisions applied to B9's lines (module docstring), in this order
REVISIONS = ["copies_apart_early", "par_shared_names"]


class Places:
    """What the harness's revisions `copies_apart_early` and `par_shared_names` read of their context: the GTF places
    of a release (`gtf`) and a cache (`_cached`, which holds Ensembl's assembly exceptions)."""

    def __init__(self, org):
        self.org, self._c = org, {}

    def _cached(self, key, make):
        if key not in self._c:
            self._c[key] = make()
        return self._c[key]

    def gtf(self, rel):
        return self._cached(("gtf", rel), lambda: H.gtf_genes(self.org, rel) or {})

    def release_done(self):
        for k in [k for k in self._c if isinstance(k, tuple) and k[0] == "gtf"]:
            del self._c[k]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def db_code(db):
    """A database name as part of a stratum name (lower case, `[a-z0-9_]`)."""
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]", "_", db.lower())).strip("_") or "none"


def era_boundary(dbs_by_release, last):
    """The first release B from which on every database filing a display name also files display names at `last`
    (module docstring); releases before B are `early`. `dbs_by_release`: {R: {db_name}}."""
    today = dbs_by_release.get(last, set())
    B = last
    for R in sorted(dbs_by_release, reverse=True):
        if dbs_by_release[R] <= today:
            B = R
        else:
            break
    return B


def synonyms_at(org, R):
    """synonym -> {gene ids} at R (Ensembl's synonyms of every cross-reference of a name database), or None."""
    return TB.core_synonyms(org, R)


class GenB9:
    def __init__(self, org, releases=None):
        self.org = org
        rels = [r for r in TB.core_display_releases(org) if r <= T and r in set(H.gtf_releases(org))]
        self.rels = [r for r in rels if r in releases] if releases else rels
        self.dbs = {}
        for R in rels:
            self.dbs[R] = {db for db, _ in TB.core_display(org, R).values()}
        self.era_from = era_boundary(self.dbs, max(rels)) if rels else None
        self.counts = collections.Counter()
        self.excluded = collections.Counter()
        self.gtf_disagrees = 0
        self.revised = collections.Counter()
        self.lrg_lines = 0
        self._lrg = {}
        self.places = Places(org)
        self.extra = {}

    def lrg_links(self, R):
        """{LRG gene id: {Ensembl gene ids}} of release R from Ensembl's core cross-references (`ENS_LRG_gene`)."""
        if R not in self._lrg:
            from idtrack_reproducibility.families.e.revise_e import lrg_genes
            self._lrg = {R: lrg_genes(self.org, R) or {}}
        return self._lrg[R]

    def era(self, R):
        return "early" if R < self.era_from else "late"

    def release_lines(self, R):
        org = self.org
        disp = TB.core_display(org, R)
        G = H.gtf_genes(org, R)
        cx = TB.core_names(org, R)
        names_r = cx[0] if cx else {}
        syn = synonyms_at(org, R) or {}
        carriers = collections.defaultdict(set)
        for g, (db, lab) in disp.items():
            if not TB.clean(lab) or lab == g:
                self.excluded[f"R{R} not a name string"] += 1
                continue
            carriers[lab].add(g)
        era = self.era(R)
        vset = sorted
        recs = []
        for v in sorted(carriers):
            allc = carriers[v]
            ing = {g for g in allc if g in G}
            if not ing:
                self.excluded["display carriers not in the GTF of the release"] += 1
                continue
            prim = {g for g in ing if G[g]["region"] != "patch"}
            E, apart = (prim, ing - prim) if prim else (ing, set())
            self.gtf_disagrees += sum(1 for g in E if G[g].get("gene_name") != v)
            dbs = sorted({disp[g][0] for g in E})
            st = f"b9_shared_{era}" if len(E) > 1 else f"b9_{era}_{'patchonly' if not prim else db_code(dbs[0])}"
            note = None
            if E and all(RE_LRG_GENE.match(g) for g in E):  # rule `lrg_answers_gene`: the Ensembl gene the record describes
                links = self.lrg_links(R)
                to = {C.base_id(x) for g in E for x in links.get(g, ())}
                if to:
                    E, note = to, LRG_NOTE
            other = (allc - ing) | set(names_r.get(v, ())) | set(syn.get(v, ()))
            B = E | apart | other
            rec = dict(qid=f"b9|{v}|{R}", family=FAMILY, kind="b9_display_at_release", stratum=st, unit=v, id=v,
                       from_release=R, to_release=R, final_database=None, strategy="both", rule="set", cmp="base",
                       expect=vset(E), broad=vset(B) if B != E else None, db=",".join(dbs), era=era,
                       name_src="core_display+core_xref" + ("+core_synonyms" if syn else ""))
            if apart:
                rec["apart"] = vset(apart)
            if note:
                rec["note"] = note
                self.lrg_lines += 1
            recs.append(rec)
        for rec in M.apply_revisions(org, recs, self.places, names=REVISIONS):
            if rec.get("revised_by"):
                self.revised["+".join(rec["revised_by"])] += 1
                E = rec["expect"]
                if len(E) > 1:
                    rec["stratum"] = f"b9_shared_{era}"
                else:
                    rec["stratum"] = f"b9_{era}_{'patchonly' if G[E[0]]['region'] == 'patch' else db_code(disp[E[0]][0])}"
                    rec["db"] = disp[E[0]][0]
            self.counts[rec["stratum"]] += 1
            yield rec
        self.places.release_done()

    def lines(self):
        for R in self.rels:
            t0, n0 = time.time(), sum(self.counts.values())
            yield from self.release_lines(R)
            log(self.org, "B9 release", R, self.era(R), "lines", sum(self.counts.values()) - n0, f"{time.time() - t0:.0f}s")
        self.extra.update(era_from=self.era_from, excluded=dict(self.excluded), gtf_name_differs=self.gtf_disagrees,
                          revised=dict(self.revised), revisions=REVISIONS, lrg_lines=self.lrg_lines,
                          display_databases={str(R): sorted(d) for R, d in sorted(self.dbs.items())})


def sources(org, rels):
    out = []
    for R in rels:
        for sub in (f"core_display/{R}.tsv.gz", f"core_xref/{R}.tsv.gz", f"core_synonyms/{R}.tsv.gz", f"gtf/genes_{R}.tsv.gz"):
            if os.path.exists(os.path.join(H.truth_dir(org), sub)):
                out.append(f"{org}/{sub}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--releases", default=None, help="only these releases (a dry run)")
    a = ap.parse_args()
    rels = [int(x) for x in a.releases.split(",")] if a.releases else None
    if rels and not a.dry:
        raise SystemExit("--releases is for dry runs: a published set holds every release")
    g = GenB9(a.organism, rels)
    log(a.organism, "releases", g.rels[0] if g.rels else None, "to", g.rels[-1] if g.rels else None, "era from", g.era_from)
    lines = list(g.lines()) if a.dry else g.lines()
    if a.dry:
        print(dict(n=len(lines), counts=dict(sorted(g.counts.items())), **g.extra))
        return
    man = L.publish(FAMILY, a.organism, SET, lines, "B9: every display name at its own release (Ensembl's display "
                    "cross-references), by the database it is filed under and the release era", __file__,
                    sources=sources(a.organism, g.rels), manifest_extra=g.extra, log=log)
    log("published", SET, man["n"], man["counts"])


if __name__ == "__main__":
    main()
