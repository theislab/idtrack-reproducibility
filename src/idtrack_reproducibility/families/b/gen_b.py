"""Family B's exhaustive sets B1 (with B2's and B3's strata) and B8, per organism, from Ensembl's tables only.

Usage: gen_b.py --organism <o> [--sets b1,b8] [--dry] [--limit-releases N]

**B1, every display name at every release** (`b1_v1`): every string that is a gene's display name in the GTF of a
release R of the window (`truth_b.clean`: no Ensembl-id-like string, no number), read at R (`from_release=R`):

- to R (kind `b1_name_at_release`): the genes carrying it at R (`set`, or `set_exact` for several; rule `support_ranks_shared_name`'s
  `support` from BioMart where the release has it); copies on patches apart (rule `judge_copies_apart`) where the GTF marks them;
  `broad` (the harness's current name kinds): the genes that carry the string at R in another name source (BioMart
  or the core xrefs of R) or as a synonym (BioMart); where R's synonyms are not read (no BioMart table), the genes that
  carry it as a synonym at any BioMart release, and where R's name sources are not read at all, the genes that carry
  it as a name at any release (the reach, `name_src` says which); a name only patch copies carry is its own stratum;
- to 115 (kind `b1_name_to_current`, R < 115): the strict walks (living ids) of the genes carrying it at R; `broad`:
  their broad walks (with ids brought back, rule `resurrection`) and those of the genes carrying it at R in another name
  source (the corrected judge `names_broad`: names, not synonyms, rule `names_before_synonyms`); where the history ends, the genes that carry
  it at 115 in any source, and the `locus` rule for a name of one gene (with the corrected judge's places across assemblies
  and exceptions); one gene that split is judged `set` (the corrected judge `split_name_set`); copies apart at 115;
- default (kind `b1_name_default`, no `from_release`, to 115): every display name of 115 (rule `start_release`: read at 115).

Strata: `b1_atR_<c>`, `b1_at115_<c>`, `b1_to115_<c>`, `b1_default_<c>` with `<c>` = `unique` (one gene), `shared`
(several: **B3**), `patchonly`; **B2 moved names**: a to-115 line whose string is carried at 115 by a gene that is no
descendant (broad walk) of its genes at R is stratum `b2_moved_from` (the walk of the old gene is the truth), and the
default line of a string that ever moved is `b2_moved_default` (the gene carrying it now). Unit: the string.

**B8, names to databases over time** (`b8_v1`): every gene with a display name at a release with a table of its
organism's name database (human `HGNC Symbol`, mouse `MGI Symbol`, others `UniProtKB Gene Name`, as the travel
matrix chooses), asked at R with that `final_database` (`cmp: exact`): the values of BioMart's table at R where the
archive answers, else the labels of Ensembl's core xrefs of R (gene-level xrefs of that database). Unit: the gene.

Every line `strategy: both`. A published set is never rewritten (F0's `lines.publish`).
"""
import argparse
import collections
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common  # noqa: F401  paths
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.b import truth_b as TB
from idtrack_reproducibility.harness import make_queries as M

FAMILY = "b"
T = TB.T
NAME_DB = {"homo_sapiens": ("HGNC Symbol", "hgnc", ("HGNC",)), "mus_musculus": ("MGI Symbol", "mgi", ("MGI",))}
OTHER_DB = ("UniProtKB Gene Name", "uniprot_gn", ("Uniprot_gn",))
LOCUS_REVISIONS = ["places", "exceptions", "exceptions_span", "core_places"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def rule_of(n):
    return "set_exact" if n > 1 else "set"


def apart_split(nb, genes, R):
    """(expected, apart): rule `judge_copies_apart`'s copies apart where the GTF of R marks patches."""
    if not nb.marks_patches(R):
        return set(genes), set()
    prim = {g for g in genes if nb.region(g, R) != "patch"}
    patch = set(genes) - prim
    return (prim, patch) if prim and patch else (set(genes), set())


class GenB1:
    def __init__(self, nb, limit_releases=None):
        self.nb, self.tr = nb, nb.tr
        self.rels = nb.rels[-limit_releases:] if limit_releases else nb.rels
        self.excluded = collections.Counter()
        self.moved = set()
        self.locus_lines = 0
        self.counts = collections.Counter()
        self.extra = {}  # the manifest's own fields, filled when the lines are all written

    def any_at(self, R, v):
        """Genes that carry v at R as display name, in another name source or as a synonym."""
        nb = self.nb
        out = set(nb.display(R).get(v, ())) | set(nb.sources(R)[0].get(v, ()))
        syn = nb.synonyms(R)
        if syn:
            out |= syn.get(v, set())
        return out & nb.gset(R)

    def projected_hgnc(self, R, v, genes):
        """Mouse: a display name that is an HGNC symbol of the genes and no name in another source at R."""
        if self.nb.org != "mus_musculus":
            return False
        src, origin, non_hgnc, _ = self.nb.sources(R)
        if origin is None or not (src.get(v, set()) & genes):
            return False
        if non_hgnc.get(v, set()) & genes:
            return False
        syn = self.nb.synonyms(R)
        return not (syn and syn.get(v, set()) & genes)

    def string_lines(self, v, R, prefix="b1", strata=None, at=True, to=True):
        """The lines of string v read at R (its genes from the GTF of R): to R (`at`) and to 115 (`to`, R < 115), as
        B1 judges them; the release's tables must be the ones loaded (releases are walked in order). `strata`: a
        function (default stratum, moved) -> stratum, for the families of lines that reuse B1's judge (B4)."""
        nb, tr = self.nb, self.tr
        disp, (src, origin), syn, G = nb.display(R), nb.sources(R)[:2], nb.synonyms(R), nb.genes(R)
        GS = nb.gset(R)
        name_src = (origin or "none+reach") + ("+synonyms_reach" if syn is None else "+synonyms")
        allg = disp.get(v, set()) & GS
        if not allg:
            return []
        prim = {g for g in allg if G[g]["region"] != "patch"}
        E, apart = (prim, (allg - prim)) if prim else (allg, set())
        cls = "patchonly" if not prim else ("shared" if len(E) > 1 else "unique")
        if self.projected_hgnc(R, v, allg):
            self.excluded[f"R{R} projected HGNC symbol (the corrected judge)"] += 1
            return []
        # the other genes that carry v at R: as a name (another name source) and as a synonym; where R's name sources
        # or its synonyms are not read, the genes that carry v so at any release (the reach)
        names_r = set(src.get(v, ())) | (set(nb.ever("names").get(v, ())) if origin is None else set())
        syn_r = set(syn.get(v, ())) if syn is not None else set(nb.ever("synonyms").get(v, ()))
        out = []
        if at:
            B = (E | apart | names_r | syn_r) & GS
            sup = None
            if len(E) > 1:
                s = {tr.vid(g, R): nb.support(R, g, v) for g in E}
                if None not in s.values() and len(set(s.values())) > 1:
                    sup = s
            st = f"{prefix}_at115_{cls}" if R == T else f"{prefix}_atR_{cls}"
            st = strata(st, False) if strata else st
            out.append(self.line(f"{prefix}|{v}|{R}|{R}", f"{prefix}_name_at_release", st, v, R, R, E, B, apart,
                                 tr.cmp_for(R), name_src, sup))
        if to and R < T:
            disp_T, G_T = nb.display(T), nb.genes(T)
            Es, Bw = nb.walks(allg, R)
            cur = {g for g in disp_T.get(v, ()) if G_T[g]["region"] != "patch"}
            moved = bool(cur - Bw - Es)  # a gene carries the name now that is no descendant of its genes at R
            Bs = set(Bw) | nb.walks(names_r & GS, R)[1]  # the corrected judge names_broad: names, not synonyms (rule `names_before_synonyms`)
            loc = None
            if not Es:
                Bs |= self.any_at(T, v)
                if len(allg) == 1:
                    loc = tr.locus(next(iter(allg)), R)
            Et, At = apart_split(nb, Es, T)
            rule = "set" if len(allg) == 1 else rule_of(len(Et))
            if moved:
                self.moved.add(v)
            st = f"{prefix}_to115_{cls}" if prefix != "b1" else ("b2_moved_from" if moved else f"b1_to115_{cls}")
            st = strata(st, moved) if strata else st
            sup = None
            if rule == "set_exact" and nb.support(R, next(iter(allg)), v) is not None:
                s = {}
                for g in allg:
                    n = nb.support(R, g, v) or 0
                    for t in nb.walk(g, R)[0]:
                        if t in Et:
                            k = tr.vid(t, T)
                            s[k] = max(s.get(k, 0), n)
                for t in Et:
                    s.setdefault(tr.vid(t, T), 0)
                if len(set(s.values())) > 1:
                    sup = s
            ln = self.line(f"{prefix}|{v}|{R}|{T}", f"{prefix}_name_to_current", st, v, R, T, Et, Bs | Et | At, At,
                           "versioned", name_src, sup, rule=rule, locus=loc)
            if loc is not None:  # the corrected judge's places across assemblies and exceptions (the harness's own revisions)
                ln = M.apply_revisions(nb.org, [ln], nb.C, names=LOCUS_REVISIONS)[0]
                self.locus_lines += 1
            out.append(ln)
        return out

    def lines(self):
        nb, tr = self.nb, self.tr
        nb.ever()
        disp_T, G_T = nb.display(T), nb.genes(T)
        for R in self.rels:
            t0 = time.time()
            n0 = sum(self.counts.values())
            disp = nb.display(R)
            for v in sorted(disp):
                if not TB.clean(v):
                    self.excluded[f"R{R} not a name string"] += 1
                    continue
                yield from self.string_lines(v, R)
            nb._walk.clear()
            log(nb.org, "B1 release", R, "lines", sum(self.counts.values()) - n0, f"{time.time() - t0:.0f}s")
        # default lines: every display name of 115, read at 115 (rule `start_release`)
        cmpT = tr.cmp_for(T)
        origin_T = nb.sources(T)[1]
        for v in sorted(disp_T):
            if not TB.clean(v):
                continue
            allg = disp_T[v] & set(G_T)
            if self.projected_hgnc(T, v, allg):
                continue
            prim = {g for g in allg if G_T[g]["region"] != "patch"}
            E, apart = (prim, allg - prim) if prim else (allg, set())
            cls = "patchonly" if not prim else ("shared" if len(E) > 1 else "unique")
            B = self.any_at(T, v) | E | apart
            st = "b2_moved_default" if v in self.moved else f"b1_default_{cls}"
            yield self.line(f"b1|{v}|default|{T}", "b1_name_default", st, v, None, T, E, B, apart, cmpT, origin_T or "none", None)
        self.extra.update(excluded=dict(self.excluded), moved_names=len(self.moved), locus_lines=self.locus_lines)

    def line(self, qid, kind, stratum, v, R, to, E, B, apart, cmp, name_src, sup, rule=None, locus=None):
        tr = self.tr
        Ev = tr.vset(E, to) if cmp == "versioned" else sorted(E)
        Bv = tr.vset(B, to) if cmp == "versioned" else sorted(B)
        rec = dict(qid=qid, family=FAMILY, kind=kind, stratum=stratum, unit=v, id=v, from_release=R, to_release=to,
                   final_database=None, strategy="both", rule=rule or rule_of(len(E)), cmp=cmp, expect=Ev,
                   broad=Bv if set(Bv) != set(Ev) else None, locus=locus, name_src=name_src)
        if apart:
            rec["apart"] = tr.vset(apart, to) if cmp == "versioned" else sorted(apart)
        if sup:
            rec["support"] = sup
        self.counts[stratum] += 1
        return rec


def b1(nb, dry, limit_releases=None):
    g = GenB1(nb, limit_releases)
    man = L.publish(FAMILY, nb.org, "b1_v1", g.lines(), "B1 every display name at every release (to R, to 115, default); "
                    "B2 moved names and B3 names of several genes as strata", __file__, judge="b_v1",
                    sources=["gtf/genes_<R>.tsv.gz", "biomart/<R>/", "core_xref/<R>.tsv.gz", "vgnc/", "vega/", "core/ (history)"],
                    manifest_extra=g.extra,
                    dry=dry, fmt="parquet", log=lambda *a: log(*a))
    log("published" if not dry else "dry", man.get("n"), man.get("counts"), "excluded", dict(g.excluded) if dry else "",
        "moved", len(g.moved), "locus", g.locus_lines)
    return man


def b8_lines(nb, counts=None):
    """Every gene with a display name at R -> its organism's name database at R (BioMart R, else the core xrefs)."""
    tr, org = nb.tr, nb.org
    db, grp, xdbs = NAME_DB.get(org, OTHER_DB)
    counts = collections.Counter() if counts is None else counts

    def gen():
        for R in nb.rels:
            G = nb.genes(R)
            named = {g for v, gs in nb.display(R).items() if TB.clean(v) for g in gs if g in G}
            if R in nb.bm_rels:
                pp = nb.C.bm(R, grp, 1)
                if not pp:
                    continue
                vals, src = pp[0], "biomart"
            else:
                src_all = nb.sources(R)
                if src_all[1] is None or "core_xref" not in src_all[1]:
                    continue
                by_gene = src_all[3]
                vals = {g: set().union(*(by_gene.get(g, {}).get(x, set()) for x in xdbs)) for g in named}
                if not any(vals.values()):
                    continue  # the release has no xref of this database (e.g. pig before HGNC's symbols): nothing to ask
                src = "core_xref"
            cmpR = tr.cmp_for(R)
            for g in sorted(named):
                E = sorted(vals.get(g, set()))
                counts[f"b8_{src}"] += 1
                by_case = collections.defaultdict(list)
                for x in E:
                    by_case[x.lower()].append(x)
                sp = [sorted(xs) for xs in by_case.values() if len(xs) > 1]
                extra = dict(spellings=sp) if sp else {}
                yield dict(**extra, **dict(qid=f"b8|{g}|{R}|{db}", family=FAMILY, kind="b8_name_db", stratum=f"b8_{src}", unit=g,
                           id=tr.vid(g, R) if cmpR == "versioned" else g, from_release=R, to_release=R, final_database=db,
                           strategy="both", rule="set", cmp="exact", expect=E, broad=None, locus=None, name_src=src))
            log(org, "B8 release", R, src)
    return gen()


def b8(nb, dry):
    db = NAME_DB.get(nb.org, OTHER_DB)[0]
    man = L.publish(FAMILY, nb.org, "b8_v1", b8_lines(nb), f"B8 every named gene to {db} at every release with a table",
                    __file__, judge="b_v1", sources=["biomart/<R>/", "core_xref/<R>.tsv.gz", "gtf/genes_<R>.tsv.gz"], dry=dry,
                    fmt="parquet", log=lambda *a: log(*a))
    log("published" if not dry else "dry", man.get("n"), man.get("counts"))
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="b1,b8")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--limit-releases", type=int, default=None, help="test: only the last N releases (with --dry)")
    a = ap.parse_args()
    if a.limit_releases and not a.dry:
        raise SystemExit("--limit-releases is for a dry test only")
    nb = TB.NamesB(a.organism, log=log)
    log(a.organism, "releases", nb.rels[0], "-", nb.rels[-1], "BioMart", sorted(nb.bm_rels), "core_xref", len(nb.xref_rels))
    for s in a.sets.split(","):
        if not a.dry and f"{s}_v1" in L.set_names(FAMILY, a.organism):
            log(a.organism, s, "published before: skipped (a rerun after an interruption)")
            continue
        {"b1": lambda: b1(nb, a.dry, a.limit_releases), "b8": lambda: b8(nb, a.dry)}[s]()


if __name__ == "__main__":
    main()
