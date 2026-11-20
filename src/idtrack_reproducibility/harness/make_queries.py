"""Build the fixed, seeded query sets of the accuracy harness from Ensembl's own files (never from idtrack).

Usage: make_queries.py --organism <name> --set <v1|v2|...> --kinds <kind,...|all> [--target 115]
       make_queries.py --organism <name> --set <new> --copies-apart | --names-first [--dry-run]   (a revision set)
       make_queries.py --organism <name> --revise <correction,...|all> [--dry-run]   (the corrected judge: one revision set
                       per correction, numbered after the last set; see REVISIONS)
       make_queries.py --organism <name> --set <new> --travel [--dry-run]   (the travel matrix: build_travel)

Writes `<work>/queries/<organism>/<set>.jsonl.gz` (one query per line) and `<set>.manifest.json` (counts per kind, sha256,
truth files used). A published set is never rewritten: the script refuses to overwrite an existing set; new kinds go
into a new set. Every kind draws from its own random generator (`random.Random(f"{SEED}:{organism}:{kind}")`) over a
sorted population, so adding kinds never changes the sample of another.

A query line: {qid, kind, id, from_release, to_release, final_database, strategy, rule, cmp, expect, broad, locus, note}
- rule `set`: right = `best` returns a non-empty subset of `expect` (or nothing when `expect` is empty); `all`
  returns exactly `expect`. `set_exact`: both strategies must return exactly `expect` (the right answer is the set).
- `broad` (optional) is a superset that Ensembl's own data also allow; an answer that is not right but within
  `broad` is `undetermined`.
- `locus` (only with empty `expect`): place of the query's gene; a returned gene that overlaps it on the same
  assembly and strand is `undetermined` (Ensembl's history has no successor, the genome says it may be the same
  gene), any other returned gene is wrong.
- cmp: `versioned` (returned Ensembl ids are compared with their version), `base` (without version), `exact`.
Fields of the corrected judge's revisions (README "The corrected judge"): `revised_by` (the corrections that changed the line),
`locus_mapped` (the locus mapped to the target's assembly, or why it cannot be), `support` (name sources per expected
gene, rule `support_ranks_shared_name`'s rule), `spellings` (spellings of one value), `excluded` (left out of the scores, with the reason).
A set generated now (`build`, `build_travel`) carries the corrections of the corrected judge (`apply_revisions`).
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import random
import re
import string
import sys

from idtrack_reproducibility.harness import harness_lib as H

SEED = 20261112
N = 3000          # sample size of a kind where the whole population is not taken
N_REL = 300       # per old release, for kinds stratified by release
GENE_DB = {       # idtrack database name -> (BioMart group, column) of the truth at that release
    "EntrezGene": ("entrez", 1),
    "NCBI gene (formerly Entrezgene)": ("entrez", 2),
    "HGNC Symbol": ("hgnc", 1),
    "MGI Symbol": ("mgi", 1),
    "UniProtKB Gene Name": ("uniprot_gn", 1),
    "Uniprot/SWISSPROT": ("swissprot", 1),
    "Uniprot/SPTREMBL": ("sptrembl", 1),
    "RefSeq_mRNA": ("refseq_mrna", 1),
    "RefSeq_peptide": ("refseq_peptide", 1),
}
NAME_GROUPS = [("hgnc", 1), ("mgi", 1), ("uniprot_gn", 1), ("synonym", 1), ("entrez", 2), ("wikigene", 1)]
NAME_SOURCES = [g for g in NAME_GROUPS if g[0] != "synonym"]  # the tables in which a string is a gene's name


def rng_for(org, kind):
    return random.Random(f"{SEED}:{org}:{kind}")


def sample(org, kind, population, k):
    pop = sorted(population)
    if len(pop) <= k:
        return pop
    return sorted(rng_for(org, kind).sample(pop, k))


class Truth:
    """Ensembl's history of one organism: the GTF genes of every release up to `target` and the mapping sessions and
    gene events of the release's core database, walked from one release to another.

    A mapping session is crossed where its effect appears (`placement="effective"`, the default; `effective_release`):
    between the last release before and the first release at which the GTF shows its new gene set. A session from
    release o to release n > o + 1 takes effect at the first release after o whose GTF holds the session's new ids
    (rhesus 50 -> 52 at 52, 52 -> 58 at 58: the GTFs of 50 and 51, and of 52 to 57, have the same gene ids); a
    mislabelled session (n not after o: pig 106 -> "104", human 114 -> "114") after o. `placement="old_release"` is
    the reading before 2026-11-18 (every session crossed between o and o + 1), kept to measure what the change moves
    (`rev_sessions_effective`).
    """

    def __init__(self, org, target, placement="effective"):
        self.org, self.T = org, target
        self.rels = [r for r in H.gtf_releases(org) if r <= target]
        self.genes = {r: H.gtf_genes(org, r) for r in self.rels}
        self.present = collections.defaultdict(set)
        for r, g in self.genes.items():
            for gid in g:
                self.present[gid].add(r)
        self.assembly = H.gtf_assemblies(org)
        self.core, ev, self.sessions = H.core_tables(org, target)
        self.succ = collections.defaultdict(lambda: collections.defaultdict(set))
        self.pred = collections.defaultdict(lambda: collections.defaultdict(set))
        for o, ov, n, nv, s, sc in ev:
            if o:
                self.succ[s][o].add(n)
            if n:
                self.pred[s][n].add(o)
        self.placement = placement
        self.eff = {s: self.effective_release(s) for s in self.sessions}
        self.order = sorted(self.sessions, key=lambda s: (self.eff[s], self.sessions[s]["created"]))
        self.ever = set(self.present) | {x for s in self.succ for x in self.succ[s]} | {x for s in self.pred for x in self.pred[s]}

    def effective_release(self, s):
        """The first release at which session s has taken effect (it is crossed between eff - 1 and eff).

        `old_release` placement: o + 1. `effective`: o + 1 for a session from o to o + 1 or a mislabelled one (new
        release not after the old); for o -> n > o + 1, the first GTF release r in (o, n] that holds at least half of
        the ids the session creates (new ids of its events absent from the GTF of o), else the first r in (o, n] whose
        GTF gene ids differ from those of o, else n (no GTF there: outside the window).
        """
        o, n = H.rel_int(self.sessions[s]["old_release"]), H.rel_int(self.sessions[s]["new_release"])
        if self.placement == "old_release" or n <= o + 1:
            return o + 1
        at_o = [r for r in self.rels if r <= o]
        between = [r for r in self.rels if o < r <= n]
        if not at_o or not between:
            return n
        base = self.genes[at_o[-1]]
        new = {x for x in self.pred[s] if x and x not in base}
        if new:
            for r in between:
                if 2 * len(new & set(self.genes[r])) >= len(new):
                    return r
        for r in between:
            if set(self.genes[r]) != set(base):
                return r
        return n

    def sessions_between(self, a, b):
        """Mapping sessions crossed going from release a to release b (a < b): a < effective release <= b."""
        return [s for s in self.order if a < self.eff[s] <= b]

    def sides(self, s):
        """The GTF releases on both sides of a session: the last one before it takes effect, the first one after."""
        if not hasattr(self, "_sides"):
            self._sides = {}
        if s not in self._sides:
            e = self.eff[s]
            before = [r for r in self.rels if r < e]
            after = [r for r in self.rels if r >= e]
            self._sides[s] = (before[-1] if before else None, after[0] if after else None)
        return self._sides[s]

    def lives(self, x, s):
        """Rule `id_lives_on` (the corrected judge): the id exists in the GTF on both sides of the session (it lives on, although
        Ensembl recorded only events of other ids into or out of it there)."""
        b0, a0 = self.sides(s)
        pr = self.present.get(x, ())
        return b0 is not None and a0 is not None and b0 in pr and a0 in pr

    def comes_back(self, x, s, forward):
        """Rule `resurrection` (the corrected judge): Ensembl retires the id in this session (forward; or creates it, backward) and
        the GTF has the same stable id again after (forward; or before, backward) the session."""
        e = self.eff[s]
        pr = self.present.get(x, ())
        return any(r >= e for r in pr) if forward else any(r < e for r in pr)

    def _step(self, x, s, strict, living, resurrect, forward):
        moves = (self.succ if forward else self.pred)[s]
        if x not in moves:
            return {x}
        others = moves[x] - {None}
        keep = x in others or (living and self.lives(x, s))
        if strict:
            if keep:
                return {x}
            if not others and resurrect and self.comes_back(x, s, forward):
                return {x}
            return others
        out = others | ({x} if keep else set())
        if not others and resurrect and self.comes_back(x, s, forward):
            out.add(x)
        return out

    def walk_fwd(self, gid, a, b, strict, living=False, resurrect=False):
        """From gid at release a to release b > a through Ensembl's mapping sessions.

        `strict`: where a session maps the id to itself only itself is followed; otherwise all its successors.
        `living` (the corrected judge after rule `id_lives_on`): an id present in the GTF on both sides of a session is kept, as if
        the session mapped it to itself. `resurrect` (the corrected judge after rule `resurrection`): an id retired in a session and
        present again later under the same stable id is kept. The published sets use neither.
        """
        fr = {gid}
        for s in self.sessions_between(a, b):
            nxt = set()
            for x in fr:
                nxt |= self._step(x, s, strict, living, resurrect, True)
            fr = nxt
        return {x for x in fr if b in self.present.get(x, ())}

    def walk_back(self, gid, a, b, strict, living=False, resurrect=False):
        """From gid at release a back to release b < a (the mirror image of `walk_fwd`)."""
        fr = {gid}
        for s in reversed(self.sessions_between(b, a)):
            nxt = set()
            for x in fr:
                nxt |= self._step(x, s, strict, living, resurrect, False)
            fr = nxt
        return {x for x in fr if b in self.present.get(x, ())}

    def walk_series(self, gid, a, targets, strict, living=False, resurrect=False):
        """`walk_fwd` (targets after a) or `walk_back` (targets before a) from gid at release a to every target in one
        pass over the sessions: {target: ids present there}; a target equal to a gives {gid}."""
        out = {}
        if a in targets:
            out[a] = {gid}
        for forward in (True, False):
            tt = sorted((t for t in targets if (t > a if forward else t < a)), reverse=not forward)
            if not tt:
                continue
            ss = self.sessions_between(a, tt[-1]) if forward else list(reversed(self.sessions_between(tt[-1], a)))
            fr, i = {gid}, 0
            for t in tt:
                if forward:
                    while i < len(ss) and self.eff[ss[i]] <= t:
                        fr = set().union(*(self._step(x, ss[i], strict, living, resurrect, True) for x in fr)) if fr else set()
                        i += 1
                else:
                    while i < len(ss) and self.eff[ss[i]] > t:
                        fr = set().union(*(self._step(x, ss[i], strict, living, resurrect, False) for x in fr)) if fr else set()
                        i += 1
                out[t] = {x for x in fr if t in self.present.get(x, ())}
        return out

    def ver(self, gid, rel):
        g = self.genes.get(rel, {}).get(gid)
        return g["gene_version"] if g else ""

    def vid(self, gid, rel):
        v = self.ver(gid, rel)
        return f"{gid}.{v}" if v else gid

    def vset(self, ids, rel):
        return sorted(self.vid(g, rel) for g in ids)

    def cmp_for(self, rel):
        g = self.genes.get(rel) or {}
        return "versioned" if g and next(iter(g.values()))["gene_version"] else "base"

    def locus(self, gid, rel):
        g = self.genes.get(rel, {}).get(gid)
        if not g:
            return None
        return dict(release=rel, assembly=self.assembly.get(rel, ""), seqname=g["seqname"], start=int(g["start"]),
                    end=int(g["end"]), strand=g["strand"])

    def last_before(self, gid, b):
        rr = [r for r in self.present.get(gid, ()) if r < b]
        return max(rr) if rr else None

    def eras(self):
        """[(assembly, first release, last release)] over the GTF releases up to the target."""
        out = []
        for r in self.rels:
            a = self.assembly.get(r, "")
            if out and out[-1][0] == a:
                out[-1][2] = r
            else:
                out.append([a, r, r])
        return [tuple(x) for x in out]


def q(kind, ident, to_release, expect, strategy="best", from_release=None, final_database=None, rule="set", cmp="versioned",
      broad=None, locus=None, note=None):
    rec = dict(kind=kind, id=ident, from_release=from_release, to_release=to_release, final_database=final_database,
               strategy=strategy, rule=rule, cmp=cmp, expect=sorted(expect), broad=(sorted(broad) if broad is not None else None),
               locus=locus, note=note)
    rec["qid"] = "|".join(str(x) for x in (kind, ident, from_release, to_release, final_database, strategy))
    return rec


def name_tables(org, rel):
    """display name -> genes, and any name (all name sources) -> genes, at one BioMart release."""
    display = collections.defaultdict(set)
    for r in H.biomart_genes(org, rel) or []:
        if r["name"]:
            display[r["name"]].add(r["gene_id"])
    anyname = collections.defaultdict(set)
    for v, gs in display.items():
        anyname[v] |= gs
    for grp, col in NAME_GROUPS:
        pp = H.biomart_pairs(org, rel, grp, col)
        if pp:
            for v, gs in pp[1].items():
                anyname[v] |= gs
    return display, anyname


def exact_names(org, rel):
    """name, exactly as written -> the genes that carry it as a name (not as a synonym) at one BioMart release.

    A name is the display name of a gene or its value in one of the name tables (``NAME_SOURCES``: HGNC and MGI
    symbol, UniProt gene name, NCBI gene label, WikiGene name). Rule `name_before_synonym`: a string that is one gene's name and
    another gene's synonym means the name's gene.
    """
    names = collections.defaultdict(set)
    for r in H.biomart_genes(org, rel) or []:
        if r["name"]:
            names[r["name"]].add(r["gene_id"])
    for grp, col in NAME_SOURCES:
        pp = H.biomart_pairs(org, rel, grp, col)
        if pp:
            for v, gs in pp[1].items():
                names[v] |= gs
    return names


def name_before_synonym(synonym_genes, name_genes):
    """Rule `name_before_synonym`: the genes a synonym query expects, given the genes that carry the string as an exact name.

    Where the genes that carry the string as an exact name differ from the synonym's genes, the string means the
    name and those genes alone are expected: the name of a gene that is not among the synonym's genes, and also the
    name of only some of the synonym's genes (integration Q16), since names are tried before synonyms at each step
    (rule `names_before_synonyms`). Where no gene carries it as a name, or exactly the synonym's genes do, the synonym's genes are
    expected, as before. Both sets are genes of the primary sequence at the release asked about.
    """
    name_genes, synonym_genes = set(name_genes), set(synonym_genes)
    return name_genes if name_genes and name_genes != synonym_genes else synonym_genes


def build(org, kinds, T):
    tr = Truth(org, T)
    cur = tr.genes[T]
    primary = {g for g, r in cur.items() if r["region"] == "primary"}
    patch = set(cur) - primary
    cmpT = tr.cmp_for(T)
    curset = set(cur)
    out = collections.defaultdict(list)
    want = lambda k: kinds == {"all"} or k in kinds  # noqa: E731
    display, anyname = name_tables(org, T)
    any_lower = collections.defaultdict(set)
    for v, gs in anyname.items():
        any_lower[v.lower()] |= gs
    ens_like = re.compile(r"^ENS[A-Z]*[GTP]\d{11}")
    entrez = H.biomart_pairs(org, T, "entrez", 1)
    entrez_v2g = entrez[1] if entrez else {}
    numeric_names = {v for v in anyname if v.isdigit()}

    def name_rule(E):
        return "set_exact" if len(E) > 1 else "set"

    # ---- current release: Ensembl ids ----
    for sub, pop in (("", primary), ("_altloci", patch)):
        if not pop:
            continue
        k = f"gene_versioned_to_self{sub}"
        if want(k):
            for g in sorted(pop):
                out[k].append(q(k, tr.vid(g, T), T, [tr.vid(g, T)], cmp=cmpT))
            for g in sample(org, k + ":all", pop, N):
                out[k].append(q(k, tr.vid(g, T), T, [tr.vid(g, T)], strategy="all", cmp=cmpT))
        k = f"gene_unversioned_to_self{sub}"
        if want(k):
            for g in sorted(pop):
                out[k].append(q(k, g, T, [tr.vid(g, T)], cmp=cmpT))
            for g in sample(org, k + ":from", pop, N):
                out[k].append(q(k, g, T, [tr.vid(g, T)], from_release=T, cmp=cmpT))
            for g in sample(org, k + ":all", pop, N):
                out[k].append(q(k, g, T, [tr.vid(g, T)], strategy="all", cmp=cmpT))
    txT = H.gtf_tx(org, T) or []
    for k, idcol, vcol in (("transcript_to_gene", "transcript_id", "transcript_version"),
                           ("translation_to_gene", "protein_id", "protein_version")):
        if not want(k):
            continue
        rows = {r[idcol]: r for r in txT if r[idcol] and r["region"] == "primary"}
        for t in sample(org, k, rows, N):
            r = rows[t]
            E = [tr.vid(r["gene_id"], T)]
            out[k].append(q(k, t, T, E, cmp=cmpT, note="unversioned"))
            if r[vcol]:
                out[k].append(q(k, f"{t}.{r[vcol]}", T, E, cmp=cmpT, note="versioned"))

    # ---- current release: gene -> external database ----
    k = "gene_to_db"
    if want(k):
        genes_bm = {r["gene_id"] for r in H.biomart_genes(org, T)} & primary
        gs = sample(org, k, genes_bm, 1500)
        for db, (grp, col) in GENE_DB.items():
            if db in H.NOT_A_DATABASE_OF.get(org, ()):
                continue
            pp = H.biomart_pairs(org, T, grp, col)
            if pp is None or not pp[1]:
                continue
            for g in gs:
                E = pp[0].get(g, set())
                for st in ("best", "all"):
                    out[k].append(q(k, tr.vid(g, T), T, E, strategy=st, from_release=T, final_database=db, cmp="exact",
                                    note=db))

    # ---- current release: names and external ids -> gene ----
    def name_queries(k, names, strict_of, broad_of, n=None, both=True, all_n=2000):
        pop = sorted(names)
        chosen = pop if n is None else sample(org, k, pop, n)
        for v in chosen:
            E = strict_of(v)
            B = broad_of(v) | E
            Ev, Bv = tr.vset(E, T), tr.vset(B, T)
            out[k].append(q(k, v, T, Ev, rule=name_rule(E), cmp=cmpT, broad=Bv if B != E else None))
            if both:
                out[k].append(q(k, v, T, Ev, from_release=T, rule=name_rule(E), cmp=cmpT, broad=Bv if B != E else None))
        for v in sample(org, k + ":all", chosen, all_n):
            E = strict_of(v)
            B = broad_of(v) | E
            out[k].append(q(k, v, T, tr.vset(E, T), strategy="all", from_release=T, rule=name_rule(E), cmp=cmpT,
                            broad=tr.vset(B, T) if B != E else None))

    clean = lambda v: not ens_like.match(v) and not v.isdigit() and v.strip() == v and v != ""  # noqa: E731
    disp_ok = {v: gs & primary for v, gs in display.items() if clean(v) and gs & primary}
    if want("name_unique_to_gene"):
        name_queries("name_unique_to_gene", [v for v, gs in disp_ok.items() if len(gs) == 1],
                     lambda v: disp_ok[v], lambda v: anyname[v] & curset)
    if want("name_ambiguous_to_set"):
        name_queries("name_ambiguous_to_set", [v for v, gs in disp_ok.items() if len(gs) > 1],
                     lambda v: disp_ok[v], lambda v: anyname[v] & curset)
    syn = H.biomart_pairs(org, T, "synonym", 1)
    if want("synonym_to_gene") and syn:
        pop = [v for v in syn[1] if clean(v) and v not in display and (syn[1][v] & primary)]
        named = exact_names(org, T)
        name_queries("synonym_to_gene", pop,
                     lambda v: name_before_synonym(syn[1][v] & primary, named.get(v, set()) & primary),
                     lambda v: anyname[v] & curset, n=N)
    if want("ncbi_to_gene") and entrez_v2g:
        pop = [v for v in entrez_v2g if v.isdigit() and v not in anyname and entrez_v2g[v] & primary]
        name_queries("ncbi_to_gene", pop, lambda v: entrez_v2g[v] & primary, lambda v: set())
    for k, groups in (("refseq_to_gene", ["refseq_mrna", "refseq_ncrna", "refseq_peptide", "refseq_mrna_predicted",
                                          "refseq_ncrna_predicted", "refseq_peptide_predicted"]),
                      ("uniprot_to_gene", ["swissprot", "sptrembl"])):
        if not want(k):
            continue
        v2g = collections.defaultdict(set)
        for grp in groups:
            pp = H.biomart_pairs(org, T, grp, 1)
            if pp:
                for v, gs in pp[1].items():
                    v2g[v] |= gs
        pop = [v for v in v2g if v2g[v] & primary and v not in anyname]
        name_queries(k, pop, lambda v: v2g[v] & primary, lambda v: set(), n=N)

    # ---- spelling variants ----
    k = "case_variant"
    if want(k):
        uniq = [v for v, gs in disp_ok.items() if len(gs) == 1]
        cands = []
        for v in uniq:
            for w, how in ((v.lower(), "lower"), (v.upper(), "upper")):
                if w != v and w not in anyname and any_lower[w.lower()] == anyname[v]:
                    cands.append((w, v, how))
        for w, v, how in sample(org, k, cands, N):
            out[k].append(q(k, w, T, tr.vset(disp_ok[v], T), cmp=cmpT, note=f"{how} of {v}"))
        for g in sample(org, k + ":ensembl", primary, 1000):
            out[k].append(q(k, g.lower(), T, [tr.vid(g, T)], cmp=cmpT, note=f"lower of {g}"))
    k = "separator_variant"
    if want(k):
        cands = []
        for v, gs in disp_ok.items():
            if len(gs) != 1 or not re.search(r"[-_]", v):
                continue
            w = v.translate(str.maketrans("-_", "_-"))
            if w.lower() not in any_lower:
                cands.append((w, v))
        for w, v in sample(org, k, cands, N):
            out[k].append(q(k, w, T, tr.vset(disp_ok[v], T), cmp=cmpT, note=f"separators of {v}"))

    # ---- ids that never existed ----
    k = "nonexistent"
    if want(k):
        rg = rng_for(org, k)
        prefix = H.ORGANISMS[org][1]
        allnames_lower = set(any_lower)
        for rel in H.biomart_releases(org):
            if rel != T:
                _, an = name_tables(org, rel)
                allnames_lower |= {v.lower() for v in an}
        made = set()
        while len(made) < 500:
            made.add((f"{prefix}{rg.randrange(10**10, 10**11):011d}", "unknown id with the organism's prefix"))
        while len(made) < 1000:
            s = "".join(rg.choice(string.ascii_uppercase) for _ in range(rg.randrange(7, 11)))
            if s.lower() not in allnames_lower:
                made.add((s, "random letters"))
        while len(made) < 1500:
            s = f"NM_{rg.randrange(10**8, 10**9)}"
            if s.lower() not in allnames_lower:
                made.add((s, "refseq-like accession that does not exist"))
        for ident, note in sorted(made):
            if ident in tr.ever:
                continue
            out[k].append(q(k, ident, T, [], cmp=cmpT, note=note))
    k = "nonexistent_version"
    if want(k):
        for g in sample(org, k, primary, 1000):
            vs = [int(tr.ver(g, r)) for r in tr.present[g] if tr.ver(g, r).isdigit()]
            if vs:
                out[k].append(q(k, f"{g}.{max(vs) + 37}", T, [], cmp=cmpT, note="version that never existed"))
    # Rule `version_newer_than_graph` (docs/rules.md) splits `nonexistent_version` into the two kinds below (published in
    # the set after the one that holds `nonexistent_version`, which stays as it is and is no longer part of the
    # acceptance: its queries are `version_newer_than_graph` queries with the opposite expectation).
    k = "version_newer_than_graph"
    if want(k):
        # a current id with a version above every version it ever had: the version may belong to a release newer
        # than the graph, so the answer is the id's gene (flagged by idtrack)
        for i, g in enumerate(sample(org, k, primary, 1000)):
            vs = [int(tr.ver(g, r)) for r in tr.present[g] if tr.ver(g, r).isdigit()]
            if vs:
                step = 1 if i % 2 == 0 else 37
                out[k].append(q(k, f"{g}.{max(vs) + step}", T, [tr.vid(g, T)], cmp=cmpT,
                                note=f"current id, version {step} above its highest: may be newer than the graph"))
    k = "version_unknown"
    if want(k):
        # versions that no release ever had and that cannot be newer than the graph: the answer is nothing
        retired, gap, lower = {}, {}, {}
        for g, rr in tr.present.items():
            given = [tr.ver(g, r) for r in sorted(rr)]
            vs = sorted({int(v) for v in given if v.isdigit()})
            if not vs:
                continue
            if T not in rr:
                retired[g] = (max(vs) + 1, "retired id, one above its last version")
            elif g in primary and all(v.isdigit() for v in given):
                # every release of the id has a version in the GTF, so the versions it had are known completely
                holes = sorted(set(range(vs[0], vs[-1])) - set(vs))
                if holes:
                    gap[g] = (holes[0], "current id, a version between two of its versions that no release had")
                if vs[0] > 1 and min(rr) > tr.rels[0]:
                    lower[g] = (vs[0] - 1, "current id born after the first release, one below its first version")
        for name, pop, n in (("retired", retired, 1000), ("gap", gap, 500), ("lower", lower, 500)):
            for g in sample(org, f"{k}:{name}", pop, n):
                version, note = pop[g]
                out[k].append(q(k, f"{g}.{version}", T, [], cmp=cmpT, note=note))

    # ---- history: old ids to the current release ----
    eras = tr.eras()
    old_rels = [r for r in tr.rels if r < T]
    hist_cls = collections.defaultdict(list)  # class -> [(gid, L)]
    if any(want(x) for x in ("retired_renamed", "retired_merged", "retired_split", "retired_none", "retired_none_branch")):
        for gid, rr in tr.present.items():
            if T in rr:
                continue
            L = max(rr)
            E = tr.walk_fwd(gid, L, T, True)
            B = tr.walk_fwd(gid, L, T, False)
            if not E:
                cls = "retired_none" if not B else "retired_none_branch"
            elif len(E) > 1:
                cls = "retired_split"
            else:
                merged = False
                for s in tr.sessions_between(L, T):
                    news = tr.succ[s].get(gid, set()) - {None}
                    if news and gid not in news:
                        merged = any(len(tr.pred[s][n] - {None}) > 1 for n in news)
                        break
                cls = "retired_merged" if merged else "retired_renamed"
            hist_cls[cls].append((gid, L, E, B))
        for cls, items in hist_cls.items():
            if not want(cls):
                continue
            per = -(-N // len(eras))
            chosen = []
            for a, lo, hi in eras:
                pool = [x for x in items if lo <= x[1] <= hi]
                picked = set(sample(org, f"{cls}:{a}", [y[0] for y in pool], per))
                chosen += [x for x in pool if x[0] in picked]
            for gid, L, E, B in sorted(chosen):
                Ev, Bv = tr.vset(E, T), tr.vset(B, T)
                loc = tr.locus(gid, L) if not E else None
                common = dict(cmp=cmpT, broad=Bv if B != E else None, locus=loc, note=f"last release {L}")
                out[cls].append(q(cls, gid, T, Ev, **common))
                out[cls].append(q(cls, gid, T, Ev, from_release=L, **common))
                out[cls].append(q(cls, gid, T, Ev, from_release=L, strategy="all", **common))

    # old releases used for the per-release kinds: both sides of every assembly handoff, the first, and evenly between
    picks = {eras[0][1]}
    for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
        picks |= {hi1, lo2}
    picks |= {old_rels[i] for i in range(0, len(old_rels), max(1, len(old_rels) // 6))}
    picks = sorted(r for r in picks if r < T)
    k = "unchanged_old_to_current"
    if want(k):
        for R in picks:
            pop = [g for g in tr.genes[R] if g in cur and tr.walk_fwd(g, R, T, True) == {g}]
            for g in sample(org, f"{k}:{R}", pop, N_REL):
                out[k].append(q(k, g, T, [tr.vid(g, T)], from_release=R, cmp=cmpT, note=f"from {R}"))
                if tr.ver(g, R):
                    out[k].append(q(k, tr.vid(g, R), T, [tr.vid(g, T)], from_release=R, cmp=cmpT, note=f"from {R}, versioned"))
    if want("reverse_present") or want("reverse_absent"):
        for R in picks:
            cR = tr.cmp_for(R)
            pres = [g for g in primary if R in tr.present[g]]
            absent = [g for g in primary if R not in tr.present[g]]
            for k, pop in (("reverse_present", pres), ("reverse_absent", absent)):
                if not want(k):
                    continue
                for g in sample(org, f"{k}:{R}", pop, N_REL):
                    E, B = tr.walk_back(g, T, R, True), tr.walk_back(g, T, R, False)
                    Ev, Bv = tr.vset(E, R), tr.vset(B, R)
                    loc = tr.locus(g, T) if not E else None
                    for st in ("best", "all"):
                        out[k].append(q(k, g, R, Ev, strategy=st, from_release=T, cmp=cR, broad=Bv if B != E else None,
                                        locus=loc, note=f"to {R}"))
    k = "cross_assembly"
    if want(k):
        for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
            for R1, R2 in sorted({(hi1, lo2), (max(lo1, hi1 - 5), min(hi2, lo2 + 5))}):
                if R2 >= T:
                    continue
                cR = tr.cmp_for(R2)
                for g in sample(org, f"{k}:{R1}:{R2}", tr.genes[R1], N_REL):
                    E, B = tr.walk_fwd(g, R1, R2, True), tr.walk_fwd(g, R1, R2, False)
                    loc = tr.locus(g, R1) if not E else None
                    for st in ("best", "all"):
                        out[k].append(q(k, g, R2, tr.vset(E, R2), strategy=st, from_release=R1, cmp=cR,
                                        broad=tr.vset(B, R2) if B != E else None, locus=loc, note=f"{a1} {R1} -> {a2} {R2}"))

    # ---- old releases: names and external databases (BioMart archives) ----
    bm_old = [r for r in H.biomart_releases(org) if r < T and r in tr.genes]
    for R in bm_old:
        dR, aR = name_tables(org, R)
        bm_genes = {r["gene_id"] for r in H.biomart_genes(org, R)} & set(tr.genes[R])
        k = "old_gene_to_db"
        if want(k):
            gs = sample(org, f"{k}:{R}", bm_genes, N_REL)
            for db, (grp, col) in GENE_DB.items():
                if db in H.NOT_A_DATABASE_OF.get(org, ()):
                    continue
                pp = H.biomart_pairs(org, R, grp, col)
                if pp is None or not pp[1]:
                    continue
                for g in gs:
                    out[k].append(q(k, tr.vid(g, R), R, pp[0].get(g, set()), from_release=R, final_database=db, cmp="exact",
                                    note=f"{db} at {R}"))
        k = "name_old_to_current"
        if want(k):
            pop = [v for v, gs in dR.items() if clean(v) and gs & bm_genes]
            for v in sample(org, f"{k}:{R}", pop, 500):
                src = dR[v] & bm_genes
                E = set().union(*(tr.walk_fwd(g, R, T, True) for g in src))
                B = set().union(*(tr.walk_fwd(g, R, T, False) for g in src)) | E
                loc = tr.locus(next(iter(src)), R) if not E and len(src) == 1 else None
                if not E:  # Ensembl's history ends: a gene that carries the name now is not refuted by it
                    B |= anyname.get(v, set()) & curset
                rule = "set_exact" if len(E) > 1 else "set"
                for st in ("best", "all"):
                    out[k].append(q(k, v, T, tr.vset(E, T), strategy=st, from_release=R, rule=rule, cmp=cmpT,
                                    broad=tr.vset(B, T) if B != E else None, locus=loc, note=f"name at {R}"))
        k = "name_changed_gene"
        if want(k):
            pop = []
            for v, gs in dR.items():
                if not clean(v) or len(gs) != 1 or v not in disp_ok or len(disp_ok[v]) != 1:
                    continue
                (a,), (b,) = gs, disp_ok[v]
                if a != b and a in tr.genes[R] and b not in tr.walk_fwd(a, R, T, False):
                    pop.append(v)
            for v in sample(org, f"{k}:{R}", pop, 500):
                (a,), (b,) = dR[v], disp_ok[v]
                E = tr.walk_fwd(a, R, T, True)
                B = tr.walk_fwd(a, R, T, False) | E
                if not E:  # Ensembl's history of the old gene ends: the gene with the name now is not refuted
                    B |= {b}
                out[k].append(q(k, v, T, [tr.vid(b, T)], cmp=cmpT, note=f"default: the gene with the name at {T}"))
                out[k].append(q(k, v, T, tr.vset(E, T), from_release=R, cmp=cmpT, broad=tr.vset(B, T) if B != E else None,
                                locus=tr.locus(a, R) if not E else None, note=f"from {R}: the fate of {a}"))
    return out, tr


def copies_apart(org, lines):
    """Rule `judge_copies_apart`: the truth of gene answers follows rule `copies_apart` where the GTF of the target release marks patches.

    Args:
        org: The organism.
        lines: The published query lines (the latest line of every qid).

    Returns:
        tuple: ``{kind: [revised lines]}`` (only the lines that change, with their published ``qid``) and the
        releases whose GTF marks genes on patches and alternative loci (the only target releases revised).

    A line changes when it expects Ensembl gene ids (no ``final_database``) and its expected set holds at least one
    gene of the primary sequence and at least one gene on a patch or alternative locus at the target release
    (``region`` of the GTF tables: the primary GTF file against the ``chr_patch_hapl_scaff`` file of the release).
    Its ``expect`` becomes the genes of the primary sequence, the others go to the new field ``apart`` (the copies
    that the answer must report apart, in ``alternative_locus_genes``); ``rule`` ``set_exact`` stays only while
    more than one gene is expected; everything else (``broad``, ``locus``, ``cmp``, the query itself) is kept.
    Lines whose expected genes are all on one side, lines with a ``final_database``, and lines whose target
    release has no patch file in the truth (the releases before the first one listed) are not touched.
    """
    regions, marked = {}, []
    out = collections.defaultdict(list)
    for rec in lines:
        if rec["final_database"] is not None or len(rec["expect"]) < 2:
            continue
        release = rec["to_release"]
        if release not in regions:
            try:
                regions[release] = {g: row["region"] for g, row in H.gtf_genes(org, release).items()}
            except FileNotFoundError:
                regions[release] = {}
            if "patch" in regions[release].values():
                marked.append(release)
        where = regions[release]
        if release not in marked:
            continue
        sides = collections.defaultdict(list)
        for gene in rec["expect"]:
            sides[where.get(gene.split(".")[0] if rec["cmp"] == "versioned" else gene)].append(gene)
        if None in sides or not sides["primary"] or not sides["patch"]:
            continue
        new = dict(rec)
        new.pop("set", None)
        new["expect"], new["apart"] = sorted(sides["primary"]), sorted(sides["patch"])
        if rec["rule"] == "set_exact" and len(new["expect"]) == 1:
            new["rule"] = "set"
        new["note"] = ((rec.get("note") or "") + "; " if rec.get("note") else "") + "copies apart (rule `judge_copies_apart`)"
        out[rec["kind"]].append(new)
    return out, sorted(marked)


def names_first(org, lines):
    """Rule `name_before_synonym`: a published ``synonym_to_gene`` line whose string is an exact name of genes expects those genes.

    Args:
        org: The organism.
        lines: The published query lines (the latest line of every qid).

    Returns:
        tuple: ``{kind: [revised lines]}`` (only the lines that change, with their published ``qid``) and the
        releases whose name tables were read.

    The lines of the kind were drawn from BioMart's synonym table without asking whether the string is also a
    name. A line changes when, at the release it asks about (``from_release``, else ``to_release``), the genes of
    the primary sequence that carry the string as an exact name (``exact_names``) are not empty and differ from
    the genes the line expects: a name of a gene that is not among them, or a name of only some of them
    (integration Q16). Its ``expect`` becomes the genes that carry the string as a name (``name_before_synonym``,
    the rule ``build`` applies when a set is generated). ``rule`` is ``set_exact`` when several genes are
    expected and ``set`` when one; ``broad`` (the genes that carry the string in any name source, the synonym
    table among them) and everything else are kept, so an answer with the synonym's genes is ``undetermined``, as
    an answer with the name's genes was before; where the line had no ``broad`` (its expected genes were all the
    genes that carry the string in any source), the genes it expected become ``broad``, as ``build`` gives it. A
    line whose string is the name of exactly the genes it expects, or of none, is not touched.
    """
    tables = {}
    out = collections.defaultdict(list)
    for rec in lines:
        if rec["kind"] != "synonym_to_gene" or rec["final_database"] is not None:
            continue
        release = rec["from_release"] or rec["to_release"]
        if release not in tables:
            genes = H.gtf_genes(org, release)
            tables[release] = (exact_names(org, release), genes, {g for g, row in genes.items() if row["region"] == "primary"})
        names, genes, primary = tables[release]
        base = (lambda g: g.split(".")[0]) if rec["cmp"] == "versioned" else (lambda g: g)
        expected = {base(g) for g in rec["expect"]}
        revised = name_before_synonym(expected, names.get(rec["id"], set()) & primary)
        if revised == expected:
            continue
        new = dict(rec)
        new.pop("set", None)
        version = (lambda g: f"{g}.{genes[g]['gene_version']}") if rec["cmp"] == "versioned" else (lambda g: g)
        new["expect"] = sorted(version(g) for g in revised)
        new["rule"] = "set_exact" if len(revised) > 1 else "set"
        if new["broad"] is None:  # the expected genes were all the genes with the string in any source
            new["broad"] = sorted(rec["expect"])
        if not set(new["expect"]) <= set(new["broad"]):
            raise SystemExit(f"{rec['qid']}: the name's genes are not within the broad set")
        if set(new["expect"]) == set(new["broad"]):
            new["broad"] = None
        new["note"] = ((rec.get("note") or "") + "; " if rec.get("note") else "") + "name before synonym (rule `name_before_synonym`)"
        out[rec["kind"]].append(new)
    return out, sorted(tables)


# ---------------- the corrected judge: the judge corrected and widened (revision sets; README "The corrected judge") ----------------

ENS_ID = re.compile(rf"^({H.STABLE_ID_FORM})(\.\d+)?$")


def gene_base(x):
    m = ENS_ID.match(str(x))
    return m.group(1) if m else str(x)


class Ctx78:
    """What the revisions of the corrected judge read: the history (``Truth``), BioMart's tables, VGNC, Vega, the GTF."""

    def __init__(self, org, T, truth=None):
        self.org, self.T, self._tr, self._c = org, T, truth, {}

    @property
    def tr(self):
        if self._tr is None:
            self._tr = Truth(self.org, self.T)
        return self._tr

    def _cached(self, key, make):
        if key not in self._c:
            self._c[key] = make()
        return self._c[key]

    def bm(self, rel, group, col=1):
        return self._cached(("bm", rel, group, col), lambda: H.biomart_pairs(self.org, rel, group, col))

    def display(self, rel):
        """(name -> {genes with it as display name}, gene -> (name, source)) at a BioMart release."""
        def make():
            n2g, g2n = collections.defaultdict(set), {}
            for r in H.biomart_genes(self.org, rel) or []:
                g2n[r["gene_id"]] = (r["name"], r["source"])
                if r["name"]:
                    n2g[r["name"]].add(r["gene_id"])
            return n2g, g2n
        return self._cached(("display", rel), make)

    def gtf(self, rel):
        return self._cached(("gtf", rel), lambda: H.gtf_genes(self.org, rel) or {})

    def primary(self, rel):
        return self._cached(("primary", rel), lambda: {g for g, r in self.gtf(rel).items() if r["region"] == "primary"})

    def vgnc(self, rel):
        return self._cached(("vgnc", rel), lambda: H.vgnc_symbols(self.org, rel))

    def vega(self, rel):
        return self._cached(("vega", rel), lambda: H.vega_names(self.org, rel))

    def name_genes(self, rel, v, synonyms=False):
        """Genes that carry string v at a release as a name: display name, the BioMart name tables (HGNC, MGI,
        UniProt gene name, NCBI label, WikiGene), the VGNC symbol (VGNC's symbol of BioMart's VGNC id) and the Vega
        gene name (core xrefs); with `synonyms` also BioMart's synonyms."""
        out = set(self.display(rel)[0].get(v, ()))
        for grp, col in NAME_SOURCES + ([("synonym", 1)] if synonyms else []):
            pp = self.bm(rel, grp, col)
            if pp:
                out |= pp[1].get(v, set())
        for pp in (self.vgnc(rel), self.vega(rel)):
            if pp:
                out |= pp[1].get(v, set())
        return out

    def support(self, rel, gene, v):
        """In how many of BioMart's name sources a gene carries string v at a release (rule `support_ranks_shared_name`): the display
        name, HGNC, MGI, UniProt gene name, NCBI label, WikiGene, synonym."""
        n = int(self.display(rel)[1].get(gene, ("", ""))[0] == v)
        for grp, col in NAME_SOURCES + [("synonym", 1)]:
            pp = self.bm(rel, grp, col)
            if pp and v in pp[0].get(gene, ()):
                n += 1
        return n

    def bm_genes(self, rel):
        return self._cached(("bm_genes", rel), lambda: {r["gene_id"] for r in H.biomart_genes(self.org, rel) or []}
                            & set(self.tr.genes.get(rel, {})))


def revised(rec, name, note, **changes):
    new = dict(rec)
    new.pop("set", None)
    new.update(changes)
    new["revised_by"] = list(rec.get("revised_by") or []) + [name]
    new["note"] = ((rec.get("note") or "") + "; " if rec.get("note") else "") + note
    if new.get("broad") is not None and set(new["broad"]) == set(new["expect"]):
        new["broad"] = None
    return new


def _vers(C, genes, rel, cmp):
    return sorted(C.tr.vid(g, rel) for g in genes) if cmp == "versioned" else sorted(genes)


def _name_of_variant(rec):
    m = re.match(r"^(lower|upper|separators) of (.*)$", rec.get("note") or "")
    return m.group(2) if m else None


def _history_walks(C, rec, living=True, resurrect=False):
    """(strict, broad, release of the answer) of a history line from Ensembl's sessions, or None when the line's
    truth does not come from a walk of the query id (or of the genes that carried a name at R)."""
    tr, k, T = C.tr, rec["kind"], C.T
    if k in ("name_old_to_current", "name_changed_gene"):
        R = rec["from_release"]
        if R is None:
            return None
        n2g = C.display(R)[0]
        src = n2g.get(rec["id"], set()) & C.bm_genes(R) if k == "name_old_to_current" else n2g.get(rec["id"], set())
        if not src:
            return None
        E = set().union(*(tr.walk_fwd(g, R, T, True, living, resurrect) for g in src))
        B = set().union(*(tr.walk_fwd(g, R, T, False, living, resurrect) for g in src)) | E
        return E, B, T
    gid = gene_base(rec["id"])
    if gid not in tr.present:
        return None
    if k in ("reverse_present", "reverse_absent"):
        a, b = T, rec["to_release"]
        return tr.walk_back(gid, a, b, True, living, resurrect), tr.walk_back(gid, a, b, False, living, resurrect), b
    if k == "cross_assembly":
        a, b = rec["from_release"], rec["to_release"]
        return tr.walk_fwd(gid, a, b, True, living, resurrect), tr.walk_fwd(gid, a, b, False, living, resurrect), b
    if k.startswith("retired_"):
        L = max(tr.present[gid])
        return tr.walk_fwd(gid, L, T, True, living, resurrect), tr.walk_fwd(gid, L, T, False, living, resurrect), T
    return None


def rev_swissprot75(org, lines, C):
    """BioMart 75's `uniprot_swissprot` is the entry name (`VATG2_PIG`); the accession is its attribute
    `uniprot_swissprot_accession` (fetched as the table `swissprot_acc`). Lines to Uniprot/SWISSPROT at a release
    that has that table expect the accessions."""
    out = collections.defaultdict(list)
    for rec in lines:
        if rec["final_database"] != "Uniprot/SWISSPROT" or rec["kind"] not in ("gene_to_db", "old_gene_to_db"):
            continue
        acc = C.bm(rec["to_release"], "swissprot_acc", 1)
        if acc is None:
            continue
        E = sorted(acc[0].get(gene_base(rec["id"]), set()))
        if E != sorted(rec["expect"]):
            out[rec["kind"]].append(revised(rec, "swissprot75", "Swiss-Prot by accession (the corrected judge)", expect=E))
    return out


def rev_living_ids(org, lines, C):
    """Rule `id_lives_on`: Ensembl records no self event beside an event of another id into or out of an id that lives on,
    so the strict walk keeps an id present in the GTF on both sides of such a session (and the broad walk adds it).
    Every line whose truth is a walk is recomputed; rule `judge_copies_apart`'s copies apart are applied to the new set."""
    out = collections.defaultdict(list)
    regions = {}
    for rec in lines:
        k = rec["kind"]
        if rec.get("final_database") or k == "unchanged_old_to_current":
            continue
        w = _history_walks(C, rec, living=True)
        if w is None:
            continue
        E, B, rel = w
        if k in ("name_old_to_current", "name_changed_gene") and not E:
            continue  # their broad then holds the genes with the name now: left as published
        Ev, Bv = _vers(C, E, rel, rec["cmp"]), _vers(C, B, rel, rec["cmp"])
        old_E, old_B = sorted(rec["expect"]), sorted(rec["broad"] or rec["expect"])
        newB = sorted(set(Bv) | set(old_B)) if (rec["broad"] is not None and not E) else Bv
        if Ev == old_E and (newB == old_B or set(newB) == set(old_E)):
            continue
        ch = dict(expect=Ev, broad=newB if set(newB) != set(Ev) else None)
        apart = None
        if k == "name_old_to_current":
            ch["rule"] = "set_exact" if len(Ev) > 1 else "set"
        if rel not in regions:
            where = {g: r["region"] for g, r in (H.gtf_genes(org, rel) or {}).items()}
            regions[rel] = where if "patch" in where.values() else {}
        where = regions[rel]
        if where and len(Ev) > 1:  # rule `judge_copies_apart`, as copies_apart applies it
            sides = collections.defaultdict(list)
            for gene in Ev:
                sides[where.get(gene_base(gene))].append(gene)
            if None not in sides and sides["primary"] and sides["patch"]:
                ch["expect"], apart = sorted(sides["primary"]), sorted(sides["patch"])
                if ch.get("rule", rec["rule"]) == "set_exact" and len(ch["expect"]) == 1:
                    ch["rule"] = "set"
        new = revised(rec, "living_ids", "strict walk keeps a living id (the corrected judge after rule `id_lives_on`)", **ch)
        new.pop("apart", None)
        if apart:
            new["apart"] = apart
        if sorted(new["expect"]) == old_E and sorted(new.get("apart") or []) == sorted(rec.get("apart") or []) and \
                (new["broad"] or []) == (rec["broad"] or []):
            continue
        out[k].append(new)
    return out


def rev_split_name_set(org, lines, C):
    """name_old_to_current (and travel_name): one gene carried the name at R and Ensembl split it; `best` may return
    any non-empty part of the successors, as `retired_split` judges the same history (rule `set`, not `set_exact`)."""
    out = collections.defaultdict(list)
    for rec in lines:
        if rec["kind"] not in ("name_old_to_current", "travel_name") or rec["rule"] != "set_exact" or \
                rec["from_release"] is None:
            continue
        R = rec["from_release"]
        src = C.display(R)[0].get(rec["id"], set()) & C.bm_genes(R)
        if len(src) == 1:
            out[rec["kind"]].append(revised(rec, "split_name_set", "one gene had the name at R: rule set (the corrected judge)",
                                            rule="set"))
    return out


def rev_names_broad(org, lines, C):
    """Names in other name sources are `broad` (the truth cannot decide between genes that carry a string as a
    name): the name kinds add the genes with the string as VGNC symbol or Vega gene name (their broad already holds
    BioMart's name tables); case and separator variants get the genes that carry the original name in any name
    source; name_old_to_current, name_changed_gene and travel_name add the broad walks of the genes that carried the
    string at R in any name source (and, read at the target, the genes that carry it there). Synonyms are no names here
    (rule `names_before_synonyms`: names before synonyms)."""
    out = collections.defaultdict(list)
    T = C.T
    for rec in lines:
        k = rec["kind"]
        if rec.get("final_database") or rec.get("excluded"):
            continue
        base = set(rec["broad"] or rec["expect"])
        cmp = rec["cmp"]
        if k in ("name_unique_to_gene", "name_ambiguous_to_set", "synonym_to_gene"):
            rel = rec["from_release"] or rec["to_release"]
            add = set()
            for pp in (C.vgnc(rel), C.vega(rel)):
                if pp:
                    add |= pp[1].get(rec["id"], set())
            add &= set(C.gtf(rec["to_release"]))
            new = base | set(_vers(C, add, rec["to_release"], cmp))
        elif k in ("case_variant", "separator_variant"):
            v = _name_of_variant(rec)
            if v is None or ENS_ID.match(v):
                continue
            add = C.name_genes(T, v) & set(C.gtf(T))
            new = base | set(_vers(C, add, T, cmp))
        elif k in ("name_old_to_current", "name_changed_gene", "travel_name"):
            R = rec["from_release"]
            if R is None:  # name_changed_gene by default: the string read at the target
                add = C.name_genes(T, rec["id"]) & set(C.gtf(T))
                new = base | set(_vers(C, add, T, cmp))
            else:
                n2g = C.display(R)[0]
                src = n2g.get(rec["id"], set()) & C.bm_genes(R) if k != "name_changed_gene" else n2g.get(rec["id"], set())
                others = (C.name_genes(R, rec["id"]) & set(C.tr.genes.get(R, {}))) - src
                tgt = rec["to_release"]  # T for the name kinds; any release from R on for travel_name
                add = set()
                for g in others:
                    add |= C.tr.walk_fwd(g, R, tgt, False, living=True) if tgt > R else ({g} if tgt == R else set())
                new = base | set(_vers(C, add, tgt, cmp))
        else:
            continue
        if new != base and new != set(rec["expect"]):
            out[k].append(revised(rec, "names_broad", "names in other name sources in broad (the corrected judge)",
                                  broad=sorted(new)))
    return out


def rev_spellings_once(org, lines, C):
    """Two spellings of one value of one database for one gene (`IFN-DELTA-8`, `IFN-delta-8`) are one value: the
    line lists them in `spellings`, and an answer with either spelling counts as that value."""
    out = collections.defaultdict(list)
    for rec in lines:
        if rec["kind"] not in ("gene_to_db", "old_gene_to_db") or len(rec["expect"]) < 2:
            continue
        groups = collections.defaultdict(list)
        for v in rec["expect"]:
            groups[v.casefold()].append(v)
        sp = sorted(sorted(g) for g in groups.values() if len(g) > 1)
        if sp and rec.get("spellings") != sp:
            out[rec["kind"]].append(revised(rec, "spellings_once", "two spellings of one value counted once "
                                            "(the corrected judge)", spellings=sp))
    return out


def rev_projected_hgnc(org, lines, C):
    """Mouse: a display name that is a human HGNC symbol projected onto a mouse gene (display source `HGNC Symbol`)
    and no MGI symbol, UniProt gene name, NCBI label, WikiGene name or synonym of any mouse gene in BioMart at the
    release read asks for a database the mouse graph does not have (as `harness_lib.EXCLUDED` says for
    `gene_to_db`): the line is left out (`excluded`)."""
    out = collections.defaultdict(list)
    if org != "mus_musculus":
        return out
    kinds = ("name_unique_to_gene", "name_ambiguous_to_set", "case_variant", "separator_variant", "name_old_to_current",
             "name_changed_gene")
    for rec in lines:
        if rec["kind"] not in kinds or rec.get("final_database") or rec.get("excluded"):
            continue
        v = _name_of_variant(rec) if rec["kind"] in ("case_variant", "separator_variant") else rec["id"]
        if v is None or ENS_ID.match(v):
            continue
        rel = rec["from_release"] or rec["to_release"]
        n2g, g2n = C.display(rel)
        genes = n2g.get(v, set())
        if not genes or any(g2n[g][1] != "HGNC Symbol" for g in genes):
            continue
        other = set()
        for grp, col in [("mgi", 1), ("uniprot_gn", 1), ("entrez", 2), ("wikigene", 1), ("synonym", 1)]:
            pp = C.bm(rel, grp, col)
            if pp:
                other |= pp[1].get(v, set())
        if other:
            continue
        out[rec["kind"]].append(revised(rec, "projected_hgnc", "left out: projected HGNC symbol (the corrected judge)",
                                        excluded=f"{v} at {rel}: display name from HGNC Symbol only (a human symbol "
                                                 f"projected onto the mouse gene); no MGI, UniProt, NCBI, WikiGene or "
                                                 f"synonym entry"))
    return out


def rev_support(org, lines, C):
    """Rule `support_ranks_shared_name`: a `best` answer that is a non-empty part of the expected set whose genes carry the name in
    strictly more name sources than every gene left out is `undetermined`. The line holds the support of every
    expected gene (`support`: BioMart's name sources at the release read; for name_old_to_current the largest
    support at R of the genes whose strict walk reaches it); score.py applies the rule. Only lines whose expected
    genes differ in support."""
    out = collections.defaultdict(list)
    T = C.T
    for rec in lines:
        k = rec["kind"]
        if k not in ("name_ambiguous_to_set", "synonym_to_gene", "name_old_to_current") or rec["strategy"] != "best" \
                or rec["rule"] != "set_exact" or len(rec["expect"]) < 2 or rec.get("final_database"):
            continue
        by_base = {gene_base(e): e for e in rec["expect"]}
        sup = {e: 0 for e in rec["expect"]}
        if k == "name_old_to_current":
            R = rec["from_release"]
            src = C.display(R)[0].get(rec["id"], set()) & C.bm_genes(R)
            for g in src:
                n = C.support(R, g, rec["id"])
                for t in C.tr.walk_fwd(g, R, T, True, living=True):
                    if t in by_base:
                        sup[by_base[t]] = max(sup[by_base[t]], n)
        else:
            rel = rec["from_release"] or rec["to_release"]
            for b, e in by_base.items():
                sup[e] = C.support(rel, b, rec["id"])
        if len(set(sup.values())) > 1 and rec.get("support") != sup:
            out[k].append(revised(rec, "support", "support per gene for rule `support_ranks_shared_name`'s rule (the corrected judge)", support=sup))
    return out


def rev_resurrection(org, lines, C):
    """Rule `resurrection`: an id that Ensembl retired and later brought back under the same stable id is one gene
    throughout; the ids that the walks reach only through it (with `resurrect`, not without) go to `broad`."""
    out = collections.defaultdict(list)
    for rec in lines:
        if rec.get("final_database") or rec["kind"] == "unchanged_old_to_current":
            continue
        w = _history_walks(C, rec, living=True, resurrect=True)
        if w is None:
            continue
        E, B, rel = w
        E0, B0, _ = _history_walks(C, rec, living=True, resurrect=False)
        add = (B | E) - (B0 | E0)
        base = set(rec["broad"] or rec["expect"])
        new = base | set(_vers(C, add, rel, rec["cmp"]))
        if new != base:
            out[rec["kind"]].append(revised(rec, "resurrection", "a returned id in broad (the corrected judge after "
                                            "rule `resurrection`)", broad=sorted(new)))
    return out


def _lower_map(v2g):
    low = collections.defaultdict(set)
    for v, gs in v2g.items():
        low[v.lower()] |= gs
    return low


def rev_capitals_mgi(org, lines, C):
    """Rule `single_case_symbols` (strings written entirely in capitals): in an organism whose authority
    writes symbols in mixed letter case (`mixed_case_symbols`: MGI's), a string written entirely in capitals that is no
    exact symbol and is, without letter case, the symbol of exactly one gene is read as that symbol first: the line
    expects that gene (read at an older release: its strict walk), and the genes it expected go to `broad`, as rule
    `name_before_synonym` made the judge follow rule `spelling_of_release`."""
    return _single_case(org, lines, C, "capitals_mgi", "upper")


def rev_places(org, lines, C):
    """Places across assemblies: a line judged by the locus rule (Ensembl's history has no successor) whose answer
    lies on another assembly than the query's place gets the place mapped through Ensembl's chain files
    (`locus_mapped`, or the reason it cannot be mapped); score.py then applies the overlap rule as on one assembly."""
    from idtrack_reproducibility.harness import assembly_map
    out = collections.defaultdict(list)
    asm = H.gtf_assemblies(org)
    mapper = C._cached("mapper", lambda: assembly_map.Mapper(org))
    for rec in lines:
        loc = rec.get("locus")
        if not loc or rec["expect"] or rec.get("locus_mapped"):
            continue
        to_asm = asm.get(rec["to_release"])
        if not to_asm or to_asm == loc["assembly"]:
            continue
        out[rec["kind"]].append(revised(rec, "places", "place mapped across assemblies (the corrected judge)",
                                        locus_mapped=mapper.map(loc, to_asm)))
    return out


def _walk_source(C, rec):
    """(gene, release read, target) of a line whose truth is a walk of its own id, or None."""
    k, tr, T = rec["kind"], C.tr, C.T
    if not (k.startswith("travel_") or k.startswith("retired_") or k.startswith("reverse_") or k == "cross_assembly") \
            or k == "travel_name":
        return None
    g = gene_base(rec["id"])
    if g not in tr.present:
        return None
    t = rec["to_release"]
    if k.startswith("reverse_"):
        return g, T, t
    if k.startswith("retired_"):
        return g, max(tr.present[g]), T
    R = rec["from_release"]
    if R is None:
        R = t if t in tr.present[g] else max(tr.present[g])
    return g, R, t


def rev_sessions_effective(org, lines, C):
    """Mapping sessions are crossed where their effect appears in the GTFs (`Truth.effective_release`), not between
    their old release and the next: rhesus 50 -> 52 and 52 -> 58 take effect at 52 and 58 (the GTFs of 50-51 and of
    52-57 hold the same gene ids). Every line whose truth is a walk of its id (travel, retired, reverse,
    cross-assembly; travel lines to a name database) is walked under both placements (strict with living ids; broad
    with living and returned ids); where they differ the line's expected genes (values), broad and place are
    computed again under the effective placement."""
    out = collections.defaultdict(list)
    old = C._cached("truth_old_release", lambda: Truth(org, C.T, placement="old_release"))
    if all(old.eff[s] == C.tr.eff[s] for s in C.tr.sessions):
        return out
    mapper = None
    for rec in lines:
        src = _walk_source(C, rec)
        if src is None:
            continue
        g, R, t = src
        walks = [(tr.walk_series(g, R, [t], True, living=True)[t],
                  tr.walk_series(g, R, [t], False, living=True, resurrect=True)[t]) for tr in (old, C.tr)]
        if walks[0] == walks[1]:
            continue
        E, B = walks[1][0], walks[1][1] | walks[1][0]
        ch = {}
        if rec.get("final_database"):
            grp, col = GENE_DB[rec["final_database"]]
            pp = C.bm(t, grp, col)
            vals = set().union(*(pp[0].get(x, set()) for x in E)) if pp else set()
            bvals = set().union(*(pp[0].get(x, set()) for x in B)) if pp else set()
            ch = dict(expect=sorted(vals), broad=sorted(bvals) if bvals != vals else None)
        else:
            Ev, Bv = _vers(C, E, t, rec["cmp"]), _vers(C, B, t, rec["cmp"])
            ch = dict(expect=Ev, broad=Bv if set(Bv) != set(Ev) else None,
                      locus=C.tr.locus(g, R) if not E else None)
            if not E and ch["locus"] and H.gtf_assemblies(org).get(t) != ch["locus"]["assembly"]:
                if mapper is None:
                    from idtrack_reproducibility.harness import assembly_map
                    mapper = assembly_map.Mapper(org)
                ch["locus_mapped"] = mapper.map(ch["locus"], H.gtf_assemblies(org).get(t))
        if all(rec.get(k) == v for k, v in ch.items()):
            continue
        new = revised(rec, "sessions_effective", "sessions crossed where their effect appears (the corrected judge, 2026-11-18)",
                      **ch)
        if E:
            new.pop("locus_mapped", None)
        out[rec["kind"]].append(new)
    return out


def rev_exceptions(org, lines, C):
    """Places on haplotypes, patches, alternative loci and the Y PARs compared on the reference sequence: the locus
    rule compared places by sequence name, so a copy on such a sequence and the reference gene it stands for (same
    chromosome coordinates, another sequence name) counted as two places. Ensembl's assembly exceptions of the release
    (`exceptions/<release>.tsv`, fetch_tables.py exceptions) put the query's place on its reference sequence
    (`locus_projected`; `locus_mapped` mapped again from it across assemblies), and score.py projects every returned
    gene's place at the target release the same way (`exceptions`: true). Every locus line of an organism whose
    releases have assembly exceptions, where the query's release or the target release has any, is revised."""
    from idtrack_reproducibility.harness import assembly_map
    out = collections.defaultdict(list)
    exc = C._cached("exceptions", lambda: assembly_map.Exceptions(org))
    mapper = None
    asm = H.gtf_assemblies(org)
    for rec in lines:
        loc = rec.get("locus")
        if not loc or rec["expect"] or rec.get("exceptions"):
            continue
        rel, t = loc["release"], rec["to_release"]
        if not (exc.has(rel) or exc.has(t)):
            continue
        proj = exc.project(loc, rel)
        ch = dict(exceptions=True)
        if proj is not loc:
            ch["locus_projected"] = proj
            to_asm = asm.get(t)
            if to_asm and to_asm != loc["assembly"]:
                if mapper is None:
                    mapper = C._cached("mapper", lambda: assembly_map.Mapper(org))
                ch["locus_mapped"] = dict(unmapped=proj["unplaced"]) if proj.get("unplaced") else mapper.map(proj, to_asm)
        out[rec["kind"]].append(revised(rec, "exceptions", "places on haplotypes, patches and PARs projected on the "
                                        "reference sequence (the corrected judge, 2026-11-18)", **ch))
    return out


def _reproject(org, lines, C, name, note):
    """Lines of `exceptions` whose query place, projected as assembly_map.Exceptions.project does now, differs from
    their `locus_projected`: the new place (and `locus_mapped` mapped again from it)."""
    from idtrack_reproducibility.harness import assembly_map
    out = collections.defaultdict(list)
    exc = C._cached("exceptions", lambda: assembly_map.Exceptions(org))
    asm = H.gtf_assemblies(org)
    mapper = None
    for rec in lines:
        if not rec.get("exceptions") or not rec.get("locus_projected"):
            continue
        loc = rec["locus"]
        proj = exc.project(loc, loc["release"])
        if proj == rec["locus_projected"]:
            continue
        ch = dict(locus_projected=proj)
        to_asm = asm.get(rec["to_release"])
        if to_asm and to_asm != loc["assembly"]:
            if mapper is None:
                mapper = C._cached("mapper", lambda: assembly_map.Mapper(org))
            ch["locus_mapped"] = dict(unmapped=proj["unplaced"]) if proj.get("unplaced") else mapper.map(proj, to_asm)
        out[rec["kind"]].append(revised(rec, name, note, **ch))
    return out


def rev_exceptions_span(org, lines, C):
    """Lines of `exceptions` whose query place, projected as assembly_map.Exceptions.project does now, differs from
    their `locus_projected` get the new place (and `locus_mapped` mapped again from it). The projection was refined
    twice after `exceptions` was published, each published under this name: the span of both segment ends instead of
    linear scaling (human v14, mouse v17), the GRC's alignments of the scaffolds (human v15, mouse v18), and no place
    where no alignment holds an end inside the segment (human v16, mouse v19)."""
    return _reproject(org, lines, C, "exceptions_span", "copy's place re-projected (the corrected judge, 2026-11-18)")


def rev_core_places(org, lines, C):
    """Genes of the core database that the GTF of the release lacks (human 48-51, cynomolgus 112) have a place
    (`core_places/<release>.tsv`, from the core's gene and seq_region tables): a locus line whose target release has
    such genes gets `core_places: true`, and score.py places a returned gene absent from the GTF there instead of
    calling it 'not at the query's place' for want of a place."""
    out = collections.defaultdict(list)
    have = {int(f[:-4]) for f in os.listdir(os.path.join(H.truth_dir(org), "core_places"))} \
        if os.path.isdir(os.path.join(H.truth_dir(org), "core_places")) else set()
    for rec in lines:
        if rec.get("locus") and not rec["expect"] and rec["to_release"] in have and not rec.get("core_places"):
            out[rec["kind"]].append(revised(rec, "core_places", "genes the GTF lacks placed from the core (decision "
                                            "78, 2026-11-19)", core_places=True))
    return out


# ---------------- the judge following idtrack's decided rules (docs/rules.md; revision sets of judge `rules`) ----------------

# Rule `missing_values`: the spellings of a missing value that the decided rule refuses (compared without letter case
# and surrounding white space; the empty string too): the strings pandas reads as missing by default (its own list,
# read here) and the two the rule adds, pandas' text of a missing date or time and UniProt's placeholder name (no
# source lists those two; the rule names them, docs/rules.md).
MISSING_ADDED = frozenset({"nat", "unknown"})


def missing_values():
    """Rule `missing_values`'s spellings, casefolded: pandas' default missing values and `MISSING_ADDED`."""
    from pandas._libs.parsers import STR_NA_VALUES
    return frozenset(v.casefold() for v in STR_NA_VALUES if v) | MISSING_ADDED
# Rule `stray_punctuation`: the ASCII members of Unicode's property Terminal_Punctuation (the property is not exposed
# by `unicodedata`), and the quotation marks: Unicode's initial and final quotes (categories Pi, Pf) and ASCII's two.
TERMINAL_PUNCTUATION = "!,.:;?"
ASCII_QUOTES = "'\""


def is_missing_value(s):
    """Rule `missing_values`: the string stands for a missing value."""
    return not s.strip() or s.strip().casefold() in missing_values()


def strip_stray_punctuation(s):
    """Rule `stray_punctuation`: the string without the punctuation that ends a sentence at its end and without a pair of
    quotation marks around it, repeatedly (a quote at one end only stays); the string itself when nothing would be left."""
    import unicodedata

    def quote(ch):
        return ch in ASCII_QUOTES or unicodedata.category(ch) in ("Pi", "Pf")
    cur = s.strip()
    while True:
        before = cur
        cur = cur.rstrip(TERMINAL_PUNCTUATION).rstrip()
        if len(cur) > 2 and quote(cur[0]) and quote(cur[-1]):
            cur = cur[1:-1].strip()
        if cur == before:
            break
    return cur or s


def mixed_case_symbols(C, rel):
    """(gene -> {symbol}, symbol -> {gene}) of the organism's nomenclature authority at a release when that authority
    writes its symbols in mixed letter case, else None. The authority is the source of most genes' display names in
    BioMart's gene table (MGI Symbol for mouse, HGNC Symbol for human and pig, VGNC Symbol for the macaques); it writes
    in mixed case when more than half of those display names hold a small letter. Rule `single_case_symbols` applies
    to organisms of that kind only (idtrack reads the same convention from its graph's symbols)."""
    def make():
        g2n = C.display(rel)[1]
        sources = collections.Counter(src for name, src in g2n.values() if name and src)
        if not sources:
            return None
        src = sources.most_common(1)[0][0]
        names = [name for name, s in g2n.values() if name and s == src]
        if sum(1 for n in names if any(ch.islower() for ch in n)) * 2 <= len(names):
            return None
        if src in GENE_DB:
            pp = C.bm(rel, *GENE_DB[src])
        elif src == "VGNC Symbol":
            pp = C.vgnc(rel)
        else:
            pp = None
        return pp if pp and pp[1] else None
    return C._cached(("mixed_case_symbols", rel), make)


def _single_case(org, lines, C, name, case):
    """Rule `single_case_symbols` for strings written entirely in one letter case (`upper` or `lower`): see
    rev_capitals_mgi and rev_lower_case_symbols."""
    out = collections.defaultdict(list)
    for rec in lines:
        s = rec["id"]
        if rec.get("final_database") or rec.get("excluded") or ENS_ID.match(s) or not any(ch.isalpha() for ch in s) or \
                s != (s.upper() if case == "upper" else s.lower()):
            continue
        rel = rec["from_release"] or rec["to_release"]
        mixed = mixed_case_symbols(C, rel)
        if not mixed or s in mixed[1]:
            continue
        low = C._cached(("mixed_low", rel), lambda: _lower_map(mixed[1]))
        genes = low.get(s.lower(), set()) & C.primary(rel)
        if len(genes) != 1:
            continue
        (g,) = genes
        tgt = rec["to_release"]
        new_E = {g} if tgt == rel else C.tr.walk_fwd(g, rel, tgt, True, living=True)
        if not new_E or new_E == {gene_base(e) for e in rec["expect"]}:
            continue
        Ev = _vers(C, new_E, tgt, rec["cmp"])
        broad = set(rec["broad"] or []) | set(rec["expect"]) | set(Ev)
        words = "capitals" if case == "upper" else "lower case"
        out[rec["kind"]].append(revised(rec, name, f"{words} read as the symbol of {g} (the corrected judge after rule "
                                        f"`single_case_symbols`)", expect=Ev, rule="set_exact" if len(Ev) > 1 else "set",
                                        broad=sorted(broad)))
    return out


def rev_lower_case_symbols(org, lines, C):
    """Rule `single_case_symbols` (decided; strings written entirely in lower case): in an organism whose
    authority writes symbols in mixed letter case (`mixed_case_symbols`), a string written entirely in lower case that
    is no exact symbol and is, without letter case, the symbol of exactly one gene is read as that symbol first: the
    line expects that gene (read at an older release: its strict walk), and the genes it expected go to `broad`, as
    `capitals_mgi` does for capitals."""
    return _single_case(org, lines, C, "lower_case_symbols", "lower")


def rev_stray_punctuation(org, lines, C):
    """Rule `stray_punctuation` (decided): punctuation that ends a sentence at the end of a string and quotes
    around it are no part of an identifier, so a line whose string carries them asks for the identifier without them
    (`strip_stray_punctuation`). For a line that expects Ensembl genes of a string read as a name or synonym at the
    release asked about (from_release, else to_release), the expected genes are those that carry the string, with or
    without the punctuation, as a name there (rule `name_before_synonym`), else as a synonym (BioMart's display names,
    name tables and synonyms, VGNC and Vega names: `Ctx78.name_genes`), genes of the primary sequence; the genes it
    expected before go to `broad`. Lines that ask for another database, or whose string has no such punctuation, are
    not touched."""
    out = collections.defaultdict(list)
    for rec in lines:
        s = rec["id"]
        if rec.get("final_database") or rec.get("excluded") or ENS_ID.match(s):
            continue
        clean = strip_stray_punctuation(s)
        if clean == s.strip():
            continue
        rel = rec["from_release"] or rec["to_release"]
        prim = C.primary(rel)
        names = (C.name_genes(rel, s) | C.name_genes(rel, clean)) & prim
        syn = (C.name_genes(rel, s, synonyms=True) | C.name_genes(rel, clean, synonyms=True)) & prim
        genes = names or syn
        tgt = rec["to_release"]
        new_E = genes if tgt == rel else set().union(*(C.tr.walk_fwd(g, rel, tgt, True, living=True) for g in genes))
        if not new_E or new_E == {gene_base(e) for e in rec["expect"]}:
            continue
        Ev = _vers(C, new_E, tgt, rec["cmp"])
        broad = set(rec["broad"] or []) | set(rec["expect"]) | set(Ev)
        out[rec["kind"]].append(revised(rec, "stray_punctuation", f"read as {clean!r} (rule `stray_punctuation`)",
                                        expect=Ev, rule="set_exact" if len(Ev) > 1 else "set", broad=sorted(broad)))
    return out


def rev_missing_values(org, lines, C):
    """Rule `missing_values` (decided): a string that stands for a missing value
    (`missing_values()`) is no identifier and is refused before any lookup, unless its source database is named. The
    harness's lines carry no source database, so such a line (drawn from a table where a placeholder stands as a gene's
    synonym, e.g. a synonym `N/A`) asks for no gene: it is left out of the scores (`excluded`, with the reason). The
    refusal itself is measured by family E (E1, judge `decided`)."""
    out = collections.defaultdict(list)
    for rec in lines:
        s = rec["id"]
        if rec.get("excluded") or rec.get("from_database") or ENS_ID.match(s) or not is_missing_value(s):
            continue
        out[rec["kind"]].append(revised(rec, "missing_values", "a missing value, no identifier (rule `missing_values`)",
                                        excluded=f"{s!r} stands for a missing value (rule `missing_values`: refused; "
                                                 "family E measures the refusal)"))
    return out


def rev_par_shared_names(org, lines, C):
    """Rule `par_shared_names` (decided): a name that the X and Y copies of a pseudoautosomal gene share is
    the X gene's, with the Y copy reported apart. Where a line asked by a name expects Ensembl genes and its expected
    set at the target release holds a gene placed in a pseudoautosomal region of a sequence (Ensembl's assembly
    exceptions of the release, type PAR: the region of Y and the segment of X it equals) together with a gene of the
    segment it equals (the place projected by its offset overlaps that gene's place, GTF places), the copy leaves
    `expect` and is expected in the list reported apart (`apart`), as rule `judge_copies_apart` does for copies on
    patches; `rule` `set_exact` stays only while more than one gene is expected; `broad` and the rest are kept."""
    from idtrack_reproducibility.harness import assembly_map
    out = collections.defaultdict(list)
    exc = C._cached("exceptions", lambda: assembly_map.Exceptions(org))
    for rec in lines:
        if rec.get("final_database") or rec.get("excluded") or ENS_ID.match(rec["id"]) or len(rec["expect"]) < 2:
            continue
        tgt = rec["to_release"]
        par = {sq: rows for sq, rows in exc.rows(tgt).items() if any(r["type"] == "PAR" for r in rows)}
        if not par:
            continue
        places = C.gtf(tgt)
        base = {e: gene_base(e) for e in rec["expect"]}
        copies = set()
        for e, g in base.items():
            p = places.get(g)
            if not p or p["seqname"] not in par:
                continue
            proj = exc.project(dict(p), tgt)
            if proj is p or proj.get("exception") != "PAR":
                continue
            if any(o != e and assembly_map.overlaps(proj, places.get(h)) for o, h in base.items()):
                copies.add(e)
        if not copies or copies == set(rec["expect"]):
            continue
        keep = sorted(set(rec["expect"]) - copies)
        apart = sorted(set(rec.get("apart") or []) | copies)
        rule = "set" if rec["rule"] == "set_exact" and len(keep) == 1 else rec["rule"]
        out[rec["kind"]].append(revised(rec, "par_shared_names", "Y copies apart (rule `par_shared_names`)", expect=keep,
                                        apart=apart, rule=rule))
    return out


def rev_copies_apart_early(org, lines, C):
    """Rule `judge_copies_apart` before release 82: Ensembl's GTFs before 82 have no file of the genes on patches and
    haplotypes, so `copies_apart` could not tell copies there and such lines still expect the copies inside the answer.
    Ensembl's assembly exceptions of every release (`exceptions/<release>.tsv`) name the sequences that are haplotypes
    and patches: a line that expects Ensembl genes at a target release before 82, holds no `apart` yet, and expects genes
    both on such sequences (GTF places) and off them expects the genes off them, the others reported apart (`apart`);
    `rule` `set_exact` stays only while more than one gene is expected. The pseudoautosomal regions are no copies here
    (rule `par_shared_names` reads them)."""
    from idtrack_reproducibility.harness import assembly_map
    out = collections.defaultdict(list)
    exc = C._cached("exceptions", lambda: assembly_map.Exceptions(org))
    for rec in lines:
        t = rec["to_release"]
        if t is None or t >= 82 or rec.get("final_database") or rec.get("excluded") or len(rec["expect"]) < 2 or \
                rec.get("apart"):
            continue
        places, rows = C.gtf(t), exc.rows(t)
        side = {}
        for e in rec["expect"]:
            g = places.get(gene_base(e)) if ENS_ID.match(e) else None
            if g is None:
                break
            side[e] = any(r["type"] != "PAR" for r in rows.get(g["seqname"], ()))
        else:
            if all(side.values()) or not any(side.values()):
                continue
            keep = sorted(e for e, v in side.items() if not v)
            rule = "set" if rec["rule"] == "set_exact" and len(keep) == 1 else rec["rule"]
            out[rec["kind"]].append(revised(rec, "copies_apart_early", "copies apart before 82 (rule `judge_copies_apart`, "
                                            "places from Ensembl's assembly exceptions)", expect=keep,
                                            apart=sorted(e for e, v in side.items() if v), rule=rule))
    return out


# name -> function (its docstring says what it revises and why), in the order they are published
REVISIONS = {
    "swissprot75": rev_swissprot75,
    "living_ids": rev_living_ids,
    "split_name_set": rev_split_name_set,
    "names_broad": rev_names_broad,
    "spellings_once": rev_spellings_once,
    "projected_hgnc": rev_projected_hgnc,
    "support": rev_support,
    "resurrection": rev_resurrection,
    "capitals_mgi": rev_capitals_mgi,
    "places": rev_places,
    "sessions_effective": rev_sessions_effective,
    "exceptions": rev_exceptions,
    "exceptions_span": rev_exceptions_span,
    "core_places": rev_core_places,
    "missing_values": rev_missing_values,
    "lower_case_symbols": rev_lower_case_symbols,
    "stray_punctuation": rev_stray_punctuation,
    "par_shared_names": rev_par_shared_names,
    "copies_apart_early": rev_copies_apart_early,
}
# The judge of each correction's revision set (its manifest's `judge`): `78`, the judge corrected where it was wrong
# (rule `judge_corrections`); `rules`, the judge following a rule idtrack decided after 71690ab (docs/rules.md). Both are left out
# by the judge `published` (harness_lib.REVISED_JUDGES).
RULES_JUDGE = {"missing_values", "lower_case_symbols", "stray_punctuation", "par_shared_names"}


def not_revised_after(name, lines):
    """The lines a correction may revise: none that a correction later in the order has revised (a correction run
    again must not undo what a later one did, e.g. `living_ids` recomputing a `broad` that `names_broad` widened)."""
    later = set(list(REVISIONS)[list(REVISIONS).index(name) + 1:])
    return [r for r in lines if not later & set(r.get("revised_by") or ())]


def apply_revisions(org, lines, C, names=None, report=None):
    """Apply the revisions of the corrected judge in order to a list of lines; returns the revised list (for new sets)."""
    by = {r["qid"]: r for r in lines}
    for name in names or REVISIONS:
        out = REVISIONS[name](org, not_revised_after(name, by.values()), C)
        for v in out.values():
            for rec in v:
                by[rec["qid"]] = rec
        if report is not None:
            report[name] = {k: len(v) for k, v in out.items()}
    return list(by.values())


def publish_revisions(a, d, names, C=None):
    """Publish one revision set per correction of the corrected judge (`--revise`), each under the next free set name and
    holding only the lines it changes; a correction that changes nothing publishes nothing. Every correction leaves
    a line it has already revised as it is, so `--revise all` run again publishes nothing (or only what an
    interrupted run left out)."""
    lines, based_on = published_lines(d)
    C = C or Ctx78(a.organism, a.target or H.ORGANISMS[a.organism][2])
    summary = {}
    for name in names:
        func = REVISIONS[name]
        out = func(a.organism, not_revised_after(name, lines.values()), C)
        counts = {k: len(v) for k, v in sorted(out.items())}
        summary[name] = counts
        print(name, json.dumps(counts), flush=True)
        for v in out.values():  # a later correction reads the lines as the earlier ones revised them
            for rec in v:
                lines[rec["qid"]] = rec
        if a.dry_run or not counts:
            continue
        n = max(int(x[1:]) for x in H.set_names(d) if re.match(r"^v\d+$", x)) + 1
        set_name = f"v{n}"
        path = os.path.join(d, f"{set_name}.jsonl.gz")
        if os.path.exists(path):
            raise SystemExit(f"{path} exists")
        tmp = path + ".tmp"
        with gzip.open(tmp, "wt") as fo:
            for k in sorted(out):
                for rec in out[k]:
                    fo.write(H.dumps(rec) + "\n")
        os.replace(tmp, path)
        sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
        man = dict(organism=a.organism, set=set_name, seed=SEED, kinds=counts, n=sum(counts.values()), sha256=sha,
                   revises=dict(based_on), judge="rules" if name in RULES_JUDGE else "78", correction=name,
                   revision=" ".join((func.__doc__ or "").split()))
        json.dump(man, open(os.path.join(d, f"{set_name}.manifest.json"), "w"), indent=1)
        based_on[set_name] = sha
        print("wrote", path, man["n"], sha, flush=True)
    return summary


# ---------------- the travel matrix (docs/harness.md, "The travel matrix") ----------------

TRAVEL_IDS = 50       # ids per history stratum
TRAVEL_ORIGINS = 5    # releases of origin per id (an even subset of the releases where it exists)
TRAVEL_NAMES = 150    # display names at an old release, travelled to every later release
NAME_DB = {"homo_sapiens": "HGNC Symbol", "mus_musculus": "MGI Symbol"}  # others: UniProtKB Gene Name


def _even(rr, k):
    """k releases spread evenly over a sorted list (its first and last among them)."""
    if len(rr) <= k:
        return list(rr)
    return sorted({rr[round(i * (len(rr) - 1) / (k - 1))] for i in range(k)})


def travel_strata(org, tr, T):
    """Gene ids of the window by their history: {stratum: [ids]} (an id in one stratum only)."""
    rels = tr.rels
    eras = tr.eras()
    strata = collections.defaultdict(list)
    taken = set()
    for gid, rr in tr.present.items():
        if T in rr:
            if len(rr) == len(rels):
                strata["unchanged"].append(gid)
            elif min(rr) > rels[0]:
                strata["born_late"].append(gid)
            continue
        L = max(rr)
        E = tr.walk_fwd(gid, L, T, True, living=True)
        if not E:
            cls = "retired_none" if not tr.walk_fwd(gid, L, T, False, living=True) else "retired_branch"
        elif len(E) > 1:
            cls = "split"
        else:
            merged = False
            for s in tr.sessions_between(L, T):
                news = tr.succ[s].get(gid, set()) - {None}
                if news and gid not in news:
                    merged = any(len(tr.pred[s][n] - {None}) > 1 for n in news)
                    break
            cls = "merged" if merged else "renamed"
        strata[cls].append(gid)
    chosen = {}
    for k in sorted(strata):
        chosen[k] = sample(org, f"travel:{k}", strata[k], TRAVEL_IDS)
        taken |= set(chosen[k])
    for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
        k = f"handoff_{hi1}_{lo2}"
        pop = [g for g, rr in tr.present.items() if hi1 in rr and lo2 in rr and g not in taken]
        chosen[k] = sample(org, f"travel:{k}", pop, TRAVEL_IDS)
        taken |= set(chosen[k])
    return chosen


def build_travel(org, T, C):
    """The travel matrix: every sampled id asked from up to TRAVEL_ORIGINS releases where it exists to every release
    of the window, forward, backward and in place; `best` for all, `all` for every other id; a final database of a
    name for every sixth id, at target releases with a BioMart table. Truth: the strict walk with the corrected judge's rules
    (living ids; broad with returned ids), rule `judge_copies_apart`'s copies apart, the locus rule (with the place mapped across
    assemblies by `apply_revisions`). Plus display names at the oldest BioMart release, travelled to every later one."""
    tr = C.tr
    rels = tr.rels
    out = collections.defaultdict(list)
    db = NAME_DB.get(org, "UniProtKB Gene Name")
    grp, col = GENE_DB[db]
    bm_rels = [r for r in H.biomart_releases(org) if r in rels]
    regions = {}

    def apart_of(Ev, rel, cmp):
        if rel not in regions:
            where = {g: r["region"] for g, r in (H.gtf_genes(org, rel) or {}).items()}
            regions[rel] = where if "patch" in where.values() else {}
        where = regions[rel]
        if not where or len(Ev) < 2:
            return Ev, None
        sides = collections.defaultdict(list)
        for gene in Ev:
            sides[where.get(gene_base(gene))].append(gene)
        if None in sides or not sides["primary"] or not sides["patch"]:
            return Ev, None
        return sorted(sides["primary"]), sorted(sides["patch"])

    def lines_for(kind, ident, src, R1, note, with_all, with_db):
        E = {}
        B = {}
        for g in src:
            for t, v in tr.walk_series(g, R1, rels, True, living=True).items():
                E.setdefault(t, set()).update(v)
            for t, v in tr.walk_series(g, R1, rels, False, living=True, resurrect=True).items():
                B.setdefault(t, set()).update(v)
        loc = tr.locus(next(iter(src)), R1) if len(src) == 1 else None
        for R2 in rels:
            e, b = E.get(R2, set()), B.get(R2, set()) | E.get(R2, set())
            cmp = tr.cmp_for(R2)
            Ev, Bv = tr.vset(e, R2), tr.vset(b, R2)
            Ev2, apart = apart_of(Ev, R2, cmp)
            rule = "set"  # an id, or a name of one gene at R (the corrected judge: `set` for an old name of one split gene)
            d = "same" if R2 == R1 else ("forward" if R2 > R1 else "backward")
            common = dict(from_release=R1, cmp=cmp, rule=rule, broad=Bv if set(Bv) != set(Ev) else None,
                          locus=loc if not e else None, note=f"{note}; {d} {R1}->{R2}")
            for st in ("best", "all") if with_all else ("best",):
                rec = q(kind, ident, R2, Ev2, strategy=st, **common)
                if apart:
                    rec["apart"] = apart
                out[kind].append(rec)
            if with_db and e and R2 in bm_rels:
                pp = C.bm(R2, grp, col)
                if pp:
                    vals = set().union(*(pp[0].get(g, set()) for g in e))
                    bvals = set().union(*(pp[0].get(g, set()) for g in b))
                    out[kind].append(q(kind, ident, R2, sorted(vals), from_release=R1, final_database=db, cmp="exact",
                                       broad=sorted(bvals) if bvals != vals else None, note=f"{note}; {d} {R1}->{R2}; {db}"))

    strata = travel_strata(org, tr, T)
    for k, ids in sorted(strata.items()):
        kind = "travel_" + ("handoff" if k.startswith("handoff_") else k)
        for i, g in enumerate(ids):
            rr = sorted(tr.present[g])
            for R1 in _even(rr, TRAVEL_ORIGINS):
                lines_for(kind, g, {g}, R1, f"{k}; origin {R1} ({tr.assembly.get(R1, '')})", with_all=(i % 2 == 0),
                          with_db=(i % 6 == 0))
    # display names of one gene at the oldest BioMart release, travelled to every later release
    R = bm_rels[0] if bm_rels else None
    if R is not None and R < T:
        n2g = C.display(R)[0]
        bmg = C.bm_genes(R)
        clean = lambda v: not ENS_ID.match(v) and not v.isdigit() and v.strip() == v and v != ""  # noqa: E731
        pop = [v for v, gs in n2g.items() if clean(v) and len(gs & bmg) == 1 and len(gs) == 1]
        for i, v in enumerate(sample(org, "travel:name", pop, TRAVEL_NAMES)):
            before = len(out["travel_name"])
            lines_for("travel_name", v, n2g[v] & bmg, R, f"name at {R}", with_all=(i % 2 == 0), with_db=False)
            out["travel_name"][before:] = [x for x in out["travel_name"][before:] if x["to_release"] >= R]
    return out, strata


def published_lines(d):
    """The latest published line of every qid of one organism, and the sets read with their sha256."""
    lines, based_on = {}, {}
    for name in H.set_names(d):
        based_on[name] = json.load(open(os.path.join(d, f"{name}.manifest.json")))["sha256"]
        with gzip.open(os.path.join(d, f"{name}.jsonl.gz"), "rt") as fh:
            for line in fh:
                rec = json.loads(line)
                lines[rec["qid"]] = rec
    return lines, based_on


def publish_names_first(a, d, path):
    """Publish the revision set of rule `name_before_synonym` for one organism (``--names-first``)."""
    lines, based_on = published_lines(d)
    out, releases = names_first(a.organism, list(lines.values()))
    counts = {k: len(v) for k, v in sorted(out.items())}
    strings = sorted({rec["id"] for v in out.values() for rec in v})
    print(json.dumps({"revised lines per kind": counts, "strings": strings, "releases read": releases}, indent=1))
    if a.dry_run or not counts:
        return
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as fo:
        for k in sorted(out):
            for rec in out[k]:
                fo.write(H.dumps(rec) + "\n")
    os.replace(tmp, path)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    man = dict(organism=a.organism, set=a.set, seed=SEED, kinds=counts, n=sum(counts.values()), sha256=sha,
               revises=based_on, revision="name before synonym (rule `name_before_synonym`): a synonym_to_gene line whose string is, "
               "at the release asked about, an exact name of genes other than exactly the genes it expects (of a gene "
               "it does not expect, or of only some of the genes it expects: integration Q16) now expects the genes "
               "that carry the string as a name", strings=strings, name_tables_of_releases=releases)
    json.dump(man, open(os.path.join(d, f"{a.set}.manifest.json"), "w"), indent=1)
    print("wrote", path, man["n"], sha)


def publish_copies_apart(a, d, path):
    """Publish the revision set of rule `judge_copies_apart` for one organism (``--copies-apart``)."""
    lines, based_on = published_lines(d)
    out, marked = copies_apart(a.organism, list(lines.values()))
    counts = {k: len(v) for k, v in sorted(out.items())}
    print(json.dumps({"revised lines per kind": counts, "target releases with patches in the truth": marked}, indent=1))
    if a.dry_run:
        return
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as fo:
        for k in sorted(out):
            for rec in out[k]:
                fo.write(H.dumps(rec) + "\n")
    os.replace(tmp, path)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    man = dict(organism=a.organism, set=a.set, seed=SEED, kinds=counts, n=sum(counts.values()), sha256=sha,
               revises=based_on, revision="copies apart (rule `judge_copies_apart`): expect = genes of the primary sequence, apart = "
               "the copies on patches and alternative loci, for lines whose expected set holds both",
               releases_with_patches=marked)
    json.dump(man, open(os.path.join(d, f"{a.set}.manifest.json"), "w"), indent=1)
    print("wrote", path, man["n"], sha)


def publish_travel(a, d, path, T):
    """Publish the travel matrix of one organism as a new set (its kinds `travel_*` must not exist yet)."""
    C = Ctx78(a.organism, T)
    out, strata = build_travel(a.organism, T, C)
    lines = [rec for k in sorted(out) for rec in out[k]]
    report = {}
    lines = apply_revisions(a.organism, lines, C, report=report)
    counts = collections.Counter(rec["kind"] for rec in lines)
    print(json.dumps({"lines per kind": dict(sorted(counts.items())), "revisions": report,
                      "ids per stratum": {k: len(v) for k, v in strata.items()}}, indent=1), flush=True)
    existing = set()
    for f in os.listdir(d):
        if f.endswith(".manifest.json"):
            existing |= set(json.load(open(os.path.join(d, f)))["kinds"])
    if existing & set(counts) and not a.dry_run:
        raise SystemExit(f"kinds already published: {sorted(existing & set(counts))}")
    if a.dry_run:
        return
    seen = set()
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as fo:
        for rec in sorted(lines, key=lambda r: (r["kind"], r["qid"])):
            if rec["qid"] in seen:
                continue
            seen.add(rec["qid"])
            fo.write(H.dumps(rec) + "\n")
    os.replace(tmp, path)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    man = dict(organism=a.organism, set=a.set, target=T, seed=SEED, kinds=dict(sorted(counts.items())), n=len(seen),
               sha256=sha, travel=dict(ids_per_stratum=TRAVEL_IDS, origins_per_id=TRAVEL_ORIGINS, names=TRAVEL_NAMES,
                                       name_database=NAME_DB.get(a.organism, "UniProtKB Gene Name"),
                                       strata={k: v for k, v in strata.items()}),
               revisions_applied=report, judge="78", gtf_releases=C.tr.rels, core=C.tr.core,
               biomart_releases=H.biomart_releases(a.organism))
    json.dump(man, open(os.path.join(d, f"{a.set}.manifest.json"), "w"), indent=1)
    print("wrote", path, len(seen), sha)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", default=None, help="the new set's name (not with --revise, which numbers its sets)")
    ap.add_argument("--kinds", default="all")
    ap.add_argument("--target", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true", help="print counts, write nothing")
    ap.add_argument("--copies-apart", action="store_true", help="publish a set that revises the truth of the published "
                    "gene-answer lines after rule `judge_copies_apart` (see copies_apart); it repeats their qids and holds only the "
                    "lines that change")
    ap.add_argument("--revise", default=None, help="the corrected judge: publish one revision set per correction named (comma-"
                    "separated, or `all` for every one in order: " + ", ".join(REVISIONS) + "), each under the next free "
                    "set name; with --dry-run only print the lines each would change (cumulatively)")
    ap.add_argument("--travel", action="store_true", help="publish the travel matrix (build_travel) as a new set")
    ap.add_argument("--names-first", action="store_true", help="publish a set that revises the truth of the published "
                    "synonym_to_gene lines after rule `name_before_synonym` (see names_first); nothing is written when no line changes")
    a = ap.parse_args()
    T = a.target or H.ORGANISMS[a.organism][2]
    d = H.query_dir(a.organism)
    os.makedirs(d, exist_ok=True)
    if a.revise:
        return publish_revisions(a, d, list(REVISIONS) if a.revise == "all" else a.revise.split(","))
    if not a.set:
        raise SystemExit("--set is required")
    path = os.path.join(d, f"{a.set}.jsonl.gz")
    if os.path.exists(path) and not a.dry_run:
        raise SystemExit(f"{path} exists: a published set is never rewritten; use a new --set")
    if a.travel:
        return publish_travel(a, d, path, T)
    if a.copies_apart:
        return publish_copies_apart(a, d, path)
    if a.names_first:
        return publish_names_first(a, d, path)
    kinds = set(a.kinds.split(","))
    out, tr = build(a.organism, kinds, T)
    # a set generated now carries the corrections of the corrected judge, as the revision sets give them to published lines
    revised_lines = apply_revisions(a.organism, [rec for k in sorted(out) for rec in out[k]], Ctx78(a.organism, T, tr))
    out = collections.defaultdict(list)
    for rec in revised_lines:
        out[rec["kind"]].append(rec)
    existing = set()
    for f in os.listdir(d):
        if f.endswith(".manifest.json"):
            manifest = json.load(open(os.path.join(d, f)))
            if not manifest.get("revises"):  # a revision set repeats kinds by design; it is no first publication
                existing |= set(manifest["kinds"])
    clash = existing & set(out)
    if clash and not a.dry_run:
        raise SystemExit(f"kinds already published in another set: {sorted(clash)}")
    counts = {k: len(v) for k, v in sorted(out.items())}
    print(json.dumps(counts, indent=1))
    if a.dry_run:
        return
    seen = set()
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as fo:
        for k in sorted(out):
            for rec in out[k]:
                if rec["qid"] in seen:
                    continue
                seen.add(rec["qid"])
                fo.write(H.dumps(rec) + "\n")
    os.replace(tmp, path)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    man = dict(organism=a.organism, set=a.set, target=T, seed=SEED, kinds=counts, n=len(seen), sha256=sha,
               gtf_releases=tr.rels, assemblies={str(k): v for k, v in sorted(tr.assembly.items())}, core=tr.core,
               biomart_releases=H.biomart_releases(a.organism))
    json.dump(man, open(os.path.join(d, f"{a.set}.manifest.json"), "w"), indent=1)
    print("wrote", path, len(seen), sha)


if __name__ == "__main__":
    main()
