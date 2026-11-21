"""Family C, external identifiers (docs/methods.md, C): the query sets of one organism, one part, published once and
never rewritten.

Usage: gen_c.py <organism> --part c1|c2|c3|c4|c5|c6|c7 --set <name> [--target 115] [--limit N] [--dry-run]

Every line has `strategy: both` (scored under `best` and `all`), `family: c`, a `stratum`, and `unit` = the external
identifier without version (C6: the gene). qid `<kind>|<id>|<from or ->|<to>[|<final database>][|db=<from database>]`.
Truth (nothing of idtrack; `truth_c.py`): Ensembl's xrefs of every release of the window (core dumps), the genes they
name at the release read; history walks of those genes to the target (family A's `GenA`: A*'s strict walk with living
ids = `expect`, broad walk = `broad`, copies apart, the locus rule with places mapped across assemblies). The graph's
scope (which databases at which releases) is the reference tree's configuration (`truth_c.Scope`). The authority's
view (each database's own file) rides along as `auth` (base gene ids, or `none`), never as the judge (principle 2).

Releases read ("table releases"): those with a BioMart table or an FTP `tsv/` table (87 on) in the window, and the
target. Rule `start_release`: without `from_release` an identifier is read at the target if the graph holds it there, else at
the latest release where it does.

Parts and strata (kinds = strata unless said):
- c1 NCBI gene ids: `c1_at` (every id at every table release where the graph holds it, read there); `c1_default` (held
  at the target); `c1_retired` (held before, not at the target: read at its last release, walked); `c1_discontinued`
  (NCBI's gene_history ids of the taxon that no release of the window has: nothing).
- c2 RefSeq (NM, NR, XM, XR, NP, XP): `c2_at`, `c2_at_v` (its version at R), `c2_default`, `c2_default_v`,
  `c2_retired`, `c2_oldver` (a version it had in the window, not at the target: read where it last had it);
  wrong versions (rule `accession_versions`, decided; sampled, 4,000 per class): `c2_ver_lower` (below every version
  of the window: nothing), `c2_ver_gap` (between two it had, never had: nothing), `c2_ver_higher` (above every version,
  held at the target: its gene), `c2_ver_higher_retired` (above every version of an accession not held at the target:
  nothing), `c2_ver_catalog` (the RefSeq catalog's version, newer than Ensembl's at the target: its gene).
- c3 UniProt: `c3_at`, `c3_default`, `c3_retired`, `c3_default_v`, `c3_ver_lower`, `c3_ver_higher` (as c2);
  `c3_isoform` (Ensembl's isoform accessions and Swiss-Prot's isoforms of the organism), `c3_entry` (entry names),
  `c3_secondary` (UniProt's secondary accessions of the organism's entries the graph never held), `c3_secondary_held`
  (those the graph held as a primary accession at some release), `c3_deleted` (accessions the graph held that UniProt
  deleted). Isoforms, entry names and secondary accessions are no database of the graph: rule `safe` (nothing or the
  gene is right, another gene wrong; the share recovered is reported).
- c4 RNA databases: `c4_rfam_acc`, `c4_rfam_label` (with `from_database` RFAM: the label is the entry; docs/methods.md:
  accession and label one value) at every table release and by default; miRBase precursor accession and name, mature
  accession and name, RNAcentral ids at the target (no database of the graph at d4c2dda: `safe`).
- c5 other gene databases: CCDS (accession and the versioned label), Havana gene, transcript, translation, Vega gene
  ids, clone-based names (at every table release where held, by default, retired); HGNC, MGI, VGNC ids (no database
  of the graph: `safe`).
- c6 genes to every database at every table release: every gene of the GTF at R x every final database of the graph
  held at R (identifier databases; the name databases are family B's B8): the values Ensembl gives the gene there
  (rule `one_kind_per_database`: `EntrezGene` the NCBI ids, `NCBI gene (formerly Entrezgene)` the NCBI symbols, at every release;
  RFAM accession and label, CCDS accession and versioned label: `spellings` of one value).
- c7 old externals to now: every value of c1 and c2 read at every table release R < target, to the target
  (exhaustive); UniProt, CCDS, Havana, Vega, clone-based, RFAM: a seeded sample of 50,000 (value, R) per organism.
"""
import argparse
import collections
import os
import random
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families.a import gen_a as GA
from idtrack_reproducibility.families.c import truth_c as TC

FAMILY = "c"
ME = os.path.abspath(__file__)
ME_SHA256 = L.sha256(ME)
SEED = 20261119
REFERENCE_TREE = "d4c2dda"  # its externals configuration is the graph's scope (truth_c.Scope)
N_VER = 4000        # sampled lines per wrong-version class
N_C7_OTHER = 50000  # c7: sampled (value, R) pairs per organism for the databases other than NCBI and RefSeq
REFSEQ_CLASS = re.compile(r"^(NM|NR|XM|XR|NP|XP|YP|WP|NG|NC)_")
# c6: the identifier databases of the graph (d4c2dda names) and how the truth reads them
C6_DBS = {
    "EntrezGene": ("acc", {"EntrezGene"}),
    "NCBI gene (formerly Entrezgene)": ("label", {"EntrezGene"}),
    "RefSeq_mRNA": ("acc", {"RefSeq_mRNA"}), "RefSeq_mRNA_predicted": ("acc", {"RefSeq_mRNA_predicted"}),
    "RefSeq_ncRNA": ("acc", {"RefSeq_ncRNA"}), "RefSeq_ncRNA_predicted": ("acc", {"RefSeq_ncRNA_predicted"}),
    "RefSeq_peptide": ("acc", {"RefSeq_peptide"}), "RefSeq_peptide_predicted": ("acc", {"RefSeq_peptide_predicted"}),
    "Uniprot/SWISSPROT": ("acc", {"Uniprot/SWISSPROT"}), "Uniprot/SPTREMBL": ("acc", {"Uniprot/SPTREMBL"}),
    "CCDS": ("acc+label", {"CCDS"}),
    "Havana gene": ("acc", {"OTTG"}), "Havana transcript": ("acc", {"OTTT"}), "Havana translation": ("acc", {"OTTP"}),
    "Vega_gene": ("acc", {"Vega_gene"}),
    "Clone_based_ensembl_gene": ("acc", {"Clone_based_ensembl_gene"}),
    "Clone_based_vega_gene": ("acc", {"Clone_based_vega_gene"}),
    "RFAM": ("acc+label", {"RFAM"}),
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rng(org, what):
    return random.Random(f"{SEED}:{org}:{what}")


def sample(org, what, pop, k):
    pop = sorted(pop)
    return pop if len(pop) <= k else sorted(rng(org, what).sample(pop, k))


class GenC:
    def __init__(self, org, T, families, authority=True):
        self.org, self.T = org, T
        self.G = GA.GenA(org, T)
        self.tr = self.G.tr
        self.counts = collections.Counter()
        asm = H.gtf_assemblies(org)
        self.scope = TC.Scope(org, CFG.get().tree(REFERENCE_TREE)["source"], asm)
        self.X = TC.Xrefs(org, self.scope, families=families, log=log) if families else None
        self.A = TC.Authority(org) if authority else None
        bm = set(H.biomart_releases(org))
        tsv = {int(d) for d in os.listdir(os.path.join(TC.TRUTH, org, "tsv"))
               if os.path.exists(os.path.join(TC.TRUTH, org, "tsv", d, "_done"))} \
            if os.path.isdir(os.path.join(TC.TRUTH, org, "tsv")) else set()
        xr = set(TC.xref_releases(org))
        self.table_rels = sorted(((bm | tsv) & set(self.tr.rels) & xr) | ({T} & xr))
        log("table releases", self.table_rels)

    # ------------------------------------------------------------ truth of one value

    def genes_ok(self, gs, R):
        """(genes present in the GTF of R, genes absent)."""
        at = self.tr.genes.get(R) or {}
        return {g for g in gs if g in at}, {g for g in gs if g not in at}

    def at_release(self, gs, R):
        """The truth of a value read at R, to R: its genes at R (copies on patches apart, rule `judge_copies_apart`)."""
        tr = self.tr
        ok, miss = self.genes_ok(gs, R)
        out = dict(cmp=tr.cmp_for(R))
        Ev = tr.vset(ok, R)
        Ev2, apart = self.G.apart_of(Ev, R)
        out["expect"] = sorted(Ev2)
        if apart:
            out["apart"] = apart
        if miss:
            out["open"] = f"Ensembl's xrefs name gene(s) absent from the GTF of {R}: {','.join(sorted(miss)[:3])}"
            self.counts["open: gene absent from the GTF"] += 1
        return out

    def walked(self, gs, R, t):
        """The truth of a value read at R, to t (`walked_truth`)."""
        return walked_truth(self.G, gs, R, t, self.counts)

    def default_truth(self, fam, v):
        """Rule `start_release`'s reading without from_release: (R read, truth to T) or (None, None) when the graph never holds v."""
        T = self.T
        h = self.X.at(fam, v, T)
        if h and h[2]:
            return T, self.at_release(h[0], T)
        R = self.X.last_in_graph(fam, v, T)
        if R is None:
            return None, None
        return R, self.walked(self.X.at(fam, v, R)[0], R, T)

    def line(self, kind, q, R, t, truth, unit, stratum=None, fdb=None, from_db=None, rule="set", **extra):
        qid = f"{kind}|{q}|{R if R is not None else '-'}|{t}" + (f"|{fdb}" if fdb else "") + \
              (f"|db={from_db}" if from_db else "")
        rec = dict(qid=qid, family=FAMILY, kind=kind, stratum=stratum or kind, unit=unit, id=q, from_release=R,
                   to_release=t, final_database=fdb, from_database=from_db, strategy="both", rule=rule)
        for k, x in truth.items():
            if x is not None:
                rec[k] = x
        if rule == "safe":
            for k in ("broad", "locus", "locus_mapped", "locus_projected", "exceptions", "apart"):
                rec.pop(k, None)
        rec.update({k: x for k, x in extra.items() if x is not None})
        self.counts[rec["stratum"]] += 1
        return rec

    def auth_genes(self, genes):
        """`auth`: the authority's genes (base ids, sorted) or ['none']."""
        return sorted(genes) if genes else ["none"]

    # ------------------------------------------------------------ the strata of an identifier database

    def held_strata(self, fam, prefix, auth=None, unit=None, vers=False, q_of=None):
        """at R, default, retired (and versions) of every value of `fam` the graph holds at some release."""
        T, X = self.T, self.X
        auth = auth or (lambda v: None)
        unit = unit or (lambda v: v)
        q_of = q_of or (lambda v: v)
        for v in sorted(X.H[fam]):
            hv = X.H[fam][v]
            ex = dict(auth=None, cls=self.cls(fam, v))
            for R in self.table_rels:
                h = hv.get(R)
                if not h:
                    continue
                if not h[2]:
                    self.counts[f"{prefix}: not held by the graph at a table release"] += 1
                    continue
                tru = self.at_release(h[0], R)
                yield self.line(f"{prefix}_at", q_of(v), R, R, tru, unit(v), **ex)
                if vers:
                    for ver in sorted((x for x in h[1] if x.isdigit()), key=int):
                        yield self.line(f"{prefix}_at_v", f"{q_of(v)}.{ver}", R, R, tru, unit(v), **ex)
            Rd, tru = self.default_truth(fam, v)
            if Rd is None:
                # Ensembl has the value, the graph's configuration never holds it (another object form, a database
                # or releases it leaves out): a stratum all the same (docs/methods.md): nothing, or its gene, is right
                Re = max(r for r in hv if r <= T) if any(r <= T for r in hv) else None
                if Re is None:
                    continue
                tru = self.at_release(hv[Re][0], T) if Re == T else self.walked(hv[Re][0], Re, T)
                yield self.line(f"{prefix}_not_held", q_of(v), None, T, tru, unit(v), rule="safe", read_at=Re,
                                ensembl_dbs=None, **ex)
                continue
            a = auth(v)
            ex["auth"] = self.auth_genes(a) if a is not None else None
            kind = f"{prefix}_default" if Rd == T else f"{prefix}_retired"
            yield self.line(kind, q_of(v), None, T, tru, unit(v), read_at=Rd, **ex)
            if vers and Rd == T:
                for ver in sorted((x for x in hv[T][1] if x.isdigit()), key=int):
                    yield self.line(f"{prefix}_default_v", f"{q_of(v)}.{ver}", None, T, tru, unit(v), read_at=Rd, **ex)

    def cls(self, fam, v):
        if fam == "refseq":
            m = REFSEQ_CLASS.match(v)
            return m.group(1) if m else "other"
        return None

    def version_strata(self, fam, prefix, catalog=None):
        """Rule `version_newer_than_graph`'s rule for accession versions (the decided answer to Q2, docs/methods.md): versions the graph's window had, and
        versions it never had (lower, gap, higher of a held or a retired accession, the database's newer version)."""
        T, X, org = self.T, self.X, self.org
        oldver, lower, gap, higher, higher_ret, cat = [], [], [], [], [], []
        for v in sorted(X.H[fam]):
            had = {int(x) for x in X.versions_in_graph(fam, v) if x.isdigit()}
            if not had:
                continue
            hT = X.at(fam, v, T)
            held_T = bool(hT and hT[2])
            atT = {int(x) for x in hT[1] if x.isdigit()} if held_T else set()
            for ver in sorted(had - atT):
                oldver.append((v, ver))
            lo, hi = min(had), max(had)
            if lo > 1:
                lower.append((v, lo - 1))
            gaps = sorted(set(range(lo, hi + 1)) - had)
            if gaps:
                gap.append((v, gaps[0]))
            (higher if held_T else higher_ret).append((v, hi + 1))
            if held_T:
                higher.append((v, hi + 5))
            if catalog is not None and held_T:
                newer = {x for x in catalog.get(v, ()) if x > hi}
                if newer:
                    cat.append((v, max(newer)))
        for v, ver in oldver:
            R = X.last_with_version(fam, v, str(ver), T)
            if R is None:
                continue
            yield self.line(f"{prefix}_oldver", f"{v}.{ver}", None, T, self.walked(X.at(fam, v, R)[0], R, T), v,
                            read_at=R, cls=self.cls(fam, v))
        nothing = dict(cmp=self.tr.cmp_for(T), expect=[])
        for kind, pop in ((f"{prefix}_ver_lower", lower), (f"{prefix}_ver_gap", gap),
                          (f"{prefix}_ver_higher_retired", higher_ret)):
            for v, ver in sample(org, kind, pop, N_VER):
                yield self.line(kind, f"{v}.{ver}", None, T, nothing, v, cls=self.cls(fam, v),
                                q2="a version the accession never had in the graph's window, not above every version "
                                   "of an accession the graph holds at its last release: nothing (rule `version_newer_than_graph`'s rule)")
        for kind, pop in ((f"{prefix}_ver_higher", higher), (f"{prefix}_ver_catalog", cat)):
            for v, ver in sample(org, kind, pop, N_VER):
                yield self.line(kind, f"{v}.{ver}", None, T, self.at_release(X.at(fam, v, T)[0], T), v,
                                cls=self.cls(fam, v), q2="above every version the accession had, held at the "
                                                         "graph's last release: its gene (rule `version_newer_than_graph`'s rule)")

    # ------------------------------------------------------------ parts

    def c1(self):
        A, X, T = self.A, self.X, self.T

        def auth(v):
            if v in A.history:
                rep = A.history[v][0]
                return A.ncbi.get(rep, set()) if rep else set()
            return A.ncbi.get(v, set())

        yield from self.held_strata("ncbi", "c1", auth=auth)
        seen = set(X.H["ncbi"])
        nothing = dict(cmp=self.tr.cmp_for(T), expect=[])
        for v, (rep, date, sym) in sorted(A.history.items()):
            if v in seen or not v.isdigit():
                continue
            yield self.line("c1_discontinued", v, None, T, nothing, v, discontinued=date, symbol=sym,
                            auth=self.auth_genes(A.ncbi.get(rep, set()) if rep else set()), replaced_by=rep)

    def c2(self):
        A = self.A

        def auth(v):
            out = set()
            for gid in A.refseq_gene.get(v, ()):
                out |= A.ncbi.get(gid, set())
            return out

        yield from self.held_strata("refseq", "c2", auth=auth, vers=True)
        yield from self.version_strata("refseq", "c2", catalog=A.catalog_ver)

    def c3(self):
        A, X, T, org = self.A, self.X, self.T, self.org
        yield from self.held_strata("uniprot", "c3", auth=lambda v: A.uniprot.get(v, set()), vers=True)
        for rec in self.version_strata("uniprot", "c3"):
            if rec["kind"] in ("c3_oldver", "c3_ver_lower", "c3_ver_higher", "c3_ver_higher_retired"):
                yield rec
        held_T = {v: h[T] for v, h in X.H["uniprot"].items() if T in h and h[T][2]}
        # isoforms: Ensembl's isoform xrefs at T and Swiss-Prot's isoforms of the organism's entries held at T
        iso = {}
        for v, h in X.H["isoform"].items():
            if T in h:
                iso[v] = ("ensembl", h[T][0])
        for v, parent in A.isoform_parent.items():
            if v not in iso and parent in held_T:
                iso[v] = ("swissprot", held_T[parent][0])
        for v, (src, gs) in sorted(iso.items()):
            if v.endswith("-1") and src == "ensembl":
                kind = "c3_isoform_canonical"
            else:
                kind = "c3_isoform"
            yield self.line(kind, v, None, T, self.at_release(gs, T), v.split("-")[0], rule="safe", src=src,
                            auth=self.auth_genes(A.uniprot.get(v.split("-")[0], set())))
        for v, h in sorted(held_T.items()):
            e = A.entry.get(v)
            if e:
                yield self.line("c3_entry", e, None, T, self.at_release(h[0], T), v, rule="safe", accession=v,
                                auth=self.auth_genes(A.uniprot.get(v, set())))
        for s, prims in sorted(A.secondary.items()):
            if s in held_T:
                continue  # also a primary accession held at T: c3_default asks it
            gs = set()
            for p in prims:
                if p in held_T:
                    gs |= held_T[p][0]
            if s in X.H["uniprot"]:
                Rd, tru = self.default_truth("uniprot", s)
                if Rd is not None:
                    yield self.line("c3_secondary_held", s, None, T, tru, s, read_at=Rd, primaries=sorted(prims),
                                    auth=self.auth_genes(gs))
                    continue
            if not gs:
                continue
            yield self.line("c3_secondary", s, None, T, self.at_release(gs, T), s, rule="safe", primaries=sorted(prims),
                            auth=self.auth_genes(gs))
        held_ever = {v for v, h in X.H["uniprot"].items() if any(x[2] for x in h.values())}
        dele = A.deleted_uniprot(held_ever)
        for v in sorted(dele):
            Rd, tru = self.default_truth("uniprot", v)
            if Rd is None:
                continue
            yield self.line("c3_deleted", v, None, T, tru, v, read_at=Rd, auth=["none"])

    def c4(self):
        A, X, T = self.A, self.X, self.T
        for fam, prefix, from_db in (("rfam", "c4_rfam_acc", None), ("rfam_label", "c4_rfam_label", "RFAM")):
            for rec in self.held_strata(fam, prefix):
                if from_db:
                    rec["from_database"] = from_db
                    rec["qid"] += f"|db={from_db}"
                    lab = rec["id"]
                    accs = {a for a, labs in X.labels["rfam"].items() if lab in labs}
                    rec["rfam_acc"] = sorted(accs)
                    rec["auth_label_of"] = sorted(a for a in accs if A.rfam_label.get(a) == lab) or ["none"]
                else:
                    rec["auth_label"] = A.rfam_label.get(rec["id"])
                yield rec
        # miRBase and RNAcentral: no database of the graph (d4c2dda): `safe`, at the target
        precursor_genes = {}
        for fam, kind in (("mirbase", "c4_mirbase_acc"), ("mirbase_label", "c4_mirbase_label")):
            for v, h in sorted(X.H[fam].items()):
                if T not in h:
                    continue
                gs = h[T][0]
                if fam == "mirbase":
                    precursor_genes[v] = gs
                auth = set()
                acc = v if fam == "mirbase" else A.mir_acc.get(v)
                for urs in A.mirbase_urs.get(acc or "", ()):
                    auth |= A.rnacentral.get(urs, set())
                yield self.line(kind, v, None, T, self.at_release(gs, T), acc or v, rule="safe",
                                auth=self.auth_genes(auth))
        for mat, precs in sorted(A.mature_of.items()):
            gs = set()
            for p in precs:
                gs |= precursor_genes.get(p, set())
            if not gs:
                continue
            for q, kind in ((mat, "c4_mature_acc"), (A.mir_name.get(mat), "c4_mature_name")):
                if q:
                    yield self.line(kind, q, None, T, self.at_release(gs, T), mat, rule="safe", precursors=sorted(precs))
        for v, h in sorted(X.H["rnacentral"].items()):
            if T in h:
                yield self.line("c4_rnacentral", v, None, T, self.at_release(h[T][0], T), v, rule="safe",
                                auth=self.auth_genes(A.rnacentral.get(v, set())))

    def c5(self):
        A, X, T = self.A, self.X, self.T

        def ccds_auth(v):
            gs = set()
            tx = self.tx_gene(T)
            for t in A.ccds_tx.get(v, ()) | A.ccds_tx.get(v.split(".")[0], ()):
                if t in tx:
                    gs.add(tx[t])
            return gs

        yield from self.held_strata("ccds", "c5_ccds", auth=ccds_auth)
        # the versioned label of every CCDS accession (Ensembl writes `CCDS11118` and `CCDS11118.1`)
        for v, h in sorted(X.H["ccds"].items()):
            labs = sorted(X.labels["ccds"].get(v, ()))
            Rd, tru = self.default_truth("ccds", v)
            if Rd is None:
                continue
            for lab in labs:
                yield self.line("c5_ccds_label" if Rd == T else "c5_ccds_label_retired", lab, None, T, tru, v,
                                read_at=Rd, auth=self.auth_genes(ccds_auth(lab)))
        for fam in ("ottg", "ottt", "ottp", "vega_gene", "clone"):
            yield from self.held_strata(fam, f"c5_{fam}")
        auth_of = dict(hgnc_id=lambda v: None, mgi_id=lambda v: None, vgnc_id=lambda v: A.vgnc.get(v, set()))
        for fam in ("hgnc_id", "mgi_id", "vgnc_id"):
            for v, h in sorted(X.H[fam].items()):
                if T in h:
                    a = auth_of[fam](v)
                    yield self.line(f"c5_{fam}", v, None, T, self.at_release(h[T][0], T), v, rule="safe",
                                    auth=self.auth_genes(a) if a is not None else None)

    def tx_gene(self, R):
        if not hasattr(self, "_tx") or self._tx[0] != R:
            m = {}
            for r in H.gtf_tx(self.org, R) or []:
                if r.get("transcript_id"):
                    m[r["transcript_id"]] = r["gene_id"]
            self._tx = (R, m)
        return self._tx[1]

    def c6(self):
        """Every gene of the GTF at a table release x every identifier database the graph holds there."""
        tr, org = self.tr, self.org
        for R in self.table_rels:
            # rule `one_kind_per_database`: the NCBI symbols carry one name at every release at which the graph holds NCBI ids
            dbs = [n for n in C6_DBS if self.scope.included(n, R) or
                   (n == "NCBI gene (formerly Entrezgene)" and self.scope.included("EntrezGene", R))]
            if not dbs:
                continue
            want = set().union(*(C6_DBS[n][1] for n in dbs))
            vals = {n: collections.defaultdict(set) for n in dbs}
            spell = {n: collections.defaultdict(set) for n in dbs}
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, R):
                if db not in want:
                    continue
                for n in dbs:
                    how, src = C6_DBS[n]
                    if db not in src:
                        continue
                    if how == "acc" and acc:
                        vals[n][gid].add(acc)
                    elif how == "label" and lab:
                        vals[n][gid].add(lab)
                    elif how == "acc+label" and acc:
                        vals[n][gid].add(acc)
                        if lab and lab != acc:
                            vals[n][gid].add(lab)
                            spell[n][gid].add((acc, lab))
            genes = sorted(tr.genes.get(R) or {})
            for g in genes:
                q = tr.vid(g, R)
                for n in dbs:
                    E = sorted(vals[n].get(g, ()))
                    rec = self.line("c6", q, R, R, dict(cmp="exact", expect=E), g, stratum=f"c6:{n}", fdb=n)
                    if spell[n].get(g):
                        rec["spellings"] = sorted([a, b] for a, b in spell[n][g])
                    yield rec
            log("c6", R, len(genes), "genes", len(dbs), "databases", dict(self.counts))

    def c7(self):
        X, T, org = self.X, self.T, self.org
        for fam, prefix in (("ncbi", "c7_ncbi"), ("refseq", "c7_refseq")):
            for v in sorted(X.H[fam]):
                for R, h in sorted(X.H[fam][v].items()):
                    if R in self.table_rels and R < T and h[2]:
                        yield self.line(prefix, v, R, T, self.walked(h[0], R, T), v, cls=self.cls(fam, v))
        pairs = []
        for fam in ("uniprot", "ccds", "ottg", "ottt", "ottp", "vega_gene", "clone", "rfam", "rfam_label"):
            if fam not in X.H:
                continue
            for v, hv in X.H[fam].items():
                for R, h in hv.items():
                    if R in self.table_rels and R < T and h[2]:
                        pairs.append((fam, v, R))
        for fam, v, R in sample(org, "c7_other", pairs, N_C7_OTHER):
            rec = self.line(f"c7_{fam}", v, R, T, self.walked(X.H[fam][v][R][0], R, T), v)
            if fam == "rfam_label":
                rec["from_database"] = "RFAM"
                rec["qid"] += "|db=RFAM"
            yield rec


PARTS = {
    "c1": ["ncbi"], "c2": ["refseq"], "c3": ["uniprot", "isoform"], "c4": ["rfam", "rfam_label", "mirbase",
                                                                         "mirbase_label", "rnacentral"],
    "c5": ["ccds", "ottg", "ottt", "ottp", "vega_gene", "clone", "hgnc_id", "mgi_id", "vgnc_id"], "c6": [],
    "c7": ["ncbi", "refseq", "uniprot", "ccds", "ottg", "ottt", "ottp", "vega_gene", "clone", "rfam", "rfam_label"],
}


# Fields of a line that explain it in words (the rule it follows, how it was made): a line whose other fields are
# unchanged is the same truth, so `--revise-from` does not revise it for the wording of its explanation.
WORDING = ("note", "q2")


def walked_truth(G, gs, R, t, counts=None):
    """The truth of a value read at R, to t: the union of its genes' walks (one gene: A's truth_at whole); G is family
    A's generator (`GenA`), `counts` a Counter of the open cases."""
    tr = G.tr
    counts = counts if counts is not None else collections.Counter()
    at = tr.genes.get(R) or {}
    ok, miss = {g for g in gs if g in at}, {g for g in gs if g not in at}
    if miss or not ok:
        counts["open: gene absent from the GTF"] += 1
        return dict(cmp=tr.cmp_for(t), expect=[], open=f"Ensembl's xrefs name gene(s) absent from the GTF of {R}: "
                                                        f"{','.join(sorted(miss)[:3])}")
    if len(ok) == 1:
        g = next(iter(ok))
        E, B = G.walks(g, R, [t])
        return G.truth_at(g, R, t, E, B)
    Es, Bs, empty = set(), set(), 0
    for g in sorted(ok):
        E, B = G.walks(g, R, [t])
        e = E.get(t, set())
        Es |= e
        Bs |= B.get(t, set()) | e
        empty += not e
    out = dict(cmp=tr.cmp_for(t))
    if not Es:
        out["expect"] = []
        out["open"] = "several genes carry the value and Ensembl's history gives none of them a successor"
        counts["open: several genes without successor"] += 1
        return out
    Ev = tr.vset(Es, t)
    Ev2, apart = G.apart_of(Ev, t)
    out["expect"] = sorted(Ev2)
    if apart:
        out["apart"] = apart
    if Bs != Es:
        out["broad"] = tr.vset(Bs, t)
    if empty:
        out["partial"] = empty
        counts["several genes, some without successor"] += 1
    return out


def truth_hash(row):
    """A hash of a typed row (lines.to_row) without its explanations in words (`WORDING`) and its shard hash."""
    import json
    r = {k: v for k, v in row.items() if k not in WORDING and k != "h"}
    if r.get("extra"):
        r["extra"] = C.dumps({k: v for k, v in json.loads(r["extra"]).items() if k not in WORDING})
    return hash(C.dumps(r))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--part", required=True, choices=sorted(PARTS))
    ap.add_argument("--set", required=True)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--limit", type=int, default=None, help="stop after N lines (a test; never for a published set)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--show", type=int, default=0, help="with --dry-run: print the first N lines of every stratum")
    ap.add_argument("--revise-from", default=None,
                    help="a published set of this part: publish the lines whose truth changed (the explanations in "
                         "words aside, `WORDING`) as a revision set "
                         "`<set>` (revises it; base lines no longer generated are left out, `excluded`) and the new "
                         "qids as the set `<set>_new` (truth that gained releases after the base set was published)")
    a = ap.parse_args()
    org = a.organism
    log("truth", org, a.part)
    Gc = GenC(org, a.target, PARTS[a.part], authority=a.part in ("c1", "c2", "c3", "c4", "c5"))
    gen = getattr(Gc, a.part)

    shown = collections.Counter()

    def all_lines():
        k = 0
        for rec in gen():
            k += 1
            if a.show and shown[rec["stratum"]] < a.show:
                shown[rec["stratum"]] += 1
                print("LINE", C.dumps(rec), flush=True)
            if a.limit and k > a.limit:
                break
            if k % 1_000_000 == 0:
                log("lines", k)
            yield rec

    sources = ["core_xref/<release>.tsv.gz (every release of the window; fetch_c.py ensembl)",
               "gtf/genes_<release>.tsv.gz", f"core/{Gc.tr.core}.stable_id_event + mapping_session",
               "chain/*.chain.gz", "exceptions/<release>.tsv"] + sorted(set(Gc.A.files) if Gc.A else [])
    extra = dict(part=a.part, generator_sha256_at_start=ME_SHA256, truth_c_sha256=L.sha256(TC.__file__),
                 gen_a_sha256=L.sha256(GA.__file__), scope_yaml=Gc.scope.yaml, table_releases=Gc.table_rels,
                 target=a.target, limit=a.limit)
    if a.revise_from:
        import pyarrow.parquet as pq
        base = {}
        for b in pq.ParquetFile(L._set_file(L.query_dir(FAMILY, org), a.revise_from)).iter_batches(batch_size=500_000):
            for r in b.to_pylist():
                base[r["qid"]] = truth_hash(r)
        changed, new_lines, seen = [], [], set()
        for rec in all_lines():
            q = rec["qid"]
            seen.add(q)
            if q not in base:
                new_lines.append(rec)
            elif base[q] != truth_hash(L.to_row(rec)):
                rec["revised_by"] = list(rec.get("revised_by") or []) + ["releases_added"]
                changed.append(rec)
        gone = [q for q in base if q not in seen]
        if gone:
            old = {r["qid"]: r for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), a.revise_from))
                   if r["qid"] in set(gone)}
            for q in gone:
                r = dict(old[q])
                r["excluded"] = "no longer generated once every release of the window was read"
                r["revised_by"] = list(r.get("revised_by") or []) + ["releases_added"]
                changed.append(r)
        log("revise-from", a.revise_from, "changed", len(changed), "new", len(new_lines), "gone", len(gone))
        at = lambda r: r.get("from_release") or r.get("read_at") or r.get("to_release")  # noqa: E731
        log("by the release read: changed", dict(sorted(collections.Counter(at(r) for r in changed).items())),
            "new", dict(sorted(collections.Counter(at(r) for r in new_lines).items())))
        man = None
        if changed:
            man = L.publish(FAMILY, org, a.set, changed,
                            f"C {a.part}: revision of {a.revise_from} (truth read with every release of the window; "
                            f"gen_c.py --revise-from)", ME, revises=[a.revise_from], judge="c-releases_added",
                            sources=sources, manifest_extra=dict(extra, correction="releases_added"), dry=a.dry_run, log=log)
        if new_lines:
            man = L.publish(FAMILY, org, a.set + "_new", new_lines,
                            f"C {a.part}: lines new beside {a.revise_from} (releases added; gen_c.py --revise-from)", ME,
                            sources=sources, manifest_extra=extra, dry=a.dry_run, log=log)
    else:
        man = L.publish(FAMILY, org, a.set, all_lines(),
                        f"C {a.part} (docs/methods.md, C; gen_c.py docstring); strategy both; truth: Ensembl's core xrefs "
                        f"of every release, A's walks; the authorities as `auth`",
                        ME, sources=sources, manifest_extra=extra, dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k not in ("counts",)})
    log("counts", (man or {}).get("counts"))
    log("generation counts", dict(Gc.counts))


if __name__ == "__main__":
    main()
