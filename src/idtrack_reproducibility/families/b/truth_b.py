"""Family B's truth: what a string names at every release, from Ensembl's own tables, and what the authorities say.

Nothing here reads idtrack. Sources (the one truth store `accuracy/truth/<organism>/`, every file in SOURCES.tsv):

- **display names**: the GTF's `gene_name` of every release (Ensembl's display name; `gtf/genes_<R>.tsv.gz`), with
  the region (primary, patch) where the release's GTF marks patches (human and mouse from 82);
- **other name sources at R** (the corrected judge `names_broad`): where an archive BioMart answers, its HGNC, MGI, UniProt
  gene name, NCBI label and WikiGene tables, VGNC's symbol of BioMart's VGNC id and the Vega names (the harness's
  `Ctx78.name_genes`); and at every release, Ensembl's core xrefs of the release (`core_xref/<R>.tsv.gz`, family C's
  extract of the MySQL dumps): the labels of HGNC, MGI, VGNC, UniProt gene name, NCBI gene, WikiGene, Vega and the
  clone-based names (RFAM and miRBase labels are no names: rules `mirbase_genes` and `rfam_display_label`);
- **synonyms at R**: BioMart's synonym table (only at BioMart releases; the core extract has none);
- **the reach** (`ever`): every gene that carries the string in any of these tables at any release read. At a release
  whose other name sources are not read, an answer within the reach is open, an answer outside it is wrong;
- **the history**: the harness's walks (`astar.truth.TruthA`: `make_queries.Truth` with sessions placed where they
  take effect, the core genes the GTF lacks), strict with living ids, broad with ids brought back (rule `resurrection`);
- **the authorities** (principle 2, measured beside idtrack, never a wrong answer): HGNC's complete set (symbol,
  prev_symbol, alias_symbol, ensembl_gene_id) and withdrawn list; MGI's MRK_List1 (withdrawn symbols and their
  current symbol, synonyms) and MRK_ENSEMBL (MGI's own Ensembl ids); VGNC's gene sets (symbol, prev_symbol,
  alias_symbol, ensembl_gene_id); NCBI's gene_info (Synonyms, dbXrefs Ensembl) where fetched.
"""
import collections
import csv
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C  # noqa: F401  puts the harness on the path
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M
from idtrack_reproducibility.families.astar.truth import TruthA

csv.field_size_limit(1 << 30)
T = 115
# name databases of the core xref extract (by `db_name`): their labels are names of the gene
NAME_DB = re.compile(r"^(HGNC|HGNC_.*|MGI|MGI_.*|VGNC|Uniprot_gn|EntrezGene|WikiGene|Vega_gene|Clone_based_.*)$")
ENS_LIKE = re.compile(r"^ENS[A-Z]*[GTPE]\d{11}")


def clean(v):
    """A string the harness asks as a name: not an Ensembl id, not a number, no white space at its ends, not empty."""
    return bool(v) and v.strip() == v and not ENS_LIKE.match(v) and not v.isdigit()


def core_synonyms(org, R):
    """synonym -> {gene ids} of Ensembl's synonyms of every xref of a name database at R (`core_synonyms/<R>.tsv.gz`,
    fetch_b.py synonyms), or None where not fetched."""
    p = os.path.join(H.truth_dir(org), "core_synonyms", f"{R}.tsv.gz")
    if not os.path.exists(p):
        return None
    out = collections.defaultdict(set)
    with gzip.open(p, "rt") as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and NAME_DB.match(c[2]) and c[0]:
                out[c[0]].add(C.base_id(c[1]))
    return out


class NamesB:
    def __init__(self, org, log=print):
        self.org, self.log = org, log
        self.tr = TruthA(org, T)
        self.C = M.Ctx78(org, T, truth=self.tr)
        self.rels = [r for r in H.gtf_releases(org) if r <= T and r in self.tr.genes]
        self.bm_rels = {r for r in H.biomart_releases(org) if r <= T}
        self.xref_rels = {int(f[:-7]) for f in os.listdir(os.path.join(H.truth_dir(org), "core_xref"))
                          if f.endswith(".tsv.gz") and f[:-7].isdigit()} if os.path.isdir(
            os.path.join(H.truth_dir(org), "core_xref")) else set()
        self._disp, self._src, self._syn, self._walk = {}, {}, {}, {}
        self._ever = None

    # ---------------- per release ----------------

    def genes(self, R):
        return self.tr.genes[R]

    def gset(self, R):
        """The set of gene ids of release R (kept: built once per release)."""
        c = self.__dict__.setdefault("_gset", {})
        if R not in c:
            c[R] = frozenset(self.tr.genes[R])
        return c[R]

    def marks_patches(self, R):
        c = self.__dict__.setdefault("_marks", {})
        if R not in c:
            c[R] = any(r["region"] == "patch" for r in self.tr.genes[R].values())
        return c[R]

    def _slot(self, cache, R, make):
        """Per-release tables: the target release's kept, of the others one at a time (releases are walked in order)."""
        if R not in cache:
            v = make()
            for k in [k for k in cache if k != T]:
                del cache[k]
            cache[R] = v
        return cache[R]

    def display(self, R):
        """name -> {gene ids} of the GTF of release R (every region)."""
        def make():
            d = collections.defaultdict(set)
            for g, r in self.tr.genes[R].items():
                if r.get("gene_name"):
                    d[r["gene_name"]].add(g)
            return d
        return self._slot(self._disp, R, make)

    def core_names(self, R):
        """(name -> {gene ids}, name -> {gene ids} without HGNC's databases, gene -> {db_name: {labels}}) from the
        core xref extract of R (name databases only), or None."""
        if R not in self.xref_rels:
            return None
        p = os.path.join(H.truth_dir(self.org), "core_xref", f"{R}.tsv.gz")
        out, non_hgnc, by_gene = collections.defaultdict(set), collections.defaultdict(set), collections.defaultdict(dict)
        with gzip.open(p, "rt") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 7 or not NAME_DB.match(c[3]):
                    continue
                lab, gid = c[6], C.base_id(c[2])
                if lab and gid and not lab.isdigit():
                    out[lab].add(gid)
                    if not c[3].startswith("HGNC"):
                        non_hgnc[lab].add(gid)
                    by_gene[gid].setdefault(c[3], set()).add(lab)
        return out, non_hgnc, by_gene

    def sources(self, R):
        """(name -> {genes} in the other name sources at R, what was read: `biomart`, `core_xref`, both, or None,
        name -> {genes} in the name sources other than HGNC's, gene -> {database: {values}} of the core extract)."""
        def make():
            out, non_hgnc, read, by_gene = collections.defaultdict(set), collections.defaultdict(set), [], {}
            if R in self.bm_rels:
                read.append("biomart")
                for grp, col in M.NAME_SOURCES:
                    pp = self.C.bm(R, grp, col)
                    if pp:
                        for v, gs in pp[1].items():
                            out[v] |= gs
                            if grp != "hgnc":
                                non_hgnc[v] |= gs
                for pp in (self.C.vgnc(R), self.C.vega(R)):
                    if pp:
                        for v, gs in pp[1].items():
                            out[v] |= gs
                            non_hgnc[v] |= gs
            cx = self.core_names(R)
            if cx is not None:
                read.append("core_xref")
                for v, gs in cx[0].items():
                    out[v] |= gs
                for v, gs in cx[1].items():
                    non_hgnc[v] |= gs
                by_gene = cx[2]
            return out, "+".join(read) or None, non_hgnc, by_gene
        return self._slot(self._src, R, make)

    def core_synonyms(self, R):
        return core_synonyms(self.org, R)

    def synonyms(self, R):
        """name -> {genes} as a synonym at R: BioMart's synonym table (the synonyms of a gene's display xref) and,
        where fetched, Ensembl's synonyms of every xref of a name database (core); None where neither is read."""
        def make():
            pp = self.C.bm(R, "synonym", 1) if R in self.bm_rels else None
            cs = self.core_synonyms(R)
            if pp is None and cs is None:
                return None
            out = collections.defaultdict(set)
            for d in ((pp[1] if pp else {}), (cs or {})):
                for v, gs in d.items():
                    out[v] |= gs
            return out
        return self._slot(self._syn, R, make)

    SUPPORT_GROUPS = (("hgnc", re.compile(r"^HGNC")), ("mgi", re.compile(r"^MGI")), ("uniprot_gn", re.compile(r"^Uniprot_gn$")),
                      ("entrez", re.compile(r"^EntrezGene$")), ("wikigene", re.compile(r"^WikiGene$")))

    def support(self, R, gene, v):
        """Rule `support_ranks_shared_name`'s support: in how many name sources the gene carries v at R. BioMart's (display name, HGNC,
        MGI, UniProt gene name, NCBI label, WikiGene, synonym) where the archive answers; elsewhere the GTF's display
        name and the same databases in the core xrefs of R, without the synonym (not read there): a support that can
        only make an answer undetermined that a full count might call wrong, never the reverse."""
        if R in self.bm_rels:
            return self.C.support(R, gene, v)
        by_gene = self.sources(R)[3]
        if not by_gene:
            return None
        n = int(self.tr.genes[R].get(gene, {}).get("gene_name") == v)
        dbs = by_gene.get(gene, {})
        for _, rx in self.SUPPORT_GROUPS:
            if any(rx.match(d) and v in labs for d, labs in dbs.items()):
                n += 1
        return n

    # ---------------- every release ----------------

    def ever(self, kind="all"):
        """string -> {gene base ids} that carry it in any table read at any release: `names` (display names and the
        other name sources), `synonyms` (BioMart's synonym tables) or `all`."""
        if self._ever is None:
            names, syns = collections.defaultdict(set), collections.defaultdict(set)
            for R in self.rels:
                for v, gs in self.display(R).items():
                    names[v] |= gs
                for v, gs in self.sources(R)[0].items():
                    names[v] |= gs
                for v, gs in (self.synonyms(R) or {}).items():
                    syns[v] |= gs
            alls = collections.defaultdict(set)
            for d in (names, syns):
                for v, gs in d.items():
                    alls[v] |= gs
            self._ever = dict(names=names, synonyms=syns, all=alls)
            self.log("ever", len(alls), "strings", len(names), "names", len(syns), "synonyms")
        return self._ever[kind]

    # ---------------- history ----------------

    def walk(self, g, R, b=T):
        """(strict walk with living ids, broad walk with living ids and returns) of gene g from R to b."""
        k = (g, R, b)
        w = self._walk.get(k)
        if w is None:
            s = self.tr.walk_fwd(g, R, b, True, living=True) if b > R else (
                self.tr.walk_back(g, R, b, True, living=True) if b < R else ({g} if g in self.tr.genes[R] else set()))
            br = self.tr.walk_fwd(g, R, b, False, living=True, resurrect=True) if b > R else (
                self.tr.walk_back(g, R, b, False, living=True, resurrect=True) if b < R else set(s))
            w = (frozenset(s), frozenset(br | s))
            self._walk[k] = w
        return w

    def walks(self, genes, R, b=T):
        s, br = set(), set()
        for g in genes:
            a, c = self.walk(g, R, b)
            s |= a
            br |= c
        return s, br

    def vset(self, genes, R):
        return sorted(self.tr.vid(g, R) for g in genes)

    def region(self, g, R):
        r = self.tr.genes[R].get(g)
        return r["region"] if r else None


# ---------------- the authorities ----------------

def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def hgnc(org):
    """HGNC's complete set: rows with symbol, prev_symbol, alias_symbol (lists), ensembl_gene_id, hgnc_id."""
    if org != "homo_sapiens":
        return None
    p = os.path.join(H.truth_dir(org), "hgnc", "hgnc_complete_set.txt")
    if not os.path.exists(p):
        return None
    out = []
    for r in read_tsv(p):
        out.append(dict(id=r["hgnc_id"], symbol=r["symbol"], prev=[x for x in r["prev_symbol"].strip('"').split("|") if x],
                        alias=[x for x in r["alias_symbol"].strip('"').split("|") if x], ensembl=r["ensembl_gene_id"] or None,
                        status=r["status"]))
    return out


def mgi(org):
    """MGI: current markers (symbol, synonyms, Ensembl ids from MRK_ENSEMBL) and withdrawn symbols with the current
    marker they lead to."""
    if org != "mus_musculus":
        return None
    d = os.path.join(H.truth_dir(org), "mgi")
    if not os.path.isdir(d):
        return None
    ens = collections.defaultdict(set)
    with open(os.path.join(d, "MRK_ENSEMBL.rpt")) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > 5 and c[5].startswith("ENSMUSG"):
                ens[c[0]].add(c[5])
    cur, withdrawn = [], []
    with open(os.path.join(d, "MRK_List1.rpt"), encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        for r in rd:
            if r["Status"] == "O":
                cur.append(dict(id=r["MGI Accession ID"], symbol=r["Marker Symbol"], ensembl=sorted(ens.get(r["MGI Accession ID"], ())),
                                alias=[x for x in (r["Marker Synonyms (pipe-separated)"] or "").split("|") if x],
                                type=r["Marker Type"]))
            elif r["Status"] == "W":
                to = r["Current MGI Accession ID (if withdrawn)"]
                withdrawn.append(dict(symbol=r["Marker Symbol"], to_id=to, to_symbol=r["Current Marker Symbol (if withdrawn)"],
                                      ensembl=sorted(ens.get(to, ())), name=r["Marker Name"]))
    return dict(current=cur, withdrawn=withdrawn)


def vgnc(org):
    d = os.path.join(H.truth_dir(org), "vgnc")
    files = [f for f in os.listdir(d) if f.endswith("_vgnc_gene_set_All.tsv")] if os.path.isdir(d) else []
    if not files:
        return None
    out = []
    for r in read_tsv(os.path.join(d, files[0])):
        out.append(dict(id=r["vgnc_id"], symbol=r["symbol"], prev=[x for x in (r.get("prev_symbol") or "").split("|") if x],
                        alias=[x for x in (r.get("alias_symbol") or "").split("|") if x],
                        ensembl=(r.get("ensembl_gene_id") or None), status=r.get("status")))
    return out


def ncbi(org):
    """NCBI gene_info of the taxon: (gene id, symbol, synonyms, Ensembl gene ids of dbXrefs), or None."""
    d = os.path.join(H.truth_dir(org), "ncbi")
    files = [f for f in os.listdir(d) if f.endswith(".gene_info.gz")] if os.path.isdir(d) else []
    if not files:
        return None
    out = []
    with gzip.open(os.path.join(d, files[0]), "rt", encoding="utf-8", errors="replace") as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            ens = sorted(x.split(":", 1)[1] for x in r.get("dbXrefs", "").split("|") if x.startswith("Ensembl:"))
            out.append(dict(id=r["GeneID"], symbol=r["Symbol"], alias=[x for x in r.get("Synonyms", "-").split("|") if x and x != "-"],
                            ensembl=ens))
    return out
