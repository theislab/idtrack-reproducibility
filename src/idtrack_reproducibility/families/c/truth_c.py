"""Truth of family C (external identifiers): Ensembl's xrefs at every release and each database's own files.

Nothing here reads idtrack's code or answers. The scope of the graph (which databases at which releases, the window)
is read from the externals configuration of the reference tree (`default_config/<organism>_externals_modified.yml`),
as the window is: it says what a graph holds, not what is right.

- **Ensembl's truth** (`Xrefs`): the core database's xrefs of every release of the window (`core_xref/<release>.tsv.gz`,
  written by `fetch_c.py ensembl` from the FTP MySQL dumps: one row per Ensembl object and xref of the databases of
  family C; the object's gene is the gene a value names). A value of a database family (`FAMILIES`) at release R names
  the genes whose gene, transcript or translation carries it at R.
- **The authorities** (`Authority`): NCBI `gene2ensembl`, `gene2refseq`, `gene_history`; UniProt `idmapping_selected`,
  `sec_ac`, `delac_sp`, `delac_tr`, Swiss-Prot isoforms; RFAM's family list; RNAcentral's Ensembl and miRBase
  mappings; miRBase's genome files; CCDS's sequence table; VGNC's gene sets. Each gives the gene(s) the database itself
  attaches to its identifier now (the file's date).
- Cross-checks (`tsv_pairs`, BioMart): the FTP `tsv/` tables (87 on) and BioMart, read for the agreement with the
  core reading.
"""
import collections
import glob
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C

TRUTH = C.TRUTH
SHARED = os.path.join(TRUTH, "_shared")
TAXON = {"homo_sapiens": "9606", "mus_musculus": "10090", "sus_scrofa": "9823", "macaca_mulatta": "9544",
         "macaca_fascicularis": "9541"}
GENE_PREFIX = {"homo_sapiens": "ENSG", "mus_musculus": "ENSMUSG", "sus_scrofa": "ENSSSCG", "macaca_mulatta": "ENSMMUG",
               "macaca_fascicularis": "ENSMFAG"}
# The YAML's assembly keys (`Assembly: "38"`) of the GTF's assembly names.
ASM_KEY = {"GRCh38": "38", "GRCh37": "37", "NCBI36": "36", "GRCm39": "39", "GRCm38": "38", "NCBIM37": "37",
           "Sscrofa11.1": "111", "Sscrofa10.2": "102", "Sscrofa9": "9", "Mmul_10": "10", "Mmul_8.0.1": "801",
           "MMUL_1": "1", "Macaca_fascicularis_6.0": "6", "Macaca_fascicularis_5.0": "5"}

# Database families: Ensembl `db_name`s whose values are one kind of identifier, the field read (acc, label or both),
# and the graph's database names (d4c2dda, the reference) that hold them.
FAMILIES = {
    "ncbi": dict(dbs={"EntrezGene"}, field="acc", graph=["EntrezGene"]),
    "ncbi_label": dict(dbs={"EntrezGene"}, field="label", graph=["NCBI gene (formerly Entrezgene)", "NCBI gene"]),
    "refseq": dict(dbs={"RefSeq_mRNA", "RefSeq_mRNA_predicted", "RefSeq_ncRNA", "RefSeq_ncRNA_predicted",
                        "RefSeq_peptide", "RefSeq_peptide_predicted"}, field="acc", versioned=True,
                   graph=["RefSeq_mRNA", "RefSeq_mRNA_predicted", "RefSeq_ncRNA", "RefSeq_ncRNA_predicted",
                          "RefSeq_peptide", "RefSeq_peptide_predicted"]),
    "uniprot": dict(dbs={"Uniprot/SWISSPROT", "Uniprot/SPTREMBL"}, field="acc", versioned=True,
                    graph=["Uniprot/SWISSPROT", "Uniprot/SPTREMBL"]),
    "isoform": dict(dbs={"Uniprot_isoform"}, field="acc", graph=[]),
    "ccds": dict(dbs={"CCDS"}, field="acc", graph=["CCDS"]),
    "ottg": dict(dbs={"OTTG"}, field="acc", graph=["Havana gene"]),
    "ottt": dict(dbs={"OTTT"}, field="acc", graph=["Havana transcript"]),
    "ottp": dict(dbs={"OTTP"}, field="acc", graph=["Havana translation"]),
    "vega_gene": dict(dbs={"Vega_gene"}, field="acc", graph=["Vega_gene"]),
    "clone": dict(dbs={"Clone_based_ensembl_gene", "Clone_based_vega_gene"}, field="acc",
                  graph=["Clone_based_ensembl_gene", "Clone_based_vega_gene"]),
    "rfam": dict(dbs={"RFAM"}, field="acc", graph=["RFAM"]),
    "rfam_label": dict(dbs={"RFAM"}, field="label", graph=["RFAM"]),
    "mirbase": dict(dbs={"miRBase"}, field="acc", graph=[]),
    "mirbase_label": dict(dbs={"miRBase"}, field="label", graph=[]),
    "rnacentral": dict(dbs={"RNAcentral", "RNACentral"}, field="acc", graph=[]),
    "hgnc_id": dict(dbs={"HGNC"}, field="acc", graph=[]),
    "mgi_id": dict(dbs={"MGI"}, field="acc", graph=[]),
    "vgnc_id": dict(dbs={"VGNC"}, field="acc", graph=[]),
}
DB_FAMILIES = collections.defaultdict(list)
for _f, _s in FAMILIES.items():
    for _d in _s["dbs"]:
        DB_FAMILIES[_d].append(_f)

VERSIONED_ACC = re.compile(r"^(.+)\.(\d+)$")


def base(x):
    """An Ensembl stable id without its version (`ENSG00000141510.20` -> `ENSG00000141510`)."""
    i = x.rfind(".")
    return x[:i] if i > 0 and x[i + 1:].isdigit() and x.startswith(("ENS", "LRG_")) else x


def xref_releases(org):
    d = os.path.join(TRUTH, org, "core_xref")
    return sorted(int(os.path.basename(p)[:-7]) for p in glob.glob(os.path.join(d, "*.tsv.gz")))


# The edit-counter bound of rule `accession_versions` (idtrack's `settings.VERSION_EDIT_COUNTER_SHARE`): a database whose stored version differs, for
# one link, in more than this share of the consecutive releases at which the link holds counts its records' edits
# (UniProt's entry version), not versions of a sequence; the two groups lie at most 14.6 % and at least 58 % apart.
EDIT_COUNTER_SHARE = 0.5


def version_change_shares(org, dbs, upto=None):
    """Rule `accession_versions`' edit counters measured from Ensembl's records, not from idtrack: per `db_name` of `dbs`, (pairs, changed): of
    the pairs of consecutive releases of the extracts (up to `upto`) at which one link (the accession, the Ensembl
    object with its version, the database) holds with a version (`xref.version` not empty or 0), the pairs and those
    in which the version differs."""
    out = collections.defaultdict(lambda: [0, 0])
    prev, last_rel = {}, None
    for rel in xref_releases(org):
        if upto is not None and rel > upto:
            continue
        cur = {}
        for ot, oid, gid, db, disp, acc, lab, ver, it in read_xref(org, rel):
            if db not in dbs or ver in ("", "0", "\\N"):
                continue
            k = (db, acc, oid)
            if k in cur:
                continue
            cur[k] = ver
            p = prev.get(k)
            if p is not None:
                out[db][0] += 1
                out[db][1] += p != ver
        prev, last_rel = cur, rel
    return {db: tuple(v) for db, v in out.items()}


def edit_counter_dbs(shares):
    """The databases of `version_change_shares` whose version changes in more than `EDIT_COUNTER_SHARE` of the pairs."""
    return {db for db, (n, c) in shares.items() if n and c / n > EDIT_COUNTER_SHARE}


def read_xref(org, rel):
    """Rows of one release's extract: (object_type, object_id, gene base id, db_name, db_display_name, acc, label,
    version, info_type)."""
    p = os.path.join(TRUTH, org, "core_xref", f"{rel}.tsv.gz")
    with gzip.open(p, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            yield f[0], f[1], base(f[2]), f[3], f[4], f[5].strip(), f[6].strip(), f[7], f[8]


class Scope:
    """Which graph database names the reference tree's configuration includes at which releases (by the assembly of
    the release in the GTF)."""

    def __init__(self, org, clone, asm_of_release):
        import yaml
        p = os.path.join(clone, "idtrack", "default_config", f"{org}_externals_modified.yml")
        d = yaml.safe_load(open(p))[org]
        self.inc = collections.defaultdict(set)   # graph name -> releases where included
        self.form = {}
        for form, dbs in d.items():
            for name, spec in dbs.items():
                for asm, a in (spec.get("Assembly") or {}).items():
                    if a.get("Include"):
                        for r in str(a["Ensembl release"]).split(","):
                            r = int(r)
                            if ASM_KEY.get(asm_of_release.get(r, ""), None) == str(asm):
                                self.inc[name].add(r)
                                self.form[name] = form
        self.yaml = p

    def included(self, name, rel):
        return rel in self.inc.get(name, ())

    def family_included(self, fam, rel):
        return any(self.included(n, rel) for n in FAMILIES[fam]["graph"])


class Xrefs:
    """Ensembl's values per family and release: H[fam][value] = {release: (genes, versions, in_graph)}.

    `genes` are base gene ids (frozensets interned); `versions` the accession versions Ensembl writes at that release;
    `in_graph`: the value is in a database the graph's configuration includes at that release.
    """

    def __init__(self, org, scope, families=None, releases=None, log=print):
        self.org = org
        self.fams = families or list(FAMILIES)
        self.rels = [r for r in xref_releases(org) if releases is None or r in releases]
        self.H = {f: collections.defaultdict(dict) for f in self.fams}
        self.dbs_seen = collections.defaultdict(collections.Counter)   # release -> db_name -> rows
        self.labels = {f: collections.defaultdict(set) for f in ("rfam", "ccds", "mirbase")}  # acc -> labels seen
        intern = {}
        prefix = GENE_PREFIX[org]
        self.non_ens = collections.Counter()
        for rel in self.rels:
            acc_genes = {f: collections.defaultdict(set) for f in self.fams}
            acc_vers = {f: collections.defaultdict(set) for f in self.fams}
            acc_dbs = {f: collections.defaultdict(set) for f in self.fams}
            for ot, oid, gid, db, disp, acc, lab, ver, it in read_xref(org, rel):
                self.dbs_seen[rel][db] += 1
                fs = DB_FAMILIES.get(db)
                if not fs:
                    continue
                if not gid.startswith(prefix):
                    self.non_ens[(rel, db)] += 1
                    continue
                for f in fs:
                    if f not in self.H:
                        continue
                    v = acc if FAMILIES[f]["field"] == "acc" else lab
                    if not v:
                        continue
                    acc_genes[f][v].add(gid)
                    # the graph's name of this row's value (d4c2dda): an accession under `db_name`, a label under
                    # `db_display_name`; where both are the same string, under either
                    # (and only from rows on objects of the form the configuration gives that name)
                    names = ((db, disp) if acc == lab else (db,)) if FAMILIES[f]["field"] == "acc" else \
                        ((disp, db) if acc == lab else (disp,))
                    acc_dbs[f][v].update((n, ot.lower()) for n in names)
                    if FAMILIES[f].get("versioned") and ver not in ("", "0", "\\N"):
                        acc_vers[f][v].add(ver)
                    if f in self.labels and lab and lab != acc:
                        self.labels[f][acc].add(lab)
            for f in self.fams:
                for v, gs in acc_genes[f].items():
                    key = frozenset(gs)
                    key = intern.setdefault(key, key)
                    vers = frozenset(acc_vers[f].get(v, ()))
                    vers = intern.setdefault(vers, vers)
                    in_graph = any(scope.included(n, rel) and scope.form.get(n) == ot for n, ot in acc_dbs[f][v])
                    self.H[f][v][rel] = (key, vers, in_graph)
            log(f"xref {org} {rel}: " + ", ".join(f"{f} {len(acc_genes[f])}" for f in self.fams if acc_genes[f]))

    def at(self, fam, v, rel):
        """(genes, versions, in_graph) of a value at a release, or None."""
        return self.H[fam].get(v, {}).get(rel)

    def last_in_graph(self, fam, v, upto):
        """The latest release <= upto at which the graph holds v (rule `start_release`'s reading without from_release)."""
        best = None
        for r, (gs, vers, ing) in self.H[fam].get(v, {}).items():
            if ing and r <= upto and (best is None or r > best):
                best = r
        return best

    def last_with_version(self, fam, v, ver, upto):
        best = None
        for r, (gs, vers, ing) in self.H[fam].get(v, {}).items():
            if ing and r <= upto and ver in vers and (best is None or r > best):
                best = r
        return best

    def versions_in_graph(self, fam, v, upto=None):
        """The versions Ensembl wrote for v at the releases where the graph's configuration holds it, up to `upto` (the
        graph's last release: the extracts reach one release beyond the graphs' window)."""
        out = set()
        for r, (gs, vers, ing) in self.H[fam].get(v, {}).items():
            if ing and (upto is None or r <= upto):
                out |= vers
        return out


# ------------------------------------------------------------------ the authorities


def mirbase_text(path):
    """A miRBase download as text: mirbase.org serves its files wrapped in HTML (`<p>`, `<br>` for new lines, entities
    such as `&gt;`); the bytes are kept as served (their md5 in SOURCES.tsv) and unwrapped here."""
    import html
    with gzip.open(path, "rt", errors="replace") as fh:
        t = fh.read()
    if t.startswith("<p>"):
        t = html.unescape(t.replace("<p>", "").replace("</p>", "").replace("<br>", "\n"))
    return t


def _rows(path, sep="\t", skip_hash=True):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if skip_hash and line.startswith("#"):
                continue
            yield line.rstrip("\n").split(sep)


class Authority:
    """Each database's own mapping of its identifiers to Ensembl genes (base ids), now (the file's date)."""

    def __init__(self, org):
        self.org = org
        t = os.path.join(TRUTH, org)
        self.files = []
        # NCBI: GeneID -> Ensembl genes; RefSeq accession (unversioned) -> GeneIDs; versions per accession
        self.ncbi = collections.defaultdict(set)
        p = os.path.join(t, "ncbi", "gene2ensembl.tsv.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p):
                if len(f) > 2 and f[2] not in ("-", ""):
                    self.ncbi[f[1]].add(base(f[2]))
        self.refseq_gene = collections.defaultdict(set)
        self.refseq_ver = collections.defaultdict(set)
        p = os.path.join(t, "ncbi", "gene2refseq.tsv.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p):
                for col in (3, 5):
                    if len(f) > col and f[col] not in ("-", ""):
                        m = VERSIONED_ACC.match(f[col])
                        a, v = (m.group(1), m.group(2)) if m else (f[col], None)
                        self.refseq_gene[a].add(f[1])
                        if v:
                            self.refseq_ver[a].add(int(v))
        # the RefSeq release catalog: accession -> versions (the release's records)
        self.catalog_ver = collections.defaultdict(set)
        for p in glob.glob(os.path.join(t, "refseq", "RefSeq-release*.catalog.tsv.gz")):
            self.files.append(p)
            for f in _rows(p):
                if len(f) > 2:
                    m = VERSIONED_ACC.match(f[2])
                    if m:
                        self.catalog_ver[m.group(1)].add(int(m.group(2)))
        # gene_history: discontinued GeneID -> replacement GeneID ('-' none) and date
        self.history = {}
        p = os.path.join(t, "ncbi", "gene_history.tsv.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p):
                if len(f) > 4:
                    self.history[f[2]] = (None if f[1] == "-" else f[1], f[4], f[3])
        # UniProt: accession -> Ensembl genes, entry name; isoforms; secondary; deleted
        self.uniprot = collections.defaultdict(set)
        self.entry = {}
        p = os.path.join(t, "uniprot", "idmapping_selected.tab.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p, skip_hash=False):
                if len(f) > 18:
                    self.entry[f[0]] = f[1]
                    for g in f[18].split(";"):
                        g = g.strip()
                        if g:
                            self.uniprot[f[0]].add(base(g))
        self.secondary = collections.defaultdict(set)   # secondary -> primaries (of this organism)
        p = os.path.join(SHARED, "uniprot", "sec_ac.txt.gz")
        if os.path.exists(p):
            self.files.append(p)
            started = False
            for line in gzip.open(p, "rt", errors="replace"):
                if not started:
                    started = line.startswith("_____")
                    continue
                f = line.split()
                if len(f) == 2 and f[1] in self.entry:
                    self.secondary[f[0]].add(f[1])
        self.isoform_parent = {}
        p = os.path.join(SHARED, "uniprot", "uniprot_sprot_varsplic.fasta.gz")
        if os.path.exists(p):
            self.files.append(p)
            ox = f"OX={TAXON[org]} "
            for line in gzip.open(p, "rt", errors="replace"):
                if line.startswith(">") and ox in line:
                    acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                    self.isoform_parent[acc] = acc.split("-")[0]
        # RNAcentral: URS -> Ensembl genes (through its Ensembl mapping's gene column) ; miRBase acc -> URS
        self.rnacentral = collections.defaultdict(set)
        self.rnacentral_tx = collections.defaultdict(set)
        p = os.path.join(t, "rnacentral", "ensembl.tsv.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p, skip_hash=False):
                if len(f) >= 6:
                    urs = f[0].split("_")[0]
                    self.rnacentral_tx[urs].add(base(f[2]))
                    if f[5].startswith(GENE_PREFIX[org]):
                        self.rnacentral[urs].add(base(f[5]))
        self.mirbase_urs = collections.defaultdict(set)
        p = os.path.join(t, "rnacentral", "mirbase.tsv.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p, skip_hash=False):
                if len(f) >= 3:
                    self.mirbase_urs[f[2]].add(f[0].split("_")[0])
        # miRBase genome file: precursor MI <-> name; mature MIMAT -> name and precursors
        self.mir_name, self.mir_acc, self.mature_of = {}, {}, collections.defaultdict(set)
        code = {"homo_sapiens": "hsa", "mus_musculus": "mmu", "sus_scrofa": "ssc", "macaca_mulatta": "mml",
                "macaca_fascicularis": "mfa"}[org]
        p = os.path.join(t, "mirbase", f"{code}.gff3.gz")
        if os.path.exists(p):
            self.files.append(p)
            for line in mirbase_text(p).split("\n"):
                f = line.split("\t")
                if line.startswith("#") or len(f) < 9:
                    continue
                att = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                if f[2] in ("miRNA_primary_transcript", "miRNA"):
                    self.mir_name[att.get("ID")] = att.get("Name")
                if f[2] == "miRNA" and att.get("Derives_from"):
                    self.mature_of[att.get("ID")].add(att["Derives_from"])
        # the fasta files of every species: accession <-> name of the organism's precursors and mature miRNAs
        for name in ("hairpin.fa.gz", "mature.fa.gz"):
            p = os.path.join(SHARED, "mirbase", "CURRENT", name)
            if os.path.exists(p):
                self.files.append(p)
                for line in mirbase_text(p).split("\n"):
                    if line.startswith(">") and line[1:].startswith(code + "-"):
                        f = line[1:].split()
                        if len(f) > 1:
                            self.mir_name.setdefault(f[1], f[0])
        self.mir_acc = {n: a for a, n in self.mir_name.items() if n}
        # RFAM: accession <-> family id (label)
        self.rfam_label = {}
        p = os.path.join(SHARED, "rfam", "family.txt.gz")
        if os.path.exists(p):
            self.files.append(p)
            for f in _rows(p, skip_hash=False):
                if len(f) > 1:
                    self.rfam_label[f[0]] = f[1]
        # CCDS: CCDS id (versioned) -> Ensembl transcripts
        self.ccds_tx = collections.defaultdict(set)
        for p in glob.glob(os.path.join(t, "ccds", "CCDS2Sequence.current.txt.gz")):
            self.files.append(p)
            for f in _rows(p):
                if len(f) > 4 and f[4].startswith("ENS"):
                    self.ccds_tx[f[0]].add(base(f[4]))
        # VGNC: id -> Ensembl gene
        self.vgnc = collections.defaultdict(set)
        for p in glob.glob(os.path.join(t, "vgnc", "*_vgnc_gene_set_All.tsv")):
            self.files.append(p)
            rows = _rows(p, skip_hash=False)
            head = next(rows)
            if "ensembl_gene_id" in head:
                i, j = head.index("vgnc_id"), head.index("ensembl_gene_id")
                for f in rows:
                    if len(f) > max(i, j):
                        for g in f[j].split("|"):
                            if g.strip():
                                self.vgnc[f[i]].add(base(g.strip()))

    def deleted_uniprot(self, candidates):
        """The accessions among `candidates` that UniProt lists as deleted (Swiss-Prot and TrEMBL lists)."""
        out = set()
        for name in ("delac_sp.txt.gz", "delac_tr.txt.gz"):
            p = os.path.join(SHARED, "uniprot", name)
            if not os.path.exists(p):
                continue
            self.files.append(p)
            for line in gzip.open(p, "rt", errors="replace"):
                a = line.strip()
                if a in candidates:
                    out.add(a)
        return out


def tsv_pairs(org, rel, kind):
    """(value -> genes) of the FTP `tsv/` table of a release (`entrez`, `refseq`, `uniprot`), or None."""
    p = os.path.join(TRUTH, org, "tsv", str(rel), f"{kind}.tsv.gz")
    if not os.path.exists(p):
        return None
    out = collections.defaultdict(set)
    rows = _rows(p, skip_hash=False)
    head = next(rows)
    gi, xi, di = head.index("gene_stable_id"), head.index("xref"), head.index("db_name")
    for f in rows:
        if len(f) > max(gi, xi, di):
            out[(f[di], f[xi])].add(base(f[gi]))
    return out
