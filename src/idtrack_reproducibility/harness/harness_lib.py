"""Shared helpers of the accuracy harness: locations, downloads, table readers, intervals. Imports nothing of idtrack.

The truth store is `<truth>/<organism>/` and the query sets `<harness>/queries/<organism>/` of the configuration
(README, "Configuration")."""
import collections
import csv
import gzip
import json
import os
import re
import time
import urllib.request

from idtrack_reproducibility import config as CFG

csv.field_size_limit(1 << 30)

ORGANISMS = {
    # name: (BioMart dataset, gene id prefix, latest release scored)
    "sus_scrofa": ("sscrofa_gene_ensembl", "ENSSSCG", 115),
    "mus_musculus": ("mmusculus_gene_ensembl", "ENSMUSG", 115),
    "homo_sapiens": ("hsapiens_gene_ensembl", "ENSG", 115),
    "macaca_mulatta": ("mmulatta_gene_ensembl", "ENSMMUG", 115),
    "macaca_fascicularis": ("mfascicularis_gene_ensembl", "ENSMFAG", 115),
}


# Queries of published sets that ask for something the organism's graphs do not offer at all. They stay in the set
# (published sets are never rewritten) and are left out when scoring: (organism, kind, final_database) -> reason.
EXCLUDED = {
    ("mus_musculus", "gene_to_db", "HGNC Symbol"): "mouse graphs have no HGNC Symbol database (BioMart lists HGNC "
                                                   "symbols of human orthologs for mouse genes)",
    ("mus_musculus", "old_gene_to_db", "HGNC Symbol"): "as above",
}
NOT_A_DATABASE_OF = {"mus_musculus": {"HGNC Symbol"}}  # used by make_queries for new sets
# A stable id of Ensembl's scheme, without its version: Ensembl's genes, transcripts, translations and exons
# (`ENS<species><G|T|P|E><11 digits>`), and LRG records with their genes, transcripts and proteins as Ensembl names
# them (`LRG_<n>`, `LRG_<n>_g<k>`, `LRG_<n>t<k>`, `LRG_<n>p<k>`). One definition for every reader of ids.
STABLE_ID_FORM = r"(?:ENS[A-Z]*[GTPE]\d{11}|LRG_\d+(?:_g\d+|[tp]\d+)?)"
# The judges of revision sets (a manifest's `judge`): `78`, the judge corrected where it was wrong (rule
# `judge_corrections`);
# `rules`, the judge following the rules idtrack decided after 71690ab (make_queries.RULES_JUDGE). The judge
# `published` leaves out every set of these.
REVISED_JUDGES = ("78", "rules")


def truth_dir(org):
    return os.path.join(CFG.get().truth, org)


def query_dir(org):
    return os.path.join(CFG.get().harness, "queries", org)


def runs_dir():
    return os.path.join(CFG.get().harness, "runs")


def set_names(d):
    """The published sets of a query folder in the order they were published (`v2` before `v10`).

    A revision set replaces the lines of earlier sets with the same qid, so the order matters; names are `v<n>`
    (a plain string sort would put `v10` before `v2`).
    """
    names = [f[: -len(".jsonl.gz")] for f in os.listdir(d) if f.endswith(".jsonl.gz")]

    def key(n):
        m = re.match(r"^v(\d+)$", n)
        return (0, int(m.group(1)), n) if m else (1, 0, n)

    return sorted(names, key=key)


def http_get(url, data=None, tries=5, timeout=600, headers=None):
    err = None
    for i in range(tries):
        try:
            h = {"User-Agent": "python-urllib"}
            h.update(headers or {})
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(5 * (i + 1))
    raise err


def write_tsv_gz(path, header, rows):
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    os.replace(tmp, path)


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def add_source(org, what, url, nbytes, md5, fetched):
    p = os.path.join(truth_dir(org), "SOURCES.tsv")
    new = not os.path.exists(p)
    with open(p, "a") as fo:
        if new:
            fo.write("what\turl\tbytes\tmd5\tfetched\n")
        fo.write(f"{what}\t{url}\t{nbytes}\t{md5}\t{fetched}\n")


def wilson(k, n, z=1.959964):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (max(0.0, c - h), min(1.0, c + h))


# ---------- truth readers ----------

def gtf_releases(org):
    d = os.path.join(truth_dir(org), "gtf")
    if not os.path.isdir(d):
        return []
    return sorted(int(f[6:-7]) for f in os.listdir(d) if f.startswith("genes_") and f.endswith(".tsv.gz"))


def gtf_genes(org, rel, region=None):
    """{gene_id: row} of one release from Ensembl's GTF (region 'primary' / 'patch' / None = both)."""
    p = os.path.join(truth_dir(org), "gtf", f"genes_{rel}.tsv.gz")
    if not os.path.exists(p):
        return None
    return {r["gene_id"]: r for r in read_tsv(p) if region is None or r["region"] == region}


def gtf_tx(org, rel):
    p = os.path.join(truth_dir(org), "gtf", f"tx_{rel}.tsv.gz")
    if not os.path.exists(p):
        return None
    return read_tsv(p)


def biomart(org, rel, group):
    p = os.path.join(truth_dir(org), "biomart", str(rel), f"{group}.tsv")
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        r = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(r)
        rows = [row + [""] * (len(header) - len(row)) for row in r if row]
    return header, rows


def biomart_pairs(org, rel, group, col=1):
    """(gene -> {value}, value -> {gene}) for non-empty values of one BioMart table column."""
    t = biomart(org, rel, group)
    if t is None:
        return None
    g2v, v2g = collections.defaultdict(set), collections.defaultdict(set)
    for row in t[1]:
        if len(row) > col and row[col] != "":
            g2v[row[0]].add(row[col])
            v2g[row[col]].add(row[0])
    return g2v, v2g


def biomart_genes(org, rel):
    t = biomart(org, rel, "genes")
    if t is None:
        return None
    header, rows = t

    def col(*names):
        for n in names:
            if n in header:
                return header.index(n)
        return None

    ci = dict(gene_id=col("ensembl_gene_id"), version=col("ensembl_gene_id_version"),
              name=col("external_gene_name", "external_gene_id"), source=col("external_gene_source", "external_gene_db"),
              biotype=col("gene_biotype"))
    return [{k: (r[i] if i is not None and i < len(r) else "") for k, i in ci.items()} for r in rows]


def vgnc_symbols(org, rel):
    """(gene -> {symbol}, symbol -> {gene}) at a BioMart release: BioMart's VGNC id of the gene (`vgnc` table) read
    through VGNC's own gene set (`vgnc/<species>_vgnc_gene_set_All.tsv`: the id's symbol and its previous symbols).
    None when either table is missing."""
    ids = biomart_pairs(org, rel, "vgnc", 1)
    d = os.path.join(truth_dir(org), "vgnc")
    files = [f for f in os.listdir(d) if f.endswith("_vgnc_gene_set_All.tsv")] if os.path.isdir(d) else []
    if ids is None or not files:
        return None
    sym = collections.defaultdict(set)
    for r in read_tsv(os.path.join(d, files[0])):
        sym[r["vgnc_id"]].add(r["symbol"])
        sym[r["vgnc_id"]] |= {x for x in r["prev_symbol"].split("|") if x}
    g2v, v2g = collections.defaultdict(set), collections.defaultdict(set)
    for g, vids in ids[0].items():
        for vid in vids:
            for s in sym.get(vid, ()):
                g2v[g].add(s)
                v2g[s].add(g)
    return g2v, v2g


def vega_names(org, rel):
    """(gene -> {Vega gene name}, name -> {gene}) of the core xref tables of a release (`vega/<release>.tsv`), or
    None when not fetched."""
    p = os.path.join(truth_dir(org), "vega", f"{rel}.tsv")
    if not os.path.exists(p):
        return None
    g2v, v2g = collections.defaultdict(set), collections.defaultdict(set)
    for r in read_tsv(p):
        if r["vega_name"]:
            g2v[r["gene_id"]].add(r["vega_name"])
            v2g[r["vega_name"]].add(r["gene_id"])
    return g2v, v2g


def biomart_releases(org):
    d = os.path.join(truth_dir(org), "biomart")
    if not os.path.isdir(d):
        return []
    return sorted(int(x) for x in os.listdir(d) if x.isdigit() and os.path.exists(os.path.join(d, x, "genes.tsv")))


def core_tables(org, release=None):
    """Ensembl's stable_id_event (gene rows) and mapping_session of the core database of `release`
    (default: the organism's latest scored release)."""
    release = release or ORGANISMS[org][2]
    raw = os.path.join(truth_dir(org), "core")
    ev_files = sorted(f for f in os.listdir(raw) if f.endswith(".stable_id_event.txt.gz") and f"_core_{release}_" in f)
    if len(ev_files) != 1:
        raise FileNotFoundError(f"no single core stable_id_event of release {release} for {org}: {ev_files}")
    core = ev_files[0][: -len(".stable_id_event.txt.gz")]
    sessions = {}
    with gzip.open(os.path.join(raw, core + ".mapping_session.txt.gz"), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sessions[int(f[0])] = dict(old_db=f[1], new_db=f[2], old_release=f[3], new_release=f[4],
                                       old_assembly=f[5], new_assembly=f[6], created=f[7])
    ev = []
    with gzip.open(os.path.join(raw, core + ".stable_id_event.txt.gz"), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[5] != "gene":
                continue
            ev.append((None if f[0] == "\\N" else f[0], f[1], None if f[2] == "\\N" else f[2], f[3], int(f[4]),
                       float(f[6])))
    return core, ev, sessions


def rel_int(x):
    import re
    m = re.match(r"^(\d+)", str(x))
    return int(m.group(1)) if m else -1


def assembly_of_gtf(url, rel):
    """Assembly name in a GTF file name: `Sus_scrofa.Sscrofa10.2.75.gtf.gz` -> `Sscrofa10.2`."""
    fn = url.rsplit("/", 1)[1]
    head = fn[: -len(f".{rel}.gtf.gz")] if fn.endswith(f".{rel}.gtf.gz") else fn
    return head.split(".", 1)[1] if "." in head else head


def gtf_assemblies(org):
    out = {}
    import re
    p = os.path.join(truth_dir(org), "SOURCES.tsv")
    for row in read_tsv(p):
        m = re.match(r"gtf/(\d+)$", row["what"])
        if m and "chr_patch" not in row["url"]:
            out[int(m.group(1))] = assembly_of_gtf(row["url"], int(m.group(1)))
    return out


def dumps(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"))
