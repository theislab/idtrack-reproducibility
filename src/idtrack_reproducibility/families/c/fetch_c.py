"""Truth of family C (external identifiers): fetch each database's own files and Ensembl's xref tables over HTTPS.

Nothing of idtrack is read. Every file fetched is listed in the one truth store's `SOURCES.tsv`
(`<truth>/<organism>/SOURCES.tsv`; files of no single organism in
`truth/_shared/SOURCES.tsv`) with URL, size and md5 of the download and the time; where only an organism's rows of a
file are kept, the row says so and gives the md5 of the kept file. Downloads stream (a file is never held in memory)
and use one connection per host: one subcommand per host, run one at a time per host.

Usage (any Python 3.10+, standard library only):
  fetch_c.py ncbi                 NCBI gene: gene2ensembl, gene2refseq, gene_history (rows of the five taxa);
                                  RefSeq release catalog (catalog, accession2geneid, removed records: five taxa);
                                  CCDS current human and mouse files                      [ftp.ncbi.nlm.nih.gov]
  fetch_c.py ebi                  UniProt idmapping_selected (five taxa), sec_ac, delac_sp, delac_tr, Swiss-Prot
                                  isoforms (varsplic); RFAM family list; RNAcentral's Ensembl and miRBase mappings
                                  (five taxa)                                              [ftp.ebi.ac.uk]
  fetch_c.py mirbase              miRBase CURRENT hairpin.fa, mature.fa, miRNA.diff, genomes/<code>.gff3;
                                  the previous release's miRNA.dat and aliases where they answer   [mirbase.org]
  fetch_c.py ensembl [--releases 87-116] [--kinds tsv,xref] [--organisms ...]
                                  Ensembl FTP `tsv/` entrez, refseq, uniprot files as published; and a compact
                                  extract of the core database's xrefs (xref, object_xref, external_db, gene,
                                  transcript, translation dumps read, not kept): `core_xref/<release>.tsv.gz`
                                  [ftp.ensembl.org]
Resumable: a kept file that exists is skipped; a file that does not exist at its URL is noted in
`<dir>/_unavailable.tsv` of its organism (URL, status, time) and skipped on later runs unless --retry.
"""
import argparse
import datetime
import gzip
import hashlib
import io
import os
import re
import sys
import time
import urllib.error
import urllib.request
import zlib

from idtrack_reproducibility import config as _CFG
TRUTH = _CFG.get().truth
TAXA = {"homo_sapiens": "9606", "mus_musculus": "10090", "sus_scrofa": "9823", "macaca_mulatta": "9544",
        "macaca_fascicularis": "9541"}
ORG_OF_TAXON = {t: o for o, t in TAXA.items()}
SHARED = "_shared"
NCBI = "https://ftp.ncbi.nlm.nih.gov"
EBI = "https://ftp.ebi.ac.uk/pub/databases"
UNIPROT = f"{EBI}/uniprot/current_release/knowledgebase"
MIRBASE = "https://mirbase.org/download"
ENSEMBL = "https://ftp.ensembl.org/pub"
MIRBASE_CODE = {"homo_sapiens": "hsa", "mus_musculus": "mmu", "sus_scrofa": "ssc", "macaca_mulatta": "mml",
                "macaca_fascicularis": "mfa"}
WINDOW = {"homo_sapiens": 48, "mus_musculus": 48, "sus_scrofa": 56, "macaca_mulatta": 48, "macaca_fascicularis": 91}
# Ensembl xref databases kept in the core extract (the rest, GO, InterPro, EMBL, PDB and the like, are not family C's)
KEEP_DB = re.compile(r"(?i)entrez|ncbi|refseq|uniprot|ccds|ott|vega|havana|clone|mirbase|rfam|rnacentral|hgnc|mgi|"
                     r"vgnc|lrg|wikigene")
UA = {"User-Agent": "python-urllib (benchmark truth fetch)"}


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def tdir(org):
    return os.path.join(TRUTH, org)


def add_source(org, what, url, nbytes, md5, fetched):
    p = os.path.join(tdir(org), "SOURCES.tsv")
    os.makedirs(tdir(org), exist_ok=True)
    new = not os.path.exists(p)
    with open(p, "a") as fo:
        if new:
            fo.write("what\turl\tbytes\tmd5\tfetched\n")
        fo.write(f"{what}\t{url}\t{nbytes}\t{md5}\t{fetched}\n")


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def unavailable(org, sub, url, why):
    d = os.path.join(tdir(org), sub)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "_unavailable.tsv"), "a") as fo:
        fo.write(f"{url}\t{why}\t{now()}\n")
    print("unavailable", url, why, flush=True)


def was_unavailable(org, sub, url):
    p = os.path.join(tdir(org), sub, "_unavailable.tsv")
    return os.path.exists(p) and any(line.split("\t")[0] == url for line in open(p))


class Download:
    """One HTTP GET streamed in chunks: md5 and size of the bytes received; `lines()` decompresses gzip on the fly."""

    def __init__(self, url, tries=5, timeout=600):
        self.url, self.md5, self.size = url, hashlib.md5(), 0
        err = None
        for i in range(tries):
            try:
                self.resp = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout)
                self.length = int(self.resp.headers.get("Content-Length") or -1)
                self.headers = dict(self.resp.headers)
                return
            except urllib.error.HTTPError as e:
                if e.code in (403, 404, 410):
                    raise
                err = e
            except Exception as e:  # noqa: BLE001
                err = e
            time.sleep(10 * (i + 1))
        raise err

    def chunks(self):
        while True:
            b = self.resp.read(1 << 20)
            if not b:
                break
            self.md5.update(b)
            self.size += len(b)
            yield b
        if self.length >= 0 and self.size != self.length:
            raise IOError(f"{self.url}: received {self.size} of {self.length} bytes")

    def lines(self, gz, bz=False):
        """Text lines (without newline) of the body; gzip (also of several members) decompressed as it arrives, after
        bzip2 when `bz` (Ensembl's dumps of some releases are `.txt.gz.bz2`)."""
        import bz2
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS) if gz else None
        bdec = bz2.BZ2Decompressor() if bz else None
        rest = b""
        for b in self.chunks():
            if bdec is not None:
                out = b""
                while b:
                    out += bdec.decompress(b)
                    b = b""
                    if bdec.eof:  # a next bzip2 stream
                        b = bdec.unused_data
                        bdec = bz2.BZ2Decompressor()
                b = out
            if dec is not None:
                out = dec.decompress(b)
                while dec.unused_data:  # a next gzip member
                    tail = dec.unused_data
                    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
                    out += dec.decompress(tail)
                b = out
            rest += b
            parts = rest.split(b"\n")
            rest = parts.pop()
            for p in parts:
                yield p.decode("utf-8", "replace")
        if dec is not None:
            rest += dec.flush()
        if rest:
            yield rest.decode("utf-8", "replace")


def open_dump(url):
    """A dump at its URL, or at the URL with `.bz2` (releases whose dumps are `.gz.bz2`): (Download, url, bz)."""
    try:
        return Download(url), url, False
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    return Download(url + ".bz2"), url + ".bz2", True


def save_whole(org, rel_path, url, gzip_it=False, what=None):
    """A file kept as published (or gzipped, `gzip_it`); one SOURCES row (md5 of the download)."""
    dst = os.path.join(tdir(org), rel_path)
    if os.path.exists(dst):
        return True
    sub = os.path.dirname(rel_path) or "."
    if was_unavailable(org, sub, url) and not ARGS.retry:
        return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        d = Download(url)
    except urllib.error.HTTPError as e:
        unavailable(org, sub, url, f"HTTP {e.code}")
        return False
    tmp = dst + ".tmp"
    with (gzip.open(tmp, "wb") if gzip_it else open(tmp, "wb")) as fo:
        for b in d.chunks():
            fo.write(b)
    os.replace(tmp, dst)
    note = what or rel_path
    if gzip_it:
        note += f" (gzipped here; kept file md5 {file_md5(dst)})"
    add_source(org, note, url, d.size, d.md5.hexdigest(), now())
    print("saved", org, rel_path, d.size, flush=True)
    return True


def save_split(url, rel_path, key, gz=True, header=None, orgs=None, what="", keep_header=True):
    """One file streamed once; each organism's rows (`key(line)` -> organism or None) kept in
    `<organism>/<rel_path>` (gzip). `header`: number of leading lines kept in every output."""
    orgs = orgs or list(TAXA)
    dsts = {o: os.path.join(tdir(o), rel_path) for o in orgs}
    if all(os.path.exists(p) for p in dsts.values()):
        return
    try:
        d = Download(url)
    except urllib.error.HTTPError as e:
        for o in orgs:
            unavailable(o, os.path.dirname(rel_path), url, f"HTTP {e.code}")
        return
    outs, counts = {}, {o: 0 for o in orgs}
    for o, p in dsts.items():
        os.makedirs(os.path.dirname(p), exist_ok=True)
        outs[o] = gzip.open(p + ".tmp", "wt")
    n = 0
    for i, line in enumerate(d.lines(gz)):
        if header and i < header:
            if keep_header:
                for fo in outs.values():
                    fo.write(line + "\n")
            continue
        o = key(line)
        if o in outs:
            outs[o].write(line + "\n")
            counts[o] += 1
        n += 1
        if n % 20_000_000 == 0:
            print(f"  {rel_path}: {n:,} lines, {d.size:,} bytes", flush=True)
    for o, fo in outs.items():
        fo.close()
        os.replace(dsts[o] + ".tmp", dsts[o])
    t = now()
    for o in orgs:
        add_source(o, f"{rel_path} (rows of {what or 'the organism'} from {url.rsplit('/', 1)[1]}: {counts[o]:,} of "
                      f"{n:,}; kept file md5 {file_md5(dsts[o])})", url, d.size, d.md5.hexdigest(), t)
    print("saved", rel_path, counts, flush=True)


def tsv_taxon(col, sep="\t"):
    if col == 0:
        return lambda line: ORG_OF_TAXON.get(line.partition(sep)[0].strip())

    def key(line):
        f = line.split(sep)
        return ORG_OF_TAXON.get(f[col].strip()) if len(f) > col else None
    return key


# ---------------------------------------------------------------- NCBI


def ncbi():
    g = f"{NCBI}/gene/DATA"
    save_split(f"{g}/gene2ensembl.gz", "ncbi/gene2ensembl.tsv.gz", tsv_taxon(0), header=1, what="taxon")
    save_split(f"{g}/gene_history.gz", "ncbi/gene_history.tsv.gz", tsv_taxon(0), header=1, what="taxon")
    save_split(f"{g}/gene2refseq.gz", "ncbi/gene2refseq.tsv.gz", tsv_taxon(0), header=1, what="taxon")
    cat = f"{NCBI}/refseq/release/release-catalog"
    rel = urllib.request.urlopen(urllib.request.Request(f"{NCBI}/refseq/release/RELEASE_NUMBER", headers=UA),
                                 timeout=120).read().decode().strip()
    for o in TAXA:
        p = os.path.join(tdir(o), "refseq", "RELEASE_NUMBER")
        if not os.path.exists(p):
            os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "w").write(rel + "\n")
            add_source(o, "refseq/RELEASE_NUMBER", f"{NCBI}/refseq/release/RELEASE_NUMBER", len(rel) + 1,
                       hashlib.md5((rel + "\n").encode()).hexdigest(), now())
    save_split(f"{cat}/release{rel}.accession2geneid.gz", f"refseq/release{rel}.accession2geneid.tsv.gz", tsv_taxon(0),
               what="taxon")
    save_split(f"{cat}/release{rel}.removed-records.gz", f"refseq/release{rel}.removed-records.tsv.gz", tsv_taxon(0),
               what="taxon")
    save_split(f"{cat}/RefSeq-release{rel}.catalog.gz", f"refseq/RefSeq-release{rel}.catalog.tsv.gz", tsv_taxon(0),
               what="taxon")
    for org, sub in (("homo_sapiens", "current_human"), ("mus_musculus", "current_mouse")):
        listing = urllib.request.urlopen(urllib.request.Request(f"{NCBI}/pub/CCDS/{sub}/", headers=UA),
                                         timeout=120).read().decode()
        for fn in sorted(set(re.findall(r'href="([^"/]+\.current\.txt)"', listing))):
            save_whole(org, f"ccds/{fn}.gz", f"{NCBI}/pub/CCDS/{sub}/{fn}", gzip_it=True, what=f"ccds/{fn}.gz")


# ---------------------------------------------------------------- EBI


def ebi():
    save_split(f"{UNIPROT}/idmapping/idmapping_selected.tab.gz", "uniprot/idmapping_selected.tab.gz", tsv_taxon(12),
               what="taxon (column 13)")
    rel = urllib.request.urlopen(urllib.request.Request(f"{UNIPROT}/complete/reldate.txt", headers=UA),
                                 timeout=120).read()
    p = os.path.join(tdir(SHARED), "uniprot", "reldate.txt")
    if not os.path.exists(p):
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "wb").write(rel)
        add_source(SHARED, "uniprot/reldate.txt", f"{UNIPROT}/complete/reldate.txt", len(rel),
                   hashlib.md5(rel).hexdigest(), now())
    save_whole(SHARED, "uniprot/sec_ac.txt.gz", f"{UNIPROT}/complete/docs/sec_ac.txt", gzip_it=True)
    save_whole(SHARED, "uniprot/delac_sp.txt.gz", f"{UNIPROT}/complete/docs/delac_sp.txt", gzip_it=True)
    save_whole(SHARED, "uniprot/delac_tr.txt.gz", f"{UNIPROT}/complete/docs/delac_tr.txt.gz")
    save_whole(SHARED, "uniprot/uniprot_sprot_varsplic.fasta.gz", f"{UNIPROT}/complete/uniprot_sprot_varsplic.fasta.gz")
    save_whole(SHARED, "rfam/family.txt.gz", f"{EBI}/Rfam/CURRENT/database_files/family.txt.gz")
    save_whole(SHARED, "rfam/README", f"{EBI}/Rfam/CURRENT/README")
    rc = f"{EBI}/RNAcentral/current_release"
    save_whole(SHARED, "rnacentral/release_notes.txt", f"{rc}/release_notes.txt")
    save_split(f"{rc}/id_mapping/database_mappings/ensembl.tsv", "rnacentral/ensembl.tsv.gz", tsv_taxon(3), gz=False,
               what="taxon (column 4)")
    save_split(f"{rc}/id_mapping/database_mappings/mirbase.tsv", "rnacentral/mirbase.tsv.gz", tsv_taxon(3), gz=False,
               what="taxon (column 4)")


# ---------------------------------------------------------------- miRBase


def mirbase():
    for fn in ("hairpin.fa", "mature.fa", "miRNA.diff", "README", "miRNA.dead"):
        save_whole(SHARED, f"mirbase/CURRENT/{fn}.gz", f"{MIRBASE}/CURRENT/{fn}", gzip_it=True)
    for org, code in MIRBASE_CODE.items():
        save_whole(org, f"mirbase/{code}.gff3.gz", f"{MIRBASE}/CURRENT/genomes/{code}.gff3", gzip_it=True)
    for fn in ("miRNA.dat.gz", "aliases.txt.gz", "miRNA.dat", "aliases.txt"):
        save_whole(SHARED, f"mirbase/22.1/{fn if fn.endswith('.gz') else fn + '.gz'}",
                   f"{MIRBASE}/PREVIOUS_RELEASES/22.1/{fn}", gzip_it=not fn.endswith(".gz"))


# ---------------------------------------------------------------- Ensembl


def http_text(url):
    d = Download(url)
    return b"".join(d.chunks()).decode("utf-8", "replace")


def core_db(org, rel):
    html = http_text(f"{ENSEMBL}/release-{rel}/mysql/")
    # folders only (release 89 also lists `<db>.gz` files); human from 112 also holds the GRCh37 archive core `_37`
    # beside the release's own `_38`: the newest assembly is the release's
    dbs = sorted(set(re.findall(rf'href="({org}_core_{rel}_[^"/.]+)/"', html)),
                 key=lambda x: [int(p) if p.isdigit() else -1 for p in x.rsplit("_", 1)[1:]])
    return dbs[-1] if dbs else None


def columns(sql, table):
    m = re.search(rf"CREATE TABLE `{table}` \((.*?)\n\)", sql, re.S)
    return [c for c in re.findall(r"^\s+`(\w+)`", m.group(1), re.M)] if m else None


def dump_rows(org, rp, base, sql, table, keep=None):
    """Rows (dicts) of one core table dump, streamed; SOURCES row for the dump (read, not kept)."""
    cols = columns(sql, table)
    if cols is None:
        return None
    d, url, bz = open_dump(f"{base}/{table}.txt.gz")
    held = ""
    for line in d.lines(True, bz):
        if line.endswith("\\") and not line.endswith("\\\\"):  # an escaped newline inside a text field
            held += line[:-1] + "\n"
            continue
        line, held = held + line, ""
        f = line.split("\t")
        if line and len(f) == len(cols):
            r = dict(zip(cols, f))
            if keep is None or keep(r):
                yield r
    add_source(org, f"{rp} (read from {url.rsplit('/', 1)[1]}; not kept)", url, d.size, d.md5.hexdigest(), now())


def ensembl_xref(org, rel):
    """Compact extract of the core xrefs of a release: one row per (Ensembl object, xref) of the kept databases."""
    rp = f"core_xref/{rel}.tsv.gz"
    dst = os.path.join(tdir(org), rp)
    marker = f"{ENSEMBL}/release-{rel}/mysql/{org}_core_{rel}"
    if os.path.exists(dst) or (was_unavailable(org, "core_xref", marker) and not ARGS.retry):
        return
    db = core_db(org, rel)
    if db is None:
        unavailable(org, "core_xref", marker, "no single core database")
        return
    base = f"{ENSEMBL}/release-{rel}/mysql/{db}"
    sd, sqlurl, sbz = open_dump(f"{base}/{db}.sql.gz")
    sql = "\n".join(sd.lines(True, sbz))
    add_source(org, f"{rp} (schema read from {sqlurl.rsplit('/', 1)[1]}; not kept)", sqlurl, sd.size, sd.md5.hexdigest(), now())
    dbs = {}
    for r in dump_rows(org, rp, base, sql, "external_db"):
        if (KEEP_DB.search(r["db_name"]) or KEEP_DB.search(r.get("db_display_name", ""))) \
                and not r["db_name"].endswith("_trans_name"):
            dbs[r["external_db_id"]] = (r["db_name"], r.get("db_display_name", ""))
    xrefs = {}
    for r in dump_rows(org, rp, base, sql, "xref", keep=lambda r: r["external_db_id"] in dbs):
        xrefs[r["xref_id"]] = (r["external_db_id"], r["dbprimary_acc"], r["display_label"], r.get("version", ""),
                               r.get("info_type", ""))
    obj = [(r["ensembl_id"], r["ensembl_object_type"], r["xref_id"])
           for r in dump_rows(org, rp, base, sql, "object_xref", keep=lambda r: r["xref_id"] in xrefs)]
    stable_in_table = "stable_id" in (columns(sql, "gene") or [])

    def sid(r):
        v = r.get("version")
        return r["stable_id"] + (f".{v}" if v not in (None, "", "\\N") else "")

    def read(kind, parent=None):
        """stable ids of one object table (`<kind>_stable_id` before release 65) and the parent column, one pass each"""
        par = {}
        if stable_in_table:
            st = {}
            for r in dump_rows(org, rp, base, sql, kind):
                st[r[f"{kind}_id"]] = sid(r)
                if parent:
                    par[r[f"{kind}_id"]] = r[parent]
            return st, par
        st = {r[f"{kind}_id"]: sid(r) for r in dump_rows(org, rp, base, sql, f"{kind}_stable_id")}
        if parent:
            par = {r[f"{kind}_id"]: r[parent] for r in dump_rows(org, rp, base, sql, kind)}
        return st, par

    gst, _ = read("gene")
    tst, t2g = read("transcript", "gene_id")
    pst, p2t = read("translation", "transcript_id")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    n = 0
    with gzip.open(dst + ".tmp", "wt") as fo:
        fo.write("object_type\tobject_id\tgene_id\tdb_name\tdb_display_name\tacc\tlabel\tversion\tinfo_type\n")
        for eid, typ, xid in obj:
            dbid, acc, lab, ver, it = xrefs[xid]
            if typ == "Gene":
                oid, gid = gst.get(eid, ""), gst.get(eid, "")
            elif typ == "Transcript":
                oid, gid = tst.get(eid, ""), gst.get(t2g.get(eid, ""), "")
            elif typ == "Translation":
                oid, gid = pst.get(eid, ""), gst.get(t2g.get(p2t.get(eid, ""), ""), "")
            else:
                continue
            fo.write(f"{typ}\t{oid}\t{gid}\t{dbs[dbid][0]}\t{dbs[dbid][1]}\t{acc}\t{lab}\t{ver}\t{it}\n")
            n += 1
    os.replace(dst + ".tmp", dst)
    add_source(org, f"{rp} (extract written here: {n:,} rows; md5 {file_md5(dst)})", base, 0, "-", now())
    print("saved", org, rp, n, flush=True)


def ensembl_tsv(org, rel):
    sub = f"tsv/{rel}"
    url = f"{ENSEMBL}/release-{rel}/tsv/{org}/"
    d = os.path.join(tdir(org), sub)
    if os.path.exists(os.path.join(d, "_done")) or (was_unavailable(org, sub, url) and not ARGS.retry):
        return
    try:
        html = http_text(url)
    except urllib.error.HTTPError as e:
        unavailable(org, sub, url, f"HTTP {e.code}")
        return
    files = sorted(set(re.findall(r'href="([^"/]+\.(entrez|refseq|uniprot)\.tsv\.gz)"', html)))
    if not files:
        unavailable(org, sub, url, "no entrez/refseq/uniprot tsv in the folder")
        return
    for fn, kind in files:
        save_whole(org, f"{sub}/{kind}.tsv.gz", url + fn, what=f"{sub}/{kind}.tsv.gz ({fn})")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "_done"), "w").write(now() + "\n")


def ensembl(rels, kinds, orgs):
    for org in orgs:
        for rel in rels:
            if rel < WINDOW[org]:
                continue
            if "tsv" in kinds and rel >= 87:
                ensembl_tsv(org, rel)
            if "xref" in kinds:
                ensembl_xref(org, rel)


def parse_rels(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def main():
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["ncbi", "ebi", "mirbase", "ensembl"])
    ap.add_argument("--releases", default="87-116")
    ap.add_argument("--kinds", default="tsv,xref")
    ap.add_argument("--organisms", default=",".join(TAXA))
    ap.add_argument("--retry", action="store_true", help="try again URLs noted as unavailable")
    ARGS = ap.parse_args()
    print(now(), "start", ARGS.what, flush=True)
    if ARGS.what == "ncbi":
        ncbi()
    elif ARGS.what == "ebi":
        ebi()
    elif ARGS.what == "mirbase":
        mirbase()
    else:
        ensembl(parse_rels(ARGS.releases), ARGS.kinds.split(","), ARGS.organisms.split(","))
    print(now(), "done", ARGS.what, flush=True)


ARGS = None
if __name__ == "__main__":
    main()
