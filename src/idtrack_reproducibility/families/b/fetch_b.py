"""Truth of family B (names over time): the authorities' own name files, fetched over HTTPS.

Nothing of idtrack is read. Every file is listed in its organism's `SOURCES.tsv` of the one truth store
(`<truth>/<organism>/`) with URL, size and md5 of the download and the time. Downloads
stream to disk (`.tmp`, then renamed), one connection per host: one subcommand per host, run one at a time.

Usage (standard library only):
  fetch_b.py hgnc     HGNC complete set (symbol, prev_symbol, alias_symbol, ensembl_gene_id, date columns) and the
                      withdrawn list; the quarterly archive where the bucket lists it      [storage.googleapis.com]
  fetch_b.py mgi      MGI MRK_List1.rpt (markers with withdrawn symbols and synonyms), MRK_List2.rpt (current),
                      MRK_ENSEMBL.rpt (MGI's own Ensembl gene ids)                         [www.informatics.jax.org]
  fetch_b.py ncbi     NCBI gene_info of the five taxa (symbol, synonyms, dbXrefs)          [ftp.ncbi.nlm.nih.gov]
  fetch_b.py synonyms [--releases 115|all|<r,r,...>] [--organisms <o>,...]
                      Ensembl's synonyms of every xref of a release (core MySQL dumps: external_synonym, xref,
                      object_xref, external_db, gene, transcript, translation; read, not kept): `core_synonyms/<R>.tsv.gz`
                      (synonym, gene, the xref's database and label). BioMart's synonym table holds only the synonyms
                      of a gene's display xref (cynomolgus 115: 113); idtrack reads those of every xref  [ftp.ensembl.org]
VGNC's gene sets (pig, rhesus, cynomolgus) are in the store already (`vgnc/`, fetched by the harness).
A kept file that exists is skipped.
"""
import datetime
import hashlib
import os
import re
import sys
import time
import urllib.request

from idtrack_reproducibility import config as _CFG
TRUTH = _CFG.get().truth
HGNC = "https://storage.googleapis.com/public-download-files/hgnc"
MGI = "https://www.informatics.jax.org/downloads/reports"
NCBI = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia"
SPECIES = {"homo_sapiens": "Homo_sapiens", "mus_musculus": "Mus_musculus", "sus_scrofa": "Sus_scrofa",
           "macaca_mulatta": "Macaca_mulatta", "macaca_fascicularis": "Macaca_fascicularis"}


def add_source(org, what, url, nbytes, md5, fetched, note=""):
    p = os.path.join(TRUTH, org, "SOURCES.tsv")
    with open(p, "a") as fo:
        fo.write(f"{what}\t{url}\t{nbytes}\t{md5}\t{fetched}{(' ' + note) if note else ''}\n")


def fetch(org, rel_path, url, tries=5):
    """Stream url to truth/<org>/<rel_path>; record it. Returns the path (None when the URL does not exist)."""
    path = os.path.join(TRUTH, org, rel_path)
    if os.path.exists(path):
        print("have", path)
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    err = None
    for i in range(tries):
        try:
            h, n = hashlib.md5(), 0
            req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
            with urllib.request.urlopen(req, timeout=600) as r, open(path + ".tmp", "wb") as fo:
                for b in iter(lambda: r.read(1 << 20), b""):
                    fo.write(b)
                    h.update(b)
                    n += len(b)
            os.replace(path + ".tmp", path)
            add_source(org, rel_path, url, n, h.hexdigest(), datetime.datetime.now().isoformat(timespec="seconds"))
            print("fetched", path, n)
            return path
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print("absent", url)
                return None
            err = e
        except Exception as e:  # noqa: BLE001
            err = e
        time.sleep(10 * (i + 1))
    raise err


def hgnc():
    o = "homo_sapiens"
    fetch(o, "hgnc/hgnc_complete_set.txt", f"{HGNC}/tsv/tsv/hgnc_complete_set.txt")
    fetch(o, "hgnc/withdrawn.txt", f"{HGNC}/tsv/tsv/withdrawn.txt")
    # the quarterly archive: listed through the bucket's XML API (no file is fetched that the listing does not name)
    try:
        req = urllib.request.Request("https://storage.googleapis.com/public-download-files?prefix=hgnc/archive/archive/quarterly/tsv/",
                                     headers={"User-Agent": "python-urllib"})
        xml = urllib.request.urlopen(req, timeout=120).read().decode()
        keys = re.findall(r"<Key>([^<]+)</Key>", xml)
        print("archive listing:", len(keys), "keys", keys[:3], keys[-3:])
        with open(os.path.join(TRUTH, o, "hgnc", "archive_listing.txt"), "w") as fo:
            fo.write("\n".join(keys) + "\n")
    except Exception as e:  # noqa: BLE001
        print("archive listing failed:", type(e).__name__, e)


def mgi():
    o = "mus_musculus"
    for f in ("MRK_List1.rpt", "MRK_List2.rpt", "MRK_ENSEMBL.rpt"):
        fetch(o, f"mgi/{f}", f"{MGI}/{f}")


def ncbi():
    missing = {}
    for o, sp in SPECIES.items():
        if fetch(o, f"ncbi/{sp}.gene_info.gz", f"{NCBI}/{sp}.gene_info.gz") is None:
            missing[o] = sp
    if missing:
        ncbi_all(missing)


def ncbi_all(missing):
    """The taxa NCBI gives no file of their own (the macaques): their rows of All_Mammalia.gene_info.gz, streamed."""
    import gzip
    import zlib
    taxa = {"macaca_mulatta": "9544", "macaca_fascicularis": "9541"}
    want = {taxa[o]: o for o in missing if o in taxa}
    url = f"{NCBI}/All_Mammalia.gene_info.gz"
    outs = {o: gzip.open(os.path.join(TRUTH, o, f"ncbi/{missing[o]}.gene_info.gz.tmp"), "wt") for o in want.values()}
    counts = dict.fromkeys(want.values(), 0)
    h, n, total = hashlib.md5(), 0, 0
    dec, buf, head = zlib.decompressobj(16 + zlib.MAX_WBITS), b"", None
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=600) as r:
        for b in iter(lambda: r.read(1 << 20), b""):
            h.update(b)
            n += len(b)
            buf += dec.decompress(b)
            *lines, buf = buf.split(b"\n")
            for ln in lines:
                t = ln.decode("utf-8", "replace")
                if head is None:
                    head = t
                    for fo in outs.values():
                        fo.write(head + "\n")
                    continue
                total += 1
                o = want.get(t.split("\t", 1)[0])
                if o:
                    outs[o].write(t + "\n")
                    counts[o] += 1
    fetched = datetime.datetime.now().isoformat(timespec="seconds")
    for o, fo in outs.items():
        fo.close()
        path = os.path.join(TRUTH, o, f"ncbi/{missing[o]}.gene_info.gz")
        os.replace(path + ".tmp", path)
        km = hashlib.md5(open(path, "rb").read()).hexdigest()
        add_source(o, f"ncbi/{missing[o]}.gene_info.gz (rows of taxon {taxa[o]} from All_Mammalia.gene_info.gz: {counts[o]:,} of "
                   f"{total:,}; kept file md5 {km})", url, n, h.hexdigest(), fetched)
        print("kept", path, counts[o])


WINDOW = {"homo_sapiens": 48, "mus_musculus": 48, "sus_scrofa": 56, "macaca_mulatta": 48, "macaca_fascicularis": 91}
ORGANISMS = ("macaca_fascicularis", "sus_scrofa", "macaca_mulatta", "mus_musculus", "homo_sapiens")  # smallest first


def synonyms(releases=(115,), organisms=ORGANISMS):
    """`core_synonyms/<R>.tsv.gz` of the organisms (default every one): Ensembl's synonyms of every xref of the
    release, by gene. The core dumps are read by family C's reader (`fetch_c`: the release's core database, its
    dumps as `.txt.gz` or `.txt.gz.bz2`, the stable-id tables before release 65), each recorded in SOURCES.tsv."""
    import collections
    import gzip
    from idtrack_reproducibility.families.c import fetch_c as FC
    for o in [x for x in ORGANISMS if x in organisms]:
        for rel in (range(WINDOW[o], 116) if releases == "all" else releases):
            rp = f"core_synonyms/{rel}.tsv.gz"
            dst = os.path.join(TRUTH, o, rp)
            if os.path.exists(dst):
                print("have", dst)
                continue
            db = FC.core_db(o, rel)
            if db is None:
                print(o, rel, "no core database")
                continue
            base = f"{FC.ENSEMBL}/release-{rel}/mysql/{db}"
            sd, sqlurl, sbz = FC.open_dump(f"{base}/{db}.sql.gz")
            sql = "\n".join(sd.lines(True, sbz))
            FC.add_source(o, f"{rp} (schema read from {sqlurl.rsplit('/', 1)[1]}; not kept)", sqlurl, sd.size,
                          sd.md5.hexdigest(), FC.now())
            syn = collections.defaultdict(set)
            for r in FC.dump_rows(o, rp, base, sql, "external_synonym"):
                syn[r["xref_id"]].add(r["synonym"])
            xr = {r["xref_id"]: (r["external_db_id"], r["display_label"])
                  for r in FC.dump_rows(o, rp, base, sql, "xref", keep=lambda r: r["xref_id"] in syn)}
            edb = {r["external_db_id"]: r["db_name"] for r in FC.dump_rows(o, rp, base, sql, "external_db")}
            obj = [(r["ensembl_object_type"], r["ensembl_id"], r["xref_id"])
                   for r in FC.dump_rows(o, rp, base, sql, "object_xref", keep=lambda r: r["xref_id"] in xr)]
            gene_table = "gene" if "stable_id" in (FC.columns(sql, "gene") or []) else "gene_stable_id"  # before 65
            gst = {r["gene_id"]: r["stable_id"] for r in FC.dump_rows(o, rp, base, sql, gene_table)}
            t2g = {r["transcript_id"]: r["gene_id"] for r in FC.dump_rows(o, rp, base, sql, "transcript")}
            p2t = {r["translation_id"]: r["transcript_id"] for r in FC.dump_rows(o, rp, base, sql, "translation")}
            n = 0
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with gzip.open(dst + ".tmp", "wt") as fo:
                fo.write("synonym\tgene_id\tdb_name\tlabel\tobject_type\n")
                for typ, eid, xid in obj:
                    g = eid if typ == "Gene" else t2g.get(eid) if typ == "Transcript" else t2g.get(p2t.get(eid))
                    gid = gst.get(g)
                    if not gid:
                        continue
                    dbid, lab = xr[xid]
                    for sy in sorted(syn[xid]):
                        fo.write(f"{sy}\t{gid}\t{edb.get(dbid, '')}\t{lab}\t{typ}\n")
                        n += 1
            os.replace(dst + ".tmp", dst)
            km = hashlib.md5(open(dst, "rb").read()).hexdigest()
            FC.add_source(o, f"{rp} (kept: {n:,} rows; md5 {km})", base, 0, km, FC.now())
            print("kept", dst, n, flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    rels = [115]
    if "--releases" in args:
        x = args[args.index("--releases") + 1]
        rels = "all" if x == "all" else [int(y) for y in x.split(",")]
    orgs = ORGANISMS
    if "--organisms" in args:
        orgs = args[args.index("--organisms") + 1].split(",")
        unknown = set(orgs) - set(ORGANISMS)
        if unknown:
            raise SystemExit(f"unknown organisms: {sorted(unknown)}")
    values = {args[i + 1] for i, x in enumerate(args[:-1]) if x in ("--releases", "--organisms")}
    for what in [x for x in args if not x.startswith("--") and x not in values]:
        {"hgnc": hgnc, "mgi": mgi, "ncbi": ncbi, "synonyms": lambda: synonyms(rels, orgs)}[what]()
