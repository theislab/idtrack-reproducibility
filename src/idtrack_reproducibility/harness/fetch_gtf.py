"""Fetch Ensembl's GTF files over HTTPS and keep the gene, transcript and translation ids with versions and places.

Usage: fetch_gtf.py --organism <name> [--releases 48-116] [--procs 6]

For every release whose FTP folder has a GTF for the organism, writes under
`<work>/truth/<organism>/gtf/`:

- `genes_<release>.tsv.gz`: gene_id, gene_version, gene_name, gene_source, gene_biotype, seqname, start, end, strand, region
  (`primary` = the main GTF file; `patch` = only in the `chr_patch_hapl_scaff` file: alternative loci, patches)
- `tx_<release>.tsv.gz`: transcript_id, transcript_version, gene_id, protein_id, protein_version, region
- one line per downloaded file in `<work>/truth/<organism>/SOURCES.tsv` (url, bytes, md5, fetched at)

Nothing here reads idtrack. Resumable: a release whose two output files exist is skipped.
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
import urllib.request
from multiprocessing import Pool

from idtrack_reproducibility.harness import harness_lib as H

FTP = "https://ftp.ensembl.org/pub"
ATTR = re.compile(r'(\S+) "([^"]*)"')


def listing(url):
    try:
        html = H.http_get(url, tries=2, timeout=120).decode("utf-8", "replace")
    except Exception:
        return None
    return re.findall(r'href="([^"?/][^"]*)"', html)


def parse(data, region, genes, txs):
    with gzip.open(io.BytesIO(data), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feat = p[2]
            if feat not in ("gene", "transcript", "CDS", "exon"):
                continue
            a = dict(ATTR.findall(p[8]))
            gid = a.get("gene_id")
            if not gid:
                continue
            if gid not in genes:
                biot = a.get("gene_biotype") or a.get("gene_type") or (p[1] if feat == "exon" else "")
                genes[gid] = [a.get("gene_version", ""), a.get("gene_name", ""), a.get("gene_source", ""), biot,
                              p[0], p[3], p[4], p[6], region]
            g = genes[gid]
            if feat == "gene":  # old GTFs have no gene lines: the place is then the span of the exons
                g[5], g[6] = p[3], p[4]
            elif feat == "exon" and g[8] == region:
                g[5] = str(min(int(g[5]), int(p[3])))
                g[6] = str(max(int(g[6]), int(p[4])))
            tid = a.get("transcript_id")
            if tid:
                t = txs.get(tid)
                if t is None:
                    t = txs[tid] = [a.get("transcript_version", ""), gid, "", "", region]
                if feat == "CDS" and a.get("protein_id"):
                    t[2], t[3] = a["protein_id"], a.get("protein_version", "")


def one_release(args):
    org, rel, out = args
    gpath, tpath = os.path.join(out, f"genes_{rel}.tsv.gz"), os.path.join(out, f"tx_{rel}.tsv.gz")
    if os.path.exists(gpath) and os.path.exists(tpath):
        return rel, "exists", []
    files = listing(f"{FTP}/release-{rel}/gtf/{org}/")
    if not files:
        return rel, "no gtf folder", []
    main = [f for f in files if re.search(rf"\.{rel}\.gtf\.gz$", f) and not re.search(r"\.(abinitio|chr|chr_patch_hapl_scaff)\.gtf", f)]
    patch = [f for f in files if f.endswith(f".{rel}.chr_patch_hapl_scaff.gtf.gz")]
    if len(main) != 1:
        return rel, f"unexpected listing {files}", []
    genes, txs, src = {}, {}, []
    for fname, region in [(main[0], "primary")] + [(f, "patch") for f in patch]:
        url = f"{FTP}/release-{rel}/gtf/{org}/{fname}"
        data = H.http_get(url, timeout=900)
        src.append((url, len(data), hashlib.md5(data).hexdigest(), datetime.datetime.now().isoformat(timespec="seconds")))
        parse(data, region, genes, txs)
    H.write_tsv_gz(gpath, ["gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype", "seqname", "start", "end",
                           "strand", "region"], ([k] + v for k, v in sorted(genes.items())))
    H.write_tsv_gz(tpath, ["transcript_id", "transcript_version", "gene_id", "protein_id", "protein_version", "region"],
                   ([k] + v for k, v in sorted(txs.items())))
    return rel, f"{len(genes)} genes, {len(txs)} transcripts", src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--releases", default="48-116")
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.releases.split("-"))
    out = os.path.join(H.truth_dir(a.organism), "gtf")
    os.makedirs(out, exist_ok=True)
    with Pool(a.procs) as pool:
        for rel, msg, src in pool.imap_unordered(one_release, [(a.organism, r, out) for r in range(hi, lo - 1, -1)]):
            for s in src:
                H.add_source(a.organism, f"gtf/{rel}", *s)
            print(rel, msg, flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
