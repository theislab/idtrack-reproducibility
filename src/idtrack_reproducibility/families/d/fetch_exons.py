"""Gene models (merged exons) of Ensembl's GTF of given releases, for family D1's locus truth. Imports nothing of idtrack.

Usage: fetch_exons.py <organism>:<release>[,<release>...] ...

For every (organism, release): the main GTF of the release (`<Species>.<Assembly>.<release>.gtf.gz`, the file Cell
Ranger and rosetta's references are made from: chromosomes and unplaced scaffolds, no patches) is fetched over HTTPS
from Ensembl's FTP, one file at a time (one connection to the host), its size and md5 taken and its BSD `sum` checked
against the release folder's CHECKSUMS; then written to the one truth store:

  accuracy/truth/<organism>/exons/<release>.tsv.gz   gene_id, gene_version, gene_name, gene_biotype, seqname, strand,
                                                     start, end (the exons' span), bases (exon bases, merged),
                                                     exons ("s-e,s-e,..": merged exons, 1-based inclusive), places
                                                     (sequence:strand pairs the exons lie on; a gene on several has
                                                     no one model)
  accuracy/truth/<organism>/SOURCES.tsv              one row `exons/<release>` (url, bytes, md5, time)

The GTF itself is not kept (its row in SOURCES.tsv names it). Resumable: a release whose table exists is skipped.
"""
import datetime
import gzip
import hashlib
import os
import re
import subprocess
import sys
import time
import urllib.request

from idtrack_reproducibility.harness import harness_lib as H

FTP = "https://ftp.ensembl.org/pub"
ATTR = re.compile(r'(\S+) "([^"]*)"')
SPECIES = {"homo_sapiens": "Homo_sapiens", "mus_musculus": "Mus_musculus", "sus_scrofa": "Sus_scrofa",
           "macaca_mulatta": "Macaca_mulatta", "macaca_fascicularis": "Macaca_fascicularis"}
from idtrack_reproducibility import config as _CFG
TMP = os.path.join(_CFG.get().cache, "d", "tmp")


def assembly(org, rel):
    """The assembly of a release's main GTF, from the harness's SOURCES rows of the GTF it read."""
    a = H.gtf_assemblies(org).get(rel)
    if not a:
        raise SystemExit(f"{org} {rel}: no GTF row in SOURCES.tsv to name its assembly")
    return a


def fetch(url, path):
    for i in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
            md5, n = hashlib.md5(), 0
            with urllib.request.urlopen(req, timeout=600) as r, open(path + ".part", "wb") as fo:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    md5.update(b)
                    n += len(b)
                    fo.write(b)
            os.replace(path + ".part", path)
            return n, md5.hexdigest()
        except Exception as e:  # noqa: BLE001
            print(f"  {url}: try {i + 1} failed: {e}", flush=True)
            time.sleep(10 * (i + 1))
    raise SystemExit(f"{url}: not fetched")


def bsd_sum(path):
    out = subprocess.run(["sum", "-r", path], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def models(path):
    genes = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            a = dict(ATTR.findall(p[8]))
            gid = a.get("gene_id")
            if not gid:
                continue
            g = genes.get(gid)
            if g is None:
                g = genes[gid] = dict(version=a.get("gene_version", ""), name=a.get("gene_name", ""),
                                      biotype=a.get("gene_biotype") or a.get("gene_type") or p[1], places={}, )
            g["places"].setdefault((p[0], p[6]), []).append((int(p[3]), int(p[4])))
    rows = []
    for gid in sorted(genes):
        g = genes[gid]
        places = sorted(g["places"])
        (seq, strand) = max(places, key=lambda k: sum(e - s + 1 for s, e in g["places"][k]))
        ex = merge(g["places"][(seq, strand)])
        rows.append([gid, g["version"], g["name"], g["biotype"], seq, strand, ex[0][0], ex[-1][1],
                     sum(e - s + 1 for s, e in ex), ",".join(f"{s}-{e}" for s, e in ex),
                     ";".join(f"{a}:{b}" for a, b in places)])
    return rows


def one(org, rel):
    out_dir = os.path.join(H.truth_dir(org), "exons")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{rel}.tsv.gz")
    if os.path.exists(out):
        print(f"{org} {rel}: exists", flush=True)
        return
    asm = assembly(org, rel)
    name = f"{SPECIES[org]}.{asm}.{rel}.gtf.gz"
    base = f"{FTP}/release-{rel}/gtf/{org}"
    sums = {}
    try:
        text = H.http_get(f"{base}/CHECKSUMS", tries=2, timeout=120).decode()
    except Exception:  # noqa: BLE001  (releases before 69 have no CHECKSUMS file)
        text = ""
    for line in text.splitlines():
        f = line.split()
        if len(f) == 3:
            sums[f[2]] = (int(f[0]), int(f[1]))
    known = {row["url"]: row for row in H.read_tsv(os.path.join(H.truth_dir(org), "SOURCES.tsv"))
             if row["what"] == f"gtf/{rel}"}
    url = f"{base}/{name}"
    if name not in sums and url not in known:
        raise SystemExit(f"{name}: neither {base}/CHECKSUMS nor the harness's SOURCES.tsv gives a checksum")
    os.makedirs(TMP, exist_ok=True)
    local = os.path.join(TMP, name)
    t0 = time.time()
    n, md5 = fetch(url, local)
    fetched = datetime.datetime.now().isoformat(timespec="seconds")
    if name in sums:
        got = bsd_sum(local)
        if got != sums[name]:
            os.remove(local)
            raise SystemExit(f"{name}: sum {got} != CHECKSUMS {sums[name]}")
        check = f"CHECKSUMS sum {got[0]} {got[1]} checked"
    else:
        if md5 != known[url]["md5"] or n != int(known[url]["bytes"]):
            os.remove(local)
            raise SystemExit(f"{name}: md5 {md5} / {n} bytes differ from the harness's fetch ({known[url]['md5']}, {known[url]['bytes']})")
        check = f"no CHECKSUMS at this release; md5 and size equal the harness's fetch of {known[url]['fetched']}"
    rows = models(local)
    H.write_tsv_gz(out, ["gene_id", "gene_version", "gene_name", "gene_biotype", "seqname", "strand", "start", "end", "bases",
                         "exons", "places"], rows)
    H.add_source(org, f"exons/{rel}", url, n, md5, f"{fetched} ({check})")
    os.remove(local)
    print(f"{org} {rel}: {name} {n} bytes, {len(rows)} genes, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        org, rels = arg.split(":")
        for r in rels.split(","):
            one(org, int(r))
