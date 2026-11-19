"""The gene ids of Ensembl's GRCh37 site at every release it was served with (76 to 115), for D2's release-mode truth.
Nothing of idtrack.

Usage: fetch_grch37.py [first last]

Ensembl keeps serving human GRCh37 beside GRCh38 (`ftp.ensembl.org/pub/grch37/release-<n>/`), and idtrack's human graph
holds that site's genes as assembly 37 at those releases; `infer_source` in release mode merges the assemblies, so a
release's gene ids are the main GTF's (GRCh38) and the GRCh37 site's. For each release: the core database's `gene`
table (`mysql/homo_sapiens_core_<n>_37/gene.txt.gz`, or `.gz.bz2`/`.bz2` at some releases; the stable id is the first column that is an Ensembl gene id, the version the next: 13 and 14 from 90 on, 14 and 15 before), its BSD `sum`
checked against the folder's CHECKSUMS, written to `accuracy/truth/homo_sapiens/grch37_site/genes_<n>.tsv.gz`
(gene_id, gene_version, biotype) with a SOURCES row; one file at a time.
"""
import datetime
import gzip
import hashlib
import io
import os
import subprocess
import sys
import tempfile

from idtrack_reproducibility.harness import harness_lib as H

ORG = "homo_sapiens"


def bsd_sum(data):
    with tempfile.NamedTemporaryFile() as fh:
        fh.write(data)
        fh.flush()
        out = subprocess.run(["sum", "-r", fh.name], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def main():
    first, last = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (76, 115)
    out_dir = os.path.join(H.truth_dir(ORG), "grch37_site")
    os.makedirs(out_dir, exist_ok=True)
    for rel in range(first, last + 1):
        out = os.path.join(out_dir, f"genes_{rel}.tsv.gz")
        if os.path.exists(out):
            print(rel, "exists", flush=True)
            continue
        import bz2
        import re
        dbs = sorted(set(re.findall(r'href="(homo_sapiens_core_\d+_37)/"', H.http_get(f"https://ftp.ensembl.org/pub/grch37/release-{rel}/mysql/",
                                                                                  timeout=120).decode())))
        if not dbs:
            raise SystemExit(f"{rel}: no homo_sapiens_core_*_37 in the release's mysql folder")
        # from 112 on the folder holds the GRCh37 core of 111 (the site's core is not rebuilt): the core it holds is used
        base = f"https://ftp.ensembl.org/pub/grch37/release-{rel}/mysql/{dbs[-1]}"
        listing = H.http_get(base + "/", timeout=120).decode()
        names = set(re.findall(r'href="([^"?/]+)"', listing))
        cks = next((n for n in ("CHECKSUMS", "CHECKSUMS.gz", "CHECKSUMS.bz2") if n in names), None)
        gname = next((n for n in ("gene.txt.gz", "gene.txt.gz.bz2", "gene.txt.bz2") if n in names), None)
        if gname is None:
            raise SystemExit(f"{rel}: no gene table in {base}/")
        sums = {}
        if cks:
            raw = H.http_get(f"{base}/{cks}", timeout=120)
            text = (gzip.decompress(raw) if cks.endswith(".gz") else bz2.decompress(raw) if cks.endswith(".bz2") else raw).decode()
            for line in text.splitlines():
                f = line.split()
                if len(f) == 3:
                    sums[f[2]] = (int(f[0]), int(f[1]))
        data = H.http_get(f"{base}/{gname}", timeout=600)
        fetched = datetime.datetime.now().isoformat(timespec="seconds")
        if gname in sums and bsd_sum(data) != sums[gname]:
            raise SystemExit(f"{rel}: {gname} sum differs from CHECKSUMS")
        if gname.endswith(".bz2"):
            data_txt = bz2.decompress(data)
            data_txt = gzip.decompress(data_txt) if gname.endswith(".gz.bz2") else data_txt
        else:
            data_txt = gzip.decompress(data)
        rows = []
        with io.TextIOWrapper(io.BytesIO(data_txt)) as fh:
            col = None
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if col is None:  # the stable id's column (12 from 90 on; 13 before, with the `status` column)
                    col = next(i for i, x in enumerate(c) if re.fullmatch(r"ENSG\d{11}", x))
                rows.append((c[col], c[col + 1], c[1]))
        H.write_tsv_gz(out, ["gene_id", "gene_version", "biotype"], sorted(rows))
        known = {r["what"] for r in H.read_tsv(os.path.join(H.truth_dir(ORG), "SOURCES.tsv"))}
        if f"grch37_site/{rel}" not in known:
            H.add_source(ORG, f"grch37_site/{rel}", f"{base}/{gname}", len(data), hashlib.md5(data).hexdigest(),
                         f"{fetched} ({cks or 'no CHECKSUMS'} sum {'checked' if gname in sums else 'not listed'})")
        print(rel, len(rows), "genes", flush=True)


if __name__ == "__main__":
    main()
