"""Family D3's inputs: published marker lists, and the dates of Ensembl's releases (for the release current at a list's
date). One file at a time, one connection per host. Nothing of idtrack.

Usage: fetch_d3.py

Writes `<families>/d/inputs/` (the lists are the queries' source, not truth; the truth is the truth
store's GTF tables and HGNC/MGI files) with `SOURCES.tsv` (what, url, bytes, md5, fetched, last_modified):

- `panglaodb/PanglaoDB_markers_27_Mar_2020.tsv.gz` (PanglaoDB, its file of 2020-03-27)
- `cellmarker/all_cell_marker.zip` (CellMarker 2.0, http://bio-bigdata.hrbmu.edu.cn/CellMarker/, Download page)
- `msigdb/c8.all.v2024.1.Hs.symbols.gmt`, `msigdb/m8.all.v2024.1.Mm.symbols.gmt` (MSigDB 2024.1, cell type
  signatures, human and mouse)
- `ensembl_releases.tsv` (release, month: Ensembl's archive index, `Ensembl <n>: <month year>`), and
  `ensembl_gtf_dates.tsv` (the date of the human GTF of each release the index does not list: its build, about a month
  before the release is published)

Resumable: a file present with its SOURCES row is skipped.
"""
import datetime
import hashlib
import os
import re
import sys
import time
import urllib.request

from idtrack_reproducibility import config as _CFG
OUT = os.path.join(_CFG.get().families, "d", "inputs")
UA = {"User-Agent": "Mozilla/5.0 (compatible; benchmark-fetch)"}
FILES = [
    ("panglaodb/PanglaoDB_markers_27_Mar_2020.tsv.gz", "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"),
    ("cellmarker/all_cell_marker.zip", "http://bio-bigdata.hrbmu.edu.cn/CellMarker/file/all_cell_marker.zip"),
    ("msigdb/c8.all.v2024.1.Hs.symbols.gmt", "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c8.all.v2024.1.Hs.symbols.gmt"),
    ("msigdb/m8.all.v2024.1.Mm.symbols.gmt", "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Mm/m8.all.v2024.1.Mm.symbols.gmt"),
]


def sources():
    p = os.path.join(OUT, "SOURCES.tsv")
    if not os.path.exists(p):
        return {}
    return {line.split("\t")[0]: line for line in open(p).read().splitlines()[1:]}


def add(what, url, n, md5, lm):
    p = os.path.join(OUT, "SOURCES.tsv")
    new = not os.path.exists(p)
    with open(p, "a") as fo:
        if new:
            fo.write("what\turl\tbytes\tmd5\tfetched\tlast_modified\n")
        fo.write(f"{what}\t{url}\t{n}\t{md5}\t{datetime.datetime.now().isoformat(timespec='seconds')}\t{lm}\n")


def get(url, path):
    for i in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                data = r.read()
                lm = r.headers.get("Last-Modified", "")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path + ".part", "wb") as fo:
                fo.write(data)
            os.replace(path + ".part", path)
            return len(data), hashlib.md5(data).hexdigest(), lm
        except Exception as e:  # noqa: BLE001
            print(f"{url}: try {i + 1}: {e}", flush=True)
            time.sleep(10 * (i + 1))
    raise SystemExit(f"{url}: not fetched")


def main():
    os.makedirs(OUT, exist_ok=True)
    have = sources()
    for what, url in FILES:
        path = os.path.join(OUT, what)
        if what in have and os.path.exists(path):
            print(what, "exists")
            continue
        n, md5, lm = get(url, path)
        add(what, url, n, md5, lm)
        print(what, n, md5, lm, flush=True)
    what = "ensembl_releases.tsv"
    if what not in have:
        url = "https://www.ensembl.org/info/website/archives/index.html"
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
            html = r.read().decode()
        rows = sorted({(int(a), b) for a, b in re.findall(r"Ensembl (\d+): ([A-Z][a-z]{2} \d{4})", html)})
        if not rows:
            sys.exit("no 'Ensembl <n>: <Mon YYYY>' rows in Ensembl's archive index")
        with open(os.path.join(OUT, what), "w") as fo:
            fo.write("release\tmonth\n")
            for rel, d in rows:
                fo.write(f"{rel}\t{datetime.datetime.strptime(d, '%b %Y').strftime('%Y-%m')}\n")
        add(what, url + " (the archive index: 'Ensembl <n>: <month year>')", len(html), hashlib.md5(html.encode()).hexdigest(), "")
        print(what, len(rows), "releases")
    what = "ensembl_gtf_dates.tsv"
    if what not in have:
        known = {int(x.split("\t")[0]) for x in open(os.path.join(OUT, "ensembl_releases.tsv")).read().splitlines()[1:]}
        rows = []
        for rel in range(55, 117):
            if rel in known:
                continue
            asm = "GRCh37" if rel <= 75 else "GRCh38"
            url = f"https://ftp.ensembl.org/pub/release-{rel}/gtf/homo_sapiens/"
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                html = r.read().decode()
            m = re.search(rf"Homo_sapiens\.{asm}\.{rel}\.gtf\.gz</a></td><td align=\"right\">(\d{{4}}-\d{{2}}-\d{{2}})", html)
            rows.append((rel, m.group(1) if m else ""))
            time.sleep(1)
        with open(os.path.join(OUT, what), "w") as fo:
            fo.write("release\tgtf_date\n")
            for rel, d in rows:
                fo.write(f"{rel}\t{d}\n")
        add(what, "https://ftp.ensembl.org/pub/release-<n>/gtf/homo_sapiens/ (listing: the date of the human GTF of each release "
                  "the archive index does not list)", len(rows), "-", "")
        print(what, len(rows), "releases")


if __name__ == "__main__":
    main()
