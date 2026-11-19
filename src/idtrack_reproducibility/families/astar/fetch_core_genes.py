"""Fetch the core database's gene stable ids of the releases whose GTF lacks genes that the core database has.

Ensembl's GTFs of human 48-51 leave out genes on unplaced scaffolds, and the GTF of cynomolgus 112 leaves out 47 genes
(measured by the weakness hunt `cleanup_hunt_history`, `mine/<organism>.gtf_vs_core.tsv`, on the tables idtrack's
builder downloaded). The graph is built from the core tables, so the history truth reads presence from both. This
script takes those genes from Ensembl's own FTP MySQL dumps (not from idtrack's downloads): `gene_stable_id.txt.gz`
(releases before 65) or `gene.txt.gz` (stable_id and version in the gene table), into the one truth store
(`accuracy/truth/<organism>/core/<core db>.<table>.txt.gz`, a line in SOURCES.tsv each).

Usage: fetch_core_genes.py            (the releases of RELEASES; one connection to ftp.ensembl.org)
Reader: `core_only(org, gtf_genes_of_release)` in truth.py.
"""
import datetime
import hashlib
import os
import sys

from idtrack_reproducibility.harness import harness_lib as H

FTP = "https://ftp.ensembl.org/pub"
RELEASES = {  # organism -> {release: (core db, table)}
    "homo_sapiens": {48: ("homo_sapiens_core_48_36j", "gene_stable_id"), 49: ("homo_sapiens_core_49_36k", "gene_stable_id"),
                     50: ("homo_sapiens_core_50_36l", "gene_stable_id"), 51: ("homo_sapiens_core_51_36m", "gene_stable_id")},
    "macaca_fascicularis": {112: ("macaca_fascicularis_core_112_6", "gene")},
}


def main():
    for org, rr in RELEASES.items():
        for rel, (db, table) in sorted(rr.items()):
            path = os.path.join(H.truth_dir(org), "core", f"{db}.{table}.txt.gz")
            if os.path.exists(path):
                continue
            url = f"{FTP}/release-{rel}/mysql/{db}/{table}.txt.gz"
            data = H.http_get(url)
            with open(path + ".tmp", "wb") as fo:
                fo.write(data)
            os.replace(path + ".tmp", path)
            H.add_source(org, f"core/{db}.{table}.txt.gz", url, len(data), hashlib.md5(data).hexdigest(),
                         datetime.datetime.now().isoformat(timespec="seconds") + " (bench A*: genes of the core absent "
                         "from the GTF)")
            print(org, rel, url, len(data), flush=True)


if __name__ == "__main__":
    main()
