"""Family D's external inputs, copied into the data root (docs/families/d.md, "Inputs from the atlas project").

Usage: idtrack-repro family d inputs_atlas [--dry]

D1 to D3 read files of the atlas project that are not Ensembl's or an authority's: its per-feature locus tables
(`[sources] rosetta_m02`: `m02/<dataset>.tsv.gz`, and the place step's `m16/<dataset>.tsv.gz` beside it), its
datasets' preparation configs and its lineage gene table (`[sources] rosetta`: `stages/data_preparation/<dataset>/
configs/prepare.toml`, `stages/reference_data/gene_sets/lineage_genes.tsv`). This command copies them from the
configured sources into `<families>/d/inputs/atlas/` (`m02/`, `m16/`, `rosetta/<path>`) and writes
`inputs/atlas/SOURCES.tsv`: per file its path here, where it came from (the source name and the path in it), the
atlas checkout's commit, bytes, sha256 and when it was copied. A file already copied with the same sha256 is left
as it is. Afterwards the configuration's `[sources]` may name the copies (`rosetta = "<families>/d/inputs/atlas/
rosetta"`, `rosetta_m02 = "<families>/d/inputs/atlas/m02"`), so that D's sets can be generated again from the data
root alone.
"""
import argparse
import datetime
import glob
import hashlib
import os
import shutil
import subprocess
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C

ROSETTA_FILES = ["stages/data_preparation/*/configs/prepare.toml", "stages/reference_data/gene_sets/lineage_genes.tsv"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wanted():
    """[(path under inputs/atlas, source name, path in the source, absolute path)] of every file D reads."""
    cfg = CFG.get()
    m02 = cfg.source("rosetta_m02")
    rosetta = cfg.source("rosetta")
    out = []
    for sub, base in (("m02", m02), ("m16", os.path.join(os.path.dirname(m02), "m16"))):
        for p in sorted(glob.glob(os.path.join(base, "*.tsv.gz"))):
            out.append((f"{sub}/{os.path.basename(p)}", "rosetta_m02", os.path.relpath(p, os.path.dirname(m02)), p))
    for pat in ROSETTA_FILES:
        for p in sorted(glob.glob(os.path.join(rosetta, pat))):
            rel = os.path.relpath(p, rosetta)
            out.append((f"rosetta/{rel}", "rosetta", rel, p))
    return out, rosetta


def main():
    ap = argparse.ArgumentParser(prog="idtrack-repro family d inputs_atlas", description=__doc__.split("\n\n")[0])
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    files, rosetta = wanted()
    commit = subprocess.run(["git", "-C", rosetta, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dest = os.path.join(C.family_dir("d"), "inputs", "atlas")
    rows, now = [], datetime.datetime.now().isoformat(timespec="seconds")
    for rel, src, in_src, p in files:
        target = os.path.join(dest, rel)
        digest = sha256(p)
        if not a.dry and not (os.path.exists(target) and sha256(target) == digest):
            os.makedirs(os.path.dirname(target), exist_ok=True)
            shutil.copy2(p, target + ".tmp")
            os.replace(target + ".tmp", target)
        rows.append([rel, src, in_src, commit if src == "rosetta" else "", str(os.path.getsize(p)), digest, now])
    print(f"{len(rows)} files ({sum(int(r[4]) for r in rows) / 1e6:.1f} MB); atlas checkout {commit[:12]}")
    if a.dry:
        return 0
    with open(os.path.join(dest, "SOURCES.tsv.tmp"), "w") as fo:
        fo.write("path\tsource\tpath_in_source\tcommit\tbytes\tsha256\tcopied\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    os.replace(os.path.join(dest, "SOURCES.tsv.tmp"), os.path.join(dest, "SOURCES.tsv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
