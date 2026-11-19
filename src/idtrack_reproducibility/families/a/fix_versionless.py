"""A revision set for the published lines the generator no longer makes after the core-version fix (gen_a.GenA: the
core genes the GTF lacks carry no version at a release whose GTF writes none): their target is such a release and
their expected, broad or apart ids carry the core's version. The version is dropped there (the comparison of those
releases is by stable id); nothing else of the line changes.

Usage: fix_versionless.py <organism> --sets <set>,<its revisions in order> --name <revision set> [--dry-run]
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.a import gen_a as GA

L, H, C = GA.L, GA.H, GA.C
ME = os.path.abspath(__file__)


def versionless_releases(org):
    out = set()
    for r in H.gtf_releases(org):
        rows = H.gtf_genes(org, r) or {}
        if rows and not any(x.get("gene_version") for x in rows.values()):
            out.add(r)
    return out


def main():
    import pyarrow.parquet as pq
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--sets", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org, sets = a.organism, a.sets.split(",")
    vl = versionless_releases(org)
    GA.log("releases without versions", sorted(vl))
    d = L.query_dir(GA.FAMILY, org)
    latest = {}
    for s in sets:
        for b in pq.ParquetFile(L._set_file(d, s)).iter_batches(batch_size=500_000):
            for row in b.to_pylist():
                if row["to_release"] in vl:
                    latest[row["qid"]] = row
                else:
                    latest.pop(row["qid"], None)
    GA.log("lines at those targets", len(latest))

    def changed():
        for qid, row in latest.items():
            line = L.from_row(row)
            new = dict(line)
            hit = False
            for k in ("expect", "broad", "apart"):
                v = line.get(k)
                if v:
                    s = sorted({C.base_id(x) for x in v})
                    if s != sorted(v):
                        new[k] = s
                        hit = True
            if hit:
                yield new

    man = L.publish(GA.FAMILY, org, a.name, changed(), f"revision of {sets}: the core's version dropped at targets whose "
                    f"GTF writes none (fix_versionless.py)", ME, revises=sets, judge="revised", sources=[f"as {sets[0]}"],
                    fmt="parquet", dry=a.dry_run, log=GA.log)
    GA.log("revision", a.name, (man or {}).get("n"), (man or {}).get("counts"))
    if man and not a.dry_run and man.get("n") == 0:
        for f in (L._set_file(d, a.name), os.path.join(d, f"{a.name}.manifest.json")):
            if f and os.path.exists(f):
                os.remove(f)
        GA.log("no line changed: no revision set")


if __name__ == "__main__":
    main()
