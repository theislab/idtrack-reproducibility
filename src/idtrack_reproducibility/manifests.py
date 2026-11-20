"""`idtrack-repro manifest`: the data manifests, one per area of the data root (docs/data.md, "Manifests").

The areas are declared in `manifests/areas.toml`: a name, the folders of the data root it covers (with `{organism}`
expanded), regular expressions of the files it leaves out (logs, scratch, superseded runs) and how its files are
rebuilt (`rebuild`). For every file the manifest has a row: path in the data root, bytes, sha256, source (the URL a
fetched file came from, read from the truth store's `SOURCES.tsv`; else the command or tree that made it), when it was
fetched or made, and how it is rebuilt (`fetch`, `generate`, `run` (answers of a tree: the same answers, other
seconds), or `none: ...` where a service's answer of a date cannot be asked again). A manifest of more than 5,000 rows
is written gzip-compressed (`.tsv.gz`).

Usage: idtrack-repro manifest <area>[,<area>...] | --all       (hashes every file: run large areas as jobs)
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import os
import re
import sys
import tomllib

from . import config as CFG

COLS = ["path", "bytes", "sha256", "source", "time", "rebuild"]
CHECK_HASH_BELOW = 50 * 1024 * 1024


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def areas(root):
    with open(os.path.join(root, "manifests", "areas.toml"), "rb") as fh:
        return tomllib.load(fh)["area"]


def expand(area):
    out = []
    for f in area["folders"]:
        if "{organism}" in f:
            out += [f.format(organism=o) for o in area.get("organisms", CFG.ORGANISMS)]
        else:
            out.append(f)
    return out


def sources_of(data):
    """{data-root path: (url, fetched)} from the truth store's SOURCES.tsv files (`what` is the path in the organism's
    folder, or a prefix of it)."""
    out = {}
    for org in os.listdir(os.path.join(data, "truth")) if os.path.isdir(os.path.join(data, "truth")) else []:
        p = os.path.join(data, "truth", org, "SOURCES.tsv")
        if not os.path.exists(p):
            continue
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                out[f"truth/{org}/{r['what']}"] = (r.get("url", ""), (r.get("fetched") or "").split(" ")[0])
    return out


def source_for(path, srcs):
    if path in srcs:
        return srcs[path]
    m = re.match(r"^(truth/[^/]+/gtf)/(?:genes|tx)_(\d+)\.tsv\.gz$", path)
    if m and f"{m.group(1)}/{m.group(2)}" in srcs:
        return srcs[f"{m.group(1)}/{m.group(2)}"]
    stem = re.sub(r"\.(tsv|txt|json)(\.gz)?$", "", path)
    return srcs.get(stem, ("", ""))


def area_sources(data, area):
    """{data-root path: origin} of an area that names a table of its files' origins (`sources`: columns `path` in the
    table's folder, `source`, `path_in_source`, `commit`), e.g. external inputs copied from another project."""
    if not area.get("sources"):
        return {}
    p = os.path.join(data, area["sources"])
    base = os.path.dirname(area["sources"])
    out = {}
    if os.path.exists(p):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                where = f"{r.get('source', '')}:{r.get('path_in_source', '')}" + (f" at {r['commit']}" if r.get("commit") else "")
                out[os.path.join(base, r["path"])] = (where, r.get("copied", ""))
    return out


def rebuild_of(rel, area):
    for pat, how in area.get("rules", []):
        if re.search(pat, rel):
            return how
    return area.get("rebuild", "")


def walk(data, area):
    skip = [re.compile(x) for x in area.get("exclude", [])]
    keep = [re.compile(x) for x in area.get("include", [])]
    for top in expand(area):
        base = os.path.join(data, top)
        if not os.path.exists(base):
            continue
        for d, dirs, files in os.walk(base, followlinks=True):
            dirs.sort()
            for f in sorted(files):
                p = os.path.join(d, f)
                rel = os.path.relpath(p, data)
                if any(s.search(rel) for s in skip) or (keep and not any(k.search(rel) for k in keep)):
                    continue
                yield rel, p


def build(root, name, area, dry=False):
    cfg = CFG.get()
    data = cfg.root
    srcs = dict(sources_of(data), **area_sources(data, area))
    rows = []
    for rel, p in walk(data, area):
        st = os.stat(p)
        url, fetched = source_for(rel, srcs)
        rows.append(dict(path=rel, bytes=st.st_size, sha256="" if dry else sha256(p),
                         source=url or area.get("source", ""),
                         time=fetched or datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
                         rebuild=rebuild_of(rel, area)))
    if dry:
        return rows
    path = os.path.join(root, "manifests", name + (".tsv.gz" if len(rows) > 5000 else ".tsv"))
    for old in (os.path.join(root, "manifests", name + ".tsv"), os.path.join(root, "manifests", name + ".tsv.gz")):
        if os.path.exists(old) and old != path:
            os.remove(old)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    buf = io.StringIO()
    buf.write("\t".join(COLS) + "\n")
    for r in rows:
        buf.write("\t".join(str(r[c]).replace("\t", " ") for c in COLS) + "\n")
    if path.endswith(".gz"):
        with gzip.GzipFile(path + ".tmp", "wb", mtime=0) as fo:
            fo.write(buf.getvalue().encode())
    else:
        open(path + ".tmp", "w").write(buf.getvalue())
    os.replace(path + ".tmp", path)
    return rows


def read_manifest(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def check(root, problems, skipped, full=False):
    """Every manifest's files against the data root: present, size, sha256 (files under 50 MB, or all with full)."""
    d = os.path.join(root, "manifests")
    files = [os.path.join(dp, f) for dp, _, fs in os.walk(d) for f in fs if f.endswith((".tsv", ".tsv.gz"))]
    if not files:
        problems.append(("manifests", "manifests/", "no manifest"))
        return
    data = CFG.get().root
    if not os.path.isdir(data):
        skipped.append(f"manifests: no data root at {data}")
        return
    declared = {a["name"] for a in areas(root)}
    for mf in sorted(files):
        name = os.path.relpath(mf, d).replace(".tsv.gz", "").replace(".tsv", "")
        if name not in declared:
            problems.append(("manifests", os.path.relpath(mf, root), "no area of this name in manifests/areas.toml"))
        for r in read_manifest(mf):
            p = os.path.join(data, r["path"])
            if not os.path.exists(p):
                problems.append(("manifests", os.path.relpath(mf, root), f"missing in the data root: {r['path']}"))
                continue
            if os.path.getsize(p) != int(r["bytes"]):
                problems.append(("manifests", os.path.relpath(mf, root), f"size differs: {r['path']}"))
                continue
            if (full or int(r["bytes"]) < CHECK_HASH_BELOW) and sha256(p) != r["sha256"]:
                problems.append(("manifests", os.path.relpath(mf, root), f"sha256 differs: {r['path']}"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro manifest", description=__doc__.split("\n\n")[0])
    ap.add_argument("areas", nargs="?", default=None, help="comma-separated area names (manifests/areas.toml)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry", action="store_true", help="count files and bytes only")
    a = ap.parse_args(argv)
    from .check import ROOT
    decl = {x["name"]: x for x in areas(ROOT)}
    names = list(decl) if a.all else (a.areas or "").split(",")
    for n in names:
        if n not in decl:
            raise SystemExit(f"no area {n!r} in manifests/areas.toml ({', '.join(decl)})")
        rows = build(ROOT, n, decl[n], a.dry)
        print(f"{n}: {len(rows)} files, {sum(int(r['bytes']) for r in rows) / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
