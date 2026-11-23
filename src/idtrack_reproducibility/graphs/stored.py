"""A copy of a graph file with every member stored, not compressed (docs/graphs.md, "A stored copy").

Usage: idtrack-repro graph-stored <graph file> --out <folder>

A graph file is a zip archive of numpy arrays (`.npz`) whose members idtrack compresses where that makes them markedly
smaller. The copy holds the same members with the same bytes and dates, each stored (`ZIP_STORED`), under the same
name in `--out`; it loads as the same graph (the same content digest). Its use: family H's `load` measure compares
the load of the two layouts when a tree is configured a second time with the stored copies as its graphs
(docs/families/h.md). The copy is written to a temporary name and renamed when complete; the sha256 of both files is
printed. Nothing imports idtrack.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import zipfile


def stored_copy(source: str, out: str) -> str:
    """Write `<out>/<name of source>` with every member of `source` stored; return its path."""
    os.makedirs(out, exist_ok=True)
    target = os.path.join(out, os.path.basename(source))
    if os.path.realpath(target) == os.path.realpath(source):
        raise SystemExit(f"{target} is the source itself: choose another --out")
    with zipfile.ZipFile(source) as a, zipfile.ZipFile(target + ".tmp", "w", zipfile.ZIP_STORED, allowZip64=True) as b:
        for info in a.infolist():
            b.writestr(zipfile.ZipInfo(info.filename, date_time=info.date_time), a.read(info.filename),
                       compress_type=zipfile.ZIP_STORED)
    os.replace(target + ".tmp", target)
    return target


def same_members(a: str, b: str) -> bool:
    """Whether two zip archives hold the same members with the same bytes, in the same order."""
    with zipfile.ZipFile(a) as x, zipfile.ZipFile(b) as y:
        if [i.filename for i in x.infolist()] != [i.filename for i in y.infolist()]:
            return False
        return all(x.read(name) == y.read(name) for name in x.namelist())


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro graph-stored", description=__doc__.split("\n\n")[0])
    ap.add_argument("graph")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    target = stored_copy(a.graph, a.out)
    if not same_members(a.graph, target):
        raise SystemExit(f"{target}: its members differ from {a.graph}'s")
    print(f"{a.graph}\t{os.path.getsize(a.graph)} bytes\tsha256 {sha256(a.graph)}")
    print(f"{target}\t{os.path.getsize(target)} bytes\tsha256 {sha256(target)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
