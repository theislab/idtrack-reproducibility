"""The store of family G: what is fetched, where it lies, and its manifest (docs/families/g.md, "The store").

Usage: python -m idtrack_reproducibility.families.g.store [--check]

Everything family G fetched lies in `<families>/g/raw/` and is never removed, overwritten or truncated: every
response of every tool (failed, refused and partial ones under their own names), the version calls of every pass
(`<tool>/_versions/`), pyensembl's GTFs, the R tools' inputs and outputs (`_r_work/`), the pilot's and the probes'
files (`_pilot/`, `_probes/`). A pass only adds what is missing.

Writes `<families>/g/raw/STORE.tsv` (one row per stored file: area, tool, version, organism, file, bytes, sha256, the
requests that returned it, their endpoints, the first and last date asked) and `STORE_SUMMARY.tsv` (per tool, version
and organism: files, bytes, requests, usable requests, lines held, endpoints, dates, and the question sets they were
asked from with their sha256: G's published sets of the organism). `--check`: every file an index names exists with
its size and sha256 (the content-addressed names are their sha256), nothing else is written.
"""
import argparse
import collections
import hashlib
import json
import os
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.g import client as K


def raw_root():
    return os.path.join(C.family_dir("g"), "raw")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stores():
    """(tool, version, organism, folder) of every store with an index."""
    root = raw_root()
    for tool in sorted(os.listdir(root)):
        td = os.path.join(root, tool)
        if tool.startswith("_") or not os.path.isdir(td):
            continue
        for v in sorted(os.listdir(td)):
            for o in sorted(os.listdir(os.path.join(td, v))) if os.path.isdir(os.path.join(td, v)) else []:
                f = os.path.join(td, v, o)
                if os.path.exists(os.path.join(f, "index.jsonl")):
                    yield tool, v, o, f


def question_sets(organism):
    return {n: L.manifest("g", organism, n)["sha256"] for n in L.set_names("g", organism)}


def run(check=False):
    rows, summary, problems = [], [], []
    for tool, v, o, folder in stores():
        idx = [json.loads(x) for x in open(os.path.join(folder, "index.jsonl")) if x.strip()]
        by_file = collections.defaultdict(list)
        for r in idx:
            if r.get("file"):
                by_file[r["file"]].append(r)
        for f, rs in sorted(by_file.items()):
            p = os.path.join(folder, f)
            if not os.path.exists(p):
                problems.append(f"{tool}/{v}/{o}/{f}: named by the index, absent")
                continue
            size = os.path.getsize(p)
            sha = rs[-1].get("sha256") or sha256(p)
            if check and not f.endswith(".part") and sha256(p) != sha:
                problems.append(f"{tool}/{v}/{o}/{f}: sha256 differs from the index")
            dates = sorted(r["asked_at"] for r in rs if r.get("asked_at"))
            rows.append(dict(area="store", tool=tool, version=v, organism=o, file=f, bytes=size, sha256=sha,
                             requests=len(rs), endpoints=",".join(sorted({r["endpoint"] for r in rs})),
                             first_asked=dates[0] if dates else "", last_asked=dates[-1] if dates else ""))
        files = [x for x in os.listdir(folder)]
        dates = sorted(r["asked_at"] for r in idx if r.get("asked_at"))
        qs = question_sets(o) if v != "_versions" and o in CFG.ORGANISMS else {}
        summary.append(dict(tool=tool, version=v, organism=o, files=len(files),
                            bytes=sum(os.path.getsize(os.path.join(folder, x)) for x in files),
                            requests=len(idx), usable=sum(1 for r in idx if K.usable(r)),
                            lines=len({q for r in idx if K.usable(r) for q in r.get("qids") or []}),
                            endpoints=",".join(sorted({r["endpoint"] for r in idx})),
                            first_asked=dates[0] if dates else "", last_asked=dates[-1] if dates else "",
                            question_sets=",".join(f"{n}:{h[:12]}" for n, h in sorted(qs.items()))))
    for area in ("_pilot", "_probes", "_r_work"):
        top = os.path.join(raw_root(), area)
        for dp, _, fs in os.walk(top):
            for f in sorted(fs):
                p = os.path.join(dp, f)
                rows.append(dict(area=area, tool="", version="", organism="", file=os.path.relpath(p, raw_root()),
                                 bytes=os.path.getsize(p), sha256=sha256(p), requests="", endpoints="",
                                 first_asked="", last_asked=""))
    if check:
        for p in problems:
            print(p)
        return len(problems)
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(raw_root(), "STORE.tsv"), sep="\t", index=False)
    pd.DataFrame(summary).to_csv(os.path.join(raw_root(), "STORE_SUMMARY.tsv"), sep="\t", index=False)
    print(f"{len(rows)} files, {sum(r['bytes'] for r in rows):,} bytes; {len(summary)} stores; {len(problems)} problems")
    for p in problems:
        print(p)
    return len(problems)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.store", description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    sys.exit(1 if run(a.check) else 0)


if __name__ == "__main__":
    main()
