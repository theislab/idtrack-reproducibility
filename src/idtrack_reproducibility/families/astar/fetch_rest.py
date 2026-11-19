"""A*'s second, independent reading of the history: Ensembl REST `archive/id` (release of the service at fetch time)
for every selected id of one organism and every id on its `possible_replacement` chains (docs/methods.md, A*).

Usage: fetch_rest.py <organism> [--selection <set name>]

One connection to rest.ensembl.org, POST of up to 1000 ids (the endpoint's documented maximum), at most 2 requests per
second, `Retry-After` honoured. Records are kept in the one truth store,
`accuracy/truth/<organism>/rest/archive_astar.json` (id -> record; ids REST has no record of -> {"_error": ...}),
resumable (ids already there are not asked again), and listed in `SOURCES.tsv` with URL, size, md5 and time.
The service's release is read from `info/data` and stored in the record file's `_meta`.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

from idtrack_reproducibility.harness import harness_lib as H

from idtrack_reproducibility.families import common as C


def _data():
    return os.path.join(C.BENCH, "astar")
REST = "https://rest.ensembl.org"
BATCH = 1000


def call(path, payload=None, tries=8):
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(REST + path, data=json.dumps(payload).encode() if payload is not None else None,
                                         headers={"Content-Type": "application/json", "Accept": "application/json",
                                                  "User-Agent": "python-urllib"})
            with urllib.request.urlopen(req, timeout=600) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            err = f"{e.code} {e.read()[:200]!r}"
            time.sleep(float(e.headers.get("Retry-After") or 5 * (i + 1)) + 1)
        except Exception as e:  # noqa: BLE001
            err = repr(e)
            time.sleep(5 * (i + 1))
    raise RuntimeError(err)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--ids", default=None, help="a file of extra ids (one per line) to fetch as well")
    a = ap.parse_args()
    org = a.organism
    ids = set()
    with gzip.open(os.path.join(_data(), "selection", f"{org}.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[6] != "0":  # present in a GTF of the window
                ids.add(f[0])
    if a.ids:
        ids |= {x.strip() for x in open(a.ids) if x.strip()}
    d = os.path.join(H.truth_dir(org), "rest")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "archive_astar.json")
    arch = json.load(open(path)) if os.path.exists(path) else {}
    meta = arch.pop("_meta", {})
    info = call("/info/data")
    meta.setdefault("releases", [])
    if info.get("releases") and info["releases"] not in meta["releases"]:
        meta["releases"].append(info["releases"])
    meta.setdefault("calls", 0)
    todo = sorted(x for x in ids if x not in arch)
    rounds = 0
    while todo:
        rounds += 1
        for i in range(0, len(todo), BATCH):
            batch = todo[i:i + BATCH]
            ans = call("/archive/id", {"id": batch})
            meta["calls"] += 1
            got = {r["id"]: r for r in ans if isinstance(r, dict) and "id" in r}
            for x in batch:
                arch[x] = got.get(x, {"_error": "no record"})
            time.sleep(0.5)
            print(org, "round", rounds, i + len(batch), "of", len(todo), flush=True)
        tmp = dict(arch)
        tmp["_meta"] = meta
        json.dump(tmp, open(path + ".tmp", "w"), sort_keys=True)
        os.replace(path + ".tmp", path)
        nxt = set()
        for x, r in arch.items():
            for p in (r.get("possible_replacement") or []):
                if p["stable_id"] not in arch:
                    nxt.add(p["stable_id"])
        todo = sorted(nxt)
    meta["fetched"] = datetime.datetime.now().isoformat(timespec="seconds")
    arch["_meta"] = meta
    json.dump(arch, open(path + ".tmp", "w"), sort_keys=True)
    os.replace(path + ".tmp", path)
    md5 = hashlib.md5(open(path, "rb").read()).hexdigest()
    H.add_source(org, "rest/archive_astar.json",
                 f"{REST}/archive/id (POST, {BATCH} ids per call; A*'s {len(ids)} selected ids and their "
                 f"possible_replacement chains; service release {meta['releases']})",
                 os.path.getsize(path), md5, meta["fetched"])
    print(org, len(arch) - 1, "records", meta)


if __name__ == "__main__":
    main()
