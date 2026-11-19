"""Fetch Ensembl REST `archive/id` records for a seeded sample of retired gene ids (the cross-check of check_truth.py).

Usage: fetch_rest.py --organism <name> [--n 1500]

The sample: `--n` gene ids that exist in some GTF release but not at the target release (sorted, seeded), and then
every id named in a `possible_replacement` that is not current, until the chains are closed. POST requests of up to
200 ids, at most 4 requests per second. Writes `<work>/truth/<organism>/rest/archive_harness.json` and a SOURCES line.
"""
import argparse
import datetime
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness.make_queries import SEED, Truth

REST = "https://rest.ensembl.org/archive/id"


def post(ids, tries=6):
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(REST, data=json.dumps({"id": ids}).encode(),
                                         headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            err = f"{e.code} {e.read()[:200]!r}"
            time.sleep(float(e.headers.get("Retry-After", 5 * (i + 1))) + 1)
        except Exception as e:  # noqa: BLE001
            err = repr(e)
            time.sleep(5 * (i + 1))
    raise RuntimeError(err)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--n", type=int, default=1500)
    a = ap.parse_args()
    org, T = a.organism, H.ORGANISMS[a.organism][2]
    tr = Truth(org, T)
    cur = set(tr.genes[T])
    pop = sorted(g for g in tr.present if g not in cur)
    ids = sorted(random.Random(f"{SEED}:{org}:rest").sample(pop, min(a.n, len(pop))))
    d = os.path.join(H.truth_dir(org), "rest")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "archive_harness.json")
    arch = json.load(open(path)) if os.path.exists(path) else {}
    todo = [x for x in ids if x not in arch]
    while todo:
        for i in range(0, len(todo), 200):
            batch = todo[i:i + 200]
            ans = post(batch)
            got = {r["id"]: r for r in ans if isinstance(r, dict) and "id" in r}
            for x in batch:
                arch[x] = got.get(x, {"_error": "no record"})
            time.sleep(0.3)
        json.dump(arch, open(path + ".tmp", "w"), sort_keys=True)
        os.replace(path + ".tmp", path)
        nxt = set()
        for r in arch.values():
            for p in r.get("possible_replacement") or []:
                if p["stable_id"] not in arch and p["stable_id"] not in cur:
                    nxt.add(p["stable_id"])
        todo = sorted(nxt)
    H.add_source(org, "rest/archive_harness.json", f"{REST} (POST; {len(ids)} sampled ids and their chains)",
                 os.path.getsize(path), "-", datetime.datetime.now().isoformat(timespec="seconds"))
    print(org, len(arch), "records")


if __name__ == "__main__":
    main()
