"""Cross-check the harness's reading of Ensembl's history (stable_id_event, walked by make_queries.Truth) against
Ensembl's REST service `archive/id` (an independent service and code), for every id with a saved REST record.

Usage: check_truth.py --organism <name>     writes `<work>/truth/<organism>/check_rest.json` and prints the summary

For each id not present at the target release: REST's answer is followed through `possible_replacement` until ids
that exist at the target (by the target's GTF); the harness's strict walk from the id's last release must give the
same set. Agreement is counted overall and per kind of disagreement.
"""
import argparse
import collections
import json
import os
import sys

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness.make_queries import Truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    a = ap.parse_args()
    org = a.organism
    T = H.ORGANISMS[org][2]
    tr = Truth(org, T)
    cur = set(tr.genes[T])
    arch = {}
    rd = os.path.join(H.truth_dir(org), "rest")
    for f in sorted(os.listdir(rd)):
        if f.startswith("archive") and f.endswith(".json"):
            arch.update(json.load(open(os.path.join(rd, f))))

    def t_rest(x0):
        out, seen, stack, missing = set(), set(), [x0], False
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            rec = arch.get(x)
            if x in cur:
                out.add(x)
                continue
            if rec is None or "_error" in rec:
                missing = True
                continue
            for p in rec.get("possible_replacement") or []:
                stack.append(p["stable_id"])
        return out, missing

    c = collections.Counter()
    ex = collections.defaultdict(list)
    for gid in sorted(arch):
        if gid in cur or gid not in tr.present or "_error" in arch[gid]:
            c["skipped (current, unknown to the GTFs, or REST error)"] += 1
            continue
        L = max(tr.present[gid])
        rest, missing = t_rest(gid)
        if missing:
            c["REST chain incomplete"] += 1
            continue
        mine = tr.walk_fwd(gid, L, T, True)
        broad = tr.walk_fwd(gid, L, T, False)
        if mine == rest:
            k = "equal"
        elif rest <= broad and mine <= rest:
            k = "REST larger, within the broad walk"
        elif rest <= broad:
            k = "differs, REST within the broad walk"
        elif mine < rest:
            k = "REST larger, outside the broad walk"
        else:
            k = "differs"
        c[k] += 1
        if k != "equal" and len(ex[k]) < 20:
            ex[k].append(dict(id=gid, last=L, rest=sorted(rest), walk=sorted(mine), broad=sorted(broad)))
    n = sum(v for k, v in c.items() if not k.startswith("skipped") and k != "REST chain incomplete")
    lo, hi = H.wilson(c["equal"], n)
    res = dict(organism=org, compared=n, equal=c["equal"], rate=c["equal"] / n if n else None, ci=[lo, hi], counts=dict(c),
               examples=ex)
    json.dump(res, open(os.path.join(H.truth_dir(org), "check_rest.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))
    for k, v in ex.items():
        print(k, v[:3])


if __name__ == "__main__":
    main()
