"""Checks of the travel matrix that need no truth, on the stored answers of a scored run (and the queries they need).

Usage: travel_checks.py --clone <clone> --graph <graph file> --organism <name> --run <scored travel run> [--n 3000]

1. identity: an id asked at a release where it exists, to that release, answers itself (every `travel_*` line of an
   id with from = to, `best` and `all`);
2. composition: where travelling R1 -> R2 answers one gene h (`best`), h travelled R2 -> R3 gives what R1 -> R3 gives
   (a seeded sample of up to --n triples, R3 another target of the same id and origin);
3. return: where R1 -> R2 answers the id itself (`best`), the id travelled R2 -> R1 answers itself (up to --n pairs;
   the stored answer when the matrix holds the query, else asked).
Queries the matrix does not hold are asked through the clone (`API.convert_identifier`, as score.py does). Every
violation is asked again with `explain=True` (the paths, cut to 20). Writes into <run>/checks/: `summary.json`,
`violations.jsonl` (check, the queries, the answers, the paths), `asked.jsonl` (the extra queries and answers).
"""
import argparse
import collections
import gc
import glob
import json
import os
import random
import re
import sys
import time

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import score

ENS = re.compile(r"^(ENS[A-Z]*[GTPE]\d{11})(\.\d+)?$")


def base(x):
    m = ENS.match(str(x))
    return m.group(1) if m else str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--release", type=int, default=115)
    a = ap.parse_args()
    out = os.path.join(a.run, "checks")
    os.makedirs(out, exist_ok=True)
    lines, _ = score.load_queries(a.organism, None, None)
    travel = {r["qid"]: r for r in lines if r["kind"].startswith("travel_") and r["kind"] != "travel_name"}
    res = {}
    for p in glob.glob(os.path.join(a.run, "results", "shard_*.jsonl")):
        for line in open(p):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            res[r["qid"]] = r
    stored = {}  # (id, from, to, strategy) -> answer (gene answers only)
    for qid, rec in travel.items():
        r = res.get(qid)
        if rec["final_database"] is None and r and "got" in r:
            stored[(rec["id"], rec["from_release"], rec["to_release"], rec["strategy"])] = r["got"]
    rg = random.Random(20261118)

    sys.path.insert(0, a.clone)
    import idtrack  # noqa: E402
    repo = os.path.join(out, "repo")
    os.makedirs(repo, exist_ok=True)
    link = os.path.join(repo, os.path.basename(a.graph))
    if not os.path.exists(link):
        os.symlink(os.path.abspath(a.graph), link)
    api = idtrack.API(local_repository=repo)
    api.build_graph(organism_name=a.organism, snapshot_release=a.release)
    gc.collect()
    gc.freeze()
    asked = {}
    fa = open(os.path.join(out, "asked.jsonl"), "w")

    def ask(ident, fr, to, strategy="best"):
        key = (ident, fr, to, strategy)
        if key in stored:
            return stored[key]
        if key not in asked:
            t = time.time()
            try:
                r = api.convert_identifier(ident, from_release=fr, to_release=to, strategy=strategy)
                asked[key] = sorted(set(score.targets_of(r, {"final_database": None})))
            except Exception as e:  # noqa: BLE001
                asked[key] = [f"error: {type(e).__name__}: {e}"[:200]]
            fa.write(json.dumps(dict(id=ident, from_release=fr, to_release=to, strategy=strategy, got=asked[key],
                                     seconds=round(time.time() - t, 4))) + "\n")
        return asked[key]

    viol = []
    summary = {}
    # 1. identity
    n = bad = 0
    for (ident, fr, to, st), got in sorted(stored.items(), key=str):
        if fr != to:
            continue
        n += 1
        if {base(x) for x in got} != {base(ident)}:
            bad += 1
            viol.append(dict(check="identity", queries=[[ident, fr, to, st]], answers=[got]))
    summary["identity"] = dict(checked=n, violations=bad)
    # 2. composition
    by_origin = collections.defaultdict(dict)
    for (ident, fr, to, st), got in stored.items():
        if st == "best":
            by_origin[(ident, fr)][to] = got
    triples = []
    for (ident, fr), tg in sorted(by_origin.items()):
        singles = [t for t, got in tg.items() if len({base(x) for x in got}) == 1 and t != fr]
        for t2 in singles:
            for t3 in tg:
                if t3 not in (t2,):
                    triples.append((ident, fr, t2, t3))
    triples = sorted(rg.sample(triples, min(a.n, len(triples)))) if triples else []
    n = bad = 0
    for ident, fr, t2, t3 in triples:
        h = base(by_origin[(ident, fr)][t2][0])
        direct = {base(x) for x in by_origin[(ident, fr)][t3]}
        via = ask(h, t2, t3)
        n += 1
        if {base(x) for x in via} != direct:
            bad += 1
            viol.append(dict(check="composition", queries=[[ident, fr, t2, "best"], [h, t2, t3, "best"],
                                                           [ident, fr, t3, "best"]],
                             answers=[by_origin[(ident, fr)][t2], via, by_origin[(ident, fr)][t3]]))
    summary["composition"] = dict(checked=n, violations=bad)
    # 3. return
    pairs = []
    for (ident, fr), tg in sorted(by_origin.items()):
        for t, got in tg.items():
            if t != fr and {base(x) for x in got} == {base(ident)}:
                pairs.append((ident, fr, t))
    pairs = sorted(rg.sample(pairs, min(a.n, len(pairs)))) if pairs else []
    n = bad = 0
    for ident, fr, t in pairs:
        back = ask(ident, t, fr)
        n += 1
        if {base(x) for x in back} != {base(ident)}:
            bad += 1
            viol.append(dict(check="return", queries=[[ident, fr, t, "best"], [ident, t, fr, "best"]],
                             answers=[by_origin[(ident, fr)][t], back]))
    summary["return"] = dict(checked=n, violations=bad)
    fa.close()
    # every violation's queries asked again with explain=True
    for v in viol:
        v["paths"] = []
        for ident, fr, to, st in v["queries"]:
            try:
                r = api.convert_identifier(ident, from_release=fr, to_release=to, strategy=st, explain=True)
                paths = r.get("the_path") or []
                v["paths"].append(dict(query=[ident, fr, to, st], match=r.get("match"),
                                       paths=[{k: p.get(k) for k in ("ensembl_gene", "edges", "external_jump")}
                                              for p in paths[:20]], n_paths=len(paths)))
            except Exception as e:  # noqa: BLE001
                v["paths"].append(dict(query=[ident, fr, to, st], error=f"{type(e).__name__}: {e}"[:300]))
    with open(os.path.join(out, "violations.jsonl"), "w") as fo:
        for v in viol:
            fo.write(json.dumps(v, default=str) + "\n")
    summary["asked"] = len(asked)
    json.dump(summary, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
