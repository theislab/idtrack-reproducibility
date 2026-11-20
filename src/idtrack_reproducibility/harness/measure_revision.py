"""Measure a correction of the judge before it is published: the lines it would revise, and a tree's stored answers to
them scored under the published truth and under the revised one (nothing is written into the query sets).

Usage: python -m idtrack_reproducibility.harness.measure_revision --organism <o> --correction <name> --answers <folder> [...]
           [--out <json>]

`--answers`: stored answer folders of the harness (`<harness>/answers/<tree>/<organism>/{kinds,travel}`, their
`results/shard_*.jsonl`). The correction reads the published lines as the corrections before it revised them
(`make_queries.not_revised_after`). Prints and writes: the lines revised per kind, the moves of the outcome per kind
(`score.score_one`, the harness's own rule), and examples of each move.
"""
import argparse
import collections
import glob
import json
import os

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as MQ
from idtrack_reproducibility.harness import score as S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--correction", required=True, choices=list(MQ.REVISIONS))
    ap.add_argument("--answers", nargs="+", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    lines, _ = MQ.published_lines(H.query_dir(a.organism))
    C = MQ.Ctx78(a.organism, H.ORGANISMS[a.organism][2])
    rev = MQ.REVISIONS[a.correction](a.organism, MQ.not_revised_after(a.correction, lines.values()), C)
    new = {r["qid"]: r for v in rev.values() for r in v}
    stored = {}
    for folder in a.answers:
        for p in sorted(glob.glob(os.path.join(folder, "results", "shard_*.jsonl"))):
            with open(p) as fh:
                for raw in fh:
                    try:
                        r = json.loads(raw)
                    except json.JSONDecodeError:
                        continue
                    if r.get("qid") in new:
                        stored[r["qid"]] = r
    places = S.Places(a.organism)
    moves, by_kind, examples = collections.Counter(), collections.Counter(), collections.defaultdict(list)
    for qid, rec in new.items():
        res = stored.get(qid)
        if res is None:
            moves["no stored answer"] += 1
            continue
        x, y = S.score_one(lines[qid], res, places)[0], S.score_one(rec, res, places)[0]
        moves[f"{x} -> {y}"] += 1
        by_kind[f"{rec['kind']}|{x} -> {y}"] += 1
        if x != y and len(examples[f"{x} -> {y}"]) < 4:
            examples[f"{x} -> {y}"].append(dict(qid=qid, expect=lines[qid]["expect"], revised=rec["expect"],
                                                apart=rec.get("apart"), got=res.get("got"), got_apart=res.get("apart")))
    out = dict(organism=a.organism, correction=a.correction, lines_revised={k: len(v) for k, v in sorted(rev.items())},
               moves=dict(moves), by_kind=dict(sorted(by_kind.items())), examples=examples,
               answers=[os.path.relpath(os.path.abspath(f), CFG.get().root) for f in a.answers])
    print(json.dumps({k: v for k, v in out.items() if k != "examples"}, indent=1))
    if a.out:
        json.dump(out, open(a.out, "w"), indent=1, default=str)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
