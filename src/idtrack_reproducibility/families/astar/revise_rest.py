"""A* revision set: the second reading for lines generated while it had no record of their ids.

The lines of ids that REST had not answered when a set was generated (gen.py counted them `rest incomplete`) are
read again with the records available now (`truth.Rest`: REST's own, else the 116 core's tables read as REST reads
them, rest_from_core.py). A forward line (target after the release read) where the reading is complete and gives
other genes than the walk gets `open`, `rest` and `rest_source`, as gen.py writes them; only those lines are
published (F0's lines.publish, under their qids).

Usage: revise_rest.py <organism> --set <new name> --revises <set>[,<set>] --ids <file of ids to read again>
"""
import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.astar import truth as TA
from idtrack_reproducibility.harness import make_queries as M

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--revises", required=True)
    ap.add_argument("--ids", required=True)
    a = ap.parse_args()
    org = a.organism
    ids = {x.strip() for x in open(a.ids) if x.strip()}
    tr = TA.TruthA(org, 115)
    rest = TA.Rest(org, tr)
    out, c = [], collections.Counter()
    for s in a.revises.split(","):
        for q in L.iter_set(L._set_file(L.query_dir("astar", org), s)):
            g = M.gene_base(q["id"])
            if g not in ids or q.get("final_database") or q.get("open"):
                continue
            read = q.get("read_at") or q.get("from_release")
            t = q["to_release"]
            if read is None or t <= read:
                continue
            r, complete = rest.at(g, t)
            if not complete:
                c["still incomplete"] += 1
                continue
            e = {M.gene_base(x) for x in q["expect"]} | {M.gene_base(x) for x in q.get("apart") or []}
            if r == e:
                c["agrees"] += 1
                continue
            c["differs"] += 1
            new = dict(q)
            new.update(open="the walk of Ensembl's history and REST archive/id differ", rest=tr.vset(r, t),
                       rest_source=sorted(rest.sources(g, t)))
            new["revised_by"] = list(q.get("revised_by") or []) + ["rest_reading"]
            out.append(new)
    print(org, dict(c))
    if out:
        man = L.publish("astar", org, a.set, out, f"the second reading (REST or the 116 core) for lines of {len(ids)} ids "
                        f"read without it in {a.revises}", os.path.abspath(__file__), revises=a.revises.split(","))
        print(man["n"], man["sha256"])


if __name__ == "__main__":
    main()
