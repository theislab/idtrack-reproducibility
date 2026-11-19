"""A* revision set: places of the core database's genes that the GTF lacks (human 48-51, cynomolgus 112), as the
harness's correction `core_places` (make_queries.rev_core_places) gives them.

1. A locus line (no expected gene) whose target release has such genes gets `core_places: true` (score.py places a
   returned gene the GTF lacks from `core_places/<release>.tsv`).
2. A line left open because the query gene had no place (gen.py: "the query's place is unknown") gets its place from
   the core (`locus`; mapped across assemblies where the target is on another one) and is judged by the locus rule.
Only the lines that change are published (F0's lines.publish), under their qids, after the latest line of each.

Usage: revise_core.py <organism> --set <new name> --revises <set>[,<set>...]   (the first set and its revisions)
"""
import argparse
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import assembly_map as A
from idtrack_reproducibility.harness import score as S

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--revises", required=True)
    a = ap.parse_args()
    org = a.organism
    sets = a.revises.split(",")
    d = L.query_dir("astar", org)
    later = {}
    for s in sets[1:]:
        for q in L.iter_set(L._set_file(d, s)):
            later[q["qid"]] = q
    cp_dir = os.path.join(H.truth_dir(org), "core_places")
    have = {int(f[:-4]) for f in os.listdir(cp_dir)} if os.path.isdir(cp_dir) else set()
    places = S.Places(org)
    asm = H.gtf_assemblies(org)
    mapper = A.Mapper(org)
    ENS = re.compile(r"^(ENS[A-Z]*[GTP]\d{11})")
    UNKNOWN = "Ensembl has no successor and the query's place is unknown"
    out, c = [], collections.Counter()
    for q in L.iter_set(L._set_file(d, sets[0])):
        q = later.get(q["qid"], q)
        new, why = dict(q), []
        t = q["to_release"]
        if not q.get("final_database") and not q["expect"] and t in have and not q.get("core_places") and \
                (q.get("locus") or (q.get("open") or "").startswith(UNKNOWN)):
            new["core_places"] = True
            why.append("core_places")
        if (q.get("open") or "").startswith(UNKNOWN):
            read = q.get("read_at") or q.get("from_release")
            gid = q.get("gene_at_origin") or ENS.match(q["id"]).group(1)
            p = places.core(read, gid) if read in have else None
            if p:
                loc = dict(release=read, assembly=asm.get(read, ""), seqname=p["seqname"], start=int(p["start"]),
                           end=int(p["end"]), strand=p["strand"])
                new["locus"] = loc
                new.pop("open", None)
                to_asm = asm.get(t)
                if to_asm and to_asm != loc["assembly"]:
                    new["locus_mapped"] = mapper.map(loc, to_asm)
                new["core_places"] = True
                why.append("core_query_place")
        if why:
            new["revised_by"] = list(q.get("revised_by") or []) + sorted(set(why))
            new.pop("set", None)
            out.append(new)
            for w in why:
                c[(w, q["stratum"])] += 1
    print(org, len(out), dict(c))
    if out:
        man = L.publish("astar", org, a.set, out, "core_places (revise_core.py) on " + a.revises, os.path.abspath(__file__),
                        revises=sets, judge="78")
        print(man["n"], man["sha256"])


if __name__ == "__main__":
    main()
