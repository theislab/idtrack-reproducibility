"""A* revision set for two faults of its first sets at releases whose GTF has no patch file or no versions.

1. `cmp_base_core`: a line compared without versions (`cmp` base: the GTF of the target writes no versions, human
   before 77) whose expected, broad, apart or REST genes carry versions: the genes of the core database that the GTF
   lacks (human 48-51, truth.py) were written with the core's version. Versions are dropped.
2. `copies_apart_early`: rule `judge_copies_apart`'s copies apart where the GTF labels copies `primary` (no `chr_patch_hapl_scaff`
   file before 82): a gene on a haplotype, patch or alternative-locus sequence of the target release's assembly
   exceptions (`exceptions/<release>.tsv`) goes from `expect` to `apart` when `expect` also holds a gene off them
   (the regions hunt's J3, decided for idtrack's next version; gen.py writes it so from now on).
Only the lines that change are published (F0's lines.publish), under their qids, after the latest line of each.

Usage: revise_early.py <organism> --set <new name> --revises <set>[,<set>...]   (the sets in their order)
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

ENS = re.compile(r"^(ENS[A-Z]*G\d{11})\.\d+$")
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
    exc = A.Exceptions(org)
    gtf = {}


    def copy_at(rel, gid):
        if rel not in gtf:
            rows = exc.rows(rel)
            gtf[rel] = {g: (r["region"] == "patch" or (r["seqname"] in rows and any(e["type"] != "PAR" for e in rows[r["seqname"]])))
                        for g, r in (H.gtf_genes(org, rel) or {}).items()}
        return gtf[rel].get(gid)


    def base(x):
        m = ENS.match(x)
        return m.group(1) if m else x


    out, c = [], collections.Counter()
    for s in sets[:1]:  # the revision sets after the first hold only its qids (read above as `later`)
        for q in L.iter_set(L._set_file(d, s)):
            q = later.get(q["qid"], q)
            new, why = dict(q), []
            if q.get("cmp") == "base":
                for k in ("expect", "broad", "apart", "rest"):
                    if q.get(k) and any(ENS.match(x) for x in q[k]):
                        new[k] = sorted({base(x) for x in q[k]})
                if new != q:
                    why.append("cmp_base_core")
            t = q["to_release"]
            if not q.get("final_database") and t is not None and t < 82 and len(new["expect"]) >= 2:
                side = {e: copy_at(t, base(e)) for e in new["expect"]}
                if None not in side.values() and any(side.values()) and not all(side.values()):
                    new["expect"] = sorted(e for e, v in side.items() if not v)
                    new["apart"] = sorted(set(new.get("apart") or []) | {e for e, v in side.items() if v})
                    why.append("copies_apart_early")
            if why:
                new["revised_by"] = list(q.get("revised_by") or []) + why
                new.pop("set", None)
                out.append(new)
                for w in why:
                    c[(w, q["stratum"])] += 1
    print(org, len(out), dict(c))
    if out:
        man = L.publish("astar", org, a.set, out, "cmp_base_core and copies_apart_early (revise_early.py) on " + a.revises,
                        os.path.abspath(__file__), revises=sets)
        print(man["n"], man["sha256"])


if __name__ == "__main__":
    main()
