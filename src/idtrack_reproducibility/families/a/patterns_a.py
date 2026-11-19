"""Family A: the pattern of every gene id of the window, from Ensembl's tables only (the strata of the hunts' findings
measured over family A's whole sets, docs/methods.md principle 8).

Usage: patterns_a.py <organism> [--target 115]   -> <families>/a/patterns/<organism>.tsv.gz (id, unit, history_class, flags)

Flags (an id may have several; none = `plain`):
- `cross_events`: the id has an event with another id (a lineage); `lineage5`: its lineage has 5 ids or more;
- `retired_none` (hh-4): retired, and the strict walk from its last release reaches nothing at the target (Ensembl's
  history gives no successor; where idtrack answers through a bridge, rule `bridge_same_place`);
- `split`, `merged`, `renamed`, `retired_branch`: the travel matrix's history class of a retired id;
- `gap` (hh-3): absent from a GTF release between two where it is present (retired and back, or a hole);
- `version_drop` (hh-3): its version falls between two consecutive GTF releases; `version_jump`: rises by more than 1;
- `moved_seq` (hh-2): on another chromosome at another release of the same assembly (alternate-haplotype and patch
  sequences left out, as A* does); `strand_flip`: the other strand at another release of the same assembly;
- `patch` (R6, J3, J4): on a patch or alternative locus at some release (the GTF marks them for human and mouse);
- `mt` (R1): on the mitochondrial sequence at some release;
- `core_only`: in the core database but not the GTF at some release (hh-5, J2);
- `handoff_dies`, `handoff_born`: present on one side of an assembly handoff only;
- `astar`: selected by family A* (`<families>/astar/selection/<organism>.tsv.gz`), so A's average and A*'s ids can be read apart.
"""
import argparse
import collections
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.a import gen_a as GA

H, C = GA.H, GA.C
MT = {"MT", "chrM", "M", "MtDNA", "mitochondrion_genome"}


def is_chrom(name):
    return bool(name) and len(name) <= 6 and not any(x in name.upper() for x in ("HSCHR", "HAP", "PATCH", "CTG", "_"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--target", type=int, default=115)
    a = ap.parse_args()
    org = a.organism
    G = GA.GenA(org, a.target)
    tr = G.tr
    comp = collections.Counter(G.unit.values())
    astar = set()
    p = os.path.join(C.BENCH, "astar", "selection", f"{org}.tsv.gz")
    if os.path.exists(p):
        astar = {r["id"] for r in H.read_tsv(p)}
    eras = tr.eras()
    out = os.path.join(C.BENCH, "a", "patterns")
    os.makedirs(out, exist_ok=True)
    n = collections.Counter()
    with gzip.open(os.path.join(out, f"{org}.tsv.gz.tmp"), "wt") as fo:
        fo.write("id\tunit\thistory_class\tflags\n")
        for g in sorted(tr.present):
            rr = sorted(tr.present[g])
            fl = set()
            if g in G.cross:
                fl.add("cross_events")
                if comp[G.unit[g]] >= 5:
                    fl.add("lineage5")
            hc = G.history_class(g)
            if hc in ("retired_none", "split", "merged", "renamed", "retired_branch"):
                fl.add(hc)
            if not G.contiguous(g):
                fl.add("gap")
            vers = [(r, tr.ver(g, r)) for r in rr if tr.ver(g, r)]
            for (r1, v1), (r2, v2) in zip(vers, vers[1:]):
                if G.idx[r2] == G.idx[r1] + 1:
                    if int(v2) < int(v1):
                        fl.add("version_drop")
                    elif int(v2) > int(v1) + 1:
                        fl.add("version_jump")
            places = collections.defaultdict(set)
            for r in rr:
                x = tr.genes[r].get(g) or {}
                if x.get("region") == "patch":
                    fl.add("patch")
                if x.get("region") == "core_only":
                    fl.add("core_only")
                if x.get("seqname") in MT:
                    fl.add("mt")
                if x.get("region") == "primary" and is_chrom(x.get("seqname", "")):
                    places[tr.assembly.get(r, "")].add((x["seqname"], x["strand"]))
            for asm, pl in places.items():
                if len({s for s, _ in pl}) > 1:
                    fl.add("moved_seq")
                for s in {s for s, _ in pl}:
                    if len({st for s2, st in pl if s2 == s}) > 1:
                        fl.add("strand_flip")
            for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
                if hi1 in tr.present[g] and lo2 not in tr.present[g]:
                    fl.add("handoff_dies")
                if lo2 in tr.present[g] and hi1 not in tr.present[g]:
                    fl.add("handoff_born")
            if g in astar:
                fl.add("astar")
            for f in fl or {"plain"}:
                n[f] += 1
            fo.write(f"{g}\t{G.unit_of(g)}\t{hc}\t{','.join(sorted(fl)) or 'plain'}\n")
    os.replace(os.path.join(out, f"{org}.tsv.gz.tmp"), os.path.join(out, f"{org}.tsv.gz"))
    GA.log("patterns", org, dict(sorted(n.items())))


if __name__ == "__main__":
    main()
