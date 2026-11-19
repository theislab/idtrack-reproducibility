"""The error bound of D1's locus truth across an assembly handoff (gen_d1.py: models mapped exon by exon through
Ensembl's chain files). Nothing of idtrack.

Usage: check_map.py <organism> <last release of the old assembly> <first release of the new one> <out tsv>

Every gene id present in the GTFs of both releases (whatever their versions: an assembly change renews most of them) is read at
the old release, its model mapped to the new assembly, and judged against the new release's models as D1 judges a
feature against 115's: the share of genes the mapping places, and of those, the share whose locus class expects the
gene itself (`one` or `dominant` with the same id), expects another gene (the truth's error), or leaves it open.
"""
import collections
import sys

from idtrack_reproducibility.families.d import gen_d1 as G

org, a, b, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
loc = G.Locus(org, target=b)
old, new = loc.source(a), loc.target
# every id of both (an assembly change gives most ids a new version: mouse 102 -> 103 keeps 35 of 52,404)
same_version = sum(1 for g in old.genes if g in new.genes and old.genes[g]["version"] == new.genes[g]["version"])
both = sorted(g for g in old.genes if g in new.genes)
c = collections.Counter()
rows = []
for g in both:
    x = loc.of(g, a)
    if x["cls"] == "unmapped":
        verdict = "unmapped"
    elif x["expect"] == [g]:
        verdict = "itself"
    elif x["expect"]:
        verdict = "another_gene"
    elif g in x["hits"]:
        verdict = "open_holding_itself"
    else:
        verdict = "open_without_itself"
    c[(x["cls"], verdict)] += 1
    if verdict in ("another_gene", "open_without_itself"):
        rows.append((g, x["cls"], verdict, ",".join(f"{k}:{v}" for k, v in x["hits"].items()), x.get("mapped") or ""))
n = len(both)
mapped = max(1, n - sum(v for (cls, verdict), v in c.items() if verdict == "unmapped"))
with open(out, "w") as fo:
    fo.write(f"# {org} {a} -> {b}: {n} genes in both (of them with the same version: {same_version}); mapped {mapped} ({mapped / n:.4f})\n")
    for (cls, verdict), v in sorted(c.items()):
        fo.write(f"# {cls}\t{verdict}\t{v}\t{v / n:.4f}\t{v / mapped:.4f} of mapped\n")
    fo.write("gene\tclass\tverdict\thits\tmapped\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(open(out).read().split("gene\tclass")[0])
