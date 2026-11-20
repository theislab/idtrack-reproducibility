"""A*: transcripts and translations of the selected genes travelled to their gene (docs/methods.md, A*, bounded).

Usage: gen_tx.py <organism> --set <name> [--target 115] [--dry-run]

Transcript and translation ids cost 6 to 18 ms per query on d4c2dda (docs/methods.md, forty times a gene id), so the
set is bounded (the whole (id, release) matrix of transcripts is family A7's, once idtrack answers them fast):
- `astar_tx`, `astar_tl` (and `_v` with the version of the release asked from): every transcript (translation) id
  whose gene is a selected gene at some release, from its first and from its last release (GTF `tx_<release>`), to
  that release itself and to the target (115);
- `astar_tx_all`, `astar_tl_all`: for one transcript and one translation per selected gene (seeded: the smallest
  crc32 of `20261118:<id>` among those of the gene), from its first and its last release to every release of the
  window, the stable id only.
Truth: the gene the GTF of the release asked from gives the transcript (translation), then that gene's walk
(gen.Gen.truth_at: strict walk with living ids, broad, copies apart, the locus rule, REST as the second reading).
"""
import argparse
import collections
import os
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.astar import gen as GN

H = GN.H
ME = os.path.abspath(__file__)
ME_SHA256 = GN.L.sha256(ME)


def collect(org, tr, genes):
    """{(kind, id): {release: (version, gene)}} for transcripts and translations of the selected genes."""
    out = collections.defaultdict(dict)
    for r in tr.rels:
        rows = H.gtf_tx(org, r) or []
        for row in rows:
            g = row["gene_id"]
            if g not in genes:
                continue
            out[("tx", row["transcript_id"])][r] = (row["transcript_version"], g)
            if row.get("protein_id"):
                out[("tl", row["protein_id"])][r] = (row["protein_version"], g)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--selection", default=None, help="selection file name (default: the organism)")
    a = ap.parse_args()
    org, T = a.organism, a.target
    G = GN.Gen(org, T)
    sel = GN.selected(org, name=a.selection)
    unit_of = {g: u for g, u, _ in sel}
    GN.log("collect transcripts of", len(sel), "genes")
    txs = collect(org, G.tr, set(unit_of))
    GN.log("ids", collections.Counter(k for k, _ in txs))
    per_gene = {}
    for (k, x), rr in txs.items():
        for r, (v, g) in rr.items():
            key = (g, k)
            c = zlib.crc32(f"20261118:{x}".encode())
            if key not in per_gene or c < per_gene[key][0]:
                per_gene[key] = (c, x)
    chosen = {(k, x) for (g, k), (c, x) in per_gene.items()}
    rels = G.tr.rels
    walks = {}

    def walk(g, R):
        if (g, R) not in walks:
            if len(walks) > 200_000:
                walks.clear()
            walks[(g, R)] = (G.tr.walk_series(g, R, rels, True, living=True),
                             G.tr.walk_series(g, R, rels, False, living=True, resurrect=True))
        return walks[(g, R)]

    def all_lines():
        for i, ((k, x), rr) in enumerate(sorted(txs.items())):
            origins = sorted({min(rr), max(rr)})
            for R in origins:
                v, g = rr[R]
                E, B = walk(g, R)
                full = (k, x) in chosen
                targets = rels if full else sorted({R, T})
                for t in targets:
                    base = G.truth_at(g, R, t, E, B)
                    st = ("tx_" if k == "tx" else "tl_") + G.stratum(R, t)
                    forms = [(f"astar_{k}", x)]
                    if t in (R, T) and v:
                        forms.append((f"astar_{k}_v", f"{x}.{v}"))
                    for kind, q in forms:
                        kk = kind + ("_all" if full and t not in (R, T) else "")
                        yield dict(qid=f"{kk}|{q}|{R}|{t}", family=GN.FAMILY, kind=kk, stratum=st, unit=unit_of[g],
                                   id=q, from_release=R, to_release=t, strategy="both", rule="set", gene_at_origin=g,
                                   **base)
            if i % 20000 == 0:
                GN.log("ids done", i, dict(G.counts))

    man = GN.publish_stream(org, a.set, all_lines(),
                            f"A* transcripts and translations of the selected genes: every id from its first and last "
                            f"release to that release and to {T} (stable and versioned); one transcript and one "
                            f"translation per gene to every release (stable)", ["gtf/tx_<release>.tsv.gz",
                                                                                 "as gen.py"], a.dry_run,
                            generator=ME, generator_sha256=ME_SHA256)
    GN.log(man)


if __name__ == "__main__":
    main()
