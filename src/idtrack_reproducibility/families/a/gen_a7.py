"""Family A, stratum A7: transcripts and translations to their genes (docs/methods.md, A, the first wave).

Usage: gen_a7.py <organism> --set <name> [--target 115] [--sample 200000] [--dry-run]

Transcript and translation ids cost 6 to 18 ms per query on d4c2dda (docs/methods.md), so the first wave is the design's:
- every transcript (`tx`) and translation (`tl`) id of the GTFs (`tx_<release>`), from its first and from its last
  release R: `a7_<k>_at` the stable id to R itself (its gene at R; strategy `best`), and `a7_<k>_v` the id with its
  version at R to the target T (strategy `both`; the stable id where the GTF of R writes no versions, kind `a7_<k>`);
- `a7_<k>_s`: a seeded sample of (id, release) pairs of the whole matrix (about `--sample` lines over both kinds,
  chosen by crc32 of `20261118:<id>:<release>` below a threshold), the stable id to T (strategy `both`).
Truth (nothing of idtrack): the gene the GTF of R gives the id (a translation: the gene of its transcript), then that
gene's walk read by A*'s `gen.Gen.truth_at` (strict walk with living ids, broad, copies apart, locus), as in A1.
"""
import argparse
import json
import collections
import os
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.a import gen_a as GA

H, L, C = GA.H, GA.L, GA.C
ME = os.path.abspath(__file__)


def collect(org, tr):
    """{(kind, id): {release: (version, gene)}} over every GTF release of the window."""
    out = collections.defaultdict(dict)
    for r in tr.rels:
        for row in H.gtf_tx(org, r) or []:
            g = row["gene_id"]
            out[("tx", row["transcript_id"])][r] = (row.get("transcript_version", ""), g)
            if row.get("protein_id"):
                out[("tl", row["protein_id"])][r] = (row.get("protein_version", ""), g)
        GA.log("tx release", r, len(out))
    return out


def lines(org, T, sample=200_000, G=None, info=None):
    """Every line of A7 for one organism (a generator); `info` (a dict) gets the pairs and the threshold."""
    G = G or GA.GenA(org, T)
    tr = G.tr
    txs = collect(org, tr)
    pairs = sum(len(v) for v in txs.values())
    thr = min(1.0, sample / max(pairs, 1)) * 2 ** 32
    if info is not None:
        info.update(pairs=pairs, sample_threshold=thr)
    GA.log("ids", collections.Counter(k for k, _ in txs), "pairs", pairs)
    n = 0
    for (k, x), rr in sorted(txs.items()):
        origins = sorted({min(rr), max(rr)})
        for R in sorted(rr):
            v, g = rr[R]
            first_last = R in origins
            sampled = zlib.crc32(f"20261118:{x}:{R}".encode()) < thr
            if not (first_last or sampled):
                continue
            if g not in tr.present or R not in tr.present[g]:
                continue  # the GTF's gene of the transcript is not a gene of that GTF (never seen; guarded)
            E, B = G.walks(g, R, sorted({R, T}))
            if first_last:
                base = G.truth_at(g, R, R, E, B)
                yield G.line(f"a7_{k}_at", f"A7_{k}_at", g, x, R, R, base, gene_at_origin=g,
                             **G.axes(R, R)) | dict(strategy="best")
                base = G.truth_at(g, R, T, E, B)
                q, kind = (f"{x}.{v}", f"a7_{k}_v") if v else (x, f"a7_{k}")
                yield G.line(kind, f"A7_{k}_{T}", g, q, R, T, base, gene_at_origin=g, **G.axes(R, T))
            if sampled:
                base = G.truth_at(g, R, T, E, B)
                yield G.line(f"a7_{k}_s", f"A7_{k}_sample", g, x, R, T, base, gene_at_origin=g, **G.axes(R, T))
            n += 1
            if n % 200_000 == 0:
                GA.log("origins", n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--sample", type=int, default=200_000)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org, T = a.organism, a.target
    info = {}
    sources = ["gtf/tx_<release>.tsv.gz", "as gen_a.py"]
    man = L.publish(GA.FAMILY, org, a.set, lines(org, T, a.sample, info=info),
                    f"A A7: transcripts and translations from their first and last release to their gene there and to "
                    f"{T}, and a seeded sample of (id, release) pairs to {T}", ME, sources=sources,
                    manifest_extra=dict(stratum="A7", target=T, astar_gen_sha256=L.sha256(GA.AG.__file__),
                                        astar_truth_sha256=L.sha256(GA.TA.__file__)),
                    dry=a.dry_run, log=GA.log)
    if man and not a.dry_run:
        mp = os.path.join(L.query_dir(GA.FAMILY, org), f"{a.set}.manifest.json")
        m = json.load(open(mp))
        m.update(info)
        json.dump(m, open(mp, "w"), indent=1)
    GA.log("published", {k: v for k, v in (man or {}).items() if k != "counts"})
    GA.log("counts", (man or {}).get("counts"))


if __name__ == "__main__":
    main()
