"""Probe A0's rows: print what the graph holds for given nodes (history edges, cross-reference edges in and out, and
for genes the names through transcripts and translations), to read a mismatch class by hand.

Usage: probe_a0.py <organism> --clone <clone> --graph <graph> --work <folder> <node or id> ...
(an id without version: every node of it)
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.a import a0_fidelity as A0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("nodes", nargs="+")
    a = ap.parse_args()
    gr = A0.GraphReader(a.clone, a.graph, a.organism, a.work)
    gr.tx, gr.tl, gr.backbone = set(gr.tx), set(gr.tl), set(gr.backbone)
    G = gr.G
    idx = G._frozen().index
    for q in a.nodes:
        names = sorted(n for n in gr.backbone | gr.tx | gr.tl if n.startswith(q + ".")) or ([q] if q in idx else [])
        print("=" * 20, q, names, flush=True)
        for n in names:
            print("--", n)
            try:
                print("  active", sorted(G.active_releases(n)))
            except Exception as e:  # noqa: BLE001
                print("  active error", e)
            for rev in (False, True):
                try:
                    print("  history", "in" if rev else "out", G.history_edges(n, reverse=rev))
                except Exception as e:  # noqa: BLE001
                    print("  history error", e)
                try:
                    for src, con in G.xref_edges(n, reverse=rev):
                        print("  xref", "in" if rev else "out", src, con)
                        if rev and (src in gr.tx or src in gr.tl) and (n in gr.backbone or n in gr.tx):
                            for s2, c2 in G.xref_edges(src, reverse=True):
                                print("     <-", s2, c2)
                                if s2 in gr.tl:
                                    for s3, c3 in G.xref_edges(s2, reverse=True):
                                        print("        <-", s3, c3)
                except Exception as e:  # noqa: BLE001
                    print("  xref error", e)


if __name__ == "__main__":
    main()
