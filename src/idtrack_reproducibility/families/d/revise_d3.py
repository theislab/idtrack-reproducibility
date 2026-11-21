"""Revision sets of D3 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d3.py --organism <o> --set <new set> [--revision par_shared_names] [--dry]

- `par_shared_names` (2026-11-20, `d3_r1`; rule `par_shared_names`, decided, built from idtrack's next version on): a
  marker symbol that the X and Y copies of a pseudoautosomal gene share (`CD99`) names the X gene, the Y copy reported
  apart: the accuracy harness's revision of the same name (`make_queries.rev_par_shared_names`: Ensembl's assembly
  exceptions and GTF places) applied to D3's latest lines (judge `rules`).
"""
import argparse
import os

from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--revision", default="par_shared_names", choices=["par_shared_names"])
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    sets = [s for s in L.set_names("d", a.organism) if s.startswith("d3")]
    lines = [dict(q) for q in L.load("d", a.organism, sets)[0] if q["qid"].startswith("d3|")]
    for q in lines:
        q.pop("set", None)
    out = [r for v in M.rev_par_shared_names(a.organism, lines, M.Ctx78(a.organism, H.ORGANISMS[a.organism][2])).values()
           for r in v]
    print(f"{a.organism}: {len(out)} lines revised ({a.revision})")
    if out and not a.dry:
        L.publish("d", a.organism, a.set, out, description=f"D3 revision {a.revision} (revise_d3.py docstring)",
                  generator=os.path.abspath(__file__), revises=sets, judge="rules")
        print("published", a.set)


if __name__ == "__main__":
    main()
