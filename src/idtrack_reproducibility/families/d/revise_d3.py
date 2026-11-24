"""Revision sets of D3 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d3.py --organism <o> --set <new set> [--revision par_shared_names|names_broad] [--dry]

- `par_shared_names` (2026-11-20, `d3_r1`; rule `par_shared_names`, decided, built from idtrack's next version on): a
  marker symbol that the X and Y copies of a pseudoautosomal gene share (`CD99`) names the X gene, the Y copy reported
  apart: the accuracy harness's revision of the same name (`make_queries.rev_par_shared_names`: Ensembl's assembly
  exceptions and GTF places) applied to D3's latest lines (judge `rules`).
- `names_broad` (the corrected judge's `names_broad`, judge `78`): a symbol that other genes carry as a name in another
  name source at the release read (`read_at`) puts the genes those carriers become at 115 into `broad`, or, on a line
  open within the genes that carry it at 115, within them: D1's correction of that name (`revise_d1.names_broad`, its
  docstring) on D3's lines.
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
    ap.add_argument("--revision", default="par_shared_names", choices=["par_shared_names", "names_broad"])
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--lines-to", default=None, help="with --dry: the revised lines, for reading before publishing")
    a = ap.parse_args()
    sets = [s for s in L.set_names("d", a.organism) if s.startswith("d3")]
    lines = [dict(q) for q in L.load("d", a.organism, sets)[0] if q["qid"].startswith("d3|")]
    for q in lines:
        q.pop("set", None)
    if a.revision == "names_broad":
        from idtrack_reproducibility.families.d import revise_d1 as R1
        out = R1.names_broad(a.organism, lines, "d3_symbol", "genes_at_read", "read_at")
        R1.report(a.organism, out)
    else:
        out = [r for v in M.rev_par_shared_names(a.organism, lines, M.Ctx78(a.organism, H.ORGANISMS[a.organism][2])).values()
               for r in v]
    if a.dry and a.lines_to:
        import gzip
        import json
        with gzip.open(a.lines_to, "wt") as fo:
            for q in out:
                fo.write(json.dumps(q, sort_keys=True, default=str) + "\n")
    print(f"{a.organism}: {len(out)} lines revised ({a.revision})")
    if out and not a.dry:
        L.publish("d", a.organism, a.set, out, description=f"D3 revision {a.revision} (revise_d3.py docstring)",
                  generator=os.path.abspath(__file__), revises=sets, judge="78" if a.revision == "names_broad" else "rules")
        print("published", a.set)


if __name__ == "__main__":
    main()
