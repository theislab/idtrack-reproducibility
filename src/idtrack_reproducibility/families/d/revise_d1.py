"""Revision sets of D1 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d1.py --organism <o> --revision <name> --set <new set name> [--dry]

Revisions:

- `name_several_set`: a symbol that several genes carry as their name in R's GTF (`truth_class` name_several) was
  published with `set_exact` (the harness's rule for a name of several genes: every view must return them all). The
  matrix's feature is one of those genes (a counts matrix holds one row per name; R's and Seurat's make.unique name
  the others `X.1`, which are the `made_unique` lines) and which one cannot be told from the name: `best` answering any
  non-empty part of the expected genes is right, `all` must still return them all: rule `set`. Found on d4c2dda's
  first run (MALAT1 at 77: the gene and a 95-base gene inside it both named MALAT1).
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L

HERE = os.path.dirname(os.path.abspath(__file__))


def name_several_set(q):
    if q.get("truth_class") == "name_several" and q.get("rule") == "set_exact":
        return dict(q, rule="set")
    return None


REVISIONS = {"name_several_set": name_several_set}
# Revisions that regenerate the lines with gen_d1.py's current truth and publish every line whose truth changed:
# - `mutual_first` (2026-11-19): the `dominant` class's `mutual` rule chose a small new gene lying inside a gene that
#   grew (CDC42 at 95: the model shares 3,598 bases with CDC42 at 115 and 2,914 with ENSG00000289694, which it holds
#   whole, while CDC42 holds under half of its own bases inside the model); the rule now also asks that the gene share
#   no fewer of the model's bases than any other gene. Lines whose truth changes are published again (dominant ->
#   open, mostly; a name of several genes whose genes' classes change). `name_several_set` is applied to them too.
REGENERATE = {"mutual_first"}
TRUTH_FIELDS = ("expect", "open", "open_within", "truth_class", "stratum", "rule")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--revision", required=True, choices=sorted(REVISIONS) + sorted(REGENERATE))
    ap.add_argument("--set", required=True)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    d1_sets = [s for s in L.set_names("d", a.organism) if s.startswith("d1")]
    lines, _, _ = L.load("d", a.organism, d1_sets)
    out = []
    if a.revision in REGENERATE:
        from idtrack_reproducibility.families.d import gen_d1
        latest = {q["qid"]: q for q in lines if q["qid"].startswith("d1|")}
        for q in gen_d1.build(a.organism, dry=True):
            q = name_several_set(q) or q
            p = latest.get(q["qid"])
            if p is None:
                raise SystemExit(f"{q['qid']}: not published; a revision only revises")
            if any(p.get(k) != q.get(k) for k in TRUTH_FIELDS):
                q["revised_by"] = list(p.get("revised_by") or []) + [a.revision]
                out.append(q)
        lines = list(latest.values())
    for q in (lines if a.revision not in REGENERATE else []):
        if not q["qid"].startswith("d1|"):
            continue
        new = REVISIONS[a.revision](q)
        if new is not None:
            new.pop("set", None)
            new["revised_by"] = list(q.get("revised_by") or []) + [a.revision]
            out.append(new)
    print(f"{a.organism}: {len(out)} lines revised by {a.revision}")
    if out and not a.dry:
        revises = sorted({q["set"] for q in lines if q["qid"].startswith("d1|")})
        L.publish("d", a.organism, a.set, out, description=f"D1 revision {a.revision} (revise_d1.py docstring)",
                  generator=os.path.abspath(__file__), revises=revises)
        print("published", a.set)


if __name__ == "__main__":
    main()
