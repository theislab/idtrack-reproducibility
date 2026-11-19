"""Revision sets of D2 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d2.py --organism <o> --set <new set name> [--dry]

Every line is generated again by gen_d2.py (its current truth: synthetic, matrices, more) and every line whose truth
(`expect`, `open`) differs from the latest published line of its qid is published in the revision set, with `revised_by`:

- `grch37_site` (2026-11-19): human releases 76 to 115 hold, beside the main GTF's genes, the genes of Ensembl's GRCh37
  site, which idtrack's graph keeps as assembly 37 and release mode merges; the first truth read only the main GTFs, so
  lists whose ids also live on GRCh37 tied fewer releases than idtrack (d4c2dda: 1,828 of 5,440 synthetic lines wrong,
  1,382 of them a superset of the expected releases).

- `grch37_site_items` (2026-11-19, `--truth-only`): `d2_r1` regenerated the lines whole, and for the lists with noise of
  another release the noise is drawn from the ids that release holds and the drawn one does not, which the GRCh37
  site changed: 1,289 of the 2,263 lines of `d2_r1` (strata `*_other5`, `*_other20`) carry other items than the lines that were
  asked. `d2_r2` gives every line of the family its first publication's items with the truth computed on them (a
  published question never changes); `d2_r1` stays published and is superseded qid by qid.
- `core_genes` (2026-11-19, `--truth-only`): human 48 to 51 hold genes of the core database the GTF lacks (unplaced
  scaffolds; A*'s `TruthA` and finding hh-5), which idtrack's graph holds: lists of 50 to 54 tied fewer releases than
  the graph tells (d4c2dda: 31 synthetic lines wrong after `d2_r2`, every one a superset with 48-51).
- `versionless_releases` (2026-11-19, `--truth-only`): lists of versioned ids (`d2_v2`) from releases near 77 were
  expected at the releases whose GTF writes the versions, but pig and human GTFs write no versions before 77, while
  idtrack's graph holds the versions there (from the core tables; finding hh-5): such lines are open with the reason.
- `symbols_fold` (2026-11-19): a mouse table of symbols written in capitals (nowotschin2019_gutendoderm: 0610007P14RIK)
  was judged by Ensembl's names exactly as written (`Erg28`-style mixed case), so almost no item was a name of any
  release (best share 0.023); rosetta reads such a table without letter case, and so does the truth now.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.d import gen_d2 as G
from idtrack_reproducibility.families import lines as L


def truth_kind(q):
    """The truth's kind of a published line: its `truth_kind`, else its kind's label (d2_ids, d2_symbols, d2_mixed)."""
    if q.get("truth_kind"):
        return q["truth_kind"]
    if q["qid"].startswith("d2v|"):
        return "ids_versioned"
    return q["kind"][len("d2_"):]


def retruth(q, P):
    """The published line `q` (its question unchanged: batch, kwargs) with the truth computed again on its own items."""
    keep = {k: v for k, v in q.items() if k not in ("expect", "open", "best_share", "drawn_share", "drawn_in_expected", "stated_share",
                                                    "stated_in_expected", "set", "revised_by")}
    E, shares = P.truth(q["batch"], truth_kind(q))
    keep["expect"] = [str(r) for r in E]
    keep["best_share"] = round(shares[E[0]], 6)
    for k in ("drawn", "stated"):
        if keep.get(k) is not None:
            keep[f"{k}_share"] = round(shares.get(keep[k], 0.0), 6)
            keep[f"{k}_in_expected"] = keep[k] in E
    if truth_kind(q) == "ids_versioned":
        why = G.versionless_open(P, q["batch"], E, shares)
        if why:
            keep["open"] = why
    if keep.get("stated") is not None and keep["stated"] not in E:
        keep["open"] = (f"the stated release {keep['stated']} fits the list worse ({keep['stated_share']:.4f}) than "
                        f"{E[0]}..{E[-1]} ({keep['best_share']:.4f}) by Ensembl's GTFs")
    return keep


def main_truth_only(a):
    """`--truth-only`: every line's truth computed again on the items of its first publication (the question never
    changes); lines whose latest published truth or items differ are published."""
    d2_sets = [s for s in L.set_names("d", a.organism) if s.startswith("d2")]
    first, latest = {}, {}
    for s in d2_sets:
        for q in L.load("d", a.organism, [s])[0]:
            if q["qid"].startswith("d2"):
                first.setdefault(q["qid"], q)
                latest[q["qid"]] = q
    P = G.Presence(a.organism)
    import itertools  # noqa: F401
    G.synthetic_more.__globals__  # the versioned ids' truth needs P.index["ids_versioned"]: built by synthetic_more
    for _ in G.synthetic_more(a.organism, P):
        break
    out = []
    for qid, q0 in first.items():
        new = retruth(q0, P)
        p = latest[qid]
        if new.get("expect") != p.get("expect") or new.get("open") != p.get("open") or new.get("batch") != p.get("batch"):
            new["revised_by"] = list(p.get("revised_by") or []) + [a.revision]
            out.append(new)
    print(f"{a.organism}: {len(out)} lines revised ({a.revision}, truth only)")
    if out and not a.dry:
        L.publish("d", a.organism, a.set, out, description=f"D2 revision {a.revision}, truth only (revise_d2.py docstring)",
                  generator=os.path.abspath(__file__), revises=d2_sets)
        print("published", a.set)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--revision", default="grch37_site")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--truth-only", action="store_true")
    a = ap.parse_args()
    if a.truth_only:
        return main_truth_only(a)
    d2_sets = [s for s in L.set_names("d", a.organism) if s.startswith("d2")]
    old = {q["qid"]: q for q in L.load("d", a.organism, d2_sets)[0] if q["qid"].startswith("d2")}
    P = G.Presence(a.organism)
    out = []
    import itertools
    regen = itertools.chain(G.synthetic(a.organism, P), G.matrices(a.organism, P), G.synthetic_more(a.organism, P))
    for q in regen:
        p = old.get(q["qid"])
        if p is None:
            continue  # a line of a set not published for this organism (d2_v2 before its publication)
        if p.get("expect") != q.get("expect") or p.get("open") != q.get("open"):
            q["revised_by"] = list(p.get("revised_by") or []) + [a.revision]
            out.append(q)
            if len(out) <= 20:
                print(q["qid"], p["expect"][:3], "->", q["expect"][:3], q.get("best_share"))
    print(f"{a.organism}: {len(out)} lines revised")
    if out and not a.dry:
        L.publish("d", a.organism, a.set, out, description=f"D2 revision {a.revision} (revise_d2.py docstring)",
                  generator=os.path.abspath(__file__), revises=d2_sets)
        print("published", a.set)


if __name__ == "__main__":
    main()
