"""Revision sets of family E (a published set is never rewritten; a change of truth is a revision set).

Usage: revise_e.py --organism <o> --revision <name> [--dry]

- `e5_r1`: the line `e|e5_database_spelling|from_database|name_db` assumed that the name chosen for E5 is a symbol of
  the organism's symbol database (`gen_e.NAME_DB`); the truth says whether it is (BioMart 115: `hgnc`, `mgi`; VGNC's
  own sets through BioMart's VGNC ids). Where it is not (cynomolgus: USP24 is a UniProt gene name, no VGNC symbol), the
  relation has no premise: the line is `open` with the reason.
- `e5_r2`: the same line as a plain line (rule `set`, open): F0's scorer applies `open` to every rule but `relation`.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C  # noqa: F401
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.e import gen_e as G


REL_AS_SET = False


def e5_r1(org):
    qid = "e|e5_database_spelling|from_database|name_db"
    lines, _, _ = L.load("e", org, ["e5"])
    q = next((x for x in lines if x["qid"] == qid), None)
    if q is None:
        return []
    db = G.NAME_DB[org]
    if db == "VGNC Symbol":
        t = H.vgnc_symbols(org, 115)
        holds = bool(t and t[1].get(q["id"]))
    else:
        t = H.biomart_pairs(org, 115, "hgnc" if db == "HGNC Symbol" else "mgi", 1)
        holds = bool(t and t[1].get(q["id"]))
    if holds:
        return []
    r = dict(q, open=f"{q['id']} is not a {db} of this organism at 115 (the line's premise)", revised_by=["e5_r1"])
    if REL_AS_SET:
        # e5_r2: F0's scorer reads `open` before every rule but `relation`: the line becomes a plain line, open
        r.update(rule="set", expect=[], revised_by=["e5_r1", "e5_r2"])
        r.pop("relation", None)
    r.pop("set", None)
    return [r]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--revision", required=True, choices=["e5_r1", "e5_r2"])
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    global REL_AS_SET
    REL_AS_SET = a.revision == "e5_r2"
    rev = e5_r1(a.organism)
    if not rev:
        print(a.organism, a.revision, "no line changes")
        return
    m = L.publish("e", a.organism, a.revision, rev, "E5: a relation whose premise the truth does not hold, open", __file__,
                  revises=["e5"] + (["e5_r1"] if a.revision == "e5_r2" else []), judge="now", sources=[os.path.join(H.truth_dir(a.organism), "SOURCES.tsv")], dry=a.dry)
    print(a.organism, a.revision, m["n"])


if __name__ == "__main__":
    main()
