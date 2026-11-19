"""Revision sets of D4 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d4.py --organism <o> --set <new set> [--revision blank_refuses_batch|decision52_open] [--dry]

- `blank_refuses_batch` (2026-11-19): `d4_v1` expected every item of a list its own truth and a missing value its
  error alone. Both trees document that a batch holding an empty or blank string is refused whole: d4c2dda's
  `convert_identifier_multiple` ("ValueError: If it holds empty strings (positions named)") and 71690ab's
  `convert_many` ("a batch that ... holds an empty string, raises with the positions of the offending items";
  `MissingValueError`). Every line of a list that holds "" or " " now expects that error (`expect_error`
  MissingValueError|d4c2dda:ValueError; the item's own truth kept as `expect_alone`); a list with only the other
  missing values (`nan`, `NA`, `None`, `null`) keeps rule `missing_values`'s item-level error for them and its genes for the rest.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L

ERR = "MissingValueError|d4c2dda:ValueError"


MISSING_STRINGS = {"nan", "na", "none", "null"}
OPEN_52 = ("rule `missing_values` (decided, not in 71690ab): a missing-value string (nan, NA, None, null) raises MissingValueError for its item; "
           "the trees before it read it as an identifier (family E's `now` judge leaves it open the same way)")


def revise_52(lines):
    """`decision52_open` (2026-11-19, `d4_r2`): `d4_v1` expected rule `missing_values`'s MissingValueError for the missing-value
    strings, a decided rule that d4c2dda and 71690ab do not have yet (family E scores its `e1_batch_d52` lines
    open under the `now` judge for this reason). The lines of such items get `open`, the expected error kept: an
    answer is undetermined on a tree without rule `missing_values`, the error is right on a tree with it. Found on 71690ab: the
    pig graph answers `nan` with SCN11A (whose synonym is NaN) and `NA` with XK."""
    out = []
    for q in lines:
        if str(q["id"]).strip().lower() in MISSING_STRINGS and q.get("expect_error") and not q.get("open"):
            new = {k: v for k, v in q.items() if k != "set"}
            new["open"] = OPEN_52
            new["revised_by"] = list(q.get("revised_by") or []) + ["decision52_open"]
            out.append(new)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--revision", default="blank_refuses_batch", choices=["blank_refuses_batch", "decision52_open"])
    a = ap.parse_args()
    sets = [s for s in L.set_names("d", a.organism) if s.startswith("d4")]
    lines = [q for q in L.load("d", a.organism, sets)[0] if q["qid"].startswith("d4|")]
    if a.revision == "decision52_open":
        out = revise_52(lines)
        print(f"{a.organism}: {len(out)} lines revised (decision52_open)")
        if out and not a.dry:
            L.publish("d", a.organism, a.set, out, description="D4 revision decision52_open (revise_d4.py revise_52)",
                      generator=os.path.abspath(__file__), revises=sets)
            print("published", a.set)
        return
    out = []
    for q in lines:
        if any(not str(x).strip() for x in q["batch"]) and q.get("expect_error") != ERR:
            new = {k: v for k, v in q.items() if k != "set"}
            new["expect_alone"] = new.get("expect") or []
            new["expect"] = []
            new["expect_error"] = ERR
            new["revised_by"] = list(q.get("revised_by") or []) + ["blank_refuses_batch"]
            out.append(new)
    print(f"{a.organism}: {len(out)} lines revised (lists holding a blank string)")
    if out and not a.dry:
        L.publish("d", a.organism, a.set, out, description="D4 revision blank_refuses_batch (revise_d4.py docstring)",
                  generator=os.path.abspath(__file__), revises=sets)
        print("published", a.set)


if __name__ == "__main__":
    main()
