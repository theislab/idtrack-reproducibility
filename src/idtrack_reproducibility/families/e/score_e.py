"""Family E's scoring: F0's scorer under both judges, then one summary.

Usage: score_e.py --organism <o> --run <run folder>/<organism> [--sets e1,e2,...] [--judges now,decided]
                  [--scores-to <folder>]   (default: the run folder; the answers are only read)

- judge `now` (71690ab's docs): the base sets; judge `decided` (decided rules not built yet, docs/rules.md): the base sets with their
  revision sets `<set>_dec` (gen_e.py) and `<set>_dec<k>` (revise_e.py, later decided revisions). Each is the common scorer into `<scores-to>/scores_<judge>/`.
- family E's two rules are the common scorer's (`families.scorer`): `or_error` (an error whose classes include one of the line's
  `or_error`, `A|B`, is `expected_error`; an answer is judged by the line's rule) and `contains` (the published field
  `e_rule: contains`: an answer is right when every expected value is among its values: source inference, the band of
  releases holds the release). They were applied here after F0's scorer until 2026-11-19; the numbers are the same.
- `<scores-to>/e_summary.tsv`: per judge and stratum the outcomes, the rate, and `answered` (lines that named a gene or
  value; for `open` lines this is the measure: e.g. missing values that became genes).
"""
import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP

FAMILY = "e"


DEC = re.compile(r"^(?P<base>.+)_dec(?P<k>\d*)$")  # a revision set of the decided judge: <set>_dec, <set>_dec2, ...


def sets_for(organism, judge, only=None):
    """The sets of a judge in publication order: `now` the base sets, `decided` with their revision sets `<set>_dec<k>`."""
    names = L.set_names(FAMILY, organism)
    base = [s for s in names if not DEC.match(s)]
    if only:
        base = [s for s in base if s in only]
    if judge == "now":
        return base
    return [s for s in names if s in base or (DEC.match(s) and DEC.match(s).group("base") in base)]


def summary(run, organism, judges):
    """`run` is the folder holding scores_<judge>/ (the scores-to folder)."""
    rows = []
    for j in judges:
        folder = os.path.join(run, f"scores_{j}")
        df = RP.read_outcomes(folder, columns=["qid", "view", "stratum", "outcome", "n_got"])
        if df.empty:
            continue
        for st, g in df.groupby("stratum"):
            c = g["outcome"].value_counts().to_dict()
            r, ee, w, er = (c.get(k, 0) for k in ("right", "expected_error", "wrong", "error"))
            d = r + ee + w + er
            lo, hi = RP.wilson(r + ee, d) if d else (float("nan"), float("nan"))
            rows.append([j, organism, st, len(g), r, ee, w, c.get("undetermined", 0), c.get("timeout", 0), er,
                         f"{(r + ee) / d:.4f}" if d else "", f"{lo:.4f}" if d else "", f"{hi:.4f}" if d else "",
                         int((g["n_got"] > 0).sum()), int(((g["outcome"] == "undetermined") & (g["n_got"] > 0)).sum())])
    hdr = ["judge", "organism", "stratum", "n", "right", "expected_error", "wrong", "undetermined", "timeout", "error", "rate",
           "wilson_low", "wilson_high", "answered", "undetermined_answered"]
    with open(os.path.join(run, "e_summary.tsv"), "w") as fo:
        fo.write("\t".join(hdr) + "\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--sets", default=None)
    ap.add_argument("--judges", default="now,decided")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--scores-to", default=None, help="the folder of scores_<judge>/ and e_summary.tsv (default: --run)")
    a = ap.parse_args()
    run = os.path.abspath(a.run)
    to = os.path.abspath(a.scores_to or a.run)
    judges = a.judges.split(",")
    only = a.sets.split(",") if a.sets else None
    for j in judges:
        sets = sets_for(a.organism, j, only)
        folder = os.path.join(to, f"scores_{j}")
        cmd = [sys.executable, "-m", "idtrack_reproducibility.families.scorer", "--family", FAMILY, "--organism", a.organism,
               "--run", run, "--sets", ",".join(sets), "--scores-to", folder] + (["--limit", str(a.limit)] if a.limit else [])
        subprocess.run(cmd, check=True)
        sys.stderr.write(f"{j}: {len(sets)} sets scored\n")
    summary(to, a.organism, judges)


if __name__ == "__main__":
    main()
