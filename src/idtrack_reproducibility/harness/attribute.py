"""The attribution of a decision's effect: the harness's answers that moved between consecutive scorings, per step.

Usage: idtrack-repro attribute --run <name>=<folder> --run <name>=<folder> [--run ...] [--organisms <o,...>]
    --out <file.tsv>

Each `--run` names a scored harness run per organism: a folder as `harness-score --merge` writes it (`scores.tsv`,
`outcomes.tsv.gz`, `run.json`), with `{organism}` where the organism's name goes (e.g.
`--run r2=<harness>/runs/r2/{organism}/kinds`). The runs are steps in the order given: a tree, then the same tree with
one part of a rule added (a variant), then the next part, and so on (docs/harness.md, "Attribution"). For each pair of
consecutive runs and each organism, every query scored in both whose outcome moved (right, wrong, undetermined, timeout,
error) is counted by its move, exactly as `harness-compare` counts it (`--moved`); so each moved answer is attributed
to the step that moved it. The table has one row per step, organism and move: `step` (`<a>-><b>`), `organism`,
`move` (`<outcome in a> -> <outcome in b>`) and `answers`. An organism whose folder is missing in one of the two runs
is skipped with a note.
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.harness import compare_scores as CS


def moves(a: str, b: str) -> collections.Counter:
    """`{move: answers}` between two scored runs: the queries scored in both whose outcome differs."""
    oa, ob = CS.read_outcomes(a), CS.read_outcomes(b)
    return collections.Counter(f"{oa[q]['outcome']} -> {ob[q]['outcome']}" for q in set(oa) & set(ob)
                               if oa[q]["outcome"] != ob[q]["outcome"])


def attribute(runs: list[tuple[str, str]], organisms: list[str]) -> list[tuple[str, str, str, int]]:
    """Rows `(step, organism, move, answers)` for consecutive runs `(name, folder with {organism})`."""
    rows = []
    for (na, fa), (nb, fb) in zip(runs, runs[1:]):
        for org in organisms:
            a, b = fa.format(organism=org), fb.format(organism=org)
            if not (os.path.isdir(a) and os.path.isdir(b)):
                print(f"NOTE: {org} skipped in {na}->{nb}: {a if not os.path.isdir(a) else b} is missing", file=sys.stderr)
                continue
            for move, n in sorted(moves(a, b).items()):
                rows.append((f"{na}->{nb}", org, move, n))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro attribute", description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", action="append", required=True, help="<name>=<scored run folder with {organism}>")
    ap.add_argument("--organisms", default=",".join(CFG.ORGANISMS))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    runs = []
    for item in a.run:
        name, sep, folder = item.partition("=")
        if not sep or not name or "{organism}" not in folder:
            raise SystemExit(f"--run {item!r}: give <name>=<folder with {{organism}}>")
        runs.append((name, folder))
    if len(runs) < 2:
        raise SystemExit("give at least two runs")
    rows = attribute(runs, a.organisms.split(","))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("step\torganism\tmove\tanswers\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    for r in rows:
        print("\t".join(map(str, r)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
