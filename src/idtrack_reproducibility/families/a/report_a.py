"""Family A's report of a scored run: the rates along the travel matrix's axes and per pattern (docs/methods.md, A:
"metrics per stratum, organism, era of origin, direction, distance in releases and assembly crossing ... with
clustered intervals"; the hunts' findings as strata, principle 8).

Usage: report_a.py <scored run folder> --organism <o> [--out <folder>]

Reads `outcomes/*.parquet` of F0's scorer partition by partition (the sets hold tens of millions of lines) and the
pattern table `<families>/a/patterns/<organism>.tsv.gz` (patterns_a.py; from Ensembl's tables only). The qid gives the
release asked from and the target (`<kind>|<query>|<from or ->|<to>`); the era is the assembly of the release asked
from (for lines without `from_release`, of the release the truth reads the id at: the extra field `read_at`, else
the target). Writes into the folder:

- `axes.tsv`: group (`A1`, `A2`, `A3`, `A4`, `A5`, `A6`, `A7`: the stratum's prefix), axis (`era`, `direction`,
  `distance`, `crossing`, `origin`), value, view, n, right, wrong, undetermined, timeout, error, rate, Wilson 95 %,
  cluster bootstrap 95 % over `unit` (2,000 resamples, seed 20261118), units, units_wrong;
- `patterns.tsv`: group, pattern flag, view, the same columns;
- `wrong_units.tsv`: per group and unit, the wrong answers and their lines (the units behind the wrong answers).
"""
import argparse
import collections
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import report as R

BUCKETS = [(0, 0, "0"), (1, 1, "1"), (2, 5, "2-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 40, "21-40"), (41, 999, "41+")]


def bucket(d):
    a = abs(d)
    for lo, hi, name in BUCKETS:
        if lo <= a <= hi:
            return name
    return "?"


def main():
    import numpy as np
    import pandas as pd
    import pyarrow.parquet as pq
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or a.run
    org = a.organism
    asm = H.gtf_assemblies(org)
    pat = {}
    pp = os.path.join(C.BENCH, "a", "patterns", f"{org}.tsv.gz")
    if os.path.exists(pp):
        for r in H.read_tsv(pp):
            pat[r["id"]] = r["flags"].split(",")
    parts = sorted(glob.glob(os.path.join(a.run, "outcomes", "*.parquet")))
    agg_axes = collections.defaultdict(lambda: np.zeros(7, dtype=np.int64))  # (group, axis, value, view, unit)
    agg_pat = collections.defaultdict(lambda: np.zeros(7, dtype=np.int64))
    OUT = ["right", "wrong", "undetermined", "timeout", "error", "expected_error"]
    for p in parts:
        df = pq.read_table(p, columns=["qid", "view", "stratum", "unit", "outcome"]).to_pandas()
        if df.empty:
            continue
        sp = df["qid"].str.split("|", expand=True)
        df["q"], df["R"], df["t"] = sp[1], sp[2], sp[3].astype(int)
        df["group"] = df["stratum"].str.split("_").str[0]
        Rn = pd.to_numeric(df["R"], errors="coerce")
        df["R0"] = Rn.fillna(df["t"]).astype(int)
        df["era"] = df["R0"].map(asm).fillna("")
        df["direction"] = np.where(Rn.isna(), "default", np.where(df["t"] > df["R0"], "forward",
                                                                  np.where(df["t"] < df["R0"], "backward", "same")))
        df["distance"] = (df["t"] - df["R0"]).map(bucket)
        df["crossing"] = np.where(df["t"].map(asm).fillna("") != df["era"], "across", "within")
        df["origin"] = df["R0"].astype(str)
        df["base"] = df["q"].map(C.base_id)
        for o_i, o in enumerate(OUT):
            df[o] = (df["outcome"] == o).astype(np.int64)
        cols = OUT + ["one"]
        df["one"] = 1
        for axis in ("era", "direction", "distance", "crossing", "origin"):
            for view_all in (False, True):
                keys = ["group", axis, "unit"] + ([] if view_all else ["view"])
                g = df.groupby(keys)[cols].sum()
                for k, row in zip(g.index, g.to_numpy()):
                    if view_all:
                        grp, val, unit = k
                        view = "*"
                    else:
                        grp, val, unit, view = k
                    agg_axes[(grp, axis, str(val), view, unit)] += row
        if pat:
            df["flags"] = df["base"].map(lambda x: pat.get(x, ["not_in_patterns"]))
            ex = df[["group", "view", "unit", "flags"] + cols].explode("flags")
            for view_all in (False, True):
                keys = ["group", "flags", "unit"] + ([] if view_all else ["view"])
                g = ex.groupby(keys)[cols].sum()
                for k, row in zip(g.index, g.to_numpy()):
                    if view_all:
                        grp, val, unit = k
                        view = "*"
                    else:
                        grp, val, unit, view = k
                    agg_pat[(grp, "pattern", val, view, unit)] += row
        print("part", p, len(df), flush=True)

    def table(agg, path):
        by = collections.defaultdict(list)
        for (grp, axis, val, view, unit), v in agg.items():
            by[(grp, axis, val, view)].append((unit, v))
        rows = []
        for key in sorted(by):
            units = by[key]
            M = np.array([v for _, v in units])
            right, wrong, und, to, err, ee, n = M.sum(0)
            r = right + ee
            den = r + wrong + err
            lo, hi = R.wilson(r, den)
            ru, du = M[:, 0] + M[:, 5], M[:, 0] + M[:, 5] + M[:, 1] + M[:, 4]
            blo, bhi = R.cluster_bootstrap(ru, du, "|".join(map(str, (org,) + key)))
            uw = int((M[:, 1] > 0).sum())
            rows.append(list(key) + [n, right, wrong, und, to, err, f"{r / den:.5f}" if den else "nan", f"{lo:.5f}",
                                     f"{hi:.5f}", f"{blo:.5f}", f"{bhi:.5f}", len(units), uw])
        R._tsv(path, ["group", "axis", "value", "view", "n", "right", "wrong", "undetermined", "timeout", "error", "rate",
                      "wilson_low", "wilson_high", "boot_low", "boot_high", "units", "units_wrong"], rows)
        return rows

    table(agg_axes, os.path.join(out, "axes.tsv"))
    if pat:
        table(agg_pat, os.path.join(out, "patterns.tsv"))
    wu = collections.Counter()
    for (grp, axis, val, view, unit), v in agg_axes.items():
        if axis == "direction" and view == "*" and v[1]:
            wu[(grp, unit)] += int(v[1])
    R._tsv(os.path.join(out, "wrong_units.tsv"), ["group", "unit", "wrong"],
           [[g, u, n] for (g, u), n in sorted(wu.items(), key=lambda kv: -kv[1])])
    print("written", out)


if __name__ == "__main__":
    main()
