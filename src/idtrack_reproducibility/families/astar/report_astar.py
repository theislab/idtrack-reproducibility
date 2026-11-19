"""A*'s own report of one scored run (docs/methods.md, A* "Report"): per pattern, per id, per direction.

Usage: report_astar.py <scored folder> --organism <o>     (the folder F0's scorer wrote: outcomes/*.parquet)

Reads the outcomes (one row per line and view) and the selection table (the patterns of every id), and writes into
`<scored>/astar/`:
  per_pattern.tsv     pattern, view (`*` = both), n, right, wrong, undetermined, timeout, error, rate, Wilson 95 %,
                      cluster bootstrap 95 % over the lineage (`unit`; 2,000 resamples, seed 20261118), ids,
                      ids_wrong, lineages, lineages_wrong
  per_pattern_stratum.tsv   the same per pattern and stratum (direction, default, db; transcripts tx_/tl_)
  hard_ids.tsv.gz     one row per id: patterns, lineage, lineage size, rows, right, wrong, undetermined, timeout,
                      error, the wrong cells (`kind R->t view`, up to 40)
  open.tsv            the lines left open by the truth (`open: ...` reasons) per reason and stratum
An id counts under every pattern it has. Nothing of idtrack is read.
"""
import argparse
import collections
import glob
import os
import sys

import numpy as np
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H

OUT = ["right", "wrong", "undetermined", "timeout", "error"]


def wilson(k, n):
    return H.wilson(k, n) if n else (float("nan"), float("nan"))


def boot(units):
    """units: array (u, 3) of right, wrong+error, all -> 95 % interval of the rate by resampling units."""
    if len(units) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(C.SEED)
    r, bad = units[:, 0].astype(float), units[:, 1].astype(float)
    w = rng.multinomial(len(units), [1 / len(units)] * len(units), size=C.N_BOOT)
    num, den = w @ r, w @ (r + bad)
    rate = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return float(np.nanpercentile(rate, 2.5)), float(np.nanpercentile(rate, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--selection", default=None, help="selection file name (default: the organism; a later "
                    "selection holds every id of the earlier ones)")
    a = ap.parse_args()
    sel = {r["id"]: r for r in H.read_tsv(os.path.join(C.BENCH, "astar", "selection",
                                                       f"{a.selection or a.organism}.tsv.gz"))}
    per_id = collections.defaultdict(collections.Counter)  # (id, view, stratum) -> outcome counts
    wrong_cells = collections.defaultdict(list)
    opened = collections.Counter()
    for p in sorted(glob.glob(os.path.join(a.scored, "outcomes", "*.parquet"))):
        t = pq.read_table(p, columns=["qid", "view", "stratum", "outcome", "reason"]).to_pydict()
        for qid, view, st, o, why in zip(t["qid"], t["view"], t["stratum"], t["outcome"], t["reason"]):
            f = qid.split("|")
            gid = C.base_id(f[1])
            if f[0].startswith("astar_tx") or f[0].startswith("astar_tl"):
                gid = None  # transcript lines: their gene is in the line (`gene_at_origin`), read from the unit below
            per_id[(gid or f[1], view, st)][o] += 1
            if o == "wrong" and len(wrong_cells[gid or f[1]]) < 40:
                wrong_cells[gid or f[1]].append(f"{f[0]} {f[2]}->{f[3]} {view}")
            if o == "undetermined" and why.startswith("open: "):
                opened[(why, st)] += 1
    out = os.path.join(a.scored, "astar")
    os.makedirs(out, exist_ok=True)

    def table(key_of, path, keyname):
        agg = collections.defaultdict(collections.Counter)
        ids = collections.defaultdict(set)
        ids_wrong = collections.defaultdict(set)
        lin = collections.defaultdict(lambda: collections.defaultdict(lambda: np.zeros(3, dtype=np.int64)))
        for (gid, view, st), cnt in per_id.items():
            r = sel.get(gid)
            if r is None:
                continue
            for k in key_of([x for x in r["patterns"].split(",") if x], st):
                for v in (view, "*"):
                    agg[(k, v)].update(cnt)
                    ids[(k, v)].add(gid)
                    if cnt["wrong"]:
                        ids_wrong[(k, v)].add(gid)
                    u = lin[(k, v)][r["lineage"]]
                    u += (cnt["right"], cnt["wrong"] + cnt["error"], sum(cnt.values()))
        with open(path, "w") as fo:
            fo.write("\t".join([keyname, "view", "n"] + OUT + ["rate", "wilson_low", "wilson_high", "boot_low", "boot_high",
                                                             "ids", "ids_wrong", "lineages", "lineages_wrong"]) + "\n")
            for (k, v) in sorted(agg):
                c = agg[(k, v)]
                n = sum(c.values())
                den = c["right"] + c["wrong"] + c["error"]
                rate = c["right"] / den if den else float("nan")
                wl, wh = wilson(c["right"], den)
                units = np.array(list(lin[(k, v)].values()))
                bl, bh = boot(units)
                lw = int((units[:, 1] > 0).sum()) if len(units) else 0
                fo.write("\t".join(str(x) for x in [k, v, n] + [c[o] for o in OUT] +
                                   [f"{rate:.6f}", f"{wl:.6f}", f"{wh:.6f}", f"{bl:.6f}", f"{bh:.6f}",
                                    len(ids[(k, v)]), len(ids_wrong[(k, v)]), len(units), lw]) + "\n")

    table(lambda pp, st: set(pp) | {"ALL"}, os.path.join(out, "per_pattern.tsv"), "pattern")
    table(lambda pp, st: {f"{p}|{st}" for p in pp} | {f"ALL|{st}"}, os.path.join(out, "per_pattern_stratum.tsv"),
          "pattern|stratum")
    rows = []
    by_id = collections.defaultdict(collections.Counter)
    for (gid, view, st), cnt in per_id.items():
        by_id[gid].update(cnt)
    for gid in sorted(by_id):
        r = sel.get(gid, {})
        c = by_id[gid]
        rows.append([gid, r.get("patterns", ""), r.get("lineage", ""), r.get("lineage_size", ""), sum(c.values())] +
                    [c[o] for o in OUT] + ["; ".join(wrong_cells.get(gid, []))])
    H.write_tsv_gz(os.path.join(out, "hard_ids.tsv.gz"), ["id", "patterns", "lineage", "lineage_size", "rows"] + OUT +
                   ["wrong_cells"], rows)
    with open(os.path.join(out, "open.tsv"), "w") as fo:
        fo.write("reason\tstratum\trows\n")
        for (why, st), n in sorted(opened.items()):
            fo.write(f"{why}\t{st}\t{n}\n")
    print(open(os.path.join(out, "per_pattern.tsv")).read())


if __name__ == "__main__":
    main()
