"""Reports of a scored run (F0): per stratum and view, with Wilson and cluster-bootstrap intervals and latency.

Usage: report.py <scored folder> [--family f --organism o]   (scorer.py calls `write` itself)

Writes into the scored folder (from `outcomes/*.parquet` and `run.json`):

  scores.tsv          the harness's table: kind, variant (`*` = the whole kind, else `<view>,<from|default>`), n, right,
                      wrong, undetermined, timeout, error, rate, ci_low, ci_high (Wilson 95 %). For the harness's
                      sets it is comparable with score.py's scores.tsv line by line.
  report.tsv          family, organism, stratum, view (`*` = both), n, right, expected_error, wrong, undetermined,
                      timeout, error, not_supported, missing, rate, wilson_low, wilson_high, boot_low, boot_high
                      (cluster bootstrap over `unit`: 2,000 resamples, seed 20261118), units, units_wrong,
                      units_error (distinct units behind the wrong answers and errors), recovered (rule `safe`: answers
                      that named a gene), s_p50, s_p95, s_p99, s_max (seconds per call)
  wrong.tsv.gz, undetermined.tsv.gz   qid, kind, stratum, variant, unit, query, expected, broad, got, reason, rule
                      (rows ordered by qid and variant, so that a table is the same under every pandas version)
  errors.tsv, timeouts.tsv            qid, kind, stratum, variant, query, reason, seconds
                                      (every list in (qid, variant) order)

Streaming: the outcome partitions are read one at a time; memory holds the counts, one row per (stratum, view, unit)
and the seconds of the calls (family A's 121 million rows of human: a few GB).

Rate = (right + expected_error) / (right + expected_error + wrong + error).
"""
import glob
import gzip
import json
import os
import sys

from idtrack_reproducibility.families import common as C

OUTCOMES = ["right", "expected_error", "wrong", "undetermined", "timeout", "error", "not_supported"]


def outcome_parts(folder):
    """The outcome files of a scored folder: only the partitions of the scoring that run.json records, and the
    relation lines."""
    parts = sorted(glob.glob(os.path.join(folder, "outcomes", "*.parquet")))
    rj = os.path.join(folder, "run.json")
    if os.path.exists(rj) and json.load(open(rj)).get("parts"):
        n = json.load(open(rj))["parts"]
        parts = [p for p in parts if p.endswith(f"_of_{n}.parquet") or p.endswith("part_relations.parquet")]
    return parts


def read_outcomes(folder, columns=None, outcomes=None):
    """Every outcome row of a scored folder as one pandas frame (with `outcomes`: only rows of those outcomes, filtered
    while reading). For families of tens of millions of rows read partition by partition (`outcome_parts`)."""
    import pandas as pd
    import pyarrow.parquet as pq
    parts = outcome_parts(folder)
    if not parts:
        return pd.DataFrame()
    flt = [("outcome", "in", list(outcomes))] if outcomes else None
    return pd.concat([pq.read_table(p, columns=columns, filters=flt).to_pandas() for p in parts], ignore_index=True)


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return max(0.0, c - h), min(1.0, c + h)


def cluster_bootstrap(r_u, d_u, key, n_boot=C.N_BOOT, seed=C.SEED):
    """95 % percentile interval of sum(r)/sum(d) over units resampled with replacement (one generator per group,
    seeded by the run's seed and the group's name, so a group's interval does not depend on the others)."""
    import numpy as np
    k = len(r_u)
    if k == 0 or d_u.sum() == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng([seed, C.crc(key)])
    chunk = max(1, min(n_boot, int(4e7 // k)))
    rates = []
    done = 0
    while done < n_boot:
        b = min(chunk, n_boot - done)
        idx = rng.integers(0, k, size=(b, k))
        rs, ds = r_u[idx].sum(1), d_u[idx].sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            rates.append(np.where(ds > 0, rs / np.maximum(ds, 1), np.nan))
        done += b
    x = np.concatenate(rates)
    x = x[~np.isnan(x)]
    return (float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))) if len(x) else (float("nan"), float("nan"))


def _tsv(path, header, rows):
    op = gzip.open if path.endswith(".gz") else open
    with op(path + ".tmp", "wt") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join("" if x is None else str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
    os.replace(path + ".tmp", path)


def counts(g):
    c = g["outcome"].value_counts()
    return {o: int(c.get(o, 0)) for o in OUTCOMES}


LIST_COLS = ["qid", "kind", "stratum", "variant", "unit", "query", "expected", "broad", "got", "reason", "rule"]
EVENT_COLS = ["qid", "kind", "stratum", "variant", "query", "reason", "seconds"]
LISTS = (("wrong", "wrong.tsv.gz", LIST_COLS), ("undetermined", "undetermined.tsv.gz", LIST_COLS),
         ("error", "errors.tsv", EVENT_COLS), ("timeout", "timeouts.tsv", EVENT_COLS))


def _sorted_rows(chunks):
    """The rows of sorted chunk files (Parquet), merged in (qid, variant) order, as tuples (a null is None and is
    written as an empty field)."""
    import heapq
    import pyarrow.parquet as pq

    def rows(path):
        for b in pq.ParquetFile(path).iter_batches(batch_size=100_000):
            names = b.schema.names
            for r in b.to_pylist():  # a null is None under every pandas version (pandas 3 reads a null string as NaN)
                yield tuple(r[c] for c in names)
    return heapq.merge(*(rows(c) for c in chunks), key=lambda r: (r[0], r[3]))


def write(folder, family=None, organism=None):
    """scores.tsv, report.tsv and the lists of answers that are not right, streaming over the outcome partitions:
    memory holds the counts, one row per (stratum, view, unit) and the seconds of the calls, never the outcome rows
    of more than one partition (family A's 121 million rows of human need a few GB, not 350)."""
    import shutil
    from collections import defaultdict

    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    info = json.load(open(os.path.join(folder, "run.json"))) if os.path.exists(os.path.join(folder, "run.json")) else {}
    family, organism = family or info.get("family"), organism or info.get("organism")
    parts = outcome_parts(folder)
    kv_n, kv_c = defaultdict(int), defaultdict(lambda: dict.fromkeys(OUTCOMES, 0))  # (kind, variant): n, outcomes
    sv_n, sv_c, sv_rec = defaultdict(int), defaultdict(lambda: dict.fromkeys(OUTCOMES, 0)), defaultdict(int)
    secs = defaultdict(list)  # (stratum, view) and (stratum, "*"): the seconds of the calls
    units = None  # (stratum, view, unit) -> r, d, w, e (rows right, rows in the rate, rows wrong, rows error)
    tmp = os.path.join(folder, ".report_tmp")  # sorted chunks of the lists, removed at the end
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    chunks = defaultdict(list)
    total = 0
    cols = sorted(set(LIST_COLS) | set(EVENT_COLS) | {"view", "outcome", "seconds"})
    for i, p in enumerate(parts):
        df = pq.read_table(p, columns=cols).to_pandas()
        if df.empty:
            continue
        total += len(df)
        for (kind, var, o), n in df.groupby(["kind", "variant", "outcome"]).size().items():
            kv_n[(kind, var)] += int(n)
            if o in OUTCOMES:
                kv_c[(kind, var)][o] += int(n)
        for (st, v, o), n in df.groupby(["stratum", "view", "outcome"]).size().items():
            sv_n[(st, v)] += int(n)
            if o in OUTCOMES:
                sv_c[(st, v)][o] += int(n)
        for (st, v), n in df[df["reason"] == "recovered"].groupby(["stratum", "view"]).size().items():
            sv_rec[(st, v)] += int(n)
        df["_r"] = df["outcome"].isin(["right", "expected_error"]).astype(np.int64)
        df["_d"] = df["outcome"].isin(["right", "expected_error", "wrong", "error"]).astype(np.int64)
        df["_w"] = (df["outcome"] == "wrong").astype(np.int64)
        df["_e"] = (df["outcome"] == "error").astype(np.int64)
        agg = df.groupby(["stratum", "view", "unit"])[["_r", "_d", "_w", "_e"]].sum()
        units = agg if units is None else pd.concat([units, agg]).groupby(level=[0, 1, 2]).sum()
        # seconds of the calls: a (qid, seconds) once per group, as a single scoring's drop_duplicates counts them
        s1 = df[["stratum", "view", "qid", "seconds"]].drop_duplicates(["stratum", "view", "qid", "seconds"]).dropna(subset=["seconds"])
        for (st, v), s in s1.groupby(["stratum", "view"])["seconds"]:
            secs[(st, v)].append(s.to_numpy())
        s2 = df[["stratum", "qid", "seconds"]].drop_duplicates(["stratum", "qid", "seconds"]).dropna(subset=["seconds"])
        for st, s in s2.groupby("stratum")["seconds"]:
            secs[(st, "*")].append(s.to_numpy())
        for o, _, lc in LISTS:
            sub = df.loc[df["outcome"] == o, lc]
            if len(sub):
                c = os.path.join(tmp, f"{o}_{i}.parquet")
                pq.write_table(pa.Table.from_pandas(sub.sort_values(["qid", "variant"], kind="stable"), preserve_index=False), c)
                chunks[o].append(c)
        del df
    if not total:
        shutil.rmtree(tmp, ignore_errors=True)
        print("no outcomes")
        return
    # scores.tsv: the harness's table (kind, variant)
    kind_n, kind_c = defaultdict(int), defaultdict(lambda: dict.fromkeys(OUTCOMES, 0))
    for (kind, var), n in kv_n.items():
        kind_n[kind] += n
        for o in OUTCOMES:
            kind_c[kind][o] += kv_c[(kind, var)][o]
    rows = []
    for kind, var in sorted(list(kv_n) + [(k, "*") for k in kind_n], key=lambda kv: (kv[0], kv[1] != "*", kv[1])):
        n, c = (kind_n[kind], kind_c[kind]) if var == "*" else (kv_n[(kind, var)], kv_c[(kind, var)])
        right = c["right"] + c["expected_error"]
        den = right + c["wrong"] + c["error"]
        lo, hi = wilson(right, den)
        rate = right / den if den else float("nan")
        rows.append([kind, var, n, right, c["wrong"], c["undetermined"], c["timeout"], c["error"], f"{rate:.4f}",
                     f"{lo:.4f}", f"{hi:.4f}"])
    _tsv(os.path.join(folder, "scores.tsv"),
         ["kind", "variant", "n", "right", "wrong", "undetermined", "timeout", "error", "rate", "ci_low", "ci_high"], rows)
    # report.tsv: per stratum and view, with the clustered interval, units and latency
    missing = info.get("missing", {})
    st_n, st_c, st_rec = defaultdict(int), defaultdict(lambda: dict.fromkeys(OUTCOMES, 0)), defaultdict(int)
    for (st, v), n in sv_n.items():
        st_n[st] += n
        st_rec[st] += sv_rec.get((st, v), 0)
        for o in OUTCOMES:
            st_c[st][o] += sv_c[(st, v)][o]
    ustar = units.groupby(level=[0, 2]).sum()
    rep = []
    keys = list(sv_n) + [(st, "*") for st in st_n]
    for st, view in sorted(keys, key=lambda kv: (kv[0], kv[1] != "*", kv[1])):
        try:  # the units of the group (sorted by unit, as a single scoring's groupby gives them)
            u = ustar.loc[st] if view == "*" else units.loc[(st, view)]
        except KeyError:  # no row of the group has a unit
            u = pd.DataFrame(columns=["_r", "_d", "_w", "_e"], dtype=np.int64)
        if view == "*":
            n, c, rec = st_n[st], st_c[st], st_rec[st]
        else:
            n, c, rec = sv_n[(st, view)], sv_c[(st, view)], sv_rec.get((st, view), 0)
        right = c["right"] + c["expected_error"]
        den = right + c["wrong"] + c["error"]
        lo, hi = wilson(right, den)
        blo, bhi = cluster_bootstrap(u["_r"].to_numpy(), u["_d"].to_numpy(), f"{family}|{organism}|{st}|{view}")
        sec = np.concatenate(secs[(st, view)]) if secs.get((st, view)) else np.array([])
        q = np.quantile(sec, [0.5, 0.95, 0.99]) if len(sec) else [float("nan")] * 3
        rep.append([family, organism, st, view, n, c["right"], c["expected_error"], c["wrong"], c["undetermined"],
                    c["timeout"], c["error"], c["not_supported"], missing.get(st, 0) if view == "*" else "",
                    f"{right / den:.4f}" if den else "nan", f"{lo:.4f}", f"{hi:.4f}", f"{blo:.4f}", f"{bhi:.4f}",
                    len(u), int((u["_w"] > 0).sum()), int((u["_e"] > 0).sum()), rec,
                    f"{q[0]:.5f}", f"{q[1]:.5f}", f"{q[2]:.5f}", f"{sec.max():.4f}" if len(sec) else "nan"])
    _tsv(os.path.join(folder, "report.tsv"),
         ["family", "organism", "stratum", "view", "n", "right", "expected_error", "wrong", "undetermined", "timeout", "error",
          "not_supported", "missing", "rate", "wilson_low", "wilson_high", "boot_low", "boot_high", "units", "units_wrong",
          "units_error", "recovered", "s_p50", "s_p95", "s_p99", "s_max"], rep)
    # the answers that are not right, in (qid, variant) order: the sorted chunks of the partitions merged
    for o, fn, lc in LISTS:
        _tsv(os.path.join(folder, fn), lc, _sorted_rows(chunks[o]))
    shutil.rmtree(tmp, ignore_errors=True)
    print(open(os.path.join(folder, "report.tsv")).read() if len(rep) < 80 else f"{len(rep)} report rows")


if __name__ == "__main__":
    args = sys.argv[1:]
    fam = args[args.index("--family") + 1] if "--family" in args else None
    org = args[args.index("--organism") + 1] if "--organism" in args else None
    write(os.path.abspath(args[0]), fam, org)
