"""The record: one row per query line, view and subject (an idtrack tree or another tool), over every part (docs/records.md).

    from idtrack_reproducibility import records
    df = records.load(parts=["a"], organisms=["sus_scrofa"])            # a pandas frame, one long table
    df = records.load(folders=["<run>/sus_scrofa/e/scores_now"])         # any scored folder (a modest run's part)
    hunts = records.load_hunts("<hunts>/reproduce/d4c2dda")              # the hunts' reproducers, one row per call

The rows are the scored folders' `outcomes/*.parquet` (written by the common scorer, `families/scorer.py`) with the
folder's own facts added from its `run.json` and the answers' load records (the subject, its version, the graph's
digest, the Python, when it was scored) and, where the folder was classified, the cause of each wrong answer
(`causes.tsv.gz`) and its row of the known-limitations table. Record version 1 (`run.json` `record_version`) holds the
query, the truth and the answer on every row; a scoring of version 0 holds them on the rows that are not right only:
`load(..., complete=True)` fills `got`, `status`, `match` and `error_type` of its other rows from the run's answers
and `query`, `expected` and `broad` from the family's sets (the latest line of the qid in the sets the scoring
read), so that every row of every stored scoring is the same record.

Usage (a quick look): python -m idtrack_reproducibility.records [--parts a,e] [--organisms ...] [--complete] [--out x.parquet]
"""
import argparse
import csv
import glob
import gzip
import json
import os
import sys

from . import config as CFG

RECORD_VERSION = 1
# Columns of the record, in order (docs/records.md has each column's meaning, type and values).
COLUMNS = ["record_version", "subject_kind", "subject", "subject_version", "graph_digest", "asked_at", "endpoint",
           "python", "part", "scoring", "organism", "family", "set", "qid", "view", "stratum", "kind", "variant",
           "unit", "query", "judge", "expected", "broad", "outcome", "reason", "got", "n_got", "status", "match",
           "error_type", "seconds", "calls", "cause", "cause_class", "limitation", "scored_at", "folder"]
OUTCOMES = ("right", "wrong", "undetermined", "timeout", "error", "expected_error", "not_supported", "missing")


def registered_folders(parts=None, organisms=None, subjects=None):
    """[(part, subject, organism, scored folder)] of the runs `results.toml` names, with E's two scorings each."""
    from .check import ROOT
    from . import results as RES
    data = CFG.get().root
    out = []
    for e in RES.registry(ROOT).get("run", []):
        part = e["family"] if "/a0/" not in e["folder"] else "a0"
        if part == "a0" or (parts and part not in parts) or (subjects and e["tree"] not in subjects):
            continue
        for o in e["organisms"]:
            if organisms and o not in organisms:
                continue
            base = os.path.join(data, e["folder"].format(organism=o))
            subs = [os.path.dirname(os.path.join(base, s)) for s in e.get("scored", ["run.json"])]
            for f in subs:
                if glob.glob(os.path.join(f, "outcomes", "*.parquet")):
                    out.append((part, e["tree"], o, f))
    return out


def _run_info(folder):
    info = json.load(open(os.path.join(folder, "run.json")))
    loads = info.get("load") or []
    subject = info.get("subject") or {}  # another tool's runs (family G) name it; an idtrack tree's runs do not
    commit = info.get("commit") or ""
    return dict(record_version=int(info.get("record_version", 0)),
                subject_kind=subject.get("kind", "tree"), subject=subject.get("name", commit[:7]),
                subject_version=subject.get("version", commit), asked_at=subject.get("asked_at"),
                endpoint=subject.get("endpoint"),
                graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                python=",".join(sorted({str(x.get("python")) for x in loads if x.get("python")})),
                organism=info.get("organism"), family=info.get("family"), scored_at=info.get("scored_at"),
                run=info.get("run"), sets=list((info.get("sets") or {}).keys()))


def _causes(folder):
    """{(qid, view): (cause, class)} of a classified folder (wrong answers; `causes_und` for the undetermined)."""
    out = {}
    for name in ("causes.tsv.gz", "causes_und.tsv.gz"):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            with gzip.open(p, "rt", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                    out[(r["qid"], r["view"])] = (r["cause"], r["class"])
    return out


def _limitations():
    """{(family, organism, cause): limitation} of the known-limitations table."""
    from .check import ROOT
    p = os.path.join(ROOT, "results", "limitations", "limitations.tsv")
    out = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t", quoting=csv.QUOTE_NONE):
            out[(r["family"], r["organism"], r["cause"])] = r["limitation"]
    return out


def _complete(df, info):
    """Fill the rows a scoring of record version 0 left without answer or truth (right answers) from the run's answers
    and the family's sets."""
    import pandas as pd
    import pyarrow.parquet as pq
    from .families import lines as L
    need = df["got"].isna()
    if need.any() and info.get("run") and os.path.isdir(os.path.join(info["run"], "answers")):
        want = set(zip(df.loc[need, "qid"], df.loc[need, "view"]))
        got = {}
        for p in sorted(glob.glob(os.path.join(info["run"], "answers", "*.parquet"))):
            t = pq.read_table(p, columns=["qid", "view", "got", "status", "match", "error_type"]).to_pandas()
            for r in t.itertuples(index=False):
                if (r.qid, r.view) in want:
                    got[(r.qid, r.view)] = r
        for c in ("got", "status", "match", "error_type"):
            vals = []
            for q, v, cur in zip(df["qid"], df["view"], df[c] if c in df else [None] * len(df)):
                r = got.get((q, v))
                if r is None or not (cur is None or (isinstance(cur, float) and pd.isna(cur))):
                    vals.append(cur)
                else:
                    x = getattr(r, c)
                    vals.append(",".join(x) if c == "got" and x is not None else x)
            df[c] = vals
    need = df["expected"].isna()
    if need.any():
        fam = info.get("family")
        lines = L.load_qids(fam, info["organism"], set(df.loc[need, "qid"]), sets=info.get("sets") or None)
        exp, brd, qry = [], [], []
        for q, v, e, b, qq in zip(df["qid"], df["view"], df["expected"], df["broad"], df["query"]):
            ln = lines.get(q)
            if ln is None or not (e is None or (isinstance(e, float) and pd.isna(e))):
                exp.append(e), brd.append(b), qry.append(qq)
                continue
            exp.append(",".join(ln.get("expect") or []))
            brd.append(None if ln.get("broad") is None else ",".join(ln["broad"]))
            qd = {k: ln.get(k) for k in ("id", "from_release", "to_release", "final_database", "from_database",
                                         "expect_error") if ln.get(k) is not None or k in ("from_release", "final_database")}
            qd["strategy"] = v
            qry.append(json.dumps(qd, sort_keys=True, separators=(",", ":"), default=str))
        df["expected"], df["broad"], df["query"] = exp, brd, qry
    return df


def load_folder(folder, part=None, complete=False, columns=None, limitations=None):
    """The record rows of one scored folder (a pandas frame with `COLUMNS`)."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    info = _run_info(folder)
    tabs = [pq.read_table(p) for p in sorted(glob.glob(os.path.join(folder, "outcomes", "*.parquet")))]
    df = pa.concat_tables(tabs, promote_options="default").to_pandas() if tabs else pd.DataFrame()
    for c in ("status", "match", "error_type", "got", "query", "expected", "broad"):
        if c not in df:
            df[c] = None
    if complete and info["record_version"] < 1 and len(df):
        df = _complete(df, info)
    for k in ("record_version", "subject_kind", "subject", "subject_version", "graph_digest", "asked_at", "endpoint",
              "python", "organism", "scored_at"):
        df[k] = info[k]
    df["part"] = part or info["family"]
    df["scoring"] = os.path.basename(folder) if os.path.basename(folder).startswith("scores_") else ""
    df["folder"] = os.path.relpath(folder, CFG.get().root) if os.path.isabs(folder) else folder
    df["judge"] = df.pop("rule") if "rule" in df else None
    cz = _causes(folder)
    pairs = list(zip(df["qid"], df["view"]))
    df["cause"] = [cz.get(k, (None, None))[0] for k in pairs]
    df["cause_class"] = [cz.get(k, (None, None))[1] for k in pairs]
    lim = _limitations() if limitations is None else limitations
    fam = part or info["family"] or ""
    # the known-limitations table measures what users get: a `decided` scoring's causes are no limitation rows
    judged_now = df["scoring"].iloc[0:1].tolist() != ["scores_decided"]
    df["limitation"] = [lim.get((fam, o, c)) if c and judged_now else None for o, c in zip(df["organism"], df["cause"])]
    df = df.reindex(columns=COLUMNS)
    return df[columns] if columns else df


def load(parts=None, organisms=None, subjects=None, folders=None, complete=False, columns=None):
    """One long table over the given scored folders, or over the runs `results.toml` registers (filtered by part,
    organism and subject: a tree's name in the registry)."""
    import pandas as pd
    lim = _limitations()
    todo = [(None, None, None, f) for f in folders] if folders else registered_folders(parts, organisms, subjects)
    frames = [load_folder(f, part=p, complete=complete, columns=columns, limitations=lim) for p, _, _, f in todo]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns or COLUMNS)


def load_hunts(folder):
    """The hunts' reproducers of one run folder (`all_all_<organism>.json` or a modest run's `hunts.json`): one row per
    call, the finding as `stratum`, `outcome` `right` where the finding holds as recorded and `wrong` where it changed."""
    import pandas as pd
    rows = []
    paths = sorted(glob.glob(os.path.join(folder, "all_all_*.json"))) + sorted(glob.glob(os.path.join(folder, "hunts.json")))
    for p in paths:
        d = json.load(open(p))
        for f in d.get("results", []):
            for c in f.get("calls", []):
                call, obs = c.get("call") or {}, c.get("observed") or {}
                rows.append(dict(record_version=RECORD_VERSION, subject_kind="tree", subject=(d.get("commit") or "")[:7],
                                 subject_version=d.get("commit"), python=d.get("python"), part="hunts",
                                 organism=call.get("organism"), family="hunts", stratum=f.get("id"), kind=call.get("kind"),
                                 qid=f"{f.get('id')}|{call.get('call')}", view=call.get("view"),
                                 query=json.dumps(call.get("kwargs") or {}, sort_keys=True) + " " + str(call.get("id")),
                                 expected=json.dumps(call.get("recorded"), sort_keys=True), judge="hunt finding",
                                 outcome="right" if not c.get("differences") else "wrong",
                                 reason="; ".join(map(str, c.get("differences") or [])),
                                 got=",".join(obs.get("got") or []), n_got=len(obs.get("got") or []),
                                 status=obs.get("status"), match=obs.get("match"), error_type=obs.get("error_type"),
                                 seconds=obs.get("seconds"), scored_at=d.get("when"), folder=p))
    return pd.DataFrame(rows).reindex(columns=COLUMNS)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.records", description=__doc__.split("\n\n")[0])
    ap.add_argument("--parts", default=None)
    ap.add_argument("--organisms", default=None)
    ap.add_argument("--subjects", default=None)
    ap.add_argument("--folders", default=None, help="comma-separated scored folders (instead of the registry's runs)")
    ap.add_argument("--complete", action="store_true")
    ap.add_argument("--out", default=None, help="write the table as Parquet")
    a = ap.parse_args(argv)
    sp = lambda x: x.split(",") if x else None  # noqa: E731
    df = load(sp(a.parts), sp(a.organisms), sp(a.subjects), sp(a.folders), a.complete)
    print(f"{len(df):,} rows; by part and outcome:")
    if len(df):
        print(df.groupby(["part", "subject", "organism", "outcome"]).size().to_string())
    if a.out:
        df.to_parquet(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
