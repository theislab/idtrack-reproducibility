"""The record: one row per query line, view and subject (an idtrack tree or another tool), over every part (docs/records.md).

    from idtrack_reproducibility import records
    df = records.load(parts=["a"], organisms=["sus_scrofa"])            # a pandas frame, one long table
    df = records.load(folders=["<run>/sus_scrofa/e/scores_now"])         # any scored folder (a modest run's part)
    hunts = records.load_hunts("<hunts>/reproduce/d4c2dda")              # the hunts' reproducers, one row per call
    corpus = records.load_corpus("<corpus>/runs/d4c2dda/sus_scrofa/run1.jsonl.gz")   # the regression corpus's calls
    h = records.load_h("<performance>/runs/<tree>/<organism>")           # family H's timed calls
    raw = records.load_raw("<families>/g/raw/<tool>/<version>/<organism>")  # another tool's raw responses (family G)

The rows are the scored folders' `outcomes/*.parquet` (written by the common scorer, `families/scorer.py`, and by the
harness's own scorer, `harness/score.py`) with the folder's own facts added from its `run.json` and the answers' load
records (the subject, its version, the graph's digest, the Python, when it was scored) and, where the folder was
classified, the cause of each wrong answer (`causes.tsv.gz`) and its row of the known-limitations table. Every row
holds the query, the truth, the answer as the subject returned it and the verdict (record version `RECORD_VERSION`;
a scoring of an earlier version is scored again from its stored answers: `idtrack-repro rescore`).

Usage (a quick look): python -m idtrack_reproducibility.records [--parts a,e] [--organisms ...] [--out x.parquet]
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
OUTCOMES = ("right", "wrong", "undetermined", "timeout", "error", "expected_error", "not_supported", "missing", "unjudged")
# The index of another tool's raw responses (family G): one row per request, beside the files that hold the responses
# as the tool returned them (docs/records.md, "Another tool").
RAW_COLUMNS = ["request_id", "tool", "tool_version", "organism", "endpoint", "method", "url", "params", "body", "qids",
               "asked_at", "http_status", "seconds", "bytes", "sha256", "file", "error"]


def registered_folders(parts=None, organisms=None, subjects=None, with_rows=True):
    """[(part, subject, organism, scored folder)] of the runs `results.toml` names, with E's two scorings each (with
    `with_rows`, only those whose rows of the record are in the data root)."""
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
                if not with_rows or glob.glob(os.path.join(f, "outcomes", "*.parquet")):
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


def load_folder(folder, part=None, columns=None, limitations=None):
    """The record rows of one scored folder (a pandas frame with `COLUMNS`)."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    info = _run_info(folder)
    if info["record_version"] < RECORD_VERSION:
        raise SystemExit(f"{folder}: a scoring of record version {info['record_version']} (before {RECORD_VERSION}); score "
                         f"it again from its stored answers: idtrack-repro rescore (or the scorer on its run)")
    tabs = [pq.read_table(p) for p in sorted(glob.glob(os.path.join(folder, "outcomes", "*.parquet")))]
    df = pa.concat_tables(tabs, promote_options="default").to_pandas() if tabs else pd.DataFrame(columns=COLUMNS)
    for k in ("record_version", "subject_kind", "subject", "subject_version", "graph_digest", "asked_at", "endpoint",
              "python", "organism", "scored_at"):
        df[k] = info[k]
    df["part"] = part or info["family"]
    df["scoring"] = os.path.basename(folder) if os.path.basename(folder).startswith("scores_") else ""
    df["folder"] = CFG.in_root(folder)
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


def load(parts=None, organisms=None, subjects=None, folders=None, columns=None):
    """One long table over the given scored folders, or over the runs `results.toml` registers (filtered by part,
    organism and subject: a tree's name in the registry)."""
    import pandas as pd
    lim = _limitations()
    todo = [(None, None, None, f) for f in folders] if folders else registered_folders(parts, organisms, subjects)
    frames = [load_folder(f, part=p, columns=columns, limitations=lim) for p, _, _, f in todo]
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


def _tree_of(path):
    """(name, commit) of the tree whose run a file of `<runs>/<tree>/<organism>/` is: the configured tree's commit."""
    name = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(path))))
    t = CFG.get().trees.get(name) or {}
    return name, t.get("commit", "")


def load_corpus(answers, corpus=None):
    """The regression corpus's answers of one run (`<corpus>/runs/<tree>/<organism>/run1.jsonl.gz` or `merged.jsonl.gz`)
    as rows of the record: one row per call, the call's arguments as `query`, the answer's target ids as `got` and its
    canonical form in `reason`, outcome `unjudged` (the corpus has no truth: docs/corpus.md); `corpus` is the corpus
    file (default `<corpus>/corpus/<organism>.jsonl`)."""
    import pandas as pd
    org = os.path.basename(os.path.dirname(os.path.abspath(answers)))
    corpus = corpus or os.path.join(CFG.get().corpus, "corpus", f"{org}.jsonl")
    calls = {}
    with open(corpus) as fh:
        for x in fh:
            c = json.loads(x)
            calls[c["qid"]] = c
    tree, commit = _tree_of(answers)
    rows, meta = [], {}
    opener = gzip.open if answers.endswith(".gz") else open
    with opener(answers, "rt") as fh:
        for x in fh:
            r = json.loads(x)
            if "meta" in r:
                meta = r["meta"]
                continue
            if "qid" not in r:
                continue
            c = calls.get(r["qid"], {})
            res = r.get("result") if isinstance(r.get("result"), dict) else {}
            got = [str(t) for t in (res.get("target_id") or []) if not isinstance(t, (dict, list))]
            err = r.get("error") or {}
            rows.append(dict(record_version=RECORD_VERSION, subject_kind="tree", subject=tree, subject_version=commit,
                             python=meta.get("python"), part="corpus", organism=org, family="corpus", qid=r["qid"],
                             stratum=c.get("stratum"), kind=r.get("entry"), variant=None, unit=None,
                             query=json.dumps(c.get("args"), sort_keys=True, default=str), judge="none (no truth)",
                             outcome="unjudged", reason=json.dumps(r.get("result"), sort_keys=True, default=str)[:2000]
                             if r.get("ok") else None, got=",".join(got), n_got=len(got),
                             status=res.get("match") and ("converted" if got else None), match=res.get("match"),
                             error_type=err.get("type") if isinstance(err, dict) else None, seconds=r.get("dt"),
                             calls=1, scored_at=meta.get("start"), folder=CFG.in_root(answers)))
    return pd.DataFrame(rows).reindex(columns=COLUMNS)


def load_h(folder):
    """Family H's timed calls of one run folder (`<performance>/runs/<tree>/<organism>/calls_rep<k>.calls.parquet`):
    one row per call and repetition, outcome `unjudged` (H measures time, not answers), `kind` the measure (`A`: one
    call per identifier with the defaults, caches cold; `F`: every corpus call with its own arguments)."""
    import pandas as pd
    import pyarrow.parquet as pq
    frames = []
    for p in sorted(glob.glob(os.path.join(folder, "*.calls.parquet"))):
        df = pq.read_table(p).to_pandas()
        df["folder"] = CFG.in_root(p)
        frames.append(df)
    return pd.concat(frames, ignore_index=True).reindex(columns=COLUMNS) if frames else pd.DataFrame(columns=COLUMNS)


def raw_schema():
    import pyarrow as pa
    types = {"http_status": pa.int32(), "seconds": pa.float64(), "bytes": pa.int64(), "qids": pa.list_(pa.string())}
    return pa.schema([(c, types.get(c, pa.string())) for c in RAW_COLUMNS])


def load_raw(folder):
    """Another tool's raw responses of one organism (`<families>/g/raw/<tool>/<version>/<organism>/index.parquet`): one
    row per request with the file that holds the response as the tool returned it (`RAW_COLUMNS`)."""
    import pyarrow.parquet as pq
    return pq.read_table(os.path.join(folder, "index.parquet"), schema=raw_schema()).to_pandas()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.records", description=__doc__.split("\n\n")[0])
    ap.add_argument("--parts", default=None)
    ap.add_argument("--organisms", default=None)
    ap.add_argument("--subjects", default=None)
    ap.add_argument("--folders", default=None, help="comma-separated scored folders (instead of the registry's runs)")
    ap.add_argument("--out", default=None, help="write the table as Parquet")
    a = ap.parse_args(argv)
    sp = lambda x: x.split(",") if x else None  # noqa: E731
    df = load(sp(a.parts), sp(a.organisms), sp(a.subjects), sp(a.folders))
    print(f"{len(df):,} rows; by part and outcome:")
    if len(df):
        print(df.groupby(["part", "subject", "organism", "outcome"]).size().to_string())
    if a.out:
        df.to_parquet(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
