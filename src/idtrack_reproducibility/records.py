"""The record: one row per query line, view and subject (an idtrack tree or another tool), over every part (docs/records.md).

    from idtrack_reproducibility import records
    df = records.load(parts=["a"], organisms=["sus_scrofa"])            # a pandas frame, one long table
    df = records.load(folders=["<run>/sus_scrofa/e/scores_now"])         # any scored folder (a modest run's part)
    df = records.load(export="<export>", parts=["a"])                    # an export (`idtrack-repro records --out`)
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

Usage: idtrack-repro records [--parts a,e|all] [--organisms ...] [--subjects ...] [--out <folder>]
    without --out, the rows counted by part, subject, organism and outcome; with it, the export: one Parquet file per
    part, organism and subject (`part=<p>/organism=<o>/subject=<s>/records.parquet`), `schema.json` (the record's
    version, every column with its type and meaning, from `SCHEMA`; the parts kept as summaries) and `index.tsv` (rows,
    bytes and sha256 per file); `records.load(export=<folder>)` reads it back.
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
# The record's columns, in order: (column, type, meaning, values). This table is the one source of the column list
# (`COLUMNS`), of the export's Parquet types (`schema`), of its `schema.json`, and of the table in docs/records.md
# ("Columns", which `schema_table()` prints and the tests compare). Types: `int32`, `float64`, `string`; `string (JSON)`
# is a string that holds a JSON object.
SCHEMA = [
    ("record_version", "int32", "the version of the record the scored folder was written as",
     "`1`: every row holds the query, the truth and the answer"),
    ("subject_kind", "string", "what answered", "`tree` (an idtrack checkout), `tool` (another tool: family G)"),
    ("subject", "string", "its name",
     "a tree: its commit's first seven characters; a tool: its name (`run.json` `subject.name`)"),
    ("subject_version", "string", "its version",
     "a tree: the full commit; a tool: its version or release as the tool reports it"),
    ("graph_digest", "string", "the content digest of the graph a tree was asked on (`,` between several)", ""),
    ("asked_at", "string",
     "for a tool: the date (ISO 8601) its service was asked; empty for a tree (the scoring's `scored_at` says when)", ""),
    ("endpoint", "string", "for a tool: the service's address or the package and its version; empty for a tree", ""),
    ("python", "string", "the Python the subject ran in", "e.g. `3.11.12`"),
    ("part", "string", "the benchmark",
     "`harness` (the harness's own scorer), `f0` (the harness's sets through the common machinery), `astar`, `a`, `b`, "
     "`c`, `d`, `e`, `hunts`, `corpus`, `h`, a later family's name"),
    ("scoring", "string", "the judge's scoring folder where a part has several",
     "`scores_now`, `scores_decided` (family E, B), or empty"),
    ("organism", "string", "the organism", "the five organisms' Ensembl names"),
    ("family", "string", "the family of the line's sets", "`harness` for the harness's own sets"),
    ("set", "string",
     "the published set whose line is the truth (the latest line of the qid: a revision set's line replaces its base's)", ""),
    ("qid", "string", "the line's stable id within the family and organism", ""),
    ("view", "string", "the view scored", "`best`, `all`"),
    ("stratum", "string", "the line's stratum (the family's design)", ""),
    ("kind", "string", "the line's kind (the family's design)", ""),
    ("variant", "string", "the view and whether a release was given", "`best,from`, `best,default`, `all,from`, `all,default`"),
    ("unit", "string", "the cluster of the bootstrap (gene id, lineage, name, dataset)", ""),
    ("query", "string (JSON)", "the call",
     "`id`, `from_release`, `to_release`, `final_database`, `from_database`, `strategy` (the view), `expect_error`"),
    ("judge", "string", "the judge of the line", "`published`, or `revised:<corrections>`"),
    ("expected", "string", "the truth: the expected ids (`,` between)", ""),
    ("broad", "string", "the wider set Ensembl's data also allow (`,` between); empty where none", ""),
    ("outcome", "string", "the verdict",
     "`right`, `wrong`, `undetermined` (the truth cannot decide), `timeout`, `error`, `expected_error` (the error the "
     "docs promise), `not_supported` (a tool that cannot ask such a query), `missing` (no answer stored), `unjudged` (a "
     "call without truth: the regression corpus, H's timed calls)"),
    ("reason", "string", "why, in words (the rule that decided it)", ""),
    ("got", "string", "the answer's ids as the subject returned them (`,` between)", ""),
    ("n_got", "int32", "how many", ""),
    ("status", "string", "the answer's own status",
     "`converted`, `tie`, `no_target`, `no_conversion`, `not_found`, `unknown_version`, as the subject reports it"),
    ("match", "string", "how the identifier was matched", "as the subject reports it"),
    ("error_type", "string", "the class of an error raised", ""),
    ("seconds", "float64", "the time of the call that gave this view", ""),
    ("calls", "int32", "the calls made for the line", ""),
    ("cause", "string", "for a classified wrong (or undetermined) answer: its cause ([F0's causes](families/f0.md#causes))", ""),
    ("cause_class", "string", "the class of that cause", "classes `a` to `e`"),
    ("limitation", "string",
     "its row of [the known-limitations table](../results/limitations/limitations.tsv) (`L-...`); empty for a decided "
     "judge's scoring", ""),
    ("scored_at", "string", "when the scoring was made", ""),
    ("folder", "string", "the scored folder in the data root", ""),
]
COLUMNS = [c for c, *_ in SCHEMA]
# The parts that keep only some of their lines as rows and count the others (docs/records.md, "What each part keeps"):
# stated in the export's schema.json as they are.
SUMMARIES = [
    ("a0", "family A's graph fidelity (the graph against Ensembl's tables)",
     "every mismatch that is not right (`<check>.tsv.gz`)", "the matches, counted per class (`summary.tsv`)"),
    ("e7", "family E's determinism and concurrency", "every difference (`e7_differences.jsonl.gz`)",
     "the agreeing calls, counted (`e7_report.tsv`)"),
    ("e9", "family E's web layer", "every difference (`differences.jsonl`)", "the agreeing requests, counted (`report.tsv`)"),
    ("f_checks", "family F's checks of the matrices, metadata, order and result cache",
     "every violation (`f_violations.jsonl.gz`)", "the items of each check, counted (`f_checks.tsv`)"),
]
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


def _frame(df, info, folder, part, causes, limitations):
    """Outcome rows of a scored folder as rows of the record (`COLUMNS`): the folder's facts, the causes and the
    limitations added."""
    for k in ("record_version", "subject_kind", "subject", "subject_version", "graph_digest", "asked_at", "endpoint",
              "python", "organism", "scored_at"):
        df[k] = info[k]
    df["part"] = part or info["family"]
    df["scoring"] = os.path.basename(folder) if os.path.basename(folder).startswith("scores_") else ""
    df["folder"] = CFG.in_root(folder)
    df["judge"] = df.pop("rule") if "rule" in df else None
    pairs = list(zip(df["qid"], df["view"]))
    df["cause"] = [causes.get(k, (None, None))[0] for k in pairs]
    df["cause_class"] = [causes.get(k, (None, None))[1] for k in pairs]
    fam = part or info["family"] or ""
    # the known-limitations table measures what users get: a `decided` scoring's causes are no limitation rows
    judged_now = os.path.basename(folder) != "scores_decided"
    df["limitation"] = [limitations.get((fam, o, c)) if c and judged_now else None
                        for o, c in zip(df["organism"], df["cause"])]
    return df.reindex(columns=COLUMNS)


def _checked_info(folder):
    info = _run_info(folder)
    if info["record_version"] < RECORD_VERSION:
        raise SystemExit(f"{folder}: a scoring of record version {info['record_version']} (before {RECORD_VERSION}); score "
                         f"it again from its stored answers: idtrack-repro rescore (or the scorer on its run)")
    return info


def load_folder(folder, part=None, columns=None, limitations=None):
    """The record rows of one scored folder (a pandas frame with `COLUMNS`)."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    info = _checked_info(folder)
    tabs = [pq.read_table(p) for p in sorted(glob.glob(os.path.join(folder, "outcomes", "*.parquet")))]
    df = pa.concat_tables(tabs, promote_options="default").to_pandas() if tabs else pd.DataFrame(columns=COLUMNS)
    df = _frame(df, info, folder, part, _causes(folder), _limitations() if limitations is None else limitations)
    return df[columns] if columns else df


def load(parts=None, organisms=None, subjects=None, folders=None, columns=None, export=None):
    """One long table over the given scored folders, over an export's files (`export`, filtered by part, organism and
    subject), or over the runs `results.toml` registers (filtered by part, organism and subject: a tree's name in the
    registry)."""
    import pandas as pd
    if export:
        return load_export(export, parts, organisms, subjects, columns)
    lim = _limitations()
    todo = [(None, None, None, f) for f in folders] if folders else registered_folders(parts, organisms, subjects)
    frames = [load_folder(f, part=p, columns=columns, limitations=lim) for p, _, _, f in todo]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns or COLUMNS)


# ---------------------------------------------------------------------------------------------------------------------
# The export: one Parquet file per part, organism and subject, with the schema and an index


ARROW = {"int32": "int32", "float64": "float64", "string": "string", "string (JSON)": "string"}
INDEX_COLS = ["file", "part", "organism", "subject", "rows", "bytes", "sha256", "folders"]


def schema():
    """The record's Parquet schema (types from `SCHEMA`)."""
    import pyarrow as pa
    return pa.schema([(c, getattr(pa, ARROW[t])()) for c, t, *_ in SCHEMA])


def schema_table():
    """The Markdown table of the columns that docs/records.md shows ("Columns"), from `SCHEMA`."""
    lines = ["| column | type | meaning | values |", "| --- | --- | --- | --- |"]
    lines += [f"| `{c}` | {t} | {m} | {v} |" for c, t, m, v in SCHEMA]
    return "\n".join(lines)


def to_table(df):
    """A frame of the record as an Arrow table of `schema()` (columns in order; missing values as nulls)."""
    import pandas as pd
    import pyarrow as pa
    cols = {}
    for c, t, *_ in SCHEMA:
        x = df[c] if c in df else pd.Series([None] * len(df), dtype=object)
        if ARROW[t] == "string":
            try:
                cols[c] = pa.array(x, type=pa.string(), from_pandas=True)
            except (pa.ArrowTypeError, pa.ArrowInvalid):  # values that are not strings: their text
                cols[c] = pa.array([None if pd.isna(v) else str(v) for v in x.tolist()], type=pa.string())
        else:
            cols[c] = pa.array(pd.to_numeric(x, errors="coerce").astype("float64"), from_pandas=True).cast(
                getattr(pa, ARROW[t])())
    return pa.table(cols, schema=schema())


def _safe(value):
    """A value as a folder name of the partitions (`[A-Za-z0-9_.-]`, others as `_`)."""
    return "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in str(value or "none"))


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def schema_document():
    """`schema.json` of an export: the record's version, every column with its type and meaning, the partitions, and
    the parts kept as summaries."""
    from .check import ROOT
    from . import results as RES
    kept = []
    try:
        reg = RES.registry(ROOT).get("run", [])
    except OSError:
        reg = []
    for part, what, rows, counted in SUMMARIES:
        folders = [e["folder"] for e in reg if e.get("rows") and (part == "a0" and "/a0/" in e["folder"]
                                                                  or part != "a0" and part in e["folder"])]
        kept.append(dict(part=part, what=what, rows=rows, kept_only_as_summary=counted, registered=folders))
    return dict(record_version=RECORD_VERSION,
                columns=[dict(name=c, type=t, arrow=ARROW[t], meaning=m, values=v) for c, t, m, v in SCHEMA],
                partitioning=["part", "organism", "subject"],
                files="part=<part>/organism=<organism>/subject=<subject>/records.parquet (folder names: the values "
                      "with characters other than letters, digits, `_`, `.`, `-` as `_`)",
                kept_as_summaries=kept,
                other_rows=dict(hunts="records.load_hunts", corpus="records.load_corpus", h="records.load_h",
                                tool_raw_responses="records.load_raw"))


def export(out, parts=None, organisms=None, subjects=None, folders=None):
    """Write the record of the registered runs (or of `folders`) to the folder `out`: one Parquet file per part,
    organism and subject (`part=<p>/organism=<o>/subject=<s>/records.parquet`, written outcome file by outcome file),
    `schema.json` and `index.tsv` (rows, bytes and sha256 per file). Returns the index rows."""
    import pyarrow.parquet as pq
    lim = _limitations()
    todo = [(None, None, None, f) for f in folders] if folders else registered_folders(parts, organisms)
    groups = {}
    for part, tree, _, folder in todo:
        info = _checked_info(folder)
        key = (part or info["family"], info["organism"], info["subject"])
        if subjects and info["subject"] not in subjects and tree not in subjects:
            continue
        groups.setdefault(key, []).append((folder, info))
    os.makedirs(out, exist_ok=True)
    index = []
    for (part, organism, subject), members in sorted(groups.items()):
        rel = os.path.join(f"part={_safe(part)}", f"organism={_safe(organism)}", f"subject={_safe(subject)}",
                           "records.parquet")
        path = os.path.join(out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        n = 0
        with pq.ParquetWriter(path + ".tmp", schema(), compression="zstd") as writer:
            for folder, info in members:
                causes = _causes(folder)
                for p in sorted(glob.glob(os.path.join(folder, "outcomes", "*.parquet"))):
                    t = to_table(_frame(pq.read_table(p).to_pandas(), info, folder, part, causes, lim))
                    writer.write_table(t)
                    n += t.num_rows
        os.replace(path + ".tmp", path)
        index.append(dict(file=rel, part=part, organism=organism, subject=subject, rows=n, bytes=os.path.getsize(path),
                          sha256=_sha256(path), folders=",".join(CFG.in_root(f) for f, _ in members)))
    with open(os.path.join(out, "schema.json"), "w") as fh:
        json.dump(schema_document(), fh, indent=1)
        fh.write("\n")
    with open(os.path.join(out, "index.tsv"), "w") as fh:
        fh.write("\t".join(INDEX_COLS) + "\n")
        for r in index:
            fh.write("\t".join(str(r[c]) for c in INDEX_COLS) + "\n")
    return index


def load_export(folder, parts=None, organisms=None, subjects=None, columns=None):
    """The rows of an export (`export`), filtered by part, organism and subject: one long table with `COLUMNS` and
    their types as the schema gives them (pandas' nullable integers for `int32`)."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    doc = json.load(open(os.path.join(folder, "schema.json")))
    if doc.get("record_version") != RECORD_VERSION:
        raise SystemExit(f"{folder}: an export of record version {doc.get('record_version')}, not {RECORD_VERSION}")
    rows = list(csv.DictReader(open(os.path.join(folder, "index.tsv")), delimiter="\t", quoting=csv.QUOTE_NONE))
    keep = [r for r in rows if (not parts or r["part"] in parts) and (not organisms or r["organism"] in organisms)
            and (not subjects or r["subject"] in subjects)]
    tabs = [pq.read_table(os.path.join(folder, r["file"]), schema=schema(), columns=columns) for r in keep]
    if not tabs:
        return pd.DataFrame(columns=columns or COLUMNS)
    t = pa.concat_tables(tabs)
    df = t.to_pandas(types_mapper={pa.int32(): pd.Int32Dtype()}.get)
    return df.reindex(columns=columns or COLUMNS)


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
    ap = argparse.ArgumentParser(prog="idtrack-repro records", description=__doc__.split("\n\n")[0])
    ap.add_argument("--parts", default=None, help="comma-separated, or `all` (default: all)")
    ap.add_argument("--organisms", default=None)
    ap.add_argument("--subjects", default=None, help="a subject's name (a tree's short commit, a tool) or a tree's name "
                    "in the registry")
    ap.add_argument("--folders", default=None, help="comma-separated scored folders (instead of the registry's runs)")
    ap.add_argument("--out", default=None, help="the export's folder (one Parquet per part, organism and subject, "
                    "schema.json, index.tsv); without it, the rows are counted only")
    ap.add_argument("--schema-table", action="store_true", help="print the columns' table of docs/records.md")
    a = ap.parse_args(argv)
    if a.schema_table:
        print(schema_table())
        return 0
    sp = lambda x: None if not x or x == "all" else x.split(",")  # noqa: E731
    if a.out:
        index = export(a.out, sp(a.parts), sp(a.organisms), sp(a.subjects), sp(a.folders))
        print(f"{len(index)} files, {sum(r['rows'] for r in index):,} rows: {os.path.join(a.out, 'index.tsv')}")
        return 0
    df = load(sp(a.parts), sp(a.organisms), sp(a.subjects), sp(a.folders))
    print(f"{len(df):,} rows; by part and outcome:")
    if len(df):
        print(df.groupby(["part", "subject", "organism", "outcome"]).size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
