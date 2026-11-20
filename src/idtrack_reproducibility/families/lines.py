"""Query lines and query sets of the benchmark families (F0; docs/methods.md).

A query line is the accuracy harness's line (`qid`, `kind`, `id`, `from_release`, `to_release`, `final_database`,
`strategy`, `rule`, `cmp`, `expect`, `broad`, `locus`, `note`, and the corrected judge's `locus_mapped`, `support`,
`spellings`, `apart`, `revised_by`, `excluded`) plus:

- `family` (`astar`, `a`, ...; `harness` for the harness's own sets), `stratum` (default: `kind`), `unit` (the cluster
  of the bootstrap: gene id, lineage component, name, dataset; default: the id without its version);
- `strategy` may be `both`: the line is scored under `best` and under `all` (one call from 71690ab on, two calls on
  d4c2dda); `best` or `all` score one view, as the harness's lines do;
- `from_database` (the API argument), `expect_error` (families E, B6/B7: the error class the docs promise; the rule
  `expect_error`), `id_codec` (how `id` is turned into the argument: `str` (default), `json` (the JSON value: null,
  numbers, lists, dicts), `bytes_hex`, `float_nan`), and any other key a family needs (kept in `extra`).

A set is `<family dir>/queries/<organism>/<set>.jsonl.gz` (under a million lines) or `<set>.parquet` (a fixed
schema, `extra` as JSON), with `<set>.manifest.json`: counts per stratum, sha256 of the file, the order number `seq`,
time, generator, description, `revises` for a revision set. **A published set is never rewritten**: `publish`
refuses an existing name; a change of truth is a revision set holding only the changed lines under their published
qid, and the line of the latest set (highest `seq`) is the truth of a qid. The harness's sets (`family="harness"`)
are read where they are (`accuracy/queries/<organism>/`), in the harness's order (`harness_lib.set_names`), with its
exclusions (`harness_lib.EXCLUDED`, the `excluded` field) and its two judges (`published`: without the revision sets
of the corrected judge and of the judge following idtrack's decided rules, `harness_lib.REVISED_JUDGES`).
"""
import datetime
import gzip
import hashlib
import json
import os

from idtrack_reproducibility.families import common as C

# Typed columns of a Parquet set (and of the frame the scorer reads); every other key goes to `extra` (JSON).
TYPED = ["qid", "family", "stratum", "unit", "kind", "id", "id_codec", "from_release", "to_release", "final_database",
         "from_database", "strategy", "rule", "cmp", "expect", "broad", "apart", "revised_by", "expect_error", "note"]
LIST_COLS = ("expect", "broad", "apart", "revised_by")
INT_COLS = ("from_release", "to_release")


def schema():
    import pyarrow as pa
    f = []
    for k in TYPED:
        if k in LIST_COLS:
            f.append(pa.field(k, pa.list_(pa.string())))
        elif k in INT_COLS:
            f.append(pa.field(k, pa.int64()))
        else:
            f.append(pa.field(k, pa.string()))
    f += [pa.field("extra", pa.string()), pa.field("h", pa.uint32())]
    return pa.schema(f)


def views(line):
    s = line.get("strategy") or "best"
    return ["best", "all"] if s == "both" else [s]


def variant(line, view):
    return f"{view},{'from' if line.get('from_release') is not None else 'default'}"


def stratum_of(line):
    return line.get("stratum") or line["kind"]


def unit_of(line):
    u = line.get("unit")
    return str(u) if u not in (None, "") else C.base_id(line["id"])


def judged_by(line):
    """`published`, or `revised:<corrections of the corrected judge>` (as the harness's score.py names the judge)."""
    rb = line.get("revised_by")
    return "revised:" + "+".join(rb) if rb else "published"


def to_row(line):
    """A line as a row of the typed schema (`extra`: the other keys as JSON, or None)."""
    row = {k: line.get(k) for k in TYPED}
    for k in ("from_release", "to_release"):
        if row[k] is not None:
            row[k] = int(row[k])
    for k in LIST_COLS:
        if row[k] is not None:
            row[k] = [str(x) for x in row[k]]
    if row["id"] is not None and not isinstance(row["id"], str):
        raise ValueError(f"{line['qid']}: `id` must be a string (use `id_codec` for other values)")
    rest = {k: v for k, v in line.items() if k not in TYPED}
    row["extra"] = C.dumps(rest) if rest else None
    row["h"] = C.crc(line["qid"])
    return row


def from_row(row):
    """A typed row back to a line (None fields dropped except those the harness's rules read)."""
    line = {k: row.get(k) for k in TYPED if row.get(k) is not None or k in ("expect", "broad", "from_release",
                                                                              "final_database")}
    if line.get("expect") is None:
        line["expect"] = []
    if row.get("extra"):
        line.update(json.loads(row["extra"]))
    return line


# ---------------- sets on disk ----------------

def query_dir(family, organism):
    if family == "harness":
        return os.path.join(C.HARNESS_QUERIES, organism)
    return os.path.join(C.family_dir(family), "queries", organism)


def _set_file(d, name):
    for ext in (".parquet", ".jsonl.gz"):
        p = os.path.join(d, name + ext)
        if os.path.exists(p):
            return p
    return None


def set_names(family, organism):
    """The published sets in their order (a later set's line replaces an earlier one with the same qid)."""
    d = query_dir(family, organism)
    if not os.path.isdir(d):
        return []
    if family == "harness":
        from idtrack_reproducibility.harness import harness_lib as H
        return H.set_names(d)
    mans = []
    for f in os.listdir(d):
        if f.endswith(".manifest.json"):
            m = json.load(open(os.path.join(d, f)))
            mans.append((m["seq"], m["name"]))
    return [n for _, n in sorted(mans)]


def manifest(family, organism, name):
    return json.load(open(os.path.join(query_dir(family, organism), f"{name}.manifest.json")))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def iter_set(path, shard=None):
    """The lines of one set file (dicts); with shard=(i, n) only those with crc32(qid) % n == i."""
    if path.endswith(".parquet"):
        import pyarrow.parquet as pq
        import pyarrow as pa
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(batch_size=200_000):
            if shard is not None:
                h = batch.column("h").to_numpy()
                batch = batch.filter(pa.array(h % shard[1] == shard[0]))
            for r in batch.to_pylist():
                yield from_row(r)
        return
    with gzip.open(path, "rt") as fh:
        for raw in fh:
            rec = json.loads(raw)
            if shard is None or C.shard_of(rec["qid"], shard[1]) == shard[0]:
                yield rec


# ---------------- panels (panel.py): published selections of qids ----------------

def panel_dir(family, organism):
    return os.path.join(C.family_dir(family), "panels", organism)


def panel_names(family, organism):
    d = panel_dir(family, organism)
    if not os.path.isdir(d):
        return []
    names = [f[: -len(".manifest.json")] for f in os.listdir(d) if f.endswith(".manifest.json")]
    return sorted(names, key=lambda n: int(n.rsplit("_v", 1)[-1]) if n.rsplit("_v", 1)[-1].isdigit() else 0)


def is_panel(family, organism, name):
    return os.path.exists(os.path.join(panel_dir(family, organism), f"{name}.manifest.json"))


def panel_qids(family, organism, name):
    import pyarrow.parquet as pq
    return set(pq.read_table(os.path.join(panel_dir(family, organism), f"{name}.parquet"), columns=["qid"]).column("qid").to_pylist())


def split_panels(family, organism, names):
    """(panel names, set names) of a list given as `sets`."""
    names = list(names or [])
    pn = [n for n in names if is_panel(family, organism, n)]
    return pn, [n for n in names if n not in pn]


def resolve_sets(family, organism, sets=None, judge="revised"):
    names = list(sets) if sets else set_names(family, organism)
    if family == "harness" and not sets and judge == "published":
        from idtrack_reproducibility.harness import harness_lib as H
        names = [s for s in names if manifest(family, organism, s).get("judge") not in H.REVISED_JUDGES]
    return names


def load(family, organism, sets=None, strata=None, kinds=None, shard=None, judge="revised", limit=None, within=None):
    """The latest line of every qid over the sets (default: every published set, in order).

    Returns (lines, sha256 per set, excluded counts). Lines left out (`excluded` field; for the harness also
    `harness_lib.EXCLUDED`) are counted, never returned. `limit`: the first N lines of each stratum (smoke runs).
    `within` (set names): only the qids these sets publish (a run that asked some sets only, scored with the revision
    sets of those sets, which may also revise lines of sets the run did not ask).
    """
    pnames, plain = split_panels(family, organism, sets) if sets else ([], None)
    sel = set().union(*(panel_qids(family, organism, p) for p in pnames)) if pnames else None
    names = resolve_sets(family, organism, plain or None, judge)
    d = query_dir(family, organism)
    by_qid, shas = {}, {}
    bound = set() if within else None
    for p in pnames:  # a panel selects qids; their lines are the latest of the sets (all, or those named beside it)
        shas[p] = json.load(open(os.path.join(panel_dir(family, organism), f"{p}.manifest.json")))["sha256"]
    for s in names:
        p = _set_file(d, s)
        if p is None:
            raise SystemExit(f"{family}/{organism}: set {s} not found in {d}")
        shas[s] = manifest(family, organism, s)["sha256"]
        for rec in iter_set(p, shard):
            if sel is not None and rec["qid"] not in sel:
                continue
            if kinds and rec["kind"] not in kinds:
                continue
            if strata and stratum_of(rec) not in strata:
                continue
            rec["set"] = s
            rec.setdefault("family", family)
            by_qid[rec["qid"]] = rec
            if bound is not None and s in within:
                bound.add(rec["qid"])
    if bound is not None:
        by_qid = {q: r for q, r in by_qid.items() if q in bound}
    excluded = {}
    out = []
    ex_h = {}
    if family == "harness":
        from idtrack_reproducibility.harness import harness_lib as H
        ex_h = H.EXCLUDED
    for q in by_qid.values():
        why = None
        if (organism, q["kind"], q.get("final_database")) in ex_h:
            why = "harness_lib.EXCLUDED"
        elif q.get("excluded"):
            why = "excluded by " + judged_by(q)
        if why:
            k = f"{stratum_of(q)}: {why}"
            excluded[k] = excluded.get(k, 0) + 1
            continue
        out.append(q)
    if limit:
        per = {}
        out = [q for q in out if per.setdefault(stratum_of(q), []).append(1) or len(per[stratum_of(q)]) <= limit]
    return out, shas, excluded


def load_qids(family, organism, qids, sets=None, judge="revised"):
    """{qid: latest line} of the given qids only, streaming over the sets (a set of tens of millions of lines is read
    in batches and filtered; nothing else is held). The latest set (highest `seq`) wins, as in `load`; exclusions are
    not applied (the qids come from a scored run, which applied them)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    want = set(qids)
    out = {}
    if not want:
        return out
    value_set = pa.array(sorted(want), type=pa.string())
    d = query_dir(family, organism)
    plain = split_panels(family, organism, sets)[1] if sets else None
    for s in resolve_sets(family, organism, plain or None, judge):
        p = _set_file(d, s)
        if p is None:
            raise SystemExit(f"{family}/{organism}: set {s} not found in {d}")
        if p.endswith(".parquet"):
            for batch in pq.ParquetFile(p).iter_batches(batch_size=500_000):
                batch = batch.filter(pc.is_in(batch.column("qid"), value_set=value_set))
                for r in batch.to_pylist():
                    rec = from_row(r)
                    rec["set"] = s
                    rec.setdefault("family", family)
                    out[rec["qid"]] = rec
        else:
            with gzip.open(p, "rt") as fh:
                for raw in fh:
                    rec = json.loads(raw)
                    if rec["qid"] in want:
                        rec["set"] = s
                        rec.setdefault("family", family)
                        out[rec["qid"]] = rec
    return out


def frame(lines):
    """Lines as a pandas frame of the typed columns (lists as Python lists) plus `extra` and `h`."""
    import pandas as pd
    return pd.DataFrame([to_row(q) for q in lines], columns=[f.name for f in schema()])


# ---------------- publishing ----------------

def _qid_column(path, parquet=None):
    """The qids of a set file (a pyarrow array)."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    if parquet or (parquet is None and path.endswith(".parquet")):
        return pq.read_table(path, columns=["qid"]).column("qid").combine_chunks()
    with gzip.open(path, "rt") as fh:
        return pa.array([json.loads(x)["qid"] for x in fh], type=pa.string())


def publish(family, organism, name, lines, description, generator, revises=None, judge=None, sources=None,
            fmt=None, manifest_extra=None, dry=False, chunk=500_000, log=None):
    """Write a new set and its manifest, streaming; refuse an existing name or a qid twice in the set.

    `lines`: any iterable of lines (a generator for sets of tens of millions: nothing is held but one chunk).
    `generator`: the script (path) that made it, recorded with its sha256. `revises`: the sets whose lines this one
    revises (every qid must be published in one of them). `sources`: the truth files used (paths under the truth store,
    whose md5 is in SOURCES.tsv, or descriptions with a sha256). `fmt`: `parquet` or `jsonl.gz`; by default Parquet,
    except for a list of fewer than a million lines. `manifest_extra`: further manifest fields of the family. `dry`:
    count only (returns n and the counts per stratum; nothing written). The checks (qids distinct; a revision set's
    qids all published before) run on the written file before it is given its name; a failed check removes it.
    """
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    if family == "harness":
        raise SystemExit("the harness's sets are published by accuracy/make_queries.py")
    d = query_dir(family, organism)
    os.makedirs(d, exist_ok=True)
    if _set_file(d, name) or os.path.exists(os.path.join(d, f"{name}.manifest.json")):
        raise SystemExit(f"{family}/{organism}/{name} is published; a published set is never rewritten")
    if fmt is None:
        fmt = "jsonl.gz" if isinstance(lines, (list, tuple)) and len(lines) < 1_000_000 else "parquet"
    path = os.path.join(d, f"{name}.{fmt}")
    tmp = path + ".tmp"
    counts, n, buf = {}, 0, []
    writer = fo = None
    if not dry:
        if fmt == "parquet":
            writer = pq.ParquetWriter(tmp, schema(), compression="zstd")
        else:
            fo = gzip.open(tmp, "wt")
    try:
        for q in lines:
            q.setdefault("family", family)
            k = stratum_of(q)
            counts[k] = counts.get(k, 0) + 1
            n += 1
            if dry:
                continue
            if writer is not None:
                buf.append(to_row(q))
                if len(buf) >= chunk:
                    writer.write_table(pa.Table.from_pylist(buf, schema=schema()))
                    buf = []
            else:
                fo.write(C.dumps(q) + "\n")
            if log and n % 2_000_000 == 0:
                log("lines", n)
        if dry:
            return dict(n=n, counts=dict(sorted(counts.items())))
        if writer is not None:
            if buf:
                writer.write_table(pa.Table.from_pylist(buf, schema=schema()))
            writer.close()
        else:
            fo.close()
        qids = _qid_column(tmp, parquet=fmt == "parquet")
        if pc.count_distinct(qids).as_py() != len(qids):
            vc = pc.value_counts(qids)
            dup = vc.field("values").filter(pc.greater(vc.field("counts"), 1))
            raise SystemExit(f"qid twice in {name}: {len(dup)} qids, e.g. {dup[:3].to_pylist()}")
        if revises:
            earlier = pa.chunked_array([_qid_column(_set_file(d, s)) for s in revises]).combine_chunks()
            new = qids.filter(pc.invert(pc.is_in(qids, value_set=earlier)))
            if len(new):
                raise SystemExit(f"a revision set holds only published qids; {len(new)} are new, e.g. {new[:3].to_pylist()}")
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.replace(tmp, path)
    seqs = [json.load(open(os.path.join(d, f)))["seq"] for f in os.listdir(d) if f.endswith(".manifest.json")]
    gen = os.path.abspath(generator)
    man = dict(name=name, family=family, organism=organism, seq=max(seqs, default=0) + 1, format=fmt, n=n,
               counts=dict(sorted(counts.items())), sha256=sha256(path),
               published=datetime.datetime.now().isoformat(timespec="seconds"), description=description,
               generator=gen, generator_sha256=sha256(gen) if os.path.exists(gen) else None,
               revises=list(revises or []), judge=judge, sources=list(sources or []))
    man.update(manifest_extra or {})
    json.dump(man, open(os.path.join(d, f"{name}.manifest.json"), "w"), indent=1)
    os.chmod(path, 0o400)
    return man
