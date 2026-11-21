"""The panel of a family (the common machinery; docs/methods.md, "Measurement design"): the lines a stage asks during
development.

Usage: panel.py --family <f> --organism <o> [--scored <scored folder> ...] [--per-stratum 50000] [--dry]

A family runs exhaustively only on the reference tree (d4c2dda) and on the final tree; every stage in between asks
the family's panel. The panel of an organism is the union of

- **every line that was not right on any tree scored so far**: wrong, undetermined, error or timeout, in any view and
  under any judge (the outcome partitions of every scored folder of the family and organism: `--scored`, default
  every folder under `<families>/<family>/runs/` that has a `run.json` of that organism and was not a smoke run);
- **every hunt reproducer of the family**: the lines of its sets named `hunt*`, the lines whose stratum starts with
  `hunt`, and the lines of the units a family names as the hunts' (`HUNT_UNITS`: A* and A, the ids A*'s selection took
  from the weakness hunts);
- **a seeded stratified random sample**: per stratum (of the family's own sets, revision sets left out: they repeat
  qids) the `--per-stratum` lines (default 50,000) with the smallest hash of (seed, family, organism, qid) -- all of a
  stratum that is smaller. The sample does not depend on the order of reading and is the same in every run;
- **the lines the others are related to**: a line of rule `relation` is judged against the line it names (its
  `relation.with`), which joins the selection with it (reason `relation`), so that every relation chosen is scored.

A panel is a **selection of qids, not copies of lines**: `<families>/<family>/panels/<organism>/panel_v<k>.parquet` (`qid`,
`h`, `stratum`, `why`: `failure`, `hunt`, `sample`, joined by `+`) with `panel_v<k>.manifest.json` (counts per stratum
and reason, sha256, the scored folders with their commits and graph digests, seed, per-stratum size). It is published
like a set: never rewritten (read-only); a new version `panel_v<k+1>` (a superset of the last) when new trees add lines
that were not right, and nothing when nothing changed. The truth of a panel line is always the latest line of its qid
in the family's sets (a revision published after the panel applies): `lines.load(family, organism, sets=["panel_v1"])`
gives the latest lines of the panel's qids over every set of the family, `sets=["panel_v1", <sets>]` over those sets
only (e.g. a judge's); the runner, scorer and explain take `--sets panel_v1` as they take a set.
"""
import argparse
import collections
import datetime
import glob
import gzip
import hashlib
import json
import os
import sys

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP

NOT_RIGHT = ["wrong", "undetermined", "error", "timeout"]


def _astar_hunt_ids(org):
    """The ids A*'s selection took from the weakness hunts (pattern `hunt_*`), latest selection file."""
    d = os.path.join(C.BENCH, "astar", "selection")
    files = sorted(glob.glob(os.path.join(d, f"{org}.r*.tsv.gz")), key=lambda p: int(p.split(".r")[-1].split(".")[0]))
    files = files or [os.path.join(d, f"{org}.tsv.gz")]
    out = set()
    for p in files:
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rt") as fh:
            next(fh)
            for line in fh:
                c = line.split("\t")
                if "hunt" in c[1]:
                    out.add(c[0])
    return out


# family -> function(organism) -> units (the line field `unit`) that are the hunts' reproducers
HUNT_UNITS = {"astar": _astar_hunt_ids, "a": _astar_hunt_ids}


def hkey(seed):
    """A 16-character key for pandas' fixed-key hash, from the seed."""
    return hashlib.sha256(str(seed).encode()).hexdigest()[:16]


def scored_folders(family, organism):
    """Every scored folder of the family and organism (run.json of that organism, outcomes, not a smoke run)."""
    out = []
    for rj in sorted(glob.glob(os.path.join(C.family_dir(family), "runs", "**", "run.json"), recursive=True)):
        d = os.path.dirname(rj)
        try:
            info = json.load(open(rj))
        except Exception:  # noqa: BLE001
            continue
        if info.get("organism") != organism or info.get("limit") or not RP.outcome_parts(d):
            continue
        out.append(d)
    return out


def failures(folders):
    """{qid: stratum} of every row that was not right in any of the scored folders, and per folder its counts."""
    import pyarrow.parquet as pq
    out, per = {}, {}
    for d in folders:
        n = 0
        for p in RP.outcome_parts(d):
            t = pq.read_table(p, columns=["qid", "stratum"], filters=[("outcome", "in", NOT_RIGHT)])
            for q, s in zip(t.column("qid").to_pylist(), t.column("stratum").to_pylist()):
                out.setdefault(q, s)
            n += t.num_rows
        per[d] = n
    return out, per


def scan_sets(family, organism, per_stratum, seed, hunt_units):
    """(sample {qid: stratum}, hunt {qid: stratum}, lines per stratum) over the family's own sets (revisions left
    out). The sample: per stratum the `per_stratum` qids with the smallest hash (seed, family, organism, qid)."""
    import numpy as np
    import pandas as pd
    import pyarrow.parquet as pq
    key = hkey(f"{seed}|{family}|{organism}")
    best = {}  # stratum -> (hashes, qids) of the current smallest
    hunt, sizes = {}, collections.Counter()
    d = L.query_dir(family, organism)

    def take(qids, strata, units, set_name):
        qids = np.asarray(qids, dtype=object)
        strata = np.asarray(strata, dtype=object)
        hv = pd.util.hash_array(qids, hash_key=key, categorize=False)
        for st in pd.unique(strata):
            m = strata == st
            sizes[st] += int(m.sum())
            h, q = hv[m], qids[m]
            if st in best:
                h, q = np.concatenate([best[st][0], h]), np.concatenate([best[st][1], q])
            if len(h) > per_stratum:
                i = np.argpartition(h, per_stratum - 1)[:per_stratum]
                h, q = h[i], q[i]
            best[st] = (h, q)
        is_hunt = set_name.startswith("hunt")
        for q, st, u in zip(qids, strata, units):
            if is_hunt or str(st).startswith("hunt") or (hunt_units and u in hunt_units):
                hunt[q] = st

    for s in L.set_names(family, organism):
        m = L.manifest(family, organism, s)
        if m.get("revises"):
            continue
        p = L._set_file(d, s)
        if p.endswith(".parquet"):
            for b in pq.ParquetFile(p).iter_batches(batch_size=1_000_000, columns=["qid", "stratum", "kind", "unit", "id"]):
                qids = b.column("qid").to_pylist()
                strata = [st or k for st, k in zip(b.column("stratum").to_pylist(), b.column("kind").to_pylist())]
                units = [u if u not in (None, "") else C.base_id(i) for u, i in zip(b.column("unit").to_pylist(), b.column("id").to_pylist())]
                take(qids, strata, units, s)
        else:
            qids, strata, units = [], [], []
            for ln in L.iter_set(p):
                qids.append(ln["qid"])
                strata.append(L.stratum_of(ln))
                units.append(L.unit_of(ln))
            if qids:
                take(qids, strata, units, s)
    sample = {}
    for st, (h, q) in best.items():
        for x in q:
            sample[x] = st
    return sample, hunt, sizes


def run_info(d):
    info = json.load(open(os.path.join(d, "run.json")))
    return dict(folder=d, commit=(info.get("commit") or "")[:12], graph_digests=info.get("graph_digests"),
                judge=info.get("judge"), scored_at=info.get("scored_at"), sets=list(info.get("sets") or []))


def relation_partners(family, organism, qids):
    """{qid: stratum} of the lines that the relation lines among `qids` are judged against and `qids` lacks."""
    lines = L.load_qids(family, organism, qids)
    want = {q["relation"]["with"] for q in lines.values() if q.get("rule") == "relation" and q.get("relation")}
    want -= set(qids)
    return {q: L.stratum_of(x) for q, x in L.load_qids(family, organism, want).items()}


def build(family, organism, folders=None, per_stratum=50_000, seed=C.SEED, dry=False, log=print, kind="panel"):
    """A new selection `<kind>_v<k>`: `panel` (failures, hunt reproducers, the sample) or `modest` (the modest set's
    part of the family: hunt reproducers and the sample only, docs/modest.md)."""
    if kind not in ("panel", "modest"):
        raise SystemExit(f"kind must be panel or modest, not {kind!r}")
    if kind == "modest":
        folders = []
    folders = scored_folders(family, organism) if folders is None else [os.path.abspath(x) for x in folders]
    fail, per_folder = failures(folders)
    hu = HUNT_UNITS.get(family, lambda o: set())(organism)
    sample, hunt, sizes = scan_sets(family, organism, per_stratum, seed, hu)
    why = collections.defaultdict(set)
    strat = {}
    for src, d in (("failure", fail), ("hunt", hunt), ("sample", sample)):
        for q, st in d.items():
            why[q].add(src)
            strat.setdefault(q, st)
    for q, st in relation_partners(family, organism, set(why)).items():
        why[q].add("relation")
        strat.setdefault(q, st)
    pdir = L.panel_dir(family, organism)
    prev = [n for n in L.panel_names(family, organism) if n.startswith(f"{kind}_v")]
    if prev:
        last = L.panel_qids(family, organism, prev[-1])
        for q in last:  # a new version holds the last one
            if q not in why:
                why[q].add("earlier_panel")
        old = set(last)
        if set(why) == old:
            log(f"{family}/{organism}: nothing new since {prev[-1]} ({len(old):,} qids)")
            return None
    name = f"{kind}_v{len(prev) + 1}"
    rows = [dict(qid=q, h=C.crc(q), stratum=strat.get(q) or "", why="+".join(sorted(w))) for q, w in why.items()]
    rows.sort(key=lambda r: r["qid"])
    counts = collections.Counter(r["stratum"] for r in rows)
    by_why = collections.Counter()
    for r in rows:
        for w in r["why"].split("+"):
            by_why[w] += 1
    man = dict(name=name, kind=kind, family=family, organism=organism, n=len(rows), counts=dict(sorted(counts.items())),
               by_reason=dict(sorted(by_why.items())), per_stratum=per_stratum, seed=seed,
               lines_per_stratum=dict(sorted(sizes.items())), failures_per_folder=per_folder,
               scored=[run_info(d) for d in folders], previous=prev[-1] if prev else None,
               published=datetime.datetime.now().isoformat(timespec="seconds"), generator=os.path.abspath(__file__),
               generator_sha256=L.sha256(os.path.abspath(__file__)))
    log(f"{family}/{organism} {name}: {len(rows):,} qids; {dict(by_why)}; {len(folders)} scored folders")
    if dry:
        return man
    import pyarrow as pa
    import pyarrow.parquet as pq
    os.makedirs(pdir, exist_ok=True)
    path = os.path.join(pdir, f"{name}.parquet")
    if os.path.exists(path):
        raise SystemExit(f"{path} exists; a panel is never rewritten")
    sch = pa.schema([("qid", pa.string()), ("h", pa.uint32()), ("stratum", pa.string()), ("why", pa.string())])
    pq.write_table(pa.Table.from_pylist(rows, schema=sch), path + ".tmp", compression="zstd")
    os.replace(path + ".tmp", path)
    man["sha256"] = L.sha256(path)
    json.dump(man, open(os.path.join(pdir, f"{name}.manifest.json"), "w"), indent=1)
    os.chmod(path, 0o400)
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--scored", nargs="*", default=None, help="scored folders (default: every one of the family)")
    ap.add_argument("--per-stratum", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    m = build(a.family, a.organism, a.scored, a.per_stratum, a.seed, a.dry)
    if m:
        print(json.dumps({k: m[k] for k in ("name", "n", "by_reason")}))


if __name__ == "__main__":
    main()
