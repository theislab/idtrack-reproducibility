"""Family A's revision sets: when the shared truth changes (A*'s `truth.py`/`gen.py`, the harness's `make_queries.py`
or `assembly_map.py`), every published line whose truth changes gets a revision line under its published qid; a
published set is never rewritten (docs/families/f0.md, "Query sets").

Usage: revise_a.py <organism> --stratum A1..A8 --set <published set>[,<its revision sets in order>] --name <revision set>
       [--snapshot S] [--dry-run]

The generator of the stratum (gen_a.py, gen_a7.py) runs again with the code as it is now; each line is compared
with the published line of its qid (its Parquet row without the hash column); a line that differs goes into the
revision set (`revises=[<published set>]`). A qid the generator makes now that the set lacks, or one of the set that
it no longer makes, is counted and reported, never published (a change of the query itself is a new set).
"""
import argparse
import collections
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.a import gen_a as GA

L, C = GA.L, GA.C
ME = os.path.abspath(__file__)


def digest(row):
    r = {k: v for k, v in row.items() if k not in ("h",)}
    if r.get("extra"):
        r["extra"] = json.dumps(json.loads(r["extra"]), sort_keys=True)
    return hashlib.blake2b(json.dumps(r, sort_keys=True, default=str).encode(), digest_size=8).digest()


def published(org, name):
    import pyarrow.parquet as pq
    p = L._set_file(L.query_dir(GA.FAMILY, org), name)
    pf = pq.ParquetFile(p)
    out = {}
    for b in pf.iter_batches(batch_size=500_000):
        for row in b.to_pylist():
            out[row["qid"]] = digest(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--snapshot", type=int, default=None)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org = a.organism
    sets = a.set.split(",")  # the set and its earlier revisions, in order: the latest line of a qid counts
    GA.log("published", sets)
    old = {}
    for s_ in sets:
        old.update(published(org, s_))
    GA.log("lines", len(old))
    counts = collections.Counter()
    changed_by = collections.Counter()
    if a.stratum == "A7":
        from idtrack_reproducibility.families.a import gen_a7 as G7
        gen_lines = G7.lines(org, a.target)
        G = None
    else:
        G = GA.GenA(org, a.target)
        ids = sorted(G.tr.present)
        if a.stratum == "A6":
            ids, G.a6_class = G.a6_ids()
        if a.stratum == "A8":
            G.snapshot = a.snapshot
            ids = [g for g in ids if min(G.tr.present[g]) <= a.snapshot]
        gen_lines = getattr(G, a.stratum.lower())(ids)
    seen = set()

    def changed():
        for q in gen_lines:
            row = L.to_row(dict(q, family=GA.FAMILY))
            qid = row["qid"]
            seen.add(qid)
            d = old.get(qid)
            if d is None:
                counts["new qid (not published)"] += 1
                continue
            if d == digest(row):
                counts["unchanged"] += 1
                continue
            counts["changed"] += 1
            changed_by[q.get("stratum")] += 1
            yield q

    man = L.publish(GA.FAMILY, org, a.name, changed(), f"revision of {a.set}: lines whose truth changed with the shared "
                    f"truth's code (revise_a.py)", ME, revises=sets, judge="revised",
                    sources=[f"as {a.set}"], fmt="parquet",
                    manifest_extra=dict(stratum=a.stratum, revises_set=a.set,
                                        astar_gen_sha256=L.sha256(GA.AG.__file__), astar_truth_sha256=L.sha256(GA.TA.__file__),
                                        make_queries_sha256=L.sha256(GA.AG.M.__file__),
                                        assembly_map_sha256=L.sha256(os.path.join(C.ACCURACY, "assembly_map.py"))),
                    dry=a.dry_run, log=GA.log)
    counts["published qids no longer made"] = len(set(old) - seen)
    if man and not a.dry_run and man.get("n") == 0:  # nothing changed: the empty set just written is removed
        d = L.query_dir(GA.FAMILY, org)
        for f in (L._set_file(d, a.name), os.path.join(d, f"{a.name}.manifest.json")):
            if f and os.path.exists(f):
                os.remove(f)
        man = None
        GA.log("no line changed: no revision set")
    GA.log("revision", a.name, dict(counts), dict(changed_by))
    if man and not a.dry_run:
        mp = os.path.join(L.query_dir(GA.FAMILY, org), f"{a.name}.manifest.json")
        m = json.load(open(mp))
        m["revision_counts"] = dict(counts)
        json.dump(m, open(mp, "w"), indent=1)


if __name__ == "__main__":
    main()
