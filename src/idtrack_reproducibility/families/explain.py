"""Explained paths of the answers that are not right (F0): the lines asked again with explain=True.

Usage:
  explain.py --family <f> --organism <o> --clone <clone> --graph <graph file> --scored <scored folder>
             [--outcomes wrong] [--per-unit 3] [--qids <file of qid<TAB>view lines>] [--out <file.jsonl.gz>]

Takes the (qid, view) of the scored folder whose outcome is in --outcomes (default `wrong`), at most --per-unit per
unit and stratum (default: all), asks each line once more through the adapter with explain=True on the tree and graph
given, and writes one JSON line per (qid, view) to --out (default `<scored>/explain_<outcomes>.jsonl.gz`): `qid`,
`view`, `rec` (the query line), `res` (`target_id`: the answer of this call, `graph_id`, `match`, `from_release`,
`alternative_locus_genes`, `paths`: per candidate `ensembl_gene`, `target_id`, `entry_gene`, `external_jump`,
`edges`, `edge_dbs`; the classifier's format for both api generations), `same_answer` (the explained call returned
the scored answer), `seconds`. Resumable (a `.part` file). Run it in the tree's Python (common.python_for) and, for
large sets, as a job.

Memory: only the rows to explain are read (each outcome partition filtered by outcome as it is read) and only their
query lines (`lines.load_qids`, streaming over the sets): a family of tens of millions of lines is explained in the
memory of its wrong answers. A scored folder with nothing to explain gets an empty output file (exit 0).
"""
import argparse
import gzip
import json
import os
import sys
import time

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report


def wanted(scored, outcomes, per_unit=0):
    """The (qid, view) to explain, in (stratum, unit, qid, view) order, and their scored answers."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    cols = ["qid", "view", "stratum", "unit", "outcome", "got"]
    tabs = []
    for p in report.outcome_parts(scored):
        t = pq.read_table(p, columns=cols, filters=[("outcome", "in", list(outcomes))])
        if t.num_rows:
            tabs.append(t)
    if not tabs:
        return [], {}
    df = pa.concat_tables(tabs).to_pandas()
    df = df.sort_values(["stratum", "unit", "qid", "view"], kind="stable")
    if per_unit:
        df = df.groupby(["stratum", "unit"], sort=False, dropna=False).head(per_unit)
    want = list(zip(df["qid"], df["view"]))
    got = {(q, v): g or "" for q, v, g in zip(df["qid"], df["view"], df["got"])}
    _ = pd
    return want, got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--scored", required=True)
    ap.add_argument("--outcomes", default="wrong")
    ap.add_argument("--per-unit", type=int, default=0)
    ap.add_argument("--qids", default=None)
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("--release", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
    a.scored = os.path.abspath(a.scored)
    out = a.out or os.path.join(a.scored, f"explain_{a.outcomes.replace(',', '_')}.jsonl.gz")
    info = json.load(open(os.path.join(a.scored, "run.json")))
    scored_got = {}
    if a.qids:
        want = [tuple(x.rstrip("\n").split("\t")[:2]) for x in open(a.qids) if x.strip()]
    else:
        want, scored_got = wanted(a.scored, a.outcomes.split(","), a.per_unit)
    by = L.load_qids(a.family, a.organism, {q for q, _ in want}, sets=list(info.get("sets") or []) or None,
                     judge=info.get("judge", "revised"))
    part = out + ".part"
    done = set()
    if os.path.exists(part):
        for x in open(part):
            try:
                d = json.loads(x)
                done.add((d["qid"], d["view"]))
            except Exception:  # noqa: BLE001  a line cut by an interruption is asked again
                pass
    todo = [k for k in dict.fromkeys(want) if k not in done and k[0] in by]
    sys.stderr.write(f"{len(want)} wanted, {len(todo)} to ask\n")
    if todo:
        from idtrack_reproducibility.families import adapter
        ad = adapter.Adapter(a.clone, a.graph, a.organism, a.release or C.release_of_graph(a.graph),
                             os.path.dirname(out), timeout=a.timeout)
        with open(part, "a") as fh:
            for i, (q, v) in enumerate(todo):
                rec = dict(by[q], strategy=v)
                t = time.time()
                (ans,) = ad.ask(rec, explain=True)
                res = None
                if ans["outcome"] == "answer":
                    res = dict(target_id=ans["got"], graph_id=ans.get("graph_id"), match=ans["match"],
                               from_release=ans["from_release"], alternative_locus_genes=ans["apart"] or [],
                               status=ans["status"], paths=json.loads(ans["path"] or "[]"))
                same = None if (q, v) not in scored_got else ",".join(ans["got"]) == scored_got[(q, v)]
                fh.write(json.dumps(dict(qid=q, view=v, rec=by[q], res=res, error=ans.get("error"), same_answer=same,
                                         seconds=round(time.time() - t, 3)), default=str) + "\n")
                if i % 200 == 0:
                    fh.flush()
                    sys.stderr.write(f"{i}/{len(todo)}\n")
    with gzip.open(out + ".tmp", "wt") as fo:
        if os.path.exists(part):  # nothing to explain: an empty file, so that the next step reads "no record"
            with open(part) as fi:
                for x in fi:
                    fo.write(x)
    os.replace(out + ".tmp", out)
    print("written", out, f"({len(want)} answers to explain)")


if __name__ == "__main__":
    main()
