"""F0's explain.py for family A's sets of tens of millions of lines: the same output, but only the lines of the
answers to explain are read (F0's explain.py loads the whole family into memory, which a 30-million-line set does not
fit in 40 GB). Everything else is F0's: the adapter, the record format, the resumable `.part` file.

Usage: explain_a.py --organism <o> --clone <clone> --graph <graph> --scored <scored folder> [--outcomes wrong]
"""
import argparse
import gzip
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report


def lines_of(org, sets, qids):
    """{qid: line} of the given qids from the family's sets (the latest set wins, as lines.load)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    d = L.query_dir("a", org)
    mans = sorted((json.load(open(os.path.join(d, f"{s}.manifest.json"))) for s in sets), key=lambda m: m["seq"])
    out, want = {}, pa.array(sorted(qids))
    for m in mans:
        p = L._set_file(d, m["name"])
        t = pq.read_table(p)
        t = t.filter(pc.is_in(t["qid"], value_set=want))
        for row in t.to_pylist():
            out[row["qid"]] = L.from_row(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--scored", required=True)
    ap.add_argument("--outcomes", default="wrong")
    ap.add_argument("--timeout", type=float, default=60.0)
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
    a.scored = os.path.abspath(a.scored)
    out = os.path.join(a.scored, f"explain_{a.outcomes.replace(',', '_')}.jsonl.gz")
    info = json.load(open(os.path.join(a.scored, "run.json")))
    df = report.read_outcomes(a.scored, ["qid", "view", "stratum", "unit", "outcome", "got"])
    df = df[df["outcome"].isin(a.outcomes.split(","))].sort_values(["stratum", "unit", "qid", "view"])
    want = list(zip(df["qid"], df["view"]))
    scored_got = {(q, v): g or "" for q, v, g in zip(df["qid"], df["view"], df["got"])}
    del df
    by = lines_of(a.organism, list(info["sets"]), {q for q, _ in want})
    part = out + ".part"
    done = set()
    if os.path.exists(part):
        for x in open(part):
            try:
                d = json.loads(x)
                done.add((d["qid"], d["view"]))
            except Exception:  # noqa: BLE001
                pass
    todo = [k for k in dict.fromkeys(want) if k not in done and k[0] in by]
    sys.stderr.write(f"{len(want)} wanted, {len(todo)} to ask\n")
    if todo:
        from idtrack_reproducibility.families import adapter
        ad = adapter.Adapter(a.clone, a.graph, a.organism, C.release_of_graph(a.graph), os.path.dirname(out),
                             timeout=a.timeout)
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
    with open(part) as fi, gzip.open(out + ".tmp", "wt") as fo:
        fo.write(fi.read())
    os.replace(out + ".tmp", out)
    print("written", out)


if __name__ == "__main__":
    main()
