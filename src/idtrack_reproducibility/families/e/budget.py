"""E6 budgets: the slowest queries of every family asked under budgets from one search state to the default.

Usage (a tree of the `convert` generation, 71690ab on: `budget_states` and `budget_seconds` are its arguments; d4c2dda
has neither):
  budget.py run --organism <o> --clone <clone> --graph <graph> --out <run>/<organism> --from <answers folder> [...]
  budget.py submit ... (one job)

The lines: the 200 slowest single calls of each run folder given with `--from` (F0's answers: the harness's sets
asked by F0's check run; family E's runs; A*'s), by their `seconds`, each with its own arguments. For each line:

- the unbounded answer U (`budget_states` 10^9, `budget_seconds` 10^6), as its canonical `to_dict()` (det.py);
- under `budget_states` b for b in 1, 2, 3, 5, 10, 20, 50, ..., 200,000 (the default) and `budget_seconds` 10^6:
  the outcome must be `SearchBudgetExceededError` or exactly U (never a partial answer), and once an answer is given
  at b every larger b gives U (a larger budget never changes a finished answer);
- under `budget_seconds` s for s in 10^-5 to 30 (states 10^9): the error or exactly U;
- in a batch (`convert_many`) of the line's identifier beside two others, under `budget_states` 1: the item's status is
  `error` with the budget error, and the other items are answered (rule `typed_api` g: a batch never stops).

Writes `<out>/e6/lines.jsonl.gz` (per line: the outcomes per budget) and `<out>/e6_report.tsv` (per family: lines,
violations by kind: partial answer, a finished answer changed by a larger budget, another error, a batch stopped).
"""
import argparse
import glob
import gzip
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

STATES = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10_000, 20_000, 50_000, 100_000, 200_000]
SECONDS = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 30.0]


def lines_of(family, organism, qids):
    """The latest line of each of the given qids (sets of tens of millions of lines read by qid, not whole)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    out = {}
    d = L.query_dir(family, organism)
    for s in L.set_names(family, organism):
        p = L._set_file(d, s)
        if p.endswith(".parquet"):
            t = pq.read_table(p, filters=[("qid", "in", sorted(qids))])
            for r in t.to_pylist():
                out[r["qid"]] = L.from_row(r)
        else:
            for q in L.iter_set(p):
                if q["qid"] in qids:
                    out[q["qid"]] = q
    return out


def slowest(run, family, organism, k=200):
    """The k slowest lines (by the seconds of their call) of a run folder, with the family of its lines."""
    import pyarrow.parquet as pq
    import pandas as pd
    parts = glob.glob(os.path.join(run, "answers", "*.parquet"))
    if not parts:
        return []
    df = pd.concat([pq.read_table(p, columns=["qid", "seconds", "outcome"]).to_pandas() for p in parts])
    df = df[df["outcome"] != "error"].groupby("qid")["seconds"].max().sort_values(ascending=False)
    want = list(df.index[: 5 * k])
    by = lines_of(family, organism, set(want))
    out = []
    for qid in want:
        q = by.get(qid)
        if q is None or (q.get("call") or "convert") != "convert" or q.get("kwargs") or (q.get("id_codec") or "str") != "str":
            continue
        out.append(dict(family=family, qid=qid, id=q["id"], seconds=float(df[qid]),
                        kw={x: q.get(x) for x in ("from_release", "to_release", "final_database", "from_database") if q.get(x) is not None}))
        if len(out) >= k:
            break
    return out


def run(a):
    sys.path.insert(0, a.clone)
    from idtrack_reproducibility.families import adapter as F0
    import idtrack
    from idtrack.errors import SearchBudgetExceededError as Budget
    from idtrack_reproducibility.families.e import det
    link = F0.link_graph(a.graph, a.out)
    api = idtrack.API.load(link)
    import gc
    gc.collect()
    gc.freeze()
    picked = []
    for spec in a.sources:
        fam, run_ = spec.split("=", 1)
        picked += slowest(run_, fam, a.organism)
    os.makedirs(os.path.join(a.out, "e6"), exist_ok=True)
    others = None
    res_rows = []
    t0 = time.time()
    for q in picked:
        def ask(**b):
            try:
                return ("answer", det.dig(api.convert(q["id"], **q["kw"], **b).to_dict() | {"path": None}))
            except Budget:
                return ("budget", None)
            except Exception as e:  # noqa: BLE001
                return ("error", f"{type(e).__name__}: {str(e)[:200]}")
        u = ask(budget_states=10 ** 9, budget_seconds=1e6)
        st = {b: ask(budget_states=b, budget_seconds=1e6) for b in STATES}
        se = {s: ask(budget_states=10 ** 9, budget_seconds=s) for s in SECONDS}
        v = []
        if u[0] != "answer":
            v.append(f"unbounded: {u[0]} {u[1]}")
        for b, r in list(st.items()) + [(f"s{s}", r) for s, r in se.items()]:
            if r[0] == "answer" and r[1] != u[1]:
                v.append(f"partial_or_other_answer@{b}")
            if r[0] == "error":
                v.append(f"other_error@{b}: {r[1]}")
        first = next((b for b in STATES if st[b][0] == "answer"), None)
        if first is not None and any(st[b][0] != "answer" for b in STATES if b > first):
            v.append(f"finished_then_stopped_after@{first}")
        # in a batch: the item stopped by the budget, the others answered
        if others is None:
            others = ["zz_nothing_zz", q["id"]]
        try:
            br = api.convert_many([q["id"], "zz_nothing_zz"], **q["kw"], budget_states=1, budget_seconds=1e6)
            if br[0].status == "error" and br[0].error.type != "SearchBudgetExceededError":
                v.append(f"batch_item_other_error: {br[0].error.type}")
            if br[1].status == "error":
                v.append("batch_other_item_error")
        except Exception as e:  # noqa: BLE001
            v.append(f"batch_raised: {type(e).__name__}: {str(e)[:200]}")
        res_rows.append(dict(q, unbounded=u, states={str(k): r[0] for k, r in st.items()}, first_finished_states=first,
                             seconds_outcomes={str(k): r[0] for k, r in se.items()}, violations=v))
    with gzip.open(os.path.join(a.out, "e6", "lines.jsonl.gz"), "wt") as fo:
        for r in res_rows:
            fo.write(json.dumps(r, default=str) + "\n")
    fams = sorted({r["family"] for r in res_rows})
    hdr = ["family", "organism", "lines", "lines_with_violation", "partial_answer", "changed_by_larger", "other_error",
           "batch", "unbounded_not_answered", "median_first_finished_states", "max_first_finished_states"]
    with open(os.path.join(a.out, "e6_report.tsv"), "w") as fo:
        fo.write("\t".join(hdr) + "\n")
        for f in fams:
            rr = [r for r in res_rows if r["family"] == f]
            ff = sorted(r["first_finished_states"] for r in rr if r["first_finished_states"] is not None)

            def cnt(p):
                return sum(any(x.startswith(p) for x in r["violations"]) for r in rr)
            fo.write("\t".join(str(x) for x in [f, a.organism, len(rr), sum(bool(r["violations"]) for r in rr),
                                                cnt("partial_or_other_answer"), cnt("finished_then_stopped"), cnt("other_error"),
                                                cnt("batch"), cnt("unbounded"), ff[len(ff) // 2] if ff else "", ff[-1] if ff else ""]) + "\n")
    json.dump(dict(organism=a.organism, lines=len(res_rows), seconds=round(time.time() - t0, 1), clone=a.clone,
                   commit=C.git_commit(a.clone)[0], sources=a.sources), open(os.path.join(a.out, "e6", "info.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["run", "submit"])
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--from", dest="sources", action="append", default=[], help="<family>=<run folder of that family>")
    ap.add_argument("--mem", type=int, default=16, help="GB of the submitted job (a mouse or human run on main took more than 16)")
    a = ap.parse_args()
    a.out, a.clone, a.graph = (os.path.abspath(x) for x in (a.out, a.clone, a.graph))
    if C.generation(a.clone) != "convert":
        raise SystemExit("E6 asks budget_states and budget_seconds: a tree of the convert generation (71690ab on)")
    if a.what == "run":
        return run(a)
    logs = os.path.join(a.out, "e6", "logs")
    os.makedirs(logs, exist_ok=True)
    src = " ".join(f"--from {s}" for s in a.sources)
    from idtrack_reproducibility import jobs
    print(jobs.submit(f"e6_{a.organism}", jobs.module(C.python_for(a.clone), "idtrack_reproducibility.families.e.budget",
                                                      f"run --organism {a.organism} --clone {a.clone} --graph {a.graph} "
                                                      f"--out {a.out} {src}"),
                      logs, cpus=2, mem_gb=a.mem, time="12:00:00"))


if __name__ == "__main__":
    main()
