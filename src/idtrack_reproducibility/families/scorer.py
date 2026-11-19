"""The scorer (F0): the harness's rules, columnar over Parquet answers, plus `expect_error`, `safe`, `relation`,
`open` (`open_within`), `or_error` and `contains`.

Usage:
  scorer.py --family <f> --organism <o> --run <run folder>/<organism> [--sets a,b] [--strata x] [--kinds k]
            [--judge revised|published] [--scores-to <folder>] [--parts P] [--limit N]
  scorer.py ... --parts P --part p     one partition (crc32(qid) % P == p), e.g. P jobs in parallel; then
  scorer.py ... --parts P --merge      the relation lines, run.json and the reports

Reads the family's lines (lines.load: the latest line of every qid) and the runner's answers
(`<run>/answers/*.parquet`, one row per line and view), and writes into `--scores-to` (default: the run folder):

  outcomes/part_<p>_of_<P>.parquet  one row per line and view: qid, view, family, stratum, kind, variant, unit, set,
                                    outcome, reason, n_got, got (non-right rows only), seconds, calls, rule (the
                                    judge: `published` or `revised:<corrections>`), and for non-right rows the query,
                                    expected and broad
  scores.tsv, report.tsv, wrong.tsv.gz, undetermined.tsv.gz, errors.tsv, timeouts.tsv, run.json   (report.py)

The rules (docs/harness.md "Scoring rules" and the corrected judge; docs/methods.md, F0):

- the harness's rules (`set`, `set_exact`, `broad`, `locus`, `locus_mapped`, `apart`, `support`, `spellings`, `cmp`)
  are the harness's own function `score.score_one`, called with the line's strategy set to the view scored. A
  columnar fast path decides only what that function surely calls right (answer equal to `expect` after the line's
  `cmp`, on lines without `apart` and `spellings`); every other row goes through the function itself. The check
  (README) is that the harness's sets scored here equal score.py's outcomes line by line.
- `expect_error` (a line field: the error class the docs promise; `A|B` for alternatives, e.g. the class of each api
  generation): that class (or a subclass: the answer's `error_mro`) raised -> `expected_error` (counted as right in the rate); another class, or an answer -> `wrong`.
- `safe` (rule): nothing, or a non-empty part of `expect`, is right (reason `nothing` / `recovered`); another gene
  is wrong. Coverage (the share recovered) is read from `n_got`.
- `or_error` (a line field, `A|B`): an error of one of these classes (or a subclass) is `expected_error`; an answer
  is judged by the line's rule (a decision that allows a refusal or an answer; family E). First of all rules.
- `open` (a line field, a string): the truth leaves the line open; an answer is `undetermined` with that reason,
  before every other rule but `or_error` (a timeout or an error stays what it is: no answer was given). With
  `open_within` (a list of ids) beside it, only an answer within those ids (or empty) is undetermined; an answer with
  an id outside them is `wrong` (reason `open: <reason>; outside the place: <ids>`; ids compared after `cmp`).
- `contains` (rule; or family E's published field `e_rule: contains`): an answer is right when it is not empty and
  every expected value is among its values (e.g. source inference: the band of releases holds the release).
- `relation` (rule; family E8): no truth; the line's `relation` = {"op": equal|subset|superset|disjoint,
  "with": <qid>, "with_view": <view>} compares its answer with another line's (scored after all parts); `or_error`
  and `open` (with `open_within`) apply to it first, as to every line.

Outcomes: right, wrong, undetermined, timeout, error, expected_error (and not_supported for other tools). Rate =
(right + expected_error) / (right + expected_error + wrong + error); undetermined and timeout apart.
"""
import argparse
import datetime
import glob
import json
import math
import os
import sys

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

SEP = "\x1f"
OUT_COLS = ["qid", "view", "family", "stratum", "kind", "variant", "unit", "set", "outcome", "reason", "n_got", "got",
            "seconds", "calls", "rule", "query", "expected", "broad"]


def out_schema():
    import pyarrow as pa
    t = {"n_got": pa.int32(), "seconds": pa.float64(), "calls": pa.int32()}
    return pa.schema([(c, t.get(c, pa.string())) for c in OUT_COLS])


def norm_key(ids, cmp):
    if ids is None:
        return ""
    if cmp == "base":
        s = {C.base_id(x) for x in ids}
    else:
        s = set(ids)
    return SEP.join(sorted(s))


def read_answers(run, part, nparts):
    """The answers of one partition (h % nparts == part) as a pandas frame; the last answer of a (qid, view) wins."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    tabs = []
    for p in sorted(glob.glob(os.path.join(run, "answers", "*.parquet"))):
        t = pq.read_table(p, columns=[c for c in pq.ParquetFile(p).schema_arrow.names if c not in ("path", "traceback")])
        if nparts > 1:
            t = t.filter(pa.array(t.column("h").to_numpy() % nparts == part))
        tabs.append(t)
    if not tabs:
        return pd.DataFrame(columns=[f.name for f in __import__("idtrack_reproducibility.families.runner", fromlist=["runner"]).answer_schema()])
    df = pa.concat_tables(tabs).to_pandas()
    return df.drop_duplicates(["qid", "view"], keep="last")


def _list(x):
    """A list column's value (None, a list, or the NumPy array pandas makes of it) as a list."""
    return [] if x is None else list(x)


def or_error_outcome(line, ans):
    """`or_error` (`A|B`): an error of one of these classes (or a subclass: the answer's `error_mro`) is
    `expected_error`; None when the line has no `or_error` or the answer is no such error (the line's rule judges)."""
    allowed = line.get("or_error")
    if not allowed or ans is None or ans["outcome"] != "error":
        return None
    if set(allowed.split("|")) & set(_list(ans.get("error_mro")) or [ans.get("error_type")]):
        return "expected_error", f"raised {ans.get('error_type')}, which the line allows"
    return None


def open_outcome(line, got):
    """`open` (with `open_within`, a list of ids): the outcome of an answer on an open line. Without `open_within` every
    answer is undetermined; with it, an answer within `open_within` (or empty) is undetermined and an answer with an
    id outside it is wrong (ids compared after the line's `cmp`)."""
    within = line.get("open_within")
    if within is None:
        return "undetermined", f"open: {line['open']}"
    g = set(norm_key(got, line.get("cmp")).split(SEP)) - {""}
    out = g - set(norm_key(list(within), line.get("cmp")).split(SEP))
    if not out:
        return "undetermined", f"open: {line['open']}"
    return "wrong", f"open: {line['open']}; outside the place: {','.join(sorted(out)[:5])}"


def is_contains(line):
    """The rule `contains` (or family E's published field `e_rule: contains`)."""
    return line.get("rule") == "contains" or line.get("e_rule") == "contains"


def score_line(line, view, ans, places):
    """(outcome, reason) of one line and view; ans is the answer row (dict) or None."""
    from idtrack_reproducibility.harness import score
    if ans is None:
        return "missing", "not run"
    o = ans["outcome"]
    oe = or_error_outcome(line, ans)
    if oe:
        return oe
    if line.get("open") and o == "answer":
        # the truth leaves the line open (e.g. A*: the history walk and REST differ; no successor and no place)
        return open_outcome(line, _list(ans["got"]))
    ee = line.get("expect_error")
    if ee:
        mro = _list(ans.get("error_mro")) or [ans.get("error_type")]
        if o == "error" and set(ee.split("|")) & set(mro):
            return "expected_error", f"raised {ans.get('error_type')}"
        if o == "timeout":
            return "timeout", ans.get("error") or ""
        if o == "error":
            return "wrong", f"raised {ans.get('error_type')}, the docs promise {ee}: {(ans.get('error') or '')[:200]}"
        return "wrong", f"answered ({','.join(_list(ans['got'])[:5]) or 'nothing'}) where {ee} was expected"
    if o == "timeout":
        return "timeout", ans.get("error") or "over the budget"
    if o == "error":
        return "error", ans.get("error") or ""
    got = _list(ans["got"])
    rule = line.get("rule")
    if is_contains(line):
        g = set(got)
        if g and set(line.get("expect") or []) <= g:
            return "right", ""
        return "wrong", f"the expected values are not among the answer's: {','.join(sorted(g)[:5])}"
    if rule == "safe":
        g = set(norm_key(got, line.get("cmp")).split(SEP)) - {""}
        E = set(line.get("expect") or [])
        if not g:
            return "right", "nothing"
        if g <= E:
            return "right", "recovered"
        return "wrong", f"another gene: {','.join(sorted(g - E)[:5])}"
    rec = dict(line, strategy=view)
    flags = ans.get("flags")
    res = {"got": got, "flags": json.loads(flags) if isinstance(flags, str) and flags else {}}
    if ans.get("apart") is not None:
        res["apart"] = list(ans["apart"])
    return score.score_one(rec, res, places)


def relation_outcome(line, view, ans, other):
    rel = line["relation"]
    if ans is None or other is None:
        return "missing", "not run" if ans is None else f"the related line {rel['with']} has no answer"
    oe = or_error_outcome(line, ans)
    if oe:
        return oe
    if line.get("open") and ans["outcome"] == "answer":
        return open_outcome(line, _list(ans["got"]))
    if ans["outcome"] != "answer" or other["outcome"] != "answer":
        o = ans["outcome"] if ans["outcome"] != "answer" else "undetermined"
        return o, f"no answer to compare: {ans['outcome']} / {other['outcome']}"
    a, b = set(_list(ans["got"])), set(_list(other["got"]))
    ok = {"equal": a == b, "subset": a <= b, "superset": a >= b, "disjoint": not (a & b)}[rel["op"]]
    return ("right", "") if ok else ("wrong", f"{rel['op']} fails against {rel['with']} ({rel.get('with_view', view)}): "
                                              f"{','.join(sorted(a)[:5])} vs {','.join(sorted(b)[:5])}")


def score_part(a, part, nparts, places):
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq

    lines, shas, excluded = L.load(a.family, a.organism, a.sets, a.strata, a.kinds, shard=(part, nparts) if nparts > 1 else None,
                                   judge=a.judge, limit=a.limit)
    ans = read_answers(a.run, part, nparts)
    ans_by = {(q, v): i for i, (q, v) in enumerate(zip(ans["qid"].tolist(), ans["view"].tolist()))}
    cols = {c: ans[c].tolist() for c in ans.columns} if len(ans) else {}

    def answer(qid, view):
        i = ans_by.get((qid, view))
        return None if i is None else {c: cols[c][i] for c in cols}

    rows, deferred, missing = [], [], {}
    for line in lines:
        st = L.stratum_of(line)
        unit = L.unit_of(line)
        rule_label = L.judged_by(line)
        for view in L.views(line):
            i = ans_by.get((line["qid"], view))
            base = dict(qid=line["qid"], view=view, family=line.get("family") or a.family, stratum=st, kind=line["kind"],
                        variant=L.variant(line, view), unit=unit, set=line.get("set"), rule=rule_label)
            if line.get("rule") == "relation":
                deferred.append((line, view, base))
                continue
            if i is None:
                missing[st] = missing.get(st, 0) + 1
                continue
            got = cols["got"][i]
            got = list(got) if got is not None else []
            sec, calls = cols["seconds"][i], cols["calls"][i]
            # the fast path: what the harness's function surely calls right
            fast = (cols["outcome"][i] == "answer" and not line.get("expect_error") and not line.get("open")
                    and line.get("rule") in ("set", "set_exact", "safe") and not is_contains(line)
                    and not line.get("apart") and not line.get("spellings")
                    and norm_key(got, line.get("cmp")) == SEP.join(sorted(set(line.get("expect") or []))))
            if fast:
                rows.append(dict(base, outcome="right", reason="" if line.get("rule") != "safe" else ("recovered" if got else "nothing"),
                                 n_got=len(got), got=None, seconds=sec, calls=calls, query=None, expected=None, broad=None))
                continue
            o, why = score_line(line, view, answer(line["qid"], view), places)
            rows.append(_row(base, line, o, why, got, sec, calls))
    out_dir = a.scores_to
    os.makedirs(os.path.join(out_dir, "outcomes"), exist_ok=True)
    path = os.path.join(out_dir, "outcomes", f"part_{part}_of_{nparts}.parquet")
    pq.write_table(pa.Table.from_pylist(rows, schema=out_schema()), path + ".tmp", compression="zstd")
    os.replace(path + ".tmp", path)
    return dict(n=len(rows), missing=missing, excluded=excluded, sets=shas, deferred=deferred)


def _row(base, line, o, why, got, sec, calls):
    right = o in ("right", "expected_error")
    q = {k: line.get(k) for k in ("id", "from_release", "to_release", "final_database", "strategy", "from_database",
                                  "expect_error") if line.get(k) is not None or k in ("from_release", "final_database")}
    q["strategy"] = base["view"]  # the view scored (a line of strategy `both` is scored twice)
    return dict(base, outcome=o, reason=why, n_got=len(got), got=None if right else ",".join(got), seconds=sec, calls=calls,
                query=None if right else C.dumps(q), expected=None if right else ",".join(line.get("expect") or []),
                broad=None if right or line.get("broad") is None else ",".join(line["broad"]))


def score_relations(a, deferred):
    """Relation lines: their own answers and those of the lines they name, read by qid from every answer file."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    path = os.path.join(a.scores_to, "outcomes", "part_relations.parquet")
    if not deferred:
        if os.path.exists(path):
            os.remove(path)
        return 0
    want = {d[0]["qid"] for d in deferred} | {d[0]["relation"]["with"] for d in deferred}
    tabs = []
    for p in sorted(glob.glob(os.path.join(a.run, "answers", "*.parquet"))):
        t = pq.read_table(p)
        tabs.append(t.filter(pc.is_in(t.column("qid"), value_set=pa.array(sorted(want)))))
    by = {}
    for r in pa.concat_tables(tabs).to_pylist():
        by[(r["qid"], r["view"])] = r
    rows = []
    for line, view, base in deferred:
        rel = line["relation"]
        own = by.get((line["qid"], view))
        other = by.get((rel["with"], rel.get("with_view", view)))
        o, why = relation_outcome(line, view, own, other)
        if o == "missing":
            continue
        got = list(own["got"] or []) if own else []
        rows.append(_row(base, line, o, why, got, own and own["seconds"], own and own["calls"]))
    pq.write_table(pa.Table.from_pylist(rows, schema=out_schema()), path + ".tmp", compression="zstd")
    os.replace(path + ".tmp", path)
    return len(rows)


def _part_stats_path(a, p, n):
    return os.path.join(a.scores_to, "outcomes", f"part_{p}_of_{n}.json")


def run_part(a, p, places):
    """Score one partition; its counts and its relation lines go beside its outcomes (for the merge)."""
    r = score_part(a, p, a.parts, places)
    deferred = [dict(line=line, view=view, base=base) for line, view, base in r.pop("deferred")]
    json.dump(dict(r, deferred=deferred), open(_part_stats_path(a, p, a.parts) + ".tmp", "w"), default=str)
    os.replace(_part_stats_path(a, p, a.parts) + ".tmp", _part_stats_path(a, p, a.parts))
    sys.stderr.write(f"part {p}/{a.parts}: {r['n']} rows\n")


def merge(a):
    """The relation lines, run.json and the reports, once every partition is scored."""
    stats = []
    for p in range(a.parts):
        sp = _part_stats_path(a, p, a.parts)
        if not os.path.exists(sp):
            raise SystemExit(f"partition {p} of {a.parts} is not scored ({sp})")
        stats.append(json.load(open(sp)))
    tot, missing, excluded, sets, deferred = 0, {}, {}, {}, []
    for r in stats:
        tot += r["n"]
        sets.update(r["sets"])
        deferred += [(d["line"], d["view"], d["base"]) for d in r["deferred"]]
        for d_, src in ((missing, r["missing"]), (excluded, r["excluded"])):
            for k, v in src.items():
                d_[k] = d_.get(k, 0) + v
    import re
    for old in glob.glob(os.path.join(a.scores_to, "outcomes", "part_*")):
        m = re.search(r"_of_(\d+)\.", os.path.basename(old))
        if m and int(m.group(1)) != a.parts:
            os.remove(old)  # a partition of an earlier scoring with another number of parts
    nrel = score_relations(a, deferred)
    loads = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(a.run, "answers", "load_*.json")))]
    clones = sorted({x.get("idtrack_file", "") for x in loads})
    digests = sorted({str(x.get("content_digest")) for x in loads})
    info = dict(family=a.family, organism=a.organism, run=a.run, judge=a.judge, sets=sets, strata=a.strata, kinds=a.kinds,
                limit=a.limit, n_rows=tot + nrel, n_relations=nrel, n_missing=sum(missing.values()), missing=missing,
                n_excluded=sum(excluded.values()), excluded=excluded, parts=a.parts,
                scored_at=datetime.datetime.now().isoformat(timespec="seconds"), idtrack_files=clones, graph_digests=digests,
                load=[{k: x.get(k) for k in ("generation", "load_seconds", "files_written_in_repo", "pythonhashseed", "graph_real",
                                             "python", "shard", "slurm_job")} for x in loads])
    if loads:
        clone = os.path.dirname(os.path.dirname(loads[0]["idtrack_file"]))
        info["clone"], (info["commit"], info["clone_modified"]) = clone, C.git_commit(clone)
    json.dump(info, open(os.path.join(a.scores_to, "run.json"), "w"), indent=1)
    from idtrack_reproducibility.families import report
    report.write(a.scores_to, a.family, a.organism)


def n_parts(run, per=3_000_000):
    import pyarrow.parquet as pq
    n = sum(pq.ParquetFile(p).metadata.num_rows for p in glob.glob(os.path.join(run, "answers", "*.parquet")))
    return max(1, math.ceil(n / per))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True, help="the run folder of this organism (holds answers/)")
    ap.add_argument("--sets", default=None)
    ap.add_argument("--strata", default=None)
    ap.add_argument("--kinds", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--judge", default="revised", choices=["revised", "published"])
    ap.add_argument("--scores-to", default=None)
    ap.add_argument("--parts", type=int, default=None, help="partitions (default: one per 3 million answers)")
    ap.add_argument("--part", type=int, default=None, help="score only this partition (jobs in parallel; then --merge)")
    ap.add_argument("--merge", action="store_true", help="only merge the scored partitions: relations, run.json, reports")
    a = ap.parse_args()
    a.run = os.path.abspath(a.run)
    a.scores_to = os.path.abspath(a.scores_to or a.run)
    for k in ("sets", "strata", "kinds"):
        setattr(a, k, getattr(a, k).split(",") if getattr(a, k) else None)
    if a.limit:
        a.parts = 1  # the first N lines of each stratum, as a runner of one shard took them
    if a.parts is None:
        if a.part is not None or a.merge:
            raise SystemExit("--part and --merge need --parts")
        a.parts = n_parts(a.run)
    if a.merge:
        return merge(a)
    from idtrack_reproducibility.harness import score
    places = score.Places(a.organism)
    if a.part is not None:
        return run_part(a, a.part, places)
    for old in glob.glob(os.path.join(a.scores_to, "outcomes", "part_*")):
        os.remove(old)  # a scoring is whole: the parts of an earlier scoring with another partition go
    for p in range(a.parts):
        run_part(a, p, places)
    merge(a)


if __name__ == "__main__":
    main()
