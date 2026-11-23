"""Rows of the known-limitations table from a scored folder's outcomes, where no cause classification is needed:
the class of a wrong answer is read from the line's stratum and the answer itself (F0's record and the run's answers).

Usage: limits_outcomes.py --family <f> --map <map.tsv> --out <tsv> <scored folder> [<scored folder> ...]

The map (one per family, `families/<family>/limitations_map_outcomes.tsv`) has the columns `limitation`, `group`,
`title`, `what`, `user_action`, `fixed_by`, `rule`, `class` (the classes of `causes.py`: a, idtrack is wrong under its
own decisions; c, the graph lacks what the truth has; ...), `strata` (a regular expression on the stratum) and
`when`, one of

- `read_elsewhere`: the string was not read at the release asked (the answer says `from_release_fallback`, or the
  string is not found or not converted at all): the class of a name the graph does not hold at that release;
- `copy`: the answer's match is a copy of another item of the batch (`copy:<step>`);
- `reason:<regular expression>`: the judge's reason for the wrong answer (the outcome's `reason`: `not expected:
  XM_006530327`, `no answer`) matches it; `status:<regular expression>`: the answer's status (`no_target`) does: a
  class the answer's own fields name;
- `any`: every wrong answer of the strata the rows above left.

Every wrong outcome row (line and view) goes to the first map row whose `strata` and `when` hold; one table row per
(limitation, organism) with its measured size, as F0's builder writes them (`limitations.HEADER`). Nothing is typed by
hand.
"""
import argparse
import collections
import csv
import glob
import json
import os
import re
import sys

from idtrack_reproducibility.families import limitations as LIM

WHEN = ("read_elsewhere", "copy", "any")
FIELD_WHEN = ("reason", "status")  # `when` = <field>:<regular expression> on the wrong outcome's field


def load_map(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            kind, _, expr = r["when"].partition(":")
            if r["when"] not in WHEN and kind not in FIELD_WHEN:
                raise SystemExit(f"{path}: `when` must be one of {WHEN} or <reason|status>:<expression>, not {r['when']!r}")
            field = (kind, re.compile(expr)) if kind in FIELD_WHEN else None
            rows.append(dict(r, rx=re.compile(r["strata"]), field=field))
    return rows


def read_elsewhere(ans):
    """The string was not read at the release asked: a fallback, or nothing found or converted."""
    return bool(ans.get("fallback")) or ans.get("status") in ("not_found", "no_conversion")


def classify(m, stratum, ans):
    """The first map row (a dict) whose strata and condition hold for a wrong answer, or None."""
    for r in m:
        if not r["rx"].search(stratum):
            continue
        if r["when"] == "any" or (r["when"] == "copy" and (ans.get("match") or "").startswith("copy:")) or \
                (r["when"] == "read_elsewhere" and read_elsewhere(ans)) or \
                (r.get("field") is not None and r["field"][1].search(ans.get(r["field"][0]) or "")):
            return r
    return None


def answers_of(run, keys):
    """{(qid, view): answer row (status, match, fallback)} of the run's answer files for the keys asked."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from idtrack_reproducibility.families import common as C
    out = {}
    if not keys:
        return out
    qids = pa.array(sorted({q for q, _ in keys}), type=pa.string())
    for p in C.answer_files(run):
        t = pq.read_table(p, columns=["qid", "view", "status", "match", "fallback"])
        t = t.filter(pc.is_in(t.column("qid"), value_set=qids))
        for r in t.to_pylist():
            if (r["qid"], r["view"]) in keys:
                out[(r["qid"], r["view"])] = r
    return out


def build(family, scored, m):
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    from idtrack_reproducibility import config as CFG
    rows = []
    for d in scored:
        info = json.load(open(os.path.join(d, "run.json")))
        org = info["organism"]
        rep = {}
        with open(os.path.join(d, "report.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["view"] == "*":
                    rep[r["stratum"]] = int(r["n"])
        wrong = []
        for p in sorted(glob.glob(os.path.join(d, "outcomes", "*.parquet"))):
            t = pq.read_table(p, columns=["qid", "view", "stratum", "unit", "outcome", "status", "match", "reason"])
            t = t.filter(pc.equal(t.column("outcome"), "wrong")).to_pylist()
            wrong += [r for r in t if any(x["rx"].search(r["stratum"]) for x in m)]
        run = os.path.join(CFG.get().root, info["run"]) if info.get("run") else os.path.dirname(d)
        ans = answers_of(run, {(r["qid"], r["view"]) for r in wrong})
        agg = collections.defaultdict(lambda: [0, set(), set()])
        for r in wrong:
            a = dict(ans.get((r["qid"], r["view"])) or {}, status=r["status"], match=r["match"], reason=r["reason"])
            c = classify(m, r["stratum"], a)
            if c is None:
                continue
            x = agg[c["limitation"]]
            x[0] += 1
            x[1].add(r["unit"])
            x[2].add(r["stratum"])
        by_id = {r["limitation"]: r for r in m}
        commit = (info.get("commit") or "")[:12]
        for lid, (n, units, strata) in sorted(agg.items()):
            lim = by_id[lid]
            of = sum(rep.get(s, 0) for s in strata)
            rows.append(dict(limitation=lid, group=lim["group"], title=lim["title"], what=lim["what"], organism=org,
                             answers=n, units=len(units), of_answers=of, share=f"{n / of:.5f}" if of else "",
                             family=family, strata=",".join(sorted(strata)), run=LIM.where(d), commit=commit,
                             graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                             **{"class": lim["class"]}, cause=lim["when"],
                             rule=lim["rule"], user_action=lim["user_action"], fixed_by=lim["fixed_by"],
                             measured_by=f"{family}: outcomes of {LIM.where(d)}", since=commit))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("scored", nargs="+")
    a = ap.parse_args()
    rows = build(a.family, a.scored, load_map(a.map))
    LIM.write(rows, a.out)
    print(f"{len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    sys.exit(main())
