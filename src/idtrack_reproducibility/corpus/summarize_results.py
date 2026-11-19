"""Summarise one answer file of the corpus: answers and exceptions by entry and stratum, seconds by entry.

Usage: python -m idtrack_reproducibility.corpus.summarize_results <run.jsonl.gz> --corpus <corpus.jsonl>
           [--json <out.json>]   (or `python -m idtrack_reproducibility.corpus summarize --tree <t> --organism <o>`)
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="summarize_results")
    ap.add_argument("results")
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)
    with open(args.corpus) as fh:
        corpus = {q["qid"]: q for q in map(json.loads, fh)}
    by_entry = collections.Counter()
    errors = collections.Counter()
    err_strata = collections.defaultdict(collections.Counter)
    err_examples = {}
    flags = collections.Counter()
    ntargets = collections.Counter()
    times = collections.defaultdict(float)
    meta = {}
    with gzip.open(args.results, "rt") as fh:
        for ln in fh:
            r = json.loads(ln)
            if "meta" in r or "meta_end" in r:
                meta.update(r)
                continue
            q = corpus[r["qid"]]
            by_entry[(r["entry"], r["ok"])] += 1
            times[r["entry"]] += r.get("dt", 0)
            if not r["ok"]:
                k = (r["entry"], r["error"]["type"], r["error"]["message_norm"][:160])
                errors[k] += 1
                err_strata[k][q["stratum"]] += 1
                err_examples.setdefault(k, {"qid": r["qid"], "args": q["args"], "where": r["error"].get("where")})
            elif r["entry"] == "convert_identifier":
                res = r["result"]
                f = "no_corresponding" if res["no_corresponding"] else "no_conversion" if res["no_conversion"] else (
                    "no_target" if res["no_target"] else "converted")
                flags[f] += 1
                ntargets[min(len(res["target_id"]), 3)] += 1
    out = {
        "meta": meta,
        "calls_by_entry_ok": {f"{e}|ok={o}": c for (e, o), c in sorted(by_entry.items())},
        "seconds_by_entry": dict(times),
        "convert_identifier_outcomes": dict(flags),
        "convert_identifier_n_targets(3=3+)": dict(sorted(ntargets.items())),
        "errors": [{"entry": e, "type": t, "message": m, "count": c, "strata": dict(err_strata[(e, t, m)]),
                    "example": err_examples[(e, t, m)]} for (e, t, m), c in errors.most_common()],
    }
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "meta"}, indent=1)[:12000])


if __name__ == "__main__":
    main()
