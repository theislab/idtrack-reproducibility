"""Compare two answer files of the corpus (`run_corpus.py` or `run_harmonize.py` output) call by call.

Usage: python -m idtrack_reproducibility.corpus.compare <a.jsonl.gz> <b.jsonl.gz> [--corpus <corpus.jsonl>]
           [--examples 5] [--json <report.json>]

Prints `IDENTICAL (<n> calls)` (exit 0) when every call has the same canonical answer, else `DIFFERENT: <k> of <n>
calls` and the differing calls grouped by entry and class of difference, with examples (exit 1). The classes:

- `missing_in_a` / `missing_in_b`: the call ran in one file only;
- `timeout_a` / `timeout_b` / `timeout_both`: the call hit the cap of `run_corpus.py --timeout` on that side (not an
  answer difference; `timeout_both` is listed so that no call is silently skipped);
- `ok_vs_error`: one side answered, the other raised;
- `error_type`: both raised, different exception types;
- `error_message`: same type, different normalised message;
- `error_message_order_only`: the messages differ only in the order of the items of a printed set or dict (the same
  normalised message, or the same multiset of characters);
- `order_only`: the answers differ only in the order of list items;
- for `convert_identifier`: `target_id:subset` (b drops targets), `target_id:superset` (b adds targets),
  `target_id:changed`; then `value:<field>`, the first top-level field that differs, read in the order `target_id`,
  `no_corresponding`, `no_conversion`, `no_target`, `final_database`, `last_node` and then alphabetically (the answer
  first, then how it was reached);
- for `track_convert`: `candidates_changed` (the set of Ensembl genes reached differs) or `candidate:<fields>` (the
  same genes; the fields that differ, `final_elements_changed` when the targets of a gene differ and
  `final_conversion_details` when only their details do);
- `value:<field>` for other dict answers, `value` for anything else.

Groups are keyed by (entry, class) and, when a corpus is given, list the strata involved. `dt` is ignored.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import sys

from idtrack_reproducibility.corpus.run_corpus import norm_message

CONVERT_FIELDS_FIRST = ["target_id", "no_corresponding", "no_conversion", "no_target", "final_database", "last_node"]


def load(path):
    """(answers by qid, meta) of an answer file; the meta lines (`meta`, `meta_end`) are merged."""
    out, meta = {}, {}
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            r = json.loads(ln)
            if "meta" in r or "meta_end" in r:
                meta.update(r)
                continue
            r.pop("dt", None)
            out[r["qid"]] = r
    return out, meta


def _key(x):
    return json.dumps(x, sort_keys=True)


def unorder(x):
    """The answer with every list sorted (recursively), for order-only detection."""
    if isinstance(x, dict):
        return {k: unorder(v) for k, v in x.items()}
    if isinstance(x, list):
        return sorted((unorder(i) for i in x), key=_key)
    return x


def classify(a, b):
    """The class of difference between two answer records of one call, or None when they are the same answer."""
    ta, tb = "timeout" in a, "timeout" in b
    if ta or tb:
        return "timeout_both" if ta and tb else ("timeout_a" if ta else "timeout_b")
    if a["ok"] != b["ok"]:
        return "ok_vs_error"
    if not a["ok"]:
        ea, eb = a["error"], b["error"]
        if ea["type"] != eb["type"]:
            return "error_type"
        if norm_message(ea["message"]) != norm_message(eb["message"]):
            if sorted(ea["message"]) == sorted(eb["message"]):
                return "error_message_order_only"  # same characters: only the order of items in a printed container
            return "error_message"
        if ea["message"] != eb["message"]:
            return "error_message_order_only"
        return None
    ra, rb = a["result"], b["result"]
    if ra == rb:
        return None
    if unorder(ra) == unorder(rb):
        return "order_only"
    if a["entry"] == "track_convert" and isinstance(ra, dict) and isinstance(rb, dict):
        if set(ra) != set(rb):
            return "candidates_changed"  # the Ensembl genes reached differ
        sub = set()
        for g in ra:
            if ra[g] != rb[g]:
                for k in set(ra[g]) | set(rb[g]):
                    if ra[g].get(k) != rb[g].get(k):
                        if k == "final_conversion":
                            fa, fb = ra[g][k], rb[g][k]
                            if set(fa.get("final_elements", {})) != set(fb.get("final_elements", {})):
                                sub.add("final_elements_changed")
                            else:
                                sub.add("final_conversion_details")
                        else:
                            sub.add(k)
        return "candidate:" + "+".join(sorted(sub))
    if isinstance(ra, dict) and isinstance(rb, dict):
        keys = sorted(set(ra) | set(rb))
        if a["entry"] == "convert_identifier":  # the answer first, then how it was reached
            keys = [k for k in CONVERT_FIELDS_FIRST if k in keys] + [k for k in keys if k not in CONVERT_FIELDS_FIRST]
        for k in keys:
            if ra.get(k) != rb.get(k):
                if k == "target_id" and a["entry"] == "convert_identifier":
                    sa, sb = {_key(i) for i in ra[k]}, {_key(i) for i in rb[k]}
                    if sb < sa:
                        return "target_id:subset"
                    if sb > sa:
                        return "target_id:superset"
                    return "target_id:changed"
                return f"value:{k}"
    return "value"


def compare_records(A: dict, B: dict, corpus: dict | None = None, examples: int = 5) -> dict:
    """The comparison report of two answer sets (qid -> record)."""
    corpus = corpus or {}
    groups = collections.defaultdict(list)
    for qid in sorted(set(A) | set(B)):
        if qid not in A:
            groups[(B[qid]["entry"], "missing_in_a")].append(qid)
            continue
        if qid not in B:
            groups[(A[qid]["entry"], "missing_in_b")].append(qid)
            continue
        kind = classify(A[qid], B[qid])
        if kind:
            groups[(A[qid]["entry"], kind)].append(qid)
    report = {"calls": len(set(A) | set(B)), "differing": sum(len(v) for v in groups.values()), "groups": []}
    for (entry, kind), qids in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        g = {"entry": entry, "kind": kind, "count": len(qids)}
        if corpus:
            g["strata"] = dict(collections.Counter(corpus[q]["stratum"] for q in qids if q in corpus).most_common())
        ex = []
        for q in qids[:examples]:
            e = {"qid": q}
            if q in corpus:
                e["args"] = corpus[q]["args"]
            for side, D in (("a", A), ("b", B)):
                r = D.get(q)
                e[side] = None if r is None else (r.get("result") if r["ok"] else r["error"])
            ex.append(e)
        g["examples"] = ex
        g["qids"] = qids
        report["groups"].append(g)
    return report


def read_corpus(path):
    out = {}
    if path:
        with open(path) as fh:
            for ln in fh:
                q = json.loads(ln)
                out[q["qid"]] = q
    return out


def compare(a: str, b: str, corpus: str | None = None, examples: int = 5) -> dict:
    A, _ = load(a)
    B, _ = load(b)
    report = compare_records(A, B, read_corpus(corpus), examples)
    return {"a": a, "b": b, **report}


def text(report: dict) -> str:
    """The printed form: the first line is `IDENTICAL (<n> calls)` or `DIFFERENT: <k> of <n> calls`."""
    if not report["groups"]:
        return f"IDENTICAL ({report['calls']} calls)\n"
    lines = [f"DIFFERENT: {report['differing']} of {report['calls']} calls"]
    for g in report["groups"]:
        lines.append(f"\n== {g['entry']} / {g['kind']}: {g['count']}" + (f"  strata={g['strata']}" if "strata" in g else ""))
        for e in g["examples"]:
            s = json.dumps(e, sort_keys=True)
            lines.append("    " + s[:1500] + (" ..." if len(s) > 1500 else ""))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="compare", description=__doc__.split("\n\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--corpus", default=None)
    ap.add_argument("--examples", type=int, default=5)
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)
    report = compare(args.a, args.b, args.corpus, args.examples)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    sys.stdout.write(text(report))
    return 0 if not report["groups"] else 1


if __name__ == "__main__":
    sys.exit(main())
