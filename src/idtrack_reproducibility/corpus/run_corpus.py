"""Run the regression corpus of one organism through one configured tree and write one canonical line per call.

Usage (in the tree's Python; `pipeline.py run` starts it so, as a job or locally):
    python -m idtrack_reproducibility.corpus.run_corpus --tree <tree> --organism <organism> --out <file.jsonl.gz>
        [--corpus <corpus.jsonl>] [--workdir <folder>] [--timeout S] [--shuffle <seed>] [--limit N]
        [--entries convert_identifier,track_convert,...] [--no-network] [--qids <file>]

The graph is the configured one of the tree and organism (`[trees.<tree>.graphs]`), linked into `<workdir>/repo/`
(default: the run folder `<corpus>/runs/<tree>/<organism>/`); loading it writes nothing and needs no network.

Each output line: `{"qid", "entry", "ok", "result" | "error", "dt"}`. `result` is canonical (`canon`): dict keys
sorted, sets as `{"__set__": [...]}` sorted, dicts with non-string keys as `{"__dict__": [[k, v], ...]}` sorted,
floats NaN/inf as strings, numpy scalars as Python values. Lists keep their order: the order of `target_id` and
`last_node` of `convert_identifier` is the engine's ranking and part of the answer. `error` holds the exception type,
its message, `message_norm` (set literals sorted, object addresses removed) and the function that raised. `dt`
(seconds) is ignored by the comparator. A header line `{"meta": ...}` comes first and `{"meta_end": ...}` last.

`--timeout S` caps each call (SIGALRM): a call that runs longer is recorded as `{"ok": false, "timeout": S, "error":
{"type": "CallTimeout"}}`; `pipeline.py reference` reruns such calls with a long cap and merges them. The comparator
reports timeouts as their own kind, never as answer differences.

Two api generations answer the same calls (`_common.generation`):

- `d4c2dda`: every call's arguments are passed as they are (`API.convert_identifier`, `convert_identifier_multiple`
  with `classify_multiple_conversion`, `infer_identifier_source`, the list accessors, `resolve_organism`, the engine's
  `Track.convert`), except `deprioritize_lrg_genes` of `track_convert` lines (rule `lrg_copies`: LRG records follow
  the one rule for copies outside the reference sequence, and the argument is gone).
- `convert` (the typed api, rule `typed_api`): each call is made through the typed api and its result is written in
  the fields of the `d4c2dda` answer (`old_conversion`, `old_engine`): a pure re-arrangement of what the result
  holds, so that the comparator compares answers, not formats. Error types are written as the nearest class that
  `d4c2dda` raised (the typed input errors derive from the built-in errors it raised); messages that name a renamed
  argument name it as `d4c2dda` did (counted in `meta_end.messages_mapped`); every result is checked to survive its
  own JSON round trip (`meta_end.round_trips_checked`).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import random
import re
import signal
import time
import traceback

from idtrack_reproducibility.corpus import _common as C


def _key(x) -> str:
    return json.dumps(x, sort_keys=True, default=str)


def canon(x):
    """Canonical, JSON-serialisable form of an api answer."""
    try:
        import numpy as np

        if isinstance(x, np.generic):
            x = x.item()
        elif isinstance(x, np.ndarray):
            x = x.tolist()
    except ImportError:  # pragma: no cover
        pass
    if isinstance(x, bool) or x is None or isinstance(x, (int, str)):
        return x
    if isinstance(x, float):
        if math.isnan(x):
            return "NaN"
        if math.isinf(x):
            return "Infinity" if x > 0 else "-Infinity"
        return x
    if isinstance(x, dict):
        items = [(canon(k), canon(v)) for k, v in x.items()]
        if all(isinstance(k, str) for k, _ in items):
            return {k: v for k, v in sorted(items, key=lambda kv: kv[0])}
        return {"__dict__": sorted(([k, v] for k, v in items), key=_key)}
    if isinstance(x, (set, frozenset)):
        return {"__set__": sorted((canon(i) for i in x), key=_key)}
    if isinstance(x, (list, tuple)):
        return [canon(i) for i in x]
    if callable(x):
        return f"<callable {getattr(x, '__name__', type(x).__name__)}>"
    return {"__repr__": repr(x)}


_SET_RE = re.compile(r"\{([^{}]*)\}")
_ADDR_RE = re.compile(r" at 0x[0-9a-f]+")


def norm_message(msg: str) -> str:
    """Make messages comparable across processes: sort the items of set literals, drop object addresses."""

    def srt(m):
        inner = m.group(1)
        if re.search(r"'\s*:\s|\d\s*:\s", inner):  # a dict literal ("'k': v"): leave it
            return m.group(0)
        parts = [p.strip() for p in inner.split(",") if p.strip()]
        return "{" + ", ".join(sorted(parts)) + "}"

    return _ADDR_RE.sub("", _SET_RE.sub(srt, msg))


# ---------------------------------------------------------------------------------------------------------------------
# Generation `d4c2dda`: the calls as they are.


def run_one(api, q):
    entry, a = q["entry"], q["args"]
    if entry == "convert_identifier":
        return canon(api.convert_identifier(**a))
    if entry == "track_convert":
        return canon(api.track.convert(**{k: v for k, v in a.items() if k != "deprioritize_lrg_genes"}))
    if entry == "convert_identifier_multiple":
        a = dict(a)
        lst = a.pop("identifier_list")
        res = api.convert_identifier_multiple(lst, verbose=False, **a)
        idx = {id(r): i for i, r in enumerate(res)}
        bins = api.classify_multiple_conversion(res)
        return {
            "results": [canon(r) for r in res],
            "classified": {k: sorted(idx[id(r)] for r in v) for k, v in sorted(bins.items())},
        }
    if entry == "infer_identifier_source":
        return canon(api.infer_identifier_source(**a))
    if entry == "list":
        return canon(getattr(api, a["method"])())
    if entry == "resolve_organism":
        return canon(api.resolve_organism(**a))
    raise ValueError(f"unknown entry {entry}")


# ---------------------------------------------------------------------------------------------------------------------
# Generation `convert`: the typed api, its results written in the fields of the d4c2dda answer.

# Error classes whose names d4c2dda used itself; any other typed error is written as its nearest built-in base.
_OLD_ERROR_NAMES = {"InactiveAtReleaseError", "UnknownVersionError", "SearchBudgetExceededError", "GraphMismatchError"}
# Arguments the typed api renamed (one name per concept); a message is compared with them named as d4c2dda named them.
_RENAMED_ARGUMENTS = {"allow_bridges": "go_external", "bridge_limit": "jump_limit", "identifiers": "identifier_list"}
_RENAMED_ENGINE = {"identifier": "from_id"}
COUNTS = {"round_trips_checked": 0, "messages_mapped": 0}


def error_type(error: BaseException) -> str:
    """The name d4c2dda gave this kind of error: the class itself where d4c2dda had a class of that name, else the
    nearest built-in base."""
    for cls in type(error).__mro__:
        if cls.__name__ in _OLD_ERROR_NAMES or cls.__module__ == "builtins":
            return cls.__name__
    return type(error).__name__


def old_message(msg: str, entry: str) -> str:
    """A message of the typed api with the renamed arguments named as d4c2dda named them."""
    out = msg
    names = dict(_RENAMED_ARGUMENTS)
    if entry == "track_convert":
        names.update(_RENAMED_ENGINE)
    for new, old in names.items():
        out = out.replace(f"`{new}`", f"`{old}`")
    if out != msg:
        COUNTS["messages_mapped"] += 1
    return out


def raw_steps(path, indices):
    """The steps of an explained path as the engine records them: [source, target, key] or [source, target, 0, release]."""
    out = []
    for index in indices:
        step = path.steps[index]
        if step.kind == "identity":
            out.append([None, step.target, None])
        elif step.kind in ("history", "self_loop"):
            out.append([step.source, step.target, step.key])
        else:
            out.append([step.source, step.target, 0, step.release])
    return out


def round_trip(r) -> None:
    """Every result is strict JSON and reads back equal to itself (rule `typed_api`); a failure is a recorded error."""
    from idtrack.results import ConversionResult

    text = json.dumps(r.to_dict(), allow_nan=False)
    if ConversionResult.from_dict(json.loads(text)) != r:
        raise AssertionError(f"to_dict/from_dict changed the result of {r.query_id!r}")
    COUNTS["round_trips_checked"] += 1


def old_conversion(r, strategy: str, explain: bool) -> dict:
    """A typed ConversionResult as d4c2dda's `convert_identifier` returned it for one strategy."""
    round_trip(r)
    cands = r.view(strategy)
    targets = list(dict.fromkeys(c.target for c in cands))
    out = {
        "query_id": r.query_id,
        "graph_id": r.graph_id,
        "match": r.match,
        "known_versions": list(r.known_versions),
        "from_release": r.from_release,
        "from_release_inferred": r.from_release_inferred,
        "from_release_given": r.from_release_given,
        "from_release_fallback": r.from_release_fallback,
        "from_databases": list(r.from_databases),
        "to_release": r.to_release,
        "strategy": strategy,
        "final_database": r.final_database,
        "target_id": targets,
        "last_node": [[c.gene, c.target] for c in cands],
        "ambiguous": len(targets) > 1,
        "tie": r.tie,
        "alternative_locus_genes": list(r.alternative_locus_genes),
        "no_corresponding": r.status in ("not_found", "unknown_version"),
        "no_conversion": r.status == "no_conversion",
        "no_target": r.no_target,
        "no_target_genes": list(r.no_target_genes),
        "error": None if r.error is None else {"type": r.error.type, "message": old_message(r.error.message, "")},
    }
    if explain:
        out["the_path"] = [
            {
                "ensembl_gene": c.gene,
                "target_id": c.target,
                "entry_gene": c.entry_gene,
                "external_jump": c.bridges,
                "external_gap": c.bridge_gap,
                "score": c.score,
                "edges": raw_steps(r.path, c.path + c.final_path),
            }
            for c in cands
        ]
    return out


def old_engine(conv, prioritize: bool, return_path: bool):
    """A typed engine Conversion as d4c2dda's `Track.convert` returned it (None without candidates)."""
    if not conv.candidates:
        return None
    view = [c for c in conv.candidates if c.best] if prioritize else list(conv.candidates)
    result = {}
    for c in view:
        if c.gene not in result:
            entry = {
                "from_id": conv.graph_id,
                "from_release": conv.from_release,
                "from_release_inferred": conv.from_release_inferred,
                "from_databases": list(conv.from_databases),
                "entry_gene": c.entry_gene,
                "entry_support": c.entry_support,
                "to_release": conv.to_release,
                "external_jump": c.bridges,
                "external_gap": c.bridge_gap,
                "external_step": c.external_steps,
                "ensembl_step": c.history_steps,
                "score": c.score,
                "tie": conv.tie,
                "alternative_locus_genes": list(conv.alternative_locus_genes),
                "final_conversion": {
                    "final_database": conv.final_database,
                    "final_conversion_confidence": c.final_confidence,
                    "no_target": conv.no_target,
                    "genes_without_target": list(conv.no_target_genes),
                    "final_elements": {},
                },
            }
            if return_path:
                entry["the_path"] = raw_steps(conv.path, c.path)
            result[c.gene] = entry
        element = {}
        if return_path:
            element["the_path"] = raw_steps(conv.path, c.final_path)
        if prioritize:
            element["filter_scores"] = {
                "initial_filter": [c.final_confidence, c.bridges, c.bridge_gap, c.entry_support, round(c.score, 12)],
                "same_as_input_filter": c.same_as_query,
                "node_importance_filter": None if c.annotation is None else list(c.annotation),
            }
        result[c.gene]["final_conversion"]["final_elements"][c.target] = element
    return result


_LIST = {
    "list_external_databases": lambda api: api.external_databases(),
    "list_genome_assemblies": lambda api: sorted(a.assembly for a in api.describe().assemblies),
    # d4c2dda keyed this by integers; the typed api keys it by strings (strict JSON): read back as d4c2dda keyed it.
    "list_external_databases_by_assembly": lambda api: {int(k): v for k, v in api.databases_by_assembly().items()},
    "external_database_forms": lambda api: api.database_forms(),
    "list_ensembl_releases": lambda api: api.releases(),
}
_MODES = {"ensembl_release": "release", "assembly_ensembl_release": "assembly_release"}


def run_one_typed(api, q):
    entry, a = q["entry"], dict(q["args"])
    if "go_external" in a:
        a["allow_bridges"] = a.pop("go_external")
    if entry == "convert_identifier":
        strategy = a.pop("strategy", "best")
        explain = a.pop("explain", False)
        r = api.convert(a.pop("identifier"), explain=explain, **a)
        return canon(old_conversion(r, strategy, explain))
    if entry == "track_convert":
        a.pop("deprioritize_lrg_genes", None)  # rule `lrg_copies`, as for d4c2dda
        prioritize = a.pop("prioritize_to_one_filter", False)
        return_path = a.pop("return_path", False)
        conv = api.engine.convert(a.pop("from_id"), explain=return_path, **a)
        return canon(old_engine(conv, prioritize, return_path))
    if entry == "convert_identifier_multiple":
        strategy = a.pop("strategy", "best")
        explain = a.pop("explain", False)
        lst = a.pop("identifier_list")
        res = api.convert_many(lst, explain=explain, **a)
        bins = [api.classify([r], strategy) for r in res]
        names = sorted(bins[0]) if bins else []
        return {
            "results": [canon(old_conversion(r, strategy, explain)) for r in res],
            "classified": {k: [i for i, b in enumerate(bins) if b[k]] for k in names},
        }
    if entry == "infer_identifier_source":
        mode = a.get("mode", "assembly_ensembl_release")
        return canon(api.infer_source(a["id_list"], mode=_MODES.get(mode, mode),
                                      report_only_winner=a.get("report_only_winner", True)))
    if entry == "list":
        return canon(_LIST[a["method"]](api))
    if entry == "resolve_organism":
        return canon(api.ensembl_organism(a["tentative_organism_name"]))
    raise ValueError(f"unknown entry {entry}")


# ---------------------------------------------------------------------------------------------------------------------


class CallTimeout(BaseException):  # not swallowed by `except Exception` inside idtrack
    pass


def load_queries(path: str, entries: str = "", no_network: bool = False, qids: str = "", limit: int | None = None,
                 shuffle: int | None = None) -> list[dict]:
    with open(path) as fh:
        queries = [json.loads(ln) for ln in fh]
    if entries:
        keep = set(entries.split(","))
        queries = [q for q in queries if q["entry"] in keep]
    if no_network:
        queries = [q for q in queries if q["entry"] != "resolve_organism"]
    if qids:
        with open(qids) as fh:
            want = {ln.strip() for ln in fh if ln.strip()}
        queries = [q for q in queries if q["qid"] in want]
    if limit:
        queries = queries[:limit]
    if shuffle is not None:
        random.Random(shuffle).shuffle(queries)
    return queries


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="run_corpus", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--corpus", default=None, help="default: <corpus>/corpus/<organism>.jsonl")
    ap.add_argument("--workdir", default=None, help="where the graph is linked (default: the run folder)")
    ap.add_argument("--shuffle", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--entries", default="")
    ap.add_argument("--no-network", action="store_true", help="leave out resolve_organism")
    ap.add_argument("--rss-every", type=int, default=1000)
    ap.add_argument("--timeout", type=int, default=0, help="seconds per call (SIGALRM); 0 = none")
    ap.add_argument("--qids", default="", help="file with one qid per line: run only these")
    args = ap.parse_args(argv)
    os.environ.setdefault("IDTRACK_ENSEMBL_ACCESS", "https")
    t = C.tree(args.tree)
    graph = C.graph_file(args.tree, args.organism)
    organism, snapshot = C.graph_of(graph)
    if organism != args.organism:
        raise SystemExit(f"the configured graph {graph} is of {organism}, not {args.organism}")
    corpus = args.corpus or C.corpus_file(args.organism)
    workdir = args.workdir or C.run_dir(args.tree, args.organism)
    link = C.link_graph(graph, workdir)
    C.use_clone(t["source"])
    queries = load_queries(corpus, args.entries, args.no_network, args.qids, args.limit, args.shuffle)

    def _alarm(signum, frame):
        raise CallTimeout(f"call exceeded {args.timeout} s")

    signal.signal(signal.SIGALRM, _alarm)
    t0 = C.now()
    import idtrack

    t_import = C.now() - t0
    t1 = C.now()
    api, gen = C.load_api(t["source"], link)
    t_load = C.now() - t1
    ask = run_one if gen == "d4c2dda" else run_one_typed
    meta = {
        "tree": args.tree, "commit": t["commit"], "generation": gen, "idtrack_file": idtrack.__file__,
        "graph": os.path.basename(graph), "graph_real": os.path.realpath(graph), "corpus": os.path.abspath(corpus),
        "n": len(queries), "shuffle": args.shuffle, "timeout": args.timeout,
        "pythonhashseed": os.environ.get("PYTHONHASHSEED"), "host": os.uname().nodename,
        "slurm_job": os.environ.get("SLURM_JOB_ID"), "start": time.strftime("%Y-%m-%d %H:%M:%S"),
        "import_s": t_import, "load_s": t_load, "rss_after_load": C.rss_bytes(),
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = args.out + ".partial"
    rss_trace = []
    with gzip.open(tmp, "wt") as out:
        out.write(json.dumps({"meta": meta}) + "\n")
        for i, q in enumerate(queries):
            rec = {"qid": q["qid"], "entry": q["entry"]}
            t_call = time.perf_counter()
            if args.timeout:
                signal.alarm(args.timeout)
            try:
                rec["result"] = ask(api, q)
                rec["ok"] = True
            except CallTimeout:
                rec["ok"] = False
                rec["timeout"] = args.timeout
                rec["error"] = {"type": "CallTimeout", "message": f"exceeded {args.timeout} s", "message_norm": "timeout",
                                "where": ""}
            except Exception as e:  # noqa: BLE001 - every failure is a recorded answer
                msg = str(e) if gen == "d4c2dda" else old_message(str(e), q["entry"])
                rec["ok"] = False
                rec["error"] = {"type": type(e).__name__ if gen == "d4c2dda" else error_type(e), "message": msg,
                                "message_norm": norm_message(msg), "where": traceback.extract_tb(e.__traceback__)[-1].name}
            finally:
                if args.timeout:
                    signal.alarm(0)
            rec["dt"] = time.perf_counter() - t_call
            out.write(json.dumps(rec, sort_keys=True) + "\n")
            if args.rss_every and i % args.rss_every == 0:
                rss_trace.append((i, C.rss_bytes()))
                print(f"progress {i}/{len(queries)} {time.strftime('%H:%M:%S')} rss={rss_trace[-1][1]}", flush=True)
        meta_end = {"end": time.strftime("%Y-%m-%d %H:%M:%S"), "rss_trace": rss_trace, "peak_rss": C.peak_rss_bytes()}
        if gen != "d4c2dda":
            meta_end.update(COUNTS)
        out.write(json.dumps({"meta_end": meta_end}) + "\n")
    os.replace(tmp, args.out)
    print(json.dumps({**meta, **meta_end}, default=str)[:3000])


if __name__ == "__main__":
    main()
