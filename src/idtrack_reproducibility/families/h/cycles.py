"""H, measure `cycles`: does a conversion leave reference cycles behind (garbage that only the cyclic collector frees)?

Usage (in the tree's Python; `run.py` starts it): python -m idtrack_reproducibility.families.h.cycles --tree <tree>
    --organism <organism> --out <file.json> [--n 2000]

`--n` corpus identifiers (as `calls`), each converted once with the defaults and once to the symbol database first,
so that the graph's per-node caches hold what they will hold (what remains is what a conversion itself leaves). Then,
with the collector disabled, each identifier is converted again in three ways (`default`; `symbol`: to the symbol
database; `all_explain`: every candidate with explained paths) and `gc.collect()` is called after each conversion:
the number it returns is the number of unreachable objects that conversion left in cycles, and the memory traced
before and after the collection (tracemalloc) is the memory those cycles held. Last, from `gc.garbage` with
`DEBUG_SAVEALL` over 50 conversions, the kinds of objects and the functions in the cycles.
"""

import argparse
import collections
import gc
import json
import os
import tracemalloc

from idtrack_reproducibility.families.h import common as H


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.cycles", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=2000)
    a = ap.parse_args(argv)
    org = a.organism
    report = {"measure": "cycles", "tree": a.tree, "organism": org, **H.environment()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    api = tree.load(tree.link(org, H.run_dir(a.tree, org)))
    report["content_digest"] = tree.manifest(api).get("content_digest")
    ids = H.distinct_ids(org, a.n)
    report["n"] = len(ids)
    symbol = {"final_database": H.SYMBOL_DB[org]}
    for i in ids:
        for kw in ({}, symbol):
            try:
                tree.convert(api, i, **kw)
            except Exception:  # noqa: BLE001 - only the caches matter here
                pass
    gc.collect()
    gc.disable()
    for label, kwargs in (("default", {}), ("symbol", symbol), ("all_explain", {"strategy": "all", "explain": True})):
        counts, freed = [], []
        tracemalloc.start()
        for i in ids:
            try:
                tree.convert(api, i, **kwargs)
            except Exception:  # noqa: BLE001
                pass
            before = tracemalloc.get_traced_memory()[0]
            counts.append(gc.collect())
            freed.append(before - tracemalloc.get_traced_memory()[0])
        tracemalloc.stop()
        counts.sort()
        report[label] = {
            "conversions_leaving_cycles": sum(1 for c in counts if c), "unreachable_objects_total": sum(counts),
            "per_conversion_median": counts[len(counts) // 2], "per_conversion_max": counts[-1],
            "bytes_freed_by_collection_total_MB": round(sum(freed) / 1e6, 2),
            "bytes_freed_max_MB": round(max(freed) / 1e6, 3),
        }
    gc.set_debug(gc.DEBUG_SAVEALL)
    for i in ids[:50]:
        try:
            tree.convert(api, i, strategy="all", explain=True)
        except Exception:  # noqa: BLE001
            pass
    gc.collect()
    kinds = collections.Counter(type(o).__name__ for o in gc.garbage)
    functions = collections.Counter(o.__qualname__ for o in gc.garbage if type(o).__name__ == "function")
    report["garbage_kinds_50_conversions"] = dict(kinds.most_common(12))
    report["garbage_functions_50_conversions"] = dict(functions.most_common(20))
    gc.set_debug(0)
    gc.garbage.clear()
    gc.enable()
    H.dump(a.out, report)
    print(json.dumps(report, default=str))


if __name__ == "__main__":
    main()
