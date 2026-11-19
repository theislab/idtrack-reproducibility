"""Structural difference of two graph files (e.g. two builds of one organism, or the builds of two commits).

Usage (in the tree's Python): python -m idtrack_reproducibility.corpus.graph_diff --tree <tree to read with>
           --a <graph file> --b <graph file> --out <report.json> [--examples 20]

The tree is a configured one of the d4c2dda generation (it reads the files with its `TheGraph.load`). When the corpus
answers otherwise on two graphs of one organism, this says where the graphs differ.

Reports the manifest fields that differ (equal content digests mean equal content: then nothing else is compared),
nodes only in a / only in b (by node type, with examples), nodes of another type, history edges only in a / only in b
and those with another release pair or weight, cross-reference edges only in a / only in b (by node-type pair) and
those whose connection differs (per database: added, removed, release sets changed), the genes outside the
reference sequence whose releases differ, and the external ids whose spelling for some database and release differs.
"""

from __future__ import annotations

import argparse
import collections

from idtrack_reproducibility.corpus import _common as C

VOLATILE = ("build_date", "package_version", "source_tables")  # manifest fields that differ between any two builds


def content(g):
    """Node types, history edges, cross-reference edges, non-reference marks and spellings of a graph, as plain dicts."""
    types, history, xref, marks, spelled = {}, {}, {}, {}, {}
    for n in g.nodes:
        types[n] = g.node_type(n)
        for e in g.history_edges(n):
            weight = None if e.weight != e.weight else e.weight
            history[(n, e.node, e.key)] = (e.old_release, e.new_release, weight)
        for other, connection in g.xref_edges(n):
            xref[(n, other)] = tuple((db, asm, releases.mask) for db, asm, releases in connection.items)
        releases = g.non_reference_releases(n)
        if releases:
            marks[n] = sorted(releases)
    store = g._frozen()
    for (number, database), (latest, otherwise) in store.spellings.items():
        spelled[(store.names[number], database)] = (latest, tuple(sorted((text, mask) for mask, text in otherwise)))
    return types, history, xref, marks, spelled


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="graph_diff", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True, help="a configured tree of the d4c2dda generation")
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--examples", type=int, default=20)
    args = ap.parse_args(argv)
    t = C.tree(args.tree)
    if C.generation(t["source"]) != "d4c2dda":
        raise SystemExit(f"tree {args.tree}: graph_diff reads graph files with the d4c2dda generation's TheGraph")
    C.use_clone(t["source"])
    from idtrack import TheGraph

    A, B = TheGraph.load(args.a), TheGraph.load(args.b)
    k = args.examples
    rep = {"a": args.a, "b": args.b, "nodes": [A.number_of_nodes(), B.number_of_nodes()],
           "edges": [A.number_of_edges(), B.number_of_edges()]}
    ma, mb = A.manifest, B.manifest
    rep["manifest_diff"] = {f: [str(ma.get(f))[:500], str(mb.get(f))[:500]] for f in sorted(set(ma) | set(mb))
                            if ma.get(f) != mb.get(f) and f not in VOLATILE}
    rep["same_content_digest"] = bool(ma.get("content_digest")) and ma.get("content_digest") == mb.get("content_digest")
    if not rep["same_content_digest"]:
        ta, ha, xa, na, sa = content(A)
        tb, hb, xb, nb, sb = content(B)
        for name, only, types in (("nodes_only_in_a", set(ta) - set(tb), ta), ("nodes_only_in_b", set(tb) - set(ta), tb)):
            by = collections.defaultdict(list)
            for n in sorted(only):
                by[types[n]].append(n)
            rep[name] = {t: {"count": len(v), "examples": v[:k]} for t, v in sorted(by.items())}
        other_type = sorted(n for n in set(ta) & set(tb) if ta[n] != tb[n])
        rep["nodes_of_another_type"] = {"count": len(other_type), "examples": [[n, ta[n], tb[n]] for n in other_type[:k]]}
        for name, da, db in (("history_edges", ha, hb), ("cross_reference_edges", xa, xb)):
            only_a, only_b = sorted(set(da) - set(db)), sorted(set(db) - set(da))
            changed = sorted(e for e in set(da) & set(db) if da[e] != db[e])
            if name == "history_edges":
                rep[name] = {"only_in_a": {"count": len(only_a), "examples": [[*e, *da[e]] for e in only_a[:k]]},
                             "only_in_b": {"count": len(only_b), "examples": [[*e, *db[e]] for e in only_b[:k]]},
                             "changed": {"count": len(changed), "examples": [[*e, da[e], db[e]] for e in changed[:k]]}}
                continue
            pairs = {}
            for side, only, types, data in (("only_in_a", only_a, ta, da), ("only_in_b", only_b, tb, db)):
                by = collections.defaultdict(list)
                for e in only:
                    by[f"{types[e[0]]} -> {types[e[1]]}"].append([*e, str(data[e])[:300]])
                pairs[side] = {t: {"count": len(v), "examples": v[:k]} for t, v in sorted(by.items())}
            kinds, examples = collections.Counter(), collections.defaultdict(list)
            for e in changed:
                ca = {(d, a): m for d, a, m in da[e]}
                cb = {(d, a): m for d, a, m in db[e]}
                for key in sorted(set(ca) | set(cb)):
                    if ca.get(key) == cb.get(key):
                        continue
                    kind = (f"db_removed:{key[0]}" if key not in cb else f"db_added:{key[0]}" if key not in ca
                            else f"db_releases_changed:{key[0]}")
                    kinds[kind] += 1
                    if len(examples[kind]) < 5:
                        examples[kind].append([*e, key[1], bin(ca.get(key, 0)), bin(cb.get(key, 0))])
            pairs["changed"] = {"count": len(changed), "kinds": {kind: {"count": c, "examples": examples[kind]} for kind, c in kinds.most_common()}}
            rep[name] = pairs
        marked = sorted(n for n in set(na) | set(nb) if na.get(n) != nb.get(n))
        rep["non_reference_marks_differ"] = {"count": len(marked), "examples": [[n, na.get(n), nb.get(n)] for n in marked[:k]]}
        spelled = sorted(key for key in set(sa) | set(sb) if sa.get(key) != sb.get(key))
        rep["spellings_differ"] = {"count": len(spelled), "examples": [[*key, str(sa.get(key))[:200], str(sb.get(key))[:200]] for key in spelled[:k]]}
    C.write_json(args.out, rep)

    def short(value):
        if isinstance(value, dict) and "count" in value:
            return value["count"]
        if isinstance(value, dict):
            return {key: short(v) for key, v in value.items()}
        return value

    print({key: short(v) for key, v in rep.items() if key not in ("a", "b")})


if __name__ == "__main__":
    main()
