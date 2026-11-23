"""What differs between two graph files, read from their arrays without idtrack (`idtrack-repro graph-judge --against`).

Both graph formats of 2026-11 (d4c2dda's `_narrow` files and the typed api's) store the same arrays for nodes, edges
and connections: `document` (JSON: node types, databases, record), `names` (NUL-separated), `ntype`, the edges
`e_kind` (0 history, 1 cross reference), `e_src`, `e_dst`, `e_ref`, the history releases `h_old`, `h_new`, and the
connections of a cross-reference edge `c_ptr`, `c_db`, `c_asm`, `c_mask` (release masks as 64-bit words). Counted:

- nodes added and removed, by node type; the external ones also by what they are (no letter or digit; UniProt's
  evidence text; the rest);
- cross-reference links `(external node, Ensembl node, database)` with their releases: per database, links in both
  with the same releases, with other releases, only in the old file, only in the new one; a link that left one database
  for another with the same releases is counted as moved, per pair of databases;
- history edges `(source, target, old release, new release)` added and removed.
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict

import numpy as np

ALNUM = re.compile(r"[0-9A-Za-z]")
HISTORY, LINK = 0, 1


def read(path: str) -> dict:
    """The node types, cross-reference links and history edges of a graph file."""
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    return from_arrays(arrays)


def from_arrays(arrays: dict) -> dict:
    document = json.loads(arrays["document"].tobytes().decode("utf-8"))
    names = arrays["names"].tobytes().decode("utf-8").split("\x00")
    node_types, databases = document["node_types"], document["databases"]
    ntype = arrays["ntype"].tolist()
    c_ptr, c_db, c_asm = arrays["c_ptr"].tolist(), arrays["c_db"].tolist(), arrays["c_asm"].tolist()
    masks = [sum(int(w) << (64 * i) for i, w in enumerate(row)) for row in arrays["c_mask"].tolist()]
    connections = [[(databases[c_db[j]], c_asm[j], masks[j]) for j in range(a, b)] for a, b in zip(c_ptr, c_ptr[1:])]
    kind, src, dst, ref = (arrays[k].tolist() for k in ("e_kind", "e_src", "e_dst", "e_ref"))
    h_old, h_new = arrays["h_old"].tolist(), arrays["h_new"].tolist()
    external = node_types.index("external") if "external" in node_types else -1
    links, history = {}, set()
    for e in range(len(kind)):
        if kind[e] == LINK and ntype[src[e]] == external:
            for database, _assembly, mask in connections[ref[e]]:
                key = (names[src[e]], names[dst[e]], database)
                links[key] = links.get(key, 0) | mask
        elif kind[e] == HISTORY:
            history.add((names[src[e]], names[dst[e]], int(h_old[ref[e]]), int(h_new[ref[e]])))
    types = {names[i]: node_types[t] for i, t in enumerate(ntype)}
    return {"types": types, "links": links, "history": history}


def what(name: str, types: dict) -> str:
    """A node's type; an external node's also by what its name is."""
    if types[name] != "external":
        return types[name]
    bare = name[len("synonym_id::"):] if name.startswith("synonym_id::") else name
    if not ALNUM.search(bare):
        return "external: no letter or digit"
    if "{ECO:" in name:
        return "external: UniProt evidence text"
    return "external"


def diff(old: dict, new: dict) -> dict:
    """What differs from the graph `old` to the graph `new` (both as `read` returns them), counted."""
    added, removed = set(new["types"]) - set(old["types"]), set(old["types"]) - set(new["types"])
    out = {"nodes": {"old": len(old["types"]), "new": len(new["types"]),
                     "added": dict(Counter(what(n, new["types"]) for n in added)),
                     "removed": dict(Counter(what(n, old["types"]) for n in removed)),
                     "added_examples": sorted(added)[:30], "removed_examples": sorted(removed)[:30]}}
    before, after = defaultdict(dict), defaultdict(dict)
    for (a, b, db), mask in old["links"].items():
        before[(a, b)][db] = mask
    for (a, b, db), mask in new["links"].items():
        after[(a, b)][db] = mask
    per_db, moved = defaultdict(Counter), Counter()
    for pair in set(before) | set(after):
        x, y = before.get(pair, {}), after.get(pair, {})
        for db in set(x) & set(y):
            per_db[db]["same" if x[db] == y[db] else "other_releases"] += 1
        only_old = {db: m for db, m in x.items() if db not in y}
        only_new = {db: m for db, m in y.items() if db not in x}
        for db_old, mask in sorted(only_old.items()):
            partner = next((db for db, m in sorted(only_new.items()) if m == mask), None)
            if partner is not None:
                moved[(db_old, partner)] += 1
                only_new.pop(partner)
                only_old.pop(db_old)
        for db in only_old:
            per_db[db]["only_old"] += 1
        for db in only_new:
            per_db[db]["only_new"] += 1
    out["links"] = {db: dict(c) for db, c in sorted(per_db.items())}
    out["links_moved_between_databases"] = {f"{a} -> {b}": n for (a, b), n in moved.most_common()}
    gained, lost = new["history"] - old["history"], old["history"] - new["history"]
    out["history"] = {"old": len(old["history"]), "new": len(new["history"]), "n_added": len(gained),
                      "n_removed": len(lost), "added": [list(h) for h in sorted(gained)[:50]],
                      "removed": [list(h) for h in sorted(lost)[:50]]}
    return out
