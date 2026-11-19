"""H, measure `load`: import and load one organism's graph through a tree's public api; time and memory.

Usage (in the tree's Python; `run.py` starts it): python -m idtrack_reproducibility.families.h.load --tree <tree>
    --organism <organism> --out <file.json>

In one fresh process:
1. `import idtrack`, then the attribute `idtrack.API` (which imports the conversion core), each timed, resident memory
   after each, the modules loaded and which of the heavy ones (pandas, anndata, scipy, h5py, networkx, ...).
2. The load (`API(local_repository=...).build_graph(...)` for the d4c2dda generation, `API.load(path)` for the typed
   api), timed; resident and peak memory after it; nodes, edges (history + cross-reference) and the content digest
   from the graph's manifest; the files the load wrote beside the graph (expected: none).
3. The objects the cyclic collector tracks and one `gc.collect()`, timed.
4. For the d4c2dda generation only (it reads the graph's internal store): the memory of the store by part (numpy
   arrays; names; the name index and its int objects; the type list; connections; spellings) and the steps of
   `TheGraph.load` replayed one by one, each timed; `versions_of` on an unknown name, first and second call.
5. The graph file's size and its npz members (size, compression).
"""

import argparse
import gc
import json
import os
import sys
import zipfile

from idtrack_reproducibility.families.h import common as H

HEAVY = ("pandas", "anndata", "scipy", "h5py", "torch", "networkx", "numpy", "requests", "pymysql", "tqdm", "yaml")


def store_accounting(graph, path, report):
    """d4c2dda's store, part by part, and the steps of its load replayed (private structures of that generation)."""
    import numpy as np

    from idtrack import _the_graph as TG
    from idtrack._the_graph import Connection, ReleaseSet

    store = graph._store
    arrays = {k: v for k, v in vars(store).items() if isinstance(v, np.ndarray)}
    report["store_arrays_total_MB"] = round(sum(v.nbytes for v in arrays.values()) / 1e6, 1)
    report["names_list_MB"] = round((sys.getsizeof(store.names) + sum(sys.getsizeof(s) for s in store.names)) / 1e6, 1)
    report["index_dict_MB"] = round(sys.getsizeof(store.index) / 1e6, 1)
    report["index_int_objects_MB"] = round(sum(sys.getsizeof(i) for i in store.index.values() if i > 256) / 1e6, 1)
    report["ntype_list_MB"] = round(sys.getsizeof(store.ntype_list) / 1e6, 1)

    def deep(obj, seen):
        if id(obj) in seen:
            return 0
        seen.add(id(obj))
        size = sys.getsizeof(obj)
        if isinstance(obj, dict):
            size += sum(deep(k, seen) + deep(v, seen) for k, v in obj.items())
        elif isinstance(obj, (list, tuple, set, frozenset)):
            size += sum(deep(x, seen) for x in obj)
        elif isinstance(obj, Connection):
            size += deep(obj.items, seen) + deep(obj.releases, seen) + deep(obj.as_dict, seen) + deep(obj.databases, seen)
        return size

    report["connections_MB"] = round(deep(store.connections, set()) / 1e6, 1)
    report["spellings_MB"] = round(deep(store.spellings, set()) / 1e6, 1)
    t = H.Timer()
    graph.versions_of("NOT_AN_ID_X")
    report["versions_of_first_call_s"] = t.lap()
    graph.versions_of("NOT_AN_ID_Y")
    report["versions_of_second_call_s"] = t.lap()
    ReleaseSet._interned.clear()
    Connection._interned.clear()
    steps = {}
    t.lap()
    with np.load(path, allow_pickle=False) as archive:
        arrs = {name: archive[name] for name in archive.files}
    steps["read_npz_members"] = t.lap()
    document = json.loads(arrs["document"].tobytes().decode("utf-8"))
    steps["json_document"] = t.lap()
    attributes = TG._from_json_graph(document["graph"])
    steps["graph_attributes"] = t.lap()
    names = arrs["names"].tobytes().decode("utf-8").split(TG._NAME_SEPARATOR)
    steps["names_decode_split"] = t.lap()
    databases = document["databases"]
    c_ptr, c_db, c_asm = arrs["c_ptr"].tolist(), arrs["c_db"].tolist(), arrs["c_asm"].tolist()
    mask_values = TG._words_to_masks(arrs["c_mask"])
    steps["words_to_masks"] = t.lap()
    connections = [
        TG.Connection.of((databases[c_db[j]], c_asm[j], TG.ReleaseSet.from_mask(mask_values[j])) for j in range(x, y))
        for x, y in zip(c_ptr, c_ptr[1:])
    ]
    steps["connections"] = t.lap()
    s2 = TG._Store(names, document["node_types"], arrs["ntype"], arrs["e_src"], arrs["e_dst"], arrs["e_key"],
                   arrs["e_kind"], arrs["e_ref"], arrs["in_order"], arrs["h_old"], arrs["h_new"], arrs["h_weight"],
                   connections)
    steps["store_init"] = t.lap()
    s2.globals = TG._from_json_globals(document["globals"])
    steps["globals"] = t.lap()
    s2.spellings = {(int(i), db): (latest, tuple((int(mask, 16), text) for mask, text in otherwise))
                    for i, db, latest, otherwise in document["spellings"]}
    steps["spellings"] = t.lap()
    marked = arrs["nonref_node"].tolist()
    s2.non_reference = dict(zip(marked, (TG.ReleaseSet.from_mask(m) for m in TG._words_to_masks(arrs["nonref_mask"]))))
    steps["non_reference"] = t.lap()
    g2 = TG.TheGraph()
    g2.graph.update(attributes)
    g2._install(s2)
    steps["install"] = t.lap()
    report["load_steps_s"] = steps
    t.lap()
    dict(zip(names, range(len(names))))
    report["index_dict_build_s"] = t.lap()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.load", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    report = {"measure": "load", "tree": a.tree, "organism": a.organism, **H.environment(), "rss_start": H.statm()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    link = tree.link(a.organism, H.run_dir(a.tree, a.organism))
    t = H.Timer()
    import idtrack

    report["import_idtrack_s"] = t.lap()
    report["rss_after_import_idtrack"] = H.statm()
    idtrack.API  # noqa: B018 - the attribute imports the conversion core
    report["import_api_s"] = t.lap()
    report["rss_after_import_api"] = H.statm()
    report["modules_loaded"] = len(sys.modules)
    report["heavy_modules_loaded"] = sorted(m for m in HEAVY if m in sys.modules)
    report["idtrack_file"] = idtrack.__file__
    repo = os.path.dirname(link)
    before = set(os.listdir(repo))
    t.lap()
    api = tree.load(link)
    report["load_s"] = t.lap()
    report["rss_after_load"] = H.statm()
    report["peak_after_load_GB"] = H.peak_gb()
    report["files_written_beside_graph"] = sorted(set(os.listdir(repo)) - before)
    report["modules_after_load"] = sorted(m for m in HEAVY if m in sys.modules)
    manifest = tree.manifest(api)
    counts = manifest.get("counts") or {}
    report["nodes"] = counts.get("nodes")
    report["edges"] = (counts["history_edges"] + counts["cross_reference_edges"]
                       if "history_edges" in counts and "cross_reference_edges" in counts else None)
    report["content_digest"] = manifest.get("content_digest")
    report["gc_tracked_objects"] = len(gc.get_objects())
    t.lap()
    gc.collect()
    report["gc_collect_s"] = t.lap()
    path = os.path.realpath(link)
    report["graph_file_MB"] = round(os.path.getsize(path) / 1e6, 1)
    with zipfile.ZipFile(path) as z:
        report["npz_members"] = {i.filename: {"MB": round(i.file_size / 1e6, 2), "method": i.compress_type}
                                 for i in z.infolist()}
    if tree.gen == "d4c2dda":
        store_accounting(tree.graph(api), path, report)
    report["rss_end"] = H.statm()
    report["peak_end_GB"] = H.peak_gb()
    H.dump(a.out, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ("npz_members",)}, default=str)[:3000])


if __name__ == "__main__":
    main()
