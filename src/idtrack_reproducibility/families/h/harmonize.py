"""H, measure `harmonize`: one joint harmonization of three datasets on a real graph, timed by step, with the time of
the cyclic garbage collector and the peak memory.

Usage (in the tree's Python, with anndata; `run.py` starts it): python -m idtrack_reproducibility.families.h.harmonize
    --tree <tree> --organism <organism> --out <file.json> [--freeze]

The datasets (50 cells, sparse counts, seed 7), made once per organism from the graph and read by every tree
(`<performance>/runs/datasets/<organism>/`): `d115`, every gene id active at the graph's last release, without
version; `dold`, every gene id active at release 90 (or the graph's first release when later), without version;
`s115`, the symbol of every gene of the last release that has one in the organism's symbol database. The run: the
constructor of `HarmonizeFeatures` over the three, to the symbol database at the last release, on an empty project
folder (no result cache), then `unify_multiple_anndatas`; then the constructor again on the result cache it wrote.
Recorded: the time of each, the collector's time and collections in each, the features in and out, the status counts
per dataset, the sizes of the result cache, the peak resident memory, and digests of every output (merged `var` and
its dtypes, `obs`, `X`, the conversion report, the cache files), so that two trees can be shown to give the same
output. `--freeze` calls `gc.freeze()` after the load (what a server may do) to show what the collector costs.
"""

import argparse
import gc
import hashlib
import json
import os
import shutil

from idtrack_reproducibility.corpus import _common as CC
from idtrack_reproducibility.families.h import common as H


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


def dataset_names(tree, api, organism, snapshot):
    """The feature names of d115, dold and s115, read from the graph through the tree."""
    graph = tree.graph(api)
    symbol_db = H.SYMBOL_DB[organism]
    if tree.gen == "d4c2dda":
        from idtrack._db import DB
        alternative = set(DB.alternative_versions)
    else:
        from idtrack import ids
        alternative = {ids.VOID, ids.RETIRED}
    old_release = max(90, graph.first_release)

    def genes_at(release):
        return sorted(n for n in graph.nodes
                      if graph.is_backbone(n) and graph.version(n) not in alternative and graph.is_active(n, release))

    genes_last = genes_at(snapshot)
    return {
        "d115": sorted({graph.base_id(n) for n in genes_last}),
        "dold": sorted({graph.base_id(n) for n in genes_at(old_release)}),
        "s115": sorted({found[0] for found in (tree.gene_to_database(api, n, symbol_db, snapshot)[0] for n in genes_last)
                        if found}),
    }


def make_datasets(tree, api, organism, snapshot, folder):
    """The three h5ad files (made by the first tree that needs them); returns their paths and whether made now."""
    paths = {name: os.path.join(folder, f"{name}.h5ad") for name in ("d115", "dold", "s115")}
    if all(os.path.exists(p) for p in paths.values()):
        return paths, False
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(7)
    for name, index in dataset_names(tree, api, organism, snapshot).items():
        x = sp.random(50, len(index), density=0.3, format="csr", random_state=rng, dtype=np.float64)
        x.data = rng.poisson(1.0, size=x.data.shape).astype(np.int32) + 1
        adata = ad.AnnData(X=x, obs=pd.DataFrame(index=[f"c{i}" for i in range(50)]), var=pd.DataFrame(index=index))
        adata.write_h5ad(paths[name] + ".tmp.h5ad")
        os.replace(paths[name] + ".tmp.h5ad", paths[name])
    return paths, True


def harmonizer(tree, idtrack, paths, project, folder, snapshot, symbol_db, api):
    if tree.gen == "d4c2dda":
        return idtrack.HarmonizeFeatures("p", dict(paths), project, folder, snapshot, final_database=symbol_db, api=api,
                                         verbose_level=0)
    from idtrack.harmonize import HarmonizeFeatures
    return HarmonizeFeatures("p", dict(paths), project, folder, snapshot, final_database=symbol_db, api=api, verbosity=0)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.harmonize", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--freeze", action="store_true")
    a = ap.parse_args(argv)
    org = a.organism
    report = {"measure": "harmonize", "tree": a.tree, "organism": org, "freeze": a.freeze, **H.environment()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    work = os.path.join(H.run_dir(a.tree, org), "harmonize")
    link = tree.link(org, work)
    folder = os.path.dirname(link)
    _, snapshot = CC.graph_of(link)
    symbol_db = H.SYMBOL_DB[org]
    import idtrack

    t = H.Timer()
    api = tree.load(link)
    report["load_s"] = t.lap()
    report["content_digest"] = tree.manifest(api).get("content_digest")
    paths, made = make_datasets(tree, api, org, snapshot, H.datasets_dir(org))
    report["datasets_made_here"] = made
    report["datasets_s"] = t.lap()
    report["datasets_sha256"] = {k: CC.sha256(p) for k, p in sorted(paths.items())}
    gc.collect()
    if a.freeze:
        gc.freeze()
    report["gc_tracked_before_run"] = len(gc.get_objects())
    report["rss_before_run"] = H.statm()
    gct = H.GCTimer()
    project = os.path.join(work, "joint")
    shutil.rmtree(project, ignore_errors=True)
    t.lap()
    h = harmonizer(tree, idtrack, paths, project, folder, snapshot, symbol_db, api)
    report["constructor_s"] = t.lap()
    report["constructor_gc"] = gct.take()
    merged = h.unify_multiple_anndatas()
    report["unify_s"] = t.lap()
    report["unify_gc"] = gct.take()
    report["features_in"] = {}
    report["features_out"] = int(merged.n_vars)
    import scipy.sparse as sp

    X = merged.X.tocsr() if sp.issparse(merged.X) else sp.csr_matrix(merged.X)
    report["digests"] = {
        "var_csv": _sha(merged.var.to_csv().encode()),
        "var_dtypes": _sha(str(merged.var.dtypes.to_dict()).encode()),
        "obs_csv": _sha(merged.obs.to_csv().encode()),
        "X": _sha(X.data.tobytes() + X.indices.tobytes() + X.indptr.tobytes() + str(X.dtype).encode()),
        "conversion_report_csv": _sha(h.conversion_report().to_csv().encode()),
        "cache_files": {f: _sha(open(os.path.join(project, f), "rb").read()) for f in sorted(os.listdir(project))},
    }
    report["statuses"] = {}
    for d in paths:
        table = h.feature_table(d)
        report["features_in"][d] = int(len(table))
        report["statuses"][d] = {str(k): int(v) for k, v in table["status"].value_counts().items()}
    report["result_cache_MB"] = {f: round(os.path.getsize(os.path.join(project, f)) / 1e6, 2)
                                 for f in sorted(os.listdir(project))}
    del merged, X
    t.lap()
    harmonizer(tree, idtrack, paths, project, folder, snapshot, symbol_db, api)
    report["constructor_from_cache_s"] = t.lap()
    report["from_cache_gc"] = gct.take()
    report["rss_end"] = H.statm()
    report["peak_GB"] = H.peak_gb()
    shutil.rmtree(project, ignore_errors=True)
    H.dump(a.out, report)
    print(json.dumps(report, default=str)[:4000])


if __name__ == "__main__":
    main()
