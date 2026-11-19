"""Run the HarmonizeFeatures part of the corpus through one configured tree and write canonical results.

Usage (in the tree's Python, with anndata; `pipeline.py harmonize` starts it once per configuration):
    python -m idtrack_reproducibility.corpus.run_harmonize --tree <tree> --organism <organism> --out <file.jsonl.gz>
        [--workdir <folder>] [--only symbol_latest|ensembl_older|ncbi_latest] [--seed 20261111] [--exclude-slow <run>]

Three small synthetic datasets (4 cells each) are made from the corpus identifiers, deterministically from the seed:
D1 Ensembl gene ids (versioned and not), D2 external ids (symbols, NCBI, RefSeq, UniProt, synonyms, ambiguous
ones), D3 a mix of both with overlaps, non-existent and malformed ids. With `--exclude-slow <run file>`, identifiers
whose `convert_identifier` call took more than 5 s (or timed out) in that run are left out. The datasets are written
once to `<workdir>/datasets/` and reused. Each configuration (target release, final database) runs in a fresh
project folder `<workdir>/run_<tree>_<pid>_<n>/`, removed after it, so no cached matching carries over.

The configurations: `symbol_latest` (the graph's last release, the organism's symbol database), `ensembl_older` (the
third corpus release before the last one, `ensembl_gene`), `ncbi_latest` (the last release, `NCBI gene`, or
`EntrezGene` where the graph has no `NCBI gene`). Results per configuration, one canonical line each (the format of
`run_corpus.py`): the conversion report, `feature_harmonizer` of each dataset with `merge="flag"` and `merge="sum"`
(features before and after, the var table, the column sums of X), and `unify_multiple_anndatas` in union and
intersect mode. Both api generations take the same arguments under their own names (the typed api: `to_release`,
`organism`, `snapshot_release`, `verbosity`, `features`).
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import os
import random
import shutil
import time
import traceback

from idtrack_reproducibility.corpus import _common as C
from idtrack_reproducibility.corpus.run_corpus import canon, error_type, norm_message

CONFIGURATIONS = ("symbol_latest", "ensembl_older", "ncbi_latest")


def slow_identifiers(results_path: str, corpus_path: str, limit_s: float = 5.0) -> set:
    """Identifiers with a convert_identifier call that timed out or took more than `limit_s` in a result file."""
    ids = {}
    with open(corpus_path) as fh:
        for ln in fh:
            q = json.loads(ln)
            if q["entry"] == "convert_identifier":
                ids[q["qid"]] = q["args"].get("identifier")
    slow = set()
    with gzip.open(results_path, "rt") as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("qid") in ids and ("timeout" in r or r.get("dt", 0) > limit_s):
                slow.add(ids[r["qid"]])
    return slow


def dataset_ids(corpus_path: str, seed: int, exclude: set = frozenset()) -> dict[str, list]:
    """The feature names of D1, D2 and D3 (unique, in drawing order)."""
    rng = random.Random(seed)
    by_stratum: dict[str, set] = {}
    with open(corpus_path) as fh:
        for ln in fh:
            q = json.loads(ln)
            if q["entry"] == "convert_identifier" and isinstance(q["args"].get("identifier"), str):
                if q["args"]["identifier"] not in exclude:
                    by_stratum.setdefault(q["stratum"], set()).add(q["args"]["identifier"])

    def pick(s, k):
        pool = sorted(by_stratum.get(s, ()))
        return rng.sample(pool, min(k, len(pool)))

    d1 = (pick("gene_graph_at_release", 500) + pick("gene_table_at_release", 400) + pick("gene_unversioned", 300)
          + pick("gene_split", 150) + pick("gene_merge", 150) + pick("gene_retired", 100))
    d2 = (pick("external_per_database", 900) + pick("synonym_stripped", 150) + pick("external_ambiguous", 200)
          + pick("external_changed_gene", 200))
    d3 = (rng.sample(sorted(set(d1)), 400) + rng.sample(sorted(set(d2)), 400) + pick("nonexistent", 50)
          + [x for x in pick("malformed", 80) if x.strip()])
    return {name: list(dict.fromkeys(ids)) for name, ids in (("D1", d1), ("D2", d2), ("D3", d3))}


def make_datasets(corpus_path: str, outdir: str, seed: int, exclude: set = frozenset()) -> dict[str, str]:
    import anndata as ad
    import numpy as np
    import pandas as pd

    os.makedirs(outdir, exist_ok=True)
    paths = {k: os.path.join(outdir, f"{k}.h5ad") for k in ("D1", "D2", "D3")}
    if all(os.path.exists(p) for p in paths.values()):
        return paths
    for name, ids in dataset_ids(corpus_path, seed, exclude).items():
        n = len(ids)
        X = np.array([[(i * 7 + j * 13) % 11 for i in range(n)] for j in range(4)], dtype=np.float32)
        obs = pd.DataFrame({"sample": [f"{name}_s{j}" for j in range(4)], "age": [20 + j for j in range(4)]},
                           index=[f"{name}_c{j}" for j in range(4)])
        var = pd.DataFrame(index=pd.Index(ids, name=None))
        ad.AnnData(X=X, obs=obs, var=var).write_h5ad(paths[name])
    return paths


def table_summary(table):
    """A data frame as plain rows (missing values as None)."""
    rows = table.astype(object).where(table.notna(), None).values.tolist()
    return {"columns": [str(c) for c in table.columns], "rows": canon(rows)}


def adata_summary(a):
    import numpy as np

    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    return {
        "var_index": [str(i) for i in a.var.index],
        "var": {c: canon([None if (isinstance(v, float) and v != v) else v for v in a.var[c].astype(object).tolist()])
                for c in sorted(a.var.columns)},
        "obs_index": [str(i) for i in a.obs.index],
        "obs": {c: canon([str(v) for v in a.obs[c].tolist()]) for c in sorted(a.obs.columns)},
        "X_colsums": canon([float(x) for x in X.sum(axis=0)]),
        "shape": list(X.shape),
    }


def configurations(organism: str, snapshot: int, summary: dict) -> list[tuple[str, int, str]]:
    older = [r for r in summary["releases"] if r < snapshot][-3]
    return [
        ("symbol_latest", snapshot, C.SYMBOL_DB[organism]),
        ("ensembl_older", older, "ensembl_gene"),
        ("ncbi_latest", snapshot, "NCBI gene" if "NCBI gene" in summary["databases_main"] else "EntrezGene"),
    ]


def harmonizer(gen, datasets, proj, repo, target, fdb, organism, snapshot):
    from idtrack import HarmonizeFeatures

    if gen == "d4c2dda":
        return HarmonizeFeatures(project_name="corpus", data_h5ad_dict=dict(datasets), project_local_repository=proj,
                                 idtrack_local_repository=repo, target_ensembl_release=target, final_database=fdb,
                                 organism_name=organism, graph_last_ensembl_release=snapshot, verbose_level=0)
    return HarmonizeFeatures(project_name="corpus", data_h5ad_dict=dict(datasets), project_local_repository=proj,
                             idtrack_local_repository=repo, to_release=target, final_database=fdb, organism=organism,
                             snapshot_release=snapshot, verbosity=0)


def unify(gen, hf, mode):
    kw = dict(obs_columns_to_keep=["sample", "age"], numeric_obs_columns={"age"})
    if gen == "d4c2dda":
        return hf.unify_multiple_anndatas(mode=mode, **kw)
    return hf.unify_multiple_anndatas(features=mode, **kw)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="run_harmonize", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--corpus", default=None, help="default: <corpus>/corpus/<organism>.jsonl")
    ap.add_argument("--workdir", default=None, help="default: <run folder>/harmonize")
    ap.add_argument("--seed", type=int, default=20261111)
    ap.add_argument("--only", default="", choices=("",) + CONFIGURATIONS, help="run only this configuration")
    ap.add_argument("--exclude-slow", default="", help="a run file: leave out identifiers slower than 5 s there")
    args = ap.parse_args(argv)
    os.environ.setdefault("IDTRACK_ENSEMBL_ACCESS", "https")
    t = C.tree(args.tree)
    graph = C.graph_file(args.tree, args.organism)
    organism, snapshot = C.graph_of(graph)
    corpus = args.corpus or C.corpus_file(args.organism)
    workdir = args.workdir or os.path.join(C.run_dir(args.tree, args.organism), "harmonize")
    link = C.link_graph(graph, workdir)
    repo = os.path.dirname(link)
    gen = C.generation(t["source"])
    C.use_clone(t["source"])
    exclude = slow_identifiers(args.exclude_slow, corpus) if args.exclude_slow else set()
    datasets = make_datasets(corpus, os.path.join(workdir, "datasets"), args.seed, exclude)
    with open(C.summary_file(corpus)) as fh:
        summary = json.load(fh)
    tmp = args.out + ".partial"
    with gzip.open(tmp, "wt") as out:
        out.write(json.dumps({"meta": {"tree": args.tree, "commit": t["commit"], "generation": gen,
                                       "graph": os.path.basename(graph), "datasets": datasets, "seed": args.seed,
                                       "excluded": len(exclude), "start": time.strftime("%Y-%m-%d %H:%M:%S")}}) + "\n")
        for i, (name, target, fdb) in enumerate(configurations(organism, snapshot, summary)):
            if args.only and name != args.only:
                continue
            proj = os.path.join(workdir, f"run_{args.tree}_{os.getpid()}_{i}")
            shutil.rmtree(proj, ignore_errors=True)
            os.makedirs(proj)
            steps = []
            t0 = time.perf_counter()
            try:
                hf = harmonizer(gen, datasets, proj, repo, target, fdb, organism, snapshot)
                steps.append(("conversion_report", lambda: table_summary(hf.conversion_report())))
                for d in sorted(datasets):
                    for merge in ("flag", "sum"):
                        def fh_step(d=d, merge=merge):
                            r = hf.feature_harmonizer(d, merge=merge)
                            return {"t0": r[1], "t1": r[2], "adata": adata_summary(r[0])}
                        steps.append((f"feature_harmonizer_{merge}_{d}", fh_step))
                for mode in ("union", "intersect"):
                    steps.append((f"unify_{mode}", lambda mode=mode: adata_summary(unify(gen, hf, mode))))
            except Exception as e:  # noqa: BLE001 - the constructor's failure is the answer of every step
                steps = [("init", lambda e=e: (_ for _ in ()).throw(e))]
            for sname, fn in steps:
                rec = {"qid": f"{organism}-harmonize-{name}-{sname}", "entry": "harmonize"}
                t1 = time.perf_counter()
                try:
                    rec["result"] = canon(fn())
                    rec["ok"] = True
                except Exception as e:  # noqa: BLE001
                    rec["ok"] = False
                    rec["error"] = {"type": type(e).__name__ if gen == "d4c2dda" else error_type(e), "message": str(e),
                                    "message_norm": norm_message(str(e)),
                                    "where": traceback.extract_tb(e.__traceback__)[-1].name}
                rec["dt"] = time.perf_counter() - t1
                out.write(json.dumps(rec, sort_keys=True) + "\n")
            print(name, "done", round(time.perf_counter() - t0, 1), flush=True)
            # Each HarmonizeFeatures instance loads its own graph: free it before the next configuration.
            steps = None
            hf = None
            gc.collect()
            shutil.rmtree(proj, ignore_errors=True)
    os.replace(tmp, args.out)


if __name__ == "__main__":
    main()
