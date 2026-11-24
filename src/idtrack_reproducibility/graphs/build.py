"""Build a tree's graphs from stored tables (docs/graphs.md, "Building").

Usage: idtrack-repro graph-build --tree <tree> --organisms <o,...|all> [--snapshot 115] [--tables <folder>]
    [--submit] [--mem <GB>] [--time <hh:mm:ss>]

For each organism, in its own job with `--submit` (the configured runner, `jobs.submit`, named
`graph_build_<tree>_<organism>`), else here one after another:

1. The local repository of the build, `<graphs>/<tree>/<organism>/build/`, is made a farm of the stored tables
   (`graphs.common.farm`): every raw Ensembl table of `--tables` linked (default: the reference graphs' folder
   `<graphs>/<organism>/`; any local repository of idtrack whose builds hold the tables a tree's builder reads), the
   derived tables left out, so that the tree's builder derives them itself. A table the farm lacks is downloaded by the
   builder into the build folder (over HTTPS): move it into the tables folder to keep it.
2. The graph is built in the tree's Python (`[trees.<tree>] python`, its checkout first on the path) through the
   tree's public api: `API(local_repository=...).build_graph(organism_name, snapshot_release)` for the d4c2dda
   generation, `API.build(organism, snapshot, repository)` from 71690ab on.
3. The graph file is linked (a hard link, else a copy) into `<graphs>/<tree>/<organism>/`, and `build.json` records the
   tree, its commit and whether its checkout was modified, the tables folder, the farm's counts, the build's wall
   time and peak memory, the file's size and sha256, the content digest, and whether a linked table changed during the
   build (it must not).

A graph already in the build folder is loaded, not built again, so a job lost to the scheduler is simply submitted
again. A graph file already in `<graphs>/<tree>/<organism>/` that is not the build's own (the same file or the same
bytes) is never replaced and never described by a record of another build: like a published set or a reference, a
stored graph is not rewritten. The command then builds nothing, records nothing, says so and exits 1;
`--keep-existing` keeps such a file the same way and exits 0. Add the printed lines to the configuration (`[trees.<tree>.graphs]`) to ask the graphs; judge them with
`idtrack-repro graph-judge`.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import resource
import shutil
import subprocess
import sys
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import jobs
from idtrack_reproducibility.graphs import common as G


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def content_digest(graph) -> str | None:
    """The content digest a graph records (its record from the typed api's graph format on, else its manifest)."""
    record = getattr(graph, "record", None)
    if record is not None and getattr(record, "content_digest", None):
        return record.content_digest
    manifest = getattr(graph, "manifest", None) or (getattr(graph, "graph", None) or {}).get("manifest") or {}
    return manifest.get("content_digest")


def same_file(a: str | None, b: str | None) -> bool:
    """Whether two paths are one file (a link) or files of the same bytes."""
    if not a or not b or not os.path.exists(a) or not os.path.exists(b):
        return False
    return os.path.samefile(a, b) or (os.path.getsize(a) == os.path.getsize(b) and sha256(a) == sha256(b))


def foreign_graph(present: str | None, own: str | None) -> str | None:
    """Why the graph file `present` of the tree's folder is not the build's own graph `own`, or None."""
    if not present or same_file(present, own):
        return None
    return (f"{present} exists and is not this build's graph"
            + (f" ({own} differs from it)" if own and os.path.exists(own) else "")
            + ": it is not replaced and no record of this build is written for it")


def keep(reason: str, keep_existing: bool) -> None:
    """Leave a foreign graph file as it is: say so; an error unless an existing graph is to be kept."""
    if not keep_existing:
        raise SystemExit(f"{reason}; move it away to build here, or pass --keep-existing to keep it")
    print(f"kept: {reason}")


def build_here(tree: str, organism: str, snapshot: int, tables: str, keep_existing: bool = False) -> dict | None:
    """Build one graph in this process (the tree's Python) and record it; returns the record, or None where a graph
    file of the tree's folder that is not this build's is kept (nothing built, nothing recorded)."""
    from idtrack_reproducibility.families import common as C

    t = CFG.get().tree(tree)
    gen = C.generation(t["source"])
    out, repository = G.out_dir(tree, organism), G.repository_of(tree, organism)
    reason = foreign_graph(G.graph_in(out, organism), G.graph_in(repository, organism))
    if reason:  # decided before hours of building: the file is kept and described by no record of this build
        keep(reason, keep_existing)
        return None
    os.makedirs(repository, exist_ok=True)
    farmed = G.farm(tables, repository)
    sys.path.insert(0, os.path.abspath(t["source"]))
    import idtrack

    if not os.path.realpath(idtrack.__file__).startswith(os.path.realpath(t["source"]) + os.sep):
        raise SystemExit(f"idtrack was imported from {idtrack.__file__}, not from the tree's checkout {t['source']}")
    existed = G.graph_in(repository, organism)
    started = time.perf_counter()
    if gen == "d4c2dda":
        api = idtrack.API(local_repository=repository)
        api.build_graph(organism_name=organism, snapshot_release=snapshot)
        graph = api.track.graph
    else:
        api = idtrack.API.build(organism, snapshot, repository)
        graph = api.graph
    seconds = time.perf_counter() - started
    built = G.graph_in(repository, organism)
    if built is None:
        raise SystemExit(f"the build left no graph file of {organism} in {repository}")
    target = os.path.join(out, os.path.basename(built))
    reason = foreign_graph(target if os.path.exists(target) else None, built)
    if reason:  # a file that appeared during the build: still never described by this build's record
        keep(reason, keep_existing)
        return None
    if not os.path.exists(target):
        try:
            os.link(built, target)
        except OSError:
            shutil.copy2(built, target)
    commit, modified = C.git_commit(t["source"])
    record = {
        "tree": tree, "organism": organism, "snapshot": snapshot, "generation": gen,
        "configured_commit": t["commit"], "commit": commit, "checkout_modified": modified,
        "python": sys.version.split()[0], "host": os.uname().nodename, "slurm_job": os.environ.get("SLURM_JOB_ID"),
        "date": time.strftime("%Y-%m-%d %H:%M:%S"), "tables": os.path.abspath(tables), "farm": farmed["counts"],
        "graph_file": os.path.relpath(target, CFG.get().graphs), "graph_existed": existed is not None,
        "build_s": round(seconds, 1), "peak_GB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024 / 1e9, 3),
        "graph_file_bytes": os.path.getsize(target), "graph_file_sha256": sha256(target),
        "content_digest": content_digest(graph),
        "linked_tables_changed": G.changed_sources(farmed["sources"]),
    }
    G.dump(os.path.join(out, "build.json"), record)
    return record


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro graph-build", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organisms", required=True, help="comma-separated, or `all`")
    ap.add_argument("--snapshot", type=int, default=115)
    ap.add_argument("--tables", default=None,
                    help="a local repository of idtrack holding the stored tables, with `{organism}` for the organism "
                         "(default: the reference graphs' folder <graphs>/<organism>)")
    ap.add_argument("--submit", action="store_true", help="one job per organism through the configured runner")
    ap.add_argument("--mem", type=int, default=None, help="GB of each job (default: graphs.common.MEMORY_GB)")
    ap.add_argument("--time", default=None, help="wall time of each job (default: graphs.common.TIME)")
    ap.add_argument("--keep-existing", action="store_true", help="where the tree's folder holds a graph file that is not "
                    "this build's, keep it, record nothing and exit 0 (default: exit 1)")
    ap.add_argument("--here", action="store_true", help=argparse.SUPPRESS)  # the job itself, in the tree's Python
    a = ap.parse_args(argv)
    organisms = CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")
    t = CFG.get().tree(a.tree)
    for org in organisms:
        tables = a.tables.format(organism=org) if a.tables else G.default_tables(org)
        args = f"--tree {a.tree} --organisms {org} --snapshot {a.snapshot} --tables {os.path.abspath(tables)} --here"
        args += " --keep-existing" if a.keep_existing else ""
        if a.here:
            r = build_here(a.tree, org, a.snapshot, tables, a.keep_existing)
            if r is None:
                continue
            print(f"{org}: {r['graph_file']} content {r['content_digest']} {r['build_s']} s {r['peak_GB']} GB; "
                  f"linked tables changed: {len(r['linked_tables_changed'])}")
            cfg = CFG.get()
            where = os.path.dirname(os.path.abspath(cfg.path)) if cfg.path else os.getcwd()
            print(f'[trees.{a.tree}.graphs] {org} = "{os.path.relpath(os.path.join(cfg.graphs, r["graph_file"]), where)}"')
            if r["linked_tables_changed"]:
                return 1
            continue
        command = jobs.module(t["python"], "idtrack_reproducibility.graphs.build", args)
        logs = os.path.join(G.out_dir(a.tree, org), "logs")
        if a.submit:
            print(org, jobs.submit(f"graph_build_{a.tree}_{org}", command, logs, cpus=2,
                                   mem_gb=a.mem or G.MEMORY_GB.get(org, 16), time=a.time or G.TIME.get(org, "06:00:00")))
        else:
            env = dict(os.environ, **jobs.job_env())
            r = subprocess.run(command, shell=True, env=env)
            if r.returncode:
                return r.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
