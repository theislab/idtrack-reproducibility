"""Shared by the graph commands: the places of a tree's graphs, the graph file names, the table farm, the jobs' memory.

    <graphs>/<tree>/<organism>/graph_<organism>_assembly<key>_min<first>_max<last>.npz    the graph (configure it in
                                                                                          [trees.<tree>.graphs])
    <graphs>/<tree>/<organism>/build.json        what built it: tree, commit, time, peak memory, size, content digest
    <graphs>/<tree>/<organism>/judge.json        its judgement (graphs.judge)
    <graphs>/<tree>/<organism>/build/            the local repository of the build (left out of the data manifests)
    <graphs>/<tree>/<organism>/logs/             the jobs' logs

Nothing here imports idtrack.
"""
from __future__ import annotations

import json
import os
import re
import shutil

from idtrack_reproducibility import config as CFG

# A graph file of idtrack: graph_<organism>_assembly<key>_min<first>_max<last>[_narrow|_sessions].npz (`_narrow` in the
# file names of d4c2dda's generation).
GRAPH_NAME = re.compile(r"^graph_(?P<organism>.+)_assembly-?(?P<assembly>\d+)_min(?P<first>\d+)_max(?P<snapshot>\d+)"
                        r"(?:_narrow|_sessions)?\.npz$")
# The files of a local repository of idtrack that are Ensembl's own (downloaded, never derived): the raw tables of the
# core databases (`idtrack_tables/<organism>_assembly-<key>/ens<release>_mysql_<table>__cols=<columns>.h5`) and the
# assembly chain files (`assembly_chain/`). Every other table file of `idtrack_tables/` is derived from them by a build.
RAW_TABLE = re.compile(r"^ens\d+_mysql_")
INDEX_FILE = re.compile(r"^ensembl_https_index_.+\.json$")
# Memory (GB) of a job that builds one graph or judges it: the measured peak of a build and its checks in one job
# (sacct MaxRSS, 2026-11-22, snapshot 115: cynomolgus 1.13, pig 2.10, rhesus 2.14, mouse 5.64, human 10.48 GB) plus a
# quarter, rounded up.
MEMORY_GB = {"macaca_fascicularis": 2, "sus_scrofa": 3, "macaca_mulatta": 3, "mus_musculus": 8, "homo_sapiens": 14}
# Wall time of such a job: the measured elapsed time of those jobs (9 min to 2 h 49 min) twice over.
TIME = {"macaca_fascicularis": "01:00:00", "sus_scrofa": "01:00:00", "macaca_mulatta": "01:00:00",
        "mus_musculus": "02:00:00", "homo_sapiens": "06:00:00"}


def graph_of(path: str) -> dict:
    """Organism, assembly, first and snapshot release of a graph file, read from its name."""
    m = GRAPH_NAME.match(os.path.basename(path))
    if not m:
        raise SystemExit(f"not a graph file name (graph_<organism>_assembly<key>_min<first>_max<last>.npz): {path}")
    return {"organism": m.group("organism"), "assembly": int(m.group("assembly")), "first": int(m.group("first")),
            "snapshot": int(m.group("snapshot"))}


def graph_in(folder: str, organism: str) -> str | None:
    """The one graph file of an organism in a folder, or None when there is none (an error when several)."""
    if not os.path.isdir(folder):
        return None
    found = sorted(f for f in os.listdir(folder) if (m := GRAPH_NAME.match(f)) and m.group("organism") == organism)
    if len(found) > 1:
        raise SystemExit(f"several graph files of {organism} in {folder}: {found}")
    return os.path.join(folder, found[0]) if found else None


def out_dir(tree: str, organism: str) -> str:
    return os.path.join(CFG.get().graphs, tree, organism)


def repository_of(tree: str, organism: str) -> str:
    return os.path.join(out_dir(tree, organism), "build")


def default_tables(organism: str) -> str:
    """The stored tables a build starts from by default: the reference graphs' folder `<graphs>/<organism>/`."""
    return os.path.join(CFG.get().graphs, organism)


def farm(source: str, repository: str) -> dict:
    """Link the raw tables of a local repository `source` into `repository`, file by file; copy its index files.

    Each raw table (RAW_TABLE under `idtrack_tables/`, every file under `assembly_chain/`) becomes a link to the
    resolved file, so that a build reads it and never writes through it (idtrack writes a table it lacks as a new file,
    and never rewrites one it has). The derived tables of `source` are left out: the build derives them again. The index
    files (`ensembl_https_index_*.json`, which a build may rewrite) are copied. A file already in `repository` is kept.
    Returns the counts and the resolved sources with their size and modification time (to check afterwards that none
    changed).
    """
    if not os.path.isdir(source):
        raise SystemExit(f"the tables folder {source} does not exist")
    counts = {"linked": 0, "copied": 0, "left_out": 0, "kept": 0}
    sources = {}
    for base, dirs, files in os.walk(source):
        rel = os.path.relpath(base, source)
        top = rel.split(os.sep)[0]
        for d in dirs:
            if os.path.islink(os.path.join(base, d)):
                raise SystemExit(f"{os.path.join(base, d)}: a linked folder; link its files instead")
        for f in sorted(files):
            path = os.path.join(base, f)
            target = os.path.join(repository, rel, f)
            if rel == "." and INDEX_FILE.match(f):
                kind = "copied"
            elif (top == "idtrack_tables" and RAW_TABLE.match(f)) or top == "assembly_chain":
                kind = "linked"
            else:
                counts["left_out"] += 1
                continue
            os.makedirs(os.path.dirname(target), exist_ok=True)
            real = os.path.realpath(path)
            if os.path.lexists(target):
                counts["kept"] += 1
            elif kind == "copied":
                shutil.copy(real, target)
                counts["copied"] += 1
            else:
                os.symlink(real, target)
                counts["linked"] += 1
            if kind == "linked":
                st = os.stat(real)
                sources[real] = [st.st_size, st.st_mtime_ns]
    return {"counts": counts, "sources": sources}


def changed_sources(sources: dict) -> list[str]:
    """The linked source files whose size or modification time is no longer the recorded one."""
    changed = []
    for real, (size, mtime) in sorted(sources.items()):
        st = os.stat(real)
        if [st.st_size, st.st_mtime_ns] != [size, mtime]:
            changed.append(real)
    return changed


def dump(path: str, obj) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=str)
        fh.write("\n")
    os.replace(path + ".tmp", path)
