"""Shared by the corpus modules: the configured locations, the trees and their graphs, importing a tree, memory.

Locations (docs/corpus.md, "Files"): everything lies under the configured `<corpus>` location (`[data] corpus`):

    corpus/<organism>.jsonl, corpus/<organism>.summary.json      the corpus (fixed once published)
    reference/<organism>.<tree>.jsonl.gz                          a tree's reference answers to the corpus
    reference/<organism>.<tree>.harmonize.jsonl.gz                its answers to the HarmonizeFeatures part
    runs/<tree>/<organism>/                                       every run, comparison and log

Nothing here imports idtrack: the tree is put on `sys.path` by `use_clone` in the process that asks it.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import resource
import sys
import time

from idtrack_reproducibility import config as CFG

# A graph file of idtrack (narrow, all forms): graph_<organism>_assembly<key>_min<first>_max<last>_narrow.npz
GRAPH_NAME = re.compile(r"^graph_(?P<organism>.+)_assembly-?\d+_min\d+_max(?P<snapshot>\d+)_narrow\.npz$")

# The symbol database of each organism (the harmonize part converts to it).
SYMBOL_DB = {
    "homo_sapiens": "HGNC Symbol",
    "mus_musculus": "MGI Symbol",
    "sus_scrofa": "VGNC Symbol",
    "macaca_mulatta": "VGNC Symbol",
    "macaca_fascicularis": "VGNC Symbol",
}


def root() -> str:
    """The configured `<corpus>` location."""
    return CFG.get().corpus


def corpus_file(organism: str) -> str:
    return os.path.join(root(), "corpus", f"{organism}.jsonl")


def summary_file(corpus_path: str) -> str:
    return corpus_path[: -len(".jsonl")] + ".summary.json" if corpus_path.endswith(".jsonl") else corpus_path + ".summary.json"


def reference_file(organism: str, tree: str, harmonize: bool = False) -> str:
    return os.path.join(root(), "reference", f"{organism}.{tree}{'.harmonize' if harmonize else ''}.jsonl.gz")


def run_dir(tree: str, organism: str) -> str:
    return os.path.join(root(), "runs", tree, organism)


def tree(name: str) -> dict:
    """The configured tree (`[trees.<name>]`: commit, source, python, graphs)."""
    t = CFG.get().tree(name)
    if not t["source"]:
        raise SystemExit(f"tree {name!r} has no `source` (a checkout of idtrack at its commit) in the configuration")
    return t


def graph_file(tree_name: str, organism: str) -> str:
    return CFG.get().graph(tree_name, organism)


def graph_of(path: str) -> tuple[str, int]:
    """Organism and snapshot release of a graph file, read from its name."""
    m = GRAPH_NAME.match(os.path.basename(path))
    if not m:
        raise SystemExit(f"not a graph file name (graph_<organism>_assembly<key>_min<first>_max<last>_narrow.npz): {path}")
    return m.group("organism"), int(m.group("snapshot"))


def graph_in(repo: str, organism: str) -> str:
    """The one graph file of an organism in a local repository folder."""
    found = sorted(f for f in os.listdir(repo) if (m := GRAPH_NAME.match(f)) and m.group("organism") == organism)
    if len(found) != 1:
        raise SystemExit(f"expected one graph file of {organism} in {repo}, found {found}")
    return os.path.join(repo, found[0])


def generation(clone: str) -> str:
    """`d4c2dda` (`API.convert_identifier` and the engine's `Track.convert`) or `convert` (the typed api: `API.convert`,
    from 71690ab on): read from the checkout's files, so that nothing of idtrack is imported to decide it."""
    if os.path.exists(os.path.join(clone, "idtrack", "api.py")):
        return "convert"
    if os.path.exists(os.path.join(clone, "idtrack", "_api.py")):
        return "d4c2dda"
    raise SystemExit(f"{clone}: neither idtrack/api.py nor idtrack/_api.py; not an idtrack checkout")


def use_clone(clone: str) -> None:
    """Put the checkout first on sys.path and make sure no other idtrack is importable before it."""
    clone = os.path.abspath(clone)
    if not os.path.isdir(os.path.join(clone, "idtrack")):
        raise SystemExit(f"not an idtrack checkout: {clone}")
    sys.path.insert(0, clone)
    for mod in [m for m in sys.modules if m == "idtrack" or m.startswith("idtrack.")]:
        del sys.modules[mod]


def link_graph(graph: str, workdir: str) -> str:
    """Link the graph file into `<workdir>/repo/` (nothing is ever written beside the original); return the link."""
    repo = os.path.join(workdir, "repo")
    os.makedirs(repo, exist_ok=True)
    link = os.path.join(repo, os.path.basename(graph))
    try:  # two processes of one run start together: the first makes the link
        os.symlink(os.path.abspath(graph), link)
    except FileExistsError:
        pass
    if os.path.realpath(link) != os.path.realpath(graph):
        raise SystemExit(f"{link} points to {os.path.realpath(link)}, not to {graph}: use another work folder")
    return link


def load_api(clone: str, graph_link: str):
    """Load the graph through the tree's public api; returns (api, generation). The tree must be on sys.path."""
    gen = generation(clone)
    from idtrack import API

    if gen == "d4c2dda":
        organism, snapshot = graph_of(graph_link)
        api = API(local_repository=os.path.dirname(graph_link))
        api.build_graph(organism_name=organism, snapshot_release=snapshot)
    else:
        api = API.load(graph_link)
    return api, gen


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rss_bytes() -> int:
    """Current resident set size of this process."""
    with open("/proc/self/statm") as fh:
        pages = int(fh.read().split()[1])
    return pages * os.sysconf("SC_PAGE_SIZE")


def peak_rss_bytes() -> int:
    """Peak resident set size of this process (ru_maxrss is in KiB on Linux)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024


def write_json(path: str, obj) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(tmp, path)


def now() -> float:
    return time.perf_counter()
