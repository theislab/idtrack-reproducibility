"""Locations, constants and small helpers shared by the benchmark families (the common machinery). Imports nothing of
idtrack.

Every family's code is a package `idtrack_reproducibility.families.<family>`; its data (query sets, runs, logs) lie in
`<families>/<family>/` of the data root and the one truth store is `<truth>/<organism>/` (README, "Configuration").
The locations are read from the configuration when they are used (`BENCH`, `TRUTH`, `HARNESS_QUERIES`, `HARNESS`).
"""
import json
import os
import re
import zlib

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import harness as _harness_pkg

ACCURACY = os.path.dirname(os.path.abspath(_harness_pkg.__file__))  # the harness's code (score.py, harness_lib.py, ...)
ORGANISMS = list(CFG.ORGANISMS)

# Memory per shard (GB): the graph of the reference tree d4c2dda as the harness measured it, with room for the run.
MEM_GB = CFG.MEMORY_GB


def __getattr__(name):
    """The configured locations, read at each use (so a process reads the configuration it was given)."""
    cfg = CFG.get()
    if name == "BENCH":
        return cfg.families
    if name == "TRUTH":
        return cfg.truth
    if name == "HARNESS":
        return cfg.harness
    if name == "HARNESS_QUERIES":
        return os.path.join(cfg.harness, "queries")
    if name == "CACHE":
        return cfg.cache
    raise AttributeError(name)


SEED = 20261118  # cluster bootstrap (docs/methods.md, principle 6)
N_BOOT = 2000

ENS_VERSIONED = re.compile(r"^(ENS[A-Z]*[GTPE]\d{11}|LRG_\d+)\.\d+$")
ENS_ID = re.compile(r"^(ENS[A-Z]*[GTPE]\d{11}|LRG_\d+)")


def generation(clone):
    """`d4c2dda` (API.convert_identifier, two calls for best and all) or `convert` (API.convert, 71690ab on)."""
    if os.path.exists(os.path.join(clone, "idtrack", "api.py")):
        return "convert"
    if os.path.exists(os.path.join(clone, "idtrack", "_api.py")):
        return "d4c2dda"
    raise SystemExit(f"{clone}: neither idtrack/api.py nor idtrack/_api.py; not an idtrack clone")


def python_for(clone):
    """The Python command that runs the tree checked out at `clone` (the configuration's `[trees.<name>] python`)."""
    return CFG.get().python_for(clone)


def shard_of(qid, n):
    return zlib.crc32(qid.encode()) % n


def crc(qid):
    return zlib.crc32(qid.encode())


def base_id(x):
    """An Ensembl id without its version; anything else unchanged."""
    x = str(x)
    m = ENS_VERSIONED.match(x)
    return m.group(1) if m else x


def dumps(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"), default=str)


def family_dir(family):
    return os.path.join(CFG.get().families, family)


def git_commit(clone):
    import subprocess
    try:
        commit = subprocess.run(["git", "-C", clone, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", clone, "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                               text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "", False
    return commit, bool(dirty)


def release_of_graph(graph):
    m = re.search(r"max(\d+)", os.path.basename(os.path.realpath(graph)))
    if not m:
        raise SystemExit(f"{graph}: no `max<release>` in the file name; pass --release")
    return int(m.group(1))
