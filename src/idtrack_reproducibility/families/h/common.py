"""Shared by family H's measurements: locations, the inputs drawn from the corpus, memory readings, the api calls of
both generations. Imports nothing of idtrack at module level (the measurement modules import the tree's idtrack after
`load` has put it on sys.path).

Locations: a measurement of tree `t` and organism `o` writes `<performance>/runs/<t>/<o>/<measure>_rep<k>.json`
(measures over all five organisms at once use the organism folder `all`); the harmonization's input datasets are made
once per organism under `<performance>/runs/datasets/<o>/` and read by every tree, so that every tree is given the same
files. The identifiers measured come from the regression corpus (`<corpus>/corpus/<o>.jsonl`).
"""

from __future__ import annotations

import gc
import json
import os
import random
import resource
import sys
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.corpus import _common as CC

ORGANISMS = ("sus_scrofa", "macaca_fascicularis", "macaca_mulatta", "mus_musculus", "homo_sapiens")
SYMBOL_DB = CC.SYMBOL_DB
SEED_IDS = 20261117  # the seed of the identifier samples (distinct_ids) and of the server's request sample


def root() -> str:
    return CFG.get().performance


def run_dir(tree: str, organism: str) -> str:
    return os.path.join(root(), "runs", tree, organism)


def out_file(tree: str, organism: str, measure: str, rep: int) -> str:
    return os.path.join(run_dir(tree, organism), f"{measure}_rep{rep}.json")


def datasets_dir(organism: str) -> str:
    return os.path.join(root(), "runs", "datasets", organism)


def corpus_calls(organism: str, entry: str = "convert_identifier") -> list[dict]:
    with open(CC.corpus_file(organism)) as fh:
        calls = [json.loads(line) for line in fh]
    return [c for c in calls if c["entry"] == entry]


def distinct_ids(organism: str, n: int, seed: int = SEED_IDS) -> list[str]:
    """`n` distinct identifier strings of the corpus's convert_identifier calls (every stratum, blank strings left
    out), seeded, in a fixed order; all of them when the corpus has fewer."""
    ids = sorted({c["args"]["identifier"] for c in corpus_calls(organism) if isinstance(c["args"].get("identifier"), str)})
    ids = [i for i in ids if i.strip()]
    rng = random.Random(seed)
    return rng.sample(ids, min(n, len(ids)))


def statm() -> dict:
    """Resident and shared memory of this process (GB) from /proc."""
    with open("/proc/self/statm") as fh:
        _size, resident, shared = (int(x) for x in fh.read().split()[:3])
    page = os.sysconf("SC_PAGE_SIZE")
    return {"rss_GB": round(resident * page / 1e9, 4), "shared_GB": round(shared * page / 1e9, 4)}


def rss_gb() -> float:
    return statm()["rss_GB"]


def peak_gb() -> float:
    """Peak resident memory of this process (GB; ru_maxrss is in KiB on Linux)."""
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024 / 1e9, 4)


def smaps_rollup(pid="self") -> dict:
    """Rss, Pss, Shared and Private memory (GB) of a process from /proc/<pid>/smaps_rollup (which counts in KiB)."""
    out = {}
    with open(f"/proc/{pid}/smaps_rollup") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 3 and parts[2] == "kB":
                out[parts[0].rstrip(":")] = round(int(parts[1]) * 1024 / 1e9, 4)
    return {k: out.get(k) for k in ("Rss", "Pss", "Shared_Clean", "Shared_Dirty", "Private_Clean", "Private_Dirty")}


class Timer:
    def __init__(self):
        self.t = time.perf_counter()

    def lap(self) -> float:
        now = time.perf_counter()
        dt, self.t = now - self.t, now
        return round(dt, 4)


class GCTimer:
    """Time spent in the cyclic garbage collector (gc.callbacks), per generation, between resets."""

    def __init__(self):
        self.reset()
        self._start = 0.0
        gc.callbacks.append(self._cb)

    def _cb(self, phase, info):
        if phase == "start":
            self._start = time.perf_counter()
        else:
            dt = time.perf_counter() - self._start
            self.s += dt
            self.n[info["generation"]] += 1
            self.s_by_gen[info["generation"]] += dt

    def reset(self):
        self.s, self.n, self.s_by_gen = 0.0, [0, 0, 0], [0.0, 0.0, 0.0]

    def take(self) -> dict:
        out = {"s": round(self.s, 3), "n": list(self.n), "s_by_gen": [round(x, 3) for x in self.s_by_gen]}
        self.reset()
        return out


def pct(times, q):
    """The q-quantile as the element at index int(n * q) of the sorted times (the measurement's convention)."""
    s = sorted(times)
    return s[min(len(s) - 1, int(len(s) * q))]


def latency(times) -> dict:
    """Latency summary of call times (seconds in, milliseconds out)."""
    return {
        "n": len(times), "total_s": round(sum(times), 2),
        "p50_ms": round(pct(times, 0.5) * 1e3, 3), "p90_ms": round(pct(times, 0.9) * 1e3, 2),
        "p99_ms": round(pct(times, 0.99) * 1e3, 1), "p999_ms": round(pct(times, 0.999) * 1e3, 1),
        "max_ms": round(max(times) * 1e3, 1), "mean_ms": round(sum(times) / len(times) * 1e3, 3),
    }


def dump(path: str, obj) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(path + ".tmp", path)


def dump_calls(path: str, rows: list) -> str:
    """The timed calls of a measure as rows of the record (docs/records.md: one row per call, outcome `unjudged`),
    beside its report: `<report without .json>.calls.parquet`. Returns the file."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    out = path[: -len(".json")] + ".calls.parquet" if path.endswith(".json") else path + ".calls.parquet"
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    types = {"record_version": pa.int32(), "seconds": pa.float64(), "calls": pa.int32(), "n_got": pa.int32()}
    cols = sorted({k for r in rows for k in r})
    schema = pa.schema([(c, types.get(c, pa.string())) for c in cols])
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), out + ".tmp", compression="zstd")
    os.replace(out + ".tmp", out)
    return out


MEMORY_UNIT = "GB of 10^9 bytes"


def cpu_model() -> str | None:
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return None


def environment() -> dict:
    """What a measurement records about where it ran (and the unit of its memory figures)."""
    return {"host": os.uname().nodename, "cpu_model": cpu_model(), "python": sys.version.split()[0],
            "slurm_job": os.environ.get("SLURM_JOB_ID"), "cpus": len(os.sched_getaffinity(0)),
            "pythonhashseed": os.environ.get("PYTHONHASHSEED"), "date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "memory_unit": MEMORY_UNIT}


# ---------------------------------------------------------------------------------------------------------------------
# The tree and its api.


class Tree:
    """A configured tree, put first on sys.path, with the calls of its api generation."""

    def __init__(self, name: str):
        self.name = name
        self.cfg = CC.tree(name)
        self.gen = CC.generation(self.cfg["source"])
        CC.use_clone(self.cfg["source"])

    def link(self, organism: str, workdir: str) -> str:
        return CC.link_graph(CC.graph_file(self.name, organism), workdir)

    def load(self, link: str):
        api, _ = CC.load_api(self.cfg["source"], link)
        return api

    def graph(self, api):
        return api.track.graph if self.gen == "d4c2dda" else api.graph

    def manifest(self, api) -> dict:
        g = self.graph(api)
        return dict(getattr(g, "manifest", None) or (getattr(g, "graph", None) or {}).get("manifest") or {})

    # One identifier with the defaults (latest release, Ensembl gene, best answer; the typed api ranks every candidate
    # in the same call).
    def convert(self, api, identifier, **kw):
        if self.gen == "d4c2dda":
            return api.convert_identifier(identifier, **kw)
        kw = dict(kw)
        kw.pop("strategy", None)
        if "go_external" in kw:
            kw["allow_bridges"] = kw.pop("go_external")
        return api.convert(identifier, **kw)

    def convert_many(self, api, identifiers, **kw):
        if self.gen == "d4c2dda":
            return api.convert_identifier_multiple(identifiers, verbose=False, **kw)
        kw = dict(kw)
        kw.pop("strategy", None)
        return api.convert_many(identifiers, **kw)

    def gene_to_database(self, api, gene, database, release):
        if self.gen == "d4c2dda":
            return api.track.gene_to_database(gene, database, release)
        return api.gene_to_database(gene, database, release)

    def plain(self, result):
        """A result as plain data (for equality checks between two ways of asking)."""
        if self.gen == "d4c2dda" or isinstance(result, dict):
            return result
        return result.to_dict()

    def cache_info(self, api) -> dict:
        g = self.graph(api)
        out = {}
        for name in sorted(dir(g)):
            if name.startswith("_c_"):
                f = getattr(g, name, None)
                if hasattr(f, "cache_info"):
                    out[name] = str(f.cache_info())
        return out
