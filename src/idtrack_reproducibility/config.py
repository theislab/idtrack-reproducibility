"""The one configuration of the benchmarks: where the data are, which idtrack trees are measured, how jobs run.

Read from a TOML file and environment variables (README, "Configuration"). The file is the first of:
`$IDTRACK_REPRO_CONFIG`, `./idtrack-repro.toml`; without one, every default below applies. Relative paths in the
file are read from the file's folder; relative paths in variables from the working directory.

Locations (each a folder; defaults under the data root):

    [data]
    root = "data"                  # IDTRACK_REPRO_DATA_ROOT
    truth = "<root>/truth"         # IDTRACK_REPRO_TRUTH        Ensembl's and the authorities' files, per organism
    harness = "<root>/harness"     # IDTRACK_REPRO_HARNESS      the accuracy harness: queries/<organism>/, runs/
    families = "<root>/families"   # IDTRACK_REPRO_FAMILIES     the benchmark families: <family>/queries, runs, ...
    corpus = "<root>/corpus"       # IDTRACK_REPRO_CORPUS       the regression corpus, its references and runs
    hunts = "<root>/hunts"         # IDTRACK_REPRO_HUNTS        the weakness hunts' outputs
    performance = "<root>/performance"  # IDTRACK_REPRO_PERFORMANCE
    graphs = "<root>/graphs"       # IDTRACK_REPRO_GRAPHS       graph files built by the trees
    cache = "<root>/cache"         # IDTRACK_REPRO_CACHE        caches of truth readers (rebuilt when missing)
    modest = "<root>/modest"       # IDTRACK_REPRO_MODEST       the modest set's runs, one folder per tree
    iterations = "<root>/iterations"  # IDTRACK_REPRO_ITERATIONS  the iterations' plans, states and reports, per tree

Trees (an idtrack checkout at a commit, the Python that runs it, its graphs per organism):

    [trees.d4c2dda]
    commit = "d4c2dda25b61cf0b33b0a178bed22f745de676d6"
    source = "trees/d4c2dda"       # a checkout of github.com/theislab/idtrack at that commit
    python = "python"              # default: the Python running this package
    [trees.d4c2dda.graphs]
    sus_scrofa = "graphs/d4c2dda/sus_scrofa/graph_sus_scrofa_assembly11_min56_max115_narrow.npz"

External inputs a family reads (family D: the atlas project's checkout and its per-feature tables):

    [sources]
    rosetta = "<checkout of the atlas project>"
    rosetta_m02 = "<its per-feature locus tables (m02)>"

The iterations (`idtrack-repro iterate`, docs/iterate.md): the trees in the order they were made; a tree's previous
iteration is the tree before it, the reference the one every iteration is also compared with:

    [iterations]
    order = ["d4c2dda", "r1", "r1b", "r2", "main"]
    reference = "d4c2dda"          # default: the first of `order`

Jobs (`runner` is `slurm` or `local`; IDTRACK_REPRO_RUNNER):

    [jobs]
    runner = "local"
    name_prefix = "idtrack_repro_"
    [jobs.slurm]
    partition = "cpu_p"            # IDTRACK_REPRO_SLURM_PARTITION
    qos = "cpu_normal"             # IDTRACK_REPRO_SLURM_QOS
    extra = []                     # further sbatch arguments
    [jobs.memory_gb]               # memory of a job that loads a graph, per organism
    homo_sapiens = 90
    [jobs.modest_memory_gb]        # memory of a modest-set job (its largest part's measured peak), per organism
    homo_sapiens = 8
    [jobs.limits]                  # read against the queue before each submission of `iterate`
    max_jobs = 30                  # the user's jobs at once (running and pending), these included
    max_memory_gb = 300            # the memory they ask at once
    names = ""                     # count only the jobs whose names start so (default: every job of the user)

The module reads no idtrack and runs on Python 3.11 too (the Python of an older tree imports it).
"""
import os
import sys
import tomllib

ORGANISMS = ["sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis"]
LOCATIONS = ["truth", "harness", "families", "corpus", "hunts", "performance", "graphs", "cache", "modest", "iterations"]
# Memory of one job that holds a graph of the reference tree (d4c2dda) and asks it, with room for the run (GB).
MEMORY_GB = {"homo_sapiens": 90, "mus_musculus": 50, "sus_scrofa": 24, "macaca_mulatta": 24, "macaca_fascicularis": 20}
# Memory of one job of the modest set (one organism; its parts run one after another, so the largest part's peak): the
# largest measured peak of its parts plus a quarter (docs/modest.md, "Size"): the scoring of the harness's sets in one
# partition (pig 3.52, mouse 5.45, human 6.14, rhesus 2.52, cynomolgus 1.65 GB: the kernel's peak of the scorer, 2026-11-20),
# the query shards (pig 3.5, mouse 4.6, human 4.8, rhesus 3.7, cynomolgus 3.2 GB) and the modest jobs of 2026-11-20 (at
# most 5.6 GB, mouse; 5.7 GB human).
MODEST_MEMORY_GB = {"homo_sapiens": 8, "mus_musculus": 7, "sus_scrofa": 5, "macaca_mulatta": 5, "macaca_fascicularis": 4}
SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the folder holding the package (src/)


class Config:
    def __init__(self, path=None):
        self.path = path
        raw = {}
        if path:
            with open(path, "rb") as fh:
                raw = tomllib.load(fh)
        self.raw = raw
        base = os.path.dirname(os.path.abspath(path)) if path else os.getcwd()
        data = raw.get("data", {})

        def where(value, rel_to):
            return os.path.normpath(os.path.join(rel_to, os.path.expanduser(value)))

        env_root = os.environ.get("IDTRACK_REPRO_DATA_ROOT")
        self.root = where(env_root, os.getcwd()) if env_root else where(data.get("root", "data"), base)
        for loc in LOCATIONS:
            env = os.environ.get(f"IDTRACK_REPRO_{loc.upper()}")
            if env:
                val = where(env, os.getcwd())
            elif loc in data:
                val = where(data[loc], base)
            else:
                val = os.path.join(self.root, loc)
            setattr(self, loc, val)
        self.trees = {}
        for name, t in raw.get("trees", {}).items():
            self.trees[name] = dict(
                name=name, commit=t.get("commit", ""),
                source=where(t["source"], base) if t.get("source") else None,
                python=t.get("python") or sys.executable,
                graphs={o: where(g, base) for o, g in (t.get("graphs") or {}).items()})
        jobs = raw.get("jobs", {})
        self.runner = os.environ.get("IDTRACK_REPRO_RUNNER") or jobs.get("runner", "local")
        if self.runner not in ("local", "slurm"):
            raise SystemExit(f"jobs.runner must be `local` or `slurm`, not {self.runner!r}")
        self.name_prefix = jobs.get("name_prefix", "idtrack_repro_")
        sl = jobs.get("slurm", {})
        self.slurm = dict(partition=os.environ.get("IDTRACK_REPRO_SLURM_PARTITION") or sl.get("partition", "cpu_p"),
                          qos=os.environ.get("IDTRACK_REPRO_SLURM_QOS") or sl.get("qos"),
                          extra=list(sl.get("extra", [])))
        self.memory_gb = dict(MEMORY_GB, **jobs.get("memory_gb", {}))
        self.modest_memory_gb = dict(MODEST_MEMORY_GB, **jobs.get("modest_memory_gb", {}))
        self.modest_shards = dict(jobs.get("modest_shards", {}))  # part -> shards asked at once in a modest job
        self.python = jobs.get("python") or sys.executable  # the Python of this package in jobs
        lim = jobs.get("limits", {})
        self.limits = dict(max_jobs=lim.get("max_jobs"), max_memory_gb=lim.get("max_memory_gb"),
                           names=str(lim.get("names", "")))
        for k in ("max_jobs", "max_memory_gb"):
            if self.limits[k] is not None and (not isinstance(self.limits[k], (int, float)) or self.limits[k] <= 0):
                raise SystemExit(f"[jobs.limits] {k} must be a positive number, not {self.limits[k]!r}")
        it = raw.get("iterations", {})
        order = list(it.get("order", []))
        if len(set(order)) != len(order):
            raise SystemExit(f"[iterations] order names a tree twice: {order}")
        self.iteration = dict(order=order, reference=it.get("reference") or (order[0] if order else None))

    def source(self, name):
        """An external input named in `[sources]` (e.g. the atlas project's checkout and its feature tables that family D
        reads), or an error that says how to configure it."""
        s = self.raw.get("sources", {})
        if name not in s:
            raise SystemExit(f"external input {name!r} is not configured: add [sources] {name} = <path> to "
                             f"{self.path or 'idtrack-repro.toml'}")
        base = os.path.dirname(os.path.abspath(self.path)) if self.path else os.getcwd()
        return os.path.normpath(os.path.join(base, os.path.expanduser(s[name])))

    def tree(self, name):
        if name not in self.trees:
            raise SystemExit(f"tree {name!r} is not in the configuration ({self.path or 'no file'}): add [trees.{name}]")
        return self.trees[name]

    def tree_of_source(self, source):
        """The configured tree whose checkout is `source`, or None."""
        real = os.path.realpath(source)
        for t in self.trees.values():
            if t["source"] and os.path.realpath(t["source"]) == real:
                return t
        return None

    def python_for(self, source):
        """The Python command that runs the tree checked out at `source` (the configured tree's, else this one)."""
        t = self.tree_of_source(source)
        return t["python"] if t else sys.executable

    def graph(self, tree, organism):
        t = self.tree(tree)
        if organism not in t["graphs"]:
            raise SystemExit(f"no graph of {organism} for tree {tree}: add [trees.{tree}.graphs] {organism} = ...")
        return t["graphs"][organism]

    def previous(self, tree):
        """The iteration before `tree` in `[iterations] order`, or None (the first, or a tree not in the order)."""
        order = self.iteration["order"]
        return order[order.index(tree) - 1] if tree in order and order.index(tree) > 0 else None

    def describe(self):
        out = dict(config=self.path, root=self.root, **{loc: getattr(self, loc) for loc in LOCATIONS},
                   runner=self.runner, slurm=self.slurm, name_prefix=self.name_prefix, limits=self.limits,
                   iteration=self.iteration,
                   trees={n: {k: v for k, v in t.items()} for n, t in self.trees.items()})
        return out


_CONFIG = None


KEEP = object()


def in_root(path, outside=KEEP):
    """A path as a record stores it: relative to the data root when it lies inside it (so that the data root can move);
    else the absolute path, or `outside` where given (e.g. None: a checkout outside the data root is no part of the
    record, whose commit names it)."""
    if not path:
        return path
    full = os.path.abspath(path)
    for base in (os.path.abspath(get().root), os.path.realpath(get().root)):
        if full == base or full.startswith(base + os.sep):
            return os.path.relpath(full, base)
    return full if outside is KEEP else outside


def from_root(path):
    """A path a record stores (`in_root`) as an absolute path: a relative one is read from the data root."""
    if not path:
        return path
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(get().root, path))


def find_file():
    p = os.environ.get("IDTRACK_REPRO_CONFIG")
    if p:
        if not os.path.exists(p):
            raise SystemExit(f"IDTRACK_REPRO_CONFIG={p}: no such file")
        return os.path.abspath(p)
    if os.path.exists("idtrack-repro.toml"):
        return os.path.abspath("idtrack-repro.toml")
    return None


def get():
    """The configuration of this process (read once)."""
    global _CONFIG
    if _CONFIG is None:
        _CONFIG = Config(find_file())
    return _CONFIG


def reset():
    """Read the configuration again at the next `get()` (tests)."""
    global _CONFIG
    _CONFIG = None
