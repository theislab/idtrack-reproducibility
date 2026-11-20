"""A8: the graph of an older snapshot, built by a tree (docs/families/a.md, "A8, an older snapshot").

Usage: idtrack-repro family a build_a8 --organism <o> --tree <tree> [--snapshot 110] [--out <folder>] [--submit]

Builds the graph of the organism at the snapshot release into `<families>/a/a8/<tree>/<organism>/` (default). The
tables the reference graph of the organism was built from (`<graphs>/<organism>/idtrack_tables/` and the index files
`ensembl_https_index_*.json`) are copied in first, so that the builder downloads only what they lack (the snapshot's
own history table, over HTTPS); nothing is written beside the reference graph's files. The build runs in the tree's
Python (the configuration's `[trees.<tree>]`: its checkout first on the path) through the tree's public api:
`API(local_repository=...).build_graph(organism_name, snapshot_release)` for d4c2dda, `API.build(organism, snapshot,
repository)` from 71690ab on. With `--submit` it is a job of the configured runner. A folder that already holds the
graph is left as it is (the builder reads it).
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C

BUILD = {
    "d4c2dda": "import idtrack; idtrack.API(local_repository={repo!r}).build_graph(organism_name={org!r}, "
               "snapshot_release={snap})",
    "convert": "import idtrack; idtrack.API.build({org!r}, {snap}, {repo!r})",
}


def prepare(organism, out):
    """Copy the reference graph's downloaded tables and index files into the build folder (once)."""
    src = os.path.join(CFG.get().graphs, organism)
    os.makedirs(out, exist_ok=True)
    if not os.path.isdir(os.path.join(out, "idtrack_tables")) and os.path.isdir(os.path.join(src, "idtrack_tables")):
        shutil.copytree(os.path.join(src, "idtrack_tables"), os.path.join(out, "idtrack_tables"))
    for p in glob.glob(os.path.join(src, "ensembl_https_index_*.json")):
        if not os.path.exists(os.path.join(out, os.path.basename(p))):
            shutil.copy2(p, out)


def command(tree, organism, snapshot, out):
    t = CFG.get().tree(tree)
    code = BUILD[C.generation(t["source"])].format(repo=os.path.abspath(out), org=organism, snap=snapshot)
    return [t["python"], "-c", f"import sys; sys.path.insert(0, {t['source']!r}); {code}"]


def main():
    ap = argparse.ArgumentParser(prog="idtrack-repro family a build_a8", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--snapshot", type=int, default=110)
    ap.add_argument("--out", default=None, help="default: <families>/a/a8/<tree>/<organism>")
    ap.add_argument("--submit", action="store_true", help="a job of the configured runner")
    ap.add_argument("--mem", type=int, default=8, help="GB of the job (a cold human build of the reference: 4.2 GB)")
    ap.add_argument("--time", default="04:00:00")
    a = ap.parse_args()
    out = a.out or os.path.join(C.family_dir("a"), "a8", a.tree, a.organism)
    if a.submit:
        from idtrack_reproducibility import jobs
        args = f"--organism {a.organism} --tree {a.tree} --snapshot {a.snapshot} --out {os.path.abspath(out)}"
        jid = jobs.submit(f"a_build_a8_{a.tree}_{a.organism}", jobs.module(CFG.get().python,
                          "idtrack_reproducibility.families.a.build_a8", args), os.path.join(out, "logs"),
                          cpus=2, mem_gb=a.mem, time=a.time)
        print(jid)
        return 0
    prepare(a.organism, out)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0", IDTRACK_ENSEMBL_ACCESS="https")
    return subprocess.run(command(a.tree, a.organism, a.snapshot, out), cwd=out, env=env).returncode


if __name__ == "__main__":
    sys.exit(main())
