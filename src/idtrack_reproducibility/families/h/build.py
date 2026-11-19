"""H, measure `build`: build one organism's graph through a tree's public api; time, peak memory, bytes received, and
whether the built graph is the configured one.

Usage (in the tree's Python; `run.py` starts it): python -m idtrack_reproducibility.families.h.build --tree <tree>
    --organism <organism> --out <file.json> (--cold | --warm-from <folder holding idtrack_tables/>) [--cprofile <file>]

- `--cold`: an empty local repository: every Ensembl table is downloaded (over HTTPS), parsed and stored, then the
  graph is built.
- `--warm-from <folder>`: the stored tables of an earlier build are copied into the repository first, so nothing is
  downloaded: the build's own work alone.

The repository is `<performance>/runs/<tree>/<organism>/build_<cold|warm>/` (emptied first). Recorded: the build's
wall time, the peak resident memory of the process, the bytes received (counted at Python's socket and TLS-socket
reads, i.e. the payload after TLS), the size of the repository, the graph file's size, and the content digest of the
built graph against the configured graph's (equal digests: the same graph). `--cprofile` runs the build under
cProfile (whose overhead inflates the time: such a run is marked `profiled` and is a profile, not a timing). The
repository is removed after the measurement unless `--keep`.
"""

import argparse
import cProfile
import glob
import json
import os
import shutil
import socket
import ssl
import time

from idtrack_reproducibility.corpus import _common as CC
from idtrack_reproducibility.families.h import common as H


class NetCounter:
    """The bytes received through Python-level socket reads (TLS payload for HTTPS, raw for plain TCP)."""

    received = 0
    installed = False

    @classmethod
    def install(cls) -> None:
        if cls.installed:
            return
        cls.installed = True

        def wrap_into(orig):
            def f(self, buffer, *a, **k):
                n = orig(self, buffer, *a, **k)
                cls.received += int(n or 0)
                return n
            return f

        def wrap_recv(orig):
            def f(self, *a, **k):
                d = orig(self, *a, **k)
                cls.received += len(d) if isinstance(d, (bytes, bytearray)) else 0
                return d
            return f

        ssl.SSLSocket.recv_into = wrap_into(ssl.SSLSocket.recv_into)
        ssl.SSLSocket.recv = wrap_recv(ssl.SSLSocket.recv)
        socket.socket.recv_into = wrap_into(socket.socket.recv_into)
        socket.socket.recv = wrap_recv(socket.socket.recv)


def dir_size(path: str) -> int:
    total = 0
    for base, _dirs, files in os.walk(path):
        for f in files:
            p = os.path.join(base, f)
            if not os.path.islink(p):
                total += os.path.getsize(p)
    return total


def main(argv=None):
    ap = argparse.ArgumentParser(prog="h.build", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    kind = ap.add_mutually_exclusive_group(required=True)
    kind.add_argument("--cold", action="store_true")
    kind.add_argument("--warm-from", default=None)
    ap.add_argument("--cprofile", default="")
    ap.add_argument("--keep", action="store_true", help="keep the repository of the build")
    a = ap.parse_args(argv)
    org = a.organism
    NetCounter.install()
    report = {"measure": "build", "tree": a.tree, "organism": org, "kind": "cold" if a.cold else "warm",
              "profiled": bool(a.cprofile), **H.environment()}
    tree = H.Tree(a.tree)
    report.update(generation=tree.gen, commit=tree.cfg["commit"])
    graph = CC.graph_file(a.tree, org)
    _, snapshot = CC.graph_of(graph)
    repo = os.path.join(H.run_dir(a.tree, org), f"build_{report['kind']}")
    shutil.rmtree(repo, ignore_errors=True)
    os.makedirs(repo)
    if a.warm_from:
        shutil.copytree(os.path.join(a.warm_from, "idtrack_tables"), os.path.join(repo, "idtrack_tables"))
        for f in glob.glob(os.path.join(a.warm_from, "ensembl_https_index_*.json")):
            shutil.copy(f, repo)
        report["tables_copied_bytes"] = dir_size(repo)
    import idtrack

    profiler = cProfile.Profile() if a.cprofile else None
    received0 = NetCounter.received
    t = time.perf_counter()
    if profiler:
        profiler.enable()
    if tree.gen == "d4c2dda":
        api = idtrack.API(local_repository=repo)
        api.build_graph(organism_name=org, snapshot_release=snapshot)
    else:
        api = idtrack.API.build(org, snapshot, repo)
    if profiler:
        profiler.disable()
        profiler.dump_stats(a.cprofile)
    report["build_s"] = round(time.perf_counter() - t, 1)
    report["peak_GB"] = H.peak_gb()
    report["bytes_received"] = NetCounter.received - received0
    report["repo_bytes"] = dir_size(repo)
    built = glob.glob(os.path.join(repo, f"graph_{org}_*.npz"))
    report["graph_file_bytes"] = os.path.getsize(built[0]) if len(built) == 1 else None
    report["content_digest"] = tree.manifest(api).get("content_digest")
    reference = tree.load(tree.link(org, H.run_dir(a.tree, org)))
    report["configured_content_digest"] = tree.manifest(reference).get("content_digest")
    report["same_content_as_configured"] = report["content_digest"] == report["configured_content_digest"]
    H.dump(a.out, report)
    if not a.keep:
        shutil.rmtree(repo, ignore_errors=True)
    print(json.dumps(report, default=str))


if __name__ == "__main__":
    main()
