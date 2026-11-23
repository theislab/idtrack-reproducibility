"""Judge a built graph of a tree (docs/graphs.md, "Judging").

Usage: idtrack-repro graph-judge --tree <tree> --organisms <o,...|all> [--against <tree>] [--repository <folder>]
    [--submit] [--mem <GB>] [--time <hh:mm:ss>]

The graph is the configured one (`[trees.<tree>.graphs]`). For each organism, in its own job with `--submit` (named
`graph_judge_<tree>_<organism>`), in the tree's Python, three checks, written to `<graphs>/<tree>/<organism>/judge.json`:

- `verify`: the tree's own `idtrack verify <graph> --repository <repository>`: the graph's format and record, its
  content digest computed again, and the package's self-tests against the tables it was built from (the repository:
  `--repository`, default the graph's build folder `<graphs>/<tree>/<organism>/build/`). Trees without the command
  (before it existed) are marked `not available`.
- `history`: the graph's gene history against Ensembl's (`<truth>/<organism>/core/`, the `stable_id_event` and
  `mapping_session` tables of the snapshot release's core database, read here, nothing of idtrack reads them). Every
  gene event of Ensembl between two different stable ids whose session's old release lies in the graph's window: is
  it a history edge from the old `<id>.<version>` to the new one? An event is `held`, `placed` (the record lists it as
  placed with the version of the release), `ids_in_no_table` (an id no release of the graph has) or `missing` (a
  fault). Every history edge between two different stable ids: is it an event of Ensembl (`edges_without_event`, a
  fault, otherwise)? And every gene version of the backbone has a place at every release at which it is active
  (`genes_without_place_at_an_active_release`, a fault otherwise). For trees whose graph has a record and places (the
  graph format of 2026-11 on); otherwise `not available`.
- `diff` (with `--against <tree>`): what differs from that tree's configured graph of the organism, read from both
  files' arrays without idtrack (`graphs.diff`): nodes added and removed by type, cross-reference links per database
  (only in one, other releases, moved to another database), history edges added and removed.

The exit status is 1 when `verify` fails or `history` finds a fault.
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import subprocess
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import jobs
from idtrack_reproducibility.graphs import common as G
from idtrack_reproducibility.graphs import diff as D

EVENT_COLUMNS = ["old_stable_id", "old_version", "new_stable_id", "new_version", "mapping_session_id", "type", "score"]
SESSION_COLUMNS = ["mapping_session_id", "old_db_name", "new_db_name", "old_release", "new_release", "old_assembly",
                   "new_assembly", "created"]
FAULTS = ("missing", "edges_without_event", "genes_without_place_at_an_active_release")


def core_table(organism: str, snapshot: int, name: str):
    """A table of the snapshot release's core database from the truth store (Ensembl's MySQL dump), as text columns."""
    import pandas as pd

    folder = os.path.join(CFG.get().truth, organism, "core")
    found = sorted(f for f in os.listdir(folder) if f.startswith(f"{organism}_core_{snapshot}_") and f.endswith(f".{name}.txt.gz"))
    if len(found) != 1:
        raise SystemExit(f"expected one {name} table of release {snapshot} in {folder}, found {found} "
                         f"(idtrack-repro fetch-tables core)")
    with gzip.open(os.path.join(folder, found[0]), "rt") as fh:
        text = fh.read()
    return pd.read_csv(io.StringIO(text), sep="\t", header=None, na_values=["\\N"], keep_default_na=False, dtype=str)


def gene_events(events, sessions, first: int, last: int):
    """Ensembl's gene events between two different stable ids (versions not 0) whose old release is in [first, last]."""
    import pandas as pd

    events = events.iloc[:, :len(EVENT_COLUMNS)].copy()
    events.columns = EVENT_COLUMNS
    sessions = sessions.iloc[:, :len(SESSION_COLUMNS)].copy()
    sessions.columns = SESSION_COLUMNS
    events = events.merge(sessions[["mapping_session_id", "old_release", "new_release"]], on="mapping_session_id")
    events = events[(events["type"] == "gene") & events["old_stable_id"].notna() & events["new_stable_id"].notna()]
    events = events[events["old_stable_id"] != events["new_stable_id"]]
    events = events[(events["old_version"] != "0") & (events["new_version"] != "0")]
    old_release = pd.to_numeric(events["old_release"], errors="coerce")
    return events[(old_release >= first) & (old_release <= last)]


def compare_history(events, pairs: set, placed: dict, names: set) -> dict:
    """Ensembl's gene events against the graph's history edges between different stable ids.

    `events`: rows with old_stable_id, old_version, new_stable_id, new_version, old_release, new_release (`gene_events`);
    `pairs`: the graph's history edges `(old <id>.<version>, new <id>.<version>)` between two stable ids; `placed`:
    `{event pair: pair placed with the version of the release}` (the graph's record); `names`: the graph's node names.
    """
    held = as_placed = absent = 0
    missing = []
    truth = set()
    for row in events.itertuples(index=False):
        a, b = f"{row.old_stable_id}.{row.old_version}", f"{row.new_stable_id}.{row.new_version}"
        truth.add((a, b))
        if (a, b) in pairs:
            held += 1
        elif (a, b) in placed and placed[(a, b)] in pairs:
            as_placed += 1
            truth.add(placed[(a, b)])
        elif a not in names or b not in names:
            absent += 1
        else:
            missing.append([a, b, row.old_release, row.new_release])
    invented = sorted(p for p in pairs if p not in truth)
    return {"events": int(len(events)), "edges_between_ids": len(pairs), "held": held,
            "placed_with_the_version_of_the_release": as_placed, "ids_in_no_table": absent, "missing": len(missing),
            "missing_examples": missing[:50], "edges_without_event": len(invented),
            "edges_without_event_examples": [list(p) for p in invented[:50]]}


def history_here(graph_path: str, organism: str) -> dict:
    """The history check on a graph loaded with the tree's own code (the tree first on the path)."""
    import numpy as np

    from idtrack.graph import TheGraph

    try:
        from idtrack.graph.types import _KIND_HISTORY
    except ImportError:
        return {"status": "not available (the tree's graph has no record and places)"}
    graph = TheGraph.load(graph_path)
    store = graph._frozen()
    names = store.names
    pairs = set()
    for src, dst in zip(store.e_src[store.e_kind == _KIND_HISTORY].tolist(), store.e_dst[store.e_kind == _KIND_HISTORY].tolist()):
        a, b = names[src], names[dst]
        if a.rsplit(".", 1)[0] != b.rsplit(".", 1)[0]:
            pairs.add((a, b))
    placed = {tuple(item["event"]): tuple(item["placed"])
              for item in graph.record.history.get("events_placed_with_the_version_of_the_release", [])}
    first, last = graph.first_release, graph.snapshot_release
    events = gene_events(core_table(organism, last, "stable_id_event"), core_table(organism, last, "mapping_session"),
                         first, last)
    result = {"status": "checked", **compare_history(events, pairs, placed, set(store.index))}
    unplaced, checked, examples = 0, 0, []
    for i in np.flatnonzero(store.ntype == store.backbone_code).tolist():
        name = names[i]
        if name.endswith(".Void") or name.endswith(".Retired"):
            continue
        checked += 1
        have = 0
        for place in store.places.of_node(i):
            have |= place.releases.mask
        lacking = graph.active_releases(name).mask & ~have
        if lacking:
            unplaced += 1
            if len(examples) < 50:
                examples.append([name, [r for r in range(lacking.bit_length()) if lacking >> r & 1]])
    result.update({"places_checked_genes": checked, "genes_without_place_at_an_active_release": unplaced,
                   "unplaced_examples": examples})
    return result


def verify_here(python: str, source: str, graph_path: str, repository: str | None) -> dict:
    """The tree's own `idtrack verify` (a fresh process of the tree's Python), or `not available`."""
    main = os.path.join(source, "idtrack", "__main__.py")
    if not os.path.exists(main) or '"verify"' not in open(main).read():
        return {"status": "not available (the tree has no `idtrack verify`)"}
    command = [python, "-m", "idtrack", "verify", graph_path]
    if repository:
        command += ["--repository", repository]
    env = dict(os.environ, PYTHONPATH=os.path.abspath(source) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run(command, capture_output=True, text=True, env=env, cwd=repository or None)
    return {"status": "holds" if r.returncode == 0 else "fails", "exit": r.returncode, "repository": repository,
            "output": r.stdout.splitlines()[-40:], "errors": r.stderr.splitlines()[-20:]}


def judge_here(tree: str, organism: str, against: str | None, repository: str | None) -> dict:
    t = CFG.get().tree(tree)
    graph_path = CFG.get().graph(tree, organism)
    repository = repository or G.repository_of(tree, organism)
    if not os.path.isdir(repository):
        repository = None
    report = {"tree": tree, "organism": organism, "graph": graph_path}
    report["verify"] = verify_here(t["python"], t["source"], graph_path, repository)
    sys.path.insert(0, os.path.abspath(t["source"]))
    report["history"] = history_here(graph_path, organism)
    if against:
        report["diff"] = {"against": against, **D.diff(D.read(CFG.get().graph(against, organism)), D.read(graph_path))}
    faults = [k for k in FAULTS if report["history"].get(k)]
    report["holds"] = report["verify"]["status"] != "fails" and not faults
    report["faults"] = faults + (["verify"] if report["verify"]["status"] == "fails" else [])
    G.dump(os.path.join(G.out_dir(tree, organism), "judge.json"), report)
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro graph-judge", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--organisms", required=True, help="comma-separated, or `all`")
    ap.add_argument("--against", default=None, help="a tree whose configured graph the graph is compared with")
    ap.add_argument("--repository", default=None, help="the tables the graph was built from (default: its build folder)")
    ap.add_argument("--submit", action="store_true", help="one job per organism through the configured runner")
    ap.add_argument("--mem", type=int, default=None, help="GB of each job (default: graphs.common.MEMORY_GB)")
    ap.add_argument("--time", default=None, help="wall time of each job (default: graphs.common.TIME)")
    ap.add_argument("--here", action="store_true", help=argparse.SUPPRESS)  # the job itself, in the tree's Python
    a = ap.parse_args(argv)
    organisms = CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")
    t = CFG.get().tree(a.tree)
    status = 0
    for org in organisms:
        if a.here:
            r = judge_here(a.tree, org, a.against, a.repository)
            h = r["history"]
            print(f"{org}: verify {r['verify']['status']}; history {h.get('status')}: held {h.get('held')}, placed "
                  f"{h.get('placed_with_the_version_of_the_release')}, missing {h.get('missing')}, edges without event "
                  f"{h.get('edges_without_event')}, genes without place {h.get('genes_without_place_at_an_active_release')}"
                  + (f"; diff against {a.against}: nodes {r['diff']['nodes']['old']} -> {r['diff']['nodes']['new']}, "
                     f"history edges +{r['diff']['history']['n_added']} -{r['diff']['history']['n_removed']}"
                     if a.against else "") + ("" if r["holds"] else f"; FAULTS {r['faults']}"))
            status = status or (0 if r["holds"] else 1)
            continue
        args = f"--tree {a.tree} --organisms {org} --here"
        args += f" --against {a.against}" if a.against else ""
        args += f" --repository {os.path.abspath(a.repository)}" if a.repository else ""
        command = jobs.module(t["python"], "idtrack_reproducibility.graphs.judge", args)
        if a.submit:
            print(org, jobs.submit(f"graph_judge_{a.tree}_{org}", command, os.path.join(G.out_dir(a.tree, org), "logs"),
                                   cpus=2, mem_gb=a.mem or G.MEMORY_GB.get(org, 16), time=a.time or G.TIME.get(org, "06:00:00")))
        else:
            r = subprocess.run(command, shell=True, env=dict(os.environ, **jobs.job_env()))
            status = status or r.returncode
    return status


if __name__ == "__main__":
    sys.exit(main())
