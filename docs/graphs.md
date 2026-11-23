# A tree's graphs

Every benchmark asks a tree on its graphs (`[trees.<tree>.graphs]` of the configuration, [README](../README.md#trees)). This page is how a tree's graphs are made and checked before they are asked: built from stored tables with the tree's own builder, judged against Ensembl's history and against another tree's graphs, and, for load measurements, copied with every member stored. Code: `idtrack_reproducibility.graphs`. The reference graphs every stored result of `d4c2dda` and `71690ab` was measured on are described in [data](data.md#layout); this page is about the graphs of the trees after them.

## Places

| path (under `<graphs>`) | what |
| --- | --- |
| `<tree>/<organism>/graph_<organism>_assembly<key>_min<first>_max<last>.npz` | the graph: name it in `[trees.<tree>.graphs]` |
| `<tree>/<organism>/build.json` | what built it: the tree, its commit and whether its checkout was modified, the tables folder, the build's wall time and peak memory, the file's size and sha256, the content digest |
| `<tree>/<organism>/judge.json` | its judgement (below) |
| `<tree>/<organism>/build/` | the local repository of the build: links to the stored tables, the tables the builder derived, the graph file; left out of the data manifests |
| `<tree>/<organism>/logs/` | the jobs' logs |

## Building

```
idtrack-repro graph-build --tree <tree> --organisms <o,...|all> [--snapshot 115] [--tables <folder with {organism}>] --submit
```

One job per organism through the configured runner (`graph_build_<tree>_<organism>`; memory and time per organism in `graphs.common`, from the measured peaks of builds of snapshot 115 plus a quarter). The build folder is first made a farm of the stored tables: every raw Ensembl table of `--tables` (a local repository of idtrack in which graphs were built; default the reference graphs' folder `<graphs>/<organism>/`) is linked file by file to its resolved file, Ensembl's assembly chain files too, the index files are copied, and the tables a build derives are left out, so that the tree's builder derives them itself. idtrack writes a table it lacks as a new file and never rewrites one it has, so nothing is written through a link; `build.json` records the linked files whose size or time changed during the build (none, else the command fails). A table the stored tables lack is downloaded by the builder into the build folder; move it into the tables folder to keep it (Ensembl's downloaded tables are never deleted). The build runs in the tree's Python with its checkout first on the path, through its public api (`API.build(organism, snapshot, repository)`; d4c2dda's generation `API(local_repository=...).build_graph(...)`), and the graph file is hard-linked beside the build folder. A graph already in the build folder is loaded, not built again: a job lost to the scheduler is submitted again as it was. The command prints the line to add to the configuration.

Equal graphs have equal content digests: a tree whose builder did not change builds, from the same tables, a graph with the content digest of the earlier tree's graph (`build.json`, `content_digest`). Example (2026-11-24): the package's release commit `05f0719` built the cynomolgus graph from the tables of the round-2 builds in PROOF_BUILD, with the content digest of the configured graph of tree `r2`, whose builder it keeps.

## Judging

```
idtrack-repro graph-judge --tree <tree> --organisms <o,...|all> [--against <tree>] [--repository <folder>] --submit
```

One job per organism (`graph_judge_<tree>_<organism>`), in the tree's Python, on the configured graph; `judge.json` holds three parts and the exit status is 1 on a fault:

- **verify**: the tree's own `idtrack verify <graph> --repository <build folder>`: the graph's format and record for the tree's organism table, its content digest computed again, and the package's self-tests (`GraphChecks`) against the tables it was built from. A tree without the command is marked `not available`.
- **history**: the graph's gene history against Ensembl's, read from the snapshot release's core database in the truth store (`<truth>/<organism>/core/`, `stable_id_event` and `mapping_session`; nothing of idtrack reads them). Every gene event between two different stable ids whose old release lies in the graph's window is `held` (a history edge between the two versions), `placed` (the graph's record lists it as placed with the version of the release), `ids_in_no_table` (an id no table of the graph's releases has), or `missing` (a fault); every history edge between two different stable ids that is no event of Ensembl is a fault (`edges_without_event`); every gene version of the backbone has its place at every release at which it is active (`genes_without_place_at_an_active_release`, a fault otherwise). Graphs without a record and places (before the graph format of 2026-11) are marked `not available`.
- **diff** (with `--against <tree>`): what differs from that tree's configured graph, read from both files' arrays without idtrack (`graphs.diff`): nodes added and removed by type (external names also by what they are: no letter or digit, UniProt's evidence text), cross-reference links per database (with the same releases, other releases, only in one; moved to another database with the same releases), history edges added and removed. A difference is no fault: it is what a change of the builder did, to be explained by its decision.

Example: PROOF_JUDGE

The answers of a tree on its graphs are then compared with an earlier tree's by the [regression corpus](corpus.md) (`idtrack-repro corpus`) and the [modest set](modest.md).

## A stored copy

```
idtrack-repro graph-stored <graph file> --out <folder>
```

A graph file is a zip archive of numpy arrays whose members idtrack compresses where that makes them markedly smaller. The copy holds the same members with the same bytes and dates, each stored uncompressed, under the same name in `--out` (written to a temporary name, renamed when complete; checked member by member). To compare the load of the two layouts, configure the tree a second time with the stored copies as its graphs (`[trees.<tree>_stored]`, the same `commit`, `source` and `python`) and run [family H](families/h.md)'s `load` measure on both: `idtrack-repro performance --organisms <o> --measures load --trees <tree>,<tree>_stored --reps 3 --submit` (fresh processes, interleaved on one node, three repetitions). H loads with the page cache as the node holds it; a load after the file was dropped from the page cache is not measured.
