# Data

The benchmarks read and write a **data root** outside git: truth files, query sets, the answers of each idtrack tree, scored runs, graphs. The repository holds, for every file of it a result rests on, a row of a manifest: enough to check a copy of the data and to rebuild it from its sources.

## Layout

| location (configuration) | holds | made by |
| --- | --- | --- |
| `truth/<organism>/` | Ensembl's and the authorities' files: `gtf/genes_<release>.tsv.gz` and `tx_<release>.tsv.gz` (from the GTF of every release), `core/` (history tables of the release-115 core database), `biomart/<release>/` (archive BioMart tables), `chain/` (assembly chain files), `exceptions/<release>.tsv` (assembly exceptions), `grc_alignments/`, `core_places/`, `vgnc/`, `vega/`, `rest/` (Ensembl REST records), and the authorities' files of the later families; `SOURCES.tsv` lists every fetched file with URL, size, md5 and time | `idtrack-repro fetch-gtf`, `fetch-tables`, `fetch-rest` |
| `harness/queries/<organism>/` | the harness's query sets `v<n>.jsonl.gz` with `v<n>.manifest.json` (fixed once published) | `idtrack-repro make-queries` |
| `harness/answers/<tree>/<organism>/{kinds,travel}/results/` | the answers of a tree to the harness's sets (`shard_<i>_of_<n>.jsonl`, `load_<i>_of_<n>.json`) | `idtrack-repro harness-score` |
| `harness/runs/<tree>/<organism>/{kinds,travel}/` | the scored tables of those answers under the current judge (`outcomes.tsv.gz`, `scores.tsv`, `wrong.tsv`, ...; the travel matrix's `travel_report.tsv`) | `idtrack-repro harness-score --merge`, `travel-report` |
| `families/<family>/queries/<organism>/` | the family's sets (`<set>.parquet` or `.jsonl.gz`, with `<set>.manifest.json`) | the family's generator |
| `families/<family>/runs/<tree>/<organism>/` | `answers/` (the tree's answers), `outcomes/` and the tables of the scoring, `explain_wrong.jsonl.gz`, `causes.tsv.gz` | `idtrack-repro run`, `score`, `explain`, `causes` |
| `families/<family>/...` | the family's other inputs (A*'s `selection/` and REST `truth/`, A's `patterns/`, `a0/`, `a8/`) | the family's modules |
| `graphs/<organism>/` | the graph of the reference tree (`graph_*.npz`) and the tables idtrack downloaded to build it (`idtrack_tables/`) | `idtrack build` at `d4c2dda` |
| `corpus/`, `hunts/`, `performance/` | the regression corpus and its references; the hunts' outputs; the performance measurements | [corpus](corpus.md), [hunts](hunts.md), [H](families/h.md) |
| `cache/` | caches of the truth readers (rebuilt when missing; never published) | the readers |

Every run folder records what produced it: `run.json` (the idtrack commit and whether its checkout was modified, the graph's content digest, the sets and their sha256, the judge, the time) and, per shard, `load_<i>_of_<n>.json` (the Python, the load time, files written beside the graph: none).

## Manifests

`manifests/areas.toml` declares the areas of the data root (`truth`, `harness/queries`, `harness/answers`, `harness/runs`, `families/<family>`, `graphs`): their folders, what they leave out (logs, scratch, superseded runs) and how their files are rebuilt. `idtrack-repro manifest <area>` writes `manifests/<area>.tsv` (`.tsv.gz` above 5,000 rows) with one row per file:

| column | meaning |
| --- | --- |
| `path` | the file's path in the data root |
| `bytes`, `sha256` | its size and checksum |
| `source` | the URL a fetched file came from (from `SOURCES.tsv`); empty for a file the pipeline made |
| `time` | when it was fetched, or made |
| `rebuild` | `fetch: ...` (from its URL), `generate: ...` (a command of the package, deterministic), `run: ...` (answers of a tree: rerunning gives the same answers, other seconds), or `none: ...` where it cannot be rebuilt as it is (Ensembl REST's answers on the date fetched: a later fetch answers for a later release) |

`idtrack-repro check` compares the manifests with the data root (sizes always, sha256 of files under 50 MB; `--data` every file).

The areas' sizes, summed from their manifests (`manifests/`; the family areas list the runs the results come from):

| area | files | GB |
| --- | ---: | ---: |
| truth | 2,420 | 4.28 |
| harness/queries | 138 | 0.03 |
| harness/answers | 90 | 0.64 |
| harness/runs | 105 | 0.04 |
| harness/checks | 3 | 0.00 |
| families/f0 | 453 | 0.12 |
| families/astar | 16,686 | 6.25 |
| families/a | 22,025 | 9.27 |
| families/b | 612 | 0.61 |
| families/c | 7,159 | 4.28 |
| families/d | 1,417 | 0.70 |
| families/e | 3,017 | 3.81 |
| families/design | 10 | 0.00 |
| graphs | 4,879 | 8.90 |
| corpus | 120 | 0.09 |
| hunts | 295 | 3.72 |
| performance | 117 | 0.11 |
| **total** | **59,546** | **42.9** |

## Publication

The data root is to be published with the manuscript in a data repository that gives a DOI (Zenodo or a similar one), as several records so that a reader can take what a benchmark needs:

| record | areas | about |
| --- | --- | --- |
| truth | `truth` | 4.3 GB: Ensembl's and the authorities' files as fetched (their own licences: Ensembl's data are free to use; the authorities' terms are at the URLs in `SOURCES.tsv`) |
| harness | `harness/*` | 0.7 GB: the fixed query sets, the reference tree's answers and their scoring |
| families | `families/*` | 25.0 GB: the families' sets, answers, scorings and causes (the answers are most of it; Parquet, compressed) |
| graphs | `graphs` | 8.9 GB: the reference tree's graphs and the tables they were built from (equal graphs have equal content digests, so a rebuilt graph is checked against the digest the runs record) |
| corpus, hunts, performance | `corpus`, `hunts`, `performance` | 3.9 GB: the regression corpus and its references, the hunts' outputs, the performance measurements |

Each record carries its manifests from this repository (the same rows); the repository then names each record's DOI in this section and `idtrack-repro check --data` verifies a downloaded copy. Graphs and answers can be rebuilt, at the cost of the builds (hours per organism over HTTPS) and the runs; the truth files from Ensembl's FTP and archives as long as they serve them; Ensembl REST's records cannot be asked again for the same release, which is why they are kept.
