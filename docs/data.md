# Data

The benchmarks read and write a **data root** outside git: truth files, query sets, the answers of each idtrack tree, scored runs, graphs. The repository holds, for every file of it a result rests on, a row of a manifest: enough to check a copy of the data and to rebuild it from its sources.

## Layout

| location (configuration) | holds | made by |
| --- | --- | --- |
| `truth/<organism>/` | Ensembl's and the authorities' files: `gtf/genes_<release>.tsv.gz` and `tx_<release>.tsv.gz` (from the GTF of every release), `core/` (history tables of the release-115 core database), `biomart/<release>/` (archive BioMart tables), `chain/` (assembly chain files), `exceptions/<release>.tsv` (assembly exceptions), `sequences/<release>.tsv` (the sequence regions of each assembly, at the last release of the assembly, with their coordinate system and Ensembl's `karyotype_rank` and `codon_table` attributes: which are whole chromosomes and which the mitochondrial genome; `fetch_tables sequences --releases eras`), `grc_alignments/`, `core_places/`, `vgnc/`, `vega/`, `rest/` (Ensembl REST records), the extracts of the core database of every release (`core_xref/<release>.tsv.gz`: cross-references of the kept databases, family C's `fetch_c ensembl`; `core_synonyms/<release>.tsv.gz`: synonyms of every cross-reference; `core_display/<release>.tsv.gz`: each gene's display cross-reference, the name and the database it is filed under; family B's `fetch_b synonyms`, `fetch_b display`), and the authorities' files of the later families; `SOURCES.tsv` lists every fetched file with URL, size, md5 and time | `idtrack-repro fetch-gtf`, `fetch-tables`, `fetch-rest` |
| `harness/queries/<organism>/` | the harness's query sets `v<n>.jsonl.gz` with `v<n>.manifest.json` (fixed once published) | `idtrack-repro make-queries` |
| `harness/answers/<tree>/<organism>/{kinds,travel}/results/` | the answers of a tree to the harness's sets (`shard_<i>_of_<n>.jsonl`, `load_<i>_of_<n>.json`) | `idtrack-repro harness-score` |
| `harness/runs/<tree>/<organism>/{kinds,travel}/` | the scoring of those answers under the current judge: the record (`outcomes/`) and the tables derived from it (`outcomes.tsv.gz`, `scores.tsv`, `wrong.tsv`, ...; the travel matrix's `travel_report.tsv`) | `idtrack-repro harness-score --merge`, `travel-report` |
| `families/<family>/queries/<organism>/` | the family's sets (`<set>.parquet` or `.jsonl.gz`, with `<set>.manifest.json`) | the family's generator |
| `families/<family>/runs/<tree>/<organism>/` | `answers/` (the tree's answers), `outcomes/` and the tables of the scoring, `explain_wrong.jsonl.gz`, `causes.tsv.gz` | `idtrack-repro run`, `score`, `explain`, `causes` |
| `families/<family>/...` | the family's other inputs (A*'s `selection/` and REST `truth/`, A's `patterns/`, `a0/`, `a8/`; D's copies of the atlas project's files, `d/inputs/atlas/`: external inputs, not published, [D](families/d.md#truth)) | the family's modules |
| `graphs/<organism>/` | the reference graphs (`graph_*.npz`) and the tables idtrack downloaded to build them (`idtrack_tables/`) | `idtrack build` at `f27e810` (2026-11-12), see below |
| `corpus/`, `hunts/`, `performance/` | the regression corpus and its references; the hunts' outputs; the performance measurements | [corpus](corpus.md), [hunts](hunts.md), [H](families/h.md) |
| `cache/` | caches of the truth readers (rebuilt when missing; never published) | the readers |

The data root of record, from which `results/` and the manifests were made, is the folder `idtrack_reproducibility_data/` of the project's storage (`idtrack_storage/`, beside the project's folder on the authors' cluster): the areas above, the trees' checkouts `trees/<tree>/` and the cluster's configuration `idtrack-repro.toml` (its paths relative to that folder). On that cluster the trees' Pythons are named environments built from `environments/` (`idtrack-tree-71690ab` from `tree_71690ab.yml`; `d4c2dda` runs in `idtrack_dev_env`, which holds the eight packages of `tree_d4c2dda.yml` at its versions). Run records made before 2026-11-24 name, as provenance, places in the project's former work folder (`idtrack_upgrade_work/<name>/...`, e.g. `integration/builds6`, `accuracy/truth`, a working clone `clones/<name>`): that work is kept whole in the storage's `work_2026_11/<name>/...`, a removed clone as its bundle `work_2026_11/clones/<name>.bundle`; no code reads these places.

Every run folder records what produced it: `run.json` (the record's version and the scoring key, the idtrack commit and whether its checkout was modified, the graph's content digest, the sets and their sha256, the judge, the time; places in the data root relative to it, none outside it) and, per shard, `load_<i>_of_<n>.json` (the Python, the load time, files written beside the graph: none).

**The reference graphs.** The five graphs every stored result of both trees was measured on were built on 2026-11-12 by idtrack `f27e810`, an intermediate commit that is not on GitHub; `d4c2dda` (2026-11-13) keeps its builder unchanged. Their content digests are in `results.toml` (`[graphs.d4c2dda]`) and in every run's records. A build by `d4c2dda` itself with the same content digest is shown for pig and cynomolgus (warm builds, [H](families/h.md#builds)); for mouse, human and rhesus it has not been made. The cold builds of `f27e810` took 15:27 (pig, 1.09 GiB peak), 6:50 (cynomolgus, 0.50 GiB), 15:09 (rhesus, 0.98 GiB), 37:02 (mouse, 2.69 GiB) and 1:56:23 (human, 4.20 GiB), one run each (the scheduler's elapsed time and sampled peak memory). Family A's A8 graphs of snapshot 110 (pig, human) were built by `d4c2dda`: `idtrack-repro family a build_a8 --organism <o> --tree d4c2dda` ([A](families/a.md#a8-an-older-snapshot)).

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
| truth | 3,213 | 5.34 |
| harness/queries | 162 | 0.03 |
| harness/answers | 90 | 0.64 |
| harness/runs | 115 | 0.07 |
| harness/checks | 8 | 0.00 |
| families/design | 10 | 0.00 |
| families/f0 | 460 | 0.16 |
| families/astar | 16,714 | 8.07 |
| families/a | 22,071 | 12.13 |
| families/e | 3,170 | 4.33 |
| graphs | 4,899 | 10.22 |
| hunts | 299 | 3.72 |
| corpus | 170 | 0.12 |
| performance | 117 | 0.11 |
| families/b | 2,627 | 1.60 |
| families/c | 7,284 | 5.36 |
| families/d | 1,465 | 1.01 |
| families/d_atlas_inputs | 151 | 0.04 |
| **total** | **62,818** | **52.6** |

## Publication

The data root is to be published with the manuscript in a data repository that gives a DOI (Zenodo or a similar one), as several records so that a reader can take what a benchmark needs:

| record | areas | about |
| --- | --- | --- |
| truth | `truth` | 5.3 GB: Ensembl's and the authorities' files as fetched (their own licences: Ensembl's data are free to use; the authorities' terms are at the URLs in `SOURCES.tsv`) |
| harness | `harness/*` | 0.7 GB: the fixed query sets, the reference tree's answers and their scoring |
| families | `families/*` | 32.7 GB: the families' sets, answers, scorings and causes (the answers are most of it; Parquet, compressed) |
| graphs | `graphs` | 10.2 GB: the trees' graphs (the reference's and the development trees') and the tables they were built from (equal graphs have equal content digests, so a rebuilt graph is checked against the digest the runs record) |
| corpus, hunts, performance | `corpus`, `hunts`, `performance` | 4.0 GB: the regression corpus and its references, the hunts' outputs, the performance measurements |

An area marked `publish = false` in `manifests/areas.toml` is no part of any record: family D's copies of the atlas project's files (`families/d_atlas_inputs`, 150 files, 42.5 MB), whose release is the atlas project's decision; their manifest stays in this repository (each file, its origin in the project with the commit, its size and sha256), so that a reader who obtains them from the project checks them, and nothing of D's published sets, answers or scorings needs them.

Each record carries its manifests from this repository (the same rows); the repository then names each record's DOI in this section and `idtrack-repro check --data` verifies a downloaded copy. Graphs and answers can be rebuilt, at the cost of the builds (hours per organism over HTTPS) and the runs; the truth files from Ensembl's FTP and archives as long as they serve them; Ensembl REST's records cannot be asked again for the same release, which is why they are kept.
