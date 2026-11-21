# idtrack-reproducibility

Reproducibility material for [idtrack](https://github.com/theislab/idtrack): the benchmarks that measure its accuracy, robustness and speed against independent truth, and the results a manuscript cites. Every expected answer is read from Ensembl's own files or from the authorities (HGNC, MGI, VGNC, NCBI, UniProt, ...), never from idtrack. The repository will grow with other material for the manuscript.

## What is here

| benchmark | what it measures | methods | code |
| --- | --- | --- | --- |
| the accuracy harness | the judge: about 2.6 million fixed queries over five organisms (ids, names, external identifiers, history kinds, the travel matrix between every pair of releases) | [docs/harness.md](docs/harness.md) | `idtrack_reproducibility.harness` |
| F0, the common machinery | query lines, the adapter for each idtrack api, the sharded runner, the scorer, reports with clustered intervals, comparison of trees, causes, limitations | [docs/families/f0.md](docs/families/f0.md) | `idtrack_reproducibility.families` |
| A*, the hardest lineages | history travel for the ids whose lineages are hardest, every origin to every target | [docs/families/astar.md](docs/families/astar.md) | `families.astar` |
| A, history travel | history travel exhaustive over every id and release, and the graph's fidelity to Ensembl's tables (A0) | [docs/families/a.md](docs/families/a.md) | `families.a` |
| B, names over time | display names, symbols of the authorities, aliases and their forms at every release (built; scored on rhesus and cynomolgus) | [docs/families/b.md](docs/families/b.md) | `families.b` |
| C, external identifiers | NCBI gene, RefSeq, UniProt and other databases' identifiers to genes and back, at every release, against Ensembl's and the databases' own tables (built; scored on pig, rhesus, cynomolgus) | [docs/families/c.md](docs/families/c.md) | `families.c` |
| D, real-world inputs | an atlas project's feature lists, source inference, marker lists, mixed lists, judged by locus | [docs/families/d.md](docs/families/d.md) | `families.d` |
| E, robustness | missing values, malformed strings, types, other organisms' ids, arguments, budgets, determinism and concurrency, metamorphic relations, the web layer | [docs/families/e.md](docs/families/e.md) | `families.e` |
| H, performance | build, load, memory, per-call latency, batch scaling, harmonization | [docs/families/h.md](docs/families/h.md) | `families.h` |
| the regression corpus | whether an answer changed between two idtrack trees (no truth) | [docs/corpus.md](docs/corpus.md) | `idtrack_reproducibility.corpus` |
| the weakness hunts | findings of the hunts (inputs, history, regions; the attack on idtrack's next version), with runnable reproducers | [docs/hunts.md](docs/hunts.md) | `idtrack_reproducibility.hunts` |

Every benchmark keeps one row per query line, view and subject (a tree or another tool): [the record](docs/records.md) (its columns, how to load it, a worked example). How the move into this package was verified (the stored scorings scored again, line by line) is [docs/verification.md](docs/verification.md). A measurement of any tree in a few hours is the [modest set](docs/modest.md).

The design of the whole programme (principles, needs, coverage, the families F and G still to come, intervals, known limitations, the comparison with other tools) is [docs/methods.md](docs/methods.md). The rules of idtrack and of the judge that the benchmarks cite by name are in [docs/rules.md](docs/rules.md).

## Layout

| path | what |
| --- | --- |
| `src/idtrack_reproducibility/` | the package: `config.py` (locations), `jobs.py` (runners), `cli.py` (the commands), `harness/`, `families/` (the common machinery and one package per family), `corpus/`, `hunts/`, `results.py`, `records.py` (the record's loader), `rescore.py` (scoring again), `manifests.py`, `check.py` |
| `docs/` | the methods of each benchmark, the rules, the data, how to add a benchmark |
| `results/` | the current results: tables per benchmark, tree and organism, `index.tsv` (one row per family, stratum, organism and tree) and `provenance.tsv` (what generated each file, from which source) |
| `results.toml` | the registry of the scored runs whose tables are the results |
| `manifests/` | one manifest per area of the data root (path, size, sha256, source, time, how to rebuild) and `areas.toml` |
| `hunts/` | the weakness hunts' findings tables |
| `tests/` | unit tests on small invented cases (no data needed) |
| `environments/` | the environment of the package (`repro.yml`) and of the trees the results were made with (`tree_<tree>.yml`) |
| `idtrack-repro.example.toml` | the configuration, commented |

The data (truth files, query sets, answers, scored runs, graphs) are not in git: they live in a data root whose layout and manifests are described in [docs/data.md](docs/data.md).

## Installing

```
mamba env create -p <prefix> -f environments/repro.yml     # Python 3.13 with numpy, pandas, pyarrow, pytest
<prefix>/bin/python -m pip install --no-deps -e .
```

or `python -m pip install -e ".[test]"` in any Python 3.12 or newer. The command `idtrack-repro` (or `python -m idtrack_reproducibility.cli`) runs every step.

## Trees

The idtrack that is measured is a **tree**: a checkout of github.com/theislab/idtrack at a commit, with the Python that runs it. The results here measure two trees: `d4c2dda` (the reference every later version is compared with) and `71690ab` (the same answers with a typed api; idtrack's main branch). To measure a tree, clone idtrack, check out the commit, make the tree's environment, and name the checkout, its Python and its graphs in the configuration (`[trees.<name>]`); the adapter puts the checkout first on the path, so idtrack need not be installed in it. The environments the stored results were made in are `environments/tree_d4c2dda.yml` (Python 3.11.12: every result of `d4c2dda`) and `environments/tree_71690ab.yml` (Python 3.13.18: `71690ab` and the trees after it); a tree's runs record its Python (`load_<i>_of_<n>.json`). `pip install -e ".[idtrack]"` installs idtrack `71690ab` itself into this package's environment.

## Configuration

One TOML file and environment variables; nothing is a path of a particular machine. The file is `$IDTRACK_REPRO_CONFIG`, else `./idtrack-repro.toml`; [idtrack-repro.example.toml](idtrack-repro.example.toml) shows every setting with its default:

- `[data]`: the data root (`root`, default `./data`) and, if they lie elsewhere, each location (`truth`, `harness`, `families`, `corpus`, `hunts`, `performance`, `graphs`, `cache`); variables `IDTRACK_REPRO_DATA_ROOT`, `IDTRACK_REPRO_<LOCATION>`;
- `[trees.<name>]`: `commit`, `source` (the checkout), `python`, `[trees.<name>.graphs]` per organism;
- `[jobs]`: `runner` = `local` (every job runs at once, here) or `slurm` (`sbatch`, with `[jobs.slurm]` `partition`, `qos`, `extra`; variables `IDTRACK_REPRO_RUNNER`, `IDTRACK_REPRO_SLURM_PARTITION`, `IDTRACK_REPRO_SLURM_QOS`), `name_prefix` of the jobs' names, `[jobs.memory_gb]` per organism (a job that holds a graph and asks it) and `[jobs.modest_memory_gb]` (a job of the modest set: its largest part's measured peak, [docs/modest.md](docs/modest.md)).

`idtrack-repro config` prints the configuration as a process reads it.

## Commands

`idtrack-repro --help` lists them; each takes `--help`. The steps, in order:

| step | command |
| --- | --- |
| fetch the truth | `idtrack-repro fetch-gtf`, `idtrack-repro fetch-tables <what>`, `idtrack-repro fetch-rest`; check the history reading: `idtrack-repro check-truth` |
| make the query sets | the harness: `idtrack-repro make-queries --organism <o> --set <name>`; a family: `idtrack-repro family <family> <generator> --organism <o>` (e.g. `family a gen_a`) |
| run a tree on its graphs | `idtrack-repro run --family <family> --organism <o> --clone <checkout> --graph <graph> --out <run folder> --submit <shards> --score` (the harness: `idtrack-repro harness-score ... --submit <shards>`) |
| score | `idtrack-repro score --family <family> --organism <o> --run <run folder>` (in partitions: `--parts P --part p`, then `--parts P --merge`) |
| score again | `idtrack-repro rescore [--parts ...] --workers N`: every registered scoring under today's judge (the revision sets published since, the judge's code), from its stored answers; resumable ([the record](docs/records.md#scoring-again)) |
| compare two trees | `idtrack-repro compare <scored a> <scored b> --moved <file>` |
| causes of the wrong answers | `idtrack-repro explain ...` on the tree, then `idtrack-repro causes --organism <o> --scored <folder> [--plugin <family module>]` |
| known limitations | `idtrack-repro known-limitations` (every family's rows from its current runs, `results/limitations/`); one family's rows: `idtrack-repro limitations --family <family> --scored <folders> --out <tsv>` |
| results | `idtrack-repro results` (reads `results.toml`; writes `results/`) |
| data manifests | `idtrack-repro manifest <area>` or `--all` |
| the modest set | `idtrack-repro modest build` once; `idtrack-repro modest run --tree <tree> --out <folder>` (one job per organism); `idtrack-repro modest compare <reference> <tree>` ([docs/modest.md](docs/modest.md)) |
| checks | `idtrack-repro check` |

Each benchmark's document says how to reproduce it from an empty data root ("Reproducing from scratch"): [harness](docs/harness.md), [F0](docs/families/f0.md#reproducing-from-scratch), [A*](docs/families/astar.md#reproducing-from-scratch), [A](docs/families/a.md#reproducing-from-scratch), [B](docs/families/b.md#reproducing-from-scratch), [C](docs/families/c.md#reproducing-from-scratch), [D](docs/families/d.md#reproducing-from-scratch), [E](docs/families/e.md#reproducing-from-scratch), [H](docs/families/h.md#reproducing-from-scratch), [the corpus](docs/corpus.md), [the hunts](docs/hunts.md#reproducing).

## Results

`results/index.tsv` lists every current result: family, stratum, organism, tree, scope (`exhaustive`; `panel`: a development stage measured by the family's panel, [docs/methods.md](docs/methods.md#measurement-design); `sets`: only the sets the registry entry names were asked on that tree, every line of them; `modest`: the [modest set](docs/modest.md)'s selection of the family, up to 500 sampled lines per stratum, the rate an estimate within about ±4.4 points), idtrack commit, the graph's content digest, n, the outcomes, the rate (right / (right + wrong + error); `undetermined` and `timeout` apart, never counted) and its 95 % interval (the cluster bootstrap over the units that generated the lines where the table has it, else Wilson), and the file. `results/provenance.tsv` names, for every file, the scored folder of the data root it was generated from with its sha256 and scoring time. Regenerating replaces: only the latest results of each tree are kept; git keeps the past.

The known limitations, measured by the families, are one table: `results/limitations/limitations.tsv`.

## Checks

`idtrack-repro check` fails on a machine path or a mention of a private working record in any tracked file, a broken link or section, a family without its methods document, tests or limitations rows, a result not generated by the pipeline or missing from the index, an index row of a missing file, a result whose source has changed since (stale), a result without the rows of [the record](docs/records.md) behind it in the data root, and a data manifest that does not match its files. Without the data root it says which checks it skipped. The CI workflow (`.github/workflows/tests.yml`) runs the tests and the checks on Python 3.12 to 3.14; while the organisation's Actions are not available for private repositories, run `python -m pytest -q` and `idtrack-repro check` before every push.

## Adding a benchmark

[docs/adding.md](docs/adding.md): the layout, the names, the rules of writing, the steps and the checks a new family must pass.
