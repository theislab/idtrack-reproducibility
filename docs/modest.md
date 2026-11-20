# The modest set

A measurement of an idtrack tree that runs in a few hours on one job per organism: the check to make on a development tree before the exhaustive pass of the final tree ([Measurement design](methods.md#measurement-design)). Code: `idtrack_reproducibility.modest` (`idtrack-repro modest`).

## What it holds

Per organism:

| part | lines | how they are chosen |
| --- | --- | --- |
| the accuracy harness | every line of every published set | all: the judge's fixed queries ([the accuracy harness](harness.md)) |
| the weakness hunts | every reproducer of the organism | all ([the weakness hunts](hunts.md#reproducing)) |
| each family (A*, A, B, C, D, E) | up to 500 lines per stratum of the family's own sets, and the family's hunt lines | a seeded stratified sample: in every stratum the lines with the smallest hash of (seed, family, organism, qid); all of a stratum that is smaller. It does not depend on the order of reading and is the same in every build |

A family's part is a published selection of qids, `<families>/<family>/panels/<organism>/modest_v<k>.parquet` with its manifest (counts per stratum, seed, sha256), built by the panel machinery (`idtrack_reproducibility.families.panel`, kind `modest`: no failure lines, only the sample and the hunt lines). The truth of a line is always the latest line of its qid in the family's sets, so a revision of the truth published later applies.

## Size

`idtrack-repro modest size` computes the samples' size from the sets' manifests (no set is read); `modest build` adds each family's hunt lines (for A* and A the lines of the ids A*'s selection took from the weakness hunts, every line of a set or stratum named `hunt*` in the others). The selections built on 2026-11-19 (`modest_v1` of every family and organism; 500 lines per stratum; one job, 13 minutes, 1.0 GB):

| part | pig | mouse | human | rhesus | cynomolgus | five organisms | of them hunt lines |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| harness (every line of every published set) | | | | | | 2,615,211 | |
| A* | 8,500 | 37,526 | 192,006 | 62,867 | 17,505 | 318,404 | 276,007 |
| A | 15,000 | 22,159 | 30,862 | 18,283 | 17,604 | 103,908 | 33,745 |
| B | 10,557 | 13,076 | 13,370 | 10,832 | 9,242 | 57,077 | 0 |
| C | 32,052 | 8,000 | 7,625 | 27,333 | 23,505 | 98,515 | 0 |
| D | 12,859 | 16,560 | 21,528 | 4,500 | 4,500 | 59,947 | 0 |
| E | 50,375 | 52,347 | 51,095 | 50,300 | 49,507 | 253,624 | 38 |
| hunt reproducers | | | | | | 29 findings | |
| **total** | | | | | | **3,506,686 lines** | |

Each selection's manifest (`<families>/<family>/panels/<organism>/modest_v1.manifest.json`) holds its counts per stratum and reason, the seed and the sha256 of its file; the data manifests list them.

Time: the harness's answers of `d4c2dda` took 116 to 372 seconds of query time per organism (`query_seconds_total` of the stored scorings: pig 204, mouse 372, human 346, rhesus 170, cynomolgus 116); the families' selections add 891,475 lines (most asked under both views: two calls each on `d4c2dda`, one from `71690ab` on); loading a graph takes seconds to a minute and is repeated per part, and scoring the harness's sets takes minutes per organism. One job per organism is expected to finish in one to two hours (an estimate from these measurements, not yet a measured run). Its parts run one after another, so a job needs the largest part's peak: the configuration's `[jobs.memory_gb]` of the organism (or `--mem`). Measured on the reference tree: the peak of `d4c2dda`'s query shards was 3.5 GB (pig), 4.6 (mouse), 4.8 (human), 3.7 (rhesus) and 3.2 (cynomolgus) (the scheduler's sampled peaks of family A*'s and A's shards, the largest of the families' runs), the scoring of the harness's human sets 3.3 GB; trees from `71690ab` on hold the graph in less memory (load peaks 0.2 to 1.2 GB, [H](families/h.md#load)). The configuration of record asks for these peaks plus a quarter: 5, 6, 6, 5 and 4 GB.

## What it can and cannot detect

- **Can**: every change of an answer among the harness's 2.6 million lines (the judge's fixed queries over ids, names, external identifiers, history kinds and the travel matrix), each judged right or wrong as the exhaustive harness judges it; every change of a known finding of the hunts (a reproducer exits 1 and names the field); in every stratum of every family, a change that makes a share p of its lines wrong is missed with probability (1 - p)^500, under 5 % for p of 0.6 % or more (for a stratum of fewer than 500 lines: every line, so every change). A stratum's rate is estimated to within about ±4.4 percentage points at worst (the Wilson 95 % interval on 500 lines at a rate of one half; much less near 100 %).
- **Cannot**: a regression confined to fewer than about 0.6 % of a family's stratum outside the harness's sets (it may be missed), a regression in a line never sampled (most of A's 400 million and A*'s 290 million rows are not), and the families' parts run outside the common runner (E's determinism and concurrency runs, metamorphic relations needing two calls, budgets, the harmonizer, the web layer; A's graph fidelity A0; H's timings). These need the panel ([Measurement design](methods.md#measurement-design)) or the exhaustive pass of the final tree.

## Commands

```
# once: build the families' selections (reads each family's sets once; built 2026-11-19: 13 minutes in one job)
idtrack-repro modest build                      # --per-stratum 500 by default; --dry counts without writing

# measure a tree: one job per organism (the configured runner); the reference too, once, for the comparison
idtrack-repro modest run --tree d4c2dda --out <folder>/d4c2dda
idtrack-repro modest run --tree <tree> --out <folder>/<tree>

# compare the tree with the reference, part by part (exit 1 on a regression: right -> wrong or right -> error)
idtrack-repro modest compare <folder>/d4c2dda <folder>/<tree> --out <folder>/<tree>/compare.tsv
```

`<tree>` is a configured tree (`[trees.<tree>]` in `idtrack-repro.toml`: its checkout, Python and graphs per organism). A run writes `<out>/<organism>/<part>/` (answers and scored tables of each part: `harness`, the families, `e/scores_now`, `e/scores_decided`, `hunts.json`) and `<out>/<organism>/modest.json` (the parts in order with their exit codes and seconds). One organism's run, without the job runner: `idtrack-repro modest organism --tree <tree> --organism <organism> --out <folder>`.
