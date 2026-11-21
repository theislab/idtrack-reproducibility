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

`idtrack-repro modest size` computes the samples' size from the sets' manifests (no set is read); `modest build` adds each family's hunt lines (for A* and A the lines of the ids A*'s selection took from the weakness hunts, every line of a set or stratum named `hunt*` in the others). The selections built on 2026-11-19 (`modest_v1` of every family and organism; 500 lines per stratum; one job, 13 minutes, 1.0 GB), with B and E built again on 2026-11-21 (`modest_v2`, below):

| part | pig | mouse | human | rhesus | cynomolgus | five organisms | of them hunt lines |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| harness (every line of every published set) | | | | | | 2,615,211 | |
| A* | 8,500 | 37,526 | 192,006 | 62,867 | 17,505 | 318,404 | 276,007 |
| A | 15,000 | 22,159 | 30,862 | 18,283 | 17,604 | 103,908 | 33,745 |
| B (`modest_v2`) | 19,349 | 23,390 | 24,159 | 20,974 | 12,242 | 100,114 | 0 |
| C | 32,052 | 8,000 | 7,625 | 27,333 | 23,505 | 98,515 | 0 |
| D | 12,859 | 16,560 | 21,528 | 4,500 | 4,500 | 59,947 | 0 |
| E (`modest_v2`) | 57,864 | 61,280 | 60,164 | 57,096 | 54,740 | 291,144 | 38 |
| hunt reproducers | | | | | | 34 findings | |
| **total** | | | | | | **3,587,243 lines** | |

`modest_v2` of B and E (2026-11-21; one job, 5 minutes) adds the strata of the sets published for the attack on idtrack's next version: B9's `b9_*` (pig 8,792 lines in 19 strata, mouse 10,314 in 27, human 10,789 in 28, rhesus 10,142 in 24, cynomolgus 3,000 in 6) and E8's `e8_batch_*` and `e8_onerelease_*` (pig 6,838, mouse 8,177, human 8,379, rhesus 6,133, cynomolgus 4,543 lines in 10 or 11 strata, their reference lines among them); every qid of `modest_v1` is in it (the sample is a hash of the qid). A relation line is judged against the line it names, so a selection now holds that line too (reason `relation`: E pig 3,269, mouse 4,031, human 4,169, rhesus 2,973, cynomolgus 2,533, for E8's earlier pairs as well; in `modest_v1` a sampled pair whose partner was not sampled had no outcome). A line of `e8_batch` names a whole-release list, converted once per list and arguments: the E part converts up to every list of the organism twice (68 lists of up to 77,000 names in human and mouse).

Each selection's manifest (`<families>/<family>/panels/<organism>/modest_v1.manifest.json`, `modest_v2`) holds its counts per stratum and reason, the seed and the sha256 of its file; the data manifests list them.

Time: the harness's answers of `d4c2dda` took 116 to 372 seconds of query time per organism (`query_seconds_total` of the stored scorings: pig 204, mouse 372, human 346, rhesus 170, cynomolgus 116); the families' selections add 891,475 lines (most asked under both views: two calls each on `d4c2dda`, one from `71690ab` on); loading a graph takes seconds to a minute and is repeated per part, and scoring the harness's sets takes minutes per organism. One job per organism is expected to finish in one to two hours (an estimate from these measurements, not yet a measured run). Its parts run one after another, so a job needs the largest part's peak: the configuration's `[jobs.modest_memory_gb]` of the organism (default `config.MODEST_MEMORY_GB`; or `--mem`). The largest part is the scoring of the harness's sets, one partition of every line: its peak, read by the kernel for the scorer's process (`wait4`), was 3.52 GB (pig), 5.45 (mouse), 6.14 (human), 2.52 (rhesus) and 1.65 (cynomolgus) on 2026-11-20; the query shards of `d4c2dda` peaked at 3.5, 4.6, 4.8, 3.7 and 3.2 GB (the scheduler's sampled peaks of family A*'s and A's shards), and trees from `71690ab` on hold the graph in less memory (load peaks 0.2 to 1.2 GB, [H](families/h.md#load)); the modest jobs of 2026-11-20 peaked at 2.2 (pig), 5.6 (mouse), 5.7 (human), 2.1 (rhesus) and 2.1 GB (cynomolgus), the scheduler's samples. A job asks for the largest of these plus a quarter: 5, 7, 8, 5 and 4 GB. (The first modest run asked for the query shards' peak alone, and the harness scoring of human was killed at 6.19 GB.)

**Answers already in the data root are taken, not asked again.** A run of the reference tree asks lines the exhaustive runs of that tree have answered (`<families>/<family>/runs/d4c2dda/<organism>/` and its other runs, F0's check runs for the harness's sets). The runners take such an answer by qid when the stored run's records match: the same commit, the same graph content digest (the stored `run.json` against the shard's load record) and the line's set read there with the same sha256 (`runner.reuse`); only the rest is asked, and the load records and `modest.json` (`reused`: lines per source run) say what was taken from where. The first modest run of `d4c2dda` (2026-11-20, before this) asked everything: its human family D took about 3.7 seconds per source-inference line (`d2_ids`, one `infer_source` call over a list of ids: `d4c2dda`'s source inference is 12 to 15 times slower than `71690ab`'s, 0.24 s, the cause of the slow part; the two views of the convert lines cost nothing) and D4's batches 0.4 to 0.7 s per line (0.01 s on `71690ab`), so pig's and mouse's D took about three hours and human's did not finish in the job's eight hours; the exhaustive runs hold all those answers. A part that is still long is asked in shards at once by the runner's own sharding (`--shards d=8` or `[jobs.modest_shards]`; each shard a process of the job, which then asks the shards' CPUs and, where more than the part's peak, the shards' graphs: `n` times the organism's `[jobs.memory_gb]`).

`modest organism ... --score-only` scores a run's stored answers again (the judge of today: revision sets published since the run, the judge's code), asking nothing and running no reproducer; it records its steps in `modest_scored_again.json`, beside the run's own `modest.json`. A run of an organism again (to complete a part) keeps the earlier run's record as `modest.earlier_<k>.json`.

`modest organism ... --again <part>=<revision set>` asks the lines of a revision set again on one tree: a revision that changes a line's question for a tree, not only its truth (family C's `rfam_label_rows_v1`: a label line is asked under the graph's name of the labels where the graph records one, `RFAM label` from idtrack's next version on), makes that tree's stored answers answers to another question. Only the named parts run: the runner asks the set's lines of the selection whatever the folder holds, into `answers/again_<set>_<i>_of_<n>.<k>.parquet` (taking nothing from other runs), and every reader takes a folder's answer files in `common.answer_files` order (the runner's parts, then the `again_` parts), so the answer to the new question wins; the steps are recorded in `modest_again_<k>.json` beside the run's own record. On 2026-11-20: C's label lines on `next` (pig, rhesus, cynomolgus; `d4c2dda` and `71690ab` hold labels and accessions under one name and ask as before), and D2's lines on `next` (`d2_r4`: asked again so that the answer keeps the order `infer_source` returned, `flags.order`, which the band's price reads; the question is the same).

A part killed by a signal (the scheduler's SIGKILL when the job exceeds its memory) stops the organism's run: `modest.json` names it (`killed`: the part, the command, the signal), the parts after it are not run, and the job ends with an error; a part that fails otherwise skips the rest of its own commands and the run goes on (`modest.json` records every exit).

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

# a revision that changed the question of its lines for one tree: those lines asked again on that tree, then scored
idtrack-repro modest organism --tree <tree> --organism <o> --out <folder>/<tree> --again c=rfam_label_rows_v1

# compare the tree with the reference, part by part (exit 1 on a regression: right -> wrong or right -> error)
idtrack-repro modest compare <folder>/d4c2dda <folder>/<tree> --out <folder>/<tree>/compare.tsv
```

`<tree>` is a configured tree (`[trees.<tree>]` in `idtrack-repro.toml`: its checkout, Python and graphs per organism). A run writes `<out>/<organism>/<part>/` (answers and scored tables of each part: `harness`, the families, `e/scores_now`, `e/scores_decided`, `hunts.json`) and `<out>/<organism>/modest.json` (the parts in order with their exit codes and seconds, the answers each runner took from other runs of the tree, a part killed by a signal). One organism's run, without the job runner: `idtrack-repro modest organism --tree <tree> --organism <organism> --out <folder>`.
