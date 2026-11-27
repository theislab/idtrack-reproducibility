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

## Where the runs lie

A tree's run is a folder `<out>` (the `--out` of `modest run`, [Commands](#commands)) that holds `<out>/<organism>/<part>/` and the comparisons. In the data root of record ([data](data.md#layout)) the runs are `modest/<tree>/`, one per tree of the configuration (`d4c2dda`, `71690ab`, `next`, `r1`, `r1b`, `r2`); this document writes that folder `<modest>` (so `<modest>/r2/compare.tsv` is the data root's `modest/r2/compare.tsv`), and the commands below write `<folder>` for it, so that `<folder>/<tree>` is `modest/<tree>/` there. `results.toml` names the runs whose tables are the results (`folder = "modest/<tree>/{organism}/<part>"`).

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
| hunt reproducers | | | | | | every finding of the hunts' table (29 in the runs of `d4c2dda`, `71690ab` and `next`, 34 in those of `r1`, `r1b` and `r2`; 51 from 2026-11-27) | |
| **total** | | | | | | **3,587,243 lines** | |

`modest_v2` of B and E (2026-11-21; one job, 5 minutes) adds the strata of the sets published for the attack on idtrack's next version: B9's `b9_*` (pig 8,792 lines in 19 strata, mouse 10,314 in 27, human 10,789 in 28, rhesus 10,142 in 24, cynomolgus 3,000 in 6) and E8's `e8_batch_*` and `e8_onerelease_*` (pig 6,838, mouse 8,177, human 8,379, rhesus 6,133, cynomolgus 4,543 lines in 10 or 11 strata, their reference lines among them); every qid of `modest_v1` is in it (the sample is a hash of the qid). A relation line is judged against the line it names, so a selection now holds that line too (reason `relation`: E pig 3,269, mouse 4,031, human 4,169, rhesus 2,973, cynomolgus 2,533, for E8's earlier pairs as well; in `modest_v1` a sampled pair whose partner was not sampled had no outcome). A line of `e8_batch` names a whole-release list, converted once per list and arguments: the E part converts up to every list of the organism twice (68 lists of up to 77,000 names in human and mouse).

Each selection's manifest (`<families>/<family>/panels/<organism>/modest_v1.manifest.json`, `modest_v2`) holds its counts per stratum and reason, the seed and the sha256 of its file; the data manifests list them.

Time: the harness's answers of `d4c2dda` took 116 to 372 seconds of query time per organism (`query_seconds_total` of the stored scorings: pig 204, mouse 372, human 346, rhesus 170, cynomolgus 116); the families' selections add 891,475 lines (most asked under both views: two calls each on `d4c2dda`, one from `71690ab` on); loading a graph takes seconds to a minute and is repeated per part, and scoring the harness's sets takes minutes per organism. One job per organism is expected to finish in one to two hours (an estimate from these measurements, not yet a measured run). Its parts run one after another, so a job needs the largest part's peak: the configuration's `[jobs.modest_memory_gb]` of the organism (default `config.MODEST_MEMORY_GB`; or `--mem`). The largest part is the scoring of the harness's sets, one partition of every line: its peak, read by the kernel for the scorer's process (`wait4`), was 3.52 GB (pig), 5.45 (mouse), 6.14 (human), 2.52 (rhesus) and 1.65 (cynomolgus) on 2026-11-20; the query shards of `d4c2dda` peaked at 3.5, 4.6, 4.8, 3.7 and 3.2 GB (the scheduler's sampled peaks of family A*'s and A's shards), and trees from `71690ab` on hold the graph in less memory (load peaks 0.2 to 1.2 GB, [H](families/h.md#load)); the modest jobs of 2026-11-20 peaked at 2.2 (pig), 5.6 (mouse), 5.7 (human), 2.1 (rhesus) and 2.1 GB (cynomolgus), the scheduler's samples. A job asks for the largest of these plus a quarter: 5, 7, 8, 5 and 4 GB. (The first modest run asked for the query shards' peak alone, and the harness scoring of human was killed at 6.19 GB.)

**Answers already in the data root are taken, not asked again.** A run of the reference tree asks lines the exhaustive runs of that tree have answered (`<families>/<family>/runs/d4c2dda/<organism>/` and its other runs, F0's check runs for the harness's sets). The runners take such an answer by qid when the stored run's records match: the same commit, the same graph content digest (the stored `run.json` against the shard's load record) and the line's set read there with the same sha256 (`runner.reuse`); only the rest is asked, and the load records and `modest.json` (`reused`: lines per source run) say what was taken from where. The first modest run of `d4c2dda` (2026-11-20, before this) asked everything: its human family D took about 3.7 seconds per source-inference line (`d2_ids`, one `infer_source` call over a list of ids: `d4c2dda`'s source inference is 12 to 15 times slower than `71690ab`'s, 0.24 s, the cause of the slow part; the two views of the convert lines cost nothing) and D4's batches 0.4 to 0.7 s per line (0.01 s on `71690ab`), so pig's and mouse's D took about three hours and human's did not finish in the job's eight hours; the exhaustive runs hold all those answers. A part that is still long is asked in shards at once by the runner's own sharding (`--shards d=8` or `[jobs.modest_shards]`; each shard a process of the job, which then asks the shards' CPUs and, where more than the part's peak, the shards' graphs: `n` times the organism's `[jobs.memory_gb]`).

`modest organism ... --score-only` scores a run's stored answers again (the judge of today: revision sets published since the run, the judge's code), asking nothing and running no reproducer; it records its steps in `modest_scored_again.json`, beside the run's own `modest.json`. A run of an organism again (to complete a part) keeps the earlier run's record as `modest.earlier_<k>.json`.

`modest organism ... --again <part>=<revision set>` asks the lines of a revision set again on one tree: a revision that changes a line's question for a tree, not only its truth (family C's `rfam_label_rows_v1`: a label line is asked under the graph's name of the labels where the graph records one, `RFAM label` from idtrack's next version on), makes that tree's stored answers answers to another question. Only the named parts run: the runner asks the set's lines of the selection whatever the folder holds, into `answers/again_<set>_<i>_of_<n>.<k>.parquet` (taking nothing from other runs), and every reader takes a folder's answer files in `common.answer_files` order (the runner's parts, then the `again_` parts), so the answer to the new question wins; the steps are recorded in `modest_again_<k>.json` beside the run's own record. On 2026-11-20: C's label lines on `next` (pig, rhesus, cynomolgus; `d4c2dda` and `71690ab` hold labels and accessions under one name and ask as before), and D2's lines on `next` (`d2_r4`: asked again so that the answer keeps the order `infer_source` returned, `flags.order`, which the band's price reads; the question is the same).

A part killed by a signal (the scheduler's SIGKILL when the job exceeds its memory) stops the organism's run: `modest.json` names it (`killed`: the part, the command, the signal), the parts after it are not run, and the job ends with an error; a part that fails otherwise skips the rest of its own commands and the run goes on (`modest.json` records every exit).

## The tree `r1`

`r1` is idtrack's next version after the fix of its first round (commit `ef81eba`, the data root's `[trees.r1]`, its own graphs `graphs/r1/<organism>/`). Its modest run (`modest_v2` of B and E, `modest_v1` of the other families, the hunt reproducers; 2026-11-21, jobs of 2, 5 to 8 GB) took 49 minutes (pig), 1 h 56 (mouse), 2 h 30 (human), 45 minutes (rhesus) and 20 minutes (cynomolgus), the kernel's peaks 2.3, 3.5, 4.5, 2.4 and 1.1 GB; nothing of other runs was reused (no stored answers of `r1`). The tree before the fix (`next`) and the reference (`d4c2dda`) were not run: their stored answers were scored again under the same judge (`modest organism --score-only` on links to their answers), `next`'s B9 and E8 lines from its exhaustive runs of those sets (`<families>/{b,e}/runs/next/`); `d4c2dda` never asked B9 or E8 (their lines have no answer of it and are left out of its comparison), nor E8's earlier pairs' partner lines that `modest_v2` added (for `next` too).

Reading the losses found the judge lagging rules built in this round, and the revision sets published for them are listed in the families' documents: [the harness](harness.md#the-judge-following-idtracks-decided-rules) (`no_mito_nuclear_step`, `mirbase_genes`), A* and A (`no_mito_nuclear_step`), B (`<set>_rules2`; the corrections `places_not_comparable`, `support_from_xrefs`, `held_in_core`), C (`display_names_at_release`, `no_mito_nuclear_step`, `mirbase_genes`), E (`e3_dec`; `mirbase_genes` on the decided judge). Every part of every tree was then scored again, and `modest compare` run against `d4c2dda` (`<modest>/r1/compare.tsv`) and against `next` (`compare_vs_next.tsv`). The answers that went from right (or an expected error) to wrong, to an error and to undetermined, and the gains:

| organism | part | vs `next`: wrong · error · undetermined | gains | vs `d4c2dda`: wrong · error · undetermined | gains |
| --- | --- | ---: | ---: | ---: | ---: |
| pig | harness | 0 · 0 · 67 | 156 | 0 · 0 · 44 | 6,230 |
| pig | astar | 0 · 0 · 0 | 10 | 0 · 0 · 0 | 180 |
| pig | a | 0 · 0 · 12 | 36 | 0 · 0 · 6 | 918 |
| pig | b | 0 · 0 · 16 | 2,713 | 2 · 0 · 16 | 2,146 |
| pig | c | 2 · 0 · 190 | 963 | 8 · 0 · 188 | 2,295 |
| pig | d | 0 · 0 · 4 | 42 | 25 · 0 · 2 | 1,666 |
| pig | e now | 7 · 0 · 0 | 383 | 7 · 806 · 2 | 2 |
| pig | e decided | 0 · 0 · 0 | 390 | 0 · 0 · 2 | 12,470 |
| mouse | harness | 9 · 0 · 935 | 727 | 6 · 0 · 851 | 1,685 |
| mouse | astar | 6 · 0 · 14 | 26 | 6 · 0 · 14 | 1,090 |
| mouse | a | 4 · 0 · 14 | 6 | 0 · 0 · 16 | 140 |
| mouse | b | 0 · 0 · 80 | 4,619 | 9 · 0 · 65 | 2,107 |
| mouse | c | 26 · 0 · 87 | 800 | 26 · 0 · 87 | 834 |
| mouse | d | 0 · 0 · 4 | 23 | 0 · 0 · 2 | 2,015 |
| mouse | e now | 9 · 0 · 0 | 232 | 142 · 783 · 2 | 8 |
| mouse | e decided | 2 · 0 · 0 | 239 | 0 · 0 · 2 | 12,548 |
| human | harness | 3 · 0 · 41 | 111 | 3 · 0 · 39 | 2,939 |
| human | astar | 6 · 0 · 20 | 7,026 | 2 · 0 · 14 | 10,824 |
| human | a | 2 · 0 · 14 | 262 | 0 · 0 · 12 | 486 |
| human | b | 2 · 0 · 34 | 4,071 | 2 · 0 · 8 | 2,053 |
| human | c | 10 · 0 · 66 | 690 | 11 · 0 · 66 | 712 |
| human | d | 0 · 0 · 0 | 11 | 6 · 0 · 0 | 2,560 |
| human | e now | 7 · 0 · 0 | 785 | 8 · 789 · 2 | 8 |
| human | e decided | 0 · 0 · 0 | 792 | 0 · 0 · 2 | 12,617 |
| rhesus | harness | 0 · 0 · 61 | 607 | 0 · 0 · 49 | 2,817 |
| rhesus | astar | 0 · 0 · 6,800 | 23,380 | 0 · 0 · 0 | 23,492 |
| rhesus | a | 0 · 0 · 1,682 | 496 | 0 · 0 · 12 | 854 |
| rhesus | b | 0 · 0 · 25 | 3,638 | 2 · 0 · 2 | 2,311 |
| rhesus | c | 0 · 0 · 144 | 794 | 0 · 0 · 140 | 1,784 |
| rhesus | d | 0 · 0 · 4 | 5 | 0 · 0 · 4 | 1,655 |
| rhesus | e now | 7 · 0 · 0 | 355 | 8 · 786 · 2 | 2 |
| rhesus | e decided | 0 · 0 · 0 | 363 | 0 · 0 · 2 | 12,408 |
| cynomolgus | harness | 0 · 0 · 0 | 21 | 0 · 0 · 0 | 351 |
| cynomolgus | astar | 0 · 0 · 672 | 3,966 | 0 · 0 · 0 | 3,966 |
| cynomolgus | a | 0 · 0 · 350 | 2,148 | 0 · 0 · 0 | 2,286 |
| cynomolgus | b | 0 · 0 · 0 | 16 | 0 · 0 · 0 | 1,765 |
| cynomolgus | c | 0 · 0 · 0 | 8 | 0 · 0 · 0 | 342 |
| cynomolgus | d | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 1,596 |
| cynomolgus | e now | 7 · 0 · 0 | 139 | 9 · 755 · 2 | 0 |
| cynomolgus | e decided | 0 · 0 · 0 | 146 | 0 · 0 · 2 | 12,140 |

Against `d4c2dda` E's judge `now` counts the expected errors of the decided rules as errors (about 790 per organism, all expected under the judge `decided`), and pig's D2 lists of 500 names (25), human's D (6), pig's `c2_oldver` accessions (row I140, the accession versions: round 2) and B8's `S100A10;` (stray punctuation in Ensembl's value) are the losses `next` had already (the reading of `next`'s modest run). The new sets of the round, before (`next`) and after (`r1`), the rate right / (right + wrong + error) over their strata: B9 (every display name at its own release) pig 0.864 -> 0.999, mouse 0.730 -> 0.951, human 0.797 -> 0.984, rhesus 0.841 -> 0.999, cynomolgus 1.000 -> 1.000; E8's batch against single pig 0.986 -> 0.997, mouse 0.983 -> 0.997, human 0.902 -> 1.000, rhesus 0.974 -> 1.000, cynomolgus 0.977 -> 1.000; E8's one-release pairs pig 0.794, mouse 0.927, human 0.891, rhesus 0.848, cynomolgus 0.777 -> 1.000 everywhere. The known limitations the round closes, measured on the modest selection of `r1` (`results/limitations/`, `b_r1`, `e_r1`): `L-display-name-not-held` 2 answers in a stratum of 1,000 (pig, human, rhesus; on `next`, exhaustive: pig 98,738, mouse 27,394, human 32,295, rhesus 89,181), `L-one-release-reading` none, `L-batch-copies` pig 16, mouse 20, human 2 answers (made-unique copies of a list whose single call reads otherwise). Mouse's B9 lines of UniProt names at early releases keep answering other genes as on `next` (`L-display-name-other`, mouse 972 answers of 7,074).

What is left, by cause (the families' known limitations):

- **Bridges where Ensembl's history gives no successor**, more of them than on `next` because the names and NCBI ids of early releases are held at their releases now ([`display_names_at_release`](rules.md#display_names_at_release)) and the mitochondrial genes of the reused ids bridge to their new ids ([`no_mito_nuclear_step`](rules.md#no_mito_nuclear_step)): across an assembly change that no chain file of Ensembl maps (pig, the macaques) undetermined (`L-no-chain`); at the query's place, on one assembly or mapped by a chain (mouse, human), undetermined (`L-bridge-at-place`); at another place by the chain (mouse: a name or NCBI id that a gene family shares, e.g. `ENSMUSG00000078334` from 60 to 115 through `Gm5923` to `ENSMUSG00000096645`, which the chain puts elsewhere), wrong (`L-bridges`: across assemblies idtrack compares chromosome names only, not places).
- **A name several genes carry**: the best answer is the gene that more sources name at that release (rule [`support_ranks_shared_name`](rules.md#support_ranks_shared_name); NCBI's labels are held at early releases now), undetermined (B's `support_from_xrefs`); and the ties row I144 counts (below).
- **E's judge `now`**: the 7 lines per organism of a batch item without a value (`e3_dec` gives the decided judge's `MissingValueError`), and on mouse one E8 relation (`Six3os1.1`: rule [`one_release_reading`](rules.md#one_release_reading) reads the single call at release 75 by its separators, so in the whole-release list it is the made-unique copy of `Six3os1`, rule [`made_unique_copies`](rules.md#made_unique_copies); E8's relation expects them equal).
- **B's walk of a mitochondrial name**: human `MT-ATP6` from 53 (B1): B's truth walks the gene's id into the nuclear gene Ensembl gave it at 55; `r1` answers the mitochondrial MT-ATP6 of 115 (`all` adds MT-CO3 and MT-ATP8, which Ensembl's cross-references give that id at 54 too); rule [`no_mito_nuclear_step`](rules.md#no_mito_nuclear_step) is not applied to B's lines; one line, not revised.
- **A later reading decides**: pig `MIC-1` (C's clone name to 89; an HGNC synonym at 91) and mouse `100040631` (NCBI's id at 66; a synonym at 87): rule [`one_release_reading`](rules.md#one_release_reading) takes the latest release at which the string has a reading, C's judge reads one database; one line each, not revised.

## The tree `r1b`

`r1b` is `r1` with the graph's format changes of the first round's second part (commit `cbf8603`, the data root's `[trees.r1b]`, its graphs `graphs/r1b/<organism>/`): the version Ensembl's cross-reference attached is kept per link and release, so an accession's versions are spellings of its links (rule [`accession_versions`](rules.md#accession_versions), row I140), and a database whose stored number counts its records' edits (UniProt's Swiss-Prot and TrEMBL on the five graphs) has its typed versions not judged; and a gene's display name at a release ranks first among the genes that carry the string as a name there (row I144). The graphs' nodes, links, history and places are `r1`'s. Its modest run (2026-11-21, the same selections and jobs as `r1`'s) took 52 minutes (pig), 1 h 48 (mouse), 2 h 34 (human), 56 minutes (rhesus) and 30 minutes (cynomolgus), the kernel's peaks 3.3, 5.1, 5.7, 2.4 and 1.5 GB; nothing of other runs was reused. `r1` and `d4c2dda` were not run: their stored answers were scored again under the same judge (`modest organism --score-only` on links to their answers).

Reading the losses found family C's judge lagging rule `accession_versions` and three faults of its own truth; the revision sets (`version_not_judged`, `version_above_latest`, `versions_in_window`) are listed in [family C](families/c.md#the-sets): UniProt's typed versions read as the accession; a version above every version had read as the latest where the accession is attached at the graph's last release, and as an unknown version where it is not (row I150; `r1b` reads the latest there too); the versions counted up to the graph's last release (the generator had counted Ensembl's release 116, so a version written only there made a gap below it); a versioned line read as the line of the version it reads, with every correction published for that line; and a versioned string that is a gene's display name (RefSeq's labels before release 64) read as that name. And the support of rule [`support_ranks_shared_name`](rules.md#support_ranks_shared_name) counted synonyms among the sources of a name, which idtrack does not: the corrections `support_names` of [the harness](harness.md#the-corrections-of-the-judge) and [family B](families/b.md) count the name sources alone. Every part of every tree was then scored again and compared, against `r1` (`<modest>/r1b/<organism>/<part>/moved_vs_r1.tsv`) and against `d4c2dda` (`moved_vs_d4c2dda.tsv`):

| organism | part | vs `r1`: wrong · error · undetermined | gains | vs `d4c2dda`: wrong · error · undetermined | gains |
| --- | --- | ---: | ---: | ---: | ---: |
| pig | harness | 0 · 0 · 0 | 128 | 0 · 0 · 44 | 6,358 |
| pig | astar | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 180 |
| pig | a | 0 · 0 · 0 | 0 | 0 · 0 · 6 | 918 |
| pig | b | 1 · 0 · 0 | 49 | 3 · 0 · 16 | 2,162 |
| pig | c | 0 · 0 · 0 | 1,969 | 2 · 0 · 188 | 4,258 |
| pig | d | 0 · 0 · 0 | 0 | 25 · 0 · 2 | 1,666 |
| pig | e now | 0 · 0 · 0 | 0 | 7 · 806 · 2 | 2 |
| pig | e decided | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 12,470 |
| mouse | harness | 0 · 0 · 0 | 373 | 6 · 0 · 851 | 2,058 |
| mouse | astar | 0 · 0 · 0 | 0 | 6 · 0 · 14 | 1,090 |
| mouse | a | 0 · 0 · 0 | 0 | 0 · 0 · 16 | 140 |
| mouse | b | 0 · 0 · 0 | 186 | 8 · 0 · 55 | 2,160 |
| mouse | c | 0 · 0 · 0 | 1,962 | 26 · 0 · 87 | 2,796 |
| mouse | d | 0 · 0 · 0 | 5 | 0 · 0 · 2 | 2,020 |
| mouse | e now | 0 · 0 · 0 | 1 | 141 · 783 · 2 | 8 |
| mouse | e decided | 0 · 0 · 0 | 1 | 0 · 0 · 2 | 12,549 |
| human | harness | 0 · 0 · 0 | 121 | 3 · 0 · 39 | 3,060 |
| human | astar | 0 · 0 · 0 | 0 | 2 · 0 · 14 | 10,824 |
| human | a | 0 · 0 · 0 | 0 | 0 · 0 · 12 | 486 |
| human | b | 0 · 0 · 0 | 215 | 2 · 0 · 6 | 2,099 |
| human | c | 0 · 0 · 0 | 1,150 | 10 · 0 · 66 | 1,864 |
| human | d | 0 · 0 · 0 | 9 | 5 · 0 · 0 | 2,568 |
| human | e now | 1 · 0 · 0 | 0 | 9 · 789 · 2 | 8 |
| human | e decided | 1 · 0 · 0 | 0 | 1 · 0 · 2 | 12,617 |
| rhesus | harness | 0 · 0 · 0 | 46 | 0 · 0 · 49 | 2,863 |
| rhesus | astar | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 23,492 |
| rhesus | a | 0 · 0 · 0 | 0 | 0 · 0 · 12 | 854 |
| rhesus | b | 0 · 0 · 0 | 66 | 0 · 0 · 2 | 2,336 |
| rhesus | c | 0 · 0 · 0 | 1,877 | 0 · 0 · 140 | 3,661 |
| rhesus | d | 0 · 0 · 0 | 0 | 0 · 0 · 4 | 1,655 |
| rhesus | e now | 0 · 0 · 0 | 0 | 8 · 786 · 2 | 2 |
| rhesus | e decided | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 12,408 |
| cynomolgus | harness | 0 · 0 · 0 | 61 | 0 · 0 · 0 | 412 |
| cynomolgus | astar | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 3,966 |
| cynomolgus | a | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 2,286 |
| cynomolgus | b | 0 · 0 · 0 | 21 | 0 · 0 · 0 | 1,780 |
| cynomolgus | c | 0 · 0 · 0 | 1,500 | 0 · 0 · 0 | 1,842 |
| cynomolgus | d | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 1,596 |
| cynomolgus | e now | 0 · 0 · 0 | 1 | 8 · 755 · 2 | 0 |
| cynomolgus | e decided | 0 · 0 · 0 | 1 | 0 · 0 · 2 | 12,141 |

Against `r1` no answer of any part went from right to wrong, to an error or to undetermined but these two, faults of `r1b`:

- **pig B1 `C/EBP` at 78 to 115** (`best`): `r1b` answers ENSSSCG00000002866 alone where the truth expects the tie with ENSSSCG00000034207 (`r1`'s answer). At 78 `C/EBP` is the display name of ENSSSCG00000002866 and of ENSSSCG00000028739, and a UniProt gene name of ENSSSCG00000024816; Ensembl's history merges 28739 and 24816 into 34207 at 90 (equal scores). `r1b` keeps one path per gene and takes the candidate's display mark (and its support) from that path's entry gene, here 24816, so 34207 ranks below 2866 although a gene whose display name the string is leads to it: a fault of `r1b` (row I156: the display mark and the support taken over all the entry genes of the paths that reach a gene), reproducer `convert("C/EBP", from_release=78, to_release=115)` on pig.
- **human E `3.14159`** (`e2_junk_num`, judges `now` and `decided`): `r1b` answers ENSG00000224777 (`version_newer_than_graph`, read at 51): Ensembl's cross-references give the Havana translation accession `3` with version 1 to ENSP00000318226 at 50 and 51, and `r1b` reads a version above every version attached as the latest also for an accession no longer attached at the graph's last release: a fault of `r1b` (row I150; the rule [`accession_versions`](rules.md#accession_versions) answers it with nothing), the class of hunt IN-9 (numbers read as old numeric Havana ids).

Against `d4c2dda` the losses are those of `r1` (the classes of its section above) but: pig's `c2_oldver` accessions (row I140) are right now; pig B `C/EBP` and human E `3.14159` (above); family C's versioned lines that `d4c2dda` answered by cutting the version are judged by the rule as built; and family C's lines of row I150's class (below). Mouse B1 `Mtcp1` at 78 (`best`), which `r1` already answered so, is undetermined: idtrack ranks ENSMUSG00000031200 first by its support (a Vega gene label besides the names both carry), and the support of the name sources alone (`b1_snames1`) says so; the judge had counted a UniProt gene-name synonym of ENSMUSG00000090110 as a source and expected the tie.

The round's faults, before (`r1`) and after (`r1b`), on the modest selection (right / (right + wrong + error); answers of both views): family C's version strata (`c2_*_v`, `c2_oldver`, `c2_ver_*` and the same of c3; 1,000 answers per stratum or the whole stratum) pig 0.794 -> 0.935 (2,855 wrong -> 904), rhesus 0.799 -> 0.934 (2,793 -> 924), cynomolgus 0.803 -> 0.911 (2,744 -> 1,244), mouse 0.629 -> 0.875 (2,962 -> 1,000), human 0.782 -> 0.993 (1,186 -> 36); what `r1b` leaves wrong there is row I150's class (a version above every version of an accession no longer attached at the graph's last release, read as the latest: `c2_ver_higher_retired`, and gap lines whose upper versions lie beyond the window: cynomolgus 414 answers) and, in mouse, 6 answers of bridges where Ensembl's history gives the accession's genes no successor (`L-bridges`). B9 (row I144) pig 0.9993 -> 0.9995, mouse 0.9506 -> 0.9515, human 0.9844 -> 0.9868, rhesus 0.9993 -> 0.9995, cynomolgus 1.000 -> 1.000. The ties round 1 created (an early display name beside NCBI's symbol on another gene's translation: answers right on `next`, undetermined on `r1` within the line's wider set): the `best` view is right again for every one (mouse 6, human 6, rhesus 9 answers); the `all` view returns every gene that carries the string and stays undetermined (mouse 14, human 16, rhesus 12). `L-accession-versions` is closed but for row I150's class (family C's [known limitations](families/c.md#known-limitations)). The hunts' reproducers hold or changed as on `r1` but AT-4 (an old RefSeq version asked without a release), which changed on pig and mouse: `XM_003129161.3` answers ENSSSCG00000026015 and `NM_001044572.1` ENSSSCG00000058470, the genes Ensembl attached those versions to.

The wrong answers round 1 left by bridges through names now held that Ensembl's chain file puts elsewhere (`L-bridges`: mouse harness 9, A* 6, A 4, C 25; human harness 3, A* 6, A 2, C 10) are all wrong on `r1b` as on `r1`.


## The tree `r2`

`r2` is idtrack after the fix of its second round, the last round of fixes (its runs at commit `978f9a0`, whose package is that of the release 0.1.0, `b5d4fb4`, which the data root's `[trees.r2]` names: [Trees](../README.md#trees); its graphs `graphs/r2/<organism>/`): Ensembl's chain files carry a gene's place across an assembly change and places are compared by overlap with strand (row I147); a copy on a patch or haplotype stands for its assembly exception's reference region (I148); a database Ensembl renamed or refiled is the selected database it continues (rule [`database_filing`](rules.md#database_filing), I149 and I157); a version above every judged one is the latest only where the accession is attached at the graph's last release (I150); a target in a final database is an id of the answered gene (I151); the display mark is ranked by the spelling typed (I152); every id form of an LRG record answers as the record (I153); the display mark and the support are taken over every entry gene of the paths that reach a gene (I156); and a link read on a transcript counts for the gene that holds the transcript at the link's releases (I158). Its modest run (2026-11-22, the selections and jobs of `r1b`'s) took 51 minutes (pig), 1 h 52 (mouse), 2 h 20 (human), 54 minutes (rhesus) and 27 minutes (cynomolgus), the kernel's peaks 3.3, 5.1, 5.7, 2.4 and 1.5 GB; nothing of other runs was reused. `r1b` and `d4c2dda` were not run: their stored answers were scored again under the same judge (`modest organism --score-only` on links to their answers; `r1b`'s scoring before the revision below equals its stored one in every part).

Reading the losses found family C's judge lagging rule `database_filing`: the versions RefSeq's databases wrote under their earlier names, and the values Ensembl filed under them, were not the accessions' (`c2_ver_lower` said `XM_001082036.1` was never had: Ensembl's `RefSeq_dna_predicted` wrote it at rhesus 55 to 62). The revision `database_continuity` ([family C](families/c.md#the-sets): pig 6 lines, mouse 1,744, human 970, rhesus 2,865) reads Ensembl's cross-references for the databases that continue an included one (`truth_c.continuations`, equal to the graphs' own record of the filing); on the modest selection `r2` answers every revised line as the revision reads it but where the locus rule judges a bridge. C was then scored again on all three trees and every part compared, against `r1b` (`<modest>/r2/<organism>/<part>/moved_vs_r1b.tsv`, `compare_vs_r1b.tsv`) and against `d4c2dda` (`moved_vs_d4c2dda.tsv`, `compare.tsv`):

| organism | part | vs `r1b`: wrong · error · undetermined | gains | vs `d4c2dda`: wrong · error · undetermined | gains |
| --- | --- | ---: | ---: | ---: | ---: |
| pig | harness | 10 · 0 · 0 | 248 | 10 · 0 · 41 | 6,603 |
| pig | astar | 6 · 0 · 0 | 2 | 6 · 0 · 0 | 182 |
| pig | a | 0 · 0 · 0 | 72 | 0 · 0 · 6 | 990 |
| pig | b | 0 · 0 · 0 | 11 | 2 · 0 · 16 | 2,172 |
| pig | c | 0 · 0 · 2 | 916 | 2 · 0 · 186 | 5,170 |
| pig | d | 0 · 0 · 0 | 2 | 25 · 0 · 0 | 1,666 |
| pig | e now | 0 · 0 · 0 | 0 | 7 · 806 · 2 | 2 |
| pig | e decided | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 12,470 |
| mouse | harness | 50 · 0 · 187 | 13 | 53 · 0 · 1,041 | 2,071 |
| mouse | astar | 0 · 0 · 2 | 16 | 0 · 0 · 18 | 1,102 |
| mouse | a | 0 · 0 · 2 | 16 | 0 · 0 · 16 | 154 |
| mouse | b | 4 · 0 · 11 | 8 | 11 · 0 · 57 | 2,164 |
| mouse | c | 0 · 0 · 6 | 1,400 | 20 · 0 · 95 | 3,396 |
| mouse | d | 1 · 0 · 0 | 0 | 0 · 0 · 2 | 2,019 |
| mouse | e now | 2 · 0 · 0 | 0 | 141 · 783 · 2 | 8 |
| mouse | e decided | 2 · 0 · 0 | 0 | 0 · 0 · 2 | 12,549 |
| human | harness | 15 · 0 · 15 | 22 | 13 · 0 · 42 | 3,065 |
| human | astar | 3 · 0 · 7 | 3,520 | 4 · 0 · 14 | 14,336 |
| human | a | 4 · 0 · 21 | 210 | 1 · 0 · 14 | 674 |
| human | b | 0 · 0 · 20 | 4 | 2 · 0 · 6 | 2,082 |
| human | c | 0 · 0 · 4 | 120 | 5 · 0 · 71 | 1,816 |
| human | d | 0 · 0 · 0 | 1 | 5 · 0 · 0 | 2,569 |
| human | e now | 0 · 0 · 0 | 1 | 8 · 789 · 2 | 8 |
| human | e decided | 0 · 0 · 0 | 1 | 0 · 0 · 2 | 12,617 |
| rhesus | harness | 0 · 0 · 0 | 0 | 0 · 0 · 49 | 2,863 |
| rhesus | astar | 10 · 0 · 0 | 0 | 10 · 0 · 0 | 23,492 |
| rhesus | a | 0 · 0 · 0 | 0 | 0 · 0 · 12 | 854 |
| rhesus | b | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 2,336 |
| rhesus | c | 0 · 0 · 114 | 1,508 | 0 · 0 · 146 | 3,933 |
| rhesus | d | 0 · 0 · 0 | 0 | 0 · 0 · 4 | 1,655 |
| rhesus | e now | 0 · 0 · 0 | 0 | 8 · 786 · 2 | 2 |
| rhesus | e decided | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 12,408 |
| cynomolgus | harness | 10 · 0 · 0 | 0 | 10 · 0 · 0 | 412 |
| cynomolgus | astar | 10 · 0 · 0 | 0 | 10 · 0 · 0 | 3,966 |
| cynomolgus | a | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 2,286 |
| cynomolgus | b | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 1,780 |
| cynomolgus | c | 0 · 0 · 0 | 1,244 | 0 · 0 · 0 | 3,086 |
| cynomolgus | d | 0 · 0 · 0 | 0 | 0 · 0 · 0 | 1,596 |
| cynomolgus | e now | 0 · 0 · 0 | 0 | 8 · 755 · 2 | 0 |
| cynomolgus | e decided | 0 · 0 · 0 | 0 | 0 · 0 · 2 | 12,141 |

Against `r1b` the answers that went from right to wrong are of these causes (answers of both views):

- **A split gene converted to a name database** (row I151 as built: the answer is first ranked as Ensembl genes, and a target is an id of that gene): where Ensembl's history splits the gene and the descendant ranked first has no id in the database, the result is `no_target` with the genes, in the view `all` too, although other descendants carry ids there; `r1b` and `d4c2dda` answer those ids, as the judges of A*'s `astar_db` and the harness's `travel_split` expect (rule [`tie_in_final_database`](rules.md#tie_in_final_database)). Pig 16 answers (harness 10, A* 6), human 8 (harness 6, A* 2), rhesus 10 (A*), cynomolgus 20 (harness 10, A* 10), mouse none; e.g. cynomolgus `ENSMFAG00000002704` from 95 to 105 to `UniProtKB Gene Name`: nine descendants, three of them carry H2BC1, H2BC11 and H2BC17, the first by annotation, `ENSMFAG00000053830`, none. Whether row I151 is meant to answer a split so is open (a question of this reading); it is the known-limitations row `L-split-no-target`.
- **RefSeq's gene-level accessions of mouse 80 and 81 as RefSeq_mRNA** (a fault of the filing's form rule, `L-filing-forms`): Ensembl's `RefSeq_dna` on mouse genes at 80 and 81 holds `NM_`, `XM_`, `XR_`, `NR_` and `NG_` accessions; the filing reads it as RefSeq_mRNA and sends the forms that never moved to that main continuation, which has no identifier scheme (`NM_` is 0.89 of its accessions over the releases, Ensembl wrote gene labels and NCBI ids there), although `XM_` is RefSeq_mRNA_predicted's scheme (1.000) and `NG_` (RefSeqGene) no mRNA. Mouse harness `old_gene_to_db` to RefSeq_mRNA at 80: 41 answers that add `XM_`, `XR_` or `NG_` accessions to the `NM_` that BioMart's table gives.
- **The judge lagging rules of the round, not revised** (one to six lines each): the harness's BioMart tables at mouse 80 (3 answers: further `NM_` accessions Ensembl's gene-level `RefSeq_dna` gives the gene, rule `database_filing`) and at 75 (mouse 6, human 3 and human 80's `4540`: `NR_` accessions and an NCBI gene id that Ensembl filed under RefSeq_mRNA, which the filing reads by their scheme as RefSeq_ncRNA's and NCBI's); family B's `BASIC` (mouse; an MGI synonym on translations at 48 to 50, held by the filing and read without letter case for a one-case string, rule [`single_case_symbols`](rules.md#single_case_symbols); B's check of the core records reads the spelling) and `OTTMUSG00000010207` at 48 (row I152: the authority's display name, MGI's symbol, ranks before the Havana label of another gene); D's list of 500 names of mouse 48 (inferred 48 to 50: MGI's names on translations are held at 49 and 50 too); A7's `LRG_3_t1` at 58 (human; row I153: the gene).
- **A support counted in another letter case** (no decision covers it): mouse `Adora3` from 62 to 115 (`best`): UniProt's gene name `ADORA3` on ENSMUSG00000000562's translations, held at 62 now, counts for the typed mixed-case `Adora3`, so that gene ranks before ENSMUSG00000074344 where B's judge expects the tie; `r1b` and `r2` count it so from 63 on (mouse's B lines of shared names wrong by a missing gene: 15 on `r1b`, 16 on `r2`).
- **Bridges where Ensembl's history gives no successor** (`L-bridges`, through rows I147 and I148): human 9 answers (harness 5, A 3, A* 1), e.g. `ENSG00000204645` (SSX4, X on GRCh37) from 75: `r2` bridges at 75 to the patch copy ENSG00000268009 at the same place, whose id Ensembl gave the reference gene at 76, where the chain maps the old place onto the other copy of the duplication (61.7 % aligned, the other strand); and mouse E8's batch item `Igh-6` of release 53, a copy of `Igh` in the batch (`L-batch-copies`).

The answers that went from right to undetermined (all parts) are bridges the judge cannot place, at the query's place (`L-bridge-at-place`: mouse 185, human 25) or across an assembly change no chain maps (`L-no-chain`: rhesus 114, human 21, mouse 12, pig 2), the supports of shared names (human 10, mouse 1), mouse's early display-name ties (6, row I144) and answers within `broad` (human 11, mouse 4) (bridges `r2` takes and `r1b` did not: places mapped by chains and copies placed by their exceptions, rows I147 and I148, and names held at their releases, row I149). Against `d4c2dda` the losses are these and `r1b`'s (its section above) but pig B `C/EBP`, now the tie (row I156), and human E `3.14159`, now nothing (row I150).

The known-limitations rows of `r2` (`results/limitations/f0_r2.tsv`, `astar_r2.tsv`, `b_r2.tsv`, `e_r2.tsv`, from the outcomes): `L-split-no-target` pig 10 of the 22,920 `travel_split` answers and 6 of 1,000 `db` answers, human 6 of 23,586 and 2 of 5,070, rhesus 10 of 3,126 `db`, cynomolgus 10 of 9,510 and 10 of 1,218; `L-filing-forms` mouse 41 of 8,700 `old_gene_to_db` answers; `L-display-name-other` and `L-display-name-not-held` (B9) as on `r1b`; `L-batch-copies` (E8) as on `r1b` but mouse 22 (`r1b` 20: `Igh-6` above).

The round's faults, before (`r1b`) and after (`r2`), on the modest selection (right / (right + wrong + error), both views): family C's version strata pig 0.935 -> 1.000 (904 wrong -> 0), mouse 0.825 -> 0.9997 (1,402 -> 2, one line of a bridge, `L-bridges`), human 0.993 -> 1.000 (36 -> 0), rhesus 0.892 -> 1.000 (1,504 -> 0), cynomolgus 0.911 -> 1.000 (1,244 -> 0); B9 unchanged (pig 0.9995, mouse 0.9515, human 0.9868, rhesus 0.9995, cynomolgus 1.000). Round 1's wrong bridges through names held now that Ensembl's chain puts elsewhere (65 lines, all wrong on `r1b`): 10 right, 15 undetermined, 40 wrong on `r2`. The hunts' reproducers hold or changed as on `r1b` but AT-5 (mouse; a copy on a patch bridged far from its place), which changed to nothing; round 2's audit reproducers R2-1 to R2-6 and H5, asked on `r2`'s graphs, give the answers Ensembl's records expect (R2-2's three symbols of genes without one: `no_target` with the gene; R2-5's batch item: the record's id form, `lrg:exact`).

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

# a family's part on another published selection than its modest one (its panel: what `iterate` asks where one exists)
idtrack-repro modest organism --tree <tree> --organism <o> --out <folder>/<tree> --selections a=panel_v1

# compare the tree with the reference, part by part (exit 1 on a regression: right -> wrong or right -> error)
idtrack-repro modest compare <folder>/d4c2dda <folder>/<tree> --out <folder>/<tree>/compare.tsv
```

`<tree>` is a configured tree (`[trees.<tree>]` in `idtrack-repro.toml`: its checkout, Python and graphs per organism); `<folder>` is the folder that holds the trees' runs, the data root's `modest/` on the data root of record ([Where the runs lie](#where-the-runs-lie)). A run writes `<out>/<organism>/<part>/` (answers and scored tables of each part: `harness`, the families, `e/scores_now`, `e/scores_decided`, `hunts.json`) and `<out>/<organism>/modest.json` (the parts in order with their exit codes and seconds, the answers each runner took from other runs of the tree, a part killed by a signal). One organism's run, without the job runner: `idtrack-repro modest organism --tree <tree> --organism <organism> --out <folder>`.
