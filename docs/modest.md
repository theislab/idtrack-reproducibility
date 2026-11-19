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

`idtrack-repro modest size` computes it from the sets' manifests (no set is read). With the sets published on 2026-11-19 and 500 lines per stratum:

| part | lines (five organisms) | strata |
| --- | ---: | ---: |
| harness | 2,615,211 | (every line) |
| A* | 42,500 | 85 |
| A | 70,176 | 142 |
| B | 57,077 | 139 |
| C | 98,515 | 221 |
| D | 59,947 | 166 |
| E | 253,624 | 1,002 |
| hunt reproducers | 29 findings | |
| **total** | **3,197,050 lines** | |

Time: the harness's answers of `d4c2dda` took 116 to 372 seconds of query time per organism (`query_seconds_total` of the stored scorings: pig 204, mouse 372, human 346, rhesus 170, cynomolgus 116); the families' samples add about 580,000 lines (most asked under both views: two calls each on `d4c2dda`, one from `71690ab` on); loading a graph takes seconds to a minute and is repeated per part, and scoring the harness's sets takes minutes per organism. One job per organism is expected to finish in one to two hours (an estimate from these measurements, not yet a measured run); its memory is the organism's from the configuration plus 16 GB.

## What it can and cannot detect

- **Can**: every change of an answer among the harness's 2.6 million lines (the judge's fixed queries over ids, names, external identifiers, history kinds and the travel matrix), each judged right or wrong as the exhaustive harness judges it; every change of a known finding of the hunts (a reproducer exits 1 and names the field); in every stratum of every family, a change that makes a share p of its lines wrong is missed with probability (1 - p)^500, under 5 % for p of 0.6 % or more (for a stratum of fewer than 500 lines: every line, so every change). A stratum's rate is estimated to within about ±4.4 percentage points at worst (the Wilson 95 % interval on 500 lines at a rate of one half; much less near 100 %).
- **Cannot**: a regression confined to fewer than about 0.6 % of a family's stratum outside the harness's sets (it may be missed), a regression in a line never sampled (most of A's 400 million and A*'s 290 million rows are not), and the families' parts run outside the common runner (E's determinism and concurrency runs, metamorphic relations needing two calls, budgets, the harmonizer, the web layer; A's graph fidelity A0; H's timings). These need the panel ([Measurement design](methods.md#measurement-design)) or the exhaustive pass of the final tree.

## Commands

```
# once: build the families' selections (reads each family's sets once; run it as a job: about an hour or two)
idtrack-repro modest build                      # --per-stratum 500 by default; --dry counts without writing

# measure a tree: one job per organism (the configured runner); the reference too, once, for the comparison
idtrack-repro modest run --tree d4c2dda --out <folder>/d4c2dda
idtrack-repro modest run --tree <tree> --out <folder>/<tree>

# compare the tree with the reference, part by part (exit 1 on a regression: right -> wrong or right -> error)
idtrack-repro modest compare <folder>/d4c2dda <folder>/<tree> --out <folder>/<tree>/compare.tsv
```

`<tree>` is a configured tree (`[trees.<tree>]` in `idtrack-repro.toml`: its checkout, Python and graphs per organism). A run writes `<out>/<organism>/<part>/` (answers and scored tables of each part: `harness`, the families, `e/scores_now`, `e/scores_decided`, `hunts.json`) and `<out>/<organism>/modest.json` (the parts in order with their exit codes and seconds). One organism's run, without the job runner: `idtrack-repro modest organism --tree <tree> --organism <organism> --out <folder>`.
