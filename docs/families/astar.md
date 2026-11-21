# A*: the hardest lineages

A* asks idtrack to travel the gene ids whose history is hardest to follow, from every release where each id exists to every release of the window, and reports them apart, per pattern and per id. The exhaustive family [A](a.md) measures history travel on average; A* measures it where it is hardest, where a regression would hide in an average. Its code is `idtrack_reproducibility.families.astar`; its data are `<families>/astar/` (`selection/`, `queries/<organism>/`, `runs/<tree>/<organism>/`); the design is in [the methods](../methods.md#a-the-hardest-lineages). It is built on [the common machinery](f0.md).

## What it measures

For every selected id, the answer of idtrack for the id asked from every release where it exists to every release of the window (forward, backward and in place, within an assembly and across an assembly handoff), as stable id and with its version, without `from_release`, and to a name database; and for the transcripts and translations of the selected genes, the gene they belong to. Each line is scored under `best` and under `all`. The rate is right / (right + wrong + error), with the cluster bootstrap over the id's lineage (the component of cross-id events it belongs to); undetermined answers (the truth cannot decide) and timeouts apart.

A new tree of idtrack is read on A* first: an answer that moves from right to wrong here is explained before anything else of the tree is read.

## Selection

The ids are selected from Ensembl's own files only (the GTFs of every release, the release-115 core history tables `stable_id_event` and `mapping_session`, the core gene tables the GTFs lack, the assembly chain files and exceptions), deterministically, by `idtrack_reproducibility.families.astar.selection`. Every gene id of the window gets its patterns; the selection is the union over the patterns of every id at or above the pattern's threshold, plus two fixed lists that only add ids, never truth (`astar/fixed_ids.tsv` of the package: the ids named by the weakness hunts of history and regions ([weakness hunts](../methods.md#weakness-hunts)), and the 65 human ids behind the wrong answers of the travel matrix's weakest cell, backward within GRCh37 on `d4c2dda`).

| pattern | how it is found | threshold |
| --- | --- | --- |
| `lineage_ge5` | connected components of the graph of cross-id gene events (old and new id differ, both given) of the release-115 `stable_id_event` | every id of a component of 5 ids or more |
| `many_events` | mapping sessions with an event of the id that is neither to itself nor its birth; distinct partner ids in cross-id events | 3 sessions or more, or 3 partners or more |
| `split_then_merged_back` | two event paths from one id, by different first steps, meeting again in a later id: the id and the ids where its branches meet | all |
| `merged_then_split` | an id with 2 or more predecessors and 2 or more successors over the cross-id events | all |
| `retired_and_back` | an event to nothing, then the same stable id again in a later GTF, or an event from nothing in a later session (rule [`resurrection`](../rules.md#resurrection)) | all |
| `gap` | absent from the GTF (and the core gene table, where fetched) of a release between two releases that have it | all |
| `gtf_lacks` | in the core database's gene table of a release whose GTF lacks it (human 48 to 51, cynomolgus 112) | all |
| `version_jump` | the GTF version rises by more than one between consecutive releases, or a self event raises it by more than one | all |
| `version_fall` | the GTF version falls between consecutive releases, or a self event lowers it | all |
| `moved` | within an assembly, on another sequence or the other strand at another release; across assemblies, on two different chromosomes (haplotype and patch sequences left out: their names change by convention) | all |
| `handoff_lineage` | a member of a lineage with a cross-id event in a session that crosses an assembly handoff | all |
| `handoff_moved` | present on both sides of a handoff whose chain file exists, and its place at the last release of one assembly, mapped through Ensembl's chain, does not overlap its place at the first release of the next | all |
| `rename_chain_ge3` | a member of a chain of ids, each with one successor and one predecessor | chains of 3 or more |
| `successor_overlaps` | a cross-id event whose new id, at the first GTF release after the session, overlaps another living gene on the same sequence and strand: the new and the old id | all |
| `event_contradicts_gtf` | an event whose old id is absent from the GTF of the last release before the session takes effect, whose new id is absent from the GTF of the first release from it, or whose version differs from that GTF's (rule [`event_version`](../rules.md#event_version)'s case) | all |
| `mislabelled_session` | an id with a cross-id event in a session whose new release is not after its old release (pig 106 -> 104, human 114 -> 114) | all |
| `non_reference_lineage` | a member of a lineage of 2 or more ids with a gene on a patch or an alternative locus at some release | all |
| `many_versions` | the number of distinct GTF versions of the id | strictly above the 99th percentile (at most 1 % of the ids) |
| `branch_while_living` | a cross-id event from the id in a session it lives across (in the GTF on both sides): Ensembl maps it to another id while it continues (hunt finding hh-1); the id and the ids it branches to | all |
| `hunt_history` | an id named by the history and regions hunts' findings | all |
| `non_reference_copy_ended` | a gene on a haplotype, patch or alternative locus at some release, absent at the target, whose strict walk from its last release reaches nothing (hunt findings R6, J1, J4) | all |
| `par_y` | a gene on Y inside a pseudoautosomal region at some release (hunt finding R4) | all |
| `travel_weak_cell` | human: an id behind a wrong answer of the travel matrix backward within GRCh37 | all |

The selection was made three times, each time with more patterns, and each is kept (`<families>/astar/selection/<organism>.json`, `<organism>.r2.json`, `<organism>.r3.json`, with one row per selected id in the `.tsv.gz` beside each): the first selection; `r2`, which adds `branch_while_living` and the history hunt's ids; `r3`, which adds `non_reference_copy_ended`, `par_y` and the regions hunt's ids. Ids in the window (and the full-matrix lines, id x release of origin x release of the window) of each:

| organism | window | first | `r2` | `r3` |
| --- | --- | ---: | ---: | ---: |
| pig | 56–115 | 8,849 (11,030,100) | 8,849 (11,030,100) | 8,849 (11,030,100) |
| mouse | 48–115 | 11,259 (30,623,188) | 11,260 (30,624,480) | 12,184 (31,635,232) |
| human | 48–115 | 13,109 (31,247,972) | 13,131 (31,310,328) | 15,638 (32,704,192) |
| rhesus | 48–115 | 4,556 (9,320,488) | 4,556 (9,320,488) | 4,560 (9,328,376) |
| cynomolgus | 91–115 | 2,265 (776,575) | 2,265 (776,575) | 2,269 (778,775) |

Ids in the window per pattern in `r3` (an id counts under every pattern it has):

| pattern | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `lineage_ge5` | 5,580 | 5,407 | 4,651 | 3,304 | 1,656 |
| `many_events` | 4,380 | 4,078 | 2,575 | 2,447 | 1,222 |
| `split_then_merged_back` | 4,164 | 1,906 | 1,270 | 1,443 | 0 |
| `merged_then_split` | 1,670 | 913 | 641 | 527 | 0 |
| `retired_and_back` | 0 | 3 | 84 | 0 | 0 |
| `gap` | 0 | 2 | 2 | 0 | 0 |
| `gtf_lacks` | 0 | 0 | 254 | 0 | 47 |
| `version_jump` | 0 | 3,179 | 16 | 0 | 0 |
| `version_fall` | 0 | 3 | 6 | 0 | 0 |
| `moved` | 392 | 337 | 1,359 | 801 | 170 |
| `handoff_lineage` | 6,714 | 5,644 | 7,300 | 3,761 | 2,048 |
| `handoff_moved` | 575 | 107 | 786 | 0 | 0 |
| `rename_chain_ge3` | 0 | 6 | 453 | 0 | 0 |
| `successor_overlaps` | 1,439 | 2,182 | 4,131 | 1,359 | 772 |
| `event_contradicts_gtf` | 0 | 3,181 | 917 | 0 | 0 |
| `mislabelled_session` | 2,204 | 0 | 15 | 0 | 0 |
| `non_reference_lineage` | 0 | 6 | 2,412 | 0 | 0 |
| `many_versions` | 294 | 930 | 727 | 0 | 0 |
| `branch_while_living` | 684 | 0 | 35 | 0 | 0 |
| `hunt_history` | 8 | 6 | 40 | 9 | 7 |
| `non_reference_copy_ended` | 0 | 930 | 3,645 | 0 | 0 |
| `par_y` | 0 | 0 | 71 | 0 | 0 |
| `travel_weak_cell` | 0 | 0 | 65 | 0 | 0 |

Each selection file also records the lineages (pig: 1,045, of which 379 of 5 ids or more, the largest of 225 ids; mouse 1,445, 377, 1,086; human 3,971, 329, 601; rhesus 335, 203, 117; cynomolgus 220, 109, 71), the cut of `many_versions`, and the handoffs checked through the chains (Ensembl publishes no chain for pig's 89 -> 90, rhesus's 85 -> 86 and 97 -> 98, and cynomolgus's 102 -> 103 handoffs).

## The sets

Every line has `strategy: both` and the qid `<kind>|<query>|<from or ->|<to>`; the unit is the id's lineage (the smallest id of its component of cross-id events; the id itself when it has none). The kinds of `gen` (gene ids):

- `astar_from`: the stable id from every release R where it exists to every release t of the window; `astar_from_v` the same with its version at R (where the GTF of R writes versions);
- `astar_default`, `astar_default_v`: without `from_release`, the stable id and each of its versions, to every t (rule [`start_release`](../rules.md#start_release));
- `astar_db`: for every sixth id (crc32 of the id modulo 6 is 0), the stable id from every R to every t with a BioMart table, to the name database of the travel matrix (`HGNC Symbol`, `MGI Symbol`, else `UniProtKB Gene Name`).

The strata are the direction of travel: `same`, `forward`, `backward`, each with `_handoff` when the two releases lie on different assemblies; `default` (no `from_release`); `db`. The transcripts and translations (`gen_tx`, bounded because they cost several times a gene id per call on `d4c2dda`): every transcript (`tx`) and translation (`tl`) of a selected gene from its first and its last release to that release and to 115 (stable and versioned), and one transcript and one translation per gene (seeded) from its first and last release to every release of the window; strata `tx_<direction>`, `tl_<direction>`.

The published sets (manifests `<families>/astar/queries/<organism>/<set>.manifest.json`; a revision set holds only the changed lines, under their qids):

| organism | set | lines | revises | what it holds; why a revision |
| --- | --- | ---: | --- | --- |
| pig | `v1` | 19,359,613 | – | the 8,849 ids of the first selection (also `r3`'s for pig) |
| pig | `v1r1` | 1,209 | `v1` | the corrected judge's `spellings_once` on the `db` lines, published before the generator wrote it |
| pig | `tx1` | 1,479,508 | – | transcripts and translations of the selected genes |
| mouse | `v1` | 54,344,863 | – | the 12,184 ids of `r3` |
| mouse | `tx1` | 2,173,151 | – | transcripts and translations |
| human | `v1` | 54,158,824 | – | the 15,638 ids of `r3` |
| human | `v1r1` | 4,846 | `v1` | the second reading for lines of 2,529 ids generated while it had no record of them |
| human | `v1r2` | 178,103 | `v1`, `v1r1` | the core genes the GTF lacks (48 to 51) written without the core's version where the GTF writes none; copies apart before release 82 read from the assembly exceptions |
| human | `v1r3` | 1,375,610 | `v1` to `v1r2` | the places of the core genes the GTF lacks, as the corrected judge's `core_places` |
| human | `tx1` | 3,408,606 | – | transcripts and translations |
| human | `tx1r1` | 94,950 | `tx1` | `core_places` on the transcript lines |
| rhesus | `v1` | 15,336,920 | – | the 4,556 ids of the first selection |
| rhesus | `tx1` | 724,365 | – | transcripts and translations |
| rhesus | `hunt1` | 14,872 | – | the 4 ids `r3` added |
| cynomolgus | `v1` | 1,688,112 | – | the 2,265 ids of the first selection |
| cynomolgus | `tx1` | 138,073 | – | transcripts and translations |
| cynomolgus | `hunt1` | 4,714 | – | the 4 ids `r3` added |
| rhesus | `v1r1` | 15,784 | `v1` | rule [`no_mito_nuclear_step`](../rules.md#no_mito_nuclear_step) (judge `rules`, `revise_rules.py`): the lines of the lineages of the 5 mitochondrial ids that Ensembl gave to nuclear genes at Mmul_10 |
| rhesus | `tx1r1` | 614 | `tx1` | as `v1r1`, their transcripts and translations |
| rhesus | `v1r2` | 372 | `v1`, `v1r1` | as `v1r1`, the `db` lines whose walk the rule leaves without a successor (the generator makes no `db` line of such a walk; `revise_rules.db_without_successor`): on the nuclear side nothing, on the mitochondrial side open |
| cynomolgus | `v1r1` | 1,298 | `v1` | as rhesus `v1r1` (5 ids, Macaca_fascicularis_6.0) |
| cynomolgus | `hunt1r1` | 1,947 | `hunt1` | as `v1r1` |
| cynomolgus | `tx1r1` | 80 | `tx1` | as `v1r1` |
| human | `v1r4` | 3,462 | `v1` to `v1r3` | as rhesus `v1r1` (5 ids at 48 to 55) |
| human | `tx1r2` | 509 | `tx1`, `tx1r1` | as `v1r4` |

Lines per stratum of each set are in its manifest (`counts`); for example human `v1`: `same` 753,551, `forward` 12,353,556, `backward` 12,707,962, `forward_handoff` 8,822,392, `backward_handoff` 16,604,007, `default` 2,553,536, `db` 363,820.

## Truth

Nothing of idtrack: the truth is the harness's reading of Ensembl's history ([the accuracy harness](../harness.md)) through `astar/truth.py` and `astar/gen.py` (`Gen.truth_at`):

- **the strict walk with living ids** (`expect`): the genes the history carries the id to at the target, where an id present in the GTFs on both sides of a session is kept (rule [`id_lives_on`](../rules.md#id_lives_on)); sessions are crossed where their effect appears in the GTFs, not at their recorded old release (the corrected judge's `sessions_effective`); a step between the mitochondrial genome and a nuclear sequence is not followed, also where the id stays the same (rule [`no_mito_nuclear_step`](../rules.md#no_mito_nuclear_step), `truth.TruthA`; the second reading is not compared on those lineages, since REST reads the steps the rule refuses);
- **the broad walk with returned ids** (`broad`): also the genes reached through ids that Ensembl retired and brought back (rule [`resurrection`](../rules.md#resurrection)); an answer in `broad` but not in `expect` is undetermined;
- **copies apart** (`apart`): genes on haplotypes, patches and alternative loci are expected in the list idtrack reports apart (rule [`judge_copies_apart`](../rules.md#judge_copies_apart)); before release 82, where the GTFs label every gene `primary`, the copies are read from Ensembl's assembly exceptions;
- **the locus rule**: where the history reaches nothing, an answer at the old place is undetermined and one elsewhere wrong; places on haplotypes, patches and the Y pseudoautosomal regions are projected onto the reference through Ensembl's assembly exceptions and the GRC's alignments, and mapped across assemblies by Ensembl's chain files (the corrected judge's `exceptions`, `exceptions_span`; [the corrections of the judge](../harness.md#the-corrections-of-the-judge));
- **genes the GTF lacks**: the core genes of human 48 to 51 and cynomolgus 112 are present (from Ensembl's core gene tables), compared by stable id where the GTF writes no versions, and placed by the core (`core_places`);
- `astar_db`: BioMart's values of the expected genes in the name database, with the corrected judge's `spellings_once`.

### The second reading

For these ids only, the history is read a second, independent way: Ensembl REST `archive/id` (release 116 of the service) for every selected id and every id on its `possible_replacement` chains (`astar/fetch_rest.py`; records in `<truth>/<organism>/rest/archive_astar.json`, POSTs of up to 1,000 ids). Where REST did not answer, the same records are read from the tables REST serves, the release-116 core dumps (`stable_id_event`, `mapping_session`, `gene`; `astar/rest_from_core.py`, `rest/archive_astar_core116.json`, every record marked with its source). A forward line (target after the release read) whose walk and second reading are both complete and give other genes is `open` (undetermined), with REST's genes in `rest` and the record's source in `rest_source`.

The generator's counts of forward lines checked against the second reading (`rest_reading` of each `v1` manifest): pig 5,973,507 agree; mouse 14,822,150 agree, 93 differ, 393,103 without a complete record; human 15,664,742 agree, 5,314 differ, 927,160 without a complete record (2,529 of these ids read again in `v1r1`); rhesus 4,905,087 agree; cynomolgus 395,285 agree. The answers left open because the two readings differ are in `open.tsv` ([Results](#results)).

## Rules

The rules of idtrack the expectations rest on ([rules](../rules.md)): [`start_release`](../rules.md#start_release) (lines without `from_release`), [`own_id_first`](../rules.md#own_id_first), [`ties`](../rules.md#ties), [`all_order`](../rules.md#all_order), [`copies_apart`](../rules.md#copies_apart), [`unknown_version`](../rules.md#unknown_version) and [`version_newer_than_graph`](../rules.md#version_newer_than_graph) (versioned lines), [`external_bridge`](../rules.md#external_bridge) (where the history ends idtrack may bridge). The judge's rules: [`id_lives_on`](../rules.md#id_lives_on), [`resurrection`](../rules.md#resurrection), [`judge_copies_apart`](../rules.md#judge_copies_apart), [`judge_corrections`](../rules.md#judge_corrections) (`sessions_effective`, `exceptions`, `exceptions_span`, `core_places`, `spellings_once`). The decided rules A* measures before they land: [`bridge_same_place`](../rules.md#bridge_same_place) (bridges only at the same place), [`main_line_first`](../rules.md#main_line_first) (the id's own continuation before a branch), [`event_version`](../rules.md#event_version), [`copies_placed_on_primary`](../rules.md#copies_placed_on_primary), [`par_y_ids`](../rules.md#par_y_ids), [`place_moves_marked`](../rules.md#place_moves_marked).

Scoring is F0's ([scoring and outcomes](f0.md#scoring-and-outcomes)), with the rule `open` for the lines the second reading leaves open.

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)); the gene lines and the transcript lines run into two folders, `runs/d4c2dda/<organism>` and `runs/d4c2dda_tx/<organism>`.

```
# selection (one per organism; --name keeps the selection's name, e.g. <organism>.r3)
python -m idtrack_reproducibility.families.astar.selection <o> [--name <o>.r3]
# the sets: genes and transcripts (published once), and revision sets when the truth changes
python -m idtrack_reproducibility.families.astar.gen <o> --set <name> [--selection <o>.r3] [--ids <file>] [--dry-run]
python -m idtrack_reproducibility.families.astar.gen_tx <o> --set <name> [--selection <o>.r3] [--dry-run]
python -m idtrack_reproducibility.families.astar.revise_spellings <o> --set <new> --revises <set>
python -m idtrack_reproducibility.families.astar.revise_rest <o> --set <new> --revises <sets> --ids <file>
python -m idtrack_reproducibility.families.astar.revise_early <o> --set <new> --revises <sets>
python -m idtrack_reproducibility.families.astar.revise_core <o> --set <new> --revises <sets>
python -m idtrack_reproducibility.families.astar.revise_rules <o> [--dry-run]   # a decided rule the truth now follows
# a tree's run (F0's runner and scorer)
python -m idtrack_reproducibility.families.runner --family astar --organism <o> --clone <source> --graph <graph> \
    --out <families>/astar/runs/<tree>/<o> --sets v1,... --submit <n> --score --score-parts <P>
python -m idtrack_reproducibility.families.runner --family astar --organism <o> --clone <source> --graph <graph> \
    --out <families>/astar/runs/<tree>_tx/<o> --sets tx1,... --submit <n> --score
# A*'s own report: per pattern, per pattern and stratum, per id, the open lines
python -m idtrack_reproducibility.families.astar.report_astar <families>/astar/runs/<tree>/<o> --organism <o> [--selection <o>.r3]
# causes: up to 3 explained answers per lineage and stratum, then F0's classifier
<python of the tree> -m idtrack_reproducibility.families.explain --family astar --organism <o> --clone <source> \
    --graph <graph> --scored <families>/astar/runs/<tree>/<o> --per-unit 3
python -m idtrack_reproducibility.families.causes --organism <o> --scored <families>/astar/runs/<tree>/<o>
# the known-limitations rows of every organism
python -m idtrack_reproducibility.families.astar.limitations_astar <tree> [--out <tsv>]
```

The truth's own fetches for A* (beyond the harness's): `python -m idtrack_reproducibility.families.astar.fetch_core_genes` (the core gene tables of human 48 to 51 and cynomolgus 112), `... astar.fetch_rest <o> [--selection <name>]` (REST, at most 2 requests per second, resumable) and `... astar.rest_from_core <o> --ids <file>` (the release-116 core's reading for ids REST did not answer). A* was run on `d4c2dda`; `71690ab` is answer-identical to it on every harness line ([F0's check](f0.md#the-check-against-the-harnesss-scorer)) and was not run on A*.

## Results

Runs `<families>/astar/runs/d4c2dda/<organism>/` (gene lines) and `runs/d4c2dda_tx/<organism>/` (transcripts and translations), on the reference graphs of `d4c2dda` (their content digests in `run.json`). Result files: `results/astar/d4c2dda/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv` and A*'s own tables (`per_pattern.tsv`, `per_pattern_stratum.tsv`, `hard_ids.tsv.gz`, `open.tsv`); `results/astar/d4c2dda_tx/<organism>/report.tsv` and `scores.tsv`. No answer timed out and none raised an error; every identity line (`same`) is right in every organism.

**Every gene line** (pattern `ALL` of `results/astar/d4c2dda/<organism>/per_pattern.tsv`; answers are lines x views; the interval is the cluster bootstrap over the lineage):

| organism | answers | right | wrong | undetermined | rate (95 % interval) | ids with a wrong answer | lineages with a wrong answer |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| pig | 38,719,226 | 36,727,229 | 41,442 | 1,950,555 | 0.998873 (0.998444–0.999224) | 92 of 8,849 | 59 of 2,200 |
| mouse | 108,689,726 | 105,403,927 | 715,746 | 2,570,053 | 0.993255 (0.992292–0.994139) | 534 of 12,184 | 385 of 6,177 |
| human | 108,317,648 | 102,858,510 | 1,356,197 | 4,102,941 | 0.986987 (0.984951–0.988788) | 2,148 of 15,638 | 1,729 of 8,438 |
| rhesus | 30,703,584 | 30,299,794 | 760 | 403,030 | 0.999975 (0.999919–1.000000) | 1 of 4,560 | 1 of 1,134 |
| cynomolgus | 3,385,652 | 3,384,980 | 0 | 672 | 1.000000 (1.000000–1.000000) | 0 of 2,269 | 0 of 441 |

**Per stratum** (row `*` of `results/astar/d4c2dda/<organism>/report.tsv`: rate and cluster-bootstrap interval; lineages behind the wrong answers of the stratum in `units_wrong`):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `same` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `forward` | 0.9974 (0.9964–0.9983) | 0.9987 (0.9979–0.9993) | 0.9976 (0.9966–0.9983) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `backward` | 0.9994 (0.9990–0.9998) | 0.9963 (0.9956–0.9970) | 0.9932 (0.9917–0.9944) | 0.9999 (0.9996–1.0000) | 1.0000 (1.0000–1.0000) |
| `forward_handoff` | 0.9999 (0.9998–1.0000) | 0.9901 (0.9880–0.9921) | 0.9845 (0.9818–0.9870) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `backward_handoff` | 0.9985 (0.9976–0.9992) | 0.9907 (0.9890–0.9922) | 0.9751 (0.9709–0.9789) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `default` | 0.9986 (0.9981–0.9990) | 0.9909 (0.9896–0.9921) | 0.9808 (0.9784–0.9829) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `db` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 0.9998 (0.9994–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |

The weakest strata are human `backward_handoff` (33,208,014 answers, 777,630 wrong, 1,964,227 undetermined, from 1,298 of 8,008 lineages) and `default` (90,950 wrong from 1,360 lineages); mouse `forward_handoff` (218,719 wrong, 186 lineages) and `backward_handoff` (348,849 wrong, 244 lineages); pig `forward` (16,072 wrong, 33 lineages).

**Per pattern** (row `*` of `results/astar/d4c2dda/<organism>/per_pattern.tsv`: rate (cluster bootstrap); ids with a wrong answer of the pattern's ids; – : no id of the pattern):

| pattern | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `branch_while_living` | 1.000000 (1.000000–1.000000); 0 of 684 | – | 1.000000 (1.000000–1.000000); 0 of 35 | – | – |
| `event_contradicts_gtf` | – | 0.998383 (0.997735–0.998947); 35 of 3,181 | 0.907650 (0.894822–0.921303); 225 of 917 | – | – |
| `gap` | – | 0.763476 (0.570128–1.000000); 1 of 2 | 1.000000 (1.000000–1.000000); 0 of 2 | – | – |
| `gtf_lacks` | – | – | 0.999491 (0.998632–1.000000); 2 of 254 | – | 1.000000 (1.000000–1.000000); 0 of 47 |
| `handoff_lineage` | 0.998902 (0.998336–0.999353); 58 of 6,714 | 0.994281 (0.992172–0.996026); 189 of 5,644 | 0.992191 (0.990219–0.993932); 1,052 of 7,300 | 1.000000 (1.000000–1.000000); 0 of 3,761 | 1.000000 (1.000000–1.000000); 0 of 2,048 |
| `handoff_moved` | 0.999651 (0.999256–0.999936); 4 of 575 | 0.978434 (0.966616–0.988593); 20 of 107 | 0.989800 (0.986877–0.992352); 101 of 786 | – | – |
| `hunt_history` | 1.000000 (1.000000–1.000000); 0 of 8 | 0.918049 (0.769809–1.000000); 2 of 6 | 0.967398 (0.938601–0.985621); 16 of 40 | 1.000000 (1.000000–1.000000); 0 of 9 | 1.000000 (1.000000–1.000000); 0 of 7 |
| `lineage_ge5` | 0.998989 (0.998185–0.999529); 44 of 5,580 | 0.994684 (0.992471–0.996532); 162 of 5,407 | 0.992523 (0.986807–0.995873); 371 of 4,651 | 1.000000 (1.000000–1.000000); 0 of 3,304 | 1.000000 (1.000000–1.000000); 0 of 1,656 |
| `many_events` | 0.999149 (0.998278–0.999724); 30 of 4,380 | 0.995286 (0.992751–0.997276); 103 of 4,078 | 0.986296 (0.969862–0.994496); 228 of 2,575 | 1.000000 (1.000000–1.000000); 0 of 2,447 | 1.000000 (1.000000–1.000000); 0 of 1,222 |
| `many_versions` | 0.998544 (0.997323–0.999552); 7 of 294 | 0.997816 (0.996546–0.998919); 20 of 930 | 0.997245 (0.995628–0.998622); 21 of 727 | – | – |
| `merged_then_split` | 1.000000 (1.000000–1.000000); 0 of 1,670 | 0.999946 (0.999791–1.000000); 3 of 913 | 0.980979 (0.968496–0.991114); 57 of 641 | 1.000000 (1.000000–1.000000); 0 of 527 | – |
| `mislabelled_session` | 0.999467 (0.998878–0.999903); 12 of 2,204 | – | 0.998500 (0.994942–1.000000); 1 of 15 | – | – |
| `moved` | 0.998477 (0.997218–0.999527); 7 of 392 | 0.976099 (0.968155–0.982722); 59 of 337 | 0.986967 (0.984429–0.989099); 259 of 1,359 | 0.999926 (0.999776–1.000000); 1 of 801 | 1.000000 (1.000000–1.000000); 0 of 170 |
| `non_reference_copy_ended` | – | 0.899676 (0.881260–0.917205); 150 of 930 | 0.932560 (0.923757–0.941406); 441 of 3,645 | – | – |
| `non_reference_lineage` | – | 1.000000 (1.000000–1.000000); 0 of 6 | 0.985983 (0.974198–0.991933); 238 of 2,412 | – | – |
| `par_y` | – | – | 0.928176 (0.859830–0.980121); 11 of 71 | – | – |
| `rename_chain_ge3` | – | 1.000000 (1.000000–1.000000); 0 of 6 | 0.988884 (0.984274–0.993035); 97 of 453 | – | – |
| `retired_and_back` | – | 0.926620 (0.570128–1.000000); 1 of 3 | 0.999273 (0.997726–1.000000); 1 of 84 | – | – |
| `split_then_merged_back` | 0.998987 (0.998136–0.999628); 27 of 4,164 | 0.998706 (0.997741–0.999462); 12 of 1,906 | 0.977730 (0.964871–0.987071); 138 of 1,270 | 1.000000 (1.000000–1.000000); 0 of 1,443 | – |
| `successor_overlaps` | 0.999312 (0.998571–0.999809); 13 of 1,439 | 0.994084 (0.991506–0.996266); 88 of 2,182 | 0.987419 (0.984864–0.989635); 599 of 4,131 | 1.000000 (1.000000–1.000000); 0 of 1,359 | 1.000000 (1.000000–1.000000); 0 of 772 |
| `travel_weak_cell` | – | – | 0.923390 (0.887379–0.963459); 22 of 65 | – | – |
| `version_fall` | – | 0.952282 (0.570128–1.000000); 1 of 3 | 1.000000 (1.000000–1.000000); 0 of 6 | – | – |
| `version_jump` | – | 0.998427 (0.997804–0.999016); 34 of 3,179 | 0.999496 (0.998460–1.000000); 1 of 16 | – | – |

The hardest patterns are those where Ensembl's history ends: copies on haplotypes and patches whose history ends (`non_reference_copy_ended`), events that contradict the GTF, the Y pseudoautosomal genes, the travel matrix's weak cell. Splits, merges, splits merged back, rename chains, version jumps, branches while the id lives on, mislabelled sessions and large lineages are travelled well in every organism: the lowest of these rates is human `split_then_merged_back`, 0.977730.

**Transcripts and translations** (row `*` of `results/astar/d4c2dda_tx/<organism>/report.tsv`):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `tx_same` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tx_forward` | 0.9951 (0.9929–0.9970) | 0.9954 (0.9931–0.9973) | 0.9847 (0.9814–0.9875) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tx_backward` | 0.9990 (0.9979–0.9997) | 0.9874 (0.9840–0.9901) | 0.9552 (0.9500–0.9599) | 0.9999 (0.9996–1.0000) | 1.0000 (1.0000–1.0000) |
| `tx_forward_handoff` | 1.0000 (1.0000–1.0000) | 0.9782 (0.9731–0.9828) | 0.9618 (0.9566–0.9668) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tx_backward_handoff` | 0.9989 (0.9981–0.9995) | 0.9890 (0.9867–0.9909) | 0.9411 (0.9351–0.9467) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tl_same` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tl_forward` | 0.9939 (0.9910–0.9964) | 0.9895 (0.9843–0.9937) | 0.9758 (0.9705–0.9807) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tl_backward` | 0.9993 (0.9985–0.9999) | 0.9744 (0.9652–0.9816) | 0.9266 (0.9185–0.9341) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tl_forward_handoff` | 1.0000 (1.0000–1.0000) | 0.9565 (0.9451–0.9659) | 0.9387 (0.9302–0.9465) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `tl_backward_handoff` | 0.9983 (0.9968–0.9993) | 0.9741 (0.9682–0.9796) | 0.9089 (0.8953–0.9210) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |

A transcript's answer is its gene's travel, so the transcript strata repeat the genes' weakness where the history ends, at more lines per gene (human `tx_backward_handoff`: 767,242 answers, 38,644 wrong; `tl_backward_handoff`: 290,290 answers, 22,597 wrong). The transcript runs are not classified by cause.

### Causes

The wrong answers of the gene lines are explained (`explain`, up to 3 per lineage and stratum) and classified by F0's classifier ([causes](f0.md#causes)); the other wrong answers of a lineage are counted `unexplained` in `causes_summary.tsv` (they are not asked again with `explain`). Answers summed over the strata of `results/astar/d4c2dda/<organism>/causes_summary.tsv` (cynomolgus has no wrong answer and no table):

| cause | class | pig | mouse | human | rhesus |
| --- | --- | ---: | ---: | ---: | ---: |
| `external_bridge`: Ensembl's history has no successor; the answer is reached through an external id and lies elsewhere | b | 384 | 3,131 | 12,803 | 6 |
| `gtf_lacks_core_gene`: the answer is a core gene the release's GTF lacks, at another place | d | 0 | 0 | 4,335 | 0 |
| `event_version_contradicts_tables`: an event naming a version the id did not have, left out by the builder (one id) | c | 0 | 0 | 704 | 0 |
| `path_through_resurrection`: the path follows an id Ensembl retired and brought back | e | 0 | 3 | 0 | 0 |
| not explained (beyond 3 per lineage and stratum) | – | 41,058 | 712,612 | 1,338,355 | 754 |

Where the history has a successor, idtrack travels it: the explained wrong answers are almost all bridges where Ensembl's history ends (class b; rule [`external_bridge`](../rules.md#external_bridge), which the decided rule [`bridge_same_place`](../rules.md#bridge_same_place) restricts to the same place). The scorer's reason is known for every wrong answer: the limitation `L-astar-1` below counts the wrong answers whose reason is that Ensembl's history has no successor and the answer is not at the query's place. Rhesus has one such id, `ENSMMUG00000036773`, born at release 58 and asked at releases before its birth, answered through its name with genes elsewhere (760 wrong answers; its wrong cells in `hard_ids.tsv.gz`).

### The hard-ids report

`hard_ids.tsv.gz` has one row per selected id: its patterns, its lineage and the lineage's size, its answers, right, wrong, undetermined, timeouts and errors, and up to 40 of its wrong cells (`<kind> <from>-><to> <view>`). `per_pattern_stratum.tsv` gives every pattern per stratum; `open.tsv` the lines the truth leaves open, per reason and stratum. The open lines are all of one reason, the two readings of Ensembl's history differing (rows of `open.tsv`):

| organism | `default` | `forward` | `forward_handoff` |
| --- | ---: | ---: | ---: |
| mouse | 14 | 176 | 168 |
| human | 1,862 | 6,704 | 12,068 |

Pig, rhesus and cynomolgus have no open line.

## Known limitations

A*'s rows of [the known-limitations table](../../results/limitations/limitations.tsv) (generated by `limitations_astar` from the `d4c2dda` runs; answers, with the distinct ids in parentheses):

| limitation | what | user action | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `L-astar-1` | no successor in Ensembl's history (retired without successor, or asked before its birth): the answer is a gene at another place, reached through names or external ids; the decided rule `bridge_same_place` restricts such bridges | check the explained path for external steps; treat answers through bridges as candidates | 41,442 (92) | 715,746 (534) | 1,354,137 (2,143) | 760 (1) | – |
| `L-astar-2` | places across an assembly handoff that no chain maps cannot be compared: no chain file (pig Sscrofa10.2 to 11.1, rhesus MMUL_1 to Mmul_8 to Mmul_10, cynomolgus 5.0 to 6.0) or a sequence outside the chain | read the answer's place yourself (e.g. a liftover with an external alignment) | 1,798,802 (1,321) | 245,006 (201) | 461,924 (581) | 403,030 (185) | 672 (1) |
| `L-astar-3` | copies whose place on the reference no GRC alignment gives (NCBI36 haplotypes, unaligned parts of patches) cannot be judged by place | none needed for reference genes; for copies, compare names | – | 46,448 (59) | 44,646 (172) | – | – |
| `L-astar-4` | Ensembl's two readings of the history (the walk of the release-115 events; REST or the release-116 core) give other genes: the line is open | ask Ensembl REST for the id; prefer the newer release's history | – | 358 (2) | 20,634 (135) | – | – |
| `L-astar-5` | the answer is a branch that left the id while it lived on (within the broad walk), not its own continuation; the decided rule `main_line_first` ranks the continuation first | use strategy `all` to see every branch | 99,717 (229) | 44,073 (175) | 876,312 (1,350) | – | – |

`L-astar-1` is wrong; `L-astar-2` to `L-astar-5` are undetermined answers. Each row also carries its share of the organism's answers, the run and the graph's digest.

## Reproducing from scratch

From an empty data root, with the configuration naming the trees and their graphs:

1. **Fetch the truth** as the harness does ([F0's steps](f0.md#reproducing-from-scratch): GTFs, core history tables, BioMart tables, chains, assembly exceptions, the GRC's alignments, core places), then A*'s own: `python -m idtrack_reproducibility.families.astar.fetch_core_genes`; after the selection, `... astar.fetch_rest <o> --selection <o>.r3` and, for ids REST does not answer, `... astar.rest_from_core <o> --ids <file>`.
2. **Select**: `python -m idtrack_reproducibility.families.astar.selection <o> --name <o>.r3` (the current code's patterns are those of `r3`; the fixed ids come with the package).
3. **Publish the sets**: `... astar.gen <o> --set v1 --selection <o>.r3` and `... astar.gen_tx <o> --set tx1 --selection <o>.r3`. The published sets were made in steps (the first selection, then `hunt1` for the ids `r3` added, and the revision sets of the table above); the current generator writes the truth of all their corrections, so a rebuild is compared with the published sets' latest lines (`lines.load`) and any difference is published as a revision set, never by rewriting a set.
4. **Run and score** on `d4c2dda`: the two runner commands of [Running it](#running-it) (`runs/d4c2dda/<o>` with the gene sets, `runs/d4c2dda_tx/<o>` with the transcript sets).
5. **Report**: `... astar.report_astar <families>/astar/runs/d4c2dda/<o> --organism <o> --selection <name>`, the selection the set `v1` names in its manifest (`sources`).
6. **Classify**: `explain` (in the tree's Python, `--per-unit 3`) and `causes` on the gene run.
7. **Limitations**: `idtrack-repro known-limitations` (A*'s rows: `... astar.limitations_astar d4c2dda --out <tsv>`).
8. **Results**: `idtrack-repro results` writes `results/astar/d4c2dda/<organism>/` and `results/astar/d4c2dda_tx/<organism>/`, the known-limitations table, the index and the provenance.
