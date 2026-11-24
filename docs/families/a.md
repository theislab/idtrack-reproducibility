# A: history travel, exhaustive

Family A asks idtrack to carry gene, transcript and translation ids through Ensembl's history in every way the history allows, exhaustively where the size permits: every id at every release it exists at to the snapshot, without a release asked from, to its neighbouring releases, across every event and every assembly handoff, the full matrix of releases, and the graph of an older snapshot. Its first stratum, A0, reads the graph itself against Ensembl's tables, without any search. Its code is `idtrack_reproducibility.families.a`; its data are `<families>/a/` (`a0/<tree>/<organism>/`, `queries/<organism>/`, `runs/<tree>/<organism>/`, `patterns/`); the design is in [the methods](../methods.md#a-history-travel). It uses A*'s truth ([A*](astar.md#truth)) and the common machinery ([F0](f0.md)).

## What it measures

- **A0, graph fidelity**: whether the graph holds what Ensembl's tables hold, release by release: every gene is present where the GTF has it and nowhere else, every history event is an edge, every transcript and translation is linked to its gene, every name and cross-reference of BioMart is linked at its release ([A0](#a0-graph-fidelity)).
- **A1 to A8, history travel**: the rate of right answers (right / (right + wrong + error); undetermined and timeouts apart) per stratum, organism, era of origin (the assembly of the release asked from), direction, distance in releases and assembly crossing (the travel matrix's axes), with clustered intervals over the id's lineage; and per pattern of an id's history (the weakness hunts' findings as strata over the whole sets).

Every line is scored under `best` and under `all` (`strategy: both`), except A7's lines to the release asked from (`best`).

## A0: graph fidelity

`idtrack_reproducibility.families.a.a0_fidelity` reads the graph through the accessors both api generations have (history edges, cross-reference edges, active releases, the store's names and types) and compares it with the truth store, over every release of the graph's window that the truth has:

- **P presence**: every GTF gene (id.version; the id where the GTF writes no versions) of every release is a backbone node active at that release, read two ways (the ranges of the history edges; the releases at which the base id and the transcripts link the gene), and no other node is.
- **H history**: every gene event of the release-115 `stable_id_event` in the window (placed where its effect appears) is the edge it should be: between ids, old.v -> new.v at the session's two releases; a version change; a birth from `<id>.Void` or `.Retired`; a retirement to `<id>.Retired`; an id mapped onto itself with the same version only by the builder's rule (score below 1, or the session maps it to or from another id). And every history edge is an event, a version change, a birth, a retirement or a return that the GTFs show.
- **T transcripts**: every GTF transcript links its gene at its release, every translation its transcript, and the converse.
- **N names**: every value of the BioMart tables of a release (display name with its source, synonym, NCBI id and label, HGNC, MGI and VGNC symbols, UniProt gene name, Swiss-Prot, TrEMBL, RefSeq, RFAM, WikiGene) is linked at that release to the gene (directly or through a transcript or translation) under the graph's database of that table, and the converse for the tables the graph includes.

Every class of rows is counted with its distinct ids and an example; a class explained by a rule of idtrack, the builder's documented rule, the externals file or the truth itself is marked `by_design` with its reason; every other class that is not `ok` is a finding. Only the assembly that serves a release is read (human: the links of GRCh37, Ensembl's frozen parallel assembly, labelled with releases 76 to 115 are counted apart; rule [`serving_assembly`](../rules.md#serving_assembly)). VGNC symbols come from VGNC's own gene sets (BioMart gives only the VGNC id), so their rows measure the authority, never the graph.

## The sets

One set per stratum and organism, `<stratum>_v1`, published once (F0's `lines.publish`); qid `<kind>|<query>|<from or ->|<to>`; unit the id's lineage (the smallest id of its component of cross-id events; else the id); extra fields `era` (the assembly of the release asked from), `dist` (target minus origin) and `cross` (whether the assembly differs).

| stratum | what is asked | truth |
| --- | --- | --- |
| A1 | every gene id at every release R where the GTF has it, from R to the target (115); the stable id and, where the GTF of R writes versions, the id with its version | the strict walk with living ids; broad; locus; copies apart |
| A2 | every gene id and every version it had, without `from_release`, to the target (rule [`start_release`](../rules.md#start_release)) | the walk from the release the rule reads the id at |
| A3 | every (id, R) to R-2, R-1, R+1, R+2 and to the first release of the window (stable id); the versioned id to R-1 and R+1 | walks |
| A4 | every cross-id gene event of `stable_id_event` in the window: the old id from the session's old release to its new release and the new id back, with the event's versions; every version change of one id, the old version forward | the walk of that session and the strict walk beyond it |
| A5 | every assembly handoff (last release L of one assembly, first F of the next): every id present at L to F and every id present at F to L; strata `A5_<L>_<F>_lives`, `_dies`, `_born` | walks; places through the chain files |
| A6 | the full matrix: every release of the id to every release of the window (stable id); every id of cynomolgus, and for the other organisms a seeded sample of 1,000 ids per history class (unchanged, born late, renamed, merged, split, retired without successor, retired to a branch) and assembly era of the id's first release | as the travel matrix |
| A7 | every transcript and translation id of the GTFs from its first and its last release R to its gene at R (`_at`, `best`) and to the target (`_115`); a seeded sample of (id, release) pairs of the whole matrix to the target (`_sample`) | the GTF's gene at R, then the gene's walk |
| A8 | A1 against the graph of an older snapshot, 110 (pig and human): every (id, R <= 110) to 110 | walks to 110 |

Lines per set and organism (manifests `<families>/a/queries/<organism>/<set>.manifest.json`):

| set | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `a1_v1` | 2,765,612 | 5,290,751 | 6,782,338 | 3,456,451 | 1,704,612 |
| `a2_v1` | 172,533 | 316,156 | 323,050 | 157,598 | 115,893 |
| `a3_v1` | 10,131,601 | 19,573,433 | 25,434,484 | 13,078,525 | 5,609,602 |
| `a4_v1` | 67,248 | 148,683 | 119,192 | 48,821 | 23,951 |
| `a5_v1` | 147,873 | 295,316 | 208,154 | 260,826 | 135,748 |
| `a6_v1` | 15,263,820 | 18,554,072 | 19,219,248 | 18,339,668 | 21,307,650 |
| `a7_v1` | 957,166 | 1,779,923 | 3,378,808 | 918,724 | 796,345 |
| `a8_s110_v1` | 2,408,840 | – | 5,981,424 | – | – |

Revision sets (only the lines whose truth changed, under their published qids; no published qid was lost):

| set | mouse | human | why a revision |
| --- | ---: | ---: | --- |
| `a1_r1` | 28,656 | 31,938 | the shared truth changed after publication: places on haplotypes, patches, alternative loci and the Y pseudoautosomal regions projected onto the reference through Ensembl's assembly exceptions and the GRC's alignments (the corrected judge's `exceptions` and `exceptions_span`); and an A4 line whose event names a version the id never had in any GTF expects nothing (rules [`unknown_version`](../rules.md#unknown_version), [`version_newer_than_graph`](../rules.md#version_newer_than_graph)) |
| `a2_r1` | 2,397 | 4,879 | as `a1_r1` |
| `a3_r1` | 21,919 | 315,871 | as `a1_r1` |
| `a4_r1` | 173 | 918 | as `a1_r1` |
| `a5_r1` | 1,634 | 4,174 | as `a1_r1` |
| `a6_r1` | 214,845 | 573,015 | as `a1_r1` |
| `a7_r1` | 7,631 | 27,643 | as `a1_r1` |
| `a8_s110_r1` | – | 31,219 | as `a1_r1` |
| `a4_r2` | – | 917 | A4 lines whose id is in no gene table at the release asked from (Ensembl's version-0 placeholder ids of human assembly patches) are open |
| `a3_r2` | – | 19,558 | the core genes the GTF lacks (human 48 to 51) are expected without the core's version at releases whose GTF writes none, as the comparison there is by stable id |
| `a6_r2` | – | 4,148 | as `a3_r2` |
| `a3_r3` | – | 1,258 | as `a3_r2`, for the published lines whose truth the generator writes differently since (`fix_versionless`) |

Pig, rhesus and cynomolgus have no revision set of these: the place correction concerns human and mouse, and regenerating their A4 sets changed no line.

Rule [`no_mito_nuclear_step`](../rules.md#no_mito_nuclear_step) (decided, built in idtrack's next version after its first round; judge `rules`): the shared truth (`astar.truth.TruthA`) cuts every step between the mitochondrial genome and a nuclear sequence, also where the id stays the same, and an id of such a step has no direct walk (`GenA.walks`). `revise_a --rule no_mito_nuclear_step` makes again only the lines of the lineages the rule changes (the units of `Truth.mito_steps`) and publishes those whose truth differs (pig and mouse have no such id):

| set | rhesus | cynomolgus | human |
| --- | ---: | ---: | ---: |
| A1 (`a1_r1`; human `a1_r2`) | 355 | 120 | 21 |
| A2 (`a2_r1`) | 8 | 5 | – |
| A3 (`a3_r1`; human `a3_r4`) | 130 | 105 | 221 |
| A4 (`a4_r1`) | 5 | 5 | – |
| A5 (`a5_r1`) | 20 | 20 | 6 |
| A6 (`a6_r1`; human `a6_r3`) | – | 1,560 | – |
| A7 (`a7_r1`; human `a7_r2`) | 22 | 19 | 18 |
| A8 (human `a8_s110_r2`) | – | – | 21 |

## Truth

Nothing of idtrack. A imports A*'s truth and never copies it (`astar.truth.TruthA` read by `astar.gen.Gen.truth_at`; [A*'s truth](astar.md#truth)), so a correction of that truth reaches A as a revision set: the strict walk with living ids (`expect`; rule [`id_lives_on`](../rules.md#id_lives_on)), sessions placed where their effect appears; the broad walk with returned ids (`broad`; rule [`resurrection`](../rules.md#resurrection)); copies on haplotypes, patches and alternative loci reported apart (rule [`judge_copies_apart`](../rules.md#judge_copies_apart)); the locus rule, with places projected through Ensembl's assembly exceptions and mapped across assemblies by the chain files; the core genes the GTF lacks (human 48 to 51, cynomolgus 112). Ensembl REST is not read: the second reading is A*'s, for A*'s ids. For A7, the gene the GTF of R gives the transcript (a translation: the gene of its transcript), then that gene's walk.

A0's truth is the truth store's tables: the GTFs of every release, the release-115 `stable_id_event` and `mapping_session`, BioMart's tables of each release and VGNC's gene sets.

The patterns of every id (`<families>/a/patterns/<organism>.tsv.gz`, `idtrack_reproducibility.families.a.patterns_a`) come from Ensembl's tables only: `cross_events` (a lineage), `lineage5` (a lineage of 5 ids or more), `retired_none` (retired and the strict walk reaches nothing; hunt finding hh-4), `split`, `merged`, `renamed`, `retired_branch` (the travel matrix's class of a retired id), `gap` and `version_drop` (hh-3), `version_jump`, `moved_seq` (another chromosome at another release of one assembly; hh-2), `strand_flip`, `patch` (R6, J3, J4), `mt` (R1), `core_only` (hh-5, J2), `handoff_dies`, `handoff_born` (present on one side of a handoff only), `astar` (selected by A*); `plain` when none applies.

## Rules

The rules of idtrack the expectations rest on are A*'s ([A*'s rules](astar.md#rules)): [`start_release`](../rules.md#start_release) (A2), [`own_id_first`](../rules.md#own_id_first), [`ties`](../rules.md#ties), [`copies_apart`](../rules.md#copies_apart), [`lrg_copies`](../rules.md#lrg_copies), [`unknown_version`](../rules.md#unknown_version), [`version_newer_than_graph`](../rules.md#version_newer_than_graph), [`external_bridge`](../rules.md#external_bridge); the judge's [`id_lives_on`](../rules.md#id_lives_on), [`resurrection`](../rules.md#resurrection), [`judge_copies_apart`](../rules.md#judge_copies_apart), [`judge_corrections`](../rules.md#judge_corrections). The decided rules A measures before they land: [`bridge_same_place`](../rules.md#bridge_same_place), [`main_line_first`](../rules.md#main_line_first), [`event_version`](../rules.md#event_version), [`two_lives`](../rules.md#two_lives), [`lrg_answers_gene`](../rules.md#lrg_answers_gene). A0's classes by design cite [`spelling_steps`](../rules.md#spelling_steps) (a node is named as the value was written at its first release), [`mirbase_genes`](../rules.md#mirbase_genes), [`rfam_display_label`](../rules.md#rfam_display_label), [`event_version`](../rules.md#event_version) and [`serving_assembly`](../rules.md#serving_assembly).

Scoring is F0's ([scoring and outcomes](f0.md#scoring-and-outcomes)).

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)). A tree variant built from an older snapshot is a tree of its own in the configuration (`[trees.d4c2dda_s110]`: the same checkout, the graphs of snapshot 110).

```
# A0 (no search; the tree's Python)
<python of the tree> -m idtrack_reproducibility.families.a.a0_fidelity <o> --clone <source> --graph <graph> \
    --out <families>/a/a0/<tree>/<o> [--checks P,H,T,N]
<python of the tree> -m idtrack_reproducibility.families.a.probe_a0 <o> --clone <source> --graph <graph> --work <folder> <node or id> ...
# the sets (published once) and the patterns
python -m idtrack_reproducibility.families.a.gen_a <o> --stratum A1|A2|A3|A4|A5|A6|A8 --set <name> [--snapshot 110] [--dry-run]
python -m idtrack_reproducibility.families.a.gen_a7 <o> --set a7_v1 [--sample 200000] [--dry-run]
python -m idtrack_reproducibility.families.a.patterns_a <o>
# revision sets when the shared truth changes (--rule: only the lineages a decided rule changes, judge `rules`)
python -m idtrack_reproducibility.families.a.revise_a <o> --stratum <S> --set <published set>[,<its revisions>] --name <new> [--rule no_mito_nuclear_step] [--dry-run]
python -m idtrack_reproducibility.families.a.fix_versionless <o> --sets <set>,<its revisions> --name <new> [--dry-run]
# a tree's run (F0's runner and scorer), then A's report along the axes and patterns
python -m idtrack_reproducibility.families.runner --family a --organism <o> --clone <source> --graph <graph> \
    --out <families>/a/runs/<tree>/<o> --submit <n> --score --score-parts <P> --time 04:00:00 [--tag a_<tree>_<o>]
python -m idtrack_reproducibility.families.a.report_a <families>/a/runs/<tree>/<o> --organism <o>
# causes: A's explain reads only the lines it explains; F0's classifier with A's plugin
<python of the tree> -m idtrack_reproducibility.families.a.explain_a --organism <o> --clone <source> --graph <graph> \
    --scored <families>/a/runs/<tree>/<o>
python -m idtrack_reproducibility.families.causes --organism <o> --scored <families>/a/runs/<tree>/<o> \
    --plugin src/idtrack_reproducibility/families/a/causes_a.py
# the known-limitations rows (wrong answers by cause, undetermined answers by the judge's reason)
python -m idtrack_reproducibility.families.a.limits_a --out <tsv> <families>/a/runs/<tree>/<o> [...]
```

The snapshot-110 graphs are built by the tree from the same downloaded tables, the snapshot's own history table fetched by idtrack's builder: `idtrack-repro family a build_a8 --organism <o> --tree d4c2dda [--submit]` (`families/a/build_a8.py`: the reference graph's tables copied into `<families>/a/a8/<tree>/<organism>/`, then the tree's public api builds snapshot 110 there), and A8 is run on them into `runs/d4c2dda_s110/<organism>` (`--sets a8_s110_v1,...`). A later tree is run with its own `--tag`, so that its job names do not repeat an earlier run's. A tree whose builder reads tables the reference graphs' folder lacks (from `r2` on: Ensembl's assembly exceptions, chain files, coordinate systems and sequence regions) builds its snapshot-110 graphs as any graph of a tree ([a tree's graphs](../graphs.md)), from the stored tables its graphs of 115 were built from: `idtrack-repro graph-build --tree main_s110 --organisms sus_scrofa,homo_sapiens --snapshot 110 --tables <its tables folder with {organism}> --submit` (`[trees.main_s110]`: main's checkout and Python, the graphs `<graphs>/main_s110/<organism>/`); for `main` the builder downloaded only the snapshot's own two history tables, `ens110_mysql_stable_id_event.h5` and `ens110_mysql_mapping_session.h5`, which lie in `<graphs>/main_s110/<organism>/build/idtrack_tables/<organism>_assembly-<key>/` (pig `sus_scrofa_assembly-111`, human `homo_sapiens_assembly-38`) and are kept there; every other table is a link to the stored tables. These graphs are not judged: `graph-judge`'s history check reads Ensembl's history from the core tables of the snapshot release in the truth store, which holds those of release 115 and not of 110 (fetched for nothing else), so A8 on `main_s110` rests on the build's record (`build.json`: digest, wall time, peak memory, no linked table changed) and on A's truth, not on a judgement of the graph (row `L-a8-graphs-not-judged` of the known limitations).

## Results

Runs `<families>/a/runs/d4c2dda/<organism>/` (the reference graphs, every set with its revision sets, judge `revised`) and `runs/d4c2dda_s110/{sus_scrofa,homo_sapiens}/` (A8 on the snapshot-110 graphs); A0 `<families>/a/a0/d4c2dda/<organism>/summary.tsv`. Result files: `results/a/d4c2dda/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv`, `axes.tsv`, `patterns.tsv`, `wrong_units.tsv` and `a0_summary.tsv`; `results/a/d4c2dda_s110/<organism>/` for A8. No answer timed out and none raised an error. Human's run on the 115 graph also holds the A8 lines (to 110, asked on the 115 graph). A was run on `d4c2dda`; `71690ab` is answer-identical to it on every harness line ([F0's check](f0.md#the-check-against-the-harnesss-scorer)) and was not run on A.

### A0

Rows `ok` / `by_design` / not explained per check, summed over the classes of `<families>/a/a0/d4c2dda/<organism>/summary.tsv` (A0's run of 2026-11-18 by the same checks of `a0_fidelity`; its reasons by design written by this package from the run's classes, `a0_fidelity <organism> --out <folder> --resummarize`; result file `results/a/d4c2dda/<organism>/a0_summary.tsv`) (in parentheses the rows counted apart: history events of sessions outside the graph's window and ids mapped onto themselves with score 1, which the builder by its rule makes no edge; human's links on GRCh37 at releases it does not serve):

| check | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| P presence | 1,621,150 / 0 / 0 | 3,157,745 / 1,859 / 0 | 4,138,229 / 138,631 / 0 | 2,157,903 / 0 / 0 | 852,306 / 0 / 0 |
| H history | 376,398 / 0 / 0 | 920,770 / 3,182 / 0 (358,475) | 1,039,514 / 3,664 / 1 (192,990) | 283,047 / 0 / 0 (89,914) | 175,974 / 1 / 0 |
| T transcripts | 4,471,790 / 0 / 0 | 11,807,276 / 0 / 0 | 21,160,774 / 96,074 / 78 | 6,322,758 / 0 / 0 | 2,719,756 / 218 / 0 |
| N names | 2,377,827 / 255,415 / 3 | 1,911,720 / 183,678 / 10 | 2,707,676 / 251,189 / 0 (1,674,399) | 1,153,800 / 173,379 / 0 | 489,353 / 55,854 / 0 |

By design (each class with its reason in A0's `summary.tsv`): LRG records in the graph that the GTF and BioMart lack (human); genes on patches at releases whose GTF has no patch file, so that every gene of the GTF reads primary (human, mouse); Ensembl events naming a version the gene tables do not have, which the builder follows the tables for (rule [`event_version`](../rules.md#event_version)'s kind), and births recorded later than the tables show; events placed at the nearest pair of releases by the builder's rule; version 0; events of ids in no table; core genes the GTF lacks and their transcripts; miRBase display names (rule [`mirbase_genes`](../rules.md#mirbase_genes)); WikiGene (not in the organism's externals); RFAM labels held as a second id (rule [`rfam_display_label`](../rules.md#rfam_display_label)); node spellings of the first release; BioMart NCBI tables without a label column; VGNC symbols that differ from Ensembl's labels.

Not explained (the classes of A0's `summary.tsv` that are neither `ok` nor by design; their rows of the known-limitations table, `limits_a --a0`: `L-id-shared-across-types` for human's transcript id, `L-event-versions` for human's event, `L-a0-unexplained` for the rest, read one by one below):

| organism | check | class | rows | what |
| --- | --- | --- | ---: | --- |
| human | history | `between_ids:no_edge` | 1 | `ENSG00000233776.3` -> `ENSG00000173575.24` (107 -> 108) is no edge: the event's old version is not the id's version at 107 (rule [`event_version`](../rules.md#event_version), decided) |
| human | transcripts | `tl:graph_extra`, `tx:not_linked_at_release` | 39 + 39 | Ensembl 107's GTF gives transcript `ENST00000703410` the protein id `ENST00000515292`, the stable id of a transcript of `ENSG00000251182` (releases 58 to 96); the graph keeps one node per name and types it as a translation, so the transcript's links at 58 to 96 are held by a translation node |
| pig | names | `display:missing`, `uniprot_gn:missing` | 2 + 1 | `MYD88`, `NPY` (BioMart 80 display names from UniProtKB Gene Name) and `PIGBRCA1` (a UniProt gene name at 77) are not linked at those releases, while the graph links them at neighbouring releases |
| mouse | names | `uniprot_gn:missing`, `uniprot_gn:under_another_database` | 8 + 2 | BioMart 75 writes Ensembl gene ids as UniProt gene names (the truth's rows) |

### A1 to A7

Rate and cluster-bootstrap interval per stratum (row `*` of `results/a/d4c2dda/<organism>/report.tsv`; – : no such stratum):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `A1` | 0.9971 (0.9968–0.9974) | 0.9966 (0.9963–0.9969) | 0.9975 (0.9973–0.9977) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A2` | 0.9949 (0.9944–0.9954) | 0.9964 (0.9961–0.9967) | 0.9955 (0.9952–0.9958) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A3_+1` | 0.9997 (0.9996–0.9997) | 0.9998 (0.9998–0.9999) | 0.9999 (0.9998–0.9999) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A3_+2` | 0.9992 (0.9991–0.9992) | 0.9996 (0.9996–0.9996) | 0.9996 (0.9996–0.9996) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A3_-1` | 0.9996 (0.9996–0.9997) | 0.9997 (0.9997–0.9997) | 0.9995 (0.9995–0.9995) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A3_-2` | 0.9989 (0.9988–0.9990) | 0.9993 (0.9992–0.9993) | 0.9986 (0.9985–0.9986) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A3_first` | 0.9863 (0.9853–0.9872) | 0.9838 (0.9828–0.9848) | 0.9661 (0.9648–0.9673) | 0.9995 (0.9993–0.9996) | 1.0000 (1.0000–1.0000) |
| `A4_fwd` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 0.9996 (0.9988–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A4_back` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A4_fwd_version_not_in_gtf` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 0.9947 (0.9923–0.9967) | – | – |
| `A4_back_version_not_in_gtf` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 0.9716 (0.9656–0.9767) | – | – |
| `A4_ver` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A6_same` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A6_forward` | 0.9980 (0.9974–0.9985) | 0.9962 (0.9949–0.9972) | 0.9978 (0.9971–0.9984) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A6_backward` | 0.9994 (0.9990–0.9996) | 0.9951 (0.9941–0.9960) | 0.9905 (0.9894–0.9915) | 0.9989 (0.9983–0.9994) | 1.0000 (1.0000–1.0000) |
| `A6_forward_handoff` | 0.9972 (0.9964–0.9980) | 0.9892 (0.9868–0.9914) | 0.9924 (0.9904–0.9941) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A6_backward_handoff` | 0.9896 (0.9879–0.9912) | 0.9859 (0.9834–0.9881) | 0.9774 (0.9746–0.9800) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tx_at` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tx_115` | 0.9966 (0.9961–0.9970) | 0.9970 (0.9967–0.9973) | 0.9965 (0.9961–0.9968) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tx_sample` | 0.9978 (0.9974–0.9982) | 0.9976 (0.9972–0.9979) | 0.9980 (0.9977–0.9983) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tl_at` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tl_115` | 0.9958 (0.9952–0.9964) | 0.9926 (0.9918–0.9933) | 0.9947 (0.9940–0.9953) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `A7_tl_sample` | 0.9974 (0.9969–0.9979) | 0.9969 (0.9964–0.9974) | 0.9978 (0.9974–0.9982) | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |

Answers, wrong and undetermined of the main strata (same table):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `A1` | 5,531,224 / 14,546 / 524,488 | 10,581,502 / 35,231 / 138,575 | 13,564,676 / 33,773 / 189,201 | 6,912,902 / 0 / 516,652 | 3,409,224 / 0 / 69,636 |
| `A2` | 345,066 / 1,612 / 28,074 | 632,312 / 2,214 / 10,520 | 646,100 / 2,791 / 20,821 | 315,196 / 0 / 22,994 | 231,786 / 0 / 5,804 |
| `A3_first` | 3,121,472 / 37,998 / 346,870 | 6,146,942 / 97,290 / 129,224 | 8,073,422 / 262,551 / 335,816 | 4,149,266 / 2,046 / 160,556 | 1,530,132 / 0 / 31,746 |
| `A6_backward_handoff` | 9,134,152 / 86,598 / 831,666 | 11,661,290 / 159,909 / 348,178 | 10,865,400 / 234,430 / 474,143 | 11,706,112 / 0 / 1,065,100 | 12,027,600 / 0 / 376,376 |
| `A7_tx_115` | 429,936 / 1,354 / 33,461 | 1,229,086 / 3,586 / 18,254 | 2,222,452 / 7,691 / 42,982 | 413,318 / 0 / 26,096 | 339,088 / 0 / 7,578 |

**A5**, every handoff (same table; the strata are named by the handoff's two releases):

| organism | stratum | answers | wrong | undetermined | rate (95 % interval) |
| --- | --- | ---: | ---: | ---: | ---: |
| pig | `A5_66_67_born` | 20,134 | 1,426 | 4,686 | 0.9077 (0.9005–0.9142) |
| pig | `A5_66_67_dies` | 11,036 | 767 | 2,375 | 0.9114 (0.9018–0.9203) |
| pig | `A5_66_67_lives` | 59,768 | 0 | 210 | 1.0000 (1.0000–1.0000) |
| pig | `A5_89_90_born` | 39,676 | 0 | 8,764 | 1.0000 (1.0000–1.0000) |
| pig | `A5_89_90_dies` | 37,444 | 0 | 17,336 | 1.0000 (1.0000–1.0000) |
| pig | `A5_89_90_lives` | 127,688 | 0 | 0 | 1.0000 (1.0000–1.0000) |
| mouse | `A5_67_68_born` | 6,254 | 419 | 865 | 0.9222 (0.9052–0.9369) |
| mouse | `A5_67_68_dies` | 7,606 | 456 | 734 | 0.9336 (0.9212–0.9449) |
| mouse | `A5_67_68_lives` | 136,752 | 0 | 0 | 1.0000 (1.0000–1.0000) |
| mouse | `A5_102_103_born` | 5,180 | 16 | 48 | 0.9969 (0.9929–0.9993) |
| mouse | `A5_102_103_dies` | 15,600 | 360 | 2,876 | 0.9717 (0.9638–0.9782) |
| mouse | `A5_102_103_lives` | 419,240 | 0 | 0 | 1.0000 (1.0000–1.0000) |
| human | `A5_54_55_born` | 31,708 | 588 | 1,181 | 0.9807 (0.9783–0.9829) |
| human | `A5_54_55_dies` | 18,008 | 138 | 1,474 | 0.9917 (0.9896–0.9935) |
| human | `A5_54_55_lives` | 113,724 | 0 | 695 | 1.0000 (1.0000–1.0000) |
| human | `A5_75_76_born` | 10,730 | 651 | 1,791 | 0.9272 (0.9198–0.9344) |
| human | `A5_75_76_dies` | 12,570 | 371 | 1,743 | 0.9657 (0.9607–0.9704) |
| human | `A5_75_76_lives` | 229,568 | 0 | 63 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_85_86_born` | 47,256 | 0 | 5,568 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_85_86_dies` | 38,696 | 0 | 7,664 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_85_86_lives` | 164,576 | 0 | 0 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_97_98_born` | 64,352 | 0 | 5,924 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_97_98_dies` | 52,316 | 0 | 7,700 | 1.0000 (1.0000–1.0000) |
| rhesus | `A5_97_98_lives` | 154,456 | 0 | 0 | 1.0000 (1.0000–1.0000) |
| cynomolgus | `A5_102_103_born` | 74,780 | 0 | 4,796 | 1.0000 (1.0000–1.0000) |
| cynomolgus | `A5_102_103_dies` | 37,876 | 0 | 5,344 | 1.0000 (1.0000–1.0000) |
| cynomolgus | `A5_102_103_lives` | 158,840 | 0 | 0 | 1.0000 (1.0000–1.0000) |

Every id that lives across a handoff is carried across it (`_lives`: no wrong answer in any organism); the wrong answers are ids that die or are born at a handoff, whose history has no successor or predecessor there.

### A8, an older snapshot

A1 against the graph of snapshot 110 (row `A8` of `results/a/d4c2dda_s110/<organism>/report.tsv`; human's A8 lines were also asked on the 115 graph, row `A8` of `results/a/d4c2dda/homo_sapiens/report.tsv`):

| graph | organism | answers | wrong | undetermined | rate (95 % interval) | lineages with a wrong answer |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| snapshot 110 | pig | 4,817,680 | 14,314 | 523,400 | 0.9967 (0.9963–0.9970) | 398 of 56,799 |
| snapshot 110 | human | 11,962,848 | 31,350 | 185,662 | 0.9973 (0.9971–0.9976) | 1,546 of 100,607 |
| snapshot 115 | human | 11,962,848 | 31,215 | 185,797 | 0.9973 (0.9971–0.9976) | 1,546 of 100,607 |

Human's A8 lines give the same rate on the graph of snapshot 110 and on the graph of 115 (31,350 and 31,215 wrong answers): a graph of an older snapshot answers as the newest graph does.

### Axes and patterns

Along the travel matrix's axes (`results/a/d4c2dda/<organism>/axes.tsv`; rate and cluster-bootstrap interval), within one assembly the rates are near 1 and the wrong answers concentrate where an assembly is crossed:

| organism | group | crossing an assembly | within one assembly |
| --- | --- | ---: | ---: |
| pig | A1 | 1.00000 (1.00000–1.00000) | 0.99549 (0.99500–0.99596) |
| pig | A6 | 0.99351 (0.99250–0.99441) | 0.99872 (0.99839–0.99905) |
| mouse | A1 | 0.99515 (0.99468–0.99562) | 0.99998 (0.99995–0.99999) |
| mouse | A6 | 0.98745 (0.98569–0.98929) | 0.99577 (0.99499–0.99648) |
| human | A1 | 0.99120 (0.99051–0.99188) | 0.99909 (0.99893–0.99924) |
| human | A6 | 0.98452 (0.98274–0.98621) | 0.99431 (0.99366–0.99492) |
| rhesus | A1 | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| rhesus | A6 | 1.00000 (1.00000–1.00000) | 0.99944 (0.99915–0.99970) |
| cynomolgus | A1 | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| cynomolgus | A6 | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |

Pig's A1 within one assembly is the current assembly, Sscrofa11.1 (era `Sscrofa11.1` 0.99549 in `axes.tsv`): its wrong answers are ids retired there without successor. By era of origin, A1 is weakest from the oldest assemblies: mouse NCBIM37 0.98789 (0.98671–0.98902), human NCBI36 0.98288 (0.98156–0.98424). The same file gives every axis per group (`era`, `direction`, `distance`, `crossing`, `origin`).

Per pattern of an id's history (`results/a/d4c2dda/<organism>/patterns.tsv`, groups A1 and A6; rate and cluster-bootstrap interval; – : no id of the pattern):

| group | pattern | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| A1 | `plain` | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A1 | `lineage5` | 0.99763 (0.99605–0.99885) | 0.99480 (0.99104–0.99741) | 0.99587 (0.99156–0.99824) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A1 | `astar` | 0.99761 (0.99660–0.99847) | 0.99662 (0.99543–0.99767) | 0.99586 (0.99452–0.99693) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A1 | `retired_none` | 0.98034 (0.97826–0.98238) | 0.92277 (0.91573–0.92991) | 0.93620 (0.93083–0.94133) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A1 | `handoff_dies` | 1.00000 (1.00000–1.00000) | 0.95927 (0.95087–0.96647) | 0.94348 (0.93772–0.94869) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A1 | `patch` | – | 0.61925 (0.54495–0.68868) | 0.99852 (0.99792–0.99902) | – | – |
| A1 | `moved_seq` | 1.00000 (1.00000–1.00000) | 0.97484 (0.96274–0.98541) | 0.97928 (0.97047–0.98610) | – | – |
| A1 | `strand_flip` | 0.99836 (0.99480–1.00000) | 0.98028 (0.95923–0.99545) | 0.98830 (0.98148–0.99389) | – | – |
| A6 | `plain` | 0.99958 (0.99933–0.99978) | 0.99558 (0.99443–0.99659) | 0.99629 (0.99547–0.99713) | 0.99963 (0.99941–0.99983) | 1.00000 (1.00000–1.00000) |
| A6 | `lineage5` | 0.99932 (0.99829–0.99990) | 0.99526 (0.99239–0.99752) | 0.98464 (0.97809–0.98953) | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) |
| A6 | `astar` | 0.99915 (0.99838–0.99969) | 0.99494 (0.99305–0.99650) | 0.97959 (0.97664–0.98239) | 0.99993 (0.99977–1.00000) | 1.00000 (1.00000–1.00000) |
| A6 | `retired_none` | 0.98520 (0.98298–0.98736) | 0.95687 (0.94976–0.96392) | 0.96250 (0.95471–0.96919) | 0.99974 (0.99942–0.99998) | 1.00000 (1.00000–1.00000) |
| A6 | `renamed` | 1.00000 (1.00000–1.00000) | 1.00000 (1.00000–1.00000) | 0.96604 (0.96029–0.97120) | 1.00000 (1.00000–1.00000) | – |
| A6 | `patch` | – | 0.90288 (0.87143–0.93515) | 0.94965 (0.93627–0.96209) | – | – |

Ids without any pattern (`plain`) are answered right in A1 in every organism; size and depth of a lineage (`lineage5`, A*'s ids) do not defeat idtrack; a missing successor does (`retired_none`, `handoff_dies`, mouse's `patch`: copies of GRCm38 patches that GRCm39 ends). `wrong_units.tsv` lists, per group and lineage, the wrong answers behind these rates.

### Causes

Every wrong answer explained on the tree that gave it (`explain_a`) and classified by F0's classifier with A's plugin; answers summed over the strata of `results/a/<run>/<organism>/causes_summary.tsv`:

| run | cause | class | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `d4c2dda` | `external_bridge`: no successor in Ensembl's history, answered through an external id with a gene that does not overlap the old place | b | 194,727 | 485,211 | 759,611 | 10,352 | 0 |
| `d4c2dda` | `path_through_resurrection`: an id retired and brought back, read through its return | e | – | 180 | – | – | – |
| `d4c2dda` | `copies_apart_before_82`: copies on patches before release 82, which the truth labels primary | d | – | – | 11,639 | – | – |
| `d4c2dda` | `gtf_lacks_core_gene`: the target's GTF lacks a gene the core has | d | – | – | 2,687 | – | – |
| `d4c2dda` | `event_version_contradicts_tables`: one event whose version contradicts the tables (rule `event_version`) | c | – | – | 602 | – | – |
| `d4c2dda` | `transcript_read_elsewhere`: a transcript id asked from a release is read at another release or through a synonym | a | – | – | 11 | – | – |
| `d4c2dda` | `core_version_on_versionless_release`: the truth's version of a core-only gene at a release without versions | d | – | – | 6 | – | – |
| `d4c2dda` | `unexplained`, read by hand: the translation id equal to a transcript id (A0 above) and a transcript read at another release | – | – | – | 2 | – | – |
| `d4c2dda_s110` | `external_bridge` | b | 14,314 | – | 31,266 | – | – |
| `d4c2dda_s110` | `event_version_contradicts_tables` | c | – | – | 84 | – | – |

The 24 human answers of A7 that the judge once called wrong because it did not reduce LRG gene ids (`LRG_<n>_g1.<v>`) are right since the judge's stable id form reads LRG ids (`harness_lib.STABLE_ID_FORM`, 2026-11-20; A was scored again from its stored answers).

**Where Ensembl's history has a successor, idtrack is right.** Every wrong answer of pig and rhesus, all of mouse's but the 180 read through a resurrection, and most of human's (cynomolgus has none) is an id whose strict walk reaches nothing, answered with a gene that a shared external id links to it and that does not overlap the old place (class b; rule [`external_bridge`](../rules.md#external_bridge), which the decided rule [`bridge_same_place`](../rules.md#bridge_same_place) restricts). The truth classes (c, d, e) are Ensembl's or the judge's. Over the whole family idtrack is wrong under its own rules (class a) in 11 answers of human transcripts read at another release or through a synonym, and in the 2 answers read by hand. A reproducer of the bridges: rhesus `convert_identifier("ENSMMUG00000036773", from_release=65, to_release=48)` on `d4c2dda` (`API.convert` from `71690ab` on) answers `ENSMMUG00000024688`, `ENSMMUG00000026611` and `ENSMMUG00000026730`, none at the id's place; Ensembl's history has no predecessor of the id, which is born at 58.

## Known limitations

A's rows of [the known-limitations table](../../results/limitations/limitations.tsv), generated by `limits_a` from the causes (wrong answers) and from the judge's reasons of the undetermined answers; answers per organism on the 115 graphs (the `d4c2dda` runs, human with A8 on the 115 graph):

| limitation | what | outcome | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `L-bridges` | no successor in Ensembl's history: a bridge through a shared id to a gene elsewhere (until the decided rule `bridge_same_place`) | wrong | 194,727 | 485,211 | 759,611 | 10,352 | – |
| `L-resurrection` | ids retired and brought back (rule `resurrection`, by design) | wrong | – | 180 | – | – | – |
| `L-copies-before-82` | copies on patches before release 82, which the GTF does not mark | wrong | – | – | 11,639 | – | – |
| `L-event-versions` | events whose version contradicts the tables (until the decided rule `event_version`) | wrong | – | – | 602 | – | – |
| `L-transcript-elsewhere` | a transcript id asked from a release is read at another release or as a synonym | wrong | – | – | 11 | – | – |
| `L-truth-gtf-short` | genes of a release missing from its GTF (the truth's) | wrong | – | – | 2,687 | – | – |
| `L-unexplained` | wrong answers read by hand (above) | wrong | – | – | 2 | – | – |
| `T-core-version` | the truth's version of a core-only gene at a release without versions (not idtrack's) | wrong | – | – | 6 | – | – |
| `L-no-chain` | no successor and no assembly chain: the answer cannot be judged (Ensembl publishes no chain) | undetermined | 2,714,444 | – | – | 2,586,456 | 974,950 |
| `L-bridge-at-place` | no successor; idtrack answers a gene at the same place (plausible, not Ensembl's statement; until `bridge_same_place`) | undetermined | 301,601 | 1,155,605 | 1,960,844 | – | – |
| `L-chain-gap` | no successor and the place does not map through the chain | undetermined | 277,712 | 112,898 | 154,782 | – | – |
| `L-branch` | answers on branches that left the id while it lived on (until the decided rule `main_line_first`) | undetermined | 50,236 | 6,235 | 73,312 | – | – |
| `L-open-truth` | the truth leaves the line open | undetermined | – | – | 17,598 | – | – |
| `L-undetermined-other` | other undetermined answers | undetermined | – | 4,654 | 10,092 | – | – |

On the snapshot-110 graphs (A8): pig `L-bridges` 14,314, `L-no-chain` 485,976, `L-bridge-at-place` 35,958, `L-branch` 1,466; human `L-bridges` 31,266, `L-event-versions` 84, `L-bridge-at-place` 164,874, `L-chain-gap` 16,744, `L-open-truth` 2,980, `L-undetermined-other` 980, `L-branch` 84. Each row also carries its units, the answers of the strata it occurs in, the share, the strata, the run and the graph's digest.

## Reproducing from scratch

From an empty data root, with the configuration naming the tree `d4c2dda` and its graphs (and `d4c2dda_s110` with the snapshot-110 graphs of pig and human):

1. **Fetch the truth** as A* does ([A*'s steps](astar.md#reproducing-from-scratch): the harness's fetches and the core gene tables of human 48 to 51 and cynomolgus 112); A0 also reads BioMart's tables of every release and VGNC's sets (`fetch_tables ... biomart`, `... vgnc`). A reads A*'s selection only for its pattern `astar`.
2. **A0**: `a0_fidelity` on each organism's graph.
3. **Publish the sets**: `gen_a <o> --stratum A1 --set a1_v1`, and so for A2 to A6 and A8 (`--snapshot 110`), `gen_a7 <o> --set a7_v1`; then `patterns_a <o>`. A rebuild with the current code writes the truth of every revision above; it is compared with the published sets' latest lines and any difference is published as a revision set (`revise_a`), never by rewriting a set.
4. **Build the snapshot-110 graphs** of pig and human with the tree (the same tables, the snapshot's history table fetched by the builder).
5. **Run and score**: the runner on `d4c2dda` for every organism, and on `d4c2dda_s110` with the A8 sets for pig and human.
6. **Report**: `report_a` on every run.
7. **Classify**: `explain_a`, then `causes` with A's plugin.
8. **Limitations**: `limits_a --out <tsv>` over the seven runs.
9. **Results**: `idtrack-repro results` writes `results/a/d4c2dda/<organism>/` and `results/a/d4c2dda_s110/<organism>/` (with `a0_summary.tsv` from A0), the known-limitations table, the index and the provenance.
