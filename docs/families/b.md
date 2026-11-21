# B: names over time

Family B asks idtrack to read gene names as Ensembl wrote them at every release: every display name of every release read at that release, carried to the target release and read without a release; names that moved to another gene; names of several genes; the authorities' previous symbols, synonyms and aliases; the forms a name takes in pasted lists (letter case, white space, separators, invisible characters, full-width letters, spreadsheet dates, Greek letters, suffixes); and every named gene to its organism's name database at every release. Its code is `idtrack_reproducibility.families.b`; its data are `<families>/b/` (`queries/<organism>/`, `runs/<tree>/<organism>/` with the scorings `scores_now/` and `scores_decided/`); the design is in [the methods](../methods.md#b-names-over-time). It reads the history as A* does ([A*'s truth](astar.md#truth)) and runs on the common machinery ([F0](f0.md)).

## State

What is built:

- the truth readers (`truth_b`), the generators of the sets (`gen_b`: B1 with B2 and B3 as strata, and B8; `gen_b2`: B4 to B7), the fetcher of the authorities' files and of Ensembl's synonyms (`fetch_b`), the builder of revision sets (`revise_b`), the report of a run (`report_b`) and the plugin of the common classifier (`causes_b`, which registers B's kinds of lines; B has no cause of its own yet);
- the truth files B fetched, each in its organism's `SOURCES.tsv`: HGNC's complete set and withdrawn list (and the listing of HGNC's quarterly archive), MGI's `MRK_List1.rpt`, `MRK_List2.rpt` and `MRK_ENSEMBL.rpt`, NCBI's `gene_info` of the five taxa, and Ensembl's synonyms of every cross-reference of release 115 (`core_synonyms/115.tsv.gz`) and, for mouse, of every release of its window (48 to 115, fetched for family E's letter-case lines; Ensembl serves the core table `external_synonym` at every release, as `.txt.gz.bz2` at 89);
- the published sets of all five organisms (version 1; [the sets](#the-sets)).

What was run and scored, on the reference tree `d4c2dda` and its graphs:

| run folder | state |
| --- | --- |
| `runs/d4c2dda/macaca_fascicularis/` | every line answered; `scores_now/` holds `run.json` and `report.tsv` (sets `b1_v1`, `b8_v1`, `b5_v1`, `b6_v1`, `b7_v1`); `scores_decided/` not scored; no causes |
| `runs/d4c2dda/macaca_mulatta/` | every line answered; `scores_now/` holds `run.json` and `report.tsv` (the same sets and `b4_v1`); `scores_decided/` not scored; no causes |
| `runs/d4c2dda/mus_musculus/` | incomplete: shards 0 to 24 of 29 complete (8,464,522 lines; `partial.json` says so), shards 25 to 28 not asked, not scored; the answers are kept in the data root (manifested), and the runner, which asks only what has no answer, finishes them with `--shard 25/29` to `28/29` |
| `runs/d4c2dda/sus_scrofa/`, `runs/d4c2dda/homo_sapiens/` | not run (no answer) |

The judge following idtrack's decided rules (judge `rules`, [`judge_follows_rules`](../rules.md#judge_follows_rules); published 2026-11-20): `<set>_rules1`, the accuracy harness's revisions `missing_values`, `capitals_mgi` and `lower_case_symbols` (rule `single_case_symbols`), `stray_punctuation` and `par_shared_names` applied to B's latest lines (`revise_b --rules`, B's truth context): pig `b5_rules1` 12 (missing values, left out of the scoring); mouse `b4_rules1` 21, `b5_rules1` 1,156, `b6_rules1` 2,312, `b1_rules1` 105 (capitals and lower case read as MGI's symbol, missing values, stray punctuation); human `b4_rules1` 35, `b5_rules1` 415, `b6_rules1` 168, `b1_rules1` 464 (the names the X and Y copies of a pseudoautosomal gene share: the Y copy apart; missing values; stray punctuation); rhesus and cynomolgus none. B6's decided judge reads its `_dec` sets after the `now` sets, so its lines get the same revisions in `b6_rules1_dec` (mouse 2,270: the decided lines of the qids `b6_rules1` holds, published again after it so that the decided truth stands in publication order too; human 23 with the Y copy apart). On the modest selection every revised line `next` answers is right but one (mouse `MUP1` from 105, `all`: the second gene its capitals reach, wrong on every tree).

What is not built:

- Ensembl's synonyms of every release for the organisms other than mouse (`fetch_b synonyms --releases all`; release 115 is fetched for all five, every release for mouse), and B's truth read from them (it reads release 115's);
- the revision sets of the truth the code writes since version 1 was published: Ensembl's synonyms of every cross-reference beside BioMart's synonyms, rule [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name)'s support from Ensembl's cross-references at releases without a BioMart table, B8's `spellings` of one value, the original name's `broad` on B6's lines judged `set`, and B6's and B7's `safe` lines accepting the original name's genes in any name source. `revise_b` was tried on cynomolgus without publishing; until these revision sets are published, the published sets are the truth of the results below;
- the runs of pig, mouse and human, the `decided` scorings, the explained paths and causes of the wrong answers, the family's limitations map and its rows of the known-limitations table, the authority measure (`report_b`, written, not run), and the family's panel ([Measurement design](../methods.md#measurement-design)).

## What it measures

- **B1, every name at every release**: whether a display name of the GTF of release R, read at R (`from_release=R`), is answered with the genes that carry it at R; carried to the target release 115, with their descendants by Ensembl's history; and, read without a release (rule [`start_release`](../rules.md#start_release)), every display name of 115.
- **B2, moved names**: a name read at R whose carrier at 115 is no descendant of its carrier at R (the old gene's history is the truth), and the default reading of a name that ever moved (its carrier now).
- **B3, names of several genes**: the B1 lines whose name several genes carry at the release (rule [`ties`](../rules.md#ties) and the support of rule [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name)).
- **B4, previous symbols** of the authorities, without a release and from the last release whose GTF displays them.
- **B5, synonyms and aliases**: every BioMart synonym of every BioMart release, read at that release and carried to 115, and every alias of the authorities without a release.
- **B6, forms** of every display name of 115: letter case, white space, separators, invisible characters, full-width letters.
- **B7, spreadsheet and typography**: the dates a spreadsheet makes of date-like symbols, Greek letters, suffixes.
- **B8, names to databases**: every gene with a display name at a release, converted to its organism's name database at that release.

The rate of right answers is right / (right + wrong + error), undetermined answers and timeouts apart, per stratum and organism, with Wilson intervals and a cluster bootstrap over the unit (the string; B8: the gene). Every line is scored under `best` and under `all` (`strategy: both`). Two judges read the same answers: `now` scores the sets without `_dec` in their name (what the reference trees document), `decided` scores every set, with `b6_dec` and `b7_dec` replacing their lines by the decided rules. Where an authority (HGNC, MGI, VGNC, NCBI) gives a line a gene, it rides along as the line's `authority` and is measured beside idtrack (how often Ensembl and the authority agree, how often idtrack gives the authority's gene), never as a wrong answer of idtrack (principle 2 of [the methods](../methods.md#principles)).

## The sets

One set per family of strata and organism, published once (F0's `lines.publish`); a change of truth is a revision set.

| set | strata | what is asked | truth |
| --- | --- | --- | --- |
| `b1_v1` | `b1_atR_<c>`, `b1_at115_<c>`, `b1_to115_<c>`, `b1_default_<c>` (`<c>`: `unique` one gene, `shared` several genes (B3), `patchonly` only copies on patches carry it); `b2_moved_from`, `b2_moved_default` (B2) | every display name of every release R (not an Ensembl-like id, not a number), read at R, to R and to 115; every display name of 115 without a release | the genes carrying it in the GTF of R (rule `set`, `set_exact` for several), then their walks; [the truth](#truth) |
| `b8_v1` | `b8_biomart`, `b8_core_xref` | every gene with a display name at R, with `final_database` the organism's name database (human `HGNC Symbol`, mouse `MGI Symbol`, the others `UniProtKB Gene Name`), `cmp: exact` | BioMart's values of that database at R where an archive answers, else the labels of Ensembl's cross-references of R |
| `b4_v1` | `b4_prev_default`, `b4_prev_from` | every previous symbol (HGNC `prev_symbol` and the withdrawn symbols merged into an approved one, MGI's withdrawn symbols, VGNC `prev_symbol`), without a release and from the last release whose GTF displays it | without a release: read at 115, names before synonyms; a string that no name or synonym of 115 holds is `open` within the walks of every gene that ever carried it (`open_within`), so a gene never linked to it is wrong; from a release: B1's truth of that string at that release |
| `b5_v1` | `b5_syn_atR`, `b5_syn_at115`, `b5_syn_to115`, `b5_alias_<authority>` | every BioMart synonym of every BioMart release R, read at R, to R and to 115; every alias of HGNC, MGI, VGNC and NCBI without a release | the synonym's genes after rule [`name_before_synonym`](../rules.md#name_before_synonym), then their walks; aliases as B4's lines without a release |
| `b6_v1`, `b6_dec` | `b6_case`, `b6_space` (every display name of 115), `b6_sep`, `b6_invisible`, `b6_fullwidth` (20,000 sampled (name, form) pairs per organism) | lower, upper and title case; the name with spaces around it and with a tab and CR/LF; dashes and underscores swapped, replaced by a dot or a space, or removed; a non-breaking space, a zero-width space, a byte-order mark; full-width letters | `set` with the original's genes for what the reference trees document (rule [`spelling_steps`](../rules.md#spelling_steps): letter case and the dash/underscore swap; white space trimmed); `safe` (nothing or the original name's genes) for forms no built rule reads; a form that is an exact name of 115 is left out (its own question), one that is an exact synonym of another gene expects that gene (rule [`names_before_synonyms`](../rules.md#names_before_synonyms)), one that another gene carries in another letter case is `open`. `b6_dec`: rule [`unicode_normalised`](../rules.md#unicode_normalised) (the gene) and, in mouse, rule [`single_case_symbols`](../rules.md#single_case_symbols) |
| `b7_v1`, `b7_dec` | `b7_date`, `b7_greek`, `b7_suffix` | the forms a spreadsheet writes of every date-like symbol (a month and a day: `SEPT7`, `MARCH1`, `DEC1`, `OCT4`) of any release or authority (`7-Sep`, `07-Sep`, `Sep-07`, `2007-09-07`, `9/7/2007`, the serial `39332`); a curated list of 51 Greek-letter spellings (`TNFα`, `IL-1β`, `IFN-γ`, `TGF-β1`, ...); `-AS1`, `_at`, `_s_at` on 2,000 sampled names of one gene | `safe`: the symbol's genes or nothing is right, another gene wrong; a form that is a name or synonym of 115 is left out. `b7_dec` (the dates and Greek letters): the one gene the form stands for, else nothing ([Q1](../methods.md#open-methodological-choices)) |

qids: `b1|<name>|<R>|<R or 115>`, `b1|<name>|default|115`, `b8|<gene>|<R>|<database>`, `b4|<symbol>|default|115`, `b4|<symbol>|<R>|115`, `b5|syn|<synonym>|<R>|<R or 115>`, `b5|<authority>|<alias>|default|115`, `b6|<form kind>|<name>`, `b7|date|<form kind>|<month>|<day>`, `b7|greek|<form>`, `b7|suffix|<kind>|<name>`.

Lines per set and organism (manifests `<families>/b/queries/<organism>/<set>.manifest.json`):

| set | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `b1_v1` | 1,772,670 | 5,993,260 | 6,782,498 | 2,367,832 | 784,658 |
| `b8_v1` | 721,484 | 3,031,531 | 3,755,866 | 978,314 | 528,174 |
| `b4_v1` | 318 | 85,366 | 22,533 | 479 | – |
| `b5_v1` | 14,547 | 368,810 | 381,523 | 7,024 | 3,840 |
| `b6_v1` | 85,648 | 333,059 | 185,817 | 92,300 | 83,608 |
| `b6_dec` | 19,537 | 21,637 | 17,515 | 19,355 | 19,815 |
| `b7_v1` | 6,201 | 6,255 | 6,203 | 6,201 | 6,171 |
| `b7_dec` | 202 | 255 | 303 | 201 | 172 |

Cynomolgus has no `b4_v1`: its authority file lists no previous symbol. The human and mouse sets were generated before family C's extracts of Ensembl's cross-references of human release 112 and mouse release 89 existed, so those two releases are read through the reach ([the truth](#truth)); their revision sets will read them.

## Truth

Nothing of idtrack; every file is in its organism's `SOURCES.tsv` (`truth_b.NamesB`):

- **Display names at every release**: the GTF's `gene_name` of every release of the window (pig 56 to 115, human, mouse and rhesus 48 to 115, cynomolgus 91 to 115), with the region where the release's GTF marks patches (human and mouse from 82): copies on patches are reported apart (rule [`judge_copies_apart`](../rules.md#judge_copies_apart)).
- **Other name sources at R**: where an archive BioMart answers, its HGNC, MGI, UniProt gene name, NCBI label and WikiGene tables, VGNC's symbol of BioMart's VGNC id and the Vega names, as the accuracy harness reads them; and at every release Ensembl's cross-references of the release (`core_xref/<R>.tsv.gz`, family C's extract of the release's MySQL dumps): the labels of HGNC, MGI, VGNC, UniProt gene name, NCBI gene, WikiGene, Vega and the clone-based names. RFAM and miRBase labels are no names (rules [`mirbase_genes`](../rules.md#mirbase_genes), [`rfam_display_label`](../rules.md#rfam_display_label)). A gene carrying the string in another name source at R goes into `broad` (the corrected judge's `names_broad`, rule [`judge_corrections`](../rules.md#judge_corrections)).
- **Synonyms at R**: BioMart's synonym table where an archive answers (it holds the synonyms of a gene's display cross-reference only), and Ensembl's synonyms of every cross-reference where `core_synonyms/<R>.tsv.gz` is fetched (B's published sets read release 115's).
- **The reach**: where the synonyms of R are not read, the genes that carry the string as a synonym at any release read; where no name source of R is read, the genes that carry it as a name at any release. An answer within the reach is open, one outside it wrong; every line says what was read (`name_src`: `biomart`, `core_xref`, `none+reach`, with `+synonyms` or `+synonyms_reach`).
- **The history**: A*'s walks (`astar.truth.TruthA`): the strict walk with living ids (rule [`id_lives_on`](../rules.md#id_lives_on)) is `expect`, the broad walk with ids brought back (rule [`resurrection`](../rules.md#resurrection)) goes into `broad`; where the history ends for a name of one gene, the locus rule with the corrected judge's places across assemblies, assembly exceptions and core places; one gene that split is judged `set` (the corrected judge's `split_name_set`); mouse display names that are projected HGNC symbols and no name in another source are left out (`projected_hgnc`).
- **The authorities** (beside the truth, never the judge): HGNC's complete set (symbol, `prev_symbol`, `alias_symbol`, Ensembl gene id) and withdrawn list; MGI's `MRK_List1.rpt` (withdrawn symbols with their current marker, synonyms) and `MRK_ENSEMBL.rpt` (MGI's own Ensembl ids); VGNC's gene sets; NCBI's `gene_info` (synonyms, Ensembl ids of `dbXrefs`).

## Rules

The rules of idtrack the expectations rest on: [`start_release`](../rules.md#start_release) (lines without a release), [`spelling_steps`](../rules.md#spelling_steps) and [`names_before_synonyms`](../rules.md#names_before_synonyms) (B4 to B6), [`spelling_of_release`](../rules.md#spelling_of_release) and [`string_at_release`](../rules.md#string_at_release) (B1), [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name) and [`ties`](../rules.md#ties) (B3), [`own_id_first`](../rules.md#own_id_first), [`copies_apart`](../rules.md#copies_apart), [`external_bridge`](../rules.md#external_bridge); the judge's [`name_before_synonym`](../rules.md#name_before_synonym), [`judge_copies_apart`](../rules.md#judge_copies_apart), [`id_lives_on`](../rules.md#id_lives_on), [`resurrection`](../rules.md#resurrection) and [`judge_corrections`](../rules.md#judge_corrections) (`names_broad`, `living_ids`, `split_name_set`, `projected_hgnc`, `places`, `exceptions`, `exceptions_span`, `core_places`). The decided rules B measures before they land (`b6_dec`): [`unicode_normalised`](../rules.md#unicode_normalised) and, in mouse, [`single_case_symbols`](../rules.md#single_case_symbols); `b7_dec` scores the spreadsheet forms and Greek letters by one answer to [Q1](../methods.md#open-methodological-choices).

Scoring is F0's ([scoring and outcomes](f0.md#scoring-and-outcomes)): `set`, `set_exact`, `safe`, `open` with `open_within`, `cmp: exact` for B8 with its `spellings`.

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)). The judge `now` scores the sets without `_dec`, the judge `decided` every set; `--sets` lists the sets in publication order (a later set replaces the lines of its qids).

```
# the authorities' files and Ensembl's synonyms (once; one host at a time)
python -m idtrack_reproducibility.families.b.fetch_b hgnc|mgi|ncbi
python -m idtrack_reproducibility.families.b.fetch_b synonyms --releases 115|all|<r,r,...> [--organisms <o>,...]
# the sets (published once): B1 and B8 read every release's tables (a job of the configured runner)
python -m idtrack_reproducibility.families.b.gen_b --organism <o> [--sets b1,b8] [--dry]
python -m idtrack_reproducibility.families.b.gen_b2 --organism <o> [--sets b4,b5,b6,b7] [--dry]
# revision sets when the truth changes: only the lines whose truth changed, under their published qids
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name r1 --why "<what changed>" [--sets b4,b5,b6,b7,b8,b1] [--dry]
# the harness's revisions of decided rules on B's lines (judge `rules`; with `_dec` sets for B6 and B7)
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name rules1 \
    --rules missing_values,capitals_mgi,lower_case_symbols,stray_punctuation,par_shared_names
# a tree's run (F0's runner), then the scoring of each judge into its own folder
python -m idtrack_reproducibility.families.runner --family b --organism <o> --clone <source> --graph <graph> \
    --out <families>/b/runs/<tree>/<o> --sets <every set> --submit <n> [--tag b_<tree>_<o>]
python -m idtrack_reproducibility.families.runner --family b --organism <o> --clone <source> --graph <graph> \
    --out <families>/b/runs/<tree>/<o> --submit 0 --score --score-parts <P> --sets <sets without _dec> \
    --scores-to <families>/b/runs/<tree>/<o>/scores_now
python -m idtrack_reproducibility.families.runner ... --submit 0 --score --score-parts <P> --sets <every set> \
    --scores-to <families>/b/runs/<tree>/<o>/scores_decided
# causes: the wrong answers explained on the tree that gave them, then F0's classifier with B's kinds
<python of the tree> -m idtrack_reproducibility.families.explain --family b --organism <o> --clone <source> --graph <graph> \
    --scored <families>/b/runs/<tree>/<o>/scores_now --per-unit 3
python -m idtrack_reproducibility.families.causes --organism <o> --scored <families>/b/runs/<tree>/<o>/scores_now \
    --plugin src/idtrack_reproducibility/families/b/causes_b.py
# both judges side by side, the authority measure, outcomes by the tables the truth read
python -m idtrack_reproducibility.families.b.report_b --organism <o> --run <families>/b/runs/<tree>/<o>
```

`report_b` writes `report_b.tsv`, `authority_b.tsv` and `name_src_b.tsv` into the run folder.

## Results

Runs `<families>/b/runs/d4c2dda/macaca_fascicularis/scores_now/` and `<families>/b/runs/d4c2dda/macaca_mulatta/scores_now/` (the reference graphs, the published sets of version 1, the judge `now`); result files `results/b/d4c2dda/<organism>/report.tsv` and `scores.tsv`. No answer timed out and none raised an error. These are the only complete scored runs of B: pig, mouse and human are not run, no `decided` scoring and no cause exists, and the truth's revisions above are not yet published, so no wrong answer below is attributed yet to idtrack or to the truth.

Rate and cluster-bootstrap interval per stratum (row `*` of `report.tsv`: both views):

| stratum | cynomolgus | rhesus |
| --- | ---: | ---: |
| `b1_atR_unique` | 1.0000 (0.9999–1.0000) | 0.9599 (0.9586–0.9612) |
| `b1_at115_unique` | 1.0000 (1.0000–1.0000) | 0.9939 (0.9928–0.9950) |
| `b1_to115_unique` | 1.0000 (1.0000–1.0000) | 0.9774 (0.9762–0.9785) |
| `b1_default_unique` | 1.0000 (1.0000–1.0000) | 0.9930 (0.9917–0.9943) |
| `b1_atR_shared` | 0.9777 (0.9687–0.9852) | 0.7733 (0.7558–0.7917) |
| `b1_at115_shared` | 1.0000 (1.0000–1.0000) | 0.9196 (0.8718–0.9623) |
| `b1_to115_shared` | 0.9586 (0.9410–0.9741) | 0.8240 (0.8032–0.8430) |
| `b1_default_shared` | 1.0000 (1.0000–1.0000) | 0.9545 (0.8864–1.0000) |
| `b2_moved_from` | 0.9676 (0.9522–0.9808) | 0.7511 (0.7303–0.7723) |
| `b2_moved_default` | 0.9968 (0.9941–0.9991) | 0.9900 (0.9873–0.9925) |
| `b4_prev_default` | – | 1.0000 (1.0000–1.0000) |
| `b4_prev_from` | – | 0.9978 (0.9934–1.0000) |
| `b5_syn_atR` | 1.0000 (1.0000–1.0000) | 0.9954 (0.9884–1.0000) |
| `b5_syn_at115` | 1.0000 (1.0000–1.0000) | 0.9964 (0.9910–1.0000) |
| `b5_syn_to115` | 1.0000 (1.0000–1.0000) | 0.9954 (0.9883–1.0000) |
| `b5_alias_ncbi` | 0.6491 (0.6294–0.6704) | 0.6895 (0.6730–0.7061) |
| `b5_alias_vgnc` | 1.0000 (1.0000–1.0000) | 0.7778 (0.4444–1.0000) |
| `b6_case` | 1.0000 (1.0000–1.0000) | 0.9938 (0.9926–0.9949) |
| `b6_space` | 0.9981 (0.9975–0.9987) | 0.9927 (0.9915–0.9939) |
| `b6_sep` | 1.0000 (1.0000–1.0000) | 0.9333 (0.9123–0.9528) |
| `b6_invisible` | 0.9995 (0.9991–0.9998) | 0.9992 (0.9988–0.9996) |
| `b6_fullwidth` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `b7_date` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `b7_greek` | 1.0000 (1.0000–1.0000) | 1.0000 (1.0000–1.0000) |
| `b7_suffix` | 0.9998 (0.9995–1.0000) | 1.0000 (1.0000–1.0000) |
| `b8_biomart` | 0.9999 (0.9999–1.0000) | 1.0000 (1.0000–1.0000) |
| `b8_core_xref` | 1.0000 (0.9999–1.0000) | 1.0000 (0.9999–1.0000) |

Answers, wrong and undetermined of the largest strata and of those with the most wrong answers (same table):

| stratum | cynomolgus | rhesus |
| --- | ---: | ---: |
| `b1_atR_unique` | 748,620 / 31 / 791 | 2,303,890 / 92,105 / 8,630 |
| `b1_to115_unique` | 724,836 / 17 / 7,751 | 2,104,244 / 45,521 / 92,229 |
| `b1_atR_shared` | 4,222 / 94 / 11 | 27,870 / 6,207 / 490 |
| `b1_to115_shared` | 3,914 / 162 / 3 | 24,480 / 4,267 / 231 |
| `b2_moved_from` | 24,092 / 104 / 20,878 | 203,036 / 4,435 / 185,218 |
| `b5_alias_ncbi` | 6,066 / 1,519 / 1,737 | 9,464 / 1,840 / 3,538 |
| `b6_space` | 63,632 / 120 / 0 | 72,144 / 524 / 0 |
| `b8_core_xref` | 924,288 / 37 / 0 | 1,722,396 / 53 / 0 |

In cynomolgus the names of one gene have 31 wrong answers in 748,620 read at their release and 17 in 724,836 carried to 115; its wrong answers lie in the names of several genes, the moved names and NCBI's aliases. In rhesus the names read at their own release are wrong 92,105 times (`b1_atR_unique`, 5,946 names behind them), against 31 in cynomolgus; which of these the graph and which the truth's reading of rhesus's older releases explains is open until the causes are classified. NCBI's aliases are the weakest stratum in both organisms; the truth reads BioMart's synonyms of display cross-references only, and Ensembl's synonyms of every cross-reference only at 115 ([Known limitations](#known-limitations)).

## Known limitations

B has no rows in the known-limitations table yet: its wrong answers are not classified. What the truth cannot judge, and what the published sets leave open:

| limitation | what | size |
| --- | --- | --- |
| synonyms of every cross-reference | BioMart's synonym table holds the synonyms of a gene's display cross-reference only, while idtrack reads Ensembl's synonyms of every cross-reference; the truth reads those only at 115 (`core_synonyms/115.tsv.gz`), so an answer through a synonym of another release that the truth did not read is wrong until the synonyms of every release are fetched and published as a revision | not measured (no causes) |
| support at releases without BioMart | the published sets give rule `support_ranks_shared_name`'s support only at BioMart releases, so a `best` answer that the support would leave undetermined elsewhere is judged by the set; the code reads it from Ensembl's cross-references (revision not published) | not measured |
| moved names read at an old release | a moved name read at an old release is often answered with its carrier now, which the truth puts in `broad` when the old release's cross-references give that gene the name in another source (the corrected judge's `names_broad` at every release): such lines are undetermined | `b2_moved_from` undetermined: cynomolgus 20,878 of 24,092 answers, rhesus 185,218 of 203,036 (`results/b/d4c2dda/<organism>/report.tsv`) |
| the reach | where a release's synonyms or name sources are not read, the truth accepts the genes that carry the string at any release read (open) and calls an answer outside them wrong | every line says what was read (`name_src`) |
| human 112 and mouse 89 | the published human and mouse sets read these releases through the reach (their cross-references were not extracted yet) | lines of those releases in `b1_v1`, `b8_v1` of human and mouse |
| the authority | where HGNC, MGI, VGNC or NCBI and Ensembl disagree, the disagreement is Ensembl's and is measured beside idtrack (`authority_b.tsv`), never a wrong answer | not measured (`report_b` not run) |

## Reproducing from scratch

From an empty data root, with the configuration naming the tree `d4c2dda` and its graphs:

1. **Fetch the truth**: the accuracy harness's fetches (GTFs of every release, BioMart's archive tables, the core tables of the history, VGNC's gene sets; [the accuracy harness](../harness.md#truth)), family C's extracts of Ensembl's cross-references of every release ([C's truth](c.md#truth)), then `fetch_b hgnc`, `fetch_b mgi`, `fetch_b ncbi` and `fetch_b synonyms --releases all`.
2. **Publish the sets**: `gen_b --organism <o>` (B1, B8), `gen_b2 --organism <o>` (B4 to B7). A rebuild with the current code writes the truth of the revisions above; it is compared with the published sets and every difference is published as a revision set (`revise_b`), never by rewriting a set.
3. **Run** every set on `d4c2dda` for every organism (F0's runner).
4. **Score** both judges: the sets without `_dec` into `scores_now/`, every set into `scores_decided/`.
5. **Classify**: `explain` on the tree, then `causes` with B's plugin.
6. **Report**: `report_b` (both judges, the authority measure, outcomes by `name_src`).
7. **Results**: `idtrack-repro results` writes `results/b/d4c2dda/<organism>/` once the runs are in `results.toml`.
