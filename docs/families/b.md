# B: names over time

Family B asks idtrack to read gene names as Ensembl wrote them at every release: every display name of every release read at that release, carried to the target release and read without a release; every display name at its own release by the database Ensembl filed it under (B9); names that moved to another gene; names of several genes; the authorities' previous symbols, synonyms and aliases; the forms a name takes in pasted lists (letter case, white space, separators, invisible characters, full-width letters, spreadsheet dates, Greek letters, suffixes); and every named gene to its organism's name database at every release. Its code is `idtrack_reproducibility.families.b`; its data are `<families>/b/` (`queries/<organism>/`, `runs/<tree>/<organism>/` with the scorings `scores_now/` and `scores_decided/`); the design is in [the methods](../methods.md#b-names-over-time). It reads the history as A* does ([A*'s truth](astar.md#truth)) and runs on the common machinery ([F0](f0.md)).

## State

What is built:

- the truth readers (`truth_b`), the generators of the sets (`gen_b`: B1 with B2 and B3 as strata, and B8; `gen_b2`: B4 to B7), the fetcher of the authorities' files and of Ensembl's synonyms (`fetch_b`), the builder of revision sets (`revise_b`), the report of a run (`report_b`) and the plugin of the common classifier (`causes_b`, which registers B's kinds of lines; B has no cause of its own yet);
- the truth files B fetched, each in its organism's `SOURCES.tsv`: HGNC's complete set and withdrawn list (and the listing of HGNC's quarterly archive), MGI's `MRK_List1.rpt`, `MRK_List2.rpt` and `MRK_ENSEMBL.rpt`, NCBI's `gene_info` of the five taxa, and Ensembl's synonyms of every cross-reference of release 115 (`core_synonyms/115.tsv.gz`) and, for mouse, of every release of its window (48 to 115, fetched for family E's letter-case lines; Ensembl serves the core table `external_synonym` at every release, as `.txt.gz.bz2` at 89);
- the published sets of all five organisms (version 1; [the sets](#the-sets));
- B9 (`gen_b9`, set `b9_v1` of all five organisms, published 2026-11-21) and its truth, each gene's display cross-reference of every release (`core_display/<R>.tsv.gz`, `fetch_b display`), with Ensembl's synonyms of every cross-reference of every release fetched in the same pass for the four organisms that had only release 115's ([B9](#b9-display-names-at-their-own-release)).

What was run and scored, on the reference tree `d4c2dda` and its graphs:

| run folder | state |
| --- | --- |
| `runs/d4c2dda/macaca_fascicularis/` | every line answered; `scores_now/` holds `run.json` and `report.tsv` (sets `b1_v1`, `b8_v1`, `b5_v1`, `b6_v1`, `b7_v1`); `scores_decided/` not scored; no causes |
| `runs/d4c2dda/macaca_mulatta/` | every line answered; `scores_now/` holds `run.json` and `report.tsv` (the same sets and `b4_v1`); `scores_decided/` not scored; no causes |
| `runs/d4c2dda/mus_musculus/` | incomplete: shards 0 to 24 of 29 complete (8,464,522 lines; `partial.json` says so), shards 25 to 28 not asked, not scored; the answers are kept in the data root (manifested), and the runner, which asks only what has no answer, finishes them with `--shard 25/29` to `28/29` |
| `runs/d4c2dda/sus_scrofa/`, `runs/d4c2dda/homo_sapiens/` | not run (no answer) |

The judge following idtrack's decided rules (judge `rules`, [`judge_follows_rules`](../rules.md#judge_follows_rules); published 2026-11-20): `<set>_rules1`, the accuracy harness's revisions `missing_values`, `capitals_mgi` and `lower_case_symbols` (rule `single_case_symbols`), `stray_punctuation` and `par_shared_names` applied to B's latest lines (`revise_b --rules`, B's truth context): pig `b5_rules1` 12 (missing values, left out of the scoring); mouse `b4_rules1` 21, `b5_rules1` 1,156, `b6_rules1` 2,312, `b1_rules1` 105 (capitals and lower case read as MGI's symbol, missing values, stray punctuation); human `b4_rules1` 35, `b5_rules1` 415, `b6_rules1` 168, `b1_rules1` 464 (the names the X and Y copies of a pseudoautosomal gene share: the Y copy apart; missing values; stray punctuation); rhesus and cynomolgus none. B6's decided judge reads its `_dec` sets after the `now` sets, so its lines get the same revisions in `b6_rules1_dec` (mouse 2,270: the decided lines of the qids `b6_rules1` holds, published again after it so that the decided truth stands in publication order too; human 23 with the Y copy apart). On the modest selection every revised line `next` answers is right but one (mouse `MUP1` from 105, `all`: the second gene its capitals reach, wrong on every tree).

After idtrack's next version's first round (the tree `r1`, published 2026-11-21; [the modest set](../modest.md#the-tree-r1)): `<set>_rules2` (judge `rules`), the harness's revision `mirbase_genes` (rule [`mirbase_genes`](../rules.md#mirbase_genes): a miRBase entry of two loci is read through Ensembl's other databases, or is not found) on B's latest lines: pig `b1` 981, `b5` 18, `b6` 98 (`b6_rules2_dec` 30), `b9` 874; mouse `b1` 1,283, `b9` 661; human `b1` 715, `b5` 263, `b6` 5 (`_dec` 1), `b9` 371; rhesus `b1` 1,555, `b5` 27, `b6` 137 (`_dec` 38), `b9` 947; cynomolgus none. Two corrections of the judge (judge `78`; `revise_b --correction`, the module's docstring): `b1_places1` (`places_not_comparable`: a name several genes carry, which Ensembl's history carries to nothing, read on an assembly that no chain file joins to the target's, is open, as the line of a single carrier is undetermined by the locus rule: pig 957, rhesus 892, cynomolgus 92, mouse and human none (their assemblies have chain files); `b1_places2` puts back the line before `b1_places1` of the names `mirbase_genes` says are not held there, which it had opened too: pig 56, rhesus 136) and `b1_support1` (`support_from_xrefs`: rule `support_ranks_shared_name` needs each expected gene's support, which the harness reads from BioMart's tables that most early releases lack; a B1 line of a name several genes carry gets it from Ensembl's records at the release read, every database of the core cross-references, the synonyms and the display name: pig 3,161, mouse 8,392, human 20,052, rhesus 2,770, cynomolgus 106; `b1_snames1`, `support_names`, gives every such line the support of the name sources alone, as idtrack counts it (`Candidate.entry_support`: the databases attaching the name read; a synonym, BioMart's in the generator's support or the core's in `support_from_xrefs`'s, is none), by the procedure that gave its support: pig 354, mouse 272, human 495, rhesus 416, cynomolgus none; mouse `Mtcp1` at 78, whose two genes the core's count made equal by a UniProt gene-name synonym of ENSMUSG00000090110, is so undetermined where idtrack ranks ENSMUSG00000031200 first by its Vega gene label), and `b4_held1`, `b5_held1` (`held_in_core`: B4's and B5's lines of a string the truth found in no Ensembl table, generated before the core extracts of every release were fetched, that those extracts hold at some release are read as rule `start_release` reads a string, at the latest release holding it, names before synonyms: pig 3 and 1,032, mouse 15,374 and 11,855, human 179 and 15,565, rhesus 2 and 879, cynomolgus none and 853; on pig's modest lines every tree's answer moves from wrong to right on most of them, `r1`'s on none from right to wrong; a line whose reading rule `mirbase_genes` decided is left alone, and `b5_held2` puts back the 2 rhesus lines `b5_held1` had read otherwise).

The revision `r1` (published 2026-11-24, judge `b_r1`; `revise_b --name r1`): every set of version 1 regenerated with the truth the code writes since it was published, and the lines whose truth changed published under their qids, after the corrections of the judge published for each set since its base set were replayed on the regenerated lines in their publication order (`revise_b.replay`), so that no correction is undone; the comparison is of the truth and the question, not of the wording (`revise_b.WORDING`: `note`, `name_src`; `open` and `excluded` by whether they are set). What it reads anew: Ensembl's synonyms of every cross-reference of every release beside BioMart's synonyms (`b1`'s `broad`; B5's synonyms of BioMart releases, whose synonyms Ensembl's records add are new questions, `b5_v2`, `name_src: core_synonyms`); rule [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name)'s support from Ensembl's name sources at releases without a BioMart table; a string asked without a release that no name or synonym of 115 holds read at the latest release whose Ensembl records hold it (`gen_b2.Ctx.read_where_held`, as the correction `held_in_core` read the strings it found; before, such a line was open within every gene that ever carried the string); B8's `spellings` of one value; the original name's `broad` on B6's lines judged `set`; B6's and B7's `safe` lines accepting the original name's genes in any name source; and rule `no_mito_nuclear_step` (decided; A*'s history), which leaves a name of a mitochondrial gene that Ensembl's history carries to a nuclear gene without that successor (`b1`'s and `b2`'s lines of the macaques' mitochondrial genes). Lines per set (sample of 50 per set, stratum and kind of change read against the truth files: every added gene of a `broad` carries the string at the release read in a name source or as a synonym, every dropped one in none; every line read where held is held there and at no later release up to 115):

| set | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `b1_r1` | 12,396 | 28,924 | 47,597 | 14,229 | 3,087 |
| `b4_r1` | 4 | 6,978 | 675 | 4 | – |
| `b5_r1` | 1,180 | 2,129 | 7,981 | 1,780 | 864 |
| `b5_v2` (new lines) | 144,286 | 285,934 | 260,392 | 68,590 | – |
| `b6_r1` | 1,248 | 3,517 | 16,945 | 183 | 211 |
| `b6_r1_dec` | 291 | 554 | 1,564 | 33 | 43 |
| `b7_r1` | 67 | 64 | 121 | 15 | 21 |
| `b7_r1_dec` | – | 1 | 2 | – | – |
| `b8_r1` | 600 | – | 11 | 25 | 45 |
| `b8_v2` (new lines) | – | 52,386 | 48,311 | – | – |

No revision changes a question: every changed line keeps its query; `b7_r1` leaves out the suffix forms that Ensembl's synonyms of 115 now hold as strings of their own (`excluded`: pig 1, cynomolgus 3), and `b1_r1` mouse's 6 names no longer asked; `b8_v2` holds mouse's genes of release 89 and human's of 112, which version 1 read through the reach. The stored answers of `d4c2dda` (cynomolgus, rhesus, mouse's 25 shards) answer every revised line; `b5_v2`'s lines were never asked.

What is not built: the runs of pig, mouse and human, the `decided` scorings, the explained paths and causes of the wrong answers, the family's limitations map and its rows of the known-limitations table, the authority measure (`report_b`, written, not run), and the family's panel ([Measurement design](../methods.md#measurement-design)).

## What it measures

- **B1, every name at every release**: whether a display name of the GTF of release R, read at R (`from_release=R`), is answered with the genes that carry it at R; carried to the target release 115, with their descendants by Ensembl's history; and, read without a release (rule [`start_release`](../rules.md#start_release)), every display name of 115.
- **B2, moved names**: a name read at R whose carrier at 115 is no descendant of its carrier at R (the old gene's history is the truth), and the default reading of a name that ever moved (its carrier now).
- **B3, names of several genes**: the B1 lines whose name several genes carry at the release (rule [`ties`](../rules.md#ties) and the support of rule [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name)).
- **B4, previous symbols** of the authorities, without a release and from the last release whose GTF displays them.
- **B5, synonyms and aliases**: every BioMart synonym of every BioMart release, read at that release and carried to 115, and every alias of the authorities without a release.
- **B6, forms** of every display name of 115: letter case, white space, separators, invisible characters, full-width letters.
- **B7, spreadsheet and typography**: the dates a spreadsheet makes of date-like symbols, Greek letters, suffixes.
- **B8, names to databases**: every gene with a display name at a release, converted to its organism's name database at that release.
- **B9, display names at their own release**: every display name of every release as Ensembl's core database records it (each gene's display cross-reference, with the database Ensembl filed the name under), read at that release; stratified by the filing database and the release's era, so that a name source the graph does not hold at a release is seen by its database (rule [`display_names_at_release`](../rules.md#display_names_at_release)).

The rate of right answers is right / (right + wrong + error), undetermined answers and timeouts apart, per stratum and organism, with Wilson intervals and a cluster bootstrap over the unit (the string; B8: the gene). Every line is scored under `best` and under `all` (`strategy: both`). Two judges read the same answers: `now` scores the sets without `_dec` in their name (what the reference trees document), `decided` scores every set, with `b6_dec` and `b7_dec` replacing their lines by the decided rules. Where an authority (HGNC, MGI, VGNC, NCBI) gives a line a gene, it rides along as the line's `authority` and is measured beside idtrack (how often Ensembl and the authority agree, how often idtrack gives the authority's gene), never as a wrong answer of idtrack (principle 2 of [the methods](../methods.md#principles)).

## The sets

One set per family of strata and organism, published once (F0's `lines.publish`); a change of truth is a revision set.

| set | strata | what is asked | truth |
| --- | --- | --- | --- |
| `b1_v1` | `b1_atR_<c>`, `b1_at115_<c>`, `b1_to115_<c>`, `b1_default_<c>` (`<c>`: `unique` one gene, `shared` several genes (B3), `patchonly` only copies on patches carry it); `b2_moved_from`, `b2_moved_default` (B2) | every display name of every release R (not an Ensembl-like id, not a number), read at R, to R and to 115; every display name of 115 without a release | the genes carrying it in the GTF of R (rule `set`, `set_exact` for several), then their walks; [the truth](#truth) |
| `b8_v1` | `b8_biomart`, `b8_core_xref` | every gene with a display name at R, with `final_database` the organism's name database (human `HGNC Symbol`, mouse `MGI Symbol`, the others `UniProtKB Gene Name`), `cmp: exact` | BioMart's values of that database at R where an archive answers, else the labels of Ensembl's cross-references of R |
| `b4_v1` | `b4_prev_default`, `b4_prev_from` | every previous symbol (HGNC `prev_symbol` and the withdrawn symbols merged into an approved one, MGI's withdrawn symbols, VGNC `prev_symbol`), without a release and from the last release whose GTF displays it | without a release: read at 115, names before synonyms; a string that no name or synonym of 115 holds is read at the latest release whose Ensembl records hold it, its genes walked to 115 (since `r1`; version 1: `open` within the walks of every gene that ever carried it, `open_within`, so a gene never linked to it is wrong; still so for a string that a gene of 115 carries in another letter case); from a release: B1's truth of that string at that release |
| `b5_v1` | `b5_syn_atR`, `b5_syn_at115`, `b5_syn_to115`, `b5_alias_<authority>` | every BioMart synonym of every BioMart release R, read at R, to R and to 115 (`b5_v2`: the synonyms of every cross-reference that Ensembl's records add at R); every alias of HGNC, MGI, VGNC and NCBI without a release | the synonym's genes after rule [`name_before_synonym`](../rules.md#name_before_synonym), then their walks; aliases as B4's lines without a release |
| `b6_v1`, `b6_dec` | `b6_case`, `b6_space` (every display name of 115), `b6_sep`, `b6_invisible`, `b6_fullwidth` (20,000 sampled (name, form) pairs per organism) | lower, upper and title case; the name with spaces around it and with a tab and CR/LF; dashes and underscores swapped, replaced by a dot or a space, or removed; a non-breaking space, a zero-width space, a byte-order mark; full-width letters | `set` with the original's genes for what the reference trees document (rule [`spelling_steps`](../rules.md#spelling_steps): letter case and the dash/underscore swap; white space trimmed); `safe` (nothing or the original name's genes) for forms no built rule reads; a form that is an exact name of 115 is left out (its own question), one that is an exact synonym of another gene expects that gene (rule [`names_before_synonyms`](../rules.md#names_before_synonyms)), one that another gene carries in another letter case is `open`. `b6_dec`: rule [`unicode_normalised`](../rules.md#unicode_normalised) (the gene) and, in mouse, rule [`single_case_symbols`](../rules.md#single_case_symbols) |
| `b9_v1`, `b9_r1`, `b9_r2` | `b9_<era>_<db>` (a name of one gene; `<db>` the database of its display cross-reference, `patchonly` where only copies on patches carry it), `b9_shared_<era>` (a name several genes carry at once) | every display name of every release R whose core database is read, as each gene's display cross-reference at R gives it, read at R (`from_release=R`, `to_release=R`) | the genes of R's GTF on the primary assembly whose display cross-reference is the name (rule `set` for one gene and for several alike: under `best` a non-empty part of the carriers is right, under `all` exactly the carriers); copies on patches and haplotypes apart (before 82 by Ensembl's assembly exceptions) and the Y copy of a pseudoautosomal gene apart (`b9_r1`); a display name that is an LRG record's expects the Ensembl gene the record describes (rule `lrg_answers_gene`, `b9_r2`); `broad`: the other genes that carry it at R (a display name off the GTF, a label of a name database, a synonym of R); compared without versions ([B9](#b9-display-names-at-their-own-release)) |
| `b7_v1`, `b7_dec` | `b7_date`, `b7_greek`, `b7_suffix` | the forms a spreadsheet writes of every date-like symbol (a month and a day: `SEPT7`, `MARCH1`, `DEC1`, `OCT4`) of any release or authority (`7-Sep`, `07-Sep`, `Sep-07`, `2007-09-07`, `9/7/2007`, the serial `39332`); a curated list of 51 Greek-letter spellings (`TNFα`, `IL-1β`, `IFN-γ`, `TGF-β1`, ...); `-AS1`, `_at`, `_s_at` on 2,000 sampled names of one gene | `safe`: the symbol's genes or nothing is right, another gene wrong; a form that is a name or synonym of 115 is left out. `b7_dec` (the dates and Greek letters): the one gene the form stands for, else nothing ([Q1](../methods.md#open-methodological-choices)) |

qids: `b9|<name>|<R>`, `b1|<name>|<R>|<R or 115>`, `b1|<name>|default|115`, `b8|<gene>|<R>|<database>`, `b4|<symbol>|default|115`, `b4|<symbol>|<R>|115`, `b5|syn|<synonym>|<R>|<R or 115>`, `b5|<authority>|<alias>|default|115`, `b6|<form kind>|<name>`, `b7|date|<form kind>|<month>|<day>`, `b7|greek|<form>`, `b7|suffix|<kind>|<name>`.

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
| `b9_v1` | 886,641 | 2,996,268 | 3,391,264 | 1,183,916 | 407,922 |
| `b9_r1` | – | 395 | 32,012 | – | – |
| `b9_r2` | – | – | 32,023 | – | – |
| `b7_dec` | 202 | 255 | 303 | 201 | 172 |

Cynomolgus has no `b4_v1`: its authority file lists no previous symbol. The human and mouse sets were generated before family C's extracts of Ensembl's cross-references of human release 112 and mouse release 89 existed, so those two releases are read through the reach ([the truth](#truth)); their revision sets will read them.

## Truth

Nothing of idtrack; every file is in its organism's `SOURCES.tsv` (`truth_b.NamesB`):

- **Display names at every release**: the GTF's `gene_name` of every release of the window (pig 56 to 115, human, mouse and rhesus 48 to 115, cynomolgus 91 to 115), with the region where the release's GTF marks patches (human and mouse from 82): copies on patches are reported apart (rule [`judge_copies_apart`](../rules.md#judge_copies_apart)).
- **Other name sources at R**: where an archive BioMart answers, its HGNC, MGI, UniProt gene name, NCBI label and WikiGene tables, VGNC's symbol of BioMart's VGNC id and the Vega names, as the accuracy harness reads them; and at every release Ensembl's cross-references of the release (`core_xref/<R>.tsv.gz`, family C's extract of the release's MySQL dumps): the labels of HGNC, MGI, VGNC, UniProt gene name, NCBI gene, WikiGene, Vega and the clone-based names. RFAM and miRBase labels are no names (rules [`mirbase_genes`](../rules.md#mirbase_genes), [`rfam_display_label`](../rules.md#rfam_display_label)). A gene carrying the string in another name source at R goes into `broad` (the corrected judge's `names_broad`, rule [`judge_corrections`](../rules.md#judge_corrections)).
- **Synonyms at R**: BioMart's synonym table where an archive answers (it holds the synonyms of a gene's display cross-reference only), and Ensembl's synonyms of every cross-reference where `core_synonyms/<R>.tsv.gz` is fetched (every release of every window; the sets of version 1 read release 115's, the revision `r1` every release's).
- **The reach**: where the synonyms of R are not read, the genes that carry the string as a synonym at any release read; where no name source of R is read, the genes that carry it as a name at any release. An answer within the reach is open, one outside it wrong; every line says what was read (`name_src`: `biomart`, `core_xref`, `none+reach`, with `+synonyms` or `+synonyms_reach`).
- **The history**: A*'s walks (`astar.truth.TruthA`): the strict walk with living ids (rule [`id_lives_on`](../rules.md#id_lives_on)) is `expect`, the broad walk with ids brought back (rule [`resurrection`](../rules.md#resurrection)) goes into `broad`; where the history ends for a name of one gene, the locus rule with the corrected judge's places across assemblies, assembly exceptions and core places; one gene that split is judged `set` (the corrected judge's `split_name_set`); mouse display names that are projected HGNC symbols and no name in another source are left out (`projected_hgnc`).
- **Display cross-references** (B9): each gene's display cross-reference of every release as the release's core database holds it (`gene.display_xref_id` -> `xref` -> `external_db`, read from the MySQL dumps, `core_display/<R>.tsv.gz`): the name, the gene and the database it is filed under; Ensembl's GTF writes this label as `gene_name`, and the two differ at a few releases only (genes of the GTF whose `gene_name` is not the label): mouse 49, where the GTF writes the MGI accession (`MGI:2687328`) of 21,182 genes for their symbol; release 78 of pig, mouse and human, where the label of a UniProt gene name ends in a stray `;` the GTF leaves out (418, 13, 26 genes); human 50, where the GTF adds a `.` to 412 HGNC projections; and where a GTF writes no name or an older one (cynomolgus 108: 15,619 genes; cynomolgus 112: 120; pig 105: 405; mouse 103 and 105: 24). B9 asks the label as the core database writes it.
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
# each gene's display cross-reference of every release (B9's truth; the synonyms of every release in the same pass)
python -m idtrack_reproducibility.families.b.fetch_b display --releases all [--organisms <o>,...]
# the sets (published once): B1 and B8 read every release's tables (a job of the configured runner)
python -m idtrack_reproducibility.families.b.gen_b --organism <o> [--sets b1,b8] [--dry]
python -m idtrack_reproducibility.families.b.gen_b2 --organism <o> [--sets b4,b5,b6,b7] [--dry]
python -m idtrack_reproducibility.families.b.gen_b9 --organism <o> [--dry]
# revision sets when the truth changes: only the lines whose truth changed, under their published qids (the published
# corrections replayed on the regenerated lines); a dry run writes them beside their published lines, and a revision
# can be published from its dry run's file (B1 regenerates release by release: hours for mouse and human)
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name r1 --why "<what changed>" [--sets b4,b5,b6,b7,b8,b1] [--dry --lines-to <file>]
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name r1 --why "<what changed>" --sets b1 --publish-from <file>
# the harness's revisions of decided rules on B's lines (judge `rules`; with `_dec` sets for B6 and B7)
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name rules1 \
    --rules missing_values,capitals_mgi,lower_case_symbols,stray_punctuation,par_shared_names
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name rules2 --rules mirbase_genes --sets b1,b4,b5,b6,b7,b8,b9
# corrections of the judge on B1 (places across assemblies without a chain; support at releases without BioMart)
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name places1 --correction places_not_comparable --sets b1
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name support1 --correction support_from_xrefs --sets b1
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name snames1 --correction support_names --sets b1
python -m idtrack_reproducibility.families.b.revise_b --organism <o> --name held1 --correction held_in_core --sets b4,b5
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

Runs `<families>/b/runs/d4c2dda/macaca_fascicularis/scores_now/` and `<families>/b/runs/d4c2dda/macaca_mulatta/scores_now/` (the reference graphs, the published sets of version 1, the judge `now`); result files `results/b/d4c2dda/<organism>/report.tsv` and `scores.tsv`. No answer timed out and none raised an error. These are the only complete scored runs of B: pig, mouse and human are not run, no `decided` scoring and no cause exists, and they are judged by the sets of version 1 and the corrections before `r1` (which was published after them), so no wrong answer below is attributed yet to idtrack or to the truth.

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

### B9: display names at their own release

Runs `<families>/b/runs/next/<organism>/` (the tree `next`, idtrack `4471396`, on its own graphs `graphs/next/<organism>/`; only B9 asked; scored by the judge `now` with the sets `b9_v1` and `b9_r1` (human also `b9_r2`) into `scores_now/`, 2026-11-21); result files `results/b/next/<organism>/now_report.tsv` and `now_scores.tsv` (scope `sets` in `results/index.tsv`). The reference tree was not run on B9 (the set is newer than its results), so these lines have no reference answers; `next` is measured before the fix of its first round, and the same sets on the fixed tree measure the fix. No answer timed out and none raised.

`b9_r1` (`revise_b --name r1 --sets b9`: B9 regenerated, the lines that changed published under their qids) holds the lines of the copies the accuracy harness's revisions put apart, which `b9_v1` expected inside the answer: before release 82 the GTFs mark no copy on a haplotype or patch, so a name that the reference gene and its copies carried at once was a name of several genes; Ensembl's assembly exceptions tell them (the corrected judge's `copies_apart_early`), and the Y copy of a pseudoautosomal gene goes apart where the X gene carries the name (rule `par_shared_names`): human 32,012 lines, mouse 395, none in the other three (`<families>/b/queries/<organism>/b9_r1.manifest.json`; the regeneration reproduced every other line of `b9_v1`).

`b9_r2` (human; `revise_b --name r2 --sets b9`) brings the judge up to the decided rule [`lrg_answers_gene`](../rules.md#lrg_answers_gene), which `next` builds: a display name that is an LRG record's (11 lines of release 58, `LRG_1` to `LRG_15`) expects the Ensembl gene the record describes, the gene Ensembl's core database of the release links to the LRG gene by its cross-reference of database `ENS_LRG_gene` ("LRG display in Ensembl", whose accession is the LRG gene id; `core_xref/58.tsv.gz`), where `b9_v1` expected the LRG gene itself. It holds 32,023 lines: the 11, and the 32,012 lines of `b9_r1` again with the same truth (the revision builder then compared the earlier revision's note with the regenerated line; it now compares the truth and the question only, so a later revision publishes only what changed).

The question of every line is which gene a display name is at the release it is a name of. By era (`early`: the releases before every database filing a display name is one that files them at 115; pig before 90, rhesus 73, mouse 108, human 104, cynomolgus none), view `best`, lines (questions), from `now_report.tsv`:

| organism | era | lines | wrong | undetermined | rate |
| --- | --- | ---: | ---: | ---: | ---: |
| pig | early | 462,826 | 47,672 | 3,256 | 0.8963 |
| pig | late | 409,768 | 104 | 428 | 0.9997 |
| pig | names of several genes, early · late | 11,784 · 2,263 | 1,342 · 76 | 93 · 0 | 0.8852 · 0.9664 |
| mouse | early | 2,511,601 | 11,771 | 10,400 | 0.9953 |
| mouse | late | 467,724 | 11 | 837 | 1.0000 |
| mouse | names of several genes, early · late | 15,119 · 1,824 | 748 · 0 | 613 · 40 | 0.9484 · 1.0000 |
| human | early | 2,882,208 | 17,959 | 8,646 | 0.9938 |
| human | late | 483,000 | 14 | 448 | 1.0000 |
| human | names of several genes, early · late | 22,443 · 3,613 | 971 · 0 | 1,588 · 39 | 0.9534 · 1.0000 |
| rhesus | early | 442,974 | 41,878 | 1,252 | 0.9052 |
| rhesus | late | 726,886 | 605 | 1,958 | 0.9992 |
| rhesus | names of several genes, early · late | 6,570 · 7,486 | 697 · 95 | 158 · 21 | 0.8913 · 0.9873 |
| cynomolgus | late | 405,748 | 0 | 205 | 1.0000 |
| cynomolgus | names of several genes, late | 2,174 | 0 | 4 | 1.0000 |

By the database the name is filed under, the strata of view `best` below 0.99 hold the fault (`now_report.tsv`): every UniProt entry name is lost (`b9_early_uniprot_sptrembl`, `_swissprot`, `_varsplic`: pig 43,017 lines, rhesus 18,759, mouse 885, human 4,726), every `RFAM_gene_name` and `miRBase_gene_name` of mouse and human 61 and 62 (mouse 6,784 lines all wrong; human 9,249 of 9,366 wrong, 5 right, 112 undetermined), and the HGNC projections of pig (`b9_early_hgnc_curated_gene` 0.9229 and `_automatic_gene` 0.9283 of 41,586 and 15,021 lines) and rhesus (0.9498 and 0.9363 of 99,994 and 58,273), rhesus's RefSeq and NCBI labels (`b9_early_refseq_dna_predicted` 0.5148 of 16,531, `b9_early_entrezgene` 0.9581 of 96,340), the miRBase names of every era (rhesus `b9_late_mirbase` 0.9675, human `b9_early_mirbase` 0.9101) and the immunoglobulin and T-cell receptor segments of IMGT (mouse 0.6498 of 1,628, human 0.8204 of 1,152). Of the wrong answers under `best`, those whose name was not read at its release (the answer says `from_release_fallback`, or the name is not found) are the class of the rule [`display_names_at_release`](../rules.md#display_names_at_release): pig 48,938 of 49,194 (46,660 lost, 2,278 another gene), rhesus 43,104 of 43,275, mouse 12,015 of 12,530, human 15,829 of 18,944; the rest are names the graph holds at the release that it answers with another gene (pig 256, mouse 515, human 3,115, rhesus 171): miRBase entries (rule [`mirbase_genes`](../rules.md#mirbase_genes)), names other genes carry in a source the truth does not read, RFAM accessions that are display labels of several genes (rule [`rfam_display_label`](../rules.md#rfam_display_label) names only those), and in human before 82 a haplotype copy answered for the reference gene, or the reference gene without its copies reported apart. Human 58's 11 display names that are LRG records' (`b9_early_lrg`) expect, by `b9_r2`, the Ensembl gene each record describes (rule [`lrg_answers_gene`](../rules.md#lrg_answers_gene)), and `next` answers all 11 so. Measured beside the audit that found the class on the same tree (every display name of the GTF of every release read at its release; lines lost or answered with another gene: pig 47,382, rhesus 42,176, mouse 30,483, human 15,391): B9 counts about as many for pig, rhesus and human (48,938, 43,104, 15,829: B9 reads the display label of the core database at every release and asks the names of several genes too), and fewer in mouse, where the GTF of release 49 writes the MGI accession for the symbol of 21,182 genes (`MGI:2687328` for Sly) that the audit read as names and B9, reading the core database's display label, does not ask.

The rows of the known-limitations table these lines fall under until the fix lands (`results/limitations/b_next.tsv`; answers of both views): `L-display-name-not-held` (the name not read at its release: pig 98,738 answers, mouse 27,394, human 32,295, rhesus 89,181) and `L-display-name-other` (held, answered with another gene: pig 920, mouse 3,095, human 10,683, rhesus 505), measured from the scored outcomes and the answers by `limits_outcomes` with B's map (`families/b/limitations_map_outcomes.tsv`).

## Known limitations

B's rows of the known-limitations table are those of B9 on `next` (above: `L-display-name-not-held`, `L-display-name-other`); the wrong answers of the sets of version 1 are not classified. What the truth cannot judge, and what the published sets leave open (the synonyms of every cross-reference, the support at releases without BioMart and human's 112 and mouse's 89, read through the reach by version 1, are read by the revision `r1`):

| limitation | what | size |
| --- | --- | --- |
| moved names read at an old release | a moved name read at an old release is often answered with its carrier now, which the truth puts in `broad` when the old release's cross-references give that gene the name in another source (the corrected judge's `names_broad` at every release): such lines are undetermined | `b2_moved_from` undetermined: cynomolgus 20,878 of 24,092 answers, rhesus 185,218 of 203,036 (`results/b/d4c2dda/<organism>/report.tsv`) |
| the reach | where a release's synonyms or name sources are not read, the truth accepts the genes that carry the string at any release read (open) and calls an answer outside them wrong; since `r1` every release's are read, and the lines of version 1 whose truth the revision left as it was keep the `name_src` they were made with | every line says what was read (`name_src`) |
| the authority | where HGNC, MGI, VGNC or NCBI and Ensembl disagree, the disagreement is Ensembl's and is measured beside idtrack (`authority_b.tsv`), never a wrong answer | not measured (`report_b` not run) |

## Reproducing from scratch

From an empty data root, with the configuration naming the tree `d4c2dda` and its graphs:

1. **Fetch the truth**: the accuracy harness's fetches (GTFs of every release, BioMart's archive tables, the core tables of the history, VGNC's gene sets; [the accuracy harness](../harness.md#truth)), family C's extracts of Ensembl's cross-references of every release ([C's truth](c.md#truth)), then `fetch_b hgnc`, `fetch_b mgi`, `fetch_b ncbi`, `fetch_b synonyms --releases all` and `fetch_b display --releases all`.
2. **Publish the sets**: `gen_b --organism <o>` (B1, B8), `gen_b2 --organism <o>` (B4 to B7), `gen_b9 --organism <o>` (B9). A rebuild with the current code writes the truth of the revisions above; it is compared with the published sets and every difference is published as a revision set (`revise_b`), never by rewriting a set.
3. **Run** every set on `d4c2dda` for every organism (F0's runner).
4. **Score** both judges: the sets without `_dec` into `scores_now/`, every set into `scores_decided/`.
5. **Classify**: `explain` on the tree, then `causes` with B's plugin.
6. **Report**: `report_b` (both judges, the authority measure, outcomes by `name_src`).
7. **Results**: `idtrack-repro results` writes `results/b/d4c2dda/<organism>/` once the runs are in `results.toml`; B9's run on a tree (`--sets b9_v1`) is registered with scope `sets` and its limitation rows come from `limits_outcomes` (`limits = "outcomes"`).
