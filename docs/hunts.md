# The weakness hunts

The benchmark families measure what was foreseen; the weakness hunts look for what was not. In a hunt, an analyst searches one angle of idtrack's use with whatever probes the angle needs, and every failure found becomes a finding: a minimal reproducer, the answer an independent truth expects, the answer idtrack gives, the cause read from the explained path and the code, a count of what the class touches, and a disposition. This page gives the method, the three hunts that have run (inputs, history, regions) with their findings and what each became, and how to reproduce them. The design of the hunts within the benchmark programme is in [the methods](methods.md#weakness-hunts).

The three hunts ran on idtrack `71690ab` with the five reference graphs of the benchmarks (built by `f27e810`, whose builder `d4c2dda` keeps: [the data](data.md#layout)). The two trees give the same answer to every reproducer below ([Reproducing](#reproducing)), and family E asked the inputs findings' lines on both. Every count below is the hunt's own, on that tree, unless another tree is named.

| angle | findings | the findings table | the hunt's outputs (data root) |
| --- | --- | --- | --- |
| inputs | IN-1 to IN-11 | [hunts/inputs/findings.tsv](../hunts/inputs/findings.tsv) | `<hunts>/inputs/` |
| history | hh-1 to hh-6 | [hunts/history/findings.tsv](../hunts/history/findings.tsv) | `<hunts>/history/` |
| regions | R1 to R9, J1 to J4 | [hunts/regions/findings.tsv](../hunts/regions/findings.tsv) | `<hunts>/regions/` |

## Method

### Angles

- **inputs**: what users type and paste: missing values in every form, white space and invisible characters, Unicode look-alikes and Greek letters, letter case and separators, names mangled by spreadsheets, versions and suffixes, ids of another organism or database, very long and very many inputs, duplicates and order in a batch, source inference on lists of known origin, every argument at and beyond its bounds, determinism (hash seeds, threads, processes, batches against single calls), the explained path against the answer.
- **history**: the genes whose histories are hardest, mined from Ensembl's own history (the release-115 `stable_id_event` and `mapping_session` tables and the GTFs of every release) and travelled every way: from each release where an id exists to every release, forward and back, `best` and `all`, with and without a target database, against the strict and broad walks of Ensembl's history.
- **regions**: the places where genomes are hard: pseudoautosomal genes, patches, alternative loci and LRG records, mitochondrial genes, segmental duplications and gene families (HLA, olfactory receptors, KRTAP, histones, amylases, immunoglobulin and T-cell receptor segments), multi-copy RNAs, readthroughs, unplaced scaffolds, Y genes and their X gametologs, genes that changed coordinates across an assembly change; for each, what idtrack answers for its ids and names over time.

Angles for later rounds (external identifiers, names over time, the explained path, concurrency and the web layer, harmonization, the docs against behaviour, slow queries) are listed in [the methods](methods.md#weakness-hunts).

### How a finding is made

1. **Search.** Probes ask idtrack's public api (calls through `API.convert`, `convert_many`, `infer_source`, `HarmonizeFeatures`) on many inputs at once, and compare each answer with a truth or with a relation that must hold (the answer of the published form of an identifier, the same node asked from two releases, the strict walk of Ensembl's history).
2. **Reproduce.** Each class of failure is cut to one call or a few, on a named graph.
3. **Expect.** The expected answer comes from a truth source, never from idtrack: Ensembl's GTFs, core tables and history, BioMart, NCBI, HGNC, MGI, UniProt, GENCODE, LRG. Where no truth can decide (an id Ensembl's history carries to another chromosome), the finding says so.
4. **Explain.** The cause is read from the explained path (`explain=True`) and the code.
5. **Count.** A probe counts what the class touches (genes, inputs or calls), never estimates it.
6. **Grade.** Severity: *high*, a wrong answer given with confidence on ordinary input, or a crash; *medium*, wrong on unusual input, or the right answer with wrong evidence; *low*, the rest (mostly answers lost as not found). *known*: the class was known before the hunt; *question*: idtrack's behaviour is a choice the truth does not settle; *judge* and *the truth*: a fault of the benchmarks' truth or judge, not of idtrack.

### The reproducer format

Each angle's `findings.tsv` has one row per finding with the columns:

| column | what |
| --- | --- |
| `id` | the finding's id (`IN-<n>`, `hh-<n>`, `R<n>`, `J<n>`) |
| `severity` | as above |
| `organisms` | where it occurs |
| `how_many` | the hunt's counts, with the tree named and the probe that counted them |
| `reproducer` | `reproduce <angle> <id>`: the calls of the runnable reproducer ([Reproducing](#reproducing)) |
| `expected` | the answer the truth expects, with the truth's source |
| `got` | the answer on `d4c2dda` and on `71690ab` |
| `cause` | from the explained path and the code |
| `disposition` | what the finding became (below) |
| `measured_by` | the family, set and stratum (or pattern) whose lines measure the class now, as `<family>/<sets>/<strata>` |

The reproducers are one table in `idtrack_reproducibility.hunts.reproduce` (the calls, the expectation with its truth's source, the answer recorded), run by one command. Where a finding's class is a pattern, its lines are published in the family it belongs to and scored there with every tree: the inputs findings as family E's set `hunt_inputs` (stratum `hunt_<id>`, judge `now`) with its revision set `hunt_inputs_dec` (judge `decided`), and as strata of E's generator ([family E](families/e.md)); the history and regions findings' ids in family A*'s selection (pattern `hunt_history`, and patterns made from the findings: `branch_while_living`, `non_reference_copy_ended`, `par_y`; [family A*](families/astar.md)) and family A's pattern flags (`retired_none`, `gap`, `version_drop`, `moved_seq`, `strand_flip`, `patch`, `mt`, `core_only`; [family A](families/a.md)); the judge's faults as corrections of the judge ([the corrections of the judge](harness.md#the-corrections-of-the-judge)).

### Triage and its five outcomes

Findings are merged where they share a cause; then each becomes one of five outcomes, written in the column `disposition`:

| outcome | `disposition` |
| --- | --- |
| a fix of idtrack under a rule it already follows | `rule:<name>` |
| a new rule of idtrack, decided for its next version | `rule:<name>` (the rule's status in [rules](rules.md) says which) |
| a known limitation, its size measured by a family | `limitation:<L-id>` (a row of the known-limitations table, `results/limitations/limitations.tsv`) |
| a truth error, corrected by a revision set with its evidence | `truth_error:<correction>` ([judge_corrections](rules.md#judge_corrections)) |
| not a defect, with the reason | `not_a_defect` |

Of the 30 findings of the three hunts (the dispositions of the findings tables), 23 became rules decided for idtrack's next version (none of them is in `71690ab`, so every reproducer of these still holds on both trees), 6 became corrections of the judge, and 1 (R5) is a limitation; none was found not to be a defect.

### Rounds and stopping

A round is one full pass of every active angle on one tree, starting from the previous round's findings. Hunting on a tree stops when a pass finds nothing that is not a duplicate by cause in any angle. The first round ran on `71690ab` (with `d4c2dda`'s graphs): the inputs hunt's two final passes (1,000 ids per organism in 34 and 50 forms, fresh seeds) and the regions hunt's third pass (on genes no earlier pass had chosen) found nothing new; the history hunt travelled every selected id exhaustively and closed with genome-wide scans of the patterns it found and truth-free checks of every selected id. The next round runs on idtrack's next version, starting from these findings; each reproducer then says whether its finding still holds.

## Inputs

**What was searched.** Probes on idtrack's public api, five organisms unless named; their outputs are in `<hunts>/inputs/outbox/<probe>.out` (and `<hunts>/inputs/det/` for determinism):

| probe | what | how much |
| --- | --- | --- |
| p01 | hand-written forms of TP53, HLA-A, MARCH1 and SEPT7: spreadsheet forms, versions, suffixes, ids of other organisms and databases, internal names, wildcards, missing values (human) | about 200 forms |
| p02 | missing-value and header tokens | 86 tokens |
| p03 | metamorphic: published identifiers of 11 kinds (names, synonyms, case and separator variants, NCBI, RefSeq, UniProt, gene ids with and without version, transcripts, translations) in up to 20 mangled forms; a mangled form must give the published form's answer or nothing | 2,000 per kind |
| p04, p05 | every current symbol with the suffixes of make-unique (`-1`, `-2`, `.1`, `.2`, `_1`); every current symbol in upper, capitalised and lower case | every symbol |
| p06, p09, p10 | source inference on lists of a known release (Ensembl's GTF of that release): ids, versioned ids, names, whole releases, mixtures, junk; the cause | human, mouse |
| p07 | every argument at and beyond its bounds, batches in every container and with missing values, very long strings (human) | |
| det | determinism and agreement: published queries asked alone, in batches (in order, reversed, doubled), in 8 threads, repeated, with and without `explain`, in two processes with different hash seeds; every explained path checked | 26,741 queries |
| p11, p13, p14 | identifiers of each organism typed into every other organism's graph; numeric ids per database | 1,500 per kind |
| p12, p20 | one batch of 516,671 identifiers; every gene id and name of 4 human and 3 mouse releases converted alone | 396,206 and 378,504 calls |
| p15 | `HarmonizeFeatures` on a 14-feature human matrix holding the forms of IN-1 and IN-3 | |
| p16, p17, p18, p22 | names equal up to separators; accession forms; prefixed ids and Greek letters | |
| p21 | real inputs: the feature names of an atlas project's 33 author matrices | 33 matrices |
| p23, p24, p25, p26 | final passes (400 symbols, 300 synonyms, 300 gene ids per organism, fresh seeds, in 34 and 50 forms, every answer that differs from the published form's explained); names with stray punctuation in the graph and their sources | 1,000 ids per organism |

**What held** (no finding): every argument at and beyond its bounds is refused with an input error naming it; answers are deterministic across hash seeds, threads, processes, batches and repeats; every explained path starts at the graph's id, ends at its candidate and is contiguous; Unicode spaces, tabs and line ends are stripped; ids of another organism are not found (except the numbers of IN-9); a batch of 516,671 identifiers runs without error.

**Findings** ([hunts/inputs/findings.tsv](../hunts/inputs/findings.tsv)):

| id | severity | finding | disposition | measured by |
| --- | --- | --- | --- | --- |
| IN-1 | high | missing values (`nan`, `NA`, `NaT`, `N/A`, `-`, `unknown`) are read as names and answer genes; 5 to 8 of 17 spellings per organism | [missing_values](rules.md#missing_values) | E `hunt_inputs` `hunt_IN-1`; E `e1` |
| IN-2 | high | mouse symbols typed in capitals or lower case reach another gene through an exact synonym: 145 of 76,955 | [single_case_symbols](rules.md#single_case_symbols) | E `hunt_IN-2`; E `e2_case` |
| IN-3 | medium | names made unique by single-cell tools (`SRF.1`, `MB.1`) reach other genes, also in the harmonizer | [made_unique_copies](rules.md#made_unique_copies) | E `hunt_IN-3`; E `e2_makeunique_*`, `e2_all_makeunique` |
| IN-4 | low | zero-width characters, byte-order marks, quotes, full-width characters and Unicode dashes make an identifier not found | [unicode_normalised](rules.md#unicode_normalised) | E `hunt_IN-4`; E `e2_zw_*`, `e2_fullwidth_*`, `e2_dash_*`, `e2_lookalike_*` |
| IN-5 | low | unversioned transcript and translation ids in another letter case are not found | [one_case_schemes](rules.md#one_case_schemes) | E `hunt_IN-5`; E `e2_txcase_*` |
| IN-6 | low | the graph's internal names (`synonym_id::TP53`, `<id>.Void`) are read as input | [no_synonym_target](rules.md#no_synonym_target) | E `hunt_IN-6`; E `e2_internal_*` |
| IN-7 | high | lists of names of a known release are placed at another release (human 5 of 9, mouse 3 of 6) or assembly | [inference_per_release](rules.md#inference_per_release) | E `hunt_IN-7` |
| IN-8 | high | a version of more than 4,300 digits raises a bare `ValueError` | [bounded_arguments](rules.md#bounded_arguments) | E `hunt_IN-8`; E `e2_longversion` |
| IN-9 | medium | bare numbers are read as old numeric Havana and RefSeq_mRNA ids: 134 of 1,500 mouse NCBI gene ids convert in the human graph | [bare_numbers](rules.md#bare_numbers) | E `hunt_IN-9`; E `e4_ncbi_foreign` |
| IN-10 | low | of names equal up to separators, a third spelling takes one silently | [one_case_schemes](rules.md#one_case_schemes) | E `e2_separators_*` (no line of the tie itself) |
| IN-11 | medium | stray `;`, `.` and quotes of Ensembl's tables are names: a symbol typed with them reaches another gene | [stray_punctuation](rules.md#stray_punctuation) | E `hunt_IN-11`; E `e2_punct_*`, `e2_all_punct` |

The counts are those of the findings table. Family E's set `hunt_inputs` holds the reproducers of IN-1 to IN-9 and IN-11 for human (23 lines), mouse (12) and pig (3), and `hunt_inputs_dec` the 23 of them whose expectation the decided rules change (`<families>/e/queries/<organism>/hunt_inputs*.manifest.json`); two reproducers have no line there because the truth's name tables do not give their expected gene as one gene (IN-4's `HLA–A`, written with an en dash, and IN-10's `AY172335_1`; the hunt's `C4B.1`, whose name several genes carry, is in neither), and IN-7's lines hold each name of a GTF once while the hunt's lists hold one name per gene, so E's line of human release 55 is answered right where the reproducer's is not. Their scores per tree are family E's results (`results/e/<tree>/<organism>/now_report.tsv` and `decided_report.tsv`, strata `hunt_IN-<n>`).

## History

**What was searched.** (1) *Mining*: every stable gene id present in a GTF of the window, with features read from Ensembl's own files (the history tables of release 115 and the GTFs of every release): gaps, retirement and return under the same id, version jumps and falls, events whose versions differ from the GTF's, moves of sequence or strand (within an assembly and across a change of assembly), splits, merges, splits and merges in one lineage, rejoins, chains of renames, fan-out, the connected components of the id graph of events (`<hunts>/history/mine/<organism>.features.tsv.gz`, `.patterns.json`). (2) *Selection*: per pattern a seeded sample of 50 ids and the 10 most complex (seed 20261118), and 20 members of each of the three largest components (`<hunts>/history/ids/<organism>.tsv`). (3) *Travel*: every id from every release where the GTF has it to every release of the window, as the stable id and as `id.version`; without `from_release`; and to the organism's name database at every release with a BioMart table; each call judged for `best` and `all` by the accuracy harness's judge with the travel matrix's truth (strict walk with living ids, broad walk with returned ids, copies apart, the locus rule with places mapped across assemblies), presence and versions read from the GTFs and the core gene tables (`<hunts>/history/runs/s1c/<organism>/`). (4) Truth-free checks on the same ids: *consistency* (one graph node asked from every release where it exists must answer the same; `<hunts>/history/runs/cons_s1/`) and *best choice* (idtrack's `best` against the genes with the largest product of Ensembl's scores along the strict walk; `<hunts>/history/runs/best_s1/`); and genome-wide scans for the patterns found (`<hunts>/history/runs/branch_s1/`, `<hunts>/history/mine/<organism>.seqmove.tsv`, `.living_branch.tsv`, `.gtf_vs_core.tsv`).

**How much** (the travel; `<hunts>/history/runs/s1c/<organism>/analysis.json`, field `totals`; calls of `best`, `all` and the name database together):

| organism | ids | calls judged | right | undetermined | wrong | ids with a wrong call |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pig | 618 | 4,274,817 | 4,005,926 | 267,325 | 1,566 | 3 |
| mouse | 648 | 7,147,229 | 6,924,417 | 160,054 | 62,758 | 46 |
| human | 649 | 5,656,131 | 5,393,661 | 187,502 | 74,968 | 102 |
| rhesus | 488 | 4,131,619 | 4,061,643 | 69,976 | 0 | 0 |
| cynomolgus | 329 | 585,343 | 581,287 | 4,056 | 0 | 0 |

No call raised and none fell back to another release of origin; the consistency check (5.6 million calls) found 6 differing cells, all of hh-3. What makes a lineage hard is not its size: the largest components (multi-copy RNA and gene families) travel without a wrong answer except through bridges; splits, merges, rejoins, chains, version jumps, events whose versions disagree with the GTF and changes of assembly travel right in every organism. Four patterns defeat idtrack: a gene that branches while it lives on (hh-1), a stable id that Ensembl's history carries to another place (hh-2), an id with no successor in the history, where bridges answer (hh-4), and an id brought back with an old version (hh-3).

**Findings** ([hunts/history/findings.tsv](../hunts/history/findings.tsv)):

| id | severity | finding | disposition | measured by |
| --- | --- | --- | --- | --- |
| hh-1 | medium | a branch that left an id while it lived on outranks the id's own line (pig, 7 genes, 1,133 calls) | [main_line_first](rules.md#main_line_first) | A* patterns `branch_while_living`, `hunt_history` |
| hh-2 | medium | an id that Ensembl moved to another chromosome is answered as the same gene, unmarked (human 609, mouse 210, pig 34 ids) | [place_moves_marked](rules.md#place_moves_marked) | A patterns `moved_seq`, `strand_flip`; A* patterns `moved`, `handoff_moved` |
| hh-3 | low | a returned id that reuses an old version answers forward and backward differently | [two_lives](rules.md#two_lives) | A patterns `gap`, `version_drop`; A* patterns `retired_and_back`, `version_fall` |
| hh-4 | medium | bridges where Ensembl's history has no successor (human 82, mouse 46, pig 3 of the hard ids) | [bridge_same_place](rules.md#bridge_same_place) | A pattern `retired_none`; A* pattern `hunt_history` |
| hh-5 | low | the truth's GTFs lack genes the core database has; sessions placed at their old release | truth error: `core_places`, `sessions_effective` | the harness's `core_places` sets; A* pattern `gtf_lacks`; A pattern `core_only` |
| hh-6 | low | the judge's copies apart before 82 and its place for ids Ensembl later moved | truth error: `copies_apart_early` | A* revision set `v1r2` |

hh-6 has two parts: the copies apart before release 82 are corrected in A* (`copies_apart_early`) and in the harness (human `v19`); the locus rule's reference place (the place at the release of origin, wrong for an id Ensembl later moved) has no correction yet, so such lines stay judged by the place of origin.

## Regions

**What was searched.** (1) *Region classes*: every gene id of the GTFs of every release tagged by what it is at any release: patch, scaffold, mitochondrial, Y, pseudoautosomal (X and Y), immunoglobulin and T-cell receptor segments, rRNA and small RNAs, families by name (HLA and its orthologues, KIR, olfactory receptors, KRTAP, histones, amylases, U6, Y RNA, 5S rRNA, 7SK, SRP, the rDNA unit, MIR, SNOR), names carried by several genes at one release, readthroughs, moves of sequence or strand, gaps (`<hunts>/regions/out/regions_<organism>.pkl`). (2) *Region travel*: up to 120 genes per tag, each id asked from up to 4 releases where it exists to every release of the graph (`best` for every id, `all` for every second, a name database for every sixth), the display name of every chosen gene at each BioMart release travelled to every later release; truth and judge as the harness's travel matrix (`<hunts>/regions/out/travel<k>/<organism>/`: the lines `lines.jsonl.gz`, the judged answers `shard_0_of_1.jsonl`, the wrong answers explained `wrong_paths.jsonl`). (3) *Probes* of named cases, each explained from its path (`<hunts>/regions/queue/p01.out` to `p30.out`). The hunt also grouped the harness's stored answers of `d4c2dda` by these tags; those tables read the judge as it was before its corrections and are replaced by the harness's current scoring ([the accuracy harness](harness.md)).

**How much** (lines of the region travel, `<hunts>/regions/out/travel<k>/<organism>/lines.jsonl.gz`):

| organism | first pass (every tag) | second pass (segmental duplications, gametologs) | third pass (genes not chosen before) |
| --- | ---: | ---: | ---: |
| human | 1,587,947 | 177,018 | 1,461,869 |
| mouse | 1,296,583 | 172,410 | 1,010,786 |
| pig | 887,476 | 53,120 | 664,168 |
| rhesus | 997,616 | 64,126 | 761,380 |
| cynomolgus | 336,194 | 14,860 | 271,322 |

The third pass, on genes no earlier pass had chosen, found no new class. Not testable with these graphs: mouse patches after release 102 (GRCm39 has none in the GTF), macaque patches (none) and macaque pseudoautosomal annotation (none).

**Findings** ([hunts/regions/findings.tsv](../hunts/regions/findings.tsv)):

| id | severity | finding | disposition | measured by |
| --- | --- | --- | --- | --- |
| R1 | high | mitochondrial gene ids that Ensembl gave to nuclear genes at an assembly change: `MT-ND1` of rhesus 97 converts to a nuclear gene (14 ids in three organisms) | [no_mito_nuclear_step](rules.md#no_mito_nuclear_step) | A pattern `mt`; A* patterns `hunt_history`, `moved` |
| R2 | high | a bridge through a synonym that another gene carries beats the gene's own name or NCBI id (`MT-CO2` to PTGS2; Y copies to their X gametologs) | [bridge_same_place](rules.md#bridge_same_place), [authority_symbol_first](rules.md#authority_symbol_first) | A pattern `retired_none`; A* pattern `hunt_history` |
| R3 | medium | GENCODE's `_PAR_Y` gene ids are not found (all 50 PAR genes) | [par_y_ids](rules.md#par_y_ids) | none (this reproducer only) |
| R4 | question | a name the X and Y copies of a PAR gene share ties the two genes from release 110 (37 of 39 names) | [par_shared_names](rules.md#par_shared_names) | A* pattern `par_y` (ids only) |
| R5 | low | UniProt entry names, Ensembl's display names at 48 and 49, are not found (2,447 names) | limitation `L-uniprot-entry-names` | none (this reproducer only) |
| R6 | medium | a patch or haplotype copy whose history ends goes to a paralog (KIR2DL4 to KIR2DL2, DEFB107B to DEFB107A) | [copies_placed_on_primary](rules.md#copies_placed_on_primary) | A* pattern `non_reference_copy_ended`; A pattern `patch` |
| R7 | known | an RFAM family label names every member (`U6` from 75 to 82: a tie of 1,279 genes) | [rfam_display_label](rules.md#rfam_display_label) | the harness's name kinds, cause `rfam_family_label` |
| R8 | medium | another gene's live synonym is read before a former name (`RNR1` to NR4A2) | [authority_symbol_first](rules.md#authority_symbol_first) | none (this reproducer only) |
| R9 | question | an LRG id answers itself, not the gene it describes | [lrg_answers_gene](rules.md#lrg_answers_gene) | none (this reproducer only) |
| J1 | judge | copies on patches, haplotypes and the Y PARs judged by the place of their own sequence | truth error: `exceptions`, `exceptions_span` | the harness's `exceptions` sets |
| J2 | judge | the GTF of cynomolgus 112 lacks 47 genes the core database has | truth error: `core_places` | the harness's `core_places` set; A* pattern `gtf_lacks` |
| J3 | judge | copies before 82 labelled primary, so a name expects them in the answer | truth error: `copies_apart_early` | A* revision set `v1r2`; in the harness the correction `copies_apart_early` (human `v19`, 429 lines; on `d4c2dda`'s answers 147 wrong -> right, none right -> wrong: [the corrections of the judge](harness.md#the-corrections-of-the-judge)) |
| J4 | judge | the human travel matrix's GRCh37-backward cell: all its wrong lines are queries on GRCh37 haplotypes and patches | truth error: `exceptions`, `exceptions_span` | the harness's travel matrix (human) |

R5 is not yet a row of the known-limitations table: no family asks names of releases 48 and 49, so its size is the hunt's count until one does. R3, R8 and R9 are measured by their reproducers only until a family publishes their lines.

## Reproducing

Every finding but J4 has a runnable reproducer: one table in `idtrack_reproducibility.hunts.reproduce`, one command. The configuration names the trees (an idtrack checkout, its Python, its graphs per organism) and the data root ([configuration](../README.md#configuration)); the reproducer restarts itself in the tree's Python with `PYTHONHASHSEED=0`, loads the graph of each organism its calls need (one at a time), asks the calls through the families' adapter (so `d4c2dda`'s api and `71690ab`'s are read alike), and prints for each call the expected answer with its truth's source, the answer recorded and the answer got.

```
# the table, without a graph
python -m idtrack_reproducibility.hunts.reproduce --list

# one finding on the reference tree (default --tree d4c2dda), or on another configured tree
python -m idtrack_reproducibility.hunts.reproduce inputs IN-2
python -m idtrack_reproducibility.hunts.reproduce regions R1 --tree 71690ab

# every finding of an angle, or of all angles; one organism only; the answers as JSON
python -m idtrack_reproducibility.hunts.reproduce history all --json hh.json
python -m idtrack_reproducibility.hunts.reproduce all all --organism mus_musculus

# as jobs through the configured runner: one per organism, each tree in turn
python -m idtrack_reproducibility.hunts.reproduce all all --submit --tree d4c2dda,71690ab
```

The exit status is 0 when every finding run still answers as recorded (the finding holds), 1 when an answer changed (the finding was fixed, or changed: the lines marked `CHANGED` say which field), 2 when the run cannot be made (no tree, graph or truth file). A call is a conversion (the ids of the view asked, the status, the match, the copies reported apart, the error class), a source inference on the gene names of a GTF of the truth store, or a fact of the truth store (a gene present in a GTF or in the core database's gene table); a recorded field `got_base` compares stable ids without versions where the hunt named none. With `--submit`, the answers go to `<hunts>/reproduce/<tree>/<angle>_<id>_<organism>.json` and the logs to `<hunts>/reproduce/logs/` (or `--out`); the memory of each job is the organism's from the configuration.

On 2026-11-19 every reproducer ran on `d4c2dda` and on `71690ab` (`all all --submit --tree d4c2dda,71690ab`; the answers in `<hunts>/reproduce/<tree>/all_all_<organism>.json`): all 29 findings hold as recorded on both trees, and the two trees give the same answer to every call. On idtrack's next version, a decided rule that changes an answer makes its finding's reproducer exit 1, naming the field it changed; a rule that only adds to the answer (hh-2's mark of a move, IN-9's database named) leaves the recorded fields as they are, and the family lines named in `measured_by` judge it. The reproducers of the judge's faults (hh-5, hh-6, J1 to J3) record idtrack's right answer, which no correction of the judge changes.

**Probes not ported.** The probes that counted a class read the graph's internals (node names per database and release, name resolution at a release, the engine's source counts, a node's active releases: inputs p09, p13, p14, p16, p19, p24, p25; regions p24's synonym count; history's best-choice and branch scans read candidates' scores and Ensembl's scores along all paths), and those internals differ between `d4c2dda` (a networkx graph under `API.track`) and `71690ab` (a typed store under `API.graph`), with no public call that both trees share; the families' adapter asks only the public calls. They are not ported: their counts stand as the hunts found them (the findings tables name the probe), their outputs are kept in the data root, and the classes they counted are measured from now on by the families' strata named in `measured_by`. As measurements, the history and regions hunts' travel runs are replaced by families A and A*, which travel the same patterns with the same judge on every tree, and the inputs hunt's determinism and metamorphic probes by family E (E2, E7, E8). J4 has no call reproducer: it is a cell of the harness's travel matrix, measured by the harness's scoring under the corrections `exceptions` and `exceptions_span`.
