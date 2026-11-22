# C: external identifiers

Family C asks idtrack to convert the identifiers of other databases: NCBI gene ids, RefSeq and UniProt accessions with and without their versions, RNA databases (RFAM, miRBase, RNAcentral), CCDS, Havana and Vega ids, clone-based names and the authorities' ids, each read at every release where Ensembl gives it a gene, without a release, and from old releases to the target; and every gene to every identifier database the graph holds at every release. Ensembl's cross-references of every release are the truth; each database's own file rides along. Its code is `idtrack_reproducibility.families.c`; its data are `<families>/c/` (`queries/<organism>/`, `runs/<tree>/<organism>/`, `truth_check/`); the design is in [the methods](../methods.md#c-external-identifiers). It reads the history as family A does ([A's truth](a.md#truth)) and runs on the common machinery ([F0](f0.md)).

## State

What is built:

- the fetcher of Ensembl's cross-references and of each database's own files (`fetch_c`), the truth readers (`truth_c`: the graph's scope, Ensembl's cross-references, the authorities), the generator of the parts c1 to c7 (`gen_c`), the check of the truth against Ensembl's other publications (`check_truth_c`), the authority measure (`authority_c`), the selection of the wrong answers to explain (`explain_qids`), the plugin of the common classifier (`causes_c`) and the corrections of the judge as revision sets (`revise_c`);
- the truth of all five organisms, each file in `SOURCES.tsv` (files of no single organism in `<truth>/_shared/SOURCES.tsv`): extracts of the cross-references of every release of each window from Ensembl's MySQL dumps (`core_xref/<release>.tsv.gz`: human, mouse and rhesus 48 to 116, pig 56 to 116, cynomolgus 91 to 116), Ensembl's FTP `tsv/` tables of 87 to 116, NCBI's `gene2ensembl`, `gene2refseq` and `gene_history`, RefSeq release 238 (catalog, accessions to gene ids, removed records), UniProt 2026_04 (`idmapping_selected`, `sec_ac`, `delac_sp`, `delac_tr`, Swiss-Prot isoforms), RFAM's family list, RNAcentral release 28 (its Ensembl and miRBase mappings), miRBase 23.1 (`hairpin.fa`, `mature.fa`, the genome files of human, mouse and rhesus; pig's and cynomolgus's do not answer), CCDS (human, mouse);
- the check of the core reading against BioMart and the FTP `tsv/` tables for cynomolgus, pig and rhesus (`<families>/c/truth_check/<organism>.tsv`);
- the published sets: c1 to c7 of cynomolgus, pig and rhesus; c1 and c2 of mouse and human, with their revision for the releases added to the truth (mouse 89, human 112) and the corrections of the judge ([the sets](#the-sets)).

What was run and scored, on the reference tree `d4c2dda` and its graphs:

| run folder | sets | `run.json`, `report.tsv` | causes | authority measure |
| --- | --- | --- | --- | --- |
| `runs/d4c2dda/macaca_fascicularis/` | `c1_v1` to `c7_v1` | yes | yes (`causes.tsv.gz`, `causes_summary.tsv`) | yes (`authority.tsv`) |
| `runs/d4c2dda/macaca_mulatta/` | `c1_v1` to `c7_v1` | yes | yes | yes |
| `runs/d4c2dda/sus_scrofa/` | `c1_v1` to `c7_v1` | yes | yes | yes |
| mouse, human | – | not run | – | – |

What is not built:

- c3 to c7 of mouse and human, and their runs; the check of the truth for mouse and human, and the family's panel ([Measurement design](../methods.md#measurement-design)); C's modest selection is built (`panels/<organism>/modest_v1`, [the modest set](../modest.md)).

Mouse's and human's `c1_v1` and `c2_v1` were published before the extracts of mouse release 89 and human release 112 existed (human 112's folder holds a GRCh37 database beside the release's own, mouse 89's dumps are `.txt.gz.bz2`); `gen_c --revise-from` published the lines a generation reading every release changes (2026-11-20; the explanations in words, `note` and `q2`, are no part of the comparison: `gen_c.WORDING`): mouse `c1_r1_new` 21,109 new lines at 89, `c2_r1` 4,777 (lines retired after 88 now read where last held, at 89) and `c2_r1_new` 309,356; human `c1_r1` 543 and `c1_r1_new` 27,648 at 112, `c2_r1` 8,368 (1,071 at 112; 7,297 lines of the sampled version strata left out, `excluded`, because the sample drawn from the population with 112 is another) and `c2_r1_new` 653,136.

The corrections of the judge (revision sets of the latest lines; the results below are judged with them, scored again on 2026-11-20): `rfam_label_accession` (`rfam_label_accession_v1`: cynomolgus 26 lines, rhesus 22, pig 13) and `other_database_broad` (`other_database_broad_v1`: pig 848, rhesus 1, mouse 58,577, human 16,522; `other_database_broad_v2`: rhesus 13, and mouse 26,684, the same correction made again (`revise_c --redo`) once it no longer counted the copies a line expects apart (a gene on a patch carrying the same value, rule `judge_copies_apart`) as other genes; cynomolgus none). In mouse most of its lines are NCBI gene ids that are also numeric accessions of `Vega_translation` of other genes at the release read. The judge following idtrack's decided rules (judge `c-rfam_label_rows`, published 2026-11-20): `rfam_label_rows` (`rfam_label_rows_v1`: pig 3,611 lines, rhesus 4,536, cynomolgus 1,898) for rules [`one_kind_per_database`](../rules.md#one_kind_per_database) and [`rfam_display_label`](../rules.md#rfam_display_label): idtrack's next version holds RFAM's labels under a name of their own (`RFAM label`, the graph record's `label_names`), so a label line asked `from_database="RFAM"` asks the accession rows and gets nothing (the 2,930 C answers of pig, rhesus and cynomolgus that went from right to wrong on `next` in the modest set); each label line says it asks the label rows (`from_database_rows: labels`) and the adapter asks the graph's name for them (a graph without the record, `d4c2dda`'s and `71690ab`'s, holds labels and accessions under `RFAM` and is asked so); a label written as the entry's accession stays with the accession rows (rule `one_kind_per_database`) and is not touched. And a label names only the gene whose display name it is: a gene whose GTF display name at the release read is another string leaves `expect`, a gene whose display name is the label while another database's cross-reference of it carries the same label goes to `broad` (the extract cannot tell which is the display cross-reference); where no gene is left the line is `excluded` (no label row is held there), or `open`. The truth changed for pig 4 lines and cynomolgus 28 (8 more `open`, 8 `excluded`); every other line only asks the label rows. `next`'s label lines of the modest set were asked again under the label rows (`modest --again c=rfam_label_rows_v1`): every judged answer right (`d4c2dda` and `71690ab`: wrong where they answer genes the label does not name, pig 2 and cynomolgus 8 answers in the modest set, cause `rule_decided_after_tree`).

The judge following the rules built in idtrack's next version after its first round (`revise_c --correction <name>`, judge `c-<name>`, published 2026-11-21; each function's docstring states it): `display_names_at_release` (`display_names_at_release_v1`: pig 687, mouse 3,088, human 2,903, rhesus 1,955; cynomolgus none; `display_names_at_release_v2`, the same correction on the values of the other sets, which the first publication did not read: pig 7,980, among them `c3_entry` 7,694 and `c4_mirbase_label` 277, rhesus 627, `c3_entry` 190 and `c4_mirbase_label` 415): a value the reference tree's configuration never held (rule `safe`: NCBI gene ids Ensembl gives a translation or a transcript at early releases, `c1_not_held`; UniProt entry names that early releases give genes as their display names, `c3_entry`; labels of RFAM and miRBase that are display names) that Ensembl gives a gene as a display name, or as an NCBI gene id of any of its forms, at some release is held there now, and the line is read as a held value (at the target if held there, else where it is held last, the genes walked; rule `set`); on the modest selection every revised line `r1` answers is right or undetermined (the undetermined: a bridge across an assembly change that no chain file maps). `no_mito_nuclear_step` (`no_mito_nuclear_step_v1`: rhesus 224, cynomolgus 103, human 2): a line whose value names, at the release read, a gene of the lineages of the mitochondrial ids that Ensembl gave to nuclear genes is walked again with the rule (family A's truth). `mirbase_genes` (`mirbase_genes_v1`: pig 70, rhesus 83; mouse and human none): the harness's revision of that name on C's lines (miRBase precursor names and accessions of two loci).

The judge following rule [`accession_versions`](../rules.md#accession_versions) as built in idtrack's next version after its first round (tree `r1b`; row I140; published 2026-11-21, judge `c-<name>`; the version strata only, which the corrections above leave out; `revise_c --correction <name>`, each function's docstring states it). A versioned line asked without a release is read as a line of the sets, the accession without its version or the version it reads, and takes that line's latest truth with every correction published for it (`read_as` names it); running a correction again publishes only what its earlier sets left wrong. `version_not_judged` (`version_not_judged_v1`: pig 15,151 lines, `c3_oldver` 7,226, `c3_ver_lower` 4,000, `c3_ver_higher_retired` 3,925; rhesus 7,947, `c3_ver_lower` 4,000 and `c3_ver_higher_retired` 3,947; cynomolgus 10,012, 2,013, 3,999 and 4,000; mouse and human have no c3 sets; run again, nothing more): of the pairs of consecutive releases at which one link (accession, Ensembl object with its version, database) holds in Ensembl's cross-references, the share in which `xref.version` differs (`truth_c.version_change_shares`, read from the extracts, not from idtrack): Swiss-Prot pig 69.2 %, rhesus 67.2, cynomolgus 62.4, mouse 91.6, human 99.7, TrEMBL 68.4, 63.3, 58.5, 91.1 and 95.3, RefSeq's six databases 0 to 3.8 % (human `RefSeq_mRNA_predicted`), the shares idtrack measures on its graphs; above one half (`truth_c.EDIT_COUNTER_SHARE`) the number counts the entry's edits, so a typed UniProt version expects the accession's reading without it (`version: not judged`; the lines read at a release, `c3_at_v`, already did). `versions_in_window` reads every other versioned line by the rule on the versions the graph's window had: a version it had where it last had it; a version above every version it had as the latest where the accession is attached at the graph's last release, and as an unknown version (nothing) where it is not (row I150); any other version, nothing; the window ends at the target (the generator counted Ensembl's release 116, one beyond the graphs, so a version written only at 116 made a version below it a gap); and a versioned string that Ensembl gives a gene as its display name (before release 64 Ensembl displays RefSeq's predicted mRNAs by their labels with the version, `RefSeq_dna_predicted`: rhesus `XM_001083595.1` at 55 to 62, while the graph's RefSeq versions begin at 64) is read as that name where that is later than the version's own reading (`read_as_name`; the version strata hold such strings: rhesus 1,447 lines, human 888, mouse 475, pig 261, cynomolgus none). Its sets, in order: `version_above_latest_v1` (`c2_ver_higher_retired` pig 3,999, mouse 3,873, human 3,705, rhesus 4,000, cynomolgus 3,464: the latest read also for an accession no longer attached at the target), `versions_in_window_v1` (cynomolgus 1,774, `c2_ver_gap` 1,692 and `c2_ver_higher_retired` 82; mouse 59; human 274; rhesus 292, the display readings), `versions_in_window_v2` (human 248: 241 lines the first set had given the generator's truth without the correction `par_shared_names`, put back, and 7 lines read as lines with that correction; pig 1; mouse 5, display readings), `versions_in_window_v3` (human 5, display readings), and the last ones, the rule with row I150's condition (published 2026-11-22: pig `versions_in_window_v2` 3,999, rhesus `versions_in_window_v2` 4,000, mouse `versions_in_window_v3` 3,932, cynomolgus `versions_in_window_v3` 5,238, human `versions_in_window_v4` 3,733): every line of `version_above_latest_v1` is an accession not attached at the target and expects nothing again, as do the gap lines the first sets had read as the latest (cynomolgus 1,692 and 82, mouse 56 and 3, human 26 and 2). On 100 sampled lines of pig's `version_above_latest_v1`, 60 per stratum of its `version_not_judged_v1`, 182 of cynomolgus's `versions_in_window_v1` and the 292 of rhesus's and 5 of mouse's and human's display readings, `r1b` read every line at the release the revision says (`read_at`); `r1b` reads the latest also where the accession is not attached at the target, so it answers the walked genes on the lines of row I150's class (a fault of `r1b`, row I150).

## What it measures

The rate of right answers, right / (right + wrong + error) (undetermined answers and timeouts apart), per stratum and organism, with Wilson intervals and a cluster bootstrap over the identifier without its version (c6: the gene); every line is scored under `best` and under `all` (`strategy: both`):

- **c1 to c5, at a release and by default**: whether an identifier read at a release where the graph holds it is answered with the genes Ensembl's cross-references give it there; read without a release (rule [`start_release`](../rules.md#start_release)); read where it was last held when the graph no longer holds it, and carried to the target 115.
- **c2, c3, versions**: an accession written with the version Ensembl gave it, with an older version it had in the window, and with versions it never had (lower, in a gap, higher; the database's newer version).
- **Databases the graph does not hold** (UniProt isoforms, entry names and secondary accessions, miRBase, RNAcentral, the authorities' ids, values of a held database at releases or on objects the configuration leaves out): asked all the same with the rule `safe` (nothing or the identifier's gene is right, another gene wrong); how often idtrack recovers the gene is the coverage.
- **c6, genes to databases**: every gene of every table release converted to every identifier database the graph holds there (`final_database`), compared exactly (`cmp: exact`).
- **c7, old identifiers to now**: every value read at an old release and carried to 115.

Where a database's own file gives an identifier other genes than Ensembl's cross-references (or none), the difference is measured beside idtrack (`authority.tsv`), never as a wrong answer of idtrack (principle 2 of [the methods](../methods.md#principles)).

## The sets

One set per part and organism, `<part>_v1`, published once (F0's `lines.publish`); a change of truth is a revision set. Releases read ("table releases"): those of the window with a BioMart table or an FTP `tsv/` table (87 on) and the target 115. qid `<kind>|<identifier>|<from or ->|<to>[|<final database>][|db=<from database>]`.

| part | strata | rule | truth |
| --- | --- | --- | --- |
| c1 NCBI gene ids | `c1_at` (every id at every table release where the graph holds it), `c1_default`, `c1_retired` (held before, not at the target: read where last held, walked), `c1_not_held` (Ensembl's ids the configuration never holds), `c1_discontinued` (NCBI's `gene_history` ids of the taxon that no release holds: nothing) | `set`; `safe` for `c1_not_held` | Ensembl's `EntrezGene` cross-references; `auth`: `gene2ensembl` (a discontinued id: its replacement's genes) |
| c2 RefSeq | `c2_at`, `c2_at_v` (with the version Ensembl wrote at R), `c2_default`, `c2_default_v`, `c2_retired`, `c2_oldver` (a version it had in the window, not at the target: read where it last had it), `c2_not_held`; versions never had, 4,000 sampled per class: `c2_ver_lower`, `c2_ver_gap`, `c2_ver_higher_retired` (nothing), `c2_ver_higher`, `c2_ver_catalog` (RefSeq's catalog version newer than Ensembl's: the gene) | `set` | versions from the core `xref.version`; `auth`: `gene2refseq` through `gene2ensembl` |
| c3 UniProt | `c3_at`, `c3_at_v`, `c3_default`, `c3_default_v`, `c3_retired`, `c3_oldver`, `c3_ver_lower`, `c3_ver_higher`, `c3_ver_higher_retired`; `c3_isoform`, `c3_isoform_canonical`, `c3_entry`, `c3_secondary`; `c3_secondary_held` (a secondary accession the graph held as a primary one), `c3_deleted` (accessions the graph held that UniProt deleted) | `set`; `safe` for isoforms, entry names, secondary accessions | Swiss-Prot and TrEMBL cross-references (Ensembl writes the entry version as the version); isoforms: Ensembl's `Uniprot_isoform` and Swiss-Prot's isoforms; `auth`: `idmapping_selected` |
| c4 RNA databases | `c4_rfam_acc_<at\|default\|retired\|not_held>`, `c4_rfam_label_<...>` (with `from_database` RFAM), `c4_mirbase_acc`, `c4_mirbase_label`, `c4_mature_acc`, `c4_mature_name`, `c4_rnacentral` | `set` (RFAM); `safe` (miRBase, RNAcentral) | RFAM, miRBase and RNAcentral cross-references; mature miRNAs through their precursors (miRBase's genome files); `auth`: RFAM's family list, RNAcentral's mappings |
| c5 other gene databases | CCDS (accession and versioned label: `c5_ccds_*`, `c5_ccds_label`), Havana gene, transcript, translation (`c5_ottg_*`, `c5_ottt_*`, `c5_ottp_*`), Vega gene, clone-based names (`c5_clone_*`), each at R, by default, retired; `c5_hgnc_id`, `c5_mgi_id`, `c5_vgnc_id` | `set`; `safe` for the authorities' ids | Ensembl's cross-references; `auth`: CCDS's sequence table, VGNC's gene sets |
| c6 genes to databases | `c6:<database>`: every gene of the GTF of every table release x every identifier database the graph holds there (the name databases are B8's) | `set`, `cmp: exact`, `spellings` (RFAM accession and label, CCDS accession and label are one value) | the values Ensembl gives the gene at R; rule [`one_kind_per_database`](../rules.md#one_kind_per_database): `EntrezGene` the NCBI ids, `NCBI gene (formerly Entrezgene)` the NCBI symbols, at every release with NCBI ids |
| c7 old identifiers to now | `c7_ncbi`, `c7_refseq` (every value at every table release before 115); `c7_uniprot`, `c7_ccds`, `c7_ott*`, `c7_vega_gene`, `c7_clone`, `c7_rfam`, `c7_rfam_label` (50,000 sampled (value, release) pairs per organism) | `set` | the walks of the value's genes |

Lines per set and organism (manifests `<families>/c/queries/<organism>/<set>.manifest.json`):

| set | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `c1_v1` | 621,490 | 897,178 | 937,388 | 622,485 | 400,893 |
| `c2_v1` | 7,018,312 | 10,994,278 | 15,841,952 | 7,112,048 | 3,814,412 |
| `c3_v1` | 3,336,924 | – | – | 2,609,067 | 2,381,019 |
| `c4_v1` | 26,707 | – | – | 31,398 | 24,215 |
| `c5_v1` | 144,871 | – | – | 15,373 | 751 |
| `c6_v1` | 11,196,372 | – | – | 11,598,840 | 8,523,060 |
| `c7_v1` | 3,774,963 | – | – | 3,790,394 | 2,091,539 |

## Truth

Nothing of idtrack (`truth_c`):

- **Ensembl's cross-references of every release** (`Xrefs`): the extract of each release's MySQL dumps (`external_db`, `xref`, `object_xref`, `gene`, `transcript`, `translation`, and the `*_stable_id` tables before 65; read, not kept): one row per Ensembl object and cross-reference of C's databases, with its gene, `db_name`, `db_display_name`, accession, label, version and info type. A value of a database at release R names the genes whose gene, transcript or translation carries it at R.
- **The graph's scope** (`Scope`): which of the graph's databases the reference tree's externals configuration includes at which release (by the assembly of the release's GTF) and from which object form, read from the configured `d4c2dda` checkout (`idtrack/default_config/<organism>_externals_modified.yml`), with `d4c2dda`'s naming (an accession under its `db_name`, a label under its `db_display_name`). It says what a graph holds, not what is right: a value Ensembl has that the configuration leaves out is asked all the same (`*_not_held`, rule `safe`).
- **The history to the target**: family A's `GenA` (A*'s strict walk with living ids is `expect`, rule [`id_lives_on`](../rules.md#id_lives_on); the broad walk with ids brought back is `broad`, rule [`resurrection`](../rules.md#resurrection); copies on patches apart, rule [`judge_copies_apart`](../rules.md#judge_copies_apart); the locus rule with places mapped across assemblies). A value of several genes: the union of their walks; `open` when none of them has a successor, `partial` counted when some have none.
- **Each database's own file** (`Authority`, beside the truth, never the judge): NCBI's `gene2ensembl`, `gene2refseq`, `gene_history`; RefSeq's catalog; UniProt's `idmapping_selected`, secondary and deleted accessions, isoforms; RFAM's family list; RNAcentral's mappings; miRBase's genome and sequence files (mirbase.org serves them wrapped in HTML; they are kept as served and unwrapped by `truth_c.mirbase_text`); CCDS's sequence table; VGNC's gene sets. A line carries the database's genes as `auth` (`none` where the database attaches no gene).
- **The check** (`check_truth_c`): BioMart agrees with the core reading at every archive release and database read, except release 75's Swiss-Prot table, which holds entry names where the core holds accessions (pig 566, rhesus 171 values); the FTP `tsv/` tables agree with it except at one release per organism (cynomolgus and rhesus 112, pig 105), where they lack values the core holds and pig's holds 217 RefSeq peptides the core does not (`<families>/c/truth_check/<organism>.tsv`). The core is what the graph reads, so it is the truth.

## Rules

The rules of idtrack the expectations rest on: [`start_release`](../rules.md#start_release) (lines without a release: read at the target, else where last held), [`accession_version_cut`](../rules.md#accession_version_cut) (an accession's version is cut), [`version_newer_than_graph`](../rules.md#version_newer_than_graph) and [`unknown_version`](../rules.md#unknown_version), [`from_database`](../rules.md#from_database) (RFAM labels asked with `from_database`), [`tie_in_final_database`](../rules.md#tie_in_final_database) (c6), [`copies_apart`](../rules.md#copies_apart), [`external_bridge`](../rules.md#external_bridge); the judge's [`id_lives_on`](../rules.md#id_lives_on), [`resurrection`](../rules.md#resurrection), [`judge_copies_apart`](../rules.md#judge_copies_apart) and [`judge_corrections`](../rules.md#judge_corrections) (`spellings_once`: one value with two spellings). The decided rules C measures before they land: [`one_kind_per_database`](../rules.md#one_kind_per_database) (c6's NCBI symbols under one name at every release), [`mirbase_genes`](../rules.md#mirbase_genes), [`rfam_display_label`](../rules.md#rfam_display_label) and [`bridge_same_place`](../rules.md#bridge_same_place).

Accession versions are judged by rule [`accession_versions`](../rules.md#accession_versions) (the answer to [Q2](../methods.md#open-methodological-choices); built in `r1b`, the condition on an accession no longer attached at the target after it, row I150): a version it had is read where it had it (at the target if held there with that version, else at the latest release with it); a version above every version it had is read as the latest where the accession is attached at the target, as is the RefSeq catalog's newer version, and is an unknown version (nothing) where it is not; a version below the latest that it never had (lower, in a gap) is answered with nothing (as rule `unknown_version` answers an Ensembl version); and where the database's number counts its records' edits (UniProt's entry version, measured from Ensembl's records by `truth_c.version_change_shares`), a typed version is not judged and the line expects the accession's reading without it. The generator wrote the first and third parts; the revisions `versions_in_window` and `version_not_judged` bring the rest ([the sets](#the-sets)). `d4c2dda` to `r1` cut every accession version (rule `accession_version_cut`), so their answers to the versions never had are wrong under this rule.

Scoring is F0's ([scoring and outcomes](f0.md#scoring-and-outcomes)).

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)); the scope of the graph is read from the configured `d4c2dda` checkout.

```
# the truth (once; one subcommand per host, one at a time)
python -m idtrack_reproducibility.families.c.fetch_c ncbi|ebi|mirbase
python -m idtrack_reproducibility.families.c.fetch_c ensembl [--releases 87-116] [--kinds tsv,xref] [--organisms <o>,...]
python -m idtrack_reproducibility.families.c.check_truth_c <o> [<o> ...]
# the sets (published once; one part per call, a job of the configured runner for the large parts)
python -m idtrack_reproducibility.families.c.gen_c <o> --part c1|c2|c3|c4|c5|c6|c7 --set <part>_v1 [--dry-run [--show N]]
# a set published before its truth gained releases: the changed lines as a revision set, the new qids beside it
python -m idtrack_reproducibility.families.c.gen_c <o> --part <part> --set <revision set> --revise-from <part>_v1
# corrections of the judge, as revision sets of the latest lines
python -m idtrack_reproducibility.families.c.revise_c <o> --correction rfam_label_accession|other_database_broad|rfam_label_rows --set <revision set> [--dry-run]
# a correction fixed after its set was published: its lines corrected anew from the lines before it, the differing ones published
python -m idtrack_reproducibility.families.c.revise_c <o> --correction <name> --set <revision set> --redo <its published set>
# a tree's run and its scoring (F0's runner)
python -m idtrack_reproducibility.families.runner --family c --organism <o> --clone <source> --graph <graph> \
    --out <families>/c/runs/<tree>/<o> --submit <n> --score --score-parts <P> [--tag c_<tree>_<o>]
# causes: at most N wrong answers per unit and stratum explained (not c6 and the version strata), then classified
python -m idtrack_reproducibility.families.c.explain_qids <families>/c/runs/<tree>/<o> [N]
<python of the tree> -m idtrack_reproducibility.families.explain --family c --organism <o> --clone <source> --graph <graph> \
    --scored <families>/c/runs/<tree>/<o> --qids <families>/c/runs/<tree>/<o>/explain_qids.tsv
python -m idtrack_reproducibility.families.causes --organism <o> --scored <families>/c/runs/<tree>/<o> \
    --plugin src/idtrack_reproducibility/families/c/causes_c.py
# the authority measure
python -m idtrack_reproducibility.families.c.authority_c --organism <o> --run <families>/c/runs/<tree>/<o>
```

A rescoring after a correction of the judge reads the stored answers (`runner ... --submit 0 --score`); nothing is asked again.

## Results

Runs `<families>/c/runs/d4c2dda/<organism>/` of cynomolgus, pig and rhesus (the reference graphs, every set of version 1); result files `results/c/d4c2dda/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv` and `authority.tsv`. No answer timed out. The results are judged with the corrections of the judge (scored again on 2026-11-20 from the stored answers); mouse and human are not run.

Every stratum read at a release and by default is right on every answer in the three organisms (`c1_at`, `c1_default`, `c2_at_v`, `c2_default`, `c2_default_v`, `c3_at`, `c3_at_v`, `c3_default`, `c3_default_v`, the RFAM accessions and labels at R and by default, pig's Havana ids), and pig's `c2_at` and clone-based names read at a release, whose 42 and 718 answers once wrong are undetermined under the correction `other_database_broad` (the string is also another database's value of the answered gene); so are NCBI's discontinued ids (nothing), the versions above every version of a held accession and the catalog's newer versions (the gene), and every database the graph does not hold but miRBase's names (`safe`: nothing or the gene). Rate, cluster-bootstrap interval and wrong answers of the other strata (row `*` of `report.tsv`; – : no such stratum): Since `rfam_label_rows_v1` (rule `rfam_display_label`), `d4c2dda`'s RFAM label lines at R are wrong where it answers genes whose display name is another string: pig 8 answers (4 labels), cynomolgus 45 (24), cause `rule_decided_after_tree`; cynomolgus's 8 lines whose label names no gene at R are left out (`excluded`) and 27 answers are undetermined (`open`).

| stratum | cynomolgus | pig | rhesus |
| --- | ---: | ---: | ---: |
| `c1_retired` | 1.0000 (1.0000–1.0000), 0 | 0.9980 (0.9971–0.9989), 35 | 1.0000 (1.0000–1.0000), 0 |
| `c2_retired` | 1.0000 (1.0000–1.0000), 0 | 0.9987 (0.9984–0.9990), 150 | 1.0000 (0.9999–1.0000), 4 |
| `c2_oldver` | 0.9989 (0.9987–0.9991), 278 | 0.9885 (0.9876–0.9893), 2,046 | 0.9935 (0.9930–0.9940), 1,464 |
| `c2_ver_lower` | 0.0737 (0.0655–0.0820), 7,410 | 0.0280 (0.0230–0.0335), 7,776 | 0.0440 (0.0377–0.0505), 7,648 |
| `c2_ver_gap` | 0.0030 (0.0015–0.0048), 7,910 | 0.0032 (0.0015–0.0050), 7,974 | 0.0016 (0.0000–0.0048), 1,246 |
| `c2_ver_higher_retired` | 0.1557 (0.1442–0.1673), 6,754 | 0.0550 (0.0480–0.0617), 7,560 | 0.0462 (0.0400–0.0527), 7,630 |
| `c3_retired` | 1.0000 (1.0000–1.0000), 0 | 0.9943 (0.9936–0.9950), 510 | 1.0000 (1.0000–1.0000), 0 |
| `c3_deleted` | 1.0000 (1.0000–1.0000), 0 | 0.9957 (0.9952–0.9962), 491 | 1.0000 (1.0000–1.0000), 0 |
| `c3_secondary_held` | 1.0000 (1.0000–1.0000), 0 | 0.9995 (0.9990–0.9998), 10 | 1.0000 (1.0000–1.0000), 0 |
| `c3_oldver` | 0.9983 (0.9981–0.9984), 1,312 | 0.9943 (0.9940–0.9946), 6,385 | 1.0000 (1.0000–1.0000), 0 |
| `c3_ver_lower` | 0.0772 (0.0688–0.0855), 7,382 | 0.0190 (0.0150–0.0230), 7,848 | 0.0000 (0.0000–0.0000), 8,000 |
| `c3_ver_higher_retired` | 0.2318 (0.2190–0.2450), 6,146 | 0.0838 (0.0752–0.0927), 7,330 | 0.0000 (0.0000–0.0000), 7,894 |
| `c4_mirbase_label` | – | 0.9905 (0.9778–1.0000), 6 | 0.9387 (0.9158–0.9595), 56 |
| `c4_rfam_label_retired` | 0.9993 (0.9980–1.0000), 1 | 0.9986 (0.9965–1.0000), 2 | 0.9987 (0.9968–1.0000), 2 |
| `c5_clone_retired` | – | 0.9990 (0.9986–0.9993), 37 | – |
| `c6:EntrezGene` | 1.0000 (1.0000–1.0000), 0 | 0.9341 (0.9335–0.9349), 127,392 | 0.9546 (0.9541–0.9550), 95,808 |
| `c6:NCBI gene (formerly Entrezgene)` | 0.8187 (0.8169–0.8204), 309,084 | 0.6563 (0.6542–0.6583), 361,076 | 0.5739 (0.5716–0.5762), 65,068 |
| `c7_ncbi` | 1.0000 (1.0000–1.0000), 0 | 0.9951 (0.9945–0.9957), 5,107 | 1.0000 (1.0000–1.0000), 0 |
| `c7_refseq` | 1.0000 (1.0000–1.0000), 0 | 0.9971 (0.9969–0.9973), 18,056 | 1.0000 (1.0000–1.0000), 0 |
| `c7_uniprot` | 1.0000 (1.0000–1.0000), 0 | 0.9974 (0.9969–0.9979), 236 | 1.0000 (1.0000–1.0000), 0 |

The only errors are in `c6:NCBI gene (formerly Entrezgene)`: pig 303,864 and rhesus 833,476 answers raise a `ValueError` naming the database (`errors.tsv` of the run), where `d4c2dda` does not hold that name at the release asked. Every other gene-to-database conversion of c6 (RefSeq, UniProt, RFAM, Havana, clone-based names) is right on every answer.

### Causes

Every wrong answer of the three organisms classified by F0's classifier with C's plugin (the strata of c6 and of the versions need no explained path; of the others at most 5 answers per unit and stratum are explained; the classifier meets the rows release by release, so a release's cross-references are read once: pig's 560,027 rows take minutes); answers summed over the strata of `results/c/d4c2dda/<organism>/causes_summary.tsv`:

| cause | class | cynomolgus | pig | rhesus |
| --- | --- | ---: | ---: | ---: |
| `ncbi_label_name_by_release`: a gene converted to `NCBI gene (formerly Entrezgene)` before 100, where the graph files the NCBI symbols under `EntrezGene` (to 87) or `NCBI gene` (88 to 99) (rule `one_kind_per_database`, decided) | a | 309,084 | 361,076 | 65,068 |
| `entrezgene_label_rows`: `EntrezGene` up to 87 returns the NCBI labels beside the numbers (rule `one_kind_per_database`) | a | – | 127,392 | 95,808 |
| `accession_version_never_had`: an accession version never had in the window is answered with the accession's gene ([Q2](../methods.md#open-methodological-choices)) | a | 35,602 | 38,488 | 32,418 |
| `accession_old_version_read_at_target`: a version had at an older release is read at the target without its version | a | 1,408 | 6,634 | 1,214 |
| `no_successor_answer_elsewhere`: no successor in Ensembl's history and the answer lies elsewhere, without an explained path (the class of `external_bridge`) | b | – | 14,903 | – |
| `external_bridge`: no successor in Ensembl's history, answered through an external id with a gene elsewhere | b | – | 8,827 | 17 |
| `value_of_another_database`: the string is also another database's value of the answered gene at the release read | e | – | 1,377 | 4 |
| `mirbase_not_in_graph`: miRBase names are not in the graph (rule `mirbase_genes`, decided) | c | – | 6 | 56 |
| `synonym_spelling_merged`: the string is carried in two letter cases | c | – | 2 | – |
| `rfam_label_is_accession`: the RFAM label line's string is the entry's accession (the truth's; what the correction `rfam_label_accession` leaves) | d | 1 | 2 | 2 |
| `unexplained` (old versions of `c2_oldver` and `c3_oldver`, and pig's 14 `c5_clone_retired`, that no test explains) | – | 182 | 1,320 | 233 |

**Two rules give almost every wrong answer.** The NCBI symbols under one name at every release (rule `one_kind_per_database`, decided for idtrack's next version) and the versions of accessions (Q2) give every classified wrong answer of cynomolgus but one RFAM label, of rhesus but 56 miRBase names, 17 bridges, 4 values of another database and 2 RFAM labels; in pig 533,590 of 560,027, the rest mostly where Ensembl's history ends: carried to a gene elsewhere (23,730, class b: `c7_refseq`, `c7_ncbi`), and values of another database (1,377). The old versions that no test explains (cynomolgus 182, pig 1,306, rhesus 233) and pig's 14 clone-based names are read by hand.

### The authority beside Ensembl

Where a database's own file and Ensembl's cross-references give an identifier at 115 other genes (`truths_differ`, `best` rows of `results/c/d4c2dda/<organism>/authority.tsv`; mostly the database attaching no gene, `auth_only_none`):

| stratum | cynomolgus | pig | rhesus |
| --- | ---: | ---: | ---: |
| `c1_default` (NCBI's `gene2ensembl`) | 7,158 of 13,563 (0.5278; no gene 6,554) | 16,477 of 16,477 (1.0000; no gene 1,747) | 2,695 of 21,446 (0.1257; no gene 2,582) |
| `c2_default` (`gene2refseq` through `gene2ensembl`) | 18,724 of 41,811 (0.4478) | 115,991 of 115,991 (1.0000) | 51,416 of 129,010 (0.3985) |
| `c3_default` (UniProt's `idmapping_selected`) | 2,786 of 51,341 (0.0543) | 27,281 of 70,701 (0.3859) | 6,924 of 48,212 (0.1436) |
| `c3_deleted` (accessions UniProt deleted that Ensembl's tables still give a gene) | 11,662 of 16,828 (0.6930) | 49,989 of 62,652 (0.7979) | 19,434 of 30,557 (0.6360) |

NCBI's `gene2ensembl` of pig names genes of another Ensembl annotation of the species (ids `ENSSSCG00220...`, not the reference assembly's `ENSSSCG000...`), so for pig no NCBI gene id and no RefSeq accession is given the reference gene by NCBI's own files: the pig rows of c1 and c2 measure that difference, not idtrack.

idtrack follows Ensembl, as it should; these shares are the size of the limitation "Ensembl against the authority".

## Known limitations

C's rows of [the known-limitations table](../../results/limitations/limitations.tsv) are generated by F0's builder with C's map (`families/c/limitations_map_c.tsv`) from the causes of the three organisms. Answers per organism:

| limitation | what | user action | until | pig | rhesus | cynomolgus |
| --- | --- | --- | --- | ---: | ---: | ---: |
| `L-ncbi-label-names` | Ensembl files the NCBI symbols under `EntrezGene` up to 87 and `NCBI gene` from 88 to 99, and under `NCBI gene (formerly Entrezgene)` from 100 on; a symbol asked under its name of 100 at an earlier release is not found | ask for the database by its name at the release asked (`final_databases(release)`) | decided: rule `one_kind_per_database` | 361,076 | 65,068 | 309,084 |
| `L-ncbi-labels` | up to release 87 Ensembl's EntrezGene holds NCBI labels beside the numbers; a gene converted to EntrezGene at those releases also returns its label | ask for `NCBI gene (formerly Entrezgene)` or filter numeric ids | decided: rule `one_kind_per_database` | 127,392 | 95,808 | – |
| `L-accession-versions` | every version of a RefSeq or UniProt accession is cut, so a version the accession never had in the window answers the accession's gene | drop accession versions you cannot vouch for, or check them against the database | decided: rule `accession_versions` | 45,122 | 33,632 | 37,010 |
| `L-bridges` | an id Ensembl's history ends is carried through an external id to a gene that does not overlap it | `allow_bridges=False`; check the answer's place (explain) | decided: rule `bridge_same_place` | 23,730 | 17 | – |
| `L-other-database-values` | the string is also a value or name of the returned gene in another database at the release read; without `from_database` it names both | pass `from_database` | none (Ensembl's data; undeterminable) | 1,377 | 4 | – |
| `L-unexplained` | answers no test of the classifier explains (docs/families/c.md, "Causes") | - | - | 1,320 | 233 | 182 |
| `L-mirbase` | miRNA genes whose display name is a miRBase name are not found by that name | use the Ensembl id or the HGNC/MGI symbol of the miRNA gene | decided: rule `mirbase_genes` | 6 | 56 | – |
| `L-spellings` | two genes with one string in two letter cases | type the exact spelling | decided: rule `spelling_per_edge` | 2 | – | – |
| `L-truth-rfam-labels` | at some releases Ensembl writes an RFAM entry's accession as its label: the label line expects the accession's string | none needed | the correction `rfam_label_accession` (published; these are the lines it leaves) | 2 | 2 | 1 |

`L-accession-versions` is closed in idtrack's next version but for one class (tree `r1b`, rule [`accession_versions`](../rules.md#accession_versions)): on the [modest set](../modest.md#the-tree-r1b)'s selection of the version strata (`c2_*_v`, `c2_oldver`, `c2_ver_*` and the same of c3; 500 sampled lines per stratum, both views), judged with the revisions `version_not_judged` and `versions_in_window`, the answers of `r1b` wrong by the reading of a version are all of row I150's class, a version above every version of an accession no longer attached at the graph's last release read as the latest (pig 904 of 14,000 answers, rhesus 924 of 14,000, cynomolgus 1,244 of 14,000, mouse 994 of 8,000, human 36 of 5,440; `r1`: 2,855, 2,793, 2,744, 2,962 and 1,186 wrong); mouse's 6 others are bridges where Ensembl's history gives the accession's genes no successor (`L-bridges`).

What the results show beyond the classified causes, and what the truth cannot judge:

| limitation | what | until | size |
| --- | --- | --- | --- |
| NCBI symbols by release | converting to `NCBI gene (formerly Entrezgene)` before 100 gives nothing, or an error where `d4c2dda` does not hold the name at the release; `EntrezGene` up to 87 returns the labels beside the numbers | the decided rule `one_kind_per_database` | cause rows above; errors pig 303,864, rhesus 833,476 (`report.tsv`) |
| accession versions | every accession version is cut: a version never had is answered with the gene, an old version is read at the target | the decided rule `accession_versions` (the answer to Q2) | `c2_ver_*`, `c3_ver_*`, `c2_oldver`, `c3_oldver` above |
| bridges | an identifier whose genes Ensembl retired without successor is answered with a gene elsewhere | the decided rule `bridge_same_place` | pig 23,730, rhesus 17 |
| databases not in the graph | UniProt isoforms, entry names and secondary accessions, miRBase, RNAcentral, the authorities' ids are answered with nothing (right under `safe`; the gene is recovered nowhere) | the decided rule `mirbase_genes` for miRBase | the `safe` strata above |
| RFAM labels written as accessions (the truth's) | Ensembl writes an entry's accession as its label at some releases; the label line expected only the label rows' genes | the correction `rfam_label_accession` (published; the entry's genes, at the release read and at the target, go to `broad`): of cynomolgus's 38, pig's 20 and rhesus's 36 wrong answers, 37, 18 and 34 became undetermined | cynomolgus 1, pig 2, rhesus 2 left |
| a value of several databases (the truth's) | a string that another database or the GTF gives other genes at the release read (pig's clone-based names such as `ABCA1`, RefSeq values written `LOC...`, mouse's NCBI ids that are numeric `Vega_translation` accessions): without `from_database` it names both | the correction `other_database_broad` (published: the other genes the string names at the release read, or at the target for a line asked without a release and without a version, go to `broad`, or to `expect` of a `safe` line; the copies a line expects apart are not other genes) | pig's `c5_clone_at` 718, `c2_at` 42, `c5_clone_retired` 46, `c7_clone` 3, `c2_retired` 5 and `c7_refseq` 12 wrong answers became undetermined, rhesus's `c2_retired` 4; classified as `value_of_another_database` still: pig 1,377, rhesus 4 |
| the authority | where a database's own file and Ensembl disagree, the disagreement is Ensembl's | none (Ensembl's data) | [the authority beside Ensembl](#the-authority-beside-ensembl) |

## Reproducing from scratch

From an empty data root, with the configuration naming the tree `d4c2dda`, its checkout and its graphs:

1. **Fetch the truth**: the accuracy harness's fetches (GTFs, BioMart's archive tables, the core tables of the history, VGNC's gene sets; [the accuracy harness](../harness.md#truth)) and A*'s ([A*'s steps](astar.md#reproducing-from-scratch)); then `fetch_c ncbi`, `fetch_c ebi`, `fetch_c mirbase`, `fetch_c ensembl` (one host at a time), and `check_truth_c` for every organism.
2. **Publish the sets**: `gen_c <o> --part <part> --set <part>_v1` for c1 to c7 and every organism. A rebuild with the current code reads every release of every window; for mouse and human it differs from the published `c1_v1` and `c2_v1` at releases 89 and 112, which are published as a revision set (`--revise-from`), never by rewriting a set; then the corrections of the judge (`revise_c`).
3. **Run and score** every set on `d4c2dda` for every organism (F0's runner).
4. **Classify**: `explain_qids`, `explain` on the tree, then `causes` with C's plugin; `authority_c`.
5. **Results**: `idtrack-repro results` writes `results/c/d4c2dda/<organism>/` once the runs are in `results.toml`.
