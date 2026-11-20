# The rules the benchmarks refer to

The benchmarks judge idtrack against independent truth, and some answers depend on a rule that idtrack follows (how a string is read, which candidate is the best answer) or on a rule of the judge (how a truth is read where Ensembl's data are ambiguous). Every such rule has a name; the code, the query sets' notes, the causes of wrong answers and the known-limitations table cite rules by these names.

Status of each idtrack rule:

- **d4c2dda**: idtrack behaves so from commit `d4c2dda` (the reference tree of the benchmarks) on;
- **71690ab**: from commit `71690ab` on (`d4c2dda`'s answers with a typed api; idtrack's main branch);
- **decided**: decided for idtrack's next version, not in `71690ab`; the families score such lines under two judges where the rule changes the expected answer (`now`: the documented behaviour of `71690ab`; `decided`: the rule), so that the change is measured when it lands.

Rules of the judge (the truth's reading) are marked **judge**.

## How a query is read

### start_release
*d4c2dda.* With `from_release` given, the identifier is read at that release; an identifier that is not active there is an error that names the releases where it is active. Without it: if the identifier is active at the target release it is read there, otherwise at the latest release where it is active. A versionless Ensembl id resolves to the version active at the starting release. One rule for Ensembl ids and external identifiers.

### from_release_fallback
*d4c2dda.* An external identifier given with a `from_release` at which the graph does not hold it is converted from the release at which the graph holds it, and the result names the release used and carries a flag; nothing is silent.

### spelling_steps
*d4c2dda.* Input matching is exact first, then without letter case, then by dash and underscore variants, each step only if the previous found nothing; a match by a variant is reported as such.

### names_before_synonyms
*d4c2dda.* Each step of `spelling_steps` is tried on names and then on synonyms before the next step begins: exact (names, synonyms), without letter case (names, synonyms), separator variants (names, synonyms), version cut (`accession_version_cut`).

### string_at_release
*d4c2dda.* Between the nodes a typed string can mean, one that exists at the release asked from (without `from_release`: the graph's last release) is taken before one that does not; an Ensembl stable id is never passed over for a synonym spelt the same.

### serving_assembly
*d4c2dda.* A string exists at a release only on the assembly that serves that release; a frozen parallel assembly (human GRCh37, listed at every release from 76) never keeps an old name alive.

### spelling_of_release
*d4c2dda.* A name typed as a database wrote it at the release asked about is an exact match.

### accession_version_cut
*d4c2dda.* A suffix `.<digits>` is cut from an external identifier only when the identifier without it belongs to a database of accessions (RefSeq, UniProt, CCDS, EMBL), never to one of the organism's symbol databases (`NM_000546.9` is `NM_000546`; `SLA.7` is not `SLA`).

### unknown_version
*d4c2dda.* An Ensembl id given with a version it never had is not answered with another version's gene: the answer is nothing, with the versions that exist.

### version_newer_than_graph
*d4c2dda.* A version higher than the last version the graph holds for an id that is still active at the graph's last release may belong to a release newer than the graph: it is answered with the id's gene and flagged (`version_newer_than_graph`). Every other unknown version (lower, in a gap, or of a retired id) is answered with nothing (`unknown_version`).

### accession_versions
*decided.* A RefSeq or UniProt accession written with a version: a version above every version it had, where the graph still holds the accession at the target, is answered with its gene and flagged (as `version_newer_than_graph` for an Ensembl version); a version it never had otherwise (lower, in a gap, above the versions of an accession the graph no longer holds) is answered with nothing (as `unknown_version`); a version it had at an older release is read at that release. `d4c2dda` and `71690ab` cut every accession version (`accession_version_cut`). Family C judges accession versions by this rule ([family C](families/c.md#rules)).

### from_database
*d4c2dda.* With `from_database`, the named database decides which genes a string enters; the support of an entry gene still counts every database that attaches the string to it.

### missing_values
*decided (beyond blank strings, which `71690ab` already refuses with `MissingValueError`).* A string that stands for a missing value (empty; the strings pandas reads as missing by default: `nan`, `NaN`, `NA`, `N/A`, `n/a`, `None`, `null`, `NULL`, `<NA>`, `#N/A`, `#NA`, `-NaN`, `1.#IND`, ...; and `NaT` and UniProt's placeholder name `unknown`; in any letter case, with white space) is refused before any lookup with an input error that names it; it is still looked up when the caller names its source database. The graph drops identifiers without any letter or digit. The harmonizer leaves such features out with that reason.

### single_case_symbols
*decided.* In an organism whose symbols are written in mixed case (mouse), a string written entirely in capitals, or entirely in lower case, that is not exactly a name is read first as a symbol without letter case; where that gives exactly one gene it is the answer and the other reading is among the alternatives.

### no_synonym_target
*decided.* Synonyms are no target database (`synonym_id::<database>` is no `final_database`); a match through a synonym is reported as `synonym:<database>`; the graph's internal names (`<id>.Void`, `<id>.Retired`) are refused as input.

### made_unique_copies
*decided.* A name made unique by a single-cell tool (`X-1`, `X.1`, `X_1`) is a copy of `X`, never another gene: a separator variant never joins such a suffix into another gene's name; in a batch or the harmonizer such a feature is reported as a copy.

### stray_punctuation
*decided.* Stray punctuation and quotes at the ends of an identifier (`A1BG;`, `'AGT'`) are no part of it; a match made after removing them has its own match kind.

### unicode_normalised
*decided.* Input is normalised (NFKC; zero-width characters and byte-order marks removed; Unicode dashes read as `-`); a match made after normalising has its own match kind.

### bounded_arguments
*decided.* Versions and identifiers longer than any real one are refused with an input error that names them (no bare `ValueError`).

### bare_numbers
*decided.* A bare number is read in every database whose identifiers are numbers; where it is an identifier of more than one, the answer is their tie with each database named.

### one_case_schemes
*decided.* Identifiers of schemes written in one letter case (Ensembl, RefSeq, UniProt accessions) are read without regard to letter case; a separator reading that reaches several names is a tie.

### inference_per_release
*decided.* Source inference reads a list's names at each candidate release and answers with the band of releases the list cannot tell apart, with the margin to the next.

### inference_distinct_ids
*decided.* Source inference counts each distinct identifier of a list once, however often the list repeats it (the documented reading of `71690ab`, which counted each occurrence).

### typed_api
*71690ab.* One error hierarchy rooted in `IdtrackError` (input errors, graph-file errors, build errors); typed results with a lossless `to_dict()`; a per-call budget of time and search states (`SearchBudgetExceededError`).

## Which answer is the best

### ties
*d4c2dda.* Candidates are ranked by the complete score; candidates still tied are all returned in ascending order of id and the result says that the answer is a tie.

### own_id_first
*d4c2dda.* An id that exists at the target release is its own answer: with `best` it is answered with itself; with `all` it stays in the answer beside the other genes its history reaches.

### all_order
*d4c2dda.* `strategy="all"` returns its answers best first by the order that chooses the best answer.

### support_counts_all
*d4c2dda.* The databases an identifier is read through, and the support of its gene, count every attachment: to the gene itself and to its transcripts and translations.

### support_ranks_shared_name
*d4c2dda (ranking); judge (scoring).* A name that several genes carry: the gene more name sources attach it to is the `best` answer; `all` returns the set. The judge: a `best` answer that is a non-empty part of the expected set whose genes carry the name in strictly more name sources than every expected gene left out is `undetermined`.

### tie_in_final_database
*d4c2dda.* A tie is judged among the genes that can be named in the final database; genes without an identifier there leave the result and are listed apart.

### authority_symbol_first
*decided.* An official symbol of a name authority (HGNC, MGI, VGNC, NCBI's official symbol) at the release ranks before any synonym of another gene.

### main_line_first
*decided.* Candidates on the strict walk of an id's history (the line the id lives on) rank before candidates reached only through a branch that left it while it lived on.

## History and places

### external_bridge
*d4c2dda.* Where Ensembl's history of an id ends without a successor, an external identifier that the id and a later gene share can carry the id to that gene (a bridge); the result says so.

### bridge_same_place
*decided.* A bridge is evidence of a lineage only when the id leaves one gene for another at the same place: no bridge through synonyms, none through an identifier both stable ids hold at one release, and on one assembly only between genes whose places overlap.

### id_lives_on
*judge.* Ensembl writes no self event beside an event of another id into or out of an id that lives on; an id present in the GTFs on both sides of a mapping session is kept by the strict walk.

### resurrection
*d4c2dda (graph); judge.* An id that Ensembl retired and later brought back under the same stable id is one gene throughout; the judge puts the ids reached only through its return into `broad`.

### two_lives
*decided.* An id brought back after its retirement keeps its two lives apart in the graph, so that travel gives the same lineage forward and backward.

### event_version
*decided.* A history event whose old version did not exist at the session's old release, while its stable id did, is placed with the version the id had there.

### ambiguous_event
*decided.* A history event that two release pairs fit equally is counted and logged; the builder's nearest-pair rule is documented.

### place_moves_marked
*decided.* Where Ensembl's history carries an id to another sequence region or far along it, idtrack follows the history and marks the move in the result and the explained path.

### no_mito_nuclear_step
*decided.* A history step between the mitochondrial genome and a nuclear sequence is not followed.

### copies_apart
*d4c2dda.* Genes on alternative loci, patches and LRG records are marked in the graph; when a gene of the primary sequence is among the candidates, such copies leave the answer and are reported apart (`alternative_locus_genes`).

### lrg_copies
*d4c2dda.* LRG records follow `copies_apart`.

### copies_placed_on_primary
*decided.* A copy on a patch or an alternative locus is placed on the primary assembly (Ensembl's assembly exceptions) before a place test.

### par_y_ids
*decided.* GENCODE's `<id>_PAR_Y` identifiers are read as the Y copy of the pseudoautosomal gene, reported apart.

### par_shared_names
*decided.* A name the X and Y copies of a pseudoautosomal gene share is the X gene's, with the Y copy reported apart.

### lrg_answers_gene
*decided.* An LRG identifier answers the gene it describes.

### primary_assembly_only
*decided.* Graphs are built on the organism's default primary assembly only.

### versioned_organisms
*d4c2dda.* Organisms whose Ensembl ids carry no versions are not supported.

## Names and databases in the graph

### one_kind_per_database
*decided.* One database name means one kind of identifier at every release: where Ensembl gives the label rows and the accession rows of a database one name (`EntrezGene` up to release 87, `RFAM` at every release), the label rows take a name of their own.

### mirbase_genes
*decided.* Gene-level miRBase names are in the graph for genes of a miRNA biotype.

### rfam_display_label
*decided.* An RFAM family label names only the gene whose display name it is.

### spelling_per_edge
*decided.* The graph keeps the spelling per gene, database and release; two spellings of one value of one database for one gene are one value.

## The harmonizer, the server, caches

### harmonize_from_database
*d4c2dda.* `HarmonizeFeatures` takes the source database per dataset (optional).

### cache_by_digest
*decided.* The harmonizer's result cache is keyed by the graph's content digest and a rules version, not by the package version.

### process_workers
*decided.* The web server runs workers as processes over a shared, memory-mapped graph, with a per-call budget.

### no_graph_caches
*d4c2dda.* A graph computes nothing after it is built or loaded (there is no separate caches step).

### no_external_mappers
*71690ab.* idtrack does not wrap other tools (MyGene, g:Profiler, gget, BioMart); the comparison with them is a benchmark family of its own.

## The judge

### truth_judges
*judge.* The accuracy harness, never an earlier version's answers, says what is right; a change of idtrack's answers is accepted only when the judge lists no answer that went from right to wrong.

### judge_copies_apart
*judge.* Where the expected set of a query holds a gene of the primary sequence at the target release, the expected answer is the genes of the primary sequence and the copies on patches and alternative loci are expected in the list reported apart.

### name_before_synonym
*judge.* A string that is a gene's name (display name, or its value in a name table at the release asked about) and another gene's synonym means the name's gene.

### judge_corrections
*judge.* The judge is corrected where it is wrong, with evidence that does not come from idtrack, by revision sets (a published set is never rewritten): each correction is named (`swissprot75`, `living_ids`, `split_name_set`, `names_broad`, `spellings_once`, `projected_hgnc`, `support`, `resurrection`, `capitals_mgi`, `places`, `sessions_effective`, `exceptions`, `exceptions_span`, `core_places`, `copies_apart_early`); [the accuracy harness](harness.md#the-corrections-of-the-judge) states each with its evidence. The families' own revision sets (A*'s `copies_apart_early`, its REST and core-gene readings; A's and E's) are listed in their documents.

### stable_id_form
*judge.* A stable id of Ensembl's scheme is compared without its version where a line compares by stable id: Ensembl's genes, transcripts, translations and exons (`ENS<species><G|T|P|E><11 digits>`) and LRG records with their genes, transcripts and proteins (`LRG_<n>`, `LRG_<n>_g<k>`, `LRG_<n>t<k>`, `LRG_<n>p<k>`); one definition (`harness_lib.STABLE_ID_FORM`) for the harness, the families and the hunts. Until 2026-11-20 LRG gene ids kept their version (family A's `T-lrg-base`, family E's `L-truth-lrg-form`).

### judge_follows_rules
*judge.* Where a rule decided for idtrack after `71690ab` changes what a line asks or expects, the judge follows the rule, by revision sets as for its corrections: the accuracy harness's sets of judge `rules` (`missing_values`, `lower_case_symbols`, `stray_punctuation`, `par_shared_names`; [the accuracy harness](harness.md#the-judge-following-idtracks-decided-rules)), and the families' `decided` judges (family E's `_dec` sets). The judge `published` leaves out every revision set.
