# E: robustness, adversarial inputs, metamorphic relations

Family E asks idtrack what users actually type and how they actually call it: missing values, malformed strings, values of the wrong type, identifiers of another organism, arguments at and beyond their bounds; under budgets; under many configurations of processes, threads, seeds and batches; through relations that must hold between answers; and through the web layer. Its code is `idtrack_reproducibility.families.e`; its data are `<families>/e/` (`queries/<organism>/`, `runs/<tree>/<organism>/`, `runs/<tree>_web/`); the design is in [the methods](../methods.md#e-robustness). It is built on [the common machinery](f0.md).

## What it measures

The specification is the documentation of `71690ab`: the docstrings of `API.convert`, `convert_many`, `infer_source`, `classify` and the other calls, the error hierarchy (rule [`typed_api`](../rules.md#typed_api)) and `HarmonizeFeatures`. An input whose documented outcome is an error is right only with that error class (the line field `expect_error`; on `d4c2dda`, whose errors are Python's built-in classes, the class `71690ab`'s derives from). The expectations of names and identifiers come from Ensembl's tables and the harness's published lines, never from idtrack. The rate is (right + expected_error) / (right + expected_error + wrong + error); undetermined lines (open: the documentation or the truth leaves them open) are apart, and for them the measure is how many named a gene (`answered`).

## The two judges

Where a rule of idtrack is decided but not yet built ([rules](../rules.md), status **decided**), every line is scored both ways:

- judge `now`: the documented behaviour of `71690ab`; the base set;
- judge `decided`: the base set with its revision set `<set>_dec`, which holds only the lines whose expectation a decided rule changes.

The decided rules E's lines touch: [`missing_values`](../rules.md#missing_values), [`single_case_symbols`](../rules.md#single_case_symbols), [`no_synonym_target`](../rules.md#no_synonym_target), [`made_unique_copies`](../rules.md#made_unique_copies), [`stray_punctuation`](../rules.md#stray_punctuation), [`unicode_normalised`](../rules.md#unicode_normalised), [`bounded_arguments`](../rules.md#bounded_arguments), [`bare_numbers`](../rules.md#bare_numbers), [`one_case_schemes`](../rules.md#one_case_schemes), [`inference_per_release`](../rules.md#inference_per_release) and [`lrg_answers_gene`](../rules.md#lrg_answers_gene) (the revision set `e8_identity_dec`, human: the 11 lines asked by an LRG gene id at release 58 expect the Ensembl gene that Ensembl's core database of that release links to the LRG record, its xref `ENS_LRG_gene`; `revise_e --revision e8_identity_dec`). The decided rules [`par_y_ids`](../rules.md#par_y_ids) and [`par_shared_names`](../rules.md#par_shared_names) change no line of E (no line is asked by a `_PAR_Y` id, and no line asks a name of the X and Y copies under the view `all`; the hunts' reproducers and the harness's revision sets of judge `rules` measure them), nor does [`inference_distinct_ids`](../rules.md#inference_distinct_ids) (E8's relation `e8_infer_duplicates` expects it under both judges). When a tree builds them, its `decided` rate is the one that must reach the base set's; until then the `decided` judge measures how far the tree is from them.

## E1 to E9

- **E1, missing values** (hunt finding IN-1; rule `missing_values`): the spellings the rule names (`nan`, `NA`, `N/A`, `None`, `null` in several letter cases and with white space: strata `*_d52`; `NaT` and `unknown`: `*_d81`), the hunt's further spellings (`<NA>`, `#N/A`, `-`, `.`, `?`, ...: `*_other`) and blanks; alone, with `from_database` named (then they are looked up), in a batch at every position, and the other items of those batches (`e1_batch_neighbour`: they keep their answer). `now`: blanks are refused with `MissingValueError`; the other spellings are open (`71690ab` reads them as names; the ones that name a gene are counted). `decided`: the rule's spellings are refused; the rule as built reads pandas' own list of missing values, so `<NA>`, `#N/A` and `#NA` of `*_other` are refused too (revision set `e1_dec2`, 12 lines per organism: `revise_e --revision e1_dec2`); the other spellings of `*_other` stay open, but for `"None."` (with quotes and a period), which rule `stray_punctuation` reads as `None` and so refuses (`e1_dec3`, 4 lines per organism).
- **E2, malformed strings** (IN-3, IN-4, IN-5, IN-6, IN-8, IN-11): forms of real identifiers of nine kinds (`gen_e.py` takes 400 names, 150 synonyms, 150 versioned and 150 unversioned gene ids, 150 transcript and 150 translation ids, 200 NCBI, 100 RefSeq and 100 UniProt ids) whose plain form `d4c2dda` answers right in the harness's sets (F0's check run), so that a wrong answer of a form is the form's own: white space; zero-width characters and byte-order marks (`now`: the gene or nothing; `decided`: the gene, rule `unicode_normalised`); bidi marks and the soft hyphen; control characters; full-width letters and Unicode dashes (`unicode_normalised`); look-alike letters of other scripts; emoji; quotes and stray punctuation (`decided`: the gene, rule `stray_punctuation`); brackets; inner spaces; database prefixes; make-unique suffixes (never another gene, rule `made_unique_copies`); the graph's internal names (`synonym_id::X`: refused, rule `no_synonym_target`); separators flipped and doubled; the letter case of ids; version forms and versions of 3 to 100,000 digits (rule `bounded_arguments`); the letter case of unversioned transcript ids (rule `one_case_schemes`); NCBI ids with leading zeros, as floats and in scientific notation; injection strings, wildcards, numbers, scripts, strings of 10^3 to 10^6 characters (`decided`: refused). `e2_all`: every display name of one gene at 115 with stray punctuation (`;`, `.`, single quotes; IN-11) and with make-unique suffixes (`-1`, `.1`, `_1`; IN-3), the two classes measured whole. A form that is itself an identifier of another gene in Ensembl's tables (any table, any release, in any case and separator form) is open.
- **E2, mouse case** (`e2_case`, IN-2): every MGI symbol of mouse in capitals and in lower case. `now`: the documented rule (an exact synonym of another gene before a match without letter case); `decided` (rule `single_case_symbols`): the MGI symbol's gene.
- **E3, types**: every non-string as an identifier (`ArgumentTypeError`), `np.str_` (a string), every container of a batch, non-string items, a string given as a batch, an empty batch.
- **E4, another organism**: every current gene id of each of the four other organisms (versioned), 5,000 unversioned ids, 3,000 transcripts, 1,000 translations, 3,000 symbols (judged where this organism's tables give one gene, exactly or by letter case; else open) and 1,500 NCBI gene ids (IN-9; open under both judges: rule `bare_numbers` reads a number in every numeric database, so the share answered is the measure), asked in this organism's graph.
- **E5, arguments at and beyond their bounds**: releases below and above the window, as strings, floats, booleans, NaN, NumPy integers; travel backward; database spellings, unknown databases, the synonym prefix; budgets 0, negative, NaN, infinite, of the wrong type, huge, tiny; switches; `infer_source` modes and switches; the other api calls; `classify`.
- **E6, budgets** (`budget.py`): the 200 slowest single calls of each family's runs (the harness's sets through F0's check run, E's own, and A*'s where its answers exist), each under `budget_states` from 1 to 200,000 (the default) and `budget_seconds` from 10^-5 to 30: the outcome must be the documented budget error or exactly the unbounded answer, never a partial answer, and a larger budget must never change a finished answer; in a batch the stopped item is an item error and the others are answered.
- **E7, determinism and concurrency** (`det.py`): the harness's published sets and E's `convert` lines of E1, E2 and E8, under 21 configurations: hash seeds 0, 1, 2 and random, each in 1, 2 and 8 forked processes over one memory-mapped graph (rule [`process_workers`](../rules.md#process_workers)); 8 threads sharing one `API`; shuffled order; batches of 1, 100 and 10,000 against single calls; duplicates; junk items interleaved; after a call that raised; with `explain`. Every canonical answer (`to_dict()` of one `API.convert` call; on `d4c2dda` the result dicts of `convert_identifier`) must be byte-identical to the reference configuration's. Searches are bounded by states, not seconds, so that no answer depends on timing.
- **E8, metamorphic relations** (no truth; rule `relation`): identity (every gene id at every release it exists at, to that release; `e8_identity`, exhaustive); the versioned and unversioned forms of 3,000 current ids give one gene; 3,000 unique names and their case variants give one gene; 2,000 identifiers with and without white space; `infer_source` of GTF lists (ids and names; 100 and 2,000) shuffled and with duplicates. With two calls in a row (`relations_e.py`, 3,000 each): composition (R1 -> R2 -> R3 against R1 -> R3, a violation read leg by leg against the history: idtrack's where the history itself composes); return (an id the history carries unchanged answers itself there and back); final database and back (a gene converted to a database X, and every id it gets read back from X, reaches the gene). And the explained path against Ensembl's tables: every history step of an explained path is an event of the release-115 `stable_id_event` between the releases it names, every bridge an identifier both genes held at its release, the last node the candidate. The harmonizer (`harm_e.py`; human, mouse, pig; three datasets each, made from the truth): every documented refusal raises its class; the union under dataset order and feature order is identical for `merge` "sum" and "flag"; make-unique copies never reach another gene; missing values (`decided`: rule `missing_values`) are left out.
- **E9, the web layer** (`web_e.py`): idtrack's web service, served from the same graphs with the idtrack it pins (`71690ab`), against the package's own answers: per organism 10,000 seeded lines of the harness's sets and 10,000 of E's `convert` lines (and A*'s 1,209 pig lines), as single requests and in batches of 100; every answer must equal `API.convert(...).to_dict()` field for field, every float bit for bit (nothing rounded or reworded); an identifier the service refuses by its documented validation (empty, longer than its limit, control characters) must be refused with `400 invalid_identifier`; and every limit (batch size, identifier length, body and response size, rate) answers its documented status at and beyond its bound. E9 needs idtrack's web service, the repository github.com/theislab/idtrack-web, which is private for now: until it is public a reader can read `results/e/71690ab/web/e9_report.tsv` but cannot run E9 again.
- **hunt_inputs**: the reproducers of the input hunt's findings IN-1 to IN-11 (human, mouse, pig), and IN-7's source inference of the name lists of human 55, 75, 80, 105, 110 and mouse 60, 75, 110 (`e_rule: contains`: the release must be within the answer's band, rule `inference_per_release`).

## Which tree answers what

| part | `d4c2dda` | `71690ab` |
| --- | --- | --- |
| E1 to E4, E5's release and database lines, E8's static and two-call relations, `hunt_inputs` | asked | asked |
| E7 | 17 configurations (no 8-process configuration) | 21 configurations |
| E5's budget, switch and api-call lines, E6, the explained path check, the harmonizer, E9 | not asked (`generations: ["convert"]`: arguments `d4c2dda` does not have) | asked; the reference a later tree is compared with |

Lines whose `generations` leave out a tree's api generation are not asked on it (`run_e.py` counts them). For everything both trees answer, a later tree is compared with `d4c2dda`; for the rest, with `71690ab`.

## The sets

Per organism (manifests `<families>/e/queries/<organism>/<set>.manifest.json`; `judge` is the judge a set belongs to; a `_dec` or `_dec<k>` set revises its base set for the judge `decided`, an `_r<k>` set for both judges). The judge `now` reads the sets that are no `_dec` set in the order they were published; the judge `decided` reads the same, then the `_dec` sets in the order they were published (`score_e.sets_for`): a revision of the truth is published for both judges (`<set>_r<k>` and the decided judge's line revised alike, `<set>_dec<k>`), so the decided line stays the decided judge's truth even where the `_r<k>` set came after it:

| set | judge | revises | pig | mouse | human | rhesus | cynomolgus | what |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `e1` | now | – | 930 | 930 | 930 | 930 | 930 | E1 missing values |
| `e1_dec` | decided | `e1` | 212 | 212 | 212 | 212 | 212 | E1 lines the decided rules change |
| `e2` | now | – | 74,025 | 73,462 | 74,351 | 73,848 | 71,329 | E2 malformed strings |
| `e2_dec` | decided | `e2` | 25,136 | 25,170 | 25,366 | 25,022 | 24,124 | E2 lines the decided rules change |
| `e2_all` | now | – | 98,202 | 463,170 | 228,492 | 107,490 | 95,076 | every display name with stray punctuation and make-unique suffixes |
| `e2_all_dec` | decided | `e2_all` | 49,101 | 231,585 | 114,187 | 53,745 | 47,538 | `e2_all` lines the decided rules change |
| `e2_case` | now | – | – | 154,460 | – | – | – | mouse symbols in another case |
| `e2_case_dec` | decided | `e2_case` | – | 145 | – | – | – | `e2_case` lines the decided rules change |
| `e3` | now | – | 75 | 75 | 75 | 75 | 75 | E3 types |
| `e4` | now | – | 292,685 | 250,033 | 241,998 | 292,935 | 289,817 | E4 another organism's identifiers |
| `e5` | now | – | 107 | 107 | 107 | 107 | 107 | E5 arguments |
| `e5_dec` | decided | `e5` | 1 | 1 | 1 | 1 | 1 | the synonym database as a target |
| `e5_r1` | now | `e5` | – | – | – | – | 1 | revision: the relation `e5_database_spelling` assumed the name chosen for E5 is a symbol of the organism's symbol database; for cynomolgus it is not (`USP24` is a UniProt gene name, no VGNC symbol): the relation has no premise and the line is open |
| `e5_r2` | now | `e5`, `e5_r1` | – | – | – | – | 1 | revision: the same line as a plain line (rule `set`, open), since the scorer applies `open` to every rule but `relation` |
| `e8` | now | – | 14,226 | 18,465 | 15,434 | 14,777 | 15,119 | E8 static relations |
| `e8_identity` | now | – | 1,621,150 | 3,157,745 | 4,137,346 | 2,157,903 | 852,259 | every gene id at every release it exists at (GTF), to that release |
| `hunt_inputs` | now | – | 3 | 12 | 23 | – | – | the input hunt's reproducers |
| `hunt_inputs_dec` | decided | `hunt_inputs` | 1 | 8 | 14 | – | – | reproducers the decided rules change |
| `e1_dec2` | decided | `e1`, `e1_dec` | 12 | 12 | 12 | 12 | 12 | pandas' spellings `<NA>`, `#N/A`, `#NA` refused (rule `missing_values` as built) |
| `e1_dec3`, `e2_dec3` | decided | `e1` and its `_dec` sets; `e2` and its `_dec` sets | 4 · – | 4 · 1 | 4 · 1 | 4 · 1 | 4 · 1 | rule `stray_punctuation` as built, composed with the others (`revise_e --revision stray_decided`): `"None."` is a missing value once the period and quotes are cut (`MissingValueError`); a junk string whose cut form is a name (`CLEC1?`, `e2_junk_wild`) expects that name's genes |
| `e8_identity_dec` | decided | `e8_identity` | – | – | 11 | – | – | an LRG gene id answers its gene (rule `lrg_answers_gene`) |
| `e2_r1`, `e2_all_r1`, `hunt_inputs_r1` | now | their base set | 2,750 · 49,101 · – | 2,750 · 231,585 · – | 2,750 · 114,246 · 3 | 2,750 · 53,745 · – | 2,500 · 47,538 · – | the truth's reach (`e2_names_every_release`): every made-unique line records the releases whose synonyms the truth reads (`names_read_at`); a form Ensembl gives exactly to other genes at a release of the window expects them too (pig 4, mouse 35, human 230, rhesus 2 forms) |
| `e2_dec2`, `e2_all_dec2`, `hunt_inputs_dec2` | decided | their base set and its `_dec` | as `_r1` | as `_r1` | as `_r1` | as `_r1` | as `_r1` | the same revision of the decided judge's lines |
| `e2_case_r1` | now | `e2_case` | – | 9 | – | – | – | revision `e2_case_r1`: the documented rule (an exact name or synonym before a match without letter case) names every gene that carries the form exactly at 115, read from Ensembl's synonyms of every cross-reference (`core_synonyms/115.tsv.gz`) beside BioMart's, the MGI symbol's own gene too (1 `e2_case_upper_competing`, `ESP1`; 8 `e2_case_*_folded` lines no longer open); the decided judge's lines are unchanged (rule `single_case_symbols`: the MGI symbol's gene) |

Each manifest gives the lines per stratum (`counts`); the strata are named `<set>_<form>_<kind of identifier>` (for example `e2_zw_name`: zero-width characters in names; `e2_punct_gv`: stray punctuation on versioned gene ids).

## Truth

Nothing of idtrack. The made-unique forms are read against every release's names too: the core cross-references and the GTF of each release (revision `e2_names_every_release`, `revise_e.py`); synonyms only where the truth holds them (BioMart's synonym tables, the core synonym extract of 115), recorded on each such line (`names_read_at`). `truth_e.py` reads the one truth store: `Names` (every string the store knows for the organism, as Ensembl's tables write it: the GTF's gene ids and names of every release and every column of every BioMart table of every release the store holds, plus VGNC's own symbols; folded as idtrack's case and separator steps fold them, to know which strings another gene owns) and `Current` (the genes, transcripts and translations of GTF 115). The bases of E2's forms are the harness's published lines whose plain form `d4c2dda` answers right. E8's explained-path check reads the release-115 `stable_id_event` and BioMart's cross-references; the harmonizer's datasets are built from GTF 115 and BioMart 115. The documentation of `71690ab` decides errors and what is left open.

## Rules

Scoring is F0's scorer under each judge ([scoring and outcomes](f0.md#scoring-and-outcomes)): the harness's rules, `expect_error`, `safe` (a form idtrack may not read: the gene or nothing is right, another gene wrong), `open` and `relation`; then family E's two rules (`score_e.py`, applied to the lines that carry `e_rule`): `or_error` (an error of the classes named is right as well as the line's rule, for a bound the documentation leaves to the tree, e.g. a 10-digit version or a 1,000-character identifier) and `contains` (the expected values must be among the answer's, e.g. a release within `infer_source`'s band).

The rules of idtrack the expectations rest on: [`spelling_steps`](../rules.md#spelling_steps), [`names_before_synonyms`](../rules.md#names_before_synonyms), [`string_at_release`](../rules.md#string_at_release), [`accession_version_cut`](../rules.md#accession_version_cut), [`from_database`](../rules.md#from_database), [`unknown_version`](../rules.md#unknown_version), [`version_newer_than_graph`](../rules.md#version_newer_than_graph), [`typed_api`](../rules.md#typed_api), [`harmonize_from_database`](../rules.md#harmonize_from_database), and the decided rules of [The two judges](#the-two-judges).

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)); `<t>` is `d4c2dda` or `71690ab` (or a later tree).

```
# the sets (once per organism; reads only the truth store and the harness's published sets)
python -m idtrack_reproducibility.families.e.gen_e --organism <o> [--sets e1,e2,...] [--dry]
python -m idtrack_reproducibility.families.e.revise_e --organism <o> --revision e5_r1|e5_r2|e1_dec2|e8_identity_dec|e2_names_every_release|e2_case_r1|stray_decided [--dry]
# E1 to E5, E8 static: F0's runner with E's adapter, then both judges and E's rules
python -m idtrack_reproducibility.families.e.run_e --organism <o> --clone <source> --graph <graph> \
    --out <families>/e/runs/<t>/<o> --submit 8
python -m idtrack_reproducibility.families.e.score_e --organism <o> --run <families>/e/runs/<t>/<o>
# E7, E8 with two calls and the explained path, E6, the harmonizer (the tree's Python)
<python of the tree> -m idtrack_reproducibility.families.e.det submit --organism <o> --clone <source> --graph <graph> --out <families>/e/runs/<t>/<o>
<python of the tree> -m idtrack_reproducibility.families.e.det compare --out <families>/e/runs/<t>/<o>
<python of the tree> -m idtrack_reproducibility.families.e.relations_e submit --organism <o> --clone <source> --graph <graph> --out <families>/e/runs/<t>/<o>
<python of the tree> -m idtrack_reproducibility.families.e.budget submit --organism <o> --clone <source> --graph <graph> \
    --out <families>/e/runs/<t>/<o> --from harness=<families>/f0/runs/check_<t>/<o> --from e=<families>/e/runs/<t>/<o> \
    [--from astar=<families>/astar/runs/<tree>/<o>]
<python of the tree> -m idtrack_reproducibility.families.e.harm_e --organism <o> --clone <source> --graph <graph> --out <families>/e/runs/<t>/<o>
# E9: the web service's own environment (its idtrack is the commit it pins)
<python of the web service> -m idtrack_reproducibility.families.e.web_e --out <families>/e/runs/<t>_web [--n 10000]
# causes, comparison with the reference tree, limitations, the table over trees
python -m idtrack_reproducibility.families.e.causes_e --organism <o> --run <families>/e/runs/<t>/<o>
python -m idtrack_reproducibility.families.compare <families>/e/runs/d4c2dda/<o>/scores_now <families>/e/runs/<t>/<o>/scores_now
python -m idtrack_reproducibility.families.e.limits_e --runs <families>/e/runs/<t>/<o> [...] --out <tsv>
python -m idtrack_reproducibility.families.e.summary_e --trees d4c2dda,71690ab --out <families>/e/e_rates.tsv
```

E6, the explained path, the harmonizer and E9 run only on trees of the `convert` generation. The harmonizer runs in an environment with idtrack's `harmonize` extra (`anndata` 0.13 or later, which writes and reads the datasets' `.h5ad` files).

## Results

Runs `<families>/e/runs/d4c2dda/<organism>/` and `runs/71690ab/<organism>/`, both on the reference graphs, and `runs/71690ab_web/e9/`. Result files per tree and organism: `results/e/<tree>/<organism>/now_report.tsv` and `decided_report.tsv` (F0's `report.tsv` under each judge), `e_summary.tsv` (per judge and stratum: outcomes, rate, Wilson interval, `answered`, `undetermined_answered`), `e7_report.tsv`, `e8_report.tsv`; for `71690ab` also `e6_report.tsv`, `compare_vs_d4c2dda_now.tsv` and `compare_vs_d4c2dda_decided.tsv`; and `results/e/71690ab/web/e9_report.tsv`. No line timed out.

**`d4c2dda` and `71690ab` give the same answers on every line both are asked**, under both judges and in all five organisms, but 32 per organism (`compare_vs_d4c2dda_<judge>.tsv`): the 30 batch lines in NumPy and pandas containers, which `d4c2dda` refuses with its documented `TypeError` and `71690ab` answers right, and the 2 iterables a batch is not documented to take (a generator, a set), refused on `d4c2dda` and open on `71690ab`. No answer moved from right to wrong or to an error.

**Totals** (sums over the strata of `results/e/<tree>/<organism>/e_summary.tsv`; `open named a gene` is `undetermined_answered`):

| tree | judge | organism | lines | right | expected error | wrong | undetermined | error | open, named a gene |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `71690ab` | now | pig | 2,101,403 | 2,086,663 | 158 | 160 | 13,972 | 450 | 3,203 |
| `71690ab` | now | mouse | 4,118,459 | 4,094,593 | 158 | 560 | 22,698 | 450 | 13,956 |
| `71690ab` | now | human | 4,698,756 | 4,686,113 | 158 | 247 | 11,787 | 451 | 2,913 |
| `71690ab` | now | rhesus | 2,648,065 | 2,634,602 | 158 | 141 | 12,714 | 450 | 2,262 |
| `71690ab` | now | cynomolgus | 1,324,712 | 1,311,642 | 158 | 115 | 12,347 | 450 | 1,474 |
| `71690ab` | decided | pig | 2,101,403 | 2,013,875 | 158 | 73,612 | 13,758 | 0 | 3,062 |
| `71690ab` | decided | mouse | 4,118,459 | 3,853,175 | 158 | 242,648 | 22,478 | 0 | 13,837 |
| `71690ab` | decided | human | 4,698,756 | 4,563,847 | 158 | 123,181 | 11,570 | 0 | 2,769 |
| `71690ab` | decided | rhesus | 2,648,065 | 2,556,782 | 158 | 78,624 | 12,501 | 0 | 2,146 |
| `71690ab` | decided | cynomolgus | 1,324,712 | 1,240,981 | 158 | 71,439 | 12,134 | 0 | 1,382 |
| `d4c2dda` | now | pig | 2,101,356 | 2,086,627 | 154 | 160 | 13,965 | 450 | 3,198 |
| `d4c2dda` | now | mouse | 4,118,409 | 4,094,556 | 154 | 558 | 22,691 | 450 | 13,951 |
| `d4c2dda` | now | human | 4,698,704 | 4,686,076 | 154 | 243 | 11,780 | 451 | 2,908 |
| `d4c2dda` | now | rhesus | 2,648,018 | 2,634,566 | 154 | 141 | 12,707 | 450 | 2,257 |
| `d4c2dda` | now | cynomolgus | 1,324,665 | 1,311,606 | 154 | 115 | 12,340 | 450 | 1,470 |
| `d4c2dda` | decided | pig | 2,101,356 | 2,013,839 | 154 | 73,612 | 13,751 | 0 | 3,057 |
| `d4c2dda` | decided | mouse | 4,118,409 | 3,853,138 | 154 | 242,646 | 22,471 | 0 | 13,832 |
| `d4c2dda` | decided | human | 4,698,704 | 4,563,810 | 154 | 123,177 | 11,563 | 0 | 2,764 |
| `d4c2dda` | decided | rhesus | 2,648,018 | 2,556,746 | 154 | 78,624 | 12,494 | 0 | 2,141 |
| `d4c2dda` | decided | cynomolgus | 1,324,665 | 1,240,945 | 154 | 71,439 | 12,127 | 0 | 1,378 |

The errors under `now` are the versions of more than 4,300 digits (Python's `ValueError`, IN-8); under `decided` they count as wrong (the rule `bounded_arguments` refuses them with an input error).

**Per stratum** on `71690ab` (rows of `results/e/71690ab/<organism>/e_summary.tsv`: the rate under `now` · the rate under `decided` (lines); for open lines, how many named a gene). The full table has every stratum; these are the ones that show each form's behaviour:

| group | stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| E1 | `e1_single_blank` | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) |
| E1 | `e1_batch_blank` | 1.0000 · 1.0000 (27) | 1.0000 · 1.0000 (27) | 1.0000 · 1.0000 (27) | 1.0000 · 1.0000 (27) | 1.0000 · 1.0000 (27) |
| E1 | `e1_single_d52` | open, 29 named a gene · 0.0000 (41) | open, 16 named a gene · 0.0000 (41) | open, 23 named a gene · 0.0000 (41) | open, 23 named a gene · 0.0000 (41) | open, 23 named a gene · 0.0000 (41) |
| E1 | `e1_single_d81` | open, 6 named a gene · 0.0000 (12) | open, 12 named a gene · 0.0000 (12) | open, 12 named a gene · 0.0000 (12) | open, 6 named a gene · 0.0000 (12) | open, 0 named a gene · 0.0000 (12) |
| E1 | `e1_single_other` | open, 1 named a gene · open, 1 named a gene (20) | open, 2 named a gene · open, 2 named a gene (20) | open, 1 named a gene · open, 1 named a gene (20) | open, 1 named a gene · open, 1 named a gene (20) | open, 1 named a gene · open, 1 named a gene (20) |
| E1 | `e1_batch_d52` | open, 87 named a gene · 0.0000 (123) | open, 48 named a gene · 0.0000 (123) | open, 69 named a gene · 0.0000 (123) | open, 69 named a gene · 0.0000 (123) | open, 69 named a gene · 0.0000 (123) |
| E1 | `e1_batch_neighbour` | 1.0000 · 1.0000 (438) | 1.0000 · 1.0000 (438) | 1.0000 · 1.0000 (438) | 1.0000 · 1.0000 (438) | 1.0000 · 1.0000 (438) |
| E1 | `e1_fromdb_d52` | open, 0 named a gene · open, 0 named a gene (82) | open, 0 named a gene · open, 0 named a gene (82) | open, 0 named a gene · open, 0 named a gene (82) | open, 0 named a gene · open, 0 named a gene (82) | open, 0 named a gene · open, 0 named a gene (82) |
| E2 | `e2_ws_name` | 1.0000 · 1.0000 (2,400) | 1.0000 · 1.0000 (2,400) | 1.0000 · 1.0000 (2,400) | 1.0000 · 1.0000 (2,400) | 1.0000 · 1.0000 (2,400) |
| E2 | `e2_idcase_gu` | 1.0000 · 1.0000 (300) | 1.0000 · 1.0000 (300) | 1.0000 · 1.0000 (300) | 1.0000 · 1.0000 (300) | 1.0000 · 1.0000 (300) |
| E2 | `e2_separators_name` | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (42) | 1.0000 · 1.0000 (147) | 1.0000 · 1.0000 (48) | 1.0000 · 1.0000 (9) |
| E2 | `e2_zw_name` | 1.0000 · 0.0000 (2,800) | 1.0000 · 0.0000 (2,800) | 1.0000 · 0.0000 (2,800) | 1.0000 · 0.0000 (2,800) | 1.0000 · 0.0000 (2,800) |
| E2 | `e2_bidi_name` | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) |
| E2 | `e2_fullwidth_name` | 1.0000 · 0.0000 (400) | 1.0000 · 0.0000 (400) | 1.0000 · 0.0000 (400) | 1.0000 · 0.0000 (400) | 1.0000 · 0.0000 (400) |
| E2 | `e2_dash_name` | 1.0000 · 0.0000 (24) | 1.0000 · 0.0000 (112) | 1.0000 · 0.0000 (392) | 1.0000 · 0.0000 (128) | 1.0000 · 0.0000 (24) |
| E2 | `e2_ctrl_name` | 1.0000 · 1.0000 (2,800) | 1.0000 · 1.0000 (2,800) | 1.0000 · 1.0000 (2,800) | 1.0000 · 1.0000 (2,800) | 1.0000 · 1.0000 (2,800) |
| E2 | `e2_lookalike_name` | 1.0000 · 1.0000 (754) | 1.0000 · 1.0000 (221) | 1.0000 · 1.0000 (769) | 1.0000 · 1.0000 (737) | 1.0000 · 1.0000 (754) |
| E2 | `e2_emoji_name` | 1.0000 · 1.0000 (800) | 1.0000 · 1.0000 (800) | 1.0000 · 1.0000 (800) | 1.0000 · 1.0000 (800) | 1.0000 · 1.0000 (800) |
| E2 | `e2_punct_name` | 1.0000 · 0.0036 (2,800) | 0.9993 · 0.0271 (2,800) | 0.9993 · 0.0646 (2,800) | 1.0000 · 0.0000 (2,800) | 1.0000 · 0.0000 (2,800) |
| E2 | `e2_brackets_name` | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) |
| E2 | `e2_innerspace_name` | 1.0000 · 1.0000 (400) | 1.0000 · 1.0000 (400) | 1.0000 · 1.0000 (400) | 1.0000 · 1.0000 (400) | 1.0000 · 1.0000 (400) |
| E2 | `e2_prefix` | 1.0000 · 1.0000 (500) | 1.0000 · 1.0000 (500) | 1.0000 · 1.0000 (500) | 1.0000 · 1.0000 (500) | 1.0000 · 1.0000 (500) |
| E2 | `e2_makeunique_name` | 0.9985 · 0.9985 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) | 1.0000 · 1.0000 (2,000) |
| E2 | `e2_internal_name` | 0.9650 · 0.8950 (400) | 0.9775 · 0.9450 (400) | 0.9850 · 0.9375 (400) | 0.9575 · 0.8150 (400) | 0.9750 · 0.9675 (400) |
| E2 | `e2_internal_syn` | 1.0000 · 0.0000 (150) | 0.9933 · 0.0000 (150) | 1.0000 · 0.0000 (150) | 1.0000 · 0.0000 (150) | 1.0000 · 0.0000 (100) |
| E2 | `e2_longversion` | 0.6250 · 0.2500 (1,200) | 0.6250 · 0.2500 (1,200) | 0.6250 · 0.2500 (1,200) | 0.6250 · 0.2500 (1,200) | 0.6250 · 0.2500 (1,200) |
| E2 | `e2_long` | 1.0000 · 0.2500 (8) | 1.0000 · 0.2500 (8) | 1.0000 · 0.2500 (8) | 1.0000 · 0.2500 (8) | 1.0000 · 0.2500 (8) |
| E2 | `e2_versionform_gv` | 1.0000 · 1.0000 (750) | 1.0000 · 1.0000 (750) | 1.0000 · 1.0000 (750) | 1.0000 · 1.0000 (750) | 1.0000 · 1.0000 (750) |
| E2 | `e2_txcase_tx_unversioned` | 1.0000 · 0.0000 (144) | 1.0000 · 0.0000 (134) | 1.0000 · 0.0000 (156) | 1.0000 · 0.0000 (162) | 1.0000 · 0.0000 (156) |
| E2 | `e2_txcase_tx_versioned` | 1.0000 · 1.0000 (156) | 1.0000 · 1.0000 (166) | 1.0000 · 1.0000 (144) | 1.0000 · 1.0000 (138) | 1.0000 · 1.0000 (144) |
| E2 | `e2_numbers_ncbi` | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) | 1.0000 · 1.0000 (1,600) |
| E2 | `e2_junk_inj` | 1.0000 · 1.0000 (17) | 1.0000 · 1.0000 (17) | 1.0000 · 1.0000 (17) | 1.0000 · 1.0000 (17) | 1.0000 · 1.0000 (17) |
| E2 | `e2_all_punct` | 0.9995 · 0.0089 (49,101) | 0.9979 · 0.0596 (231,585) | 0.9985 · 0.1396 (114,246) | 1.0000 · 0.0000 (53,745) | 1.0000 · 0.0000 (47,538) |
| E2 | `e2_all_makeunique` | 0.9978 · 0.9978 (49,101) | 0.9999 · 0.9999 (231,585) | 0.9997 · 0.9997 (114,246) | 0.9979 · 0.9979 (53,745) | 0.9979 · 0.9979 (47,538) |
| E2 | `e2_case_upper_alone` | – | 1.0000 · 1.0000 (70,980) | – | – | – |
| E2 | `e2_case_upper_competing` | – | 0.9237 · 0.0000 (131) | – | – | – |
| E2 | `e2_case_lower_competing` | – | 0.8571 · 0.1429 (14) | – | – | – |
| E3 | `e3_single_nonstring` | 1.0000 · 1.0000 (18) | 1.0000 · 1.0000 (18) | 1.0000 · 1.0000 (18) | 1.0000 · 1.0000 (18) | 1.0000 · 1.0000 (18) |
| E3 | `e3_single_npstr` | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) |
| E3 | `e3_batch_container` | 1.0000 · 1.0000 (40) | 1.0000 · 1.0000 (40) | 1.0000 · 1.0000 (40) | 1.0000 · 1.0000 (40) | 1.0000 · 1.0000 (40) |
| E3 | `e3_batch_nonstring_item` | 1.0000 · 1.0000 (10) | 1.0000 · 1.0000 (10) | 1.0000 · 1.0000 (10) | 1.0000 · 1.0000 (10) | 1.0000 · 1.0000 (10) |
| E3 | `e3_batch_str` | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) | 1.0000 · 1.0000 (1) |
| E4 | `e4_gene_versioned` | 1.0000 · 1.0000 (238,685) | 1.0000 · 1.0000 (196,033) | 1.0000 · 1.0000 (187,998) | 1.0000 · 1.0000 (238,935) | 1.0000 · 1.0000 (235,817) |
| E4 | `e4_gene_unversioned` | 1.0000 · 1.0000 (20,000) | 1.0000 · 1.0000 (20,000) | 1.0000 · 1.0000 (20,000) | 1.0000 · 1.0000 (20,000) | 1.0000 · 1.0000 (20,000) |
| E4 | `e4_transcript` | 1.0000 · 1.0000 (12,000) | 1.0000 · 1.0000 (12,000) | 1.0000 · 1.0000 (12,000) | 1.0000 · 1.0000 (12,000) | 1.0000 · 1.0000 (12,000) |
| E4 | `e4_translation` | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) |
| E4 | `e4_symbol_exact` | 1.0000 · 1.0000 (4,180) | 0.5000 · 0.5000 (4) | 1.0000 · 1.0000 (6,939) | 1.0000 · 1.0000 (5,346) | 1.0000 · 1.0000 (5,657) |
| E4 | `e4_symbol_case` | 1.0000 · 1.0000 (375) | 1.0000 · 1.0000 (8,363) | 1.0000 · 1.0000 (541) | 1.0000 · 1.0000 (478) | 1.0000 · 1.0000 (502) |
| E4 | `e4_symbol_shared` | open, 2,800 named a gene · open, 2,800 named a gene (2,817) | open, 977 named a gene · open, 977 named a gene (978) | open, 2,003 named a gene · open, 2,003 named a gene (2,003) | open, 1,911 named a gene · open, 1,911 named a gene (1,932) | open, 743 named a gene · open, 743 named a gene (743) |
| E4 | `e4_ncbi_foreign` | open, 0 named a gene · open, 0 named a gene (6,000) | open, 149 named a gene · open, 149 named a gene (6,000) | open, 124 named a gene · open, 124 named a gene (6,000) | open, 0 named a gene · open, 0 named a gene (6,000) | open, 0 named a gene · open, 0 named a gene (6,000) |
| E5 | `e5_release_below` | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) |
| E5 | `e5_release_str` | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) |
| E5 | `e5_release_float` | 1.0000 · 1.0000 (3) | 1.0000 · 1.0000 (3) | 1.0000 · 1.0000 (3) | 1.0000 · 1.0000 (3) | 1.0000 · 1.0000 (3) |
| E5 | `e5_database_unknown` | 1.0000 · 1.0000 (6) | 1.0000 · 1.0000 (6) | 1.0000 · 1.0000 (6) | 1.0000 · 1.0000 (6) | 1.0000 · 1.0000 (6) |
| E5 | `e5_database_synonym` | open, 0 named a gene · 0.0000 (1) | open, 0 named a gene · 0.0000 (1) | open, 0 named a gene · 0.0000 (1) | open, 0 named a gene · 0.0000 (1) | open, 0 named a gene · 0.0000 (1) |
| E5 | `e5_budget_invalid` | 1.0000 · 1.0000 (7) | 1.0000 · 1.0000 (7) | 1.0000 · 1.0000 (7) | 1.0000 · 1.0000 (7) | 1.0000 · 1.0000 (7) |
| E5 | `e5_budget_tiny` | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) | 1.0000 · 1.0000 (2) |
| E5 | `e5_api_error` | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) | 1.0000 · 1.0000 (9) |
| E8 | `e8_identity_versioned_primary` | 1.0000 · 1.0000 (1,144,462) | 1.0000 · 1.0000 (2,120,132) | 1.0000 · 1.0000 (2,428,765) | 1.0000 · 1.0000 (1,298,548) | 1.0000 · 1.0000 (852,259) |
| E8 | `e8_identity_unversioned_primary` | 1.0000 · 1.0000 (476,688) | 1.0000 · 1.0000 (1,024,739) | 1.0000 · 1.0000 (1,494,120) | 1.0000 · 1.0000 (859,355) | – |
| E8 | `e8_identity_versioned_patch` | – | 1.0000 · 1.0000 (12,874) | 1.0000 · 1.0000 (214,461) | – | – |
| E8 | `e8_versioned_pair` | 1.0000 · 1.0000 (6,000) | 1.0000 · 1.0000 (6,000) | 1.0000 · 1.0000 (6,000) | 1.0000 · 1.0000 (6,000) | 1.0000 · 1.0000 (6,000) |
| E8 | `e8_case_pair` | 1.0000 · 1.0000 (4,106) | 1.0000 · 1.0000 (8,345) | 1.0000 · 1.0000 (5,314) | 1.0000 · 1.0000 (4,657) | 1.0000 · 1.0000 (5,599) |
| E8 | `e8_whitespace_pair` | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (4,000) | 1.0000 · 1.0000 (3,400) |
| E8 | `e8_infer_order` | 1.0000 · 1.0000 (96) | 1.0000 · 1.0000 (96) | 1.0000 · 1.0000 (96) | 1.0000 · 1.0000 (96) | 1.0000 · 1.0000 (96) |
| E8 | `e8_infer_duplicates` | 1.0000 · 1.0000 (24) | 0.9583 · 0.9583 (24) | 0.9583 · 0.9583 (24) | 0.9583 · 0.9583 (24) | 1.0000 · 1.0000 (24) |
| hunt | `hunt_IN-1` | open, 1 named a gene · 0.0000 (1) | open, 3 named a gene · 0.0000 (3) | open, 4 named a gene · 0.0000 (4) | – | – |
| hunt | `hunt_IN-2` | 0.0000 · 0.0000 (2) | open, 5 named a gene · 0.0000 (5) | – | – | – |
| hunt | `hunt_IN-3` | – | – | 0.0000 · 0.0000 (3) | – | – |
| hunt | `hunt_IN-7` | – | 0.3333 · 0.3333 (3) | 0.2000 · 0.2000 (5) | – | – |
| hunt | `hunt_IN-8` | – | – | 0.0000 · 0.0000 (1) | – | – |
| hunt | `hunt_IN-11` | – | 0.0000 · 0.0000 (1) | 0.0000 · 0.0000 (1) | – | – |

Read by form: white space, letter case of ids, separators, versioned transcript ids and NCBI numbers are found; zero-width characters, full-width letters, Unicode dashes, stray punctuation and unversioned transcript ids in another letter case never give another gene today (`now` 1.0000, the rule `safe`) and are not yet read as the decided rules read them (`decided`); control characters, bidi marks, look-alikes, emoji, brackets, inner spaces and prefixes are refused or not found, under both judges; the graph's internal names and the stray punctuation and make-unique suffixes over every name are where idtrack answers another gene today. E3 and E5: every documented error class is raised on `71690ab`. E4: no gene id, transcript or translation of another organism is found in any organism; symbols judged exactly or by case are right but 2 mouse genes whose display name is a projected human HGNC symbol; NCBI ids of another organism are read through numeric ids of other databases in mouse (149 of 6,000) and human (124 of 6,000). E8: identity holds on every line but 11 human LRG ids (the NCBI36 GTFs write `LRG_<n>_g1`, the graph `LRG_<n>_g1.1`: the judge's form); `infer_source` gives the same answer to a shuffled list always, and to a list with repeats in 23 of 24 lists in mouse, human and rhesus. On `d4c2dda` the same table holds but for the lines it is not asked and the 32 lines above (`results/e/d4c2dda/<organism>/e_summary.tsv`).

**E6, budgets** (`results/e/71690ab/<organism>/e6_report.tsv`): 200 lines of each family's runs per organism (the harness's, E's, and A*'s where its answers existed: pig, rhesus, cynomolgus), 2,600 lines; every line under every budget gave either the documented budget error or exactly the unbounded answer: 0 lines with a violation, 0 partial answers, 0 answers changed by a larger budget, 0 other errors, 0 batch violations. The states needed to finish (median, maximum) are in the table, for example the harness's lines of human 1 and 2,000, of rhesus 200 and 5,000.

**E7, determinism** (`results/e/<tree>/<organism>/e7_report.tsv`): `71690ab`, 21 configurations; `d4c2dda`, 17; every canonical answer identical to the reference in every configuration, organism and tree (0 different, 0 missing). Lines per configuration: pig 568,193, mouse 798,614, human 751,798, rhesus 553,623, cynomolgus 381,700.

**E8 with two calls and the explained path** (`results/e/<tree>/<organism>/e8_report.tsv`; the same on both trees but the path, which needs `71690ab`):

| organism | composition: holds / by the history / refuted (of 3,000) | return: holds (of 3,000) | final database and back: holds (of ids) | path: history steps that are events (of steps) | bridge steps held / undetermined | final steps held (of steps) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pig | 2,798 / 200 / 2 | 3,000 | 2,432 of 2,432 | 7,821 of 7,821 | 27 / 449 | 31 of 31 |
| mouse | 2,916 / 82 / 2 | 3,000 | 2,678 of 2,678 | 12,919 of 12,919 | 70 / 306 | 163 of 163 |
| human | 2,863 / 134 / 3 | 3,000 | 2,057 of 2,057 | 9,515 of 9,515 | 74 / 428 | 17 of 17 |
| rhesus | 2,915 / 83 / 2 | 3,000 | 2,383 of 2,383 | 6,198 of 6,198 | 0 / 322 | 19 of 19 |
| cynomolgus | 2,954 / 46 / 0 | 3,000 | 1,117 of 1,117 | 2,765 of 2,765 | 0 / 190 | 5 of 5 |

A composition that does not hold where the history itself composes (`refuted`, 9 over the five organisms) is a bridge where the history has no successor (rule [`external_bridge`](../rules.md#external_bridge)); the others do not compose because the history does not (splits and merges). In the path check, 40 mouse history steps (session 102 -> 103) carry the GTF's version where Ensembl's event names another (`history_event_versions_differ`; rule [`event_version`](../rules.md#event_version)'s class); no bridge step is contradicted by the tables, and those undetermined are at releases without the table that would decide them. (Bridge steps are the near and far ends of each bridge, as `e8_report.tsv` counts them.)

**The harmonizer** (`<families>/e/runs/71690ab/<organism>/harm/report.tsv`; human, mouse, pig): every documented refusal raised its class (`args`), the union under dataset order and under feature order is identical for `merge` "sum" and "flag" (`order_*`), and make-unique copies never reached another gene (`makeunique_copy`), with no violation in any organism; under `decided`, the missing values not yet left out are the rows `missing_*` with violations (human `NA`, `NaT`, `nan`, `unknown`; mouse `NaT`; pig `NA`, `NaT`, `nan`).

**E9, the web layer** (`results/e/71690ab/web/e9_report.tsv`): every answer equal to the package's, field for field:

| organism | harness lines: equal (asked) | E's lines: equal + refused by validation (asked) | batches of 100: equal | A* lines: equal (asked) |
| --- | ---: | ---: | ---: | ---: |
| pig | 10,000 (10,000) | 9,591 + 409 (10,000) | 2,140 | 1,209 (1,209) |
| mouse | 10,000 (10,000) | 9,733 + 267 (10,000) | 1,879 | – |
| human | 10,000 (10,000) | 9,568 + 432 (10,000) | 1,869 | – |
| rhesus | 10,000 (10,000) | 9,616 + 384 (10,000) | 1,618 | – |
| cynomolgus | 10,000 (10,000) | 9,622 + 378 (10,000) | 482 | – |

Every limit answered its documented status (19 checks, each `ok`): batch size at and beyond its maximum (1,000), identifier length at and beyond its maximum (256), empty and control-character identifiers, body and response beyond their maxima and a small response within, the rate burst accepted, then refused with `Retry-After`, health never rate-limited, an unknown parameter, a release given as a JSON string, an unknown organism, a body that is not JSON, an unknown release and an unknown database.

### Causes

Every wrong answer and error of both judges on both trees has a cause (`causes_e`; `<families>/e/runs/<tree>/<organism>/scores_<judge>/causes_summary.tsv`, no `unexplained` row). Answers on `71690ab`, `now` · `decided`, summed over strata:

| cause | class | rule | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `long_version_value_error`: a version of more than 4,300 digits raises Python's `ValueError` (IN-8) | a | `bounded_arguments` | 450 · 450 | 450 · 450 | 451 · 451 | 450 · 450 | 450 · 450 |
| `stray_punctuation_artifact_name`: `X;`, `X.`, `'X'` reach another gene through names of Ensembl's tables that carry the same stray character (IN-11) | a | `stray_punctuation` | 27 · 27 | 499 · 499 | 180 · 180 | 0 · 0 | 0 · 0 |
| `separator_joins_suffix`: `X-1`, `X.1`, `X_1` joined by the case or separator step into another gene's name (IN-3) | b | `made_unique_copies` | 89 · 89 | 28 · 28 | 34 · 34 | 93 · 93 | 79 · 79 |
| `internal_name_accepted`: `synonym_id::X` read as an identifier (IN-6) | a | `no_synonym_target` | 14 · 192 | 10 · 172 | 7 · 178 | 17 · 224 | 10 · 113 |
| `synonym_spelling_merged`: the asked string is carried in two letter cases at the release read (Ensembl's synonyms of every cross-reference); the graph keeps one spelling per database and release, so the exact match holds the genes of both cases, or the case step answers where an exact synonym exists | c | `spelling_steps`; `spelling_per_edge` (decided) | 0 · 0 | 13 · 0 | 0 · 0 | 0 · 0 | 0 · 0 |
| `lrg_record_not_its_gene`: an LRG gene id answers the LRG record, not the Ensembl gene it describes | a | `lrg_answers_gene` | 0 · 0 | 0 · 0 | 0 · 11 | 0 · 0 | 0 · 0 |
| `source_inference_reads_names_at_snapshot`: name lists placed at a later release (IN-7) | a | `inference_per_release` | 0 · 0 | 2 · 2 | 4 · 4 | 0 · 0 | 0 · 0 |
| `infer_duplicates_counted`: a list with repeats placed at another release than without them | a | `inference_distinct_ids` | 0 · 0 | 1 · 1 | 1 · 1 | 1 · 1 | 0 · 0 |
| `lowercase_exact_synonym_other_organism`: a lower-case symbol that is exactly another gene's synonym (pig `afp`, `Has2`) | b | – (`single_case_symbols` covers mixed-case organisms only) | 2 · 2 | 0 · 0 | 0 · 0 | 0 · 0 | 0 · 0 |
| `mouse_hgnc_projection`: a mouse gene whose display name is a projected human HGNC symbol | d | – | 0 · 0 | 2 · 2 | 0 · 0 | 0 · 0 | 0 · 0 |
| `nonalnum_name_in_graph`: mouse `_` reaches a gene through an NCBI synonym `.` by the separator step | a | `missing_values` (the graph part) | 0 · 0 | 1 · 1 | 0 · 0 | 0 · 0 | 0 · 0 |
| `missing_value_read_as_name`, `missing_value_not_refused` (IN-1) | a | `missing_values` | 0 · 225 | 0 · 226 | 0 · 228 | 0 · 224 | 0 · 224 |
| `stray_punctuation_not_read` (IN-11) | a | `stray_punctuation` | 0 · 59,128 | 0 · 227,609 | 0 · 108,378 | 0 · 64,245 | 0 · 57,688 |
| `unicode_not_normalised` (IN-4) | a | `unicode_normalised` | 0 · 12,720 | 0 · 12,776 | 0 · 12,964 | 0 · 12,616 | 0 · 12,096 |
| `long_input_not_refused` (IN-8) | a | `bounded_arguments` | 0 · 456 | 0 · 456 | 0 · 456 | 0 · 456 | 0 · 456 |
| `unversioned_tx_case_not_found` (IN-5) | a | `one_case_schemes` | 0 · 306 | 0 · 284 | 0 · 297 | 0 · 296 | 0 · 318 |
| `case_exact_synonym_first` (IN-2) | a | `single_case_symbols` | 0 · 0 | 0 · 150 | 0 · 0 | 0 · 0 | 0 · 0 |
| `synonym_database_offered`: the synonym database accepted as a target | a | `no_synonym_target` | 0 · 1 | 0 · 1 | 0 · 1 | 0 · 1 | 0 · 1 |

Scored again on 2026-11-20 under the judge of today: the suffixed forms that are exactly another gene's name or synonym at a release the truth read only through the revision `e2_names_every_release` are undetermined or right (the cause `suffix_form_is_name_at_unread_release`, 26 · 5 · 10 · 30 · 26 answers before, is gone; a form read through the case or separator step stays `separator_joins_suffix`); the human LRG identity lines are right under `now` (the judge reads LRG ids) and wrong under `decided` (rule `lrg_answers_gene`, not built); mouse's letter-case lines are judged with `e2_case_r1`, their 13 wrong answers the graph's merged spellings; pandas' missing-value spellings (`e1_dec2`) add 12 wrong answers per organism under `decided`. On `d4c2dda` the causes and their counts are the same but `source_inference_reads_names_at_snapshot`, whose lines its api does not have. Every class a cause under `now` is one that a decided rule removes (`infer_duplicates_counted`: rule `inference_distinct_ids`, which the line's relation already expects under both judges, since `71690ab`'s documentation states it); under `decided` the classes a are, by construction, the decided rules not yet built.

## Known limitations

E's rows of [the known-limitations table](../../results/limitations/limitations.tsv), generated by `limits_e` from the `71690ab` runs under the judge `now` (what users get): the causes of the wrong answers and errors, the forms idtrack does not read (rule `safe`: not found), and the open lines that named a gene. Answers of the strata's answers per organism:

| limitation | what | user action | until | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `L-internal-names` | `synonym_id::X` and `<id>.Void` are read as identifiers | do not type them | the decided rule `no_synonym_target` | 14 of 400 | 10 of 550 | 7 of 403 | 17 of 400 | 10 of 400 |
| `L-long-inputs` | a version of more than 4,300 digits raises Python's `ValueError` | none needed (no real input) | the decided rule `bounded_arguments` | 450 of 1,200 | 450 of 1,200 | 451 of 1,201 | 450 of 1,200 | 450 of 1,200 |
| `L-make-unique` | `X-1`, `X.1`, `X_1` are read through the case and separator steps as another gene's name | convert the names before making them unique | the decided rule `made_unique_copies` | 89 of 51,851 | 28 of 231,585 | 34 of 114,999 | 93 of 54,495 | 79 of 48,038 |
| `L-stray-punctuation` | `X;`, `X.`, `'X'` match names of Ensembl's tables that carry the same stray character | strip quotes and trailing punctuation | the decided rule `stray_punctuation` | 27 of 50,151 | 499 of 235,436 | 180 of 118,097 | – | – |
| `L-missing-values` | the rule's spellings of a missing value that answer a gene (single calls) | drop missing values before converting | the decided rule `missing_values` | 35 of 53 | 28 of 53 (and 1 of 7: `_`) | 35 of 53 | 29 of 53 | 23 of 53 |
| `L-infer-names` | lists of names placed at a later release | infer the release from Ensembl ids | the decided rule `inference_per_release` | – | 2 of 3 | 4 of 5 | – | – |
| `L-case-synonym` | in organisms whose symbols are written in one case, a symbol typed in another case that is exactly another gene's synonym answers that gene | type symbols as the organism writes them, or pass `from_database` | none | 2 of 2 | – | – | – | – |
| `L-bare-numbers` | an NCBI gene id of another organism answers a gene through numeric ids of other databases | pass `from_database` for NCBI gene ids | none (rule `bare_numbers`: a number is read in every numeric database, each named) | 0 of 6,000 | 149 of 6,000 | 124 of 6,000 | 0 of 6,000 | 0 of 6,000 |
| `L-foreign-symbols` | a symbol of another organism that no table of this organism names answers a gene | convert symbols in their own organism's graph | none (the truth's reach) | 20 of 4,384 | 26 of 2,413 | 2 of 2,509 | 10 of 4,015 | 31 of 4,211 |
| `L-infer-duplicates` | a list with repeated identifiers can be placed at another release than the same list without them | remove repeated identifiers before inferring | the decided rule `inference_distinct_ids` | – | 1 of 24 | 1 of 24 | 1 of 24 | – |
| `L-format-characters` | bidi marks, isolates and soft hyphens: not found | strip format characters | none | 7,750 of 7,750 | 7,750 of 7,750 | 7,750 of 7,750 | 7,750 of 7,750 | 7,500 of 7,500 |
| `L-control-characters` | NUL, escape sequences, other control characters: not found | strip control characters | none | 10,850 of 10,850 | 10,850 of 10,850 | 10,850 of 10,850 | 10,850 of 10,850 | 10,500 of 10,500 |
| `L-lookalikes` | a Cyrillic or Greek look-alike letter: not found | none | none (by design) | 2,624 of 2,624 | 1,985 of 1,985 | 2,599 of 2,599 | 2,593 of 2,593 | 2,539 of 2,539 |
| `L-emoji-brackets` | an emoji, brackets or a space inside: not found | clean the strings | none (by design) | 11,150 of 11,150 | 11,149 of 11,150 | 11,150 of 11,150 | 11,150 of 11,150 | 10,799 of 10,800 |
| `L-prefixes` | `ensembl:ENSG...`, `UniProtKB:P04637`, `RefSeq:NM_...`: not found | strip the prefix | none | 500 of 500 | 500 of 500 | 500 of 500 | 500 of 500 | 500 of 500 |
| `L-number-forms` | `07157`, `7157.0`, `7.157e+03`, `GeneID:7157`: not found | write NCBI ids as integers | none | 1,600 of 1,600 | 1,600 of 1,600 | 1,600 of 1,600 | 1,600 of 1,600 | 1,600 of 1,600 |
| `L-version-forms` | `ENSG...-5`, `ENSG..._5`, `ENSG....v5`: not found | write `<id>.<version>` | none | 450 of 750 | 450 of 750 | 450 of 750 | 450 of 750 | 450 of 750 |
| `L-mouse-hgnc-names` | a mouse gene whose display name is a human HGNC symbol is not found by that name | use the Ensembl id | none (the graph holds no HGNC symbols for mouse) | – | 2 of 4 | – | – | – |
| `L-spellings` | a string carried in two letter cases (mouse `GC` and `gc`): the graph keeps one spelling per database and release, so one case answers the genes of both, or the symbol matched without letter case | type the exact spelling | decided: rule `spelling_per_edge` | – | 13 of 6,143 | – | – | – |

The rows generated from the `d4c2dda` runs are the same, but `L-infer-names`, whose lines `d4c2dda`'s api does not have.

## Reproducing from scratch

From an empty data root, with the configuration naming the trees `d4c2dda` and `71690ab` and their graphs:

1. **Fetch the truth** as the harness does ([F0's steps](f0.md#reproducing-from-scratch): GTFs, core history tables, BioMart tables of every release, VGNC's sets), and **publish the harness's sets**: E2's bases are the harness's published lines.
2. **Run F0's check run** of the harness's sets on each tree (`<families>/f0/runs/check_<tree>/<o>`): E2 takes its bases from the lines `d4c2dda` answers right there, and E6 the slowest of them.
3. **Publish E's sets**: `gen_e --organism <o>` (every set and its `_dec` revision); `revise_e --organism macaca_fascicularis --revision e5_r1`, then `e5_r2`.
4. **Run and score** on both trees: `run_e`, then `score_e` (both judges, E's rules, `e_summary.tsv`).
5. **The other parts**: `det submit` and `det compare` (E7), `relations_e submit` (E8), and on `71690ab` `budget submit` (E6, after A*'s runs if its lines are to be included), `harm_e` (human, mouse, pig) and `web_e` in the web service's environment (E9).
6. **Classify**: `causes_e` on every run; **compare** each `71690ab` run with `d4c2dda`'s (`compare_vs_d4c2dda_<judge>.tsv`, under both judges).
7. **Limitations**: `limits_e --runs <families>/e/runs/71690ab/<o> ... --out <tsv>` (and the same for `d4c2dda`).
8. **Results**: `idtrack-repro results` writes `results/e/<tree>/<organism>/` and `results/e/71690ab/web/`, the known-limitations table, the index and the provenance.
