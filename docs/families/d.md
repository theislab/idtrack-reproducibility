# D: real-world inputs

Family D asks idtrack what its users ask: every feature of the atlas project's 33 authors' matrices, as the project asks it (D1); the release a list of identifiers comes from (D2); the symbols of published marker lists (D3); and mixed lists of every kind of identifier in one batch call (D4). Its code is `idtrack_reproducibility.families.d`; its data are `<families>/d/` (`queries/<organism>/`; `runs/<run>/<organism>/` with one run folder per part and tree: `d4c2dda` (D1), `d4c2dda_d2`, `d4c2dda_d3`, `d4c2dda_d4`, `71690ab_d2`; `inputs/`: D3's lists with their `SOURCES.tsv`; `truth/`: the check of the chain mapping; `checks/`; `limitations/`); the design is in [the methods](../methods.md#d-real-world-inputs). D1 and D3 are plain conversions run by F0's runner; D2 (`infer_source`) and D4 (the batch call) are asked through family E's adapter ([E](e.md)) by `run_d`. The atlas project's feature tables are read, never changed, from the configuration's `[sources] rosetta` (its checkout) and `rosetta_m02` (its per-feature tables).

## State

What is built:

- the generators (`gen_d1` to `gen_d4`), the fetchers of D's truth and inputs (`fetch_exons`: the merged exons per gene of Ensembl's main GTF of a release; `fetch_grch37`: the gene ids of Ensembl's GRCh37 site at 76 to 115; `fetch_d3`: the marker lists and the dates of Ensembl's releases), the survey of the matrices (`survey_m02`), the error bound of the chain mapping (`check_map`), the revision sets (`revise_d1`, `revise_d2`, `revise_d4`), the runner of D2 and D4 (`run_d`), the reports (`report_d1` to `report_d4`, `bridges_d1`, `causes_d4_from_single`), the direct check of missing values on `71690ab` (`check_missing_71690ab`), the plugin of the common classifier (`causes_d`) and the known-limitations builder (`limits_d` with its map `limitations_map_d.tsv`);
- the truth files D fetched, each in its organism's `SOURCES.tsv`: gene models `exons/<R>.tsv.gz` of human 68, 75, 76, 77, 84, 87, 89, 93, 95, 96, 98, 101 and 115, mouse 71, 77, 84, 87, 92, 97, 98, 99, 102, 103 and 115, pig 101 and 115; the gene tables of Ensembl's GRCh37 site of releases 76 to 115 (`grch37_site/genes_<R>.tsv.gz`); D3's inputs (`<families>/d/inputs/SOURCES.tsv`): PanglaoDB (its file of 2020-03-27), CellMarker 2.0, MSigDB 2024.1 (`c8`, `m8`), Ensembl's release dates;
- the published sets ([the sets](#the-sets)).

What was run and scored (the graphs of the reference tree `d4c2dda`, which `71690ab` also reads):

| run folder | part | organisms and sets scored | causes | D's own tables |
| --- | --- | --- | --- | --- |
| `runs/d4c2dda/` | D1 | pig, mouse, human: `d1_v1`, `d1_r1`, `d1_r2` | yes | `d1/` |
| `runs/d4c2dda_d2/` | D2 | pig and mouse: `d2_v1`; human: `d2_v1`, `d2_r1`, `d2_r2`, `d2_r3` | mouse, human (pig: no wrong answer) | `d2/` |
| `runs/d4c2dda_d3/` | D3 | mouse, human: `d3_v1` | yes | `d3/` |
| `runs/d4c2dda_d4/` | D4 | all five: `d4_v1`, `d4_r1`, `d4_r2` | yes (cynomolgus: no wrong answer) | `d4/` |
| `runs/71690ab_d2/` | D2 on `71690ab` | pig and mouse: `d2_v1`, `d2_v2`, `d2_r1`; human: `d2_v2` | yes | – |

Every one of these folders holds `run.json` and `report.tsv`. Not run: `runs/71690ab_d4/` (D4 on `71690ab`, stopped before it started; the question it would answer for missing values is checked directly: `<families>/d/checks/missing_71690ab_pig.json`).

What is not built: `d2_v2` on `d4c2dda` and human's `d2_v1` on `71690ab`; `infer_source` on D4's mixed lists (a list of several origins has no one expected origin; D2 measures the function on lists of one origin); D3 for pig and the macaques (the lists hold human and mouse only); D4's items from families B and C (added as `d4_v2` once their sets are run); the family's panel ([Measurement design](../methods.md#measurement-design)).

## What it measures

- **D1, the authors' matrices**: whether idtrack's answer for a feature, asked as the atlas project asks it (the feature reduced as the project reduces it, from the source release the project uses, to 115), is the gene of 115 at the place of the feature's gene model; by lines and weighted by the counts the features hold in the authors' matrices.
- **D2, source inference**: whether `infer_source` in release mode answers the releases whose GTFs fit a list best, on synthetic lists of every release and on the matrices' lists.
- **D3, marker lists**: every symbol of PanglaoDB, CellMarker 2.0, MSigDB's cell type signatures and the atlas project's lineage table, to 115 from the release of the list's date and without a release.
- **D4, mixed lists**: every item of a list of 5 to 100 identifiers of every kind, asked in one batch call, judged as its single call is judged, and compared with that single call.

The rate of right answers is right / (right + wrong + error) (an expected error counts as right; undetermined answers and timeouts apart), per stratum and organism, with Wilson intervals and a cluster bootstrap over the unit (D1: the feature's gene or symbol; D2: the release drawn or the dataset; D3: the symbol; D4: the item's identifier). D1, D3 and D4 lines are scored under `best` and under `all` (`strategy: both`), D2 lines under `best`.

## The sets

One set per part and organism, `<part>_v<n>`, published once (F0's `lines.publish`); a change of truth is a revision set under the published qids.

| set | strata | what is asked | truth |
| --- | --- | --- | --- |
| `d1_v1` | `d1_<id\|symbol>_<same\|cross>_<one\|dominant\|open>` (`same`: the source release lies on the target's assembly; `cross`: GRCh37 or GRCm38), `d1_id_versioned_<one\|dominant\|open>` | every feature of the 33 matrices (1,209,403 author features of 37 sources; `<families>/d/survey_m02.tsv`), reduced as the project reduces it (its suffix cut, a version removed, an Ensembl id or a symbol by the declared type), from the source release R (stated by the authors or measured by the project) to 115; an Ensembl id the GTF of R does not hold is asked without a release; the matrix whose ids carry versions (`xiang2020_humanpregastrula`) is also asked with them | the place ([the truth](#truth)) |
| `d2_v1` | `d2_synth_<size>_<noise>`, `d2_matrix_<ids\|symbols\|mixed>_<stated\|measured>` | `infer_source` (release mode) on synthetic lists of the kept gene ids (Cell Ranger's biotypes) of every release: 50, 500, 5,000 and all ids, with no noise and with 5 % and 20 % of another release's ids or of junk, 5 seeds (lists of all ids: one draw); and on each matrix source's list (all features; those with counts) | the releases whose GTF holds the largest share of the list (a tie is right only as the whole set, rule `set_exact`); a stated release outside that set is `open` |
| `d2_v2` | `d2_synth_names_<500\|all>`, `d2_synth_versioned_500` | the same on lists of one release's display names (500 names, 5 seeds; all names) and of versioned ids (500, 5 seeds; releases whose GTF writes versions) | the releases whose GTF holds the largest share of the names (of the id.version pairs) |
| `d3_v1` | `d3_<panglaodb\|cellmarker\|msigdb\|rosetta_lineage>_<from\|default>` | every distinct symbol of PanglaoDB (2020-03-27), CellMarker 2.0 (each row at its paper's year), MSigDB 2024.1 `c8`/`m8` and the atlas project's lineage table (115), from the release current at the list's date and without a release | Ensembl's names at the release read, walked to 115; HGNC and MGI beside it |
| `d4_v1` | `d4_<category>` (`id_current`, `id_old`, `transcript`, `symbol`, `ncbi`, `refseq`, `uniprot`, `nothing`, `missing`), `d4_issue85` | 2,000 seeded lists per organism of 5, 10, 20, 50 or 100 items drawn from the accuracy harness's published lines (with probability 0.3 one to three items repeated, and with probability 0.3 one or two missing values inserted), each item a line of one batch call to 115; human also idtrack issue #85's list (`TP53`, `ENSG00000141510`, `NOT_A_REAL_GENE`) | each item keeps its harness line's truth whole; a missing value expects rule [`missing_values`](../rules.md#missing_values)'s error |

Lines per set and organism (manifests `<families>/d/queries/<organism>/<set>.manifest.json`):

| set | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `d1_v1` | 58,735 | 306,473 | 816,718 | – | – |
| `d1_r1` | 37 | 134 | 3,268 | – | – |
| `d1_r2` | 2 | 223 | 1,692 | – | – |
| `d2_v1` | 4,804 | 5,457 | 5,490 | – | – |
| `d2_v2` | 555 | 603 | 603 | – | – |
| `d2_r1` | 65 | 5 | 2,263 | – | – |
| `d2_r2` | – | – | 1,338 | – | – |
| `d2_r3` | – | – | 31 | – | – |
| `d3_v1` | – | 259,104 | 414,692 | – | – |
| `d4_v1` | 76,164 | 76,265 | 78,141 | 76,391 | 76,803 |
| `d4_r1` | 9,366 | 10,338 | 9,143 | 10,356 | 10,009 |
| `d4_r2` | 610 | 566 | 583 | 646 | 553 |

The revision sets (only the lines whose truth changed, under their published qids):

| set | why |
| --- | --- |
| `d1_r1` | a symbol that several genes carry as their name at R is judged `set` (`best` may answer any of them, `all` must return them all), not `set_exact`: the matrix's feature is one of those genes and which one cannot be told from the name |
| `d1_r2` | the `mutual` branch of the dominance rule chose a small new gene inside a gene that had grown; the chosen gene must also share the most of the model's bases |
| `d2_r1` | human: release mode also counts Ensembl's GRCh37 site served with releases 76 to 115, which idtrack's graph holds as assembly 37; pig and mouse (lines of `d2_v2`): a list of versioned ids whose stable ids the releases before 77 hold as well as the expected releases hold its versions is `open` (those GTFs write no versions) |
| `d2_r2` | human: every line of `d2_r1` with the items of its first publication (the revision had redrawn the noise of 1,289 lines) and the truth computed on them |
| `d2_r3` | human: the core genes the GTFs of 48 to 51 lack (A*'s truth) count for those releases |
| `d4_r1` | both trees document that a batch holding an empty or blank string is refused whole: every item of such a list expects that error (its own truth kept as `expect_alone`) |
| `d4_r2` | the missing-value strings `nan`, `NA`, `None`, `null` expect rule `missing_values`'s error, which `d4c2dda` and `71690ab` do not have: their lines are `open`, the expected error kept, as family E's `now` judge scores them |

## Truth

Nothing of idtrack.

- **D1, the place** (`gen_d1`): the feature's gene model in Ensembl's GTF of R (merged exons; a symbol's model: the gene or genes whose `gene_name` it is at R, exactly, else without letter case) against the models of every gene of 115 (main GTF) by exon bases shared on the strand. A model of another assembly (human GRCh37 at 68 and 75, mouse GRCm38 at 71 to 99) is mapped exon by exon through Ensembl's chain files, the aligned pieces on the sequence and strand that take most bases kept, unmapped below half of the model's bases. The atlas project's thresholds: a model holds a gene from 1 shared base; the one gene it holds is its gene from half of the model's bases. Classes: `one` (one gene, at least half of the model: the strict truth); `dominant` (several genes, one told by shared bases, or by holding half of the model and half of its own bases with no fewer of the model's bases than any other; its own stratum); open (`several`, `below_share`, `none`, `unmapped`, an id only on a patch at R, an id or a name R does not hold, a made-unique name `X.1`), judged with F0's `open_within` (nothing or a gene at the place is undetermined, a gene elsewhere wrong). Every id line also carries Ensembl's history walk R to 115 (the strict walk with living ids) as the other reading of Ensembl's own data, never the judge. The error bound of the chain mapping (`check_map`, every id present on both sides of a handoff judged as a feature): human 75 to 76, 57,392 genes, 55,110 mapped (0.9602), the truth expects another gene for 37 (17 `dominant`, 20 `one`); mouse 102 to 103, 52,404 genes, 52,035 mapped (0.9930), another gene for 5 (1 and 4) (`<families>/d/truth/check_map_homo_sapiens_75_76.tsv`, `check_map_mus_musculus_102_103.tsv`).
- **D2** (`gen_d2`): per release r of the window, the share of the list's items that are genes of r (an id: a gene id of r's GTF, both regions; for human from 76 also of Ensembl's GRCh37 site, and the core genes of human 48 to 51 and cynomolgus 112; a symbol: a display name of r, without letter case for a mouse table written in capitals); the expected answer is the set of releases of the largest share.
- **D3** (`gen_d3`): the genes whose GTF display name the symbol is at the release read (primary region; exactly, else without letter case), walked to 115 by the strict walk with living ids; without a release, read as rule [`start_release`](../rules.md#start_release) reads it. A symbol that is no display name at the release read is `open`; one whose genes Ensembl retired without successor is `open` within the genes that carry it at 115. Beside it, HGNC (approved, previous and alias symbols, withdrawn symbols merged into an entry) and MGI (current and withdrawn symbols, MGI's Ensembl ids) as the line's `authority`.
- **D4** (`gen_d4`): each item is a published line of the accuracy harness (the latest line of its qid, the corrected judge) asked without a release to 115, and keeps that line's truth whole (`expect`, `broad`, `rule`, `cmp`, `apart`, `locus`, `support`, `spellings`, ...), so F0's scorer judges the item exactly as the harness judges the single call (`src_qid` names it).

### Inputs from the atlas project

D1 to D3 read files of the atlas project, which are neither Ensembl's nor an authority's: its per-feature locus tables of the 33 matrices (`m02/<dataset>.tsv.gz`: every author feature, its counts, the source release the project used; `m16/<dataset>.tsv.gz`: its place step, where it applies), its datasets' preparation configs (`stages/data_preparation/<dataset>/configs/prepare.toml`) and its lineage gene table (`stages/reference_data/gene_sets/lineage_genes.tsv`). They are external inputs: `idtrack-repro family d inputs_atlas` copies them from the configured `[sources]` (`rosetta`: the project's checkout; `rosetta_m02`: its tables) into `<families>/d/inputs/atlas/` (`m02/`, `m16/`, `rosetta/<path>`) and writes `inputs/atlas/SOURCES.tsv` with each file's origin, the checkout's commit, bytes and sha256; the data manifest of `families/d` lists them. The copies of record (150 files, 42.5 MB) were taken on 2026-11-19 from the project's checkout at `d95b5c8` (the configs and the lineage table last changed in its commit `b31eefa` of 2026-11-18) and from its tables written on 2026-11-18, both before D's sets were published (2026-11-19), so they are the files D read. A reader obtains them from the data root's `families/d` record; with `[sources]` naming the copies (`rosetta = "<families>/d/inputs/atlas/rosetta"`, `rosetta_m02 = "<families>/d/inputs/atlas/m02"`), D's sets are generated again from the data root alone.

## Rules

The rules of idtrack the expectations rest on: [`start_release`](../rules.md#start_release) (D1's ids that R does not hold, D3's default lines, D4), [`own_id_first`](../rules.md#own_id_first), [`ties`](../rules.md#ties), [`names_before_synonyms`](../rules.md#names_before_synonyms), [`external_bridge`](../rules.md#external_bridge), [`copies_apart`](../rules.md#copies_apart), [`typed_api`](../rules.md#typed_api) (`71690ab`'s `MissingValueError` for a blank item); the judge's [`id_lives_on`](../rules.md#id_lives_on), [`resurrection`](../rules.md#resurrection), [`name_before_synonym`](../rules.md#name_before_synonym) and [`judge_corrections`](../rules.md#judge_corrections) (D4's items). The decided rules D measures before they land: [`bridge_same_place`](../rules.md#bridge_same_place) (D1's bridges measured by place), [`place_moves_marked`](../rules.md#place_moves_marked) (answers that follow the history to another place), [`main_line_first`](../rules.md#main_line_first), [`made_unique_copies`](../rules.md#made_unique_copies) (D1's `X.1` features), [`inference_per_release`](../rules.md#inference_per_release) (D2), [`missing_values`](../rules.md#missing_values) (D4). Scoring is F0's ([scoring and outcomes](f0.md#scoring-and-outcomes)), with `open_within` for D1's open lines.

## Running it

Locations, trees, graphs and the job runner come from the configuration ([running F0](f0.md#running-it)), and the atlas project's inputs from `[sources] rosetta` and `rosetta_m02` ([their copies](#inputs-from-the-atlas-project)); `<t>` is `d4c2dda` or `71690ab` (or a later tree). `--sets` lists the sets in publication order: a later set replaces the lines of its qids.

```
# truth and inputs (once; one connection per host)
python -m idtrack_reproducibility.families.d.fetch_exons <o>:<release>[,<release>...] ...
python -m idtrack_reproducibility.families.d.fetch_grch37 [first last]
python -m idtrack_reproducibility.families.d.fetch_d3
python -m idtrack_reproducibility.families.d.survey_m02 <families>/d/survey_m02.tsv
python -m idtrack_reproducibility.families.d.check_map <o> <last release of the old assembly> <first of the new> <out tsv>
# the sets (published once) and their revisions
python -m idtrack_reproducibility.families.d.gen_d1 --organism <o> --publish d1_v1 [--dry]
python -m idtrack_reproducibility.families.d.gen_d2 --organism <o> --publish d2_v1 [--parts synthetic,matrices] [--dry]
python -m idtrack_reproducibility.families.d.gen_d2 --organism <o> --publish d2_v2 --parts more
python -m idtrack_reproducibility.families.d.gen_d3 --organism <o> --publish d3_v1 [--dry]
python -m idtrack_reproducibility.families.d.gen_d4 --organism <o> --publish d4_v1 [--lists 2000] [--dry]
python -m idtrack_reproducibility.families.d.revise_d1 --organism <o> --revision <revision> --set <new set> [--dry]
python -m idtrack_reproducibility.families.d.revise_d2 --organism <o> --set <new set> [--dry]
python -m idtrack_reproducibility.families.d.revise_d4 --organism <o> --set <new set> --revision <revision> [--dry]
# D1, D3: F0's runner and scorer
python -m idtrack_reproducibility.families.runner --family d --organism <o> --clone <source> --graph <graph> \
    --out <families>/d/runs/<t>/<o> --sets d1_v1,d1_r1,d1_r2 --submit <n> --score [--tag d1_<t>_<o>]
python -m idtrack_reproducibility.families.runner --family d ... --out <families>/d/runs/<t>_d3/<o> --sets d3_v1 --submit <n> --score
# D2, D4: F0's runner with family E's adapter, then F0's scoring of the stored answers
python -m idtrack_reproducibility.families.d.run_d --organism <o> --clone <source> --graph <graph> \
    --out <families>/d/runs/<t>_d2/<o> --sets d2_v1,d2_v2 --submit <n> [--tag d2_<t>_<o>]
python -m idtrack_reproducibility.families.runner --family d --organism <o> --clone <source> --graph <graph> \
    --out <families>/d/runs/<t>_d2/<o> --sets d2_v1,d2_v2,d2_r1,d2_r2,d2_r3 --submit 0 --score
python -m idtrack_reproducibility.families.d.run_d ... --out <families>/d/runs/<t>_d4/<o> --sets d4_v1 --submit <n>
python -m idtrack_reproducibility.families.runner --family d ... --out <families>/d/runs/<t>_d4/<o> --sets d4_v1,d4_r1,d4_r2 --submit 0 --score
# D's own tables
python -m idtrack_reproducibility.families.d.report_d1 <families>/d/runs/<t>/<o> --organism <o> --sets d1_v1,d1_r1,d1_r2
python -m idtrack_reproducibility.families.d.bridges_d1 <families>/d/runs/<t>/<o> --organism <o> --sets d1_v1,d1_r1,d1_r2
python -m idtrack_reproducibility.families.d.report_d2 <families>/d/runs/<t>_d2/<o> --organism <o> --sets <sets>
python -m idtrack_reproducibility.families.d.report_d3 <families>/d/runs/<t>_d3/<o> --organism <o>
python -m idtrack_reproducibility.families.d.report_d4 <families>/d/runs/<t>_d4/<o> --organism <o> \
    --single <families>/f0/runs/check_<t>/<o> --sets d4_v1,d4_r1,d4_r2
# causes: explained on the tree that gave them, classified with D's plugin; D4's also by their single call's cause
<python of the tree> -m idtrack_reproducibility.families.explain --family d --organism <o> --clone <source> --graph <graph> \
    --scored <run folder> --per-unit 3
python -m idtrack_reproducibility.families.causes --organism <o> --scored <run folder> \
    --plugin src/idtrack_reproducibility/families/d/causes_d.py
python -m idtrack_reproducibility.families.d.causes_d4_from_single <families>/d/runs/<t>_d4/<o> --organism <o> \
    --single <families>/f0/runs/check_<t>/<o>
# missing values on a tree from 71690ab, item by item (the tree's Python)
<python of the tree> -m idtrack_reproducibility.families.d.check_missing_71690ab <source> <graph> <out json>
# a later tree against the reference, and the known-limitations rows
python -m idtrack_reproducibility.families.compare <families>/d/runs/d4c2dda<part>/<o> <families>/d/runs/<t><part>/<o>
python -m idtrack_reproducibility.families.d.limits_d --out <tsv> <run folder> [...]
```

## Results

Every number below is copied from the run's `report.tsv` (row `*`: both views; row `best`), its `causes_summary.tsv`, D's own tables in `<families>/d/runs/<run>/<organism>/` or the file named; the result files are `results/d/<run>/<organism>/` with the same names (D's own tables without their subfolder). No answer timed out, and no answer raised an error.

### D1: the authors' matrices on `d4c2dda`

Runs `<families>/d/runs/d4c2dda/<organism>/` (pig, mouse, human); result files `results/d/d4c2dda/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv`, `rates.tsv`, `datasets.tsv`, `history_vs_locus.tsv`, `wrong_classes.tsv`, `bridges_no_successor.tsv`. The graphs of `d4c2dda` are the graphs the atlas project maps with, so D1 measures exactly what the project receives. Rate with the cluster-bootstrap interval by lines (row `best` of `report.tsv`) and the rate weighted by the features' counts (`rate_counts`, row `best` of `rates.tsv`); lines in parentheses:

| stratum | pig | mouse | human |
| --- | ---: | ---: | ---: |
| ids, one gene at the place, same assembly | 0.9511 (0.9480–0.9541); 0.99763 (27,251) | – | 0.9948 (0.9942–0.9954); 0.99985 (374,482) |
| ids, one gene, GRCh37 or GRCm38 | – | 0.9710 (0.9698–0.9722); 0.99565 (185,872) | 0.9783 (0.9765–0.9800); 0.99814 (27,798) |
| ids with their version, one gene | – | – | 0.9951 (0.9945–0.9957); 0.99984 (49,555) |
| ids, a dominant gene, same assembly | 0.8930 (0.8746–0.9095); 0.96804 (1,608) | – | 0.9600 (0.9563–0.9635); 0.99054 (51,647) |
| ids, a dominant gene, GRCh37 or GRCm38 | – | 0.9312 (0.9244–0.9379); 0.99516 (13,297) | 0.9339 (0.9260–0.9412); 0.99551 (3,689) |
| ids with their version, a dominant gene | – | – | 0.9111 (0.9047–0.9177); 0.98784 (7,798) |
| symbols, one gene, same assembly | 0.9931 (0.9917–0.9945); 0.99931 (12,608) | – | 0.9953 (0.9948–0.9958); 0.99985 (122,245) |
| symbols, one gene, GRCh37 or GRCm38 | – | 0.9931 (0.9925–0.9937); 0.99558 (88,793) | 0.9807 (0.9792–0.9821); 0.99776 (93,440) |
| symbols, a dominant gene, same assembly | 0.9022 (0.8804–0.9217); 0.89916 (920) | – | 0.9263 (0.9209–0.9312); 0.99034 (19,167) |
| symbols, a dominant gene, GRCh37 or GRCm38 | – | 0.9538 (0.9468–0.9604); 0.99352 (7,041) | 0.9038 (0.8958–0.9115); 0.98487 (12,959) |

The open lines (no one gene at the place) answered by a gene outside the place, `best` (wrong / lines): pig ids 679 / 13,368, symbols 24 / 2,980; mouse ids 285 / 6,819, symbols 60 / 4,651; human ids 141 / 19,037 (same assembly), 45 / 1,251 (GRCh37), 28 / 7,439 (versioned), symbols 57 / 7,222 (same assembly), 324 / 18,989 (GRCh37). The `all` view is lower on symbols, as it adds the genes that carry the name in another source (`report.tsv`, row `all`). `datasets.tsv` gives the same per matrix.

**The place against Ensembl's history** (`history_vs_locus.tsv`, the `one` lines of ids): the strict walk names exactly the place's gene for pig 25,596 of 27,251, mouse 179,383 of 185,872, human 370,844 of 374,482 (same assembly) and 26,648 of 27,798 (GRCh37); it has no successor for 419, 4,346, 3,342 and 1,070; it names several genes, the place's among them, for 1,128, 1,791, 148 and 45; it names another gene for 108, 352, 148 and 35.

**Wrong answers by cause** (`causes_summary.tsv`, F0's classifier with D's plugin; answers summed over strata, both views):

| cause | class | pig | mouse | human |
| --- | --- | ---: | ---: | ---: |
| `d1_history_to_another_place`: the answer is Ensembl's history walk, which leads to genes other than the place's | b | 437 | 1,434 | 5,790 |
| `d1_name_history_to_another_place`: the same for a symbol's genes | b | 161 | 613 | 4,256 |
| `d1_open_follows_history`: an open line answered by the history to a gene outside the place | b | 1,224 | 166 | 44 |
| `d1_no_successor_nothing`: Ensembl's history retires the id, idtrack answers nothing, the place holds a gene of 115 | b | 516 | 7,098 | 4,752 |
| `d1_name_no_successor_nothing`: the same for a symbol's genes | b | 8 | 422 | 3,442 |
| `d1_history_split_place_one`: the history splits, the place holds one of the successors (`all` returns them all) | b | 1,498 | 2,157 | 209 |
| `d1_history_part_not_place`: `best` answers a successor of a split away from the place | b | 658 | 1,459 | 57 |
| `external_bridge`: no successor in Ensembl's history, a bridge to a gene elsewhere | b | 209 | 1,140 | 1,403 |
| `rfam_family_label`: an RFAM family label names the gene | b | – | – | 78 |
| `d1_name_in_other_source`: the answered gene carries the symbol in another name source or as a synonym at R | e | 124 | 1,338 | 3,838 |
| `name_in_other_source`: the same, by F0's generic test | e | – | – | 2,236 |
| `d1_name_not_a_name_at_source`: a gene whose name the symbol is not at R, while the GTF names another gene so | e | 38 | 245 | 86 |
| `path_through_resurrection`: an id retired and brought back, read through its return | e | – | 7 | – |
| `d1_name_not_found`: the symbol is not found | c | 128 | 4 | – |
| `d1_name_several_part`: a part of the genes a name of several genes names | c | – | 16 | 53 |
| `synonym_spelling_merged`: two genes carry the string in two letter cases | c | 1 | 89 | 12 |
| `mirbase_not_in_graph`: miRBase names are not in the graph | c | – | – | 20 |
| `name_in_db_not_in_truth`: a name in a database the truth does not read | d | – | 25 | 4 |
| `unexplained` | – | 1 | 8 | 37 |

**No wrong answer leaves out a successor of Ensembl's history** (the cause `d1_successor_dropped`, class a, has no row): every wrong answer either follows Ensembl's history where the gene's model moved to another place or id (class b, measured by the decided rules `place_moves_marked` and `bridge_same_place`), or reads a name that another source gives the gene.

**Bridges measured by place** (`bridges_no_successor.tsv`, row `all strata`, `best`): where Ensembl's history has no successor and the old model's place holds one gene of 115, `d4c2dda`'s `best` answer is that gene (a bridge that lands right) for pig 168, mouse 1,025, human 2,491 lines, another gene for 34, 324 and 445, nothing for 258, 3,549 and 2,376. A bridge rule that keeps only bridges to a gene at the old place (the decided rule `bridge_same_place`) keeps the first and drops the second.

### D2: source inference

Runs `<families>/d/runs/d4c2dda_d2/<organism>/` and `<families>/d/runs/71690ab_d2/<organism>/`; result files `results/d/d4c2dda_d2/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv` (mouse, human), `synthetic.tsv`, `matrices.tsv`, and `results/d/71690ab_d2/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv`.

- **Synthetic lists of ids** (`d2_v1`): every line right on `d4c2dda` (`synthetic.tsv`, row `*`: pig 4,800 of 4,800, mouse 5,437 of 5,437, human 5,440 of 5,440) and on `71690ab` for pig and mouse (every `d2_synth_*` stratum of `report.tsv`): the inferred releases equal the GTFs' tie set in every list, with and without noise.
- **The matrices' lists** (`matrices.tsv`): every id list right but mouse's `pijuansala2019_gastrulationatlas`, which is undetermined twice (stated 92; its ids are genes of 89 only); the symbol lists are wrong but for two human matrices with a stated release (`d2_matrix_symbols_measured`: human 18 wrong of 18, mouse 6 of 6; `d2_matrix_symbols_stated`: human 6 right of 6, mouse 2 wrong of 2): `infer_source` reads every name at the graph's last release, so a name that changed is not counted and a name new since counts only at 115, and the answer is mostly 115.
- **Lists of display names and of versioned ids** (`d2_v2`, on `71690ab`): lists of one release's names are right for pig 162 of 300 (500 names) and 48 of 60 (all names), mouse 14 of 340 and 2 of 68, human 76 of 340 and 19 of 68; lists of 500 versioned ids are right for every one the truth decides (pig 130, mouse 190, human 195; undetermined pig 65, mouse 5).

Wrong answers by cause: every one is `d2_names_read_at_last_release` (class b: `infer_source`'s documented reading of names): `d4c2dda` mouse 8, human 18; `71690ab` pig 150, mouse 400, human 313.

### D3: marker lists on `d4c2dda`

Runs `<families>/d/runs/d4c2dda_d3/<organism>/` (mouse, human); result files `results/d/d4c2dda_d3/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv`, `authority.tsv`. Rate and cluster-bootstrap interval, wrong and undetermined answers (row `*`):

| list | human, from the list's release | human, without a release | mouse, from the list's release | mouse, without a release |
| --- | ---: | ---: | ---: | ---: |
| PanglaoDB | 0.9974 (0.9963–0.9985); 22 / 94 | 0.9971 (0.9959–0.9983); 25 / 2 | 0.9859 (0.9824–0.9892); 123 / 74 | 0.9848 (0.9811–0.9880); 134 / 0 |
| CellMarker 2.0 | 0.9955 (0.9950–0.9960); 1,181 / 100,903 | 0.9878 (0.9868–0.9887); 4,205 / 19,620 | 0.9959 (0.9953–0.9965); 808 / 30,819 | 0.9932 (0.9923–0.9940); 1,427 / 19,736 |
| MSigDB `c8` / `m8` | 0.9953 (0.9945–0.9960); 190 / 542 | 0.9953 (0.9946–0.9960); 189 / 534 | 0.9957 (0.9948–0.9966); 92 / 0 | 0.9953 (0.9943–0.9963); 101 / 0 |
| the atlas project's lineage table | 0.9977 (0.9931–1.0000); 1 / 0 | 0.9977 (0.9931–1.0000); 1 / 0 | – | – |

The undetermined answers are strings that are no Ensembl display name at the release read (CellMarker's aliases and papers' spellings). Wrong answers by cause (human, mouse): the answered gene carries the symbol in another name source or as a synonym (class e; `name_in_other_source` 2,382 and 1,357, `d3_name_in_other_source` 2,063 and 1,093), read without a release at another release than Ensembl's names (`d3_read_elsewhere`, e) 1,082 and 57, spellings merged (`synonym_spelling_merged`, c) 125 and 135, RFAM labels (`rfam_family_label`, b) 71 and 6, not found (`d3_nothing`, c) 44 and 6, names in a database the truth does not read (d) 32 and 28, part of a name of several genes (`d3_part_of_several`, c) 15 and –, `path_through_resurrection` (e) – and 2, unexplained – and 1. Beside the truth (`authority.tsv`, `best`): on the CellMarker lines Ensembl's names leave open from the list's release, human, idtrack returns HGNC's gene for 7,713, another gene for 312, nothing for 25.

### D4: mixed lists on `d4c2dda`

Runs `<families>/d/runs/d4c2dda_d4/<organism>/` (all five); result files `results/d/d4c2dda_d4/<organism>/report.tsv`, `scores.tsv`, `causes_summary.tsv` (all but cynomolgus), `batch_vs_single.tsv`, `causes_from_single.tsv`. The items of the lists refused for a blank string count as their documented error (`expected_error`). Rate and cluster-bootstrap interval, wrong answers (row `*`):

| category | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `d4_symbol` | 0.9876 (0.9847–0.9904); 230 | 0.9999 (0.9997–1.0000); 2 | 0.9996 (0.9992–0.9999); 7 | 0.9910 (0.9887–0.9931); 175 | 1.0000 (1.0000–1.0000); 0 |
| `d4_id_old` | 0.9971 (0.9957–0.9984); 49 | 0.9852 (0.9816–0.9887); 265 | 0.9912 (0.9888–0.9933); 165 | 1.0000 (1.0000–1.0000); 0 | 1.0000 (1.0000–1.0000); 0 |

Every other category (`d4_id_current`, `d4_transcript`, `d4_ncbi`, `d4_refseq`, `d4_uniprot`, `d4_nothing`) is right on every answer in every organism, as is issue #85's list (human, 6 of 6). The missing values (`d4_missing`) are their expected error where their list holds a blank (pig 856, mouse 898, human 760, rhesus 862, cynomolgus 878 answers) and undetermined otherwise (944, 860, 946, 1,000, 846: `d4c2dda` reads `nan`, `NA`, `None`, `null` as identifiers and answers nothing). Wrong answers by cause: bridges where Ensembl has no successor (`external_bridge`, b) pig 49, mouse 265, human 165; miRBase names not in the graph (`mirbase_not_in_graph`, `d4_mirbase_form`, c) pig 118 + 108, rhesus 84 + 68; a separator variant read as another gene's synonym (`d4_separator_synonym`, e) rhesus 21; RFAM labels (b) human 7; shared names (`support_ranks_shared_name`, b) rhesus 2; spellings (c) pig 2, mouse 2; another name source (e) pig 1; unexplained pig 1.

**A list answers as its items answer alone** (`batch_vs_single.tsv`: every item against its harness line asked alone on the same tree, F0's check run): no item is right alone and not right in the list (the table has no `right_alone_only` row); the items whose outcome differs with the same answer (`outcome_differs_same_answer`, `best`: human 55, mouse 9) are judged by later revisions of the harness's truth in D4 than in the check run.

On `71690ab` (`<families>/d/checks/missing_71690ab_pig.json`, pig): a batch call answers `nan` with `ENSSSCG00000029419.3` and `NA` with `ENSSSCG00000031868.3` (genes with those strings as synonyms), `None` and `null` with nothing, and refuses a list holding `""` or `" "` whole with `MissingValueError`.

## Known limitations

D's rows of the known-limitations table come from its causes (the wrong answers, by `limitations_map_d.tsv`) and from the judge's reasons of the undetermined answers (`limits_d`); the stored tables are `<families>/d/limitations/limitations_<part>_<tree>.tsv`, and the rows become D's rows of `results/limitations/limitations.tsv` once D's builder is registered. The wrong-answer limitations D adds, with their sizes in the cause tables above: `L-history-vs-place` (idtrack follows Ensembl's history where the gene's model moved to another id), `L-new-id-no-event` (Ensembl gave the place a new id without an event), `L-split-place` and `L-split-best` (splits whose old place holds one successor), `L-infer-names` (`infer_source` cannot tell a symbol list's release), `L-default-release-names` (a marker read without a release at another release than Ensembl's names). The undetermined answers by the judge's reason, answers per organism:

| limitation | what | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `L-place-several` | D1: the feature's model holds several genes of 115, none dominant | 2,695 | 6,243 | 25,118 | – | – |
| `L-place-empty` | D1: no gene of 115 at the old model's place | 18,992 | 4,952 | 3,874 | – | – |
| `L-place-shrunk` | D1: the one gene at the place holds under half of the old model | 6,318 | 1,095 | 2,097 | – | – |
| `L-place-unmapped` | D1: the old model does not map to the target assembly | – | 4,624 | 2,088 | – | – |
| `L-patch-source` | D1: an id of the source release only on a patch or alternative locus | – | – | 34,531 | – | – |
| `L-id-not-at-source` | D1: an id the source release does not hold | – | 130 | – | – | – |
| `L-name-not-at-source` | D1: a symbol no gene carries as its Ensembl name at the source release (another annotation's names) | 782 | 4,884 | 30,494 | – | – |
| `L-made-unique` | D1: a name made unique by the matrix's maker (`X.1`, `X-1`; the decided rule `made_unique_copies`) | 2,496 | 166 | 1,446 | – | – |
| `L-name-several-place` | D1: a name of several genes, not all with one gene at 115 | – | 122 | 6,810 | – | – |
| `L-stated-release` | D2: a stated release that fits the list worse than another (`d4c2dda` and `71690ab` each) | – | 2 | – | – | – |
| `T-versionless-releases` | D2 on `71690ab`: versions at releases whose GTF writes none (the truth's) | 65 | 5 | – | – | – |
| `L-marker-not-a-name` | D3: a marker string that is no Ensembl name at the list's release | – | 49,198 | 119,650 | – | – |
| `L-marker-retired` | D3: a marker whose gene Ensembl retired | – | 1,431 | 2,045 | – | – |
| `L-missing-values` | D4: missing-value strings read as identifiers (until the decided rule `missing_values`) | 944 | 860 | 946 | 1,000 | 846 |
| `L-no-chain` | D4: no successor and no assembly chain | 1,968 | 68 | 78 | 1,114 | 1,012 |
| `L-bridge-at-place` | D4: no successor; a gene at the same place (until `bridge_same_place`) | 149 | 1,217 | 1,451 | – | – |
| `L-branch` | D4: answers on branches that left the id (until `main_line_first`) | 132 | 71 | 77 | 14 | 33 |
| `L-shared-names` | D4: a name of several genes, `best` keeps the better-supported | 1 | 8 | 36 | 15 | 2 |
| `L-unplaced-copy` | D4: a copy on a haplotype or patch without a place on the reference | – | – | 2 | – | – |

## Reproducing from scratch

From an empty data root, with the configuration naming the trees `d4c2dda` and `71690ab`, their graphs, and the atlas project's `[sources]` (the copies of [the inputs](#inputs-from-the-atlas-project), or the project itself followed by `inputs_atlas`):

1. **Fetch the truth**: the accuracy harness's and A*'s fetches ([A*'s steps](astar.md#reproducing-from-scratch)), then `fetch_exons` for every source release of the matrices and 115, `fetch_grch37`, `fetch_d3`; `survey_m02`; `check_map` at human 75 to 76 and mouse 102 to 103. The accuracy harness's sets are published and F0's check runs exist (`<families>/f0/runs/check_<tree>/<organism>`), as D4 draws its items from those sets and compares with those runs.
2. **Publish the sets**: `gen_d1`, `gen_d2` (and `--parts more` for `d2_v2`), `gen_d3`, `gen_d4` for their organisms. A rebuild with the current code writes the truth of every revision above; it is compared with the published sets and any difference is published as a revision set (`revise_d1`, `revise_d2`, `revise_d4`), never by rewriting a set.
3. **Run and score**: D1 and D3 with F0's runner on `d4c2dda`; D2 and D4 with `run_d` on `d4c2dda` and `71690ab`, scored with every set in publication order.
4. **Report**: `report_d1`, `bridges_d1`, `report_d2`, `report_d3`, `report_d4`.
5. **Classify**: `explain` on the tree, `causes` with D's plugin; `causes_d4_from_single`.
6. **Limitations**: `limits_d --out <tsv>` over the runs of each part and tree.
7. **Results**: `idtrack-repro results` writes `results/d/<run>/<organism>/` once the runs are in `results.toml`.
