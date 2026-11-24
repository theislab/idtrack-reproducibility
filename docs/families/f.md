# F: harmonization of count matrices

Family F asks idtrack's harmonizer (`idtrack.harmonize.HarmonizeFeatures`) to bring several count matrices to one feature space at release 115, and judges what comes out against a known truth: which gene each feature became, where its counts went, whether the matrices and the metadata of the inputs come out as the documentation says (GitHub issue #85: `obs` and other metadata lost between input and output), whether the order of datasets and features changes anything, and whether the result cache is keyed by the graph's content digest. F1, synthetic datasets with a known truth, is built and measured at the modest scale; F2 (the atlas project's matrices) and F3 (scale) are designed in [the methods](../methods.md#f-harmonization) and not built. Its code is `idtrack_reproducibility.families.f`; its data are `<families>/f/` (`queries/<organism>/`, `datasets/<organism>/<set>/`, `runs/<tree>/<organism>/`); it runs on the common machinery ([F0](f0.md)) and draws its truth from families [A](a.md) and [B](b.md).

## What it measures

The specification is the documentation of the tree's harmonizer (on idtrack's main `cabfd9a`: the class docstring of `HarmonizeFeatures`, `feature_harmonizer`, `unify_multiple_anndatas` and the pure core's `output_plan`). Per organism:

- **The assignment of every feature** (one row of [the record](../records.md) per feature and dataset, part `f`): the gene the harmonizer gives it (its `feature_table`: one gene for a converted feature, the genes of an ambiguous one, nothing otherwise), judged by the common scorer with the rule of the A or B line it was drawn from, under the best view (the harmonizer reads `best`). Rate = right / (right + wrong + error), undetermined apart, with Wilson's and the cluster bootstrap's intervals over the unit (the source id's lineage, or the name).
- **Count conservation by class**: for the features whose truth names one destination (decided, [Truth](#truth)), whether their counts are in the output feature the documented rules imply (`kept`), in another one (`moved`) or in none (`dropped`), by stratum, in features and in counts; the undetermined ones apart. The expected is everything kept.
- **The matrices**: each harmonized dataset (`feature_harmonizer`, `merge="sum"` and `"flag"`) and each merged object (`unify_multiple_anndatas`, union and intersection), against the documented rules in two ways: from the harmonizer's own assignment (every output feature, its inputs, status and values, with no truth), and from the truth on the decided features.
- **The metadata** (issue #85): `obs` (every column, every cell's value, the dtypes), `obsm`, `obsp`, `uns`, the layers, `varm` and the input `var` columns, kept as the documentation says, per dataset and in the merged objects.
- **Order** (E8's relation): the union of the datasets in another order, each with its features permuted, gives the same features with the same inputs and values in every dataset's block, the same merged columns and the same answer per feature.
- **The result cache** (rule [`cache_by_digest`](../rules.md#cache_by_digest)): its headers name the graph's content digest and the rules' version and nothing of the package's version; constructing again over it gives the same conversions, faster; a stored digest that is not the graph's is computed again and written with the graph's; without `api` the graph is not loaded; datasets whose feature names changed are converted again.

Two configurations of the harmonizer are run on the same datasets: `gene` (`final_database="ensembl_gene"`: the features become Ensembl genes, which the truth names) and `symbol` (the organism's default final database, what users get: HGNC symbols for human, MGI for mouse, VGNC for pig and the macaques).

## The sets

One set per organism, `f1_v1`, published once (F0's `lines.publish`; `gen_f.py`), with its datasets `<families>/f/datasets/<organism>/f1_v1/<dataset>.h5ad`; `f1_symbol_v1`, the lines the `symbol` configuration judges otherwise (a revision set of `f1_v1` of judge `symbol`; `truth_f.py`); and `f1_r1`, the revision `copies_taken` ([Truth](#truth)). A line is a feature of a dataset: qid `f1|<dataset>|<feature>`, `kind` `f1_<form>`, `stratum` `f1_<class>`, the source line's truth, `from_release` the dataset's release, `to_release` 115, `strategy` `best`, and the fields `dataset`, `form`, `source` (`<family>:<set>:<qid>` of the line it was drawn from), `history_class`, `group` (an n-to-1 group made on purpose: merge partners, copies, duplicates), `copy_of`, `core`.

**Datasets.** Per organism five source releases spread over the window (`gen_f.choose_releases`: one per assembly era first, every further release to the era with the most releases per slot, evenly spaced within it; 115 is the target, no source) and four forms of feature names per release: `ids` (Ensembl gene ids without version), `versioned` (the id with the version it had at the release), `symbols` (the GTF's display names of the release) and `mix` (ids, versioned ids and names, with duplicates). A release whose GTF writes no versions (pig and human before 77, mouse before 76) has no `versioned` dataset, as family A has no versioned line there. Each dataset holds 2,000 features and 200 cells; its name is `r<release>_<form>`.

**Features** (shares per form in `gen_f.QUOTAS`; what was drawn in the manifest, `drawn` and `datasets.<name>.strata`):

| stratum | drawn from | what it is |
| --- | --- | --- |
| `f1_unchanged`, `f1_born_late`, `f1_renamed`, `f1_merged`, `f1_split`, `f1_retired` | A's A1 lines from the release to 115 (`a1_*` sets, the latest line of each qid), unversioned or versioned | an id by its history class (`<families>/a/patterns/<organism>.tsv.gz`; `retired` is `retired_none`); a merged id comes with its merge partner (another id of the release whose line expects the same one gene) |
| `f1_resurrected` | A1 lines whose broad walk reaches ids only through a return (rule [`resurrection`](../rules.md#resurrection)) | an id brought back |
| `f1_name_unique`, `f1_name_shared`, `f1_name_moved` | B's B1 lines read at the release to 115 (`b1_*` sets): `b1_to115_unique`, `b1_to115_shared`, `b2_moved_from` | a display name of one gene, of several genes, a name that moved to another gene |
| `f1_made_unique` | a unique name X of the dataset | `X.1` or `X-1` beside X (rule [`made_unique_copies`](../rules.md#made_unique_copies)), only where no string of the truth store folds to the form (every GTF's ids and names, every BioMart table and VGNC gene set, every release's core cross-references and synonyms); its truth is X's line |
| `f1_absent` | the accuracy harness's `nonexistent` lines | an id with the organism's prefix that never existed (the id forms), random letters and RefSeq-like accessions (the names); nothing expected |
| `f1_duplicate` | an A1 and a B1 line of the release expecting the same one gene; the unversioned and versioned A1 lines of one unchanged id | one gene asked twice in one dataset (`mix` only) |

A name holding `;` is left out: the harmonizer joins a merged feature's input names with `;` in `Query ID`, so such a name could not be told apart there. Literally repeated names are no stratum: the harmonizer refuses them with `InvalidArgumentError`, as documented (family E's argument checks). Two hundred **core** genes (unchanged at every chosen release, with a unique display name at each) are in every dataset in its form, so that the intersection over the datasets is not empty.

**Counts and metadata** (`gen_f.counts`, `gen_f.metadata`; seeds `gen_f.seed_of(organism, dataset, ...)`, in the manifest): `X` drawn per feature and cell from a negative binomial with mean m_f x s_c (m_f log-normal with mu 0 and sigma 1.5, s_c the cell's size factor, log-normal with sigma 0.3), dispersion 2 (float32, sparse); the layer `unspliced` drawn alike with a third of the mean (int32); `obs` with a categorical, text, float, integer and boolean column, and one column of some forms only (`stage` in `symbols` and `mix`, `plate` in `ids`); `obsm` `X_pca` and `X_umap`; `obsp` `connectivities`; `varm` `PCs`; `var` `feature_type`, `highly_variable`, `mean_input`; `uns` with entries equal in every dataset (`organism`, `f1_set`), entries that differ (`dataset`, `source_release`) and a nested dict (`params`). The files are written with h5py in AnnData's on-disk format (`h5ad.py`), so that the generator needs no anndata; the run reads them with anndata and checks that it read what was written (`input_*`).

SETS_TABLE

## Truth

Nothing of idtrack. A feature's expected gene(s) at 115 are the `expect` of the A or B line it was drawn from, copied with its whole truth (`broad`, `rule`, `cmp`, `locus`, `apart`, ...), the latest line of the qid at publication; a revision of the source line published later becomes a revision set of F (`truth_f --check` lists the lines it concerns). A made-unique copy's truth is the line of the name it copies; an absent string expects nothing.

**Decided features** (`truth_f.decided`): a feature's truth names one destination when it expects one gene (and its broad walk allows no other, and it is not open) or nothing (and no place rule could accept a gene: no `locus`, no `broad`). Every other feature (several genes, a wider set, a place rule, open) is undetermined for the matrices and kept apart, never forced; its assignment is still judged by its line's rule.

**The expected outputs** (`truth_f.plan`, `merged_names`, `uns_kept`: the documented rules written again, not imported, so that the package is checked against its own words): a gene reached by one input feature is one output feature named by the gene (its versioned id at 115); several input features on one gene: with `merge="sum"` one feature whose values are their sum, with `merge="flag"` the input that names the gene itself (its stable id, or under the default final database one of the gene's names there) is the gene's feature and the others keep their input names as merge partners, and where none or several name it all keep their names and the gene has no feature; every other input keeps its name (`keep_unconverted`); an unconverted name equal to a converted feature's name, in the dataset or in any dataset when merging, becomes `unconverted:<name>`. The union holds every feature of any dataset, the intersection those of every dataset; `uns` keeps the entries equal in every dataset. Whether an input names its gene is read for Ensembl ids (the id's stable id is the gene's) and, in the configuration `gene`, for every input (a gene's only name there is its id); for a name under the default final database the truth does not read the gene's names there, and a `flag` group holding a name and other inputs is undetermined.

**The revision `copies_taken`** (`f1_r1`): the generator of `f1_v1` mapped each folded string of the truth store to only one of the made-unique forms that fold to it, so a form folding like another (`X-1` and `X.1`) was not seen where only its own spelling is a string of the store, and a few forms that are names of other genes were drawn as copies (mouse `Ighv5-9-1`, a GTF name of its own beside `Ighv5-9`). The generator is corrected (`gen_f.fold_index`); the revision reads the store again for every `f1_made_unique` line and opens each whose form is a string of the store (such a form is read as itself by rule `made_unique_copies`, and the copy's truth does not hold). Its lines per organism are in its manifest.

**The `symbol` judge** (`f1_symbol_v1`): under the default final database the harmonizer judges a tie among the genes that have an id there (its documentation); which of a name's genes have one is not read by the truth, so a feature expecting exactly several genes (rule `set_exact`, or the support of rule [`support_ranks_shared_name`](../rules.md#support_ranks_shared_name)) is open within them in that configuration: an answer within them is undetermined, one outside them wrong.

## Rules

The rules the expectations rest on are A's and B's ([A's rules](a.md#rules), [B's rules](b.md#rules)) under the best view, and the harmonizer's: [`made_unique_copies`](../rules.md#made_unique_copies), [`batch_agrees_with_single`](../rules.md#batch_agrees_with_single), [`start_release`](../rules.md#start_release) and [`from_release_fallback`](../rules.md#from_release_fallback) (a feature the graph does not hold at its dataset's release is read where it holds it), [`missing_values`](../rules.md#missing_values), [`harmonize_from_database`](../rules.md#harmonize_from_database) and [`cache_by_digest`](../rules.md#cache_by_digest).

## Scoring

`score_f.py` turns the run's answers into the runner's answer files and scores them with the common scorer (`scorer.py --family f`; `scores_gene/` with `f1_v1` and its revision sets, `scores_symbol/` with `f1_symbol_v1` too): `outcomes/*.parquet` (the record: `status` is the harmonizer's status of the feature, `converted`, `converted_no_final_id`, `ambiguous`, `no_conversion`, `not_in_graph`, `missing_value` or `error`; `match` the conversion's; `seconds` the dataset's conversion time divided by its features, as the harmonizer converts a dataset in one batch) and `report.tsv` per stratum. Then F's own checks, written as summary rows, stated as such (`f_checks.tsv`: per check, configuration, merge and scope the items checked, those that hold, those that do not and an example; every violation in `f_violations.jsonl.gz`):

| checks | what holds when they hold |
| --- | --- |
| `input_file`, `input_features`, `input_x`, `input_layer` | each dataset as anndata read it in the run is the file the manifest names, with its features in order and its X and layer as the generator wrote them |
| `plan_names_in_order`, `plan_feature`, `plan_values_x`, `plan_values_layer`, `plan_var_inputs`, `plan_every_input_once` | each harmonized dataset is what the documented rules make of the harmonizer's own assignment: the output features in the order of their first input, each with its inputs, status and gene, its values in X and in the layer the sums of its inputs', its input `var` values kept where its inputs agree; every input in exactly one output |
| `plan_obs`, `plan_obs_dtypes`, `plan_obsm`, `plan_obsp`, `plan_uns`, `plan_layers`, `plan_varm_subset`, `plan_varm_left_out_when_summed` | `obs`, `obsm`, `obsp`, `uns` and the layers kept as they are (issue #85); `varm` subset to the output features, or left out where features were summed |
| `truth_feature` (and `truth_feature_apart`) | each output feature the truth implies on the decided features is there, with the same inputs, status and values; judged where every input on its gene (expected or as assigned) and every input of the output holding it is decided, else apart |
| `unify_features`, `unify_block`, `unify_every_block`, `unify_absent_blocks_zero`, `unify_intersection_column` | the merged object's features are the union or the intersection of the datasets' outputs after the documented renaming; each dataset's block of a feature holds its output's inputs and values; a feature a dataset lacks is zeros there; `intersection` is 1 for the features of every dataset |
| `unify_obs_names`, `unify_obs_values`, `unify_obs_dtypes`, `unify_obsm_blocks`, `unify_obsp_blocks`, `unify_layers`, `unify_uns` | every input `obs` column with each cell's value (missing where a dataset lacks the column), the documented dtypes (text categorical, a partial integer or boolean column nullable), the dataset column and the cell names `<dataset>_<cell>`; `obsm`, `obsp` (block-diagonal) and the layers block by block; `uns` the entries equal in every dataset (issue #85) |
| `unify_truth_union`, `unify_truth_intersect` | every output feature that matched the truth in a dataset is in the union, and in the intersection where it matched in every dataset |
| `shuffle_features`, `shuffle_merged_columns`, `shuffle_blocks`, `shuffle_obs`, `shuffle_answers` | the order relation (E8's) |
| `cache_header_digest`, `cache_header_rules_version`, `cache_header_no_package_version`, `cache_again_equal`, `cache_again_faster`, `cache_changed_digest_recomputed`, `cache_without_api_no_graph_load`, `cache_changed_features_recomputed` | rule `cache_by_digest` |

`f_counts.tsv` is the count conservation by class, `f_assign.tsv` the assignment rate per configuration, form and stratum with both intervals. The wrong assignments are explained on the tree (`explain.py`) and classified by F0's classifier with F's plugin (`causes_f.py`: F's kinds, and the causes only a dataset's batch has: `copy_of_another_feature`, `batch_differs_from_single`, and `best_keeps_part_of_shared_name` for the best view of a name of several genes).

## Running it

Locations, the tree and its graphs come from the configuration ([running F0](f0.md#running-it)). The harmonizer runs in an environment with idtrack's `harmonize` extra (anndata 0.13 or later), with the tree's checkout first on the path; the generator and the scorer run in the package's environment (pyarrow, h5py; no anndata).

```
# the set (once per organism; reads A's and B's published lines, the harness's, the truth store) and the symbol judge
python -m idtrack_reproducibility.families.f.gen_f --organism <o> [--set f1_v1] [--dry]
python -m idtrack_reproducibility.families.f.truth_f --organism <o> --symbol-judge --name f1_symbol_v1
python -m idtrack_reproducibility.families.f.truth_f --organism <o> --copies --name f1_r1  # the revision copies_taken
python -m idtrack_reproducibility.families.f.truth_f --organism <o> --check     # source lines revised since
# the run (one job per organism; the environment with the harmonize extra)
<python with anndata> -m idtrack_reproducibility.families.f.run_f --organism <o> --clone <source of the tree> \
    --graph <graph of the tree> --out <families>/f/runs/<tree>/<o>
# the judge, the causes, the limitations (the package's environment; explain in the tree's Python)
python -m idtrack_reproducibility.families.f.score_f --organism <o> --run <families>/f/runs/<tree>/<o>
<python of the tree> -m idtrack_reproducibility.families.explain --family f --organism <o> --clone <source> \
    --graph <graph> --scored <families>/f/runs/<tree>/<o>/scores_<config>
python -m idtrack_reproducibility.families.causes --organism <o> --scored <families>/f/runs/<tree>/<o>/scores_<config> \
    --plugin src/idtrack_reproducibility/families/f/causes_f.py
python -m idtrack_reproducibility.families.f.limits_f --out <tsv> <families>/f/runs/<tree>/<o> ...
```

A run folder holds, per configuration, `answers.jsonl.gz` (the harmonizer's answer per feature), `outputs/<merge>/<dataset>.tsv.gz` and `.meta.json` (every output feature with the digests of its values, and the object's metadata), `unify/<features>_<merge>.*`, `shuffle/union_<merge>.*`, `cache.json`; and `load.json`, `input_check.json`, `run_f.json`, the permuted datasets `shuffle_datasets/`. Values are compared by digests of each feature's column (its values as float64 in cell order; counts are exact in float32), so that the run need not store its matrices.

RESULTS

## Known limitations

KNOWN_LIMITATIONS

## Reproducing from scratch

From an empty data root, with the configuration naming the tree and its graphs:

1. **The sources**: the truth store, the accuracy harness's published sets, family A's A1 sets and patterns, family B's B1 sets ([A](a.md#reproducing-from-scratch), [B](b.md#reproducing-from-scratch)).
2. **Publish the sets**: `gen_f --organism <o>` for each organism, then `truth_f --symbol-judge --name f1_symbol_v1` and `truth_f --copies --name f1_r1`. The draw is seeded (`gen_f.SEED`, `seed_of`): the same sources give the same lines and the same files (their sha256 in the manifest), but for the made-unique forms the corrected `fold_index` no longer draws (those `f1_r1` opens in the published set), which another form takes.
3. **Run**: `run_f` on the tree, one job per organism.
4. **Score**: `score_f`; `explain` and `causes` with F's plugin on each scoring.
5. **Limitations**: `limits_f --out <tsv> <run folders>`.
6. **Results**: `idtrack-repro results` writes `results/f/<tree>/<organism>/` from the runs registered in `results.toml`.
