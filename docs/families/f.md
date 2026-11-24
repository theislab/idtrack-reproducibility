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

The published sets (manifests `<families>/f/queries/<organism>/f1_v1.manifest.json`, `f1_symbol_v1.manifest.json`, `f1_r1.manifest.json`; 99 MB of datasets and lines in all): `f1_r1` opens pig 2, mouse 4, human 85, rhesus 1 and cynomolgus 4 made-unique lines.

| organism | source releases (assembly) | datasets | lines | `f1_symbol_v1` |
| --- | --- | ---: | ---: | ---: |
| pig | 61 (Sscrofa9), 72 (Sscrofa10.2), 84 (Sscrofa10.2), 96 (Sscrofa11.1), 108 (Sscrofa11.1) | 18 | 36,000 | 655 |
| mouse | 58 (NCBIM37), 73 (GRCm38), 85 (GRCm38), 97 (GRCm38), 109 (GRCm39) | 18 | 36,000 | 694 |
| human | 51 (NCBI36), 65 (GRCh37), 82 (GRCh38), 95 (GRCh38), 108 (GRCh38) | 18 | 36,000 | 1,127 |
| rhesus | 54 (MMUL_1), 67 (MMUL_1), 79 (MMUL_1), 92 (Mmul_8.0.1), 106 (Mmul_10) | 18 | 36,000 | 1,027 |
| cynomolgus | 93 (Macaca_fascicularis_5.0), 97 (Macaca_fascicularis_5.0), 101 (Macaca_fascicularis_5.0), 106 (Macaca_fascicularis_6.0), 112 (Macaca_fascicularis_6.0) | 20 | 40,000 | 729 |

Lines per stratum of `f1_v1` (every form together):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `f1_absent` | 900 | 900 | 900 | 900 | 1,000 |
| `f1_born_late` | 2,781 | 3,186 | 3,164 | 2,900 | 1,895 |
| `f1_duplicate` | 650 | 650 | 650 | 650 | 750 |
| `f1_made_unique` | 1,000 | 1,000 | 1,000 | 1,000 | 1,000 |
| `f1_merged` | 1,820 | 2,584 | 2,472 | 1,360 | 1,053 |
| `f1_name_moved` | 1,812 | 2,147 | 2,175 | 1,862 | 1,394 |
| `f1_name_shared` | 1,143 | 1,227 | 1,663 | 1,304 | 784 |
| `f1_name_unique` | 11,553 | 10,964 | 10,512 | 11,460 | 12,232 |
| `f1_renamed` | 172 | 194 | 431 | 57 | – |
| `f1_resurrected` | 125 | 8 | 2 | – | – |
| `f1_retired` | 1,896 | 2,308 | 2,400 | 1,800 | 1,800 |
| `f1_split` | 1,800 | 1,800 | 731 | 1,800 | 1,800 |
| `f1_unchanged` | 10,348 | 9,032 | 9,900 | 10,907 | 16,292 |

`f1_resurrected` is small or empty where Ensembl brought few ids back (rhesus and cynomolgus none, mouse 8 and human 2 lines at the chosen releases), `f1_renamed` where few ids were renamed (cynomolgus none); a class with fewer lines at a release than its quota takes them all and the rest is filled with unchanged ids or unique names (`drawn` in the manifest).

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

## Results

The tree `main` (idtrack `cabfd9a`, version 0.1.1, the data root's checkout `trees/main`) on r2's graphs (content digests as `[graphs.r2]` of `results.toml`; `load.json` of each run), the harmonizer in Python 3.13.18 with anndata 0.13.4; runs `<families>/f/runs/main/<organism>/`, 2026-11-24. The result files, copied from the runs by `idtrack-repro results`: `results/f/main/<organism>/gene_report.tsv` and `symbol_report.tsv` (the common report per stratum), `gene_causes_summary.tsv`, `symbol_causes_summary.tsv`, `f_assign.tsv`, `f_counts.tsv` and `f_checks.tsv` (the summary rows of the checks). No conversion timed out or raised an error.

**The assignment** (`f_assign.tsv`, every form and stratum together; unit: the source id's lineage or the name):

| organism | configuration | features | right | wrong | undetermined | error | rate | Wilson 95 % | cluster 95 % | units wrong |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| pig | `gene` | 36,000 | 33,675 | 97 | 2,228 | 0 | 0.9971 | 0.9965–0.9976 | 0.9965–0.9978 | 74 |
| pig | `symbol` | 36,000 | 33,070 | 67 | 2,863 | 0 | 0.9980 | 0.9974–0.9984 | 0.9974–0.9985 | 56 |
| mouse | `gene` | 36,000 | 34,227 | 126 | 1,647 | 0 | 0.9963 | 0.9956–0.9969 | 0.9954–0.9972 | 86 |
| mouse | `symbol` | 36,000 | 33,798 | 66 | 2,136 | 0 | 0.9981 | 0.9975–0.9985 | 0.9974–0.9986 | 49 |
| human | `gene` | 36,000 | 33,573 | 115 | 2,312 | 0 | 0.9966 | 0.9959–0.9972 | 0.9959–0.9973 | 96 |
| human | `symbol` | 36,000 | 32,859 | 92 | 3,049 | 0 | 0.9972 | 0.9966–0.9977 | 0.9965–0.9979 | 76 |
| rhesus | `gene` | 36,000 | 33,790 | 83 | 2,127 | 0 | 0.9975 | 0.9970–0.9980 | 0.9967–0.9983 | 45 |
| rhesus | `symbol` | 36,000 | 32,859 | 4 | 3,137 | 0 | 0.9999 | 0.9997–1.0000 | 0.9997–1.0000 | 3 |
| cynomolgus | `gene` | 40,000 | 38,574 | 65 | 1,361 | 0 | 0.9983 | 0.9979–0.9987 | 0.9976–0.9990 | 24 |
| cynomolgus | `symbol` | 40,000 | 37,906 | 11 | 2,083 | 0 | 0.9997 | 0.9995–0.9998 | 0.9994–0.9999 | 6 |

Per stratum, configuration `gene` (rate with the cluster bootstrap's 95 % interval, wrong answers in brackets; `nan`: every line undetermined):

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `f1_absent` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_born_late` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_duplicate` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_made_unique` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 0.9990 (0.9970–1.0000) [1] | 1.0000 (1.0000–1.0000) [0] |
| `f1_merged` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_name_moved` | 0.8727 (0.8317–0.9078) [55] | 0.9769 (0.9639–0.9878) [31] | 0.9702 (0.9574–0.9804) [42] | 0.9823 (0.9664–0.9949) [7] | 0.9578 (0.9241–0.9863) [13] |
| `f1_name_shared` | 0.9786 (0.9662–0.9895) [24] | 0.9285 (0.9042–0.9495) [72] | 0.9617 (0.9497–0.9732) [49] | 0.9430 (0.9202–0.9627) [73] | 0.9331 (0.8983–0.9625) [52] |
| `f1_name_unique` | 1.0000 (1.0000–1.0000) [0] | 0.9993 (0.9986–0.9997) [8] | 0.9993 (0.9987–0.9998) [7] | 0.9998 (0.9995–1.0000) [2] | 1.0000 (1.0000–1.0000) [0] |
| `f1_renamed` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | – |
| `f1_resurrected` | 1.0000 (1.0000–1.0000) [0] | nan (nan–nan) [0] | nan (nan–nan) [0] | – | – |
| `f1_retired` | 0.9864 (0.9789–0.9927) [18] | 0.9917 (0.9867–0.9957) [15] | 0.9895 (0.9822–0.9956) [17] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_split` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_unchanged` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |

Configuration `symbol` (the default final database), the same:

| stratum | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `f1_absent` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_born_late` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_duplicate` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_made_unique` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 0.9990 (0.9970–1.0000) [1] | 1.0000 (1.0000–1.0000) [0] |
| `f1_merged` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_name_moved` | 0.8828 (0.8441–0.9171) [49] | 0.9765 (0.9633–0.9877) [31] | 0.9694 (0.9575–0.9800) [42] | 0.9974 (0.9916–1.0000) [1] | 0.9632 (0.9324–0.9906) [11] |
| `f1_name_shared` | 1.0000 (1.0000–1.0000) [0] | 0.9777 (0.9612–0.9909) [12] | 0.9552 (0.9368–0.9717) [26] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_name_unique` | 1.0000 (1.0000–1.0000) [0] | 0.9993 (0.9986–0.9997) [8] | 0.9993 (0.9987–0.9998) [7] | 0.9998 (0.9996–1.0000) [2] | 1.0000 (1.0000–1.0000) [0] |
| `f1_renamed` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | – |
| `f1_resurrected` | 1.0000 (1.0000–1.0000) [0] | nan (nan–nan) [0] | nan (nan–nan) [0] | – | – |
| `f1_retired` | 0.9864 (0.9791–0.9928) [18] | 0.9917 (0.9869–0.9961) [15] | 0.9895 (0.9817–0.9952) [17] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_split` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |
| `f1_unchanged` | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] | 1.0000 (1.0000–1.0000) [0] |

No id class but `f1_retired` has a wrong answer in any organism (the history travels of A and the versions of the source release, through the harmonizer's batch), as are the absent strings, the duplicates and the made-unique copies (but one rhesus copy whose feature is itself wrong, below). The wrong assignments are names (moved, shared, a few unique) and retired ids, for the causes below; under `symbol` the names of several genes are open where the harmonizer restricts the tie to the genes with an id in the final database ([Truth](#truth)), so `f1_name_shared` has fewer wrong answers there.

**Count conservation** (`f_counts.tsv`, summed over the strata; decided features only, the others apart; counts as the share of the decided features' counts):

| organism | configuration | merge | decided features | kept | moved | dropped | counts kept | counts moved | counts dropped | apart |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| pig | `gene` | `flag` | 29,171 | 29,171 | 0 | 0 | 100.000 % | 0.000 % | 0.000 % | 6,829 |
| pig | `gene` | `sum` | 29,171 | 29,171 | 0 | 0 | 100.000 % | 0.000 % | 0.000 % | 6,829 |
| pig | `symbol` | `flag` | 26,417 | 26,417 | 0 | 0 | 100.000 % | 0.000 % | 0.000 % | 9,583 |
| pig | `symbol` | `sum` | 29,171 | 29,171 | 0 | 0 | 100.000 % | 0.000 % | 0.000 % | 6,829 |
| mouse | `gene` | `flag` | 29,161 | 29,156 | 5 | 0 | 99.982 % | 0.018 % | 0.000 % | 6,839 |
| mouse | `gene` | `sum` | 29,161 | 29,156 | 5 | 0 | 99.982 % | 0.018 % | 0.000 % | 6,839 |
| mouse | `symbol` | `flag` | 26,534 | 26,529 | 5 | 0 | 99.981 % | 0.019 % | 0.000 % | 9,466 |
| mouse | `symbol` | `sum` | 29,161 | 29,156 | 5 | 0 | 99.982 % | 0.018 % | 0.000 % | 6,839 |
| human | `gene` | `flag` | 29,162 | 29,137 | 25 | 0 | 99.933 % | 0.067 % | 0.000 % | 6,838 |
| human | `gene` | `sum` | 29,162 | 29,137 | 25 | 0 | 99.941 % | 0.059 % | 0.000 % | 6,838 |
| human | `symbol` | `flag` | 26,793 | 26,769 | 24 | 0 | 99.928 % | 0.072 % | 0.000 % | 9,207 |
| human | `symbol` | `sum` | 29,161 | 29,136 | 25 | 0 | 99.941 % | 0.059 % | 0.000 % | 6,839 |
| rhesus | `gene` | `flag` | 28,775 | 28,774 | 1 | 0 | 99.999 % | 0.001 % | 0.000 % | 7,225 |
| rhesus | `gene` | `sum` | 28,775 | 28,772 | 3 | 0 | 99.994 % | 0.006 % | 0.000 % | 7,225 |
| rhesus | `symbol` | `flag` | 26,061 | 26,060 | 1 | 0 | 99.999 % | 0.001 % | 0.000 % | 9,939 |
| rhesus | `symbol` | `sum` | 28,775 | 28,772 | 3 | 0 | 99.994 % | 0.006 % | 0.000 % | 7,225 |
| cynomolgus | `gene` | `flag` | 33,850 | 33,840 | 10 | 0 | 99.993 % | 0.007 % | 0.000 % | 6,150 |
| cynomolgus | `gene` | `sum` | 33,850 | 33,840 | 10 | 0 | 99.993 % | 0.007 % | 0.000 % | 6,150 |
| cynomolgus | `symbol` | `flag` | 31,123 | 31,113 | 10 | 0 | 99.993 % | 0.007 % | 0.000 % | 8,877 |
| cynomolgus | `symbol` | `sum` | 33,850 | 33,840 | 10 | 0 | 99.993 % | 0.007 % | 0.000 % | 6,150 |

No count is dropped in any organism, configuration or merge: every input feature's counts are in exactly one output feature. The moved ones are the wrong assignments of decided features (pig none; mouse 5, human 24 or 25, rhesus 1 to 3, cynomolgus 10 features).

**The checks** (`f_checks.tsv`; items checked, and in brackets those that do not hold; `truth_apart_*`: the expected outputs kept apart because an undetermined feature is on their gene):

| checks | pig | mouse | human | rhesus | cynomolgus |
| --- | ---: | ---: | ---: | ---: | ---: |
| `input_*` | 72 | 72 | 72 | 72 | 80 |
| `plan_*` | 556,024 | 555,144 | 558,508 | 557,080 | 621,152 |
| `truth_*` | 109,490 | 109,405 (20 not) | 110,590 (111 not) | 108,050 (6 not) | 128,603 (40 not) |
| `truth_apart_*` | 210 | 366 | 596 | 308 | 352 |
| `unify_*` | 152,662 | 152,784 | 153,733 | 152,854 | 170,872 |
| `unify_truth_*` | 58,342 | 71,677 | 74,055 | 59,875 | 69,436 |
| `shuffle_*` | 18 | 18 | 18 | 18 | 18 |
| `cache_*` | 118 | 118 | 118 | 118 | 130 |

Every mechanical check holds in all five organisms: each harmonized dataset is what the documented rules make of the harmonizer's own assignment (the output features in order, their inputs, status and values in X and in the layer, the input `var` values), every input in exactly one output, `obs` (values and dtypes), `obsm`, `obsp`, `uns` and the layers kept as they are, `varm` subset or left out where features were summed; every merged object (union and intersection, sum and flag) holds the datasets' outputs block by block with the documented renaming, zeros for absent features, the right `intersection` column, every input `obs` column with each cell's value and the documented dtypes, `obsm`, `obsp` and the layers block by block, and `uns` the entries equal in every dataset. Issue #85's loss of `obs` and other metadata does not occur on `main`. The truth's differences (`truth_feature`) are the wrong assignments of decided features, no output assembled otherwise; every output feature that matched the truth in a dataset is in the merged objects as it should be (`unify_truth_*`).

**Order** (E8's relation; `shuffle_*`): the union of the datasets in another order with each dataset's features permuted gives the same features, blocks, merged columns, cells and answers per feature, for both merges, both configurations and all five organisms.

**The result cache** (`cache_*`, from `<run>/<config>/cache.json`):

| organism | configuration | constructed (s) | again on its cache (s) | without `api` (s) | graph loaded without `api` | changed digest rewritten | changed features converted again |
| --- | --- | ---: | ---: | ---: | --- | --- | ---: |
| pig | `gene` | 34.7 | 31.99 | 32.79 | 0 | yes | 18 of 18 |
| pig | `symbol` | 35.55 | 31.95 | 32.68 | 0 | yes | 18 of 18 |
| mouse | `gene` | 47.89 | 36.49 | 36.94 | 0 | yes | 18 of 18 |
| mouse | `symbol` | 40.9 | 34.86 | 36.51 | 0 | yes | 18 of 18 |
| human | `gene` | 41.76 | 26.37 | 27.31 | 0 | yes | 18 of 18 |
| human | `symbol` | 40.69 | 28.19 | 27.26 | 0 | yes | 18 of 18 |
| rhesus | `gene` | 54.3 | 32.74 | 32.71 | 0 | yes | 18 of 18 |
| rhesus | `symbol` | 54.55 | 33.74 | 33.73 | 0 | yes | 18 of 18 |
| cynomolgus | `gene` | 58.08 | 37.52 | 37.84 | 0 | yes | 20 of 20 |
| cynomolgus | `symbol` | 61.56 | 37.81 | 39.43 | 0 | yes | 20 of 20 |

The headers name the graph's content digest and the rules' version (`rules_version` 4) and nothing of the package's version; the constructor again over the cache gives the same conversions without rewriting a file, a stored digest that is not the graph's is computed again and written with the graph's, the cache is read without loading the graph when no `api` is given, and every dataset whose feature names changed (the permuted ones) is converted again. Reading the cache is only 1.1 to 1.7 times as fast as converting again (the documented limitation of the release: the harmonizer reads its own result cache 20 to 45 times as slowly as the version of 13 November, [H](h.md)).

**Causes of the wrong assignments** (`gene_causes_summary.tsv`, `symbol_causes_summary.tsv`; answers):

| configuration | cause | class | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `gene` | `best_keeps_part_of_shared_name` | b | 27 | 61 | 21 | 74 | 48 |
| `gene` | `external_bridge` | b | 66 | 60 | 83 | 3 | 7 |
| `gene` | `support_ranks_shared_name` | b | 4 | – | 4 | 5 | 6 |
| `symbol` | `best_keeps_part_of_shared_name` | b | 1 | 1 | 2 | – | – |
| `symbol` | `external_bridge` | b | 66 | 60 | 83 | 3 | 7 |
| `gene` | `synonym_spelling_merged` | c | – | 2 | – | – | – |
| `gene` | `unexplained` | ? | – | 3 | 5 | – | – |
| `symbol` | `synonym_spelling_merged` | c | – | 2 | – | – | – |
| `symbol` | `unexplained` | ? | – | 3 | 5 | – | – |
| `gene` | `lrg_record_as_gene` | a | – | – | 2 | – | – |
| `symbol` | `lrg_record_as_gene` | a | – | – | 2 | – | – |
| `gene` | `copy_of_another_feature` | b | – | – | – | 1 | – |
| `symbol` | `copy_of_another_feature` | b | – | – | – | 1 | – |
| `gene` | `name_in_other_source` | e | – | – | – | – | 4 |
| `symbol` | `name_in_other_source` | e | – | – | – | – | 4 |


The wrong assignments are idtrack's conversions, which the harmonizer passes on unchanged: no explained single call answered otherwise than the dataset's batch but the copies, which take their feature's answer by design (rule `made_unique_copies`). Most are the bridges where Ensembl's history ends (class b) and the best view of a name of several genes (class b: the harmonizer reads `best`, which keeps the top-ranked carriers; the truth wants every carrier). One cause is idtrack wrong under its own rules (class a), `lrg_record_as_gene`: human `IKBKG` and `ORAI1` read at 65 are answered with the LRG records `LRG_70.1` and `LRG_93.1` as their genes (rule `lrg_copies`: LRG records follow `copies_apart`), and the harmonized datasets get output features named by the records. `API.convert` answers so alone too: the reproducer is `API.convert("IKBKG", from_release=65, to_release=115)` on the human graph (best genes `LRG_70.1`, nothing apart; the same name at 115 is `ENSG00000269335.8`). Eight wrong answers are `unexplained` (read one by one): human `DUX4` and `POM121C` read at 65 (Ensembl's history ends them; the answer is not at the old place, and no bridge test explains it), `FAM138D`, `HLA-DRB7` and `VARSL` at 65 (a gene outside the carriers), mouse `mmu-let-7d` and `mmu-mir-29b-2` at 73 (a miRNA's stem-loop name answered with another locus's gene).

## Known limitations

F's rows of the known-limitations table come from its causes (`limitations_map_f.tsv`), from the judge's reasons of the undetermined assignments and from the checks of the matrices, metadata, order and cache that do not hold (none on `main`; `limits_f.py`). They are stored in `<families>/f/limitations/limitations_main.tsv`; F's builder is registered in `known_limitations.BUILDERS` (`f`), so they become F's rows of `results/limitations/limitations.tsv` when the table is generated. Answers per organism:

| limitation | configuration | what | class | until | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `L-unexplained` | `gene` | wrong answers that no cause test explains | ? | – | – | 3 | 5 | – | – |
| `L-unexplained` | `symbol` | wrong answers that no cause test explains | ? | – | – | 3 | 5 | – | – |
| `L-lrg-as-gene` | `gene` | an LRG record answered as a feature's gene | a | rule `lrg_copies` (LRG records follow `copies_apart`); rule `lrg_answers_gene` (decided) | – | – | 2 | – | – |
| `L-lrg-as-gene` | `symbol` | an LRG record answered as a feature's gene | a | rule `lrg_copies` (LRG records follow `copies_apart`); rule `lrg_answers_gene` (decided) | – | – | 2 | – | – |
| `L-bridges` | `gene` | no successor in Ensembl's history: a bridge through a shared id | b | decided: rule `bridge_same_place` | 66 | 60 | 83 | 3 | 7 |
| `L-bridges` | `symbol` | no successor in Ensembl's history: a bridge through a shared id | b | decided: rule `bridge_same_place` | 66 | 60 | 83 | 3 | 7 |
| `L-harmonize-best-view` | `gene` | the harmonizer reads the best view of a name of several genes | b | none (the harmonizer's design: one gene per feature) | 27 | 61 | 21 | 74 | 48 |
| `L-harmonize-best-view` | `symbol` | the harmonizer reads the best view of a name of several genes | b | none (the harmonizer's design: one gene per feature) | 1 | 1 | 2 | – | – |
| `L-harmonize-copies` | `gene` | a made-unique copy takes its feature's gene | b | none (rule `made_unique_copies`: by design) | – | – | – | 1 | – |
| `L-harmonize-copies` | `symbol` | a made-unique copy takes its feature's gene | b | none (rule `made_unique_copies`: by design) | – | – | – | 1 | – |
| `L-shared-names` | `gene` | a name of several genes: `best` keeps the better-supported | b | none (rule `support_ranks_shared_name`: by design) | 23 | 204 | 390 | 22 | 13 |
| `L-spellings` | `gene` | two genes with one string in two letter cases | c | decided: rule `spelling_per_edge` | – | 2 | – | – | – |
| `L-spellings` | `symbol` | two genes with one string in two letter cases | c | decided: rule `spelling_per_edge` | – | 2 | – | – | – |
| `L-branch` | `gene` | ids that branched off while the old id lived on | e | decided: rule `main_line_first` | 1,436 | 814 | 494 | 1,495 | 1,094 |
| `L-branch` | `symbol` | ids that branched off while the old id lived on | e | decided: rule `main_line_first` | 1,436 | 814 | 494 | 1,496 | 1,094 |
| `L-bridge-at-place` | `gene` | no successor in Ensembl's history; idtrack answers a gene at the same place | e | decided: rule `bridge_same_place` | 128 | 595 | 1,295 | – | – |
| `L-bridge-at-place` | `symbol` | no successor in Ensembl's history; idtrack answers a gene at the same place | e | decided: rule `bridge_same_place` | 128 | 595 | 1,295 | – | – |
| `L-harmonize-tie-final-db` | `symbol` | a feature of several genes under the default final database | e | none (by design) | 654 | 693 | 1,123 | 1,027 | 729 |
| `L-no-chain` | `gene` | no successor and no assembly chain: the answer cannot be judged | e | – | 643 | 30 | 28 | 614 | 256 |
| `L-no-chain` | `symbol` | no successor and no assembly chain: the answer cannot be judged | e | – | 643 | 30 | 28 | 613 | 256 |
| `L-other-name-sources` | `gene` | names from other name sources | e | none (Ensembl's data; undeterminable) | – | – | – | – | 4 |
| `L-other-name-sources` | `symbol` | names from other name sources | e | none (Ensembl's data; undeterminable) | – | – | – | – | 4 |
| `L-unplaced-copy` | `gene` | a copy on a haplotype or patch without a place on the reference | e | – | – | – | 24 | – | – |
| `L-unplaced-copy` | `symbol` | a copy on a haplotype or patch without a place on the reference | e | – | – | – | 24 | – | – |
| `T-made-unique-taken` | `gene` | a made-unique form that is a name of its own (not idtrack's) | e | the truth: the set's revision `copies_taken` | 2 | 4 | 85 | 1 | 4 |
| `T-made-unique-taken` | `symbol` | a made-unique form that is a name of its own (not idtrack's) | e | the truth: the set's revision `copies_taken` | 2 | 4 | 85 | 1 | 4 |

Beside these, `L-harmonize-cache-reuse` (group performance; one row per organism and configuration, from `cache.json`): constructing the harmonizer again over its valid result cache takes 0.60 to 0.92 of the time of converting every feature (the row's `share`; pig 0.90 and 0.92, mouse 0.76 and 0.85, human 0.63 and 0.69, rhesus 0.60 and 0.62, cynomolgus 0.61 and 0.65), the cost idtrack's documentation states (usage, "Known Costs of This Version"); keep the harmonizer, or the merged AnnData written to a file, rather than constructing it again.

What a user of the harmonizer should know, from these rows and the design: a feature of a name several genes carry is converted to its best-ranked genes only (`L-harmonize-best-view`), and under the default final database to those of them with an id there (`L-harmonize-tie-final-db`); a merge with `merge="flag"` keeps no feature for a gene that several inputs name, or that no input names (the documented rule: a dataset holding a gene's versioned and unversioned id keeps both as merge partners, and its partner named like the gene's id becomes `unconverted:<id>` in a merge with other datasets); made-unique copies are summed into, or kept beside, the gene of the feature they copy (`L-harmonize-copies`); a feature answered with an LRG record is named by it (`L-lrg-as-gene`); reading the result cache is barely faster than converting again. What F does not measure: real matrices (F2), scale (F3), `keep_unconverted=False`, `from_database_by_dataset`, final databases other than the Ensembl genes and the default, and features of other forms (external accessions, transcripts).

## Reproducing from scratch

From an empty data root, with the configuration naming the tree and its graphs:

1. **The sources**: the truth store, the accuracy harness's published sets, family A's A1 sets and patterns, family B's B1 sets ([A](a.md#reproducing-from-scratch), [B](b.md#reproducing-from-scratch)).
2. **Publish the sets**: `gen_f --organism <o>` for each organism, then `truth_f --symbol-judge --name f1_symbol_v1` and `truth_f --copies --name f1_r1`. The draw is seeded (`gen_f.SEED`, `seed_of`): the same sources give the same lines and the same files (their sha256 in the manifest), but for the made-unique forms the corrected `fold_index` no longer draws (those `f1_r1` opens in the published set), which another form takes.
3. **Run**: `run_f` on the tree, one job per organism.
4. **Score**: `score_f`; `explain` and `causes` with F's plugin on each scoring.
5. **Limitations**: `limits_f --out <tsv> <run folders>`.
6. **Results**: `idtrack-repro known-limitations --only f` writes F's rows of `results/limitations/limitations.tsv` (F's builder, `limits_f`, over the registered runs), then `idtrack-repro results --only f` writes `results/f/<tree>/<organism>/` from the runs registered in `results.toml`, with their rows of the index and the provenance.
