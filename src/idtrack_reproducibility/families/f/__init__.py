"""Family F: harmonization of count matrices (`idtrack.harmonize.HarmonizeFeatures`) against a known truth
(docs/families/f.md).

F1, synthetic datasets with a known truth: per organism, datasets of five source releases and four forms of feature
names drawn from family A's and B's published lines (and the harness's absent strings) with their truths
(`gen_f`), harmonized by the tree in an environment with idtrack's `harmonize` extra (`run_f`), judged per feature by
the common scorer and, for the matrices, the metadata, the order of features and datasets and the result cache, by
F's own checks (`score_f`, the expected outputs of `truth_f`); causes (`causes_f`) and known-limitation rows
(`limits_f`).
"""
