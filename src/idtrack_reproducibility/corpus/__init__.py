"""The regression corpus: a seeded corpus of api calls per organism, the answers of each tree, and their comparison.

It has no truth: it says whether an answer changed between two trees, not whether it is right (docs/corpus.md).
Modules: `make_corpus` (draw a corpus), `run_corpus` (ask it through a tree; the canonical answer format),
`run_harmonize` (its HarmonizeFeatures part), `compare` (the comparator and its classes of difference),
`summarize_results`, `graph_diff` (two graph files of the d4c2dda generation, node by node), `pipeline` (the steps as
subcommands: `python -m idtrack_reproducibility.corpus make|run|harmonize|compare|reference|summarize`).
"""
