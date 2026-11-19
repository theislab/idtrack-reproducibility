"""Family H: performance. What idtrack costs in time and memory, per tree and organism, measured by the repository's
own code on the configured graphs (docs/families/h.md).

Measures (each a module run in the tree's Python, writing `<performance>/runs/<tree>/<organism>/<measure>_rep<k>.json`):
`load`, `calls` (per-call latency, batches, threads, the corpus's calls, batch scaling), `harmonize`, `server` (the five
organisms in one process with forked workers), `cycles` (reference cycles left by a conversion), `build`. `run` runs
them interleaved and repeated as jobs; `tables` writes the result tables `results/h/<tree>/<organism>/*.tsv` from the
stored measurements and from the regression corpus's runs (the seconds of every call).
"""

# H measures time and memory, not answers: no causes and no rows in the known-limitations table (its limits are in
# docs/families/h.md).
MEASURES_ANSWERS = False
