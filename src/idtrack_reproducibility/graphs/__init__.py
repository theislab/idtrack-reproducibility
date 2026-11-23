"""A tree's graphs: building them from stored tables, judging a built graph, a stored copy for load measurements.

The commands (docs/graphs.md): `idtrack-repro graph-build` (`graphs.build`), `idtrack-repro graph-judge`
(`graphs.judge`), `idtrack-repro graph-stored` (`graphs.stored`). Their places (`graphs.common`): a tree's graph of an
organism lies in `<graphs>/<tree>/<organism>/` with its build record `build.json` and its judgement `judge.json`; the
local repository it was built in (links to the stored tables, the tables the build derived, the graph file) is the
folder `build/` inside it.
"""
