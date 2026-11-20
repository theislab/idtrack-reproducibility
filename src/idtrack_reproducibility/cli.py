"""`idtrack-repro`: one command per step of the benchmarks (README, "Commands").

    idtrack-repro <command> [arguments of the step]      e.g. idtrack-repro score --family a --organism sus_scrofa ...
    idtrack-repro <command> --help                       the step's own arguments

Every command is a module of the package (`python -m <module>` runs the same); the table below is the map.
"""
import runpy
import sys

# command -> (module, one line)
COMMANDS = {
    # the truth
    "fetch-gtf": ("idtrack_reproducibility.harness.fetch_gtf", "Ensembl's GTF of every release, as gene and transcript tables"),
    "fetch-tables": ("idtrack_reproducibility.harness.fetch_tables",
                     "core history tables, BioMart archives, chain files, assembly exceptions, VGNC, Vega, GRC alignments"),
    "fetch-rest": ("idtrack_reproducibility.harness.fetch_rest", "Ensembl REST archive/id records (a check of the history reading)"),
    "check-truth": ("idtrack_reproducibility.harness.check_truth", "the history reading against REST"),
    # the accuracy harness
    "make-queries": ("idtrack_reproducibility.harness.make_queries", "publish the harness's query sets and revision sets"),
    "harness-score": ("idtrack_reproducibility.harness.score", "ask and score the harness's sets on a tree (sharded jobs)"),
    "harness-compare": ("idtrack_reproducibility.harness.compare_scores", "two harness scorings: rates, regressions, moved answers"),
    "travel-checks": ("idtrack_reproducibility.harness.travel_checks", "the travel matrix's consistency checks"),
    "travel-report": ("idtrack_reproducibility.harness.travel_report", "the travel matrix's tables"),
    "assembly-map": ("idtrack_reproducibility.harness.assembly_map", "the check of places mapped across assemblies"),
    # the families (common machinery)
    "run": ("idtrack_reproducibility.families.runner", "ask a family's sets on a tree and its graph (sharded, resumable jobs)"),
    "score": ("idtrack_reproducibility.families.scorer", "score a family's run (partitions in parallel, then a merge)"),
    "rescore": ("idtrack_reproducibility.rescore", "every registered scoring again under today's judge, from its stored answers"),
    "report": ("idtrack_reproducibility.families.report", "a scored run's report.tsv and tables"),
    "compare": ("idtrack_reproducibility.families.compare", "two scored runs: rates, regressions, every moved answer"),
    "explain": ("idtrack_reproducibility.families.explain", "the non-right answers asked again with their explained paths"),
    "causes": ("idtrack_reproducibility.families.causes", "a cause for every wrong answer"),
    "limitations": ("idtrack_reproducibility.families.limitations", "rows of the known-limitations table"),
    "panel": ("idtrack_reproducibility.families.panel", "a family's panel: the lines a development stage asks"),
    "modest": ("idtrack_reproducibility.modest", "the modest set: build it, run it on a tree (one job per organism), compare"),
    "check-harness": ("idtrack_reproducibility.families.check_harness", "the common scorer against the harness's scorer, line by line"),
    "known-limitations": ("idtrack_reproducibility.known_limitations",
                          "the known-limitations table of every family, from the current runs"),
    # the regression corpus, the hunts, the performance family
    "corpus": ("idtrack_reproducibility.corpus.pipeline", "the regression corpus: make, run, harmonize, compare, references"),
    "hunts": ("idtrack_reproducibility.hunts.reproduce", "the weakness hunts' reproducers on a tree"),
    "performance": ("idtrack_reproducibility.families.h.run", "family H: load, calls, harmonization, server, builds"),
    # results, data, tidiness
    "results": ("idtrack_reproducibility.results", "the results/ tables, their index and provenance, from the scored runs"),
    "manifest": ("idtrack_reproducibility.manifests", "the data manifests of manifests/areas.toml"),
    "check": ("idtrack_reproducibility.check", "the repository's checks (CI runs them; run before every push)"),
}


def usage():
    w = max(len(c) for c in COMMANDS)
    lines = [__doc__.strip(), "", "Commands:"]
    lines += [f"  {c.ljust(w)}  {d}" for c, (_, d) in COMMANDS.items()]
    lines += ["", "  family <family> <module> [...]   a family's own module (idtrack_reproducibility.families.<family>.<module>)",
              "  config                           the configuration as this process reads it"]
    return "\n".join(lines)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(usage())
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "config":
        import json
        from . import config as CFG
        print(json.dumps(CFG.get().describe(), indent=1))
        return 0
    if cmd == "family":
        if len(rest) < 2:
            raise SystemExit("usage: idtrack-repro family <family> <module> [arguments]")
        mod = f"idtrack_reproducibility.families.{rest[0]}.{rest[1]}"
        rest = rest[2:]
    elif cmd in COMMANDS:
        mod = COMMANDS[cmd][0]
    else:
        raise SystemExit(f"unknown command {cmd!r}\n\n{usage()}")
    sys.argv = [f"idtrack-repro {cmd}"] + rest
    try:
        runpy.run_module(mod, run_name="__main__", alter_sys=True)
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
