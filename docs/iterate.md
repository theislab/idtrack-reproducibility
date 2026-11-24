# An iteration

An **iteration** is one measurement of an idtrack tree under development: its graphs built and judged, its answers to the [modest set](modest.md) (or to the families' panels) asked and compared with the iteration before it and with the reference, the [regression corpus](corpus.md) asked, and one short report. `idtrack-repro iterate <tree>` makes it as one resumable plan of jobs; calling it again continues the plan where it stands. Code: `idtrack_reproducibility.iterate`.

```
idtrack-repro iterate <tree> [--previous <tree>] [--reference <tree>] [--organisms <o,...|all>] [--families <f,...>]
                             [--tables <folder with {organism}>] [--snapshot 115] [--retry <step,...|failed>] [--dry]
idtrack-repro iterate <tree> --report                  # the report again from what the steps left (no job)
idtrack-repro iterate <tree> --complete [--dry]         # the record of a complete run (below)
```

## The order of the iterations

The configuration names the trees in the order they were made (`idtrack-repro.example.toml`):

```
[iterations]
order = ["d4c2dda", "r1", "r1b", "r2", "main"]
reference = "d4c2dda"            # default: the first of the order
```

A tree's **previous iteration** is the tree before it in `order` (`--previous` names another, or one for a tree not in the order); the **reference** is the tree every iteration is also compared with (`--reference`). Each tree is a configured tree (`[trees.<tree>]`: its checkout, its Python, its graph per organism; [README](../README.md#trees)). The organisms are those the tree's graphs are configured for (`--organisms` takes some of them); a graph the tree has not built yet is configured where `graph-build` writes it, `<graphs>/<tree>/<organism>/graph_<organism>_assembly<key>_min<first>_max<snapshot>.npz` ([graphs](graphs.md#places)).

## The steps

| step | what runs (the existing command, in the job) | after | done when |
| --- | --- | --- | --- |
| `graph_build_<o>` | `graph-build --tree <tree> --organisms <o>` (`--tables`, `--snapshot` passed on) | – | the configured graph exists and its build record (`build.json` beside it) names that file and records its content digest |
| `graph_judge_<o>` | `graph-judge --tree <tree> --organisms <o> --against <previous>` | its build | `judge.json` judges the configured graph against the previous tree and is newer than the build record |
| `modest_<o>` | `modest organism --tree <tree> --organism <o> --out <modest>/<tree>`; a family whose panel of the organism is published (`<families>/<family>/panels/<o>/panel_v<k>`) asks its latest panel (`--selections <family>=panel_v<k>`, [Measurement design](methods.md#measurement-design)), the others their modest selection; the hunts' reproducers are the run's last part (`hunts.json`) | its build | `<modest>/<tree>/<o>/modest.json` is of the configured commit, no part was killed and none failed |
| `compare_previous`, `compare_reference` | `modest compare <modest>/<other> <modest>/<tree> --out <modest>/<tree>/compare_vs_<other>.tsv` (each part's `moved_vs_<other>.tsv` beside its scoring) | every `modest_<o>` | the comparison is newer than both trees' run records |
| `corpus_<o>` | `corpus reference --tree <tree> --organism <o> --compare-with <nearest>` (the latest earlier tree of the order with reference answers), and `corpus compare --tree <tree> --organism <o> --tag reference --against <reference>` | its build | the tree's reference answers exist (a reference is never overwritten) and so does the comparison with the reference |
| `report` | `iterate <tree> --report` | every other step, whatever its end | – (written by every call too) |

`<modest>` and `<iterations>` are locations of the data root (`[data] modest`, `iterations`; default `<root>/modest`, `<root>/iterations`). A step whose other tree has no run (`compare_*` without the other tree's modest run, `corpus_<o>` without a corpus of the organism) is skipped and says why. The hunts' reproducers run once per tree, as the modest run's last part: the report reads them from there. Memory and time: the build and judge jobs as `graphs.common.MEMORY_GB` and `TIME` (measured peaks plus a quarter), the modest jobs as `[jobs.modest_memory_gb]`, the corpus jobs as `[jobs.memory_gb]`.

The jobs go through the configured runner (`jobs.submit`), named `<name_prefix>iterate_<tree>_<step>`, their logs in `<iterations>/<tree>/logs/`. Each step runs `afterok` its predecessors (none runs after a failure), the report `afterany`. A judge that finds a fault and a comparison that finds a regression exit 1 with their verdict written: the job counts that as its result and exits 0, so that a failed job is always a failure (the wrapper `iterate.verdict_command`). A step whose predecessor runs again in a call runs again too (its output would be stale).

## The state and resuming

`<iterations>/<tree>/iterate.json` holds the plan: the tree, its commit, the previous iteration and the reference, the organisms, and per step its command, its dependencies, its memory and time, its status with the reason, and every attempt (job id, when submitted, the scheduler's last state, elapsed seconds, peak memory). Every call computes each step's status again, from its outputs first and from its last job (squeue, then sacct):

| status | meaning | what a call does |
| --- | --- | --- |
| `skipped`, `done` | its output is current (before the plan, or by its job) | nothing |
| `active` | its job is pending or running | nothing |
| `todo` | never submitted | submits it |
| `lost` | its job was preempted or its node failed (the scheduler's), or the scheduler does not know it ten minutes after its submission | submits it again under its name |
| `blocked` | its job waits for ever on a predecessor that did not succeed | cancels that job of the plan and submits it again after its predecessor's new job |
| `failed` | its job failed, ran out of memory or time, or completed without its output | nothing until `--retry <step>` (or `--retry failed`) |
| `cancelled` | its job was cancelled | nothing until `--retry` |
| `waits` | a predecessor failed, was cancelled or is held | nothing |
| `held` | submitting it would exceed `[jobs.limits]` | nothing: the next call tries again |

Nothing resubmits by itself: a lost or held step is submitted when the command is called again. `--dry` prints the plan with every step's status and what a call would submit, and writes and submits nothing.

**Limits.** Before each submission the queue of the user is read (`squeue -u $USER`) against `[jobs.limits]`:

```
[jobs.limits]
max_jobs = 30            # jobs at once, running and pending, the new one included
max_memory_gb = 300      # the memory they ask at once
names = "idtrack_"       # count only the jobs whose names start so (default: every job of the user)
```

A step that would exceed either is `held`, and the steps after it wait.

## The report

`<iterations>/<tree>/iterate_report.md` (to read) and `iterate_report.tsv` (to analyse: one row per fact, columns `section`, `organism`, `part`, `against`, `key`, `value`, `cause`, `detail`), written by every call and by the report's job when the steps end:

- **steps**: each step's status and why;
- **graphs**: per organism the content digest, whether it equals the previous tree's (its build record, else the registry's `[graphs.<tree>]`), `verify`, the history (events held of Ensembl's), the faults, and the differences from the previous graph (nodes added and removed, cross-reference links only in the old or the new graph, with other releases, moved between databases, history edges added and removed), the build's time and peak memory;
- **parts**: per organism and part (family E's two scorings apart), against the previous iteration and against the reference, the verdict of `modest compare` (`ok`, `regression`) and the answers that moved from right (or the error the docs promise) to wrong, to an error and to undetermined, and the gains;
- **right to wrong or error**: every such answer (in the Markdown the first hundred per comparison, in the TSV all) with its cause where the tree's part was classified (`causes.tsv.gz`: `idtrack-repro explain` and `causes`), else `unexplained`;
- **hunts**: per organism the findings that hold as recorded and those that changed;
- **regression corpus**: the reference answers' determinism, their comparison with the nearest earlier tree and with the reference;
- **jobs**: every attempt with its scheduler state, elapsed time, peak memory and the memory asked (sacct).

## A complete run (`--complete`)

A tree whose families were run exhaustively (the final tree: [Measurement design](methods.md#measurement-design)) is recorded with `iterate <tree> --complete`. It reads the families' run folders of the tree, `<families>/<family>/runs/<tree>`, its variants `<tree>_<what>` (a family's runs of other graphs or parts: `<tree>_tx`, `<tree>_s110`, `<tree>_d2`) and F0's `check_<tree>`, each organism's scored folders (`run.json`, family E's and B's scorings), and per family and organism:

- the scored folders with their `run.json` (commit, record version, rows, missing and excluded lines, when scored, graph digests, sets);
- the family's published sets of the organism (revision sets left out) that no folder asked;
- every stratum with lines `missing` (no answer stored), and why: `asked in <folder>` (another folder of the tree asked it, as A*'s gene and transcript runs split the sets), `missing on the reference <tree> as well` (no folder of the reference asked it either), or `gap`;
- `complete`: no set left out, no gap, every folder of the tree's commit and of the record's version.

Family A's graph fidelity (A0, `<families>/a/a0/<tree>*/`) is listed as kept as a summary ([the record](records.md#what-the-record-keeps-as-summaries-only)). The result goes into the state file (`complete`) and the report ("The complete run"). `--complete` submits nothing: the complete run is made by each family's own commands (each family's "Reproducing from scratch"). Submitting it from here would need each family's runner to return its jobs (command, memory, time, the jobs it runs after) instead of submitting them, and the family-specific steps (A0, A8 with its graphs, D's and E's own runners, E7, E9) in the same form; then they would join the plan as steps under the same limits.

## Adding a tree

1. Check out idtrack at the tree's commit, make its environment, and add `[trees.<tree>]` with its graphs where `graph-build` writes them ([README](../README.md#trees), [graphs](graphs.md#places)).
2. Add the tree at the end of `[iterations] order`.
3. `idtrack-repro iterate <tree> --dry`, then `idtrack-repro iterate <tree>`; call it again after the jobs end (or are lost) until every step is `done`, and read `iterate_report.md`.
4. A run that a result rests on is registered in `results.toml` as the families' documents say.
