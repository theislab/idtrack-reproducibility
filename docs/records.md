# The record

Every benchmark keeps, for every query line it asks, one row per line, view and **subject** (an idtrack tree, or another tool): the query, the truth with its source and revision, the answer as the subject returned it, the verdict, the cause of a wrong answer and its known limitation, the time, and what made the answer (the tree's commit or the tool's version, the graph's digest, the Python). These rows are the record: the tables of `results/` are derived from them, and any other table or figure can be made from them. The rows live in the data root ([data](data.md)); the repository holds the code that writes and reads them, this dictionary, their manifests (with checksums) and the derived tables. `idtrack-repro check` fails where a registered result has no rows behind it in the data root ([checks](#checks)).

Record version: **1** (`idtrack_reproducibility.records.RECORD_VERSION`; a scored folder's `run.json` says `record_version`). Every stored scoring is of this version: the scorings made before 2026-11-20 (version 0, which held the query, the truth and the answer on the rows that were not right only) were scored again from their stored answers ([scoring again](#scoring-again)); the loader refuses a folder of an earlier version and says how to score it again.

## Loading

```python
from idtrack_reproducibility import records

# every registered run of family A on pig (a pandas DataFrame, one row per line, view and subject)
df = records.load(parts=["a"], organisms=["sus_scrofa"])

# the harness's sets through the common machinery, both trees; the harness's own scorer's rows are part `harness`
f0 = records.load(parts=["f0"], organisms=["homo_sapiens"])

# a scored folder that is not in the registry: a part of a modest run, a tree under development
part = records.load(folders=["<modest run>/sus_scrofa/e/scores_decided"])

# the weakness hunts' reproducers of a tree (one row per call)
hunts = records.load_hunts("<hunts>/reproduce/d4c2dda")

# the regression corpus's calls of a tree, and family H's timed calls (no truth: outcome `unjudged`)
corpus = records.load_corpus("<corpus>/runs/d4c2dda/sus_scrofa/run1.jsonl.gz")
h = records.load_h("<performance>/runs/<tree>/sus_scrofa")
```

`records.load` reads the runs that `results.toml` registers (filtered by part, organism and subject), the scored folders given, or an export (`export=<folder>`, below), and returns one long table with the columns below in this order. A family of hundreds of millions of rows (A, A*) is best loaded one organism at a time.

### The export

```
idtrack-repro records --parts <part,...|all> --organisms <o,...> --subjects <s,...> --out <folder>
```

writes the record of the registered runs (or of `--folders`) as files for any tool that reads Parquet:

| file | what |
| --- | --- |
| `part=<part>/organism=<organism>/subject=<subject>/records.parquet` | every row of that part, organism and subject (a part's several scored folders together: family E's two scorings, A*'s gene and transcript runs; `folder` and `scoring` tell them apart), the columns below in their order and with their types (zstd-compressed; written outcome file by outcome file, so a family of hundreds of millions of rows is exported in the memory of one partition). A folder name is the value with every character but letters, digits, `_`, `.` and `-` written `_` (`subject=g_Profiler`); the rows hold the value itself |
| `schema.json` | the record's version, every column with its type, its Arrow type, its meaning and its values (from `records.SCHEMA`), the partitioning, the parts kept as summaries (A0, E7, E9: what they keep per line, what only as counts, their registered folders), and the loaders of the rows kept elsewhere (the hunts, the regression corpus, family H, another tool's raw responses) |
| `index.tsv` | one row per file: `file`, `part`, `organism`, `subject`, `rows`, `bytes`, `sha256`, `folders` (the scored folders it was made from) |

`--subjects` takes a subject's name (`d4c2dda`, a tool's name) or a tree's name in the registry (`r2`, `d4c2dda_tx`). Another tool's rows (family G, `subject_kind` `tool`) pass through as its scorer wrote them. Read back:

```python
df = records.load(export="<folder>", parts=["a"], organisms=["sus_scrofa"])   # the columns and types of the schema
```

The integer columns come back as pandas' nullable integers (`Int32`), so a missing `n_got` stays missing. Without `--out` the command counts the rows by part, subject, organism and outcome.

### A worked example: the rate per stratum and tree, with an interval

(The figure needs matplotlib, which the package's environment does not hold; the table needs only the package.)

```python
import matplotlib.pyplot as plt
from idtrack_reproducibility import records
from idtrack_reproducibility.families import report

df = records.load(parts=["f0"], organisms=["mus_musculus"])
df = df[df.outcome.isin(["right", "expected_error", "wrong", "error"])]      # the rate's denominator
g = df.assign(ok=df.outcome.isin(["right", "expected_error"])).groupby(["stratum", "subject"])["ok"]
rate = g.agg(["sum", "count"]).reset_index()
rate[["low", "high"]] = [report.wilson(k, n) for k, n in zip(rate["sum"], rate["count"])]
rate["rate"] = rate["sum"] / rate["count"]

top = rate[rate.subject == rate.subject.iloc[0]].nsmallest(12, "rate")       # the twelve hardest strata
fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(top["rate"], top["stratum"], xerr=[top["rate"] - top["low"], top["high"] - top["rate"]], fmt="o")
ax.set_xlabel("right / (right + wrong + error), 95 % Wilson interval")
fig.tight_layout()
fig.savefig("mouse_hardest_strata.pdf")
```

The causes and limitations of the wrong answers are columns of the same rows: `df[df.outcome == "wrong"].groupby(["cause_class", "cause", "limitation"]).size()`.

## Columns

The table is `records.SCHEMA` (the one source of the column list, the export's Parquet types and its `schema.json`; `idtrack-repro records --schema-table` prints it, and a test compares it with this one). Types: `int32`, `float64`, `string`; `string (JSON)` holds a JSON object as a string; a missing value is a null.

| column | type | meaning | values |
| --- | --- | --- | --- |
| `record_version` | int32 | the version of the record the scored folder was written as | `1`: every row holds the query, the truth and the answer |
| `subject_kind` | string | what answered | `tree` (an idtrack checkout), `tool` (another tool: family G) |
| `subject` | string | its name | a tree: its commit's first seven characters; a tool: its name (`run.json` `subject.name`) |
| `subject_version` | string | its version | a tree: the full commit; a tool: its version or release as the tool reports it |
| `graph_digest` | string | the content digest of the graph a tree was asked on (`,` between several) |  |
| `asked_at` | string | for a tool: the date (ISO 8601) its service was asked; empty for a tree (the scoring's `scored_at` says when) |  |
| `endpoint` | string | for a tool: the service's address or the package and its version; empty for a tree |  |
| `python` | string | the Python the subject ran in | e.g. `3.11.12` |
| `part` | string | the benchmark | `harness` (the harness's own scorer), `f0` (the harness's sets through the common machinery), `astar`, `a`, `b`, `c`, `d`, `e`, `hunts`, `corpus`, `h`, a later family's name |
| `scoring` | string | the judge's scoring folder where a part has several | `scores_now`, `scores_decided` (family E, B), or empty |
| `organism` | string | the organism | the five organisms' Ensembl names |
| `family` | string | the family of the line's sets | `harness` for the harness's own sets |
| `set` | string | the published set whose line is the truth (the latest line of the qid: a revision set's line replaces its base's) |  |
| `qid` | string | the line's stable id within the family and organism |  |
| `view` | string | the view scored | `best`, `all` |
| `stratum` | string | the line's stratum (the family's design) |  |
| `kind` | string | the line's kind (the family's design) |  |
| `variant` | string | the view and whether a release was given | `best,from`, `best,default`, `all,from`, `all,default` |
| `unit` | string | the cluster of the bootstrap (gene id, lineage, name, dataset) |  |
| `query` | string (JSON) | the call | `id`, `from_release`, `to_release`, `final_database`, `from_database`, `strategy` (the view), `expect_error` |
| `judge` | string | the judge of the line | `published`, or `revised:<corrections>` |
| `expected` | string | the truth: the expected ids (`,` between) |  |
| `broad` | string | the wider set Ensembl's data also allow (`,` between); empty where none |  |
| `outcome` | string | the verdict | `right`, `wrong`, `undetermined` (the truth cannot decide), `timeout`, `error`, `expected_error` (the error the docs promise), `not_supported` (a tool that cannot ask such a query), `missing` (no answer stored), `unjudged` (a call without truth: the regression corpus, H's timed calls) |
| `reason` | string | why, in words (the rule that decided it) |  |
| `got` | string | the answer's ids as the subject returned them (`,` between) |  |
| `n_got` | int32 | how many |  |
| `status` | string | the answer's own status | `converted`, `tie`, `no_target`, `no_conversion`, `not_found`, `unknown_version`, as the subject reports it |
| `match` | string | how the identifier was matched | as the subject reports it |
| `error_type` | string | the class of an error raised |  |
| `seconds` | float64 | the time of the call that gave this view |  |
| `calls` | int32 | the calls made for the line |  |
| `cause` | string | for a classified wrong (or undetermined) answer: its cause ([F0's causes](families/f0.md#causes)) |  |
| `cause_class` | string | the class of that cause | classes `a` to `e` |
| `limitation` | string | its row of [the known-limitations table](../results/limitations/limitations.tsv) (`L-...`); empty for a decided judge's scoring |  |
| `scored_at` | string | when the scoring was made |  |
| `folder` | string | the scored folder in the data root |  |

The truth's source is the set's manifest (`<set>.manifest.json`: generator, sources with their checksums, `revises`, `judge`, `revision`); seeds are in the sets' and selections' manifests.

## Another tool (family G)

A tool is a subject of its own: its runs are scored by the same scorer into the same rows, with `run.json` naming it (`"subject": {"kind": "tool", "name": ..., "version": ..., "asked_at": ..., "endpoint": ...}`). A query the tool cannot ask (another organism, a release it does not serve, a database it does not map) is `not_supported`, kept apart from `wrong` and from the rate's denominator.

The tool's raw responses are kept as it returned them, so that its answers can be read again from what it returned: one file per request in `<families>/g/raw/<tool>/<version>/<organism>/` and beside them `index.parquet`, one row per request (`records.RAW_COLUMNS`; `records.load_raw(<folder>)` reads it):

| column | type | meaning |
| --- | --- | --- |
| `request_id` | str | the request's id within the folder |
| `tool`, `tool_version` | str | the tool and its version or release as it reports them |
| `organism` | str | |
| `endpoint` | str | the service's address, or the package and its version |
| `method`, `url` | str | the HTTP method and the address asked (or the function called) |
| `params`, `body` | JSON str | the request's parameters and body as sent |
| `qids` | list of str | the query lines the request asked (a batch asks several) |
| `asked_at` | str | when (ISO 8601) |
| `http_status` | int | the response's status (empty for a package) |
| `seconds` | float | the time of the request |
| `bytes`, `sha256` | int, str | the size and checksum of the response file |
| `file` | str | the response file, in the folder |
| `error` | str | a failed request's error |

The answers the scorer reads are made from these files by the tool's adapter, never by hand; the data manifest of `families/g` lists the files with their checksums.

## What each part keeps

| part | per line and view | kept only as a summary |
| --- | --- | --- |
| harness (`harness-score`) | the record (`outcomes/part_0_of_1.parquet`, part `harness`; `outcomes.tsv.gz` derived from it for `harness-compare` and `travel-report`), the answers `results/shard_*.jsonl` and their load records (the checkout that answered: `run.json` `commit`) | – |
| F0, A*, A, B, C, D, E (the common scorer) | the record (`outcomes/*.parquet`), the answers (`answers/*.parquet`: every field the subject returned, with the explained path where asked), the causes (`causes.tsv.gz`), the load records (`answers/load_*.json`: Python, load time, graph) | – |
| A0 (graph fidelity) | every mismatch that is not right (`<check>.tsv.gz`) | the matches, counted per class (`summary.tsv`) |
| E7 (determinism, concurrency) | every difference (`e7_differences.jsonl.gz`) | the agreeing calls, counted (`e7_report.tsv`) |
| E9 (the web layer) | every difference (`differences.jsonl`) | the agreeing requests, counted (`report.tsv`) |
| F's checks (matrices, metadata, order, result cache) | every violation (`f_violations.jsonl.gz`); F's assignments are rows of the common scorer, one per feature, dataset and configuration (`scores_gene/`, `scores_symbol/`) | the items of each check, counted (`f_checks.tsv`) |
| hunts | every reproducer's call with the recorded and the observed answer (`reproduce/<tree>/all_all_<organism>.json`; `records.load_hunts`) | – |
| the modest set | its parts are scored folders of the record (`<out>/<organism>/<part>/`); its hunts as above | `compare.tsv` (the verdict per part) |
| tree comparisons (`compare`) | every moved answer (`--moved`) | the per-stratum table |
| the regression corpus | every call's answer of each tree (`run1.jsonl.gz`, `merged.jsonl.gz`; `records.load_corpus`: one row per call, outcome `unjudged`, the canonical answer in `reason`) | the determinism and comparison reports |
| H (performance) | every call timed by the measure `calls` (A: one call per identifier, caches cold; F: every corpus call with its arguments): `calls_rep<k>.calls.parquet` beside its report (`records.load_h`, outcome `unjudged`); the families' answers carry every call's seconds too | its tables (load, latency percentiles, batches, builds); the stored runs of `d4c2dda` and `71690ab` (2026-11-19) predate the per-call rows and keep the summaries only |

## Scoring again

`idtrack-repro rescore [--parts ...] [--organisms ...] [--trees ...] --workers N` scores every scored folder that `results.toml` registers (the common scorer's and the harness's own) again from its stored answers, nothing asked again, as the record's version of today: with the sets it read and every set published after it that revises them (bounded to the qids of the sets it read; a decided judge's revision set only in a decided scoring), the same selection, and the judge's code of today. Each partition and each scoring is stored with its scoring key (the sets with their checksums, the selection, the record's version, the checksums of the judge's modules): what is current is left as it is, so the command resumes after an interruption and does nothing when nothing changed. The causes, the families' own tables and the results are then made again by their commands, as each family's document says ("Reproducing from scratch"). All stored scorings were scored again on 2026-11-20 (the judge with the revision sets published on 2026-11-19 and 2026-11-20).

## Checks

`idtrack-repro check` (results): every table of `results/` comes from a registered scored folder whose rows are in the data root: its `outcomes/*.parquet` hold as many rows as its `run.json` counts, of the record's version; a registry entry whose part keeps other rows says what it keeps (`rows = "..."`: A0, E9). The data manifests list every file of the record with its size and sha256 (`manifests/`); `check --data` verifies them all.

## What the record keeps as summaries only

- A0 (the graph against Ensembl's tables), E7 (determinism) and E9 (the web layer) keep every row that is not right or agreeing, and the others counted ([what each part keeps](#what-each-part-keeps)).
- H's stored runs of `d4c2dda` and `71690ab` keep latency summaries; every later run of the measure `calls` keeps every call.

The export's `schema.json` states the first three as they are (`kept_as_summaries`), with the registry's folders that hold them.
