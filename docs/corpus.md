# The regression corpus

The regression corpus is a fixed, seeded set of calls of idtrack's api per organism, about 23,000 each, asked through every tree. It has no truth: comparing the answers of two trees call by call says *whether* an answer changed, never whether it is right. Whether an answer is right is the [accuracy harness](harness.md)'s question (rule [`truth_judges`](rules.md#truth_judges)): a change the corpus finds is a change of behaviour to explain, and the harness and the families judge it. The corpus reaches what the harness's query lines do not: every public entry point (single and batch conversions with their classification, source inference, the engine's own search with its paths and scores, the list accessors, the organism lookup, the harmonizer), unusual argument combinations, malformed and non-string input, and the exact form of every answer and every error message.

## The corpus

One file per organism, `<corpus>/corpus/<organism>.jsonl`, one call per line: `{"qid", "entry", "stratum", "args"}`, with a summary beside it (`<organism>.summary.json`: the seed, the releases, the databases, identifiers per stratum, calls per entry). The published corpus was drawn on 2026-11-11 with seed 20261111 through the api of an earlier idtrack (commit `ba37a8b`) on its graphs, and is fixed: every reference answer belongs to it, so it is never drawn again (a new corpus is a new file with references of its own). The calls name no version of idtrack; an argument that later trees removed is left out when the call is asked (`deprioritize_lrg_genes` of the engine's search: rule [`lrg_copies`](rules.md#lrg_copies)).

Releases: twelve per organism over the graph's window, always the first and the last and both sides of every assembly handoff, the rest filled evenly (pig: 56 61 66 67 72 78 89 90 96 102 108 115).

The strata of identifiers (`make_corpus.BUDGET`; every candidate list sorted before a seeded sample):

| stratum | identifiers | budget |
| --- | --- | --- |
| `gene_graph_at_release` | Ensembl gene nodes active at each corpus release | 60 per release |
| `gene_table_at_release` | `stable_id.version` from Ensembl's raw gene table of each release, per assembly that holds it | 60 per release and assembly |
| `gene_unversioned` | gene ids without version | 600 |
| `gene_retired`, `gene_split`, `gene_merge` | genes that were retired; with history edges to two or more other ids; from two or more | 400 each |
| `gene_latest` | genes active at the graph's last release | 300 |
| `gene_internal_node_name` | the graph's internal names (`<id>.Void`, `<id>.Retired`) typed as queries | 40 |
| `transcript_versioned`, `transcript_unversioned`, `translation_versioned`, `translation_unversioned` | transcript and translation ids | 400, 150, 400, 150 |
| `external_per_database` | external identifiers of each database in the graph | 80 per database |
| `synonym_stripped`, `synonym_prefixed` | synonyms without and with the graph's `synonym_id::<database>` prefix | 300, 40 |
| `external_ambiguous` | external identifiers attached to genes of two or more ids | 300 |
| `external_changed_gene` | of those, the ones whose genes hold them at disjoint releases (the name moved) | 300 |
| `nonexistent` | ids with the organism's gene prefix and random digits, random strings, plausible foreign ids, `NaN`, `None`, `null` | up to 200 |
| `malformed` | empty and white-space strings, surrounding white space and new lines, letter case changed, versions that do not exist (`.999`, `_<n>`, `-<n>`, `.1.1`), separators swapped, 12 and 16 separators in a row, non-ASCII and full-width forms | up to 450 |
| `malformed_non_string` | `None`, `0`, `12345`, `1.5`, `["ACTB"]` | 5 |

Every identifier gives a call with the defaults, a call with a random configuration (the target release: the last or a corpus release; the final database: none, `ensembl_gene`, `base_ensembl_gene` or a database of the primary assembly, and in 2 % of the configurations a database not on the primary assembly or one that does not exist; `best` or `all`; explained paths in 10 %; a start release in 25 %), and with probability 0.4 a third with every candidate and explained paths. Beside them: 50 batches of 40 identifiers (`convert_identifier_multiple`, with the classification of the batch's results), 150 identifier lists drawn from one release with noise added, asked for their source in four modes with and without `report_only_winner`, plus a bad mode and an empty list (`infer_identifier_source`), up to 3,000 calls of the engine's own search with random settings (`track_convert`: bridges, the best-only filter, paths, the target release and database, a start release), the five list accessors, and 23 organism names (`resolve_organism`, which needs the network).

| organism | calls | `convert_identifier` | `track_convert` | `infer_identifier_source` | batches | lists | `resolve_organism` | identifiers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sus_scrofa | 23,644 | 19,366 | 2,998 | 1,202 | 50 | 5 | 23 | 8,083 |
| mus_musculus | 23,402 | 19,124 | 2,998 | 1,202 | 50 | 5 | 23 | 7,975 |
| homo_sapiens | 23,623 | 19,343 | 3,000 | 1,202 | 50 | 5 | 23 | 8,083 |
| macaca_mulatta | 22,762 | 18,486 | 2,996 | 1,202 | 50 | 5 | 23 | 7,723 |
| macaca_fascicularis | 22,481 | 18,202 | 2,999 | 1,202 | 50 | 5 | 23 | 7,586 |

(From `<corpus>/corpus/<organism>.summary.json`.)

The HarmonizeFeatures part: three small datasets of 4 cells made from the corpus's identifiers with seed 20261111 (D1 Ensembl gene ids, versioned and not; D2 external identifiers: names, NCBI, RefSeq, UniProt, synonyms, ambiguous ones; D3 a mix of both with overlaps, nonexistent and malformed ids), harmonized in three configurations (`symbol_latest`: the last release, the organism's symbol database; `ensembl_older`: the third corpus release before the last, `ensembl_gene`; `ncbi_latest`: the last release, `NCBI gene`, or `EntrezGene` where the graph has no `NCBI gene`), each in a fresh project folder. Nine answers per configuration: the conversion report, `feature_harmonizer` of each dataset with `merge="flag"` and `merge="sum"`, and `unify_multiple_anndatas` in union and intersect mode: 27 per organism.

## The canonical answer format

A run writes one line per call (`run_corpus.py`): `{"qid", "entry", "ok", "result" | "error", "dt"}`, between a header `{"meta": ...}` (the tree, its commit and api generation, the graph, the corpus, the hash seed, the cap, the load time and memory) and a closing `{"meta_end": ...}` (resident memory during the run, the peak).

- `result` is the answer in canonical form (`run_corpus.canon`): dict keys sorted; sets as `{"__set__": [...]}` sorted; dicts with keys that are not all strings as `{"__dict__": [[key, value], ...]}` sorted by their JSON text; NaN and infinities as the strings `NaN`, `Infinity`, `-Infinity`; numpy scalars and arrays as Python values; tuples as lists; any other object as `{"__repr__": ...}`. **Lists keep their order**: the order of `target_id` and `last_node` is the engine's ranking (rule [`all_order`](rules.md#all_order)), part of the answer.
- `error` is an exception: its type, its message, `message_norm` (the items of every printed set sorted, object addresses removed, so that the message does not depend on the hash seed or the process) and the function that raised.
- `dt` is the seconds of the call, including writing the answer in canonical form; the comparator ignores it. Latency is measured by [family H](families/h.md), not here.

Two api generations answer the same calls. A tree of the d4c2dda generation is asked with every call's arguments as they are. A tree with the typed api (rule [`typed_api`](rules.md#typed_api), from `71690ab` on) is asked through that api, and each result is written in the fields of the d4c2dda answer (`run_corpus.old_conversion`, `old_engine`): a pure re-arrangement of what the typed result holds (its candidates for the asked view, the status as the three flags of the d4c2dda answer, the explained path's steps), so that the comparator compares answers, not formats. Error types are written as the nearest class d4c2dda raised (the typed input errors derive from the built-in errors it raised), messages that name a renamed argument name it as d4c2dda did, and every typed result is checked to survive its own JSON round trip (`to_dict`, `from_dict`); the run's `meta_end` counts both.

## The comparator

`compare.py` compares two answer files call by call and prints `IDENTICAL (<n> calls)`, or `DIFFERENT: <k> of <n> calls` with the differing calls grouped by entry and class of difference, the strata involved and examples. The classes:

| class | the two answers |
| --- | --- |
| `missing_in_a`, `missing_in_b` | the call is in one file only |
| `timeout_a`, `timeout_b`, `timeout_both` | the call hit the run's cap on that side: never counted as an answer difference, listed so that no call is skipped silently |
| `ok_vs_error` | one side answered, the other raised |
| `error_type` | both raised, different exception types |
| `error_message` | the same type, different normalised messages |
| `error_message_order_only` | the messages differ only in the order of the items of a printed set or dict |
| `order_only` | the answers differ only in the order of list items |
| `target_id:subset`, `target_id:superset`, `target_id:changed` | `convert_identifier`: the second side drops targets, adds targets, or has others |
| `value:<field>` | the first top-level field that differs; for `convert_identifier` read in the order `target_id`, `no_corresponding`, `no_conversion`, `no_target`, `final_database`, `last_node`, then alphabetically (the answer first, then how it was reached) |
| `candidates_changed` | `track_convert`: the Ensembl genes reached differ |
| `candidate:<fields>` | `track_convert`: the same genes; the fields that differ (`final_elements_changed` when a gene's targets differ, `final_conversion_details` when only their details do) |
| `value` | any other answer that differs |

## References of the current trees

A reference is a tree's answers to the whole corpus, named by its commit: `<corpus>/reference/<organism>.<tree>.jsonl.gz` and `<organism>.<tree>.harmonize.jsonl.gz`. A reference is never overwritten. The reference of `d4c2dda` is the base every later tree is compared with; a call whose answer differs from it is a change of behaviour, which the harness and the families judge.

| tree | made | determinism (hash seeds 0 and 1) | calls over the cap of 60 s | against the reference of `d4c2dda` |
| --- | --- | --- | --- | --- |
| `d4c2dda` | 2026-11-17, every call of the five corpora, by the steps of `reference` | `IDENTICAL` in all five organisms | none | the base |
| `71690ab` | 2026-11-19, the five corpora with `--no-network` (Ensembl's REST service did not answer from the cluster that day) | `IDENTICAL` in all five organisms | none | only `error_message`, in two messages: the one for a final database the graph does not have at the release (the same error type, other words) and the one for an unknown mode of source inference (the typed api names the modes `release` and `assembly_release`); no answer, no error type and no output of the harmonizer differs (its 27 calls `IDENTICAL` in all five organisms); the 23 `resolve_organism` calls are `missing_in_b` |

The files: `<corpus>/runs/<tree>/<organism>/determinism.compare.txt` (and `.json`), `reference.json`, and for `71690ab` `compare_vs_d4c2dda.txt` and `compare_vs_d4c2dda.harmonize.txt` (with the differing calls, their strata and examples in the `.json`). A run of `d4c2dda` through this package's runner on the configured graphs (`run --tag check --no-network`, 2026-11-19) answers every call as its reference does: `compare_check_vs_d4c2dda.txt` lists only the 23 `resolve_organism` calls it left out.

## Making and comparing

```bash
# a tree's reference answers (one job per organism; the checkout must be clean and at the configured commit)
python -m idtrack_reproducibility.corpus reference --tree 71690ab --organism sus_scrofa --compare-with d4c2dda --submit
# a run of a tree (hash seed 0, cap 60 s per call), then its comparison with a reference
python -m idtrack_reproducibility.corpus run --tree 71690ab --organism sus_scrofa --tag run1 --submit
python -m idtrack_reproducibility.corpus compare --tree 71690ab --organism sus_scrofa --tag run1 --against d4c2dda
# the HarmonizeFeatures part, and its comparison
python -m idtrack_reproducibility.corpus harmonize --tree 71690ab --organism sus_scrofa
python -m idtrack_reproducibility.corpus compare --tree 71690ab --organism sus_scrofa --harmonize --against d4c2dda
# answers and exceptions of a run by entry and stratum
python -m idtrack_reproducibility.corpus summarize --tree 71690ab --organism sus_scrofa --tag run1
# draw a new corpus (a tree of the d4c2dda generation; --tables: also from Ensembl's raw gene tables)
python -m idtrack_reproducibility.corpus make --tree d4c2dda --organism sus_scrofa --out <new file.jsonl> --tables
# two graph files of one organism, node by node (when a corpus differs on two graphs)
python -m idtrack_reproducibility.corpus.graph_diff --tree d4c2dda --a <graph file> --b <graph file> --out <report.json>
```

The tree, its checkout, its Python and its graphs come from the configuration (`[trees.<tree>]`); the steps that ask idtrack run in the tree's Python, each in its own process; `--submit` runs the subcommand as one job of the configured runner (`jobs.submit`). A run's files lie in `<corpus>/runs/<tree>/<organism>/`: `<tag>.jsonl.gz`, `<tag>.harmonize.jsonl.gz`, `compare_<tag>_vs_<other>.{txt,json}`, `logs/`, and the link to the graph (`repo/`).

`reference` does, in one job: the corpus twice in two processes under the hash seeds 0 and 1 (`run1`, `run2`, cap 60 s per call), their comparison (`determinism.compare.{txt,json}`), the calls of `run1` over the cap asked again with a cap of 900 s and merged in (`timeouts.txt`, `slow.jsonl.gz`, `merged.jsonl.gz`), the HarmonizeFeatures part (`harmonize.jsonl.gz`), the copies into `<corpus>/reference/`, and a summary `reference.json` (the commit, the corpus's sha256, the determinism line, the calls over the cap); with `--compare-with`, both files compared with that tree's reference (`compare_vs_<tree>.{txt,json}`, `compare_vs_<tree>.harmonize.{txt,json}`). `--no-network` leaves out the 23 `resolve_organism` calls, which ask Ensembl over the network (`reference.json` records it): a reference made so holds the other calls, and its comparison lists those 23 as `missing_in_b`.

## Determinism

A tree must give the same answers in every process: the hash seed changes the iteration order of sets and dicts, which may not reach an answer. `reference` checks it on every organism (two processes, hash seeds 0 and 1, compared call by call); every run sets its hash seed, and the canonical form removes what legitimately varies between processes (the order of a printed set, object addresses). The family E measures determinism more widely (seeds, processes, threads, batches, order: [methods](methods.md#e-robustness)).

## Reproducing from scratch

1. The graphs of each tree (`[trees.<tree>.graphs]`), built by the tree itself.
2. The corpus: the published files (`<corpus>/corpus/`, listed in the data manifest); a corpus drawn anew with `make` is another corpus, with references of its own.
3. `reference` for each tree and organism (one graph in memory at a time), `--compare-with d4c2dda` for every later tree.

## What it cannot say

- Whether an answer is right: the corpus has no truth. A difference is a change; the harness and the families judge it.
- `resolve_organism` asks Ensembl over the network, and its answer follows Ensembl's latest release and the network of the day; its calls can change without a change of idtrack (`run --no-network` and `reference --no-network` leave them out).
- The published corpus was drawn from the graphs of an earlier tree: identifiers that later graphs lack are still asked (they are answers too: not found), and strata that later graphs hold more of (new databases, new history) are not widened.
- `dt` includes the canonical form (and, for the typed api, the mapping and the round trip): it is no latency (family H measures latency).
