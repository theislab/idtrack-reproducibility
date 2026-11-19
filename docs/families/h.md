# H: performance

Family H measures what idtrack costs in time and memory, per tree and organism, with the package's own code on the configured graphs: importing and loading a graph, single calls and batches, every conversion call of the [regression corpus](../corpus.md) with its own arguments, batch scaling, a joint harmonization, a server of the five organisms with forked workers, the reference cycles a conversion leaves to the garbage collector, and builds. It has no truth and judges no answer: where it asks one question two ways (one call per identifier and one batch; one thread and four; two trees harmonizing the same files) it checks that the answers are the same, and the accuracy of the answers is the [harness](../harness.md)'s and the families' question. This is the family's seed: the measurements of the reference tree `d4c2dda` and of `71690ab`; what the [methods](../methods.md#h-performance) ask of H beyond them is listed under [known limitations](#known-limitations).

## What it measures

| measure | module | what |
| --- | --- | --- |
| `load` | `families.h.load` | in a fresh process: `import idtrack` and its `API` (time, resident memory, the heavy modules imported), the load of one graph through the tree's public api (time, resident and peak memory, nodes and edges, files written beside the graph: none), the objects the cyclic collector tracks and one collection; for the d4c2dda generation also the memory of the graph's store by part and the steps of its load one by one |
| `calls` | `families.h.calls` | per-call latency and batch throughput on the corpus's distinct identifiers: one call each with the defaults on cold caches, one batch on warm caches (its answers equal to the single calls'), one batch to the symbol database, the best and all answers with explained paths, the batch on four threads against one; every `convert_identifier` call of the corpus with its own arguments (latency overall and per stratum); batch scaling |
| `harmonize` | `families.h.harmonize` | one joint harmonization of three datasets of 50 cells to the symbol database: the constructor (every feature converted), `unify_multiple_anndatas`, the constructor again on the result cache, with the time of the collector in each, the peak memory, and digests of every output |
| `server` | `families.h.server` | the five graphs loaded in one process, then four forked workers serving 4,000 corpus calls of all five organisms each: start time, the memory of the parent and the workers read at one moment (Pss: shared pages divided among the processes sharing them), the workers' latency |
| `cycles` | `families.h.cycles` | how many conversions leave reference cycles (objects only the cyclic collector frees) and how much memory those hold, per kind of call |
| `build` | `families.h.build` | a build through the tree's public api, warm (the tables of an earlier build copied in: the build's own work) or cold (an empty repository: every table downloaded): time, peak memory, bytes received, and whether the graph built has the configured graph's content digest |

## Methods

**Trees and graphs.** Both trees are measured on the same graph files, those `d4c2dda` built (their content digests: `results.toml`, `[graphs.d4c2dda]`), with the same Python (3.13) and the same versions of numpy, pandas, networkx and anndata. Each tree is called through its own public api: the d4c2dda generation through `API.convert_identifier` (one call per view: `best`, or `all`), the typed api (rule [`typed_api`](../rules.md#typed_api)) through `API.convert`, whose one call ranks every candidate and marks the best, so one call gives both views (`best_and_all_explained` is two batches for `d4c2dda` and one for `71690ab`).

**Inputs.** The identifiers are the distinct strings of the organism's corpus `convert_identifier` calls (every stratum, blank strings left out, sampled with seed 20261117 up to 10,000: all of them, from 7,468 for cynomolgus to 8,051 for human, `calls.tsv`, column `n` of `one_call_cold`); the corpus calls are the corpus's own (from 18,190 to 19,353 per organism, `n` of `corpus_calls`, with their release, database, view and explained paths). The harmonization's datasets are made once per organism from the graph (seed 7: `d115`, every gene id active at the last release without version; `dold`, every gene id of release 90 or the graph's first release when later; `s115`, the symbol of every gene of the last release that has one) and given to both trees (their sha256 are recorded with each repetition). The server's requests: a sample (seed 20261117) of the five organisms' corpus calls.

**Interleaving and repetitions.** All measurements of one organism run in one job on one node: for repetition 1, 2 and 3 the trees in turn, the order reversed every other repetition (`d4c2dda`, `71690ab`; `71690ab`, `d4c2dda`; `d4c2dda`, `71690ab`), each measure in a fresh process of the tree's Python with `PYTHONHASHSEED=0`. So the two trees of one organism are compared on the same machine at nearly the same time; across organisms the machines may differ (every table names the processor of each repetition). The server is measured the same way on the organism `all`; `cycles` and the builds once per tree.

**What each number is.** Times are wall-clock seconds (`time.perf_counter`) around the operation named; latency percentiles are the element at index ⌊n·q⌋ of the sorted call times (p50, p90, p99, p99.9, max), in milliseconds; the time in the collector is summed from `gc.callbacks`. Memory is in GB of 10^9 bytes: resident memory from `/proc/self/statm`, peak memory from `ru_maxrss` (the process's maximum resident set), the server's Rss, Pss and private memory from `/proc/<pid>/smaps_rollup`. Every table has one row per repetition and a row `median` per quantity: the lower median of the three repetitions, a measured value. The results below quote the median rows.

## Results

Every number below is the `median` row of the named table (the lower median of three repetitions, unless the table has one repetition); `d4c2dda` and `71690ab` of one organism were measured on one machine, interleaved.

In short (every number from the tables below):

- **Identifiers that are no node of the graph** (unversioned transcript and translation ids, unknown and malformed strings) are where the trees differ: looking up the versions of such a string, `d4c2dda` walks the graph's names up to the string's place (`itertools.islice` from the start of the list), `71690ab` starts at its place (a binary search in the sorted names). On human the median unversioned transcript id takes 48.361 ms at `d4c2dda` and 0.149 ms at `71690ab` (`results/h/<tree>/homo_sapiens/calls_by_stratum.tsv`). That is most of the difference in the corpus calls: human 46.57 s at `d4c2dda` against 4.64 s at `71690ab` for 19,330 calls, p99 53.0 ms against 1.9 ms (`calls.tsv`, `corpus_calls`). The median call is about the same in both trees (0.122 and 0.115 ms on human).
- **Loading** is faster at `71690ab` for every organism and imports neither pandas nor the build's modules (human 3.6085 s against 3.452 s, the api's import 1.3146 s against 0.5028 s; `load.tsv`); the resident memory of a loaded graph is lower for every organism (human 1.1373 against 1.0872 GB).
- **Threads** do not run conversions in parallel in either tree: the batch dealt over four threads takes longer than the same batch in one thread for every organism and tree (human `71690ab` 1.6324 s against 0.8316 s; `calls.tsv`) except cynomolgus at `71690ab` (1.0128 s against 1.0976 s).
- **Harmonization**: both trees give identical outputs (equal digests) in about the same time on an empty project (human 156.1865 s and 161.4337 s, of which the collector 45.709 s and 45.295 s), but **the constructor that reads a valid result cache is much slower at `71690ab`** for every organism (human 6.0477 s against 121.8452 s, pig 2.5047 s against 47.3983 s; `harmonize.tsv`, `from_cache_s`).
- **A server of the five organisms** starts in 9.026 s and 8.256 s and holds the five graphs in 2.3999 and 2.3496 GB; each forked worker holds 0.808 to 0.819 GB of private memory at `d4c2dda` and 0.879 to 0.89 GB at `71690ab` after 4,000 requests; the workers' p99 falls from 32.6 to 36.2 ms to 4.0 to 6.5 ms (`results/h/<tree>/all/server.tsv`).
- **Reference cycles**: in both trees 1110 of 2000 pig conversions with the defaults leave objects that only the cyclic collector frees (`cycles.tsv`).
- **Builds**: a warm build of pig takes 122.4 s at `d4c2dda` and 124.7 s at `71690ab`, with a peak of 1.1739 and 1.1708 GB, and both build a graph with the configured graph's content digest (`build.tsv`).

### Load

From `results/h/<tree>/<organism>/load.tsv`: the load of the graph through the api (`load_s`), resident and peak memory after it, the import of the api before it (`import_api_s`).

| organism | nodes | edges | graph file MB | tree | load s | resident GB | peak GB | import of the api s | heavy modules after the load |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- |
| sus_scrofa | 1073214 | 1757774 | 85.2 | d4c2dda | 0.9961 | 0.3826 | 0.401 | 1.0706 | h5py, networkx, numpy, pandas, pymysql, requests, tqdm, yaml |
| sus_scrofa | 1073214 | 1757774 | 85.2 | 71690ab | 0.7707 | 0.3292 | 0.3472 | 0.5269 | networkx, numpy |
| macaca_fascicularis | 683319 | 835071 | 47.5 | d4c2dda | 0.5417 | 0.2571 | 0.2685 | 1.0083 | h5py, networkx, numpy, pandas, pymysql, requests, tqdm, yaml |
| macaca_fascicularis | 683319 | 835071 | 47.5 | 71690ab | 0.3965 | 0.2049 | 0.2162 | 0.4738 | networkx, numpy |
| macaca_mulatta | 931551 | 1297068 | 68.6 | d4c2dda | 0.8431 | 0.3393 | 0.355 | 0.996 | h5py, networkx, numpy, pandas, pymysql, requests, tqdm, yaml |
| macaca_mulatta | 931551 | 1297068 | 68.6 | 71690ab | 0.6101 | 0.2881 | 0.3038 | 0.4719 | networkx, numpy |
| mus_musculus | 2476278 | 5367306 | 228.5 | d4c2dda | 2.8405 | 0.7641 | 0.8392 | 1.3576 | h5py, networkx, numpy, pandas, pymysql, requests, tqdm, yaml |
| mus_musculus | 2476278 | 5367306 | 228.5 | 71690ab | 2.3996 | 0.7159 | 0.7912 | 0.5864 | networkx, numpy |
| homo_sapiens | 3682041 | 8630251 | 345.2 | d4c2dda | 3.6085 | 1.1373 | 1.2614 | 1.3146 | h5py, networkx, numpy, pandas, pymysql, requests, tqdm, yaml |
| homo_sapiens | 3682041 | 8630251 | 345.2 | 71690ab | 3.452 | 1.0872 | 1.2106 | 0.5028 | networkx, numpy |

### Per-call latency

From `results/h/<tree>/<organism>/calls.tsv`: `one_call_cold`, one call per distinct identifier of the corpus with the defaults on cold caches; `corpus_calls`, every `convert_identifier` call of the corpus with its own arguments (calls that raise, e.g. for a database that is not in the graph, are counted with their time).

| organism | tree | identifiers | one call, cold: p50 / p99 / max ms | corpus calls | corpus calls: total s, p50 / p99 / max ms |
| --- | --- | ---: | --- | ---: | --- |
| sus_scrofa | d4c2dda | 8003 | 0.21 / 14.6 / 204.7 | 19353 | 17.09, 0.124 / 14.7 / 443.8 |
| sus_scrofa | 71690ab | 8003 | 0.195 / 2.0 / 209.0 | 19353 | 4.46, 0.124 / 1.6 / 25.1 |
| macaca_fascicularis | d4c2dda | 7468 | 0.141 / 9.7 / 150.5 | 18190 | 11.37, 0.091 / 9.7 / 425.6 |
| macaca_fascicularis | 71690ab | 7468 | 0.141 / 0.9 / 228.4 | 18190 | 3.38, 0.094 / 1.0 / 357.3 |
| macaca_mulatta | d4c2dda | 7627 | 0.167 / 12.7 / 195.1 | 18473 | 14.46, 0.104 / 12.8 / 469.6 |
| macaca_mulatta | 71690ab | 7627 | 0.159 / 1.7 / 166.5 | 18473 | 4.03, 0.102 / 1.4 / 349.5 |
| mus_musculus | d4c2dda | 7924 | 0.211 / 36.5 / 455.7 | 19112 | 41.54, 0.129 / 37.4 / 793.5 |
| mus_musculus | 71690ab | 7924 | 0.189 / 4.8 / 434.3 | 19112 | 13.46, 0.126 / 5.3 / 691.0 |
| homo_sapiens | d4c2dda | 8051 | 0.192 / 51.3 / 362.4 | 19330 | 46.57, 0.122 / 53.0 / 766.6 |
| homo_sapiens | 71690ab | 8051 | 0.175 / 1.5 / 299.2 | 19330 | 4.64, 0.115 / 1.9 / 20.7 |

Per stratum of the corpus (`calls_by_stratum.tsv`), the strata of identifiers that are no node of the graph (unversioned transcript and translation ids, nonexistent and malformed strings) and the slowest stratum of each tree:

| organism | tree | `transcript_unversioned` p50 / p99 ms | `translation_unversioned` p50 / p99 ms | `nonexistent` p50 / p99 ms | `malformed` p50 / p99 ms | slowest stratum by p99 (ms) |
| --- | --- | --- | --- | --- | --- | --- |
| sus_scrofa | d4c2dda | 13.869 / 29.8 | 9.732 / 21.2 | 3.898 / 16.6 | 0.195 / 26.6 | `transcript_unversioned` (29.8) |
| sus_scrofa | 71690ab | 0.138 / 1.3 | 0.152 / 1.5 | 0.013 / 0.1 | 0.097 / 1.4 | `gene_merge` (2.9) |
| macaca_fascicularis | d4c2dda | 9.031 / 19.8 | 6.354 / 14.2 | 2.323 / 14.4 | 0.139 / 20.0 | `malformed` (20.0) |
| macaca_fascicularis | 71690ab | 0.117 / 0.7 | 0.128 / 0.7 | 0.012 / 0.1 | 0.084 / 0.4 | `external_ambiguous` (4.0) |
| macaca_mulatta | d4c2dda | 11.597 / 25.7 | 8.133 / 18.3 | 2.995 / 14.4 | 0.144 / 27.0 | `malformed` (27.0) |
| macaca_mulatta | 71690ab | 0.127 / 0.9 | 0.144 / 0.8 | 0.012 / 0.1 | 0.088 / 0.7 | `gene_split` (6.9) |
| mus_musculus | d4c2dda | 30.849 / 69.8 | 17.279 / 37.6 | 7.244 / 45.7 | 0.205 / 72.1 | `gene_retired` (96.4) |
| mus_musculus | 71690ab | 0.15 / 2.4 | 0.17 / 1.3 | 0.013 / 0.1 | 0.098 / 1.0 | `gene_retired` (96.7) |
| homo_sapiens | d4c2dda | 48.361 / 108.2 | 25.701 / 57.8 | 9.473 / 80.4 | 0.214 / 103.3 | `transcript_unversioned` (108.2) |
| homo_sapiens | 71690ab | 0.149 / 2.8 | 0.173 / 2.3 | 0.011 / 0.1 | 0.094 / 1.8 | `gene_merge` (3.0) |

### Batches and threads

From `calls.tsv`, wall seconds for all distinct identifiers of the corpus: `batch_warm` (one batch, warm caches), `batch_symbols` (the same to the symbol database), `best_and_all_explained` (the best and all answers with explained paths), `threads_4` (the batch dealt over four threads) and `serial_again` (the batch once more in one thread, right after).

| organism | tree | batch warm | to symbols | best and all, explained | 4 threads | 1 thread |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| sus_scrofa | d4c2dda | 6.6549 | 7.4823 | 13.0918 | 8.0283 | 6.28 |
| sus_scrofa | 71690ab | 1.299 | 2.3406 | 2.177 | 1.5002 | 1.2416 |
| macaca_fascicularis | d4c2dda | 4.6066 | 5.6964 | 9.0003 | 4.9595 | 4.4953 |
| macaca_fascicularis | 71690ab | 0.9825 | 2.2032 | 1.8645 | 1.0128 | 1.0976 |
| macaca_mulatta | d4c2dda | 5.7795 | 7.0143 | 11.7302 | 6.4075 | 5.9433 |
| macaca_mulatta | 71690ab | 1.3476 | 2.4984 | 2.2165 | 1.7061 | 1.1359 |
| mus_musculus | d4c2dda | 15.488 | 17.2989 | 31.3553 | 18.1882 | 15.3438 |
| mus_musculus | 71690ab | 4.6219 | 6.4872 | 6.1477 | 5.3303 | 4.3494 |
| homo_sapiens | d4c2dda | 16.351 | 21.131 | 33.0805 | 19.408 | 16.3715 |
| homo_sapiens | 71690ab | 0.8736 | 4.8809 | 1.3701 | 1.6324 | 0.8316 |

Batch scaling (`batch_scaling.tsv`: ms per identifier for a batch of the first 1, 10, 100 and 1,000 identifiers, warm caches):

| organism | tree | 1 | 10 | 100 | 1,000 |
| --- | --- | ---: | ---: | ---: | ---: |
| sus_scrofa | d4c2dda | 0.3 | 0.18 | 0.89 | 0.8033 |
| sus_scrofa | 71690ab | 0.2 | 0.19 | 0.147 | 0.1846 |
| macaca_fascicularis | d4c2dda | 7.0 | 0.8 | 0.558 | 0.5964 |
| macaca_fascicularis | 71690ab | 0.1 | 0.09 | 0.092 | 0.1076 |
| macaca_mulatta | d4c2dda | 0.3 | 0.21 | 0.764 | 0.767 |
| macaca_mulatta | 71690ab | 0.2 | 0.2 | 0.118 | 0.2223 |
| mus_musculus | d4c2dda | 0.4 | 0.11 | 1.264 | 1.7061 |
| mus_musculus | 71690ab | 0.3 | 0.11 | 0.146 | 0.3137 |
| homo_sapiens | d4c2dda | 0.6 | 0.15 | 2.952 | 1.9613 |
| homo_sapiens | 71690ab | 0.5 | 0.16 | 0.113 | 0.1051 |

### Harmonization

From `results/h/<tree>/<organism>/harmonize.tsv`: the constructor (every feature of the three datasets converted, no result cache), `unify_multiple_anndatas`, the constructor again on the result cache, each with the seconds of the cyclic collector in it, and the peak memory. `output digest` covers the merged var table and its types, obs, X and the conversion report: equal digests, equal output.

| organism | features in (d115, dold, s115) | features out | tree | constructor s (collector) | unify s | from the cache s (collector) | peak GB | output digest |
| --- | --- | ---: | --- | --- | ---: | --- | ---: | --- |
| sus_scrofa | 35682, 25880, 14363 | 41308 | d4c2dda | 26.0975 (4.824) | 7.9093 | 2.5047 (1.678) | 0.9588 | 584c8878c1acd112 |
| sus_scrofa | 35682, 25880, 14363 | 41308 | 71690ab | 27.5984 (5.132) | 7.9543 | 47.3983 (2.744) | 0.893 | 584c8878c1acd112 |
| macaca_fascicularis | 38550, 28592, 47 | 46380 | d4c2dda | 22.0279 (4.756) | 7.1373 | 2.1876 (1.514) | 0.7617 | a0371579290e89b4 |
| macaca_fascicularis | 38550, 28592, 47 | 46380 | 71690ab | 23.3178 (4.686) | 7.1942 | 37.9801 (1.528) | 0.7132 | a0371579290e89b4 |
| macaca_mulatta | 35432, 32386, 15012 | 47366 | d4c2dda | 26.7222 (4.98) | 9.4665 | 1.1958 (0.357) | 0.9318 | 7fb494ccbe801120 |
| macaca_mulatta | 35432, 32386, 15012 | 47366 | 71690ab | 29.4903 (6.414) | 8.3194 | 46.8078 (3.021) | 0.874 | 7fb494ccbe801120 |
| mus_musculus | 78334, 53110, 77474 | 82852 | d4c2dda | 118.0144 (39.92) | 21.4504 | 7.1251 (4.476) | 1.9881 | f58411aed88a5f18 |
| mus_musculus | 78334, 53110, 77474 | 82852 | 71690ab | 118.6545 (39.989) | 23.592 | 128.6122 (9.952) | 1.794 | f58411aed88a5f18 |
| homo_sapiens | 87693, 64661, 41364 | 89447 | d4c2dda | 156.1865 (45.709) | 21.2909 | 6.0477 (3.892) | 2.3173 | 93050a6f2ca13664 |
| homo_sapiens | 87693, 64661, 41364 | 89447 | 71690ab | 161.4337 (45.295) | 21.1561 | 121.8452 (6.148) | 2.0976 | 93050a6f2ca13664 |

### A server of the five organisms

From `results/h/<tree>/all/server.tsv`: one process loads the five graphs (`start_s`, from the import on), then four workers forked from it serve 4,000 corpus calls each; memory read with every worker alive.

| tree | start s | resident GB after the loads | Pss of parent and workers GB | private GB per worker | worker p50 ms | worker p99 ms | slowest call ms |
| --- | ---: | ---: | ---: | --- | --- | --- | ---: |
| d4c2dda | 9.026 | 2.3999 | 5.713 | 0.808 to 0.819 | 0.265 to 0.359 | 32.6 to 36.2 | 756.4 |
| 71690ab | 8.256 | 2.3496 | 5.946 | 0.879 to 0.89 | 0.243 to 0.363 | 4.0 to 6.5 | 310.7 |

### Reference cycles

From `results/h/<tree>/sus_scrofa/cycles.tsv`: of the corpus's distinct identifiers (the first `n`), how many conversions leave objects that only the cyclic collector frees, with warm caches.

| tree | kind of call | conversions | leaving cycles | unreachable objects | median / max per conversion | MB they held |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| d4c2dda | `default` | 2000 | 1110 | 92079 | 33 / 8021 | 27.91 |
| d4c2dda | `symbol` | 2000 | 1109 | 134196 | 48 / 15029 | 33.1 |
| d4c2dda | `all_explain` | 2000 | 1110 | 92079 | 33 / 8021 | 29.21 |
| 71690ab | `default` | 2000 | 1110 | 89574 | 31 / 8019 | 23.73 |
| 71690ab | `symbol` | 2000 | 1109 | 131695 | 46 / 15027 | 28.36 |
| 71690ab | `all_explain` | 2000 | 1110 | 89574 | 31 / 8019 | 27.14 |

### Builds

From `results/h/<tree>/<organism>/build.tsv`: a warm build reads the stored tables of the configured graph's build (nothing downloaded but the release list); a cold build starts from an empty repository.

| organism | tree | kind | build s | peak GB | MB received | repository MB | graph file MB | the configured graph |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| sus_scrofa | d4c2dda | warm | 122.4 | 1.1739 | 0.1 | 649.4 | 85.2 | the same content digest |
| sus_scrofa | 71690ab | warm | 124.7 | 1.1708 | 0.1 | 649.4 | 85.2 | the same content digest |
| macaca_fascicularis | d4c2dda | warm | 55.1 | 0.6504 | 0.0 | 272.6 | 47.5 | the same content digest |
| macaca_fascicularis | 71690ab | warm | 53.9 | 0.6519 | 0.0 | 272.6 | 47.5 | the same content digest |

## Known limitations

- **Machines.** The jobs ran on shared nodes of one cluster without exclusive use: other work on a node adds noise, which the interleaving spreads over both trees and the median of three repetitions dampens. The two trees of one organism are compared on one machine; numbers of different organisms may come from different processors (named in each table).
- **What is not measured yet.** Of what the [methods](../methods.md#h-performance) ask of H: the latency of every family per stratum read from the families' own runs (each answer carries its seconds; H's tables do not read them yet); memory over 10^6 mixed queries; processes over a shared memory-mapped graph for 1 to 16 workers (rule [`process_workers`](../rules.md#process_workers), decided, not in `71690ab`: the server here forks workers over graphs each process holds in its own memory); the web container under load; harmonization over the sizes of real datasets (the datasets here have 50 cells); cold builds (every table downloaded: the cold builds of cynomolgus failed for both trees because Ensembl's REST service did not answer from the cluster on the day of the measurements); Python versions other than 3.13.
- **Batch scaling** stops at the number of distinct identifiers of the corpus (under 10,000): the batch of 10,000 is not measured.
- **Cold builds** depend on Ensembl's servers and the network on the day: their time is not the package's alone (the warm builds' is).
- **`cycles`** is measured on pig only, once per tree.
- **Threads.** The four-thread batch shows whether conversions run in parallel in one process (they hold the interpreter lock); it is no measure of a server's concurrency.

## Reproducing from scratch

```bash
# the graphs of each tree (`[trees.<tree>.graphs]`) and the corpus (`<corpus>/corpus/`, see docs/corpus.md) first
# load, calls and harmonize: one job per organism, both trees interleaved, three repetitions
python -m idtrack_reproducibility.families.h.run --organisms sus_scrofa,macaca_fascicularis,macaca_mulatta,mus_musculus,homo_sapiens \
    --measures load,calls,harmonize --trees d4c2dda,71690ab --reps 3 --submit --mem 16 --cpus 4
# the server of the five organisms (32 GB)
python -m idtrack_reproducibility.families.h.run --organisms all --measures server --reps 3 --submit --mem 32 --cpus 6
# reference cycles (pig), warm builds (pig, cynomolgus: the tables of the configured graph's folder), a cold build
python -m idtrack_reproducibility.families.h.run --organisms sus_scrofa --measures cycles --reps 1 --submit
python -m idtrack_reproducibility.families.h.run --organisms sus_scrofa --measures build --build warm --reps 1 --submit
python -m idtrack_reproducibility.families.h.run --organisms macaca_fascicularis --measures build --build warm --reps 1 --submit
python -m idtrack_reproducibility.families.h.run --organisms macaca_fascicularis --measures build --build cold --reps 1 --submit
# the result tables (results/h/<tree>/<organism>/*.tsv with their .provenance), then the repository's results
python -m idtrack_reproducibility.families.h.tables
idtrack-repro results
```

A measurement writes `<performance>/runs/<tree>/<organism>/<measure>_rep<k>.json` and is not repeated when the file exists, so a job that was stopped is submitted again as it was. A cold build needs the network to Ensembl (`[jobs.slurm] extra` can place the job on a machine that reaches it).
