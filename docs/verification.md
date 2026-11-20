# Verification of the move

The benchmarks' code was built in another repository and moved into this package: imports made package imports, every location read from the configuration, jobs submitted through the configured runner, decision numbers of the earlier working notes replaced by the rules' names ([rules](rules.md)). To show that the move changed no result, the stored scorings of the reference tree `d4c2dda` (and of `71690ab` where a family measured it) were scored again by this package from the same stored answers and compared with the stored outcomes **line by line** (every column of every outcome row: outcome, reason, answer, seconds, the judge of the line) and the tables written from them (`scores.tsv`, `report.tsv`, the lists of wrong and undetermined answers, errors, timeouts).

**What was compared, and what was not, in short.** Compared line by line, every column: the harness's ten scorings (five organisms, kinds and travel), F0's ten check runs (both trees), A* on `d4c2dda` (its ten runs, the merge of mouse not compared), and family A on pig; in part: A on mouse (31 of 32 partitions) and on human (3 of 45). Not compared line by line: family A on rhesus, cynomolgus and the snapshot-110 graphs, family E (all ten runs), and families B, C and D (imported as they stood, never scored again by this package). For E, A* and part of A an earlier check exists, of tables only (below). Every result of an uncompared unit in `results/` is the stored scoring of the earlier version of the same code.

The comparison stopped before it was complete (the jobs that were still queued were cancelled on 2026-11-19): **30 of its 47 units were compared completely, three more in part, and 14 not at all; no compared line differs.** The units are the harness's scorings per organism (`kinds`: every set but the travel matrix; `travel`: the travel matrix) and each family's scoring per tree and organism.

| unit | rows compared | identical | reason renamed | different | tables |
| --- | ---: | ---: | ---: | ---: | --- |
| harness_d4c2dda_sus_scrofa_kinds | 288,176 | 288,174 | 2 | 0 | equal |
| harness_d4c2dda_sus_scrofa_travel | 192,251 | 192,251 | 0 | 0 | equal |
| harness_d4c2dda_mus_musculus_kinds | 502,615 | 502,472 | 143 | 0 | equal |
| harness_d4c2dda_mus_musculus_travel | 204,557 | 204,557 | 0 | 0 | equal |
| harness_d4c2dda_homo_sapiens_kinds | 461,049 | 460,602 | 447 | 0 | equal |
| harness_d4c2dda_homo_sapiens_travel | 201,449 | 201,449 | 0 | 0 | equal |
| harness_d4c2dda_macaca_mulatta_kinds | 271,886 | 271,802 | 84 | 0 | equal |
| harness_d4c2dda_macaca_mulatta_travel | 193,597 | 193,597 | 0 | 0 | equal |
| harness_d4c2dda_macaca_fascicularis_kinds | 236,336 | 236,322 | 14 | 0 | equal |
| harness_d4c2dda_macaca_fascicularis_travel | 59,401 | 59,401 | 0 | 0 | equal |
| f0_check_d4c2dda_sus_scrofa | 480,425 | 480,423 | 2 | 0 | equal |
| f0_check_d4c2dda_mus_musculus | 707,170 | 707,027 | 143 | 0 | equal |
| f0_check_d4c2dda_homo_sapiens | 662,496 | 662,049 | 447 | 0 | equal |
| f0_check_d4c2dda_macaca_mulatta | 465,481 | 465,397 | 84 | 0 | equal |
| f0_check_d4c2dda_macaca_fascicularis | 295,735 | 295,721 | 14 | 0 | equal |
| f0_check_71690ab_sus_scrofa | 480,425 | 480,423 | 2 | 0 | equal |
| f0_check_71690ab_mus_musculus | 707,170 | 707,027 | 143 | 0 | equal |
| f0_check_71690ab_homo_sapiens | 662,496 | 662,049 | 447 | 0 | equal |
| f0_check_71690ab_macaca_mulatta | 465,481 | 465,397 | 84 | 0 | equal |
| f0_check_71690ab_macaca_fascicularis | 295,735 | 295,721 | 14 | 0 | equal |
| astar_d4c2dda_sus_scrofa | 38,719,226 | 38,719,226 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_mus_musculus (37 of 37 partitions) | 108,689,726 | 108,689,726 | 0 | 0 | not compared (merge not run) |
| astar_d4c2dda_homo_sapiens | 108,317,648 | 108,317,648 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_macaca_mulatta | 30,703,584 | 30,703,584 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_macaca_fascicularis | 3,385,652 | 3,385,652 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_tx_sus_scrofa | 2,959,016 | 2,959,016 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_tx_mus_musculus | 4,346,302 | 4,346,302 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_tx_homo_sapiens | 6,817,212 | 6,817,212 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_tx_macaca_mulatta | 1,448,730 | 1,448,730 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| astar_d4c2dda_tx_macaca_fascicularis | 276,146 | 276,146 | 0 | 0 | equal |
| a_d4c2dda_sus_scrofa | 58,633,038 | 58,633,038 | 0 | 0 | equal (same lines, rows of one qid reordered) |
| a_d4c2dda_mus_musculus (31 of 32 partitions) | 88,275,698 | 88,275,698 | 0 | 0 | not compared (merge not run) |
| a_d4c2dda_homo_sapiens (3 of 45 partitions) | 8,090,885 | 8,090,885 | 0 | 0 | not compared (merge not run) |
| a_d4c2dda_macaca_mulatta | not compared | | | | |
| a_d4c2dda_macaca_fascicularis | not compared | | | | |
| a_d4c2dda_s110_sus_scrofa | not compared | | | | |
| a_d4c2dda_s110_homo_sapiens | not compared | | | | |
| e_d4c2dda_sus_scrofa | not compared | | | | |
| e_d4c2dda_mus_musculus | not compared | | | | |
| e_d4c2dda_homo_sapiens | not compared | | | | |
| e_d4c2dda_macaca_mulatta | not compared | | | | |
| e_d4c2dda_macaca_fascicularis | not compared | | | | |
| e_71690ab_sus_scrofa | not compared | | | | |
| e_71690ab_mus_musculus | not compared | | | | |
| e_71690ab_homo_sapiens | not compared | | | | |
| e_71690ab_macaca_mulatta | not compared | | | | |
| e_71690ab_macaca_fascicularis | not compared | | | | |
| **all compared** | **468,496,794** | **468,494,724** | **2,070** | **0** | |


What "identical" excludes, and nothing else:

- **reason renamed**: 2,070 reasons cite a rule; the stored text names it by the number the earlier working notes gave it, the package by its name (rule `support_ranks_shared_name`: the returned genes carry the name in more name sources than every expected gene left out); the rest of each line is identical.
- **rows of one qid reordered**: in the lists of wrong and undetermined answers, the two rows of a line asked under both views (`best`, `all`) were written in the order of an unstable sort before; the package orders them by qid and variant, so the same lines appear in a fixed order.

Two faults of the moved code were found by the comparison and fixed before the units above were compared: a missing text was written `nan` under pandas 3 (the stored tables were written under pandas 2, which reads it as None; the package writes an empty field under both), and the row order just named.

**Not compared**: family A on rhesus, cynomolgus and the snapshot-110 graphs (pig, human), the merges of A* and A on mouse and of A on human (their partitions compared equal as far as the table shows), and family E (every unit: `d4c2dda` and `71690ab`, five organisms). Their results in `results/` are the stored scorings, made by the earlier version of the same code.

Further checks, without the comparison's jobs: F0's scoring of the harness's sets under the current judge (every published set) equals the harness's own current scoring kind by kind and variant by variant for all five organisms (two scorers, one judge); the travel matrix's tables written by the package from the rescored travel runs are byte-identical to the stored ones (pig and human compared); and the package ran both trees (`d4c2dda`, `71690ab`) through its runner on the cynomolgus graph (two lines per kind; nothing written beside the graph, the graph's content digest recorded and equal to the registry's).

**An earlier check, of tables only.** Before the move, the scorer this package took over was changed (streaming reports and comparisons), and the families' stored scorings were scored again by it and compared table by table (`scores.tsv` and `report.tsv` byte for byte; the lists of wrong, undetermined, error and timeout answers as sorted multisets): family E's 20 scorings (both trees, both judges, five organisms) identical but the reason text of 6 wrong rows (human 4, mouse 2: the rule `contains` names its own reason; the outcome is the same), A*'s 10 identical, A's 4 of 7 identical (pig, rhesus, cynomolgus, pig's snapshot 110; mouse, human and human's snapshot 110 not run). Its tables: `<families>/f0/recheck_2026_11_19/recheck_{e,astar,a}.tsv`. It predates the move's edits (the package's imports, the configured locations, the rules' names in reasons), so it shows the scorer, not the moved code; it is no line-by-line comparison.

To score a stored run again with the package: `idtrack-repro score --family <f> --organism <o> --run <stored run> --scores-to <folder>` (the harness: `idtrack-repro harness-score ... --out <stored answers> --merge --scores-to <folder>`), then compare the two folders' `outcomes`.
