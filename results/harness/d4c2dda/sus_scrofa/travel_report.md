# Travel matrix: sus_scrofa

rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout

| origin era | direction | within assembly | across a handoff |
| --- | --- | --- | --- |
| Sscrofa9 | backward | 100.00 (99.9-100.0); 3917, 3917/0/0/0 |  |
| Sscrofa9 | same | 100.00 (99.6-100.0); 975, 975/0/0/0 |  |
| Sscrofa9 | forward | 100.00 (99.9-100.0); 5833, 5833/0/0/0 | 100.00 (100.0-100.0); 48797, 46070/0/2727/0 |
| Sscrofa10.2 | backward | 100.00 (100.0-100.0); 9482, 9482/0/0/0 | 98.62 (98.4-98.8); 9251, 8657/121/473/0 |
| Sscrofa10.2 | same | 100.00 (99.6-100.0); 1066, 1063/0/3/0 |  |
| Sscrofa10.2 | forward | 100.00 (100.0-100.0); 12473, 12431/0/42/0 | 100.00 (100.0-100.0); 28283, 25173/0/3110/0 |
| Sscrofa11.1 | backward | 99.40 (99.3-99.5); 16669, 16429/100/140/0 | 100.00 (100.0-100.0); 40493, 37012/0/3481/0 |
| Sscrofa11.1 | same | 100.00 (99.7-100.0); 1228, 1228/0/0/0 |  |
| Sscrofa11.1 | forward | 100.00 (100.0-100.0); 13783, 13612/0/171/0 |  |

| distance (releases) | backward | forward | seconds per query: median, 99th percentile |
| --- | --- | --- | --- |
| 0 | 100.00 (99.9-100.0); 3269 | 100.00 (99.9-100.0); 3269 | 0.0001, 0.000 |
| 1-2 | 99.77 (99.6-99.9); 5357 | 100.00 (99.9-100.0); 5955 | 0.0001, 0.000 |
| 3-5 | 99.73 (99.6-99.8); 7588 | 100.00 (100.0-100.0); 8858 | 0.0001, 0.000 |
| 6-10 | 99.61 (99.5-99.7); 11576 | 100.00 (100.0-100.0); 14416 | 0.0001, 0.001 |
| 11-20 | 99.45 (99.3-99.6); 18256 | 100.00 (100.0-100.0); 25212 | 0.0001, 0.001 |
| 21-40 | 99.80 (99.7-99.9); 26978 | 100.00 (100.0-100.0); 37596 | 0.0001, 0.001 |
| >40 | 100.00 (100.0-100.0); 10057 | 100.00 (100.0-100.0); 17132 | 0.0001, 0.001 |

| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |
| --- | ---: | --- | --- |
| all | 192250 | 181882/221/10147/0 | 99.88 (99.9-99.9) |
| id, all | 60000 | 56949/82/2969/0 | 99.86 (99.8-99.9) |
| id, best | 120000 | 113833/139/6028/0 | 99.88 (99.9-99.9) |
| name, all | 3075 | 2625/0/450/0 | 100.00 (99.9-100.0) |
| name, best | 6150 | 5450/0/700/0 | 100.00 (99.9-100.0) |
| name db, best | 3025 | 3025/0/0/0 | 100.00 (99.9-100.0) |
| travel_born_late | 22830 | 20385/122/2323/0 | 99.41 (99.3-99.5) |
| travel_handoff | 45865 | 44033/44/1788/0 | 99.90 (99.9-99.9) |
| travel_merged | 22920 | 22920/0/0/0 | 100.00 (100.0-100.0) |
| travel_name | 9225 | 8075/0/1150/0 | 100.00 (100.0-100.0) |
| travel_renamed | 22950 | 22950/0/0/0 | 100.00 (100.0-100.0) |
| travel_retired_none | 22590 | 18430/55/4105/0 | 99.70 (99.6-99.8) |
| travel_split | 22920 | 22139/0/781/0 | 100.00 (100.0-100.0) |
| travel_unchanged | 22950 | 22950/0/0/0 | 100.00 (100.0-100.0) |

| check (no truth needed) | checked | violations |
| --- | ---: | ---: |
| identity | 3000 | 0 |
| composition | 3000 | 204 |
| return | 3000 | 0 |

| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness outcome of the direct query and the first leg | n |
| --- | --- | --- | --- | --- | ---: |
| composition | forward then backward | handoff between | disjoint | direct right, first leg right | 54 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg right | 31 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg undetermined | 29 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg right | 20 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg undetermined | 14 |
| composition | forward then backward | one assembly | disjoint | direct right, first leg right | 8 |
| composition | backward then forward | handoff between | composed has more | direct right, first leg right | 6 |
| composition | forward then backward | handoff between | composed has fewer | direct right, first leg right | 5 |
| composition | backward then backward | handoff between | disjoint | direct right, first leg right | 5 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg undetermined | 4 |
| composition | forward then forward | handoff between | composed has fewer | direct right, first leg right | 3 |
| composition | backward then forward | handoff between | disjoint | direct undetermined, first leg undetermined | 3 |
| composition | forward then forward | handoff between | disjoint | direct right, first leg right | 2 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg wrong | 2 |
| composition | backward then backward | handoff between | composed has fewer | direct right, first leg right | 2 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg wrong | 2 |
| composition | backward then forward | handoff between | disjoint | direct wrong, first leg undetermined | 2 |
| composition | forward then backward | one assembly | disjoint | direct right, first leg undetermined | 2 |
| composition | forward then backward | handoff between | disjoint | direct undetermined, first leg undetermined | 1 |
| composition | forward then backward | handoff between | disjoint | direct wrong, first leg undetermined | 1 |
| composition | forward then backward | one assembly | composed has more | direct right, first leg right | 1 |
| composition | backward then forward | handoff between | composed has more | direct right, first leg undetermined | 1 |
| composition | backward then forward | handoff between | composed has fewer | direct right, first leg right | 1 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg right | 1 |
| composition | backward then backward | handoff between | disjoint | direct undetermined, first leg undetermined | 1 |
| composition | backward then backward | handoff between | composed has fewer | direct undetermined, first leg right | 1 |
| composition | forward then backward | one assembly | composed has fewer | direct right, first leg right | 1 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg undetermined | 1 |
