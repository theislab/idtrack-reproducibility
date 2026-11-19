# Travel matrix: macaca_mulatta

rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout

| origin era | direction | within assembly | across a handoff |
| --- | --- | --- | --- |
| MMUL_1 | backward | 100.00 (100.0-100.0); 26159, 26159/0/0/0 |  |
| MMUL_1 | same | 100.00 (99.8-100.0); 1652, 1649/0/3/0 |  |
| MMUL_1 | forward | 100.00 (100.0-100.0); 29238, 29208/0/30/0 | 100.00 (100.0-100.0); 49977, 46625/0/3352/0 |
| Mmul_8.0.1 | backward | 100.00 (99.9-100.0); 2821, 2821/0/0/0 | 100.00 (100.0-100.0); 19038, 18354/0/684/0 |
| Mmul_8.0.1 | same | 100.00 (99.2-100.0); 499, 499/0/0/0 |  |
| Mmul_8.0.1 | forward | 100.00 (99.9-100.0); 2668, 2668/0/0/0 | 100.00 (100.0-100.0); 9093, 8589/0/504/0 |
| Mmul_10 | backward | 100.00 (99.9-100.0); 6659, 6659/0/0/0 | 100.00 (100.0-100.0); 38396, 36870/0/1526/0 |
| Mmul_10 | same | 100.00 (99.5-100.0); 799, 799/0/0/0 |  |
| Mmul_10 | forward | 100.00 (99.9-100.0); 6597, 6597/0/0/0 |  |

| distance (releases) | backward | forward | seconds per query: median, 99th percentile |
| --- | --- | --- | --- |
| 0 | 100.00 (99.9-100.0); 2950 | 100.00 (99.9-100.0); 2950 | 0.0001, 0.000 |
| 1-2 | 100.00 (99.9-100.0); 4788 | 100.00 (99.9-100.0); 5265 | 0.0001, 0.000 |
| 3-5 | 100.00 (99.9-100.0); 7118 | 100.00 (100.0-100.0); 7876 | 0.0001, 0.000 |
| 6-10 | 100.00 (100.0-100.0); 11730 | 100.00 (100.0-100.0); 12853 | 0.0001, 0.001 |
| 11-20 | 100.00 (100.0-100.0); 21010 | 100.00 (100.0-100.0); 23615 | 0.0001, 0.001 |
| 21-40 | 100.00 (100.0-100.0); 32147 | 100.00 (100.0-100.0); 32037 | 0.0001, 0.002 |
| >40 | 100.00 (100.0-100.0); 16280 | 100.00 (100.0-100.0); 15927 | 0.0001, 0.002 |

| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |
| --- | ---: | --- | --- |
| all | 193596 | 187497/0/6099/0 | 100.00 (100.0-100.0) |
| id, all | 61200 | 59740/0/1460/0 | 100.00 (100.0-100.0) |
| id, best | 121856 | 118192/0/3664/0 | 100.00 (100.0-100.0) |
| name, all | 3075 | 2717/0/358/0 | 100.00 (99.9-100.0) |
| name, best | 6150 | 5533/0/617/0 | 100.00 (99.9-100.0) |
| name db, best | 1315 | 1315/0/0/0 | 100.00 (99.7-100.0) |
| travel_born_late | 25645 | 23815/0/1830/0 | 100.00 (100.0-100.0) |
| travel_handoff | 51395 | 50855/0/540/0 | 100.00 (100.0-100.0) |
| travel_merged | 25725 | 25725/0/0/0 | 100.00 (100.0-100.0) |
| travel_name | 9225 | 8250/0/975/0 | 100.00 (100.0-100.0) |
| travel_renamed | 4810 | 4810/0/0/0 | 100.00 (99.9-100.0) |
| travel_retired_none | 25346 | 22592/0/2754/0 | 100.00 (100.0-100.0) |
| travel_split | 25725 | 25725/0/0/0 | 100.00 (100.0-100.0) |
| travel_unchanged | 25725 | 25725/0/0/0 | 100.00 (100.0-100.0) |

| check (no truth needed) | checked | violations |
| --- | ---: | ---: |
| identity | 2692 | 0 |
| composition | 3000 | 167 |
| return | 3000 | 0 |

| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness outcome of the direct query and the first leg | n |
| --- | --- | --- | --- | --- | ---: |
| composition | forward then backward | handoff between | disjoint | direct right, first leg right | 73 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg right | 40 |
| composition | backward then forward | handoff between | composed has more | direct right, first leg right | 20 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg undetermined | 12 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg right | 7 |
| composition | forward then forward | handoff between | composed has fewer | direct right, first leg right | 3 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg undetermined | 3 |
| composition | backward then backward | handoff between | composed has fewer | direct undetermined, first leg undetermined | 3 |
| composition | forward then backward | handoff between | disjoint | direct undetermined, first leg undetermined | 3 |
| composition | forward then forward | handoff between | disjoint | direct right, first leg right | 2 |
| composition | forward then forward | handoff between | composed has fewer | direct undetermined, first leg right | 1 |
