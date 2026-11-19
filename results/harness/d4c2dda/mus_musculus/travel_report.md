# Travel matrix: mus_musculus

rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout

| origin era | direction | within assembly | across a handoff |
| --- | --- | --- | --- |
| NCBIM37 | backward | 98.56 (98.3-98.8); 10150, 9771/143/236/0 |  |
| NCBIM37 | same | 100.00 (99.7-100.0); 1153, 1153/0/0/0 |  |
| NCBIM37 | forward | 100.00 (100.0-100.0); 11757, 11394/0/363/0 | 99.70 (99.7-99.7); 55901, 51981/154/3766/0 |
| GRCm38 | backward | 98.32 (98.1-98.5); 21221, 20577/352/292/0 | 97.08 (96.9-97.3); 22900, 21447/645/808/0 |
| GRCm38 | same | 100.00 (99.7-100.0); 1378, 1372/0/6/0 |  |
| GRCm38 | forward | 99.17 (99.0-99.3); 23972, 22817/192/963/0 | 97.88 (97.7-98.1); 18158, 17093/370/695/0 |
| GRCm39 | backward | 100.00 (99.9-100.0); 5281, 5250/0/31/0 | 99.43 (99.3-99.5); 30579, 29913/170/496/0 |
| GRCm39 | same | 100.00 (99.4-100.0); 591, 591/0/0/0 |  |
| GRCm39 | forward | 100.00 (99.7-100.0); 1515, 1503/0/12/0 |  |

| distance (releases) | backward | forward | seconds per query: median, 99th percentile |
| --- | --- | --- | --- |
| 0 | 100.00 (99.9-100.0); 3122 | 100.00 (99.9-100.0); 3122 | 0.0001, 0.000 |
| 1-2 | 99.75 (99.6-99.9); 5185 | 99.69 (99.5-99.8); 5560 | 0.0001, 0.001 |
| 3-5 | 99.43 (99.2-99.6); 7579 | 99.58 (99.4-99.7); 8083 | 0.0001, 0.001 |
| 6-10 | 99.11 (98.9-99.3); 11720 | 99.47 (99.3-99.6); 13176 | 0.0001, 0.001 |
| 11-20 | 98.63 (98.5-98.8); 19855 | 99.22 (99.1-99.3); 24087 | 0.0001, 0.003 |
| 21-40 | 97.93 (97.8-98.1); 27762 | 99.19 (99.1-99.3); 38184 | 0.0001, 0.007 |
| >40 | 98.16 (98.0-98.4); 18030 | 99.37 (99.3-99.5); 22213 | 0.0002, 0.025 |

| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |
| --- | ---: | --- | --- |
| all | 204556 | 194862/2026/7668/0 | 98.97 (98.9-99.0) |
| id, all | 64668 | 61197/752/2719/0 | 98.79 (98.7-98.9) |
| id, best | 129336 | 123597/1188/4551/0 | 99.05 (99.0-99.1) |
| name, all | 3075 | 2795/40/240/0 | 98.59 (98.1-99.0) |
| name, best | 6150 | 5946/46/158/0 | 99.23 (99.0-99.4) |
| name db, best | 1327 | 1327/0/0/0 | 100.00 (99.7-100.0) |
| travel_born_late | 20706 | 20561/0/145/0 | 100.00 (100.0-100.0) |
| travel_handoff | 51370 | 50310/119/941/0 | 99.76 (99.7-99.8) |
| travel_merged | 25695 | 25150/65/480/0 | 99.74 (99.7-99.8) |
| travel_name | 9225 | 8741/86/398/0 | 99.03 (98.8-99.2) |
| travel_renamed | 22635 | 22392/0/243/0 | 100.00 (100.0-100.0) |
| travel_retired_none | 24228 | 17067/1700/5461/0 | 90.94 (90.5-91.3) |
| travel_split | 24972 | 24916/56/0/0 | 99.78 (99.7-99.8) |
| travel_unchanged | 25725 | 25725/0/0/0 | 100.00 (100.0-100.0) |

| check (no truth needed) | checked | violations |
| --- | ---: | ---: |
| identity | 2853 | 0 |
| composition | 3000 | 143 |
| return | 3000 | 0 |

| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness outcome of the direct query and the first leg | n |
| --- | --- | --- | --- | --- | ---: |
| composition | forward then backward | handoff between | disjoint | direct right, first leg right | 43 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg right | 28 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg right | 10 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg wrong | 8 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg wrong | 8 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg undetermined | 8 |
| composition | forward then forward | handoff between | disjoint | direct right, first leg right | 7 |
| composition | forward then backward | handoff between | composed has fewer | direct right, first leg right | 4 |
| composition | forward then backward | one assembly | disjoint | direct right, first leg right | 4 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg undetermined | 3 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg undetermined | 3 |
| composition | forward then forward | one assembly | composed has fewer | direct right, first leg right | 2 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg undetermined | 2 |
| composition | forward then backward | handoff between | disjoint | direct undetermined, first leg right | 2 |
| composition | backward then forward | handoff between | composed has more | direct right, first leg right | 2 |
| composition | backward then backward | handoff between | composed has fewer | direct right, first leg right | 2 |
| composition | forward then backward | one assembly | disjoint | direct right, first leg undetermined | 1 |
| composition | forward then backward | one assembly | composed has more | direct right, first leg right | 1 |
| composition | backward then forward | handoff between | disjoint | direct undetermined, first leg right | 1 |
| composition | backward then backward | one assembly | disjoint | direct wrong, first leg wrong | 1 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg wrong | 1 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg right | 1 |
| composition | forward then backward | handoff between | composed has fewer | direct right, first leg undetermined | 1 |
