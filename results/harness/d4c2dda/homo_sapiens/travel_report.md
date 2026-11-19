# Travel matrix: homo_sapiens

rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout

| origin era | direction | within assembly | across a handoff |
| --- | --- | --- | --- |
| NCBI36 | backward | 100.00 (99.6-100.0); 876, 866/0/10/0 |  |
| NCBI36 | same | 100.00 (99.2-100.0); 473, 473/0/0/0 |  |
| NCBI36 | forward | 100.00 (99.8-100.0); 1962, 1962/0/0/0 | 99.15 (99.0-99.3); 29095, 25697/221/3177/0 |
| GRCh37 | backward | 92.35 (91.8-92.9); 14306, 8297/687/5322/0 | 92.59 (91.9-93.2); 8169, 6124/490/1555/0 |
| GRCh37 | same | 100.00 (99.7-100.0); 1408, 1408/0/0/0 |  |
| GRCh37 | forward | 99.79 (99.7-99.9); 9145, 9076/19/50/0 | 99.71 (99.7-99.8); 56156, 55838/160/158/0 |
| GRCh38 | backward | 99.89 (99.8-99.9); 27095, 26248/30/817/0 | 99.20 (99.1-99.3); 32753, 30994/250/1509/0 |
| GRCh38 | same | 100.00 (99.7-100.0); 1205, 1205/0/0/0 |  |
| GRCh38 | forward | 100.00 (100.0-100.0); 18805, 18805/0/0/0 |  |

| distance (releases) | backward | forward | seconds per query: median, 99th percentile |
| --- | --- | --- | --- |
| 0 | 100.00 (99.9-100.0); 3086 | 100.00 (99.9-100.0); 3086 | 0.0001, 0.000 |
| 1-2 | 99.08 (98.8-99.3); 5154 | 99.85 (99.7-99.9); 5486 | 0.0001, 0.001 |
| 3-5 | 98.08 (97.7-98.4); 7409 | 99.75 (99.6-99.8); 8060 | 0.0001, 0.001 |
| 6-10 | 97.42 (97.1-97.7); 11466 | 99.68 (99.6-99.8); 13262 | 0.0001, 0.001 |
| 11-20 | 96.91 (96.6-97.2); 20264 | 99.63 (99.5-99.7); 24787 | 0.0001, 0.001 |
| 21-40 | 98.36 (98.2-98.5); 23993 | 99.60 (99.5-99.7); 42317 | 0.0001, 0.001 |
| >40 | 98.87 (98.7-99.0); 14913 | 99.63 (99.5-99.7); 21251 | 0.0001, 0.002 |

| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |
| --- | ---: | --- | --- |
| all | 201448 | 186993/1857/12598/0 | 99.02 (99.0-99.1) |
| id, all | 63512 | 58520/649/4343/0 | 98.90 (98.8-99.0) |
| id, best | 127296 | 118031/1168/8097/0 | 99.02 (99.0-99.1) |
| name, all | 3075 | 2996/0/79/0 | 100.00 (99.9-100.0) |
| name, best | 6150 | 6031/40/79/0 | 99.34 (99.1-99.5) |
| name db, best | 1415 | 1415/0/0/0 | 100.00 (99.7-100.0) |
| travel_born_late | 23067 | 21339/131/1597/0 | 99.39 (99.3-99.5) |
| travel_handoff | 51435 | 50781/105/549/0 | 99.79 (99.8-99.8) |
| travel_merged | 23305 | 20535/584/2186/0 | 97.23 (97.0-97.4) |
| travel_name | 9225 | 9027/40/158/0 | 99.56 (99.4-99.7) |
| travel_renamed | 24010 | 19014/484/4512/0 | 97.52 (97.3-97.7) |
| travel_retired_none | 21095 | 17613/215/3267/0 | 98.79 (98.6-98.9) |
| travel_split | 23586 | 22959/298/329/0 | 98.72 (98.6-98.9) |
| travel_unchanged | 25725 | 25725/0/0/0 | 100.00 (100.0-100.0) |

| check (no truth needed) | checked | violations |
| --- | ---: | ---: |
| identity | 2806 | 0 |
| composition | 3000 | 204 |
| return | 3000 | 0 |

| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness outcome of the direct query and the first leg | n |
| --- | --- | --- | --- | --- | ---: |
| composition | forward then backward | handoff between | disjoint | direct right, first leg right | 88 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg undetermined | 21 |
| composition | forward then forward | handoff between | disjoint | direct right, first leg right | 20 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg undetermined | 13 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg wrong | 11 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg right | 10 |
| composition | forward then backward | handoff between | composed has fewer | direct right, first leg right | 7 |
| composition | forward then backward | one assembly | disjoint | direct right, first leg right | 7 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg undetermined | 6 |
| composition | forward then backward | handoff between | disjoint | direct wrong, first leg right | 4 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg wrong | 3 |
| composition | forward then backward | handoff between | disjoint | direct undetermined, first leg right | 2 |
| composition | backward then forward | handoff between | disjoint | direct undetermined, first leg undetermined | 2 |
| composition | backward then backward | one assembly | disjoint | direct undetermined, first leg wrong | 2 |
| composition | backward then forward | handoff between | disjoint | direct undetermined, first leg wrong | 2 |
| composition | forward then forward | handoff between | composed has fewer | direct undetermined, first leg right | 1 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg undetermined | 1 |
| composition | backward then forward | one assembly | disjoint | direct right, first leg wrong | 1 |
| composition | backward then forward | handoff between | disjoint | direct wrong, first leg undetermined | 1 |
| composition | forward then backward | one assembly | composed has more | direct right, first leg right | 1 |
| composition | backward then backward | handoff between | disjoint | direct wrong, first leg undetermined | 1 |
