# Travel matrix: macaca_fascicularis

rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout

| origin era | direction | within assembly | across a handoff |
| --- | --- | --- | --- |
| Macaca_fascicularis_5.0 | backward | 100.00 (99.9-100.0); 7218, 7218/0/0/0 |  |
| Macaca_fascicularis_5.0 | same | 100.00 (99.7-100.0); 1425, 1425/0/0/0 |  |
| Macaca_fascicularis_5.0 | forward | 100.00 (100.0-100.0); 8457, 8457/0/0/0 | 100.00 (100.0-100.0); 18903, 18188/0/715/0 |
| Macaca_fascicularis_6.0 | backward | 100.00 (99.9-100.0); 5058, 5058/0/0/0 | 100.00 (100.0-100.0); 9900, 9720/0/180/0 |
| Macaca_fascicularis_6.0 | same | 100.00 (99.6-100.0); 1077, 1077/0/0/0 |  |
| Macaca_fascicularis_6.0 | forward | 100.00 (99.9-100.0); 7362, 7362/0/0/0 |  |

| distance (releases) | backward | forward | seconds per query: median, 99th percentile |
| --- | --- | --- | --- |
| 0 | 100.00 (99.8-100.0); 2502 | 100.00 (99.8-100.0); 2502 | 0.0001, 0.000 |
| 1-2 | 100.00 (99.9-100.0); 3777 | 100.00 (99.9-100.0); 4554 | 0.0001, 0.001 |
| 3-5 | 100.00 (99.9-100.0); 5229 | 100.00 (99.9-100.0); 6636 | 0.0001, 0.001 |
| 6-10 | 100.00 (99.9-100.0); 6345 | 100.00 (100.0-100.0); 10035 | 0.0001, 0.001 |
| 11-20 | 100.00 (99.9-100.0); 5850 | 100.00 (100.0-100.0); 11727 | 0.0001, 0.001 |
| 21-40 | 100.00 (99.6-100.0); 975 | 100.00 (99.8-100.0); 1770 | 0.0001, 0.001 |
| >40 |  |  |  |

| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |
| --- | ---: | --- | --- |
| all | 59400 | 58505/0/895/0 | 100.00 (100.0-100.0) |
| id, all | 18750 | 18430/0/320/0 | 100.00 (100.0-100.0) |
| id, best | 37500 | 36925/0/575/0 | 100.00 (100.0-100.0) |
| name, all | 825 | 825/0/0/0 | 100.00 (99.5-100.0) |
| name, best | 1650 | 1650/0/0/0 | 100.00 (99.8-100.0) |
| name db, best | 675 | 675/0/0/0 | 100.00 (99.4-100.0) |
| travel_born_late | 9510 | 9330/0/180/0 | 100.00 (100.0-100.0) |
| travel_handoff | 9510 | 9510/0/0/0 | 100.00 (100.0-100.0) |
| travel_merged | 9510 | 9510/0/0/0 | 100.00 (100.0-100.0) |
| travel_name | 2475 | 2475/0/0/0 | 100.00 (99.8-100.0) |
| travel_retired_none | 9375 | 8660/0/715/0 | 100.00 (100.0-100.0) |
| travel_split | 9510 | 9510/0/0/0 | 100.00 (100.0-100.0) |
| travel_unchanged | 9510 | 9510/0/0/0 | 100.00 (100.0-100.0) |

| check (no truth needed) | checked | violations |
| --- | ---: | ---: |
| identity | 2250 | 0 |
| composition | 3000 | 218 |
| return | 3000 | 0 |

| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness outcome of the direct query and the first leg | n |
| --- | --- | --- | --- | --- | ---: |
| composition | forward then backward | handoff between | disjoint | direct right, first leg right | 117 |
| composition | forward then backward | handoff between | composed has more | direct right, first leg right | 70 |
| composition | forward then backward | handoff between | disjoint | direct right, first leg undetermined | 16 |
| composition | forward then backward | handoff between | composed has fewer | direct right, first leg right | 6 |
| composition | backward then forward | handoff between | disjoint | direct right, first leg right | 4 |
| composition | backward then forward | handoff between | composed has more | direct right, first leg right | 3 |
| composition | forward then forward | handoff between | composed has fewer | direct right, first leg right | 2 |
