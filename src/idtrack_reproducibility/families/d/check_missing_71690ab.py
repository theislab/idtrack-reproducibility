"""A small check of missing values in lists (rule `missing_values`) on 71690ab (pig graph, memory-mapped): convert_many on lists with missing-value strings and with
a blank string, item by item. Usage: check_missing_71690ab.py <clone> <graph file> <out json>"""
import json
import os
import sys

from idtrack_reproducibility.families import adapter as F0

clone, graph, out = sys.argv[1:4]
ad = F0.Adapter(clone, graph, "sus_scrofa", 115, os.path.dirname(out))
api = ad.api
cases = {
    "nan_inside": ["TP53", "nan", "ENSSSCG00000016117"],
    "NA_None_null": ["NA", "INS", "None", "null", "GCG"],
    "blank_inside": ["INS", "", "GCG"],
    "space_inside": ["INS", " ", "GCG"],
}
res = {}
for k, batch in cases.items():
    try:
        r = api.convert_many(batch, to_release=115)
        res[k] = [dict(item=x, status=y.status, targets=list(y.targets("best")), error=(y.error or {}).get("type") if isinstance(y.error, dict) else str(y.error) if y.error else None)
                  for x, y in zip(batch, r)]
    except Exception as e:  # noqa: BLE001
        res[k] = dict(raised=type(e).__name__, message=str(e)[:200])
json.dump(res, open(out, "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
