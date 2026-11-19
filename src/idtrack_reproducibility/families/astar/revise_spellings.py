"""A* revision set: the corrected judge's `spellings_once` for the name-database lines (`astar_db`) of the sets published
before gen.py wrote it (pig, rhesus and cynomolgus `v1`). Holds only the lines that change, under their qids.
Usage: revise_spellings.py <organism> --set <new name> --revises <set>"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.astar import gen as GN
from idtrack_reproducibility.families import lines as L

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--revises", required=True)
    a = ap.parse_args()
    out = []
    for q in L.iter_set(L._set_file(L.query_dir("astar", a.organism), a.revises)):
        if q["kind"] != "astar_db":
            continue
        sp = GN.spellings_of(q["expect"])
        if sp and q.get("spellings") != sp:
            q = dict(q)
            q["spellings"] = sp
            q["revised_by"] = list(q.get("revised_by") or []) + ["spellings_once"]
            out.append(q)
    print(a.organism, len(out), "lines")
    if out:
        man = L.publish("astar", a.organism, a.set, out, "the corrected judge spellings_once on the astar_db lines of " + a.revises,
                        os.path.abspath(__file__), revises=[a.revises], judge="78")
        print(man["n"], man["sha256"])


if __name__ == "__main__":
    main()
