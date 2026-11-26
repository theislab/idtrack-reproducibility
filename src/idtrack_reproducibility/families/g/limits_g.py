"""G's rows of the known-limitations table (docs/methods.md, "Known limitations"), from the reading of its sample.

Usage: python -m idtrack_reproducibility.families.g.limits_g --out <tsv> <scored folder of idtrack> [...]

Family G reads a seeded sample of its disagreements with idtrack (`families/g/reading.py`): every sampled line where
idtrack is wrong has a class, either a known limitation (`L-...: <cause>`, the limitation the families' causes or
G's reading name) or a new cause (`new: <class> (...)`). One row per limitation and organism: `answers` the sampled
lines where idtrack is wrong in that class (each line once, whichever tool it disagreed with), `units` their inputs,
`of_answers` every sampled line where idtrack is wrong in that organism, `share` the one over the other: a measured
size in G's sample, not in G's every line. A new cause becomes the row of G's map (`limitations_map_g.tsv`: the
first row whose `strata` matches the line's stratum and whose `when` is the reading's class, `*` any class), stated
as what it is today (cause read or not read, with its backlog row). The texts of a known limitation are its
families' (F0's map, the families' maps, else another family's row of the table); the scored folders give the run,
commit and graph digest.
"""
import argparse
import collections
import csv
import glob
import json
import os
import re

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import limitations as LIM

HERE = os.path.dirname(os.path.abspath(__file__))
MAP = os.path.join(HERE, "limitations_map_g.tsv")
TEXTS = ("group", "title", "what", "user_action", "fixed_by")


def read_map(path=MAP):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    for r in rows:
        r["rx"] = re.compile(r["strata"])
    return rows


def known_texts():
    """{limitation: its texts} from F0's map and every family's map of limitations (the first that states it)."""
    out = {}
    for lid, *texts in LIM.F0_MAP.values():
        out.setdefault(lid, dict(zip(TEXTS, texts)))
    from idtrack_reproducibility.check import ROOT
    table = os.path.join(ROOT, "results", "limitations", "limitations.tsv")  # the rows other families' builders state
    for p in sorted(glob.glob(os.path.join(os.path.dirname(HERE), "*", "limitations_map*.tsv"))) + [table]:
        if not os.path.exists(p):
            continue
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                if r.get("limitation") and r.get("family") != "g":
                    out.setdefault(r["limitation"], {k: r.get(k, "") for k in TEXTS})
    return out


def classify(cls, stratum, gmap):
    """(limitation, cause) of a sampled idtrack-wrong line's class: a known limitation as the reading names it, a new
    cause by G's map."""
    if cls.startswith("L-") or (":" in cls and not cls.startswith("new: ")):
        lid, _, cause = cls.partition(": ")
        return lid, cause
    k = re.sub(r" \([^()]*in the families' runs\)$", "", cls[len("new: "):] if cls.startswith("new: ") else cls)
    for r in gmap:
        if r["rx"].search(stratum or "") and r["when"] in (k, "*"):
            return r["limitation"], k
    return "UNMAPPED", k


def rows(reading, scored, gmap=None, texts=None):
    gmap = read_map() if gmap is None else gmap
    texts = known_texts() if texts is None else texts
    by_map = {r["limitation"]: r for r in gmap}
    info = {}
    for d in scored:
        j = json.load(open(os.path.join(d, "run.json")))
        info[j["organism"]] = (d, j)
    seen, lines = set(), []
    with open(reading, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if r["side"] != "idtrack" or (r["organism"], r["qid"]) in seen:
                continue
            seen.add((r["organism"], r["qid"]))
            lines.append(r)
    of = collections.Counter(r["organism"] for r in lines)
    agg = collections.defaultdict(lambda: [set(), set(), set()])
    for r in lines:
        lid, cause = classify(r["cls"], r["stratum"], gmap)
        a = agg[(lid, r["organism"], cause)]
        a[0].add(r["qid"])
        a[1].add(r["qid"].split("|")[1] if "|" in r["qid"] else r["qid"])
        a[2].add(r["stratum"])
    out = []
    for (lid, org, cause), (qids, units, strata) in sorted(agg.items()):
        t = by_map.get(lid) or texts.get(lid) or dict(title="UNMAPPED")
        d, j = info.get(org, ("", {}))
        commit = (j.get("commit") or "cabfd9a")[:12]
        out.append(dict(limitation=lid, group=t.get("group", ""), title=t.get("title", ""), what=t.get("what", ""),
                        organism=org, answers=len(qids), units=len(units), of_answers=of[org],
                        share=f"{len(qids) / of[org]:.5f}", family="g", strata=",".join(sorted(strata)),
                        run=LIM.where(d) if d else "", commit=commit,
                        graph_digest=",".join(x for x in (j.get("graph_digests") or []) if x and x != "None"),
                        **{"class": by_map.get(lid, {}).get("class", "")}, cause=cause,
                        rule=by_map.get(lid, {}).get("rule", ""), user_action=t.get("user_action", ""),
                        fixed_by=t.get("fixed_by", ""),
                        measured_by=f"g: the reading of its sample (report/reading.tsv), {of[org]} sampled lines "
                                    f"where idtrack is wrong", since=commit))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.limits_g", description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("scored", nargs="+", help="idtrack's scored folders of G (their run.json)")
    a = ap.parse_args(argv)
    reading = os.path.join(C.family_dir("g"), "report", "reading.tsv")
    out = rows(reading, [d for d in a.scored if os.path.exists(os.path.join(d, "run.json"))])
    LIM.write(out, a.out)
    print(f"{len(out)} rows -> {a.out}")


if __name__ == "__main__":
    main()
