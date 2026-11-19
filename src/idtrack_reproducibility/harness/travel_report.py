"""Tables of the travel matrix: rates by the era of the release of origin, the direction and distance of travel, and
whether an assembly handoff lies between; seconds per query by distance.

Usage: travel_report.py <scored travel run> [<organism>]   (writes <run>/travel_report.md and .tsv)

Reads the run's `outcomes.tsv.gz` and the published `travel_*` lines. Rate = right / (right + wrong + error) with a
95 % Wilson interval; undetermined and timeouts apart, as everywhere in the harness.
"""
import collections
import csv
import gzip
import json
import os
import statistics
import sys

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import score

BINS = [(0, 0, "0"), (1, 2, "1-2"), (3, 5, "3-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 40, "21-40"),
        (41, 999, ">40")]


def dbin(d):
    return next(lab for lo, hi, lab in BINS if lo <= d <= hi)


def fmt(c):
    den = c["right"] + c["wrong"] + c["error"]
    if not den:
        return "", ""
    lo, hi = H.wilson(c["right"], den)
    return f"{100 * c['right'] / den:.2f}", f"{100 * lo:.1f}-{100 * hi:.1f}"


def checks_summary(run, by_query, asm):
    """The consistency checks (travel_checks.py): counts, and the violations by the directions of their legs, with the
    harness's outcome of the direct query (R1 -> R3) and of the first leg (R1 -> R2)."""
    p = os.path.join(run, "checks", "summary.json")
    if not os.path.exists(p):
        return []
    sm = json.load(open(p))
    md = ["", "| check (no truth needed) | checked | violations |", "| --- | ---: | ---: |"]
    for k in ("identity", "composition", "return"):
        if k in sm:
            md.append(f"| {k} | {sm[k]['checked']} | {sm[k]['violations']} |")
    groups = collections.Counter()
    for line in open(os.path.join(run, "checks", "violations.jsonl")):
        v = json.loads(line)
        if v["check"] != "composition":
            groups[(v["check"], "", "", "", "")] += 1
            continue
        (i1, r1, r2, _), _, (_, _, r3, _) = v["queries"]
        legs = ("forward" if r2 > r1 else "backward") + " then " + ("forward" if r3 > r2 else "backward")
        cross = "handoff between" if len({asm.get(r1), asm.get(r2), asm.get(r3)}) > 1 else "one assembly"
        via, direct = ({x.split(".")[0] for x in a} for a in v["answers"][1:])
        rel = ("composed has more" if via > direct else "composed has fewer" if via < direct else
               "partly shared" if via & direct else "disjoint")
        groups[("composition", legs, cross, rel, f"direct {by_query.get((i1, r1, r3, 'best'), '?')}, "
                                                  f"first leg {by_query.get((i1, r1, r2, 'best'), '?')}")] += 1
    if groups:
        md += ["", "| violation | legs | assemblies | composed (R1 -> R2 -> R3) against direct (R1 -> R3) | harness "
               "outcome of the direct query and the first leg | n |", "| --- | --- | --- | --- | --- | ---: |"]
        for k, n in groups.most_common():
            md.append(f"| {' | '.join(k)} | {n} |")
    return md


def main():
    run = sys.argv[1]
    org = sys.argv[2] if len(sys.argv) > 2 else json.load(open(os.path.join(run, "run.json")))["organism"]
    lines, _ = score.load_queries(org, None, None)
    rec = {r["qid"]: r for r in lines if r["kind"].startswith("travel_")}
    asm = H.gtf_assemblies(org)
    tab = collections.defaultdict(lambda: collections.Counter())
    by_query = {}
    secs = collections.defaultdict(list)
    with gzip.open(os.path.join(run, "outcomes.tsv.gz"), "rt") as fh:
        for o in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            r = rec.get(o["qid"])
            if r is None:
                continue
            if r["final_database"] is None:
                by_query[(r["id"], r["from_release"], r["to_release"], r["strategy"])] = o["outcome"]
            a, b = r["from_release"], r["to_release"]
            d = abs(b - a)
            direction = "same" if a == b else ("forward" if b > a else "backward")
            handoff = "across" if asm.get(a) != asm.get(b) else "within"
            era = asm.get(a, "?")
            what = "name db" if r["final_database"] else ("name" if r["kind"] == "travel_name" else "id")
            for key in [("all", "", "", "", ""), ("era", era, "", "", ""), ("direction", direction, "", "", ""),
                        ("distance", direction, dbin(d), "", ""), ("handoff", direction, handoff, "", ""),
                        ("era x direction x handoff", era, direction, handoff, ""),
                        ("stratum", r["kind"], "", "", ""), ("query", what, r["strategy"], "", "")]:
                tab[key][o["outcome"]] += 1
            if o["seconds"]:
                secs[dbin(d)].append(float(o["seconds"]))
    rows = []
    for key, c in sorted(tab.items(), key=lambda kv: [str(x) for x in kv[0]]):
        rate, ci = fmt(c)
        rows.append(list(key) + [sum(c.values()), c["right"], c["wrong"], c["undetermined"], c["timeout"], c["error"], rate, ci])
    with open(os.path.join(run, "travel_report.tsv"), "w") as fo:
        fo.write("table\ta\tb\tc\td\tn\tright\twrong\tundetermined\ttimeout\terror\trate\tci95\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    md = [f"# Travel matrix: {org}", "", "rate in % (95 % Wilson interval); n, right/wrong/undetermined/timeout", ""]
    eras = [e for e in dict.fromkeys(asm[r] for r in sorted(asm) if r <= H.ORGANISMS[org][2])]
    md += ["| origin era | direction | within assembly | across a handoff |", "| --- | --- | --- | --- |"]
    for e in eras:
        for dr in ("backward", "same", "forward"):
            cells = []
            for h in ("within", "across"):
                c = tab.get(("era x direction x handoff", e, dr, h, ""))
                if not c:
                    cells.append("")
                    continue
                rate, ci = fmt(c)
                cells.append(f"{rate} ({ci}); {sum(c.values())}, {c['right']}/{c['wrong']}/{c['undetermined']}/{c['timeout']}")
            if any(cells):
                md.append(f"| {e} | {dr} | {cells[0]} | {cells[1]} |")
    md += ["", "| distance (releases) | backward | forward | seconds per query: median, 99th percentile |",
           "| --- | --- | --- | --- |"]
    for lo, hi, lab in BINS:
        cells = []
        for dr in ("backward", "forward"):
            c = tab.get(("distance", "same" if lab == "0" else dr, lab, "", ""))
            rate, ci = fmt(c) if c else ("", "")
            cells.append(f"{rate} ({ci}); {sum(c.values())}" if c else "")
        s = sorted(secs.get(lab, []))
        st = f"{statistics.median(s):.4f}, {s[int(0.99 * (len(s) - 1))]:.3f}" if s else ""
        md.append(f"| {lab} | {cells[0]} | {cells[1]} | {st} |")
    md += ["", "| stratum | n | right/wrong/undetermined/timeout | rate (95 % CI) |", "| --- | ---: | --- | --- |"]
    for key, c in sorted(tab.items()):
        if key[0] in ("stratum", "query", "all"):
            rate, ci = fmt(c)
            name = key[1] if key[0] != "query" else f"{key[1]}, {key[2]}"
            if key[0] == "all":
                name = "all"
            md.append(f"| {name} | {sum(c.values())} | {c['right']}/{c['wrong']}/{c['undetermined']}/{c['timeout']} | {rate} ({ci}) |")
    md += checks_summary(run, by_query, asm)
    open(os.path.join(run, "travel_report.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
