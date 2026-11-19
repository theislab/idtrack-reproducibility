"""D1's own report of a scored run (beside F0's report.tsv): rates by features and by counts, per dataset, the place
against Ensembl's history, the wrong answers by what they are, and rosetta's own mapping judged by the same truth.

Usage: report_d1.py <scored folder> --organism <o> [--sets d1_v1]

Reads F0's outcomes (`outcomes/*.parquet`) and the family's lines (lines.load), writes into `<scored folder>/d1/`:

  rates.tsv            stratum, view, n, right, wrong, undetermined, error, timeout, rate, and the same weighted by the
                       counts the features hold in the authors' matrices (`*_counts`, rate_counts)
  datasets.tsv         per dataset and view: the same, over every stratum (scored lines) and over the strict strata
                       (`one` only)
  history_vs_locus.tsv per kind and truth class of the lines with an expected gene: Ensembl's strict history walk
                       R -> 115 equals the place's gene / has no successor / names another gene
  wrong_classes.tsv    per stratum and view, the wrong answers by what they are (below), lines, units, counts
  wrong_examples.tsv   up to 20 lines per class and stratum, the highest counts first
  rosetta.tsv          rosetta's own mapping by its current rule (gen_d1.rosetta_current: the r3 place step's run
                       `out/m16` where the source shares the target's assembly, else the id route `out/m02`) against
                       the same truth, per stratum; a feature rosetta leaves unmapped because its gene is outside the
                       annotation rosetta counts (`not_in_annotation`) is apart (`unmapped_not_counted`: rosetta's
                       design)
  rosetta_lost.tsv     the features rosetta leaves unmapped where the place names one gene it counts, by rosetta's class
  rosetta_by_dataset.tsv  the same per dataset (xiang2020's versioned lines left out: rosetta asks them unversioned)

The classes of a wrong answer (`wclass`): `follows_history_elsewhere` (the answer is Ensembl's strict history walk,
which leads to genes other than the place's), `follows_history_split` (the answer is the history's walk, several genes
of which one is the place's: a split the place does not show), `history_member_not_place` (a part of the history's
genes, not the place's), `history_part` (a part of them with the place's gene and another), `history_plus_other` (the
history's genes and more), `bridge_where_no_successor` (Ensembl's history has no successor; the answer is a gene
elsewhere: rule `bridge_same_place`), `nothing` (no answer where the place names a gene), `outside_open` (an open line answered by
a gene outside the place), `name_other_gene` (a symbol line answered by a gene other than the place's), `other`.
"""
import argparse
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families import report as RP


def base_set(ids):
    return sorted({C.base_id(x) for x in (ids or [])})


def wclass(line, got, reason):
    E = set(line.get("expect") or [])
    H = set(line.get("history") or [])
    g = set(got)
    if line.get("open"):
        return "outside_open"
    if not g:
        return "nothing"
    if line["kind"] == "d1_symbol":
        return "name_other_gene"
    if line.get("history") is None:
        return "other"
    if not H:
        return "bridge_where_no_successor"
    if g == H:
        return "follows_history_split" if E <= H else "follows_history_elsewhere"
    if g < H:
        return "history_member_not_place" if not (g & E) else "history_part"
    if H < g:
        return "history_plus_other"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scored")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="d1_v1")
    a = ap.parse_args()
    out = os.path.join(a.scored, "d1")
    os.makedirs(out, exist_ok=True)
    lines, _, _ = L.load("d", a.organism, a.sets.split(","))
    by = {q["qid"]: q for q in lines if q["qid"].startswith("d1|")}
    from idtrack_reproducibility.families.d import gen_d1
    cur = gen_d1.rosetta_current(a.organism)
    for q in by.values():
        kind = "symbol" if q["kind"] == "d1_symbol" else ("id_versioned" if q["qid"].split("|")[2] == "id_versioned" else "id")
        tg, cl, which = cur.get((q["dataset"], kind, q["id"], q["source_release"]), (set(), set(), "none"))
        q["rosetta_now"], q["rosetta_now_class"], q["rosetta_now_from"] = sorted(tg), sorted(cl - {""}), which
    df = RP.read_outcomes(a.scored)
    df = df[df["qid"].isin(by.keys())]
    rates = collections.defaultdict(collections.Counter)
    ds = collections.defaultdict(collections.Counter)
    wc = collections.defaultdict(collections.Counter)
    wunits = collections.defaultdict(set)
    ex = collections.defaultdict(list)
    ros = collections.defaultdict(collections.Counter)
    hist = collections.Counter()
    ros_lost = collections.defaultdict(collections.Counter)
    ros_ds = collections.defaultdict(collections.Counter)
    ros_ex = collections.defaultdict(list)
    for qid, view, outcome, got, reason in zip(df["qid"], df["view"], df["outcome"], df["got"], df["reason"]):
        q = by[qid]
        w = int(q.get("counts") or 0)
        st = q["stratum"]
        for key in ((st, view), (st, "*")):
            rates[key][outcome] += 1
            rates[key][outcome + "_counts"] += w
        strict = st.endswith("_one")
        for key in ((q["dataset"], view, "scored"), (q["dataset"], view, "strict")):
            if key[2] == "strict" and not strict:
                continue
            ds[key][outcome] += 1
            ds[key][outcome + "_counts"] += w
        if outcome == "wrong":
            g = base_set(got.split(",") if isinstance(got, str) and got else list(got) if got is not None and not isinstance(got, str) else [])
            c = wclass(q, g, reason)
            wc[(st, view, c)]["lines"] += 1
            wc[(st, view, c)]["counts"] += w
            wunits[(st, view, c)].add(q.get("unit"))
            ex[(st, c)].append((w, view, q, g, reason))
        if view == "best":
            if q["expect"] and q["kind"] == "d1_id" and q.get("history") is not None:
                H, E = set(q["history"]), set(q["expect"])
                hist[(st, q["truth_class"], "same" if H == E else "no_successor" if not H else
                      "holds_place_gene" if E <= H else "another_gene")] += 1
            rt = set(q.get("rosetta_now") or [])
            E = set(q.get("expect") or [])
            outside = "not_in_annotation" in q.get("rosetta_now_class", [])
            if q.get("open"):
                rv = "open_mapped" if rt else "open_unmapped"
                if rt and q.get("open_within") is not None and not rt <= set(q["open_within"]):
                    rv = "open_mapped_outside"
            elif not rt:
                rv = "unmapped_not_counted" if outside else "unmapped_where_place_names_a_gene"
            else:
                # a name of several genes at R: any non-empty part of them is right (the `set` rule of `best`)
                rv = "right" if rt == E or (q.get("truth_class") == "name_several" and rt <= E) else "wrong"
            if "versioned" not in st:
                rd = ros_ds[q["dataset"]]
                rd[rv] += 1
                rd[rv + "_counts"] += w
            if rv == "unmapped_where_place_names_a_gene":
                for k in q.get("rosetta_now_class") or ["?"]:
                    ros_lost[(st, k)]["lines"] += 1
                    ros_lost[(st, k)]["counts"] += w
                    ros_ex[(st, k)].append((w, q))
            ros[st][rv] += 1
            ros[st][rv + "_counts"] += w
    cols = ["right", "wrong", "undetermined", "error", "timeout"]

    def row(c):
        den = c["right"] + c["wrong"] + c["error"]
        denw = c["right_counts"] + c["wrong_counts"] + c["error_counts"]
        return ([sum(c[k] for k in cols)] + [c[k] for k in cols] + [f"{c['right'] / den:.5f}" if den else ""]
                + [c[k + "_counts"] for k in cols] + [f"{c['right_counts'] / denw:.5f}" if denw else ""])
    head = ["n"] + cols + ["rate"] + [k + "_counts" for k in cols] + ["rate_counts"]
    with open(os.path.join(out, "rates.tsv"), "w") as fo:
        fo.write("\t".join(["stratum", "view"] + head) + "\n")
        for (st, v), c in sorted(rates.items()):
            fo.write("\t".join(map(str, [st, v] + row(c))) + "\n")
    with open(os.path.join(out, "datasets.tsv"), "w") as fo:
        fo.write("\t".join(["dataset", "view", "lines"] + head) + "\n")
        for (d, v, which), c in sorted(ds.items()):
            fo.write("\t".join(map(str, [d, v, which] + row(c))) + "\n")
    with open(os.path.join(out, "history_vs_locus.tsv"), "w") as fo:
        fo.write("stratum\ttruth_class\thistory\tlines\n")
        for k, n in sorted(hist.items()):
            fo.write("\t".join(map(str, list(k) + [n])) + "\n")
    with open(os.path.join(out, "wrong_classes.tsv"), "w") as fo:
        fo.write("stratum\tview\twclass\tlines\tunits\tcounts\n")
        for k, c in sorted(wc.items()):
            fo.write("\t".join(map(str, list(k) + [c["lines"], len(wunits[k]), c["counts"]])) + "\n")
    with open(os.path.join(out, "wrong_examples.tsv"), "w") as fo:
        fo.write("stratum\twclass\tview\tqid\tcounts\texpect\tgot\thistory\ttruth_class\tlocus_hits\trosetta_target\treason\n")
        for (st, c), rows in sorted(ex.items()):
            for w, v, q, g, reason in sorted(rows, key=lambda r: -r[0])[:20]:
                fo.write("\t".join(map(str, [st, c, v, q["qid"], w, ",".join(q["expect"]), ",".join(g), ",".join(q.get("history") or []),
                                             q["truth_class"], json.dumps(q.get("locus_hits")), ",".join(q.get("rosetta_target") or []),
                                             reason])) + "\n")
    with open(os.path.join(out, "rosetta_by_dataset.tsv"), "w") as fo:
        keys = ["right", "wrong", "unmapped_where_place_names_a_gene", "unmapped_not_counted", "open_mapped", "open_mapped_outside", "open_unmapped"]
        fo.write("\t".join(["dataset"] + keys + [k + "_counts" for k in keys]) + "\n")
        for d, c in sorted(ros_ds.items()):
            fo.write("\t".join(map(str, [d] + [c[k] for k in keys] + [c[k + "_counts"] for k in keys])) + "\n")
    with open(os.path.join(out, "rosetta_lost.tsv"), "w") as fo:
        fo.write("stratum\trosetta_class\tlines\tcounts\texamples (counts: dataset query R -> place)\n")
        for (st, k), c in sorted(ros_lost.items()):
            ex = "; ".join(f"{w}: {q['dataset']} {q['id']} {q['source_release']} -> {','.join(q['expect'])}"
                           for w, q in sorted(ros_ex[(st, k)], key=lambda r: -r[0])[:5])
            fo.write(f"{st}\t{k}\t{c['lines']}\t{c['counts']}\t{ex}\n")
    with open(os.path.join(out, "rosetta.tsv"), "w") as fo:
        keys = ["right", "wrong", "unmapped_not_counted", "unmapped_where_place_names_a_gene", "open_mapped", "open_mapped_outside",
                "open_unmapped"]
        fo.write("\t".join(["stratum"] + keys + [k + "_counts" for k in keys]) + "\n")
        for st, c in sorted(ros.items()):
            fo.write("\t".join(map(str, [st] + [c[k] for k in keys] + [c[k + "_counts"] for k in keys])) + "\n")
    print(open(os.path.join(out, "rates.tsv")).read())
    print(open(os.path.join(out, "wrong_classes.tsv")).read())


if __name__ == "__main__":
    main()
