"""B's rows of the known-limitations table (docs/methods.md, "Known limitations"), from its runs (judge `now`).

Usage: limits_b.py --out <tsv> <run folder> [<run folder> ...]   (`<run>/scores_now` scored, `<run>/report_b.tsv`)

- the wrong answers, by B's map of their class (`limitations_map_outcomes.tsv`; `limits_outcomes`);
- `L-b2-history-ends`: B2's moved names read at an old release whose old gene Ensembl retired without a successor
  (`b2_moved_from: history ends` of `report_b.tsv`): every answer is a bridge, which the truth judges by place where
  the two assemblies can be compared, and otherwise leaves undetermined; its size and its outcomes, not a rate;
- `L-mgi-havana-symbols`: B8's wrong answers of a gene converted to `MGI Symbol` at a release where the symbol MGI gave
  it was a Havana id (`OTTMUSG…`; MGI's `MRK_List1.rpt` holds each as a withdrawn MGI symbol): idtrack files such a
  label under `Havana gene` (its rule of identifier schemes, `idtrack.build.continuity`, rule 3), so the gene has no
  MGI symbol there (row I170).
"""
import argparse
import csv
import json
import os

from idtrack_reproducibility.families import limitations as LIM
from idtrack_reproducibility.families import limits_outcomes as LO

HERE = os.path.dirname(os.path.abspath(__file__))
ENDS = ("L-b2-history-ends", "history", "a moved name whose old gene Ensembl retired without a successor",
        "a name read at an old release whose carrier then Ensembl retired with no successor, and that another gene "
        "carries at 115: the answer key cannot judge a name whose gene ended (every answer is a bridge: judged by place "
        "where the assemblies can be compared, else undetermined); reported apart from B2's measured part",
        "read `result.path.bridges`; ask at the target release for the name's carrier there",
        "none (Ensembl's history ends the gene)")


def ends_rows(run, scored):
    info = json.load(open(os.path.join(scored, "run.json")))
    rep = {r["stratum"]: r for r in csv.DictReader(open(os.path.join(run, "report_b.tsv")), delimiter="\t")
           if r["judge"] == "now"}
    r, whole = rep.get("b2_moved_from: history ends"), rep.get("b2_moved_from")
    if not r or not int(r["n"]):
        return []
    n, of = int(r["n"]), int(whole["n"]) if whole else 0
    commit = (info.get("commit") or "")[:12]
    return [dict(limitation=ENDS[0], group=ENDS[1], title=ENDS[2], what=ENDS[3], organism=info["organism"], answers=n,
                 units=r["units"], of_answers=of, share=f"{n / of:.5f}" if of else "", family="b", strata="b2_moved_from",
                 run=LIM.where(scored), commit=commit,
                 graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                 **{"class": "d"}, cause=f"right {r['right']}, wrong {r['wrong']}, undetermined {r['undetermined']}",
                 rule="-", user_action=ENDS[4], fixed_by=ENDS[5],
                 measured_by=f"b: report_b.tsv of {LIM.where(run)} (b2_moved_from: history ends)", since=commit)]


HAVANA = ("L-mgi-havana-symbols", "external identifiers", "MGI's symbols that were Havana ids are filed as Havana ids",
          "at mouse releases 48 to 54 MGI gave genes without a symbol their Havana id (`OTTMUSG…`) as their symbol "
          "(MGI's records hold each as a withdrawn symbol: `OTTMUSG00000000990` of MGI:3650286, now Heatr9); idtrack "
          "files these labels under `Havana gene` (rule 3 of `idtrack.build.continuity`), so converted to `MGI Symbol` "
          "such a gene answers no symbol there (row I170)",
          "convert to `Havana gene` for these releases, or to `MGI Symbol` at a later release",
          "none (the package's filing of identifier schemes; no round 3)")


def havana_rows(run, scored):
    """B8's wrong answers whose expected MGI symbol is a Havana id (`OTTMUSG…`)."""
    import glob
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    info = json.load(open(os.path.join(scored, "run.json")))
    n, units = 0, set()
    for p in sorted(glob.glob(os.path.join(scored, "outcomes", "*.parquet"))):
        t = pq.read_table(p, columns=["stratum", "unit", "outcome", "expected"])
        t = t.filter(pc.and_(pc.equal(t.column("stratum"), "b8_core_xref"), pc.equal(t.column("outcome"), "wrong")))
        for u, e in zip(t.column("unit").to_pylist(), t.column("expected").to_pylist()):
            if any(x.startswith("OTTMUSG") for x in (e or "").split(",")):
                n += 1
                units.add(u)
    if not n:
        return []
    of = 0
    for r in csv.DictReader(open(os.path.join(scored, "report.tsv")), delimiter="\t"):
        if r["stratum"] == "b8_core_xref" and r["view"] == "*":
            of = int(r["n"])
    commit = (info.get("commit") or "")[:12]
    return [dict(limitation=HAVANA[0], group=HAVANA[1], title=HAVANA[2], what=HAVANA[3], organism=info["organism"],
                 answers=n, units=len(units), of_answers=of, share=f"{n / of:.5f}" if of else "", family="b",
                 strata="b8_core_xref", run=LIM.where(scored), commit=commit,
                 graph_digest=",".join(x for x in (info.get("graph_digests") or []) if x and x != "None"),
                 **{"class": "a"}, cause="expected an OTTMUSG id as MGI Symbol", rule="idtrack.build.continuity rule 3",
                 user_action=HAVANA[4], fixed_by=HAVANA[5], measured_by=f"b: outcomes of {LIM.where(scored)} (b8_core_xref)",
                 since=commit)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    scored = [os.path.join(r, "scores_now") for r in a.runs]
    rows = LO.build("b", scored, LO.load_map(os.path.join(HERE, "limitations_map_outcomes.tsv")))
    for r, s in zip(a.runs, scored):
        rows += ends_rows(r, s) + havana_rows(r, s)
    LIM.write(rows, a.out)
    print(f"{len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
