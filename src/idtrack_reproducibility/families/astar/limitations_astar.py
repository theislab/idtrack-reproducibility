"""A*'s rows of the known-limitations table (docs/methods.md, "Known limitations"), generated from the scored runs (no number
typed by hand).

Usage: limitations_astar.py <run name> [--tx <tx run name>]  -> <data>/limitations_astar_<run>.tsv

Columns: id, title, what, organism, rows (scored answers: line x view), lines_ids (distinct ids), rate_of_family (share
of the family's scored answers of that organism), cause, user_action, measured_by (family, stratum, run, graph
digest), since. One row per limitation and organism. Classes are read from the answers' reasons (the scorer's), so
a stage that changes them changes the table when the family is run again.
"""
import argparse
import collections
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C

DATA = os.path.join(C.BENCH, "astar")
ORGS = ["sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis"]

# limitation id -> (title, what, cause, what a user can do, test on (outcome, reason))
LIMITS = {
    "L-astar-1": ("No successor in Ensembl's history: the answer is a bridge elsewhere",
                  "an id whose history ends (retired without successor, or asked before its birth) is answered with "
                  "genes at another place, reached through names or external ids",
                  "idtrack's bridges (rule `external_bridge`); the decided rule `bridge_same_place` restricts them by place",
                  "check `ExplainedPath` for external steps; treat answers with bridges as candidates",
                  lambda o, r: o == "wrong" and r.startswith("Ensembl has no successor; not at the query's place")),
    "L-astar-2": ("Places across an assembly handoff that no chain maps cannot be compared",
                  "for an id without successor asked across an assembly handoff, a bridged answer is neither right nor "
                  "wrong when no chain maps the place: no chain file at all (pig Sscrofa10.2/11.1, rhesus MMUL_1/"
                  "Mmul_8/Mmul_10, cynomolgus 5.0/6.0), or the place's sequence is absent from the chain or under "
                  "half of it aligns (human, mouse)",
                  "Ensembl publishes no chain file for these handoffs; chains leave some sequences out",
                  "read the answer's place yourself (e.g. by liftover with an external alignment)",
                  lambda o, r: o == "undetermined" and "places not comparable" in r),
    "L-astar-3": ("Copies on haplotypes or patches without an aligned place",
                  "a copy whose place on the reference no GRC alignment gives (NCBI36 haplotypes, unaligned parts of "
                  "patches) cannot be judged by place",
                  "the GRC aligns no NCBI36 haplotype; patches carry novel sequence",
                  "none needed for reference genes; for copies, compare names",
                  lambda o, r: o == "undetermined" and "no place on the reference sequence" in r),
    "L-astar-4": ("Ensembl's two readings of the history differ",
                  "REST archive/id (or the 116 core read as REST reads it) and the walk of the release-115 events give "
                  "other genes; the line is open",
                  "Ensembl's own records disagree (e.g. events added at 116, ids absent from GTFs)",
                  "ask Ensembl REST for the id; prefer the newer release's history",
                  lambda o, r: o == "undetermined" and "REST archive/id differ" in r),
    "L-astar-5": ("Within the wider set Ensembl's history allows",
                  "the answer is a branch that left the id while it lived on (broad walk), not its own continuation",
                  "Ensembl maps an id to new ids while it continues (hunt finding hh-1)",
                  "use strategy `all` to see every branch",
                  lambda o, r: o == "undetermined" and r.startswith("within the wider set")),
}


UNDETERMINED = {"L-astar-2", "L-astar-3", "L-astar-4", "L-astar-5"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", help="the run name (runs/<run>/<organism> of family astar)")
    ap.add_argument("--out", default=None, help="the table (default: <families>/astar/limitations_astar_<run>.tsv)")
    a = ap.parse_args()
    from idtrack_reproducibility.families import limitations as LIM
    rows = []
    for org in ORGS:
        R = os.path.join(DATA, "runs", a.run, org)
        if not os.path.exists(os.path.join(R, "report.tsv")):
            continue
        total = 0
        with open(os.path.join(R, "report.tsv")) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(head, line.rstrip("\n").split("\t")))
                if f["view"] == "*":
                    total += int(f["n"])
        run = json.load(open(os.path.join(R, "run.json")))
        digest = ",".join(x for x in (run.get("graph_digests") or []) if x and x != "None")
        commit = (run.get("commit") or "")[:12]
        cnt, ids, strata = collections.Counter(), collections.defaultdict(set), collections.defaultdict(set)
        for name in ("wrong.tsv.gz", "undetermined.tsv.gz"):
            outcome = name.split(".")[0]
            with gzip.open(os.path.join(R, name), "rt") as fh:
                head = fh.readline().rstrip("\n").split("\t")
                ir, iq, ist = head.index("reason"), head.index("qid"), head.index("stratum")
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    for lid, (*_, test) in LIMITS.items():
                        if test(outcome, f[ir]):
                            cnt[lid] += 1
                            ids[lid].add(C.base_id(f[iq].split("|")[1]))
                            strata[lid].add(f[ist])
                            break
        for lid, (title, what, cause, action, _) in LIMITS.items():
            if cnt[lid]:
                rows.append(dict(limitation=lid, group="history", title=title, what=what, organism=org,
                                 answers=cnt[lid], units=len(ids[lid]), of_answers=total,
                                 share=f"{cnt[lid] / total:.6f}", family="astar", strata=",".join(sorted(strata[lid])),
                                 run=LIM.where(R), commit=commit, graph_digest=digest, **{"class": "e" if
                                 lid in UNDETERMINED else "a"}, cause=f"astar:{lid}", rule=cause, user_action=action,
                                 fixed_by="", measured_by=f"astar: wrong and undetermined answers of {LIM.where(R)}",
                                 since=commit))
    out = a.out or os.path.join(DATA, f"limitations_astar_{a.run}.tsv")
    LIM.write(rows, out)
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
