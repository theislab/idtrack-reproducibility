"""A*: select the hardest lineages of one organism from Ensembl's own files (docs/methods.md, A*).

Reads only the truth store (GTFs of every release up to the target, the release-115 core history tables, the assembly
chain files) through the harness's readers (`make_queries.Truth`, `harness_lib`, `assembly_map`), plus two id lists
that only add ids (never truth): the ids the weakness hunts named and, for human, the ids behind the travel matrix's
weakest cell (backward within GRCh37), as recorded for each published selection in `fixed_ids.tsv`.

Every gene id of the window gets its patterns; the selection is the union over patterns of every id at or above the
pattern's threshold. Patterns (thresholds of docs/methods.md; how each is read is in `PATTERNS`).

Usage: selection.py <organism> [--target 115]
Writes <data>/selection/<organism>.tsv.gz (one row per selected id) and <organism>.json (counts per pattern, union,
full-matrix lines, the parameters). Deterministic: no sampling.
"""
import argparse
import bisect
import collections
import gzip
import json
import os
import re
import sys
import time

from idtrack_reproducibility.harness import assembly_map
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families.astar import truth as TA
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

HERE = os.path.dirname(os.path.abspath(__file__))
# The ids a selection takes from outside Ensembl's tables (they only add ids, never truth): the stable gene ids named
# in the weakness hunts' findings (pattern `hunt_history`) and, for human, the ids behind the wrong answers of the
# travel matrix's weakest cell, backward within GRCh37 (`travel_weak_cell`). `fixed_ids.tsv` records them for every
# published selection, so that a selection is drawn again exactly; a new selection reads them from the hunts'
# findings tables of this repository (`hunts/<angle>/findings.tsv`) and the harness's current travel scoring.
FIXED_IDS = os.path.join(HERE, "fixed_ids.tsv")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def _data():
    return os.path.join(C.BENCH, "astar")

PATTERNS = {
    "lineage_ge5": "every id of a connected component of 5 or more ids of the graph of cross-id gene events (old != new, "
                   "both not null) of the release-115 stable_id_event",
    "many_events": "3 or more mapping sessions with an event of the id that is not to itself and not its birth "
                   "(to another id, from another id, or to nothing), or 3 or more distinct partner ids in cross-id "
                   "events (the second reading added: few ids have events in 3 sessions, since most die at their "
                   "first cross-id event)",
    "split_then_merged_back": "two event paths from one id (by different first steps) meet again in a later id: the id "
                              "and the ids where its branches meet",
    "merged_then_split": "an id with 2 or more predecessors and 2 or more successors over the cross-id events",
    "retired_and_back": "an event to nothing in one session, then the same stable id again in a GTF after it, or an "
                        "event from nothing in a later session (rule `resurrection`)",
    "gap": "absent from the GTF (and the core gene table, where fetched) of a release between two releases that have it",
    "gtf_lacks": "present in the core database's gene table of a release whose GTF lacks it (human 48-51, "
                 "cynomolgus 112; truth.CORE_GENES)",
    "version_jump": "the GTF version rises by more than one between consecutive releases where the id exists, or a "
                    "self event raises it by more than one",
    "version_fall": "the GTF version falls between consecutive releases, or a self event lowers it",
    "moved": "within one assembly: on another sequence or the other strand at another release; across assemblies: on "
             "two different chromosomes (GTF; sequences of alternative haplotypes and patches left out, their names "
             "change by convention; scaffold names change with every assembly, so a change across a handoff counts "
             "only between chromosomes; places across a handoff with a chain file are handoff_moved)",
    "handoff_lineage": "a member of a lineage (component) with a cross-id event in a session that crosses an assembly "
                       "handoff (the GTF releases on its two sides are on different assemblies)",
    "handoff_moved": "present on both sides of a handoff whose chain file exists, and its place at the last release "
                     "of one assembly, mapped through Ensembl's chain, does not overlap its own place at the first "
                     "release of the next (a place that does not map is not counted: the chain cannot tell)",
    "rename_chain_ge3": "a member of a chain of 3 or more ids, each with one successor and one predecessor",
    "successor_overlaps": "a cross-id event whose new id, at the first GTF release after the session, overlaps "
                          "another living gene on the same sequence and strand: the new id and the old id",
    "event_contradicts_gtf": "an event whose old id is absent from the GTF of the last release before the session takes "
                             "effect, or whose new id is absent from the GTF of the first release from it (truth.py: "
                             "a session takes effect at its new release; a mislabelled one after its old release), or "
                             "whose version differs from that GTF's version (rule `event_version`'s case); the event's ids",
    "mislabelled_session": "an id with a cross-id event (old != new, both not null) in a session whose new release is "
                           "not after its old release (pig 106 -> 104, human 114 -> 114)",
    "non_reference_lineage": "a member of a lineage (component of 2 or more ids) with a gene on a patch or an "
                             "alternative locus at some release (GTF marks them for human and mouse)",
    "many_versions": "more distinct GTF versions than the id at the 99th percentile (at most 1 % of the ids)",
    "branch_while_living": "a cross-id event from the id in a session it lives across (present in the GTF on both "
                           "sides): Ensembl maps it to another id while it continues (hunt finding hh-1); the id and "
                           "the ids it branches to",
    "hunt_history": "an id named by the weakness hunts of history and of regions (their findings.tsv and "
                    "reproducers, docs/hunts.md)",
    "non_reference_copy_ended": "a gene on a haplotype, patch or alternative locus (an assembly-exception sequence "
                                "of the release, or GTF region patch) at some release, absent at the target, whose "
                                "strict walk from its last release reaches nothing (hunt_regions R6, J1, J4)",
    "par_y": "a gene on Y inside a pseudoautosomal region (assembly exception PAR) at some release (hunt_regions R4)",
    "travel_weak_cell": "human: an id behind a wrong answer of the travel matrix backward within GRCh37 (d4c2dda)",
}


# sequences of alternative haplotypes and patches, whose names change by convention between releases (c6_QBL,
# HSCHR6_MHC_QBL, CHR_HSCHR6_MHC_QBL_CTG1): a gene on them is a copy, not a gene that moved (non_reference_lineage)
ALT_SEQ = re.compile(r"^(CHR_|HSCHR|MSCHR|HG\d|HS.*_PATCH|c\d+_|.*_PATCH$|.*_HAP\d|.*MHC)")
CHROM = re.compile(r"^(\d+[A-Z]?|[XYWZ]|MT)$")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def components(cross):
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in cross:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for x in list(parent):
        comp[find(x)].add(x)
    # a stable name per component: its smallest id
    name = {}
    for members in comp.values():
        nm = min(members)
        for x in members:
            name[x] = nm
    return name, {min(m): len(m) for m in comp.values()}


class Overlaps:
    """Genes of one release on (sequence, strand), for 'overlaps another gene' queries."""

    def __init__(self, genes):
        g = collections.defaultdict(list)
        for gid, r in genes.items():
            if r["region"] == "core_only":
                continue
            g[(r["seqname"], r["strand"])].append((int(r["start"]), int(r["end"]), gid))
        self.idx = {}
        for k, v in g.items():
            v.sort()
            self.idx[k] = ([x[0] for x in v], v)

    def others(self, seq, strand, s, e, me):
        if (seq, strand) not in self.idx:
            return []
        starts, v = self.idx[(seq, strand)]
        i = bisect.bisect_right(starts, e)
        return [gid for (a, b, gid) in v[:i] if b >= s and gid != me]


def _recorded(name):
    """The ids recorded for a published selection (pattern `none`: a selection that took none), or None for a new one."""
    rows = [r for r in H.read_tsv(FIXED_IDS) if r["selection"] == name]
    return [r for r in rows if r["pattern"] != "none"] if rows else None


def hunt_ids(org, prefix, name):
    """The stable gene ids of this organism that the weakness hunts named: as recorded for the selection `name`, else
    every id of the organism's prefix in the hunts' findings tables."""
    rec = _recorded(name)
    if rec is not None:
        return {r["id"] for r in rec if r["pattern"] == "hunt_history"}, {"fixed_ids.tsv": name}
    out, seen = set(), {}
    pat = re.compile(rf"\b{prefix}\d{{11}}\b")
    for angle in ("history", "regions"):
        f = os.path.join(REPO, "hunts", angle, "findings.tsv")
        if os.path.exists(f):
            out |= set(pat.findall(open(f).read()))
            seen[f"hunts/{angle}/findings.tsv"] = L.sha256(f)
    return out, seen


def weak_cell_ids(org, name):
    """The human ids behind the travel matrix's weakest cell (backward within GRCh37): as recorded for the selection
    `name`, else read from the harness's current travel scoring of the reference tree."""
    if org != "homo_sapiens":
        return set()
    rec = _recorded(name)
    if rec is not None:
        return {r["id"] for r in rec if r["pattern"] == "travel_weak_cell"}
    p = os.path.join(C.HARNESS, "runs", "d4c2dda", org, "travel", "outcomes.tsv.gz")
    out = set()
    with gzip.open(p, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            if f[3] != "wrong":
                continue
            q = f[0].split("|")
            if q[4] != "None":
                continue
            a, b = int(q[2]), int(q[3])
            if 55 <= b < a <= 75:
                out.add(q[1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--name", default=None, help="file name of the selection (default: the organism); a published "
                    "set names the selection file it was generated from (its sha256 in the manifest)")
    a = ap.parse_args()
    name = a.name or a.organism
    org, T = a.organism, a.target
    prefix = H.ORGANISMS[org][1]
    log("loading truth", org)
    tr = TA.TruthA(org, T)
    rels = tr.rels
    ridx = {r: i for i, r in enumerate(rels)}
    _, ev, sessions = H.core_tables(org, T)
    log("releases", rels[0], rels[-1], "events", len(ev))
    pat = collections.defaultdict(set)

    cross = [(o, n) for o, ov, n, nv, s, sc in ev if o and n and o != n]
    comp_of, comp_size = components(cross)
    size_of = {x: comp_size[c] for x, c in comp_of.items()}
    pat["lineage_ge5"] = {x for x, n in size_of.items() if n >= 5}

    sess_of = collections.defaultdict(set)
    for o, ov, n, nv, s, sc in ev:
        if o and o != n:
            sess_of[o].add(s)
        if n and o and o != n:
            sess_of[n].add(s)
    partners = collections.defaultdict(set)
    for o, n in cross:
        partners[o].add(n)
        partners[n].add(o)
    pat["many_events"] = {x for x, ss in sess_of.items() if len(ss) >= 3} | {x for x, p in partners.items() if len(p) >= 3}

    succ, pred = collections.defaultdict(set), collections.defaultdict(set)
    for o, n in cross:
        succ[o].add(n)
        pred[n].add(o)
    for x in list(succ):
        if len(succ[x]) < 2:
            continue
        seen = {}
        for first in succ[x]:
            stack, vis = [first], set()
            while stack:
                y = stack.pop()
                if y in vis:
                    continue
                vis.add(y)
                stack.extend(succ.get(y, ()))
            for y in vis:
                seen.setdefault(y, set()).add(first)
        meet = {y for y, v in seen.items() if len(v) >= 2}
        if meet:
            pat["split_then_merged_back"] |= {x} | meet
    pat["merged_then_split"] = {x for x in succ if len(succ[x]) >= 2 and len(pred.get(x, ())) >= 2}

    order = {s: i for i, s in enumerate(tr.order)}
    death = collections.defaultdict(list)
    birth = collections.defaultdict(list)
    for o, ov, n, nv, s, sc in ev:
        if o and n is None:
            death[o].append(s)
        if n and o is None:
            birth[n].append(s)
    for x, ss in death.items():
        first_death = min(ss, key=order.get)
        if any(r >= tr.eff[first_death] for r in tr.present.get(x, ())) or any(order[b] > order[first_death] for b in birth.get(x, ())):
            pat["retired_and_back"].add(x)

    log("gtf features")
    eras = tr.eras()
    era_of = {}
    for asm, lo, hi in eras:
        for r in rels:
            if lo <= r <= hi:
                era_of[r] = asm
    nver = {}
    for x, rr in tr.present.items():
        rr = sorted(rr)
        if ridx[rr[-1]] - ridx[rr[0]] + 1 > len(rr):
            pat["gap"].add(x)
        vs = [(r, tr.genes[r][x]["gene_version"]) for r in rr]
        vv = [(r, int(v)) for r, v in vs if v]
        nver[x] = len({v for _, v in vv})
        for (r1, v1), (r2, v2) in zip(vv, vv[1:]):
            if v2 - v1 > 1:
                pat["version_jump"].add(x)
            if v2 < v1:
                pat["version_fall"].add(x)
        placed = [r for r in rr if tr.genes[r][x]["region"] != "core_only" and not ALT_SEQ.match(tr.genes[r][x]["seqname"])]
        per_era = collections.defaultdict(set)
        chroms = set()
        for r in placed:
            g = tr.genes[r][x]
            per_era[era_of[r]].add((g["seqname"], g["strand"]))
            if CHROM.match(g["seqname"]):
                chroms.add(g["seqname"])
        if any(len(v) > 1 for v in per_era.values()) or len(chroms) > 1:
            pat["moved"].add(x)
    for o, ov, n, nv, s, sc in ev:
        if o and n and o == n and ov and nv and ov.isdigit() and nv.isdigit():
            if int(nv) - int(ov) > 1:
                pat["version_jump"].add(o)
            if int(nv) < int(ov):
                pat["version_fall"].add(o)

    # handoffs
    handoff_sessions = set()
    for s in tr.order:
        b0, a0 = tr.sides(s)
        if b0 is not None and a0 is not None and era_of[b0] != era_of[a0]:
            handoff_sessions.add(s)
    touched = set()
    for o, ov, n, nv, s, sc in ev:
        if s in handoff_sessions and o and n and o != n:
            touched.add(comp_of[o])
    pat["handoff_lineage"] = {x for x, c in comp_of.items() if c in touched}
    log("handoff sessions", sorted(handoff_sessions), "lineages", len(touched))
    mapper = assembly_map.Mapper(org)
    handoff_checked = {}
    for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
        if mapper.path(a1, a2) is None:
            handoff_checked[f"{hi1}->{lo2}"] = "no chain"
            continue
        n_both = n_moved = n_unmapped = 0
        for x, rr in tr.present.items():
            if hi1 in rr and lo2 in rr:
                n_both += 1
                m = mapper.map(tr.locus(x, hi1), a2)
                if "unmapped" in m:
                    n_unmapped += 1
                elif not assembly_map.overlaps(m, tr.genes[lo2][x]):
                    pat["handoff_moved"].add(x)
                    n_moved += 1
        handoff_checked[f"{hi1}->{lo2}"] = dict(living_across=n_both, moved=n_moved, unmapped=n_unmapped)
    log("handoff moves", handoff_checked)

    n_chains = 0
    for x in succ:
        if not pred.get(x) and len(succ[x]) == 1:
            members, y = {x}, next(iter(succ[x]))
            while True:
                members.add(y)
                if len(succ.get(y, ())) != 1 or len(pred.get(y, ())) != 1:
                    break
                y = next(iter(succ[y]))
                if y in members:
                    break
            if len(members) >= 3:
                pat["rename_chain_ge3"] |= members
                n_chains += 1

    log("overlaps and contradictions")
    ovl = {}
    n_contra = collections.Counter()
    mislabelled = {s for s in sessions if H.rel_int(sessions[s]["new_release"]) <= H.rel_int(sessions[s]["old_release"])}
    for o, ov, n, nv, s, sc in ev:
        if s not in tr.sessions:
            continue
        b0, a0 = tr.sides(s)
        if o and n and o != n and a0 is not None and n in tr.genes[a0] and tr.genes[a0][n]["region"] != "core_only":
            if a0 not in ovl:
                ovl[a0] = Overlaps(tr.genes[a0])
            g = tr.genes[a0][n]
            if ovl[a0].others(g["seqname"], g["strand"], int(g["start"]), int(g["end"]), n):
                pat["successor_overlaps"] |= {o, n}
        # contradictions with the GTF on the two sides of the session (only where both sides are in the window)
        if b0 is not None and a0 is not None:
            if o:
                g = tr.genes[b0].get(o)
                if g is None and o != n:
                    n_contra["old id absent before"] += 1
                    pat["event_contradicts_gtf"] |= {o} | ({n} if n else set())
                elif g is not None and g["gene_version"] and ov and g["gene_version"] != ov:
                    n_contra["old version differs"] += 1
                    pat["event_contradicts_gtf"] |= {o} | ({n} if n else set())
            if n:
                g = tr.genes[a0].get(n)
                if g is None and o != n:
                    n_contra["new id absent after"] += 1
                    pat["event_contradicts_gtf"] |= {n} | ({o} if o else set())
                elif g is not None and g["gene_version"] and nv and g["gene_version"] != nv:
                    n_contra["new version differs"] += 1
                    pat["event_contradicts_gtf"] |= {n} | ({o} if o else set())
        if s in mislabelled and o and n and o != n:
            pat["mislabelled_session"] |= {o, n}

    pat["gtf_lacks"] = set().union(*tr.core_only.values()) if tr.core_only else set()
    patch = {x for r in rels for x, g in tr.genes[r].items() if g["region"] != "primary"}
    bad = {comp_of[x] for x in patch if x in comp_of}
    pat["non_reference_lineage"] = {x for x, c in comp_of.items() if c in bad}

    # the ids above the 99th percentile of distinct versions (strictly: at most 1 %; none where no id stands out,
    # e.g. cynomolgus, whose ids have one or two versions)
    ranked = sorted(nver.values(), reverse=True)
    cut = ranked[max(0, round(0.01 * len(ranked)) - 1)]
    pat["many_versions"] = {x for x, v in nver.items() if v > cut}

    hunt, hunt_time = hunt_ids(org, prefix, name)
    for o, ov, n, nv, s, sc in ev:
        if o and n and o != n and s in tr.sessions and tr.lives(o, s):
            pat["branch_while_living"] |= {o, n}
    pat["hunt_history"] = hunt
    exc = assembly_map.Exceptions(org)
    for r in rels:
        rows = exc.rows(r)
        if not rows:
            continue
        for x, g in tr.genes[r].items():
            if g["region"] == "core_only" or g["seqname"] not in rows:
                continue
            kinds = {e["type"] for e in rows[g["seqname"]]}
            if kinds - {"PAR"}:
                pat["_on_exception"].add(x)
            elif exc.project(dict(seqname=g["seqname"], start=g["start"], end=g["end"], strand=g["strand"]), r) \
                    .get("exception") == "PAR":
                pat["par_y"].add(x)
    for r in rels:
        for x, g in tr.genes[r].items():
            if g["region"] == "patch":
                pat["_on_exception"].add(x)
    for x in pat.pop("_on_exception", set()):
        rr = tr.present.get(x)
        if rr and T not in rr and not tr.walk_fwd(x, max(rr), T, True, living=True):
            pat["non_reference_copy_ended"].add(x)
    pat["travel_weak_cell"] = weak_cell_ids(org, name)

    union = set().union(*pat.values())
    in_window = {x for x in union if x in tr.present}
    W = len(rels)
    lines = sum(len(tr.present[x]) * W for x in in_window)
    rows = []
    for x in sorted(union):
        pp = sorted(p for p in PATTERNS if x in pat.get(p, ()))
        rr = sorted(tr.present.get(x, ()))
        rows.append([x, ",".join(pp), comp_of.get(x, x), size_of.get(x, 1), len(sess_of.get(x, ())), nver.get(x, 0),
                     len(rr), rr[0] if rr else "", rr[-1] if rr else ""])
    os.makedirs(os.path.join(_data(), "selection"), exist_ok=True)
    H.write_tsv_gz(os.path.join(_data(), "selection", f"{name}.tsv.gz"),
                   ["id", "patterns", "lineage", "lineage_size", "event_sessions", "versions", "releases", "first",
                    "last"], rows)
    res = dict(organism=org, target=T, window=[rels[0], rels[-1]], window_releases=W, core=tr.core,
               patterns={p: dict(ids=len(pat.get(p, ())), in_window=len(pat.get(p, set()) & set(tr.present)),
                                 how=PATTERNS[p]) for p in PATTERNS},
               only_this_pattern={p: len({x for x in pat.get(p, ()) if x in in_window and
                                          all(x not in pat.get(q, ()) for q in PATTERNS if q != p)}) for p in PATTERNS},
               union_ids=len(union), union_in_window=len(in_window), full_matrix_lines=lines,
               rename_chains_ge3=n_chains, many_versions_cut=cut, contradictions=dict(n_contra),
               mislabelled_sessions=sorted(mislabelled), handoff_sessions=sorted(handoff_sessions),
               handoff_moves=handoff_checked, lineages=dict(all=len(comp_size), ge5=sum(1 for v in comp_size.values() if v >= 5),
                                                             largest=sorted(comp_size.values(), reverse=True)[:5]),
               hunt_findings=dict(files=hunt_time, ids=len(hunt)))
    json.dump(res, open(os.path.join(_data(), "selection", f"{name}.json"), "w"), indent=1)
    log(json.dumps({k: v for k, v in res.items() if k != "patterns"}))
    log(json.dumps({p: v["in_window"] for p, v in res["patterns"].items()}))


if __name__ == "__main__":
    main()
