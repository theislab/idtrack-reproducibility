"""Size of the hardest-lineage family A* per organism: ids per pattern and the full-matrix lines of their union.
Read from the harness's truth files only (GTFs, release-115 core history). Usage: python -m idtrack_reproducibility.design.hard_sizes <organism> [<truth dir> [<out json>]]
"""

import collections
import glob
import json
import re
import sys

import pandas as pd

def main():
    """Usage: python -m idtrack_reproducibility.design.hard_sizes <organism> [<truth dir> [<out json>]]
    (defaults: the configuration's truth store, <families>/design/hard_<organism>.json)."""
    from idtrack_reproducibility import config as CFG
    organism = sys.argv[1]
    truth = sys.argv[2] if len(sys.argv) > 2 else os.path.join(CFG.get().truth, organism)
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(CFG.get().families, "design", "hard_" + organism + ".json")
    T = 115
    gtf = f"{truth}/{organism}/gtf"
    releases = sorted(int(re.search(r"genes_(\d+)", f).group(1)) for f in glob.glob(f"{gtf}/genes_*.tsv.gz"))
    releases = [r for r in releases if r <= T]
    idx = {r: i for i, r in enumerate(releases)}
    present = collections.defaultdict(list)
    ver = collections.defaultdict(dict)
    place = collections.defaultdict(set)
    for r in releases:
        g = pd.read_csv(f"{gtf}/genes_{r}.tsv.gz", sep="\t", dtype=str, keep_default_na=False, usecols=["gene_id", "gene_version", "seqname", "strand"])
        for gid, v, s, st in zip(g["gene_id"], g["gene_version"], g["seqname"], g["strand"]):
            present[gid].append(r)
            if v:
                ver[gid][r] = int(v)
            place[gid].add((s, st))

    core = glob.glob(f"{truth}/{organism}/core/*_core_115_*.stable_id_event.txt.gz")[0]
    ev = pd.read_csv(core, sep="\t", header=None, dtype=str, keep_default_na=False, names=["old", "old_v", "new", "new_v", "session", "type", "score"])
    ev = ev[ev["type"] == "gene"]
    cross = ev[(ev["old"] != "\\N") & (ev["new"] != "\\N") & (ev["old"] != ev["new"])]
    parent = {}


    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x


    for a, b in zip(cross["old"], cross["new"]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for x in list(parent):
        comp[find(x)].add(x)
    size = {x: len(comp[find(x)]) for x in parent}

    pat = collections.defaultdict(set)
    pat["component_ge5"] = {x for x, s in size.items() if s >= 5}
    pat["component_ge20"] = {x for x, s in size.items() if s >= 20}
    # split then merged back: a descendant reachable by two different first steps
    succ = collections.defaultdict(set)
    for a, b in zip(cross["old"], cross["new"]):
        succ[a].add(b)
    for a in list(succ):
        if len(succ[a]) < 2:
            continue
        seen = {}
        for first in succ[a]:
            stack, vis = [first], set()
            while stack:
                y = stack.pop()
                if y in vis:
                    continue
                vis.add(y)
                stack.extend(succ.get(y, ()))
            for y in vis:
                seen.setdefault(y, set()).add(first)
        if any(len(v) >= 2 for v in seen.values()):
            pat["split_then_merged_back"] |= {a} | {y for y, v in seen.items() if len(v) >= 2}
    # merged then split: an id with >=2 predecessors and >=2 successors
    pred = collections.Counter(cross["new"])
    pat["merged_then_split"] = {x for x in succ if len(succ[x]) >= 2 and pred[x] >= 2}
    # retired and back; gaps
    ms = pd.read_csv(core.replace("stable_id_event", "mapping_session"), sep="\t", header=None, dtype=str, keep_default_na=False,
                     names=["session", "old_db", "new_db", "old_release", "new_release", "old_assembly", "new_assembly", "created"])
    order = {s: (float(re.match(r"[0-9.]+", o).group(0).rstrip(".") or 0), c) for s, o, c in zip(ms["session"], ms["old_release"], ms["created"])}
    death = ev.loc[ev["new"] == "\\N"].groupby("old")["session"].apply(lambda x: min(order[s] for s in x))
    birth = ev.loc[ev["old"] == "\\N"].groupby("new")["session"].apply(lambda x: max(order[s] for s in x))
    pat["retired_and_back"] = {x for x in set(death.index) & set(birth.index) if birth[x] > death[x]}
    pat["gaps"] = {x for x, v in present.items() if idx[max(v)] - idx[min(v)] + 1 > len(v)}
    # version jumps and falls in consecutive GTF releases
    for x, d in ver.items():
        rs = sorted(d)
        for a, b in zip(rs, rs[1:]):
            if d[b] - d[a] > 1:
                pat["version_jump"].add(x)
            if d[b] < d[a]:
                pat["version_fall"].add(x)
    pat["moved_sequence"] = {x for x, v in place.items() if len({s for s, _ in v}) > 1}
    # chains of renames: one successor, one predecessor, length >= 3
    chains = 0
    for a in succ:
        if pred[a] == 0 and len(succ[a]) == 1:
            n, y = 1, next(iter(succ[a]))
            members = {a}
            while True:
                members.add(y)
                n += 1
                if len(succ.get(y, ())) != 1 or pred[y] != 1:
                    break
                y = next(iter(succ[y]))
            if n >= 3:
                pat["rename_chain_ge3"] |= members
                chains += 1

    union = set().union(*pat.values())
    window = len(releases)
    lines = sum(len(present.get(x, [])) * window for x in union)
    res = {
        "organism": organism,
        "window": window,
        "patterns": {k: len(v) for k, v in pat.items()},
        "rename_chains_ge3": chains,
        "union_ids": len(union),
        "union_full_matrix_lines": lines,
        "largest_components": sorted((len(v) for v in comp.values()), reverse=True)[:5],
    }
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
