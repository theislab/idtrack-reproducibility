"""Sizes of the populations the benchmark families would draw from, per organism, read from the harness's truth files
(GTFs of every release, the release-115 core history tables, BioMart archives). Nothing of idtrack is read.

Usage: python -m idtrack_reproducibility.design.sizes <organism> [<truth dir> [<out json>]]
"""

import collections
import glob
import gzip
import json
import os
import re
import sys

import pandas as pd

def main():
    """Usage: python -m idtrack_reproducibility.design.sizes <organism> [<truth dir> [<out json>]]
    (defaults: the configuration's truth store, <families>/design/sizes_<organism>.json)."""
    from idtrack_reproducibility import config as CFG
    organism = sys.argv[1]
    truth = sys.argv[2] if len(sys.argv) > 2 else os.path.join(CFG.get().truth, organism)
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(CFG.get().families, "design", "sizes_" + organism + ".json")
    T = 115

    gtf = os.path.join(truth, organism, "gtf")
    releases = sorted(int(re.search(r"genes_(\d+)", f).group(1)) for f in glob.glob(f"{gtf}/genes_*.tsv.gz"))
    releases = [r for r in releases if r <= T]

    present = collections.defaultdict(list)  # id -> releases
    versions = collections.defaultdict(set)  # id -> versions
    place = collections.defaultdict(set)  # id -> (seqname, strand) over releases (primary only)
    per_release = {}
    nonprimary_ids = set()
    tx_ids, tl_ids, tx_pairs, tl_pairs = set(), set(), 0, 0
    for r in releases:
        g = pd.read_csv(f"{gtf}/genes_{r}.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        prim = g["region"] == "primary" if "region" in g else pd.Series(True, index=g.index)
        per_release[r] = {"genes": int(len(g)), "primary": int(prim.sum()), "names": int((g["gene_name"] != "").sum())}
        for gid, ver, seq, strand, reg in zip(g["gene_id"], g["gene_version"], g["seqname"], g["strand"], g.get("region", ["primary"] * len(g))):
            present[gid].append(r)
            if ver:
                versions[gid].add(ver)
            if reg == "primary":
                place[gid].add((seq, strand))
            else:
                nonprimary_ids.add(gid)
        tf = f"{gtf}/tx_{r}.tsv.gz"
        if os.path.exists(tf):
            t = pd.read_csv(tf, sep="\t", dtype=str, keep_default_na=False, usecols=["transcript_id", "protein_id"])
            tx_ids.update(t["transcript_id"])
            tx_pairs += len(t)
            p = t["protein_id"][t["protein_id"] != ""]
            tl_ids.update(p)
            tl_pairs += int(p.nunique())

    window = len(releases)
    idx = {r: i for i, r in enumerate(releases)}
    n_presence = sum(len(v) for v in present.values())
    gaps = sum(1 for v in present.values() if idx[max(v)] - idx[min(v)] + 1 > len(v))
    current = sum(1 for v in present.values() if T in v)
    retired = len(present) - current
    moved = sum(1 for v in place.values() if len({s for s, _ in v}) > 1)
    flipped = sum(1 for v in place.values() if len({(s, st) for s, st in v}) > len({s for s, _ in v}))

    # assemblies by release (from mapping_session) and the history events
    core = glob.glob(f"{truth}/{organism}/core/*_core_115_*.stable_id_event.txt.gz")[0]
    ev = pd.read_csv(core, sep="\t", header=None, dtype=str, keep_default_na=False,
                     names=["old", "old_v", "new", "new_v", "session", "type", "score"])
    ev = ev[ev["type"] == "gene"]
    ms = pd.read_csv(core.replace("stable_id_event", "mapping_session"), sep="\t", header=None, dtype=str, keep_default_na=False,
                     names=["session", "old_db", "new_db", "old_release", "new_release", "old_assembly", "new_assembly", "created"])
    handoffs = ms[(ms["old_assembly"] != ms["new_assembly"]) & (ms["old_assembly"] != "") & (ms["new_assembly"] != "")]
    births = int((ev["old"] == "\\N").sum())
    deaths = int((ev["new"] == "\\N").sum())
    both = ev[(ev["old"] != "\\N") & (ev["new"] != "\\N")]
    selfs = int((both["old"] == both["new"]).sum())
    cross = both[both["old"] != both["new"]]
    pairs = cross[["old", "new"]].drop_duplicates()
    per_session_out = cross.groupby(["session", "old"])["new"].nunique()
    per_session_in = cross.groupby(["session", "new"])["old"].nunique()
    splits = int((per_session_out >= 2).sum())
    merges = int((per_session_in >= 2).sum())
    # lineage complexity: connected components of the cross-id event graph
    parent = {}


    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x


    for a, b in zip(pairs["old"], pairs["new"]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.Counter(find(x) for x in list(parent))
    sizes = collections.Counter(comp.values())
    deg = collections.Counter(list(pairs["old"]) + list(pairs["new"]))
    events_per_id = cross.groupby("old")["session"].nunique().add(cross.groupby("new")["session"].nunique(), fill_value=0)
    version_jumps = 0
    for o, ov, n, nv in zip(both["old"], both["old_v"], both["new"], both["new_v"]):
        if o == n and ov.isdigit() and nv.isdigit() and int(nv) - int(ov) > 1:
            version_jumps += 1

    # names and external ids by BioMart release
    names = {}
    all_strings = set()
    for d in sorted(glob.glob(f"{truth}/{organism}/biomart/*/")):
        r = os.path.basename(d.rstrip("/"))
        entry = {}
        for f, col in (("genes.tsv", "external_gene_name"), ("synonym.tsv", "external_synonym"), ("hgnc.tsv", None), ("mgi.tsv", None),
                       ("uniprot_gn.tsv", None), ("entrez.tsv", None), ("wikigene.tsv", None), ("refseq_mrna.tsv", None),
                       ("refseq_peptide.tsv", None), ("swissprot.tsv", None), ("sptrembl.tsv", None), ("rfam.tsv", None), ("vgnc.tsv", None)):
            p = d + f
            if not os.path.exists(p):
                continue
            try:
                t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            except Exception as e:  # noqa: BLE001
                entry[f] = f"unreadable {type(e).__name__}"
                continue
            c = col if col in t.columns else [x for x in t.columns if x != "ensembl_gene_id"][-1]
            if f == "genes.tsv" and col not in t.columns:
                c = [x for x in t.columns if x in ("external_gene_id", "external_gene_name", "display_label")][0]
            vals = t[c][t[c] != ""]
            entry[f] = {"rows": int(len(vals)), "distinct": int(vals.nunique()), "column": c}
            if f in ("genes.tsv", "synonym.tsv", "hgnc.tsv", "mgi.tsv", "uniprot_gn.tsv", "wikigene.tsv"):
                all_strings.update(vals)
            if f == "genes.tsv":
                per_gene = t[t[c] != ""].groupby(c)["ensembl_gene_id"].nunique()
                entry["names_of_several_genes"] = int((per_gene > 1).sum())
        names[r] = entry

    res = {
        "organism": organism,
        "releases": [releases[0], releases[-1], window],
        "per_release": per_release,
        "gene_ids_ever": len(present),
        "gene_ids_current_115": current,
        "gene_ids_retired_before_115": retired,
        "gene_ids_with_gaps": gaps,
        "gene_ids_nonprimary_somewhere": len(nonprimary_ids),
        "gene_id_versions_ever": sum(max(1, len(v)) for v in versions.values()),
        "id_release_pairs": n_presence,
        "ids_moved_sequence": moved,
        "ids_flipped_strand": flipped,
        "transcript_ids_ever": len(tx_ids),
        "translation_ids_ever": len(tl_ids),
        "transcript_id_release_pairs": tx_pairs,
        "events_gene_rows": int(len(ev)),
        "events_births": births,
        "events_retirements": deaths,
        "events_self": selfs,
        "events_cross_rows": int(len(cross)),
        "events_cross_pairs_distinct": int(len(pairs)),
        "split_events": splits,
        "merge_events": merges,
        "self_version_jumps_gt1": version_jumps,
        "mapping_sessions": int(len(ms)),
        "assembly_handoffs": [[r.old_release, r.new_release, r.old_assembly, r.new_assembly] for r in handoffs.itertuples()],
        "lineage_components": int(len(comp)),
        "lineage_component_sizes_top": sorted(comp.values(), reverse=True)[:15],
        "lineage_components_ge5": int(sum(c for s, c in sizes.items() if s >= 5)),
        "lineage_components_ge10": int(sum(c for s, c in sizes.items() if s >= 10)),
        "ids_degree_ge3": int(sum(1 for v in deg.values() if v >= 3)),
        "ids_in_ge3_sessions": int((events_per_id >= 3).sum()),
        "names": names,
        "name_strings_all_biomart_releases": len(all_strings),
        "estimates": {
            "to_latest": n_presence,
            "to_neighbours_pm1_pm2": 4 * n_presence,
            "full_matrix": sum(len(v) for v in present.values()) * window,
            "event_pairs_both_directions": 2 * int(len(pairs)),
        },
    }
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_release", "names")}, indent=1))


if __name__ == "__main__":
    main()
