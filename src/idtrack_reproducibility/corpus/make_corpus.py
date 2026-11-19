"""Draw a regression corpus of one organism: a seeded, stratified set of api calls (docs/corpus.md, "The corpus").

Usage (in the tree's Python; `pipeline.py make` starts it, as a job or locally):
    python -m idtrack_reproducibility.corpus.make_corpus --tree <tree> --organism <organism> --out <file.jsonl>
        [--seed 20261111] [--tables] [--workdir <folder>]

The graph is the configured one of the tree and organism, loaded through the tree's api of the `d4c2dda` generation
only to sample identifiers: the corpus itself names no version of idtrack. With `--tables` the stratum
`gene_table_at_release` is drawn from the raw Ensembl gene table of every assembly and release, read from the
graph's own folder (the local repository the graph was built in, holding `idtrack_tables/`); without it the stratum
is left out and no table is read or downloaded.

Every random choice comes from one `random.Random(seed)` and every candidate list is sorted before sampling, so the
same graph, tables and seed give the same file. `--out` must not exist: a published corpus is never redrawn, since no
reference answer fits another corpus.

Each line is one call: `{"qid", "entry", "stratum", "args"}`; `entry` is one of convert_identifier,
convert_identifier_multiple (with classify_multiple_conversion on its result), infer_identifier_source,
track_convert (the engine's `Track.convert` with non-default parameters), list (the list accessors and
external_database_forms), resolve_organism (needs the network). A summary (`<out>.summary.json`: seed, releases,
databases, identifiers per stratum, calls per entry, candidate counts) is written beside it.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import random
import string

from idtrack_reproducibility.corpus import _common as C

# Budgets of identifiers per stratum (convert_identifier); each identifier gives 2-3 calls.
BUDGET = {
    "gene_graph_at_release": 60,  # per chosen release
    "gene_table_at_release": 60,  # per chosen release and assembly table holding it
    "gene_unversioned": 600,
    "gene_retired": 400,
    "gene_split": 400,
    "gene_merge": 400,
    "gene_latest": 300,
    "gene_internal_node_name": 40,
    "transcript_versioned": 400,
    "transcript_unversioned": 150,
    "translation_versioned": 400,
    "translation_unversioned": 150,
    "external_per_database": 80,  # per database enabled in the graph
    "synonym_stripped": 300,
    "synonym_prefixed": 40,
    "external_ambiguous": 300,
    "external_changed_gene": 300,
    "nonexistent": 200,
    "malformed_from_valid": 450,
}
N_RELEASES = 12


def choose_releases(g, DB, organism: str) -> list[int]:
    """About twelve releases over the graph's window, including both sides of every assembly handoff."""
    rels = sorted(int(r) for r in g.graph["confident_for_release"])
    lo, hi = rels[0], rels[-1]
    must = {lo, hi}
    for asm, entry in DB.assembly_mysqlport_priority[organism].items():
        first, last = entry["ReleaseRange"]
        if lo < first <= hi:
            must.update({first - 1, first})
        if last is not None and lo <= last < hi:
            must.update({last, last + 1})
    must = {r for r in must if lo <= r <= hi}
    chosen = set(must)
    # Fill evenly: repeatedly add the release farthest from the chosen ones.
    while len(chosen) < N_RELEASES:
        best = max((r for r in range(lo, hi + 1) if r not in chosen), key=lambda r: (min(abs(r - c) for c in chosen), -r))
        chosen.add(best)
    return sorted(chosen)


def table_ids(api, organism: str, release: int, snapshot: int) -> dict[int, list[str]]:
    """Gene `stable_id.version` of a release, per assembly that has the release (the stored raw Ensembl tables)."""
    from idtrack._graph_maker import GraphMaker

    maker = GraphMaker(api.get_database_manager(organism_name=organism, snapshot_release=snapshot))
    out = {}
    for asm, releases in sorted(maker._assembly_releases().items()):
        if release in releases:
            table = maker._assembly_manager("gene", asm, release).get_db("ids")
            out[asm] = sorted(f"{s}.{v}" for s, v in table.values)
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="make_corpus", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True, help="a configured tree of the d4c2dda generation")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True, help="the corpus file to write (must not exist)")
    ap.add_argument("--seed", type=int, default=20261111)
    ap.add_argument("--tables", action="store_true", help="draw gene_table_at_release from the graph folder's tables")
    ap.add_argument("--workdir", default=None, help="where the graph is linked without --tables (default: beside --out)")
    args = ap.parse_args(argv)
    if os.path.exists(args.out):
        raise SystemExit(f"{args.out} exists: a corpus is never redrawn (choose another file)")
    os.environ.setdefault("IDTRACK_ENSEMBL_ACCESS", "https")
    t = C.tree(args.tree)
    if C.generation(t["source"]) != "d4c2dda":
        raise SystemExit(f"tree {args.tree}: the corpus is drawn through the api of the d4c2dda generation")
    graph = C.graph_file(args.tree, args.organism)
    organism, snapshot = C.graph_of(graph)
    if args.tables:
        repo = os.path.dirname(os.path.realpath(graph))
    else:
        repo = os.path.dirname(C.link_graph(graph, args.workdir or os.path.join(os.path.dirname(os.path.abspath(args.out)),
                                                                              f"make_{args.organism}")))
    C.use_clone(t["source"])
    from idtrack import API, DB

    rng = random.Random(args.seed)
    api = API(local_repository=repo)
    api.build_graph(organism_name=organism, snapshot_release=snapshot)
    g = api.track.graph
    NT = DB.node_type_str

    def sample(lst, k):
        lst = sorted(set(lst))
        return lst if len(lst) <= k else rng.sample(lst, k)

    by_type = collections.defaultdict(list)
    for n, d in g.nodes(data=True):
        by_type[d[NT]].append(n)
    for k in by_type:
        by_type[k].sort()
    genes = [n for n in by_type["ensembl_gene"] if not (n.endswith(".Void") or n.endswith(".Retired"))]
    internal = [n for n in by_type["ensembl_gene"] if n.endswith(".Void") or n.endswith(".Retired")]
    releases = choose_releases(g, DB, args.organism)
    main_asm = int(g.graph["genome_assembly"])
    dbs_main = sorted(api.list_external_databases_by_assembly()[main_asm])
    dbs_all = sorted(api.list_external_databases())
    dbs_not_main = sorted(set(dbs_all) - set(dbs_main))

    ids: list[tuple[str, object, dict]] = []  # (stratum, identifier, extra)

    # --- Ensembl genes present at each chosen release (graph) and in the raw Ensembl tables.
    for r in releases:
        present = [n for n in genes if g.is_active(n, r)]
        for n in sample(present, BUDGET["gene_graph_at_release"]):
            ids.append(("gene_graph_at_release", n, {"release": r}))
        for asm, tids in (table_ids(api, args.organism, r, snapshot) if args.tables else {}).items():
            for n in sample(tids, BUDGET["gene_table_at_release"]):
                ids.append(("gene_table_at_release", n, {"release": r, "assembly": asm}))

    for n in sample(by_type["base_ensembl_gene"], BUDGET["gene_unversioned"]):
        ids.append(("gene_unversioned", n, {}))

    def gid(n):
        return g.nodes[n].get("ID")

    split, merge, retired = [], [], []
    for n in genes:
        outs = {gid(w) for _, w, d in g.out_edges(n, data=True) if "old_release" in d and g.nodes[w][NT] == "ensembl_gene"}
        ins = {gid(w) for w, _, d in g.in_edges(n, data=True) if "old_release" in d and g.nodes[w][NT] == "ensembl_gene"}
        if len(outs - {None}) >= 2:
            split.append(n)
        if len(ins - {None}) >= 2:
            merge.append(n)
        if any(w.endswith(".Retired") for w in g.successors(n)):
            retired.append(n)
    for n in sample(retired, BUDGET["gene_retired"]):
        ids.append(("gene_retired", n, {}))
    for n in sample(split, BUDGET["gene_split"]):
        ids.append(("gene_split", n, {}))
    for n in sample(merge, BUDGET["gene_merge"]):
        ids.append(("gene_merge", n, {}))
    latest = [n for n in genes if g.is_active(n, g.last_release)]
    for n in sample(latest, BUDGET["gene_latest"]):
        ids.append(("gene_latest", n, {}))
    for n in sample(internal, BUDGET["gene_internal_node_name"]):
        ids.append(("gene_internal_node_name", n, {}))

    for form in ("transcript", "translation"):
        nodes = by_type[f"ensembl_{form}"]
        for n in sample(nodes, BUDGET[f"{form}_versioned"]):
            ids.append((f"{form}_versioned", n, {}))
        bases = sorted({g.nodes[n].get("ID") or n.split(".")[0] for n in nodes})
        for n in sample(bases, BUDGET[f"{form}_unversioned"]):
            ids.append((f"{form}_unversioned", n, {}))

    # --- External identifiers: per database, synonyms, ambiguous ones, ones that changed their gene.
    ext_by_db = collections.defaultdict(set)
    ambiguous, changed = [], []
    prefix = DB.synonym_id_nodes_prefix
    for n in by_type["external"]:
        base_genes = collections.defaultdict(set)
        for _, w, d in g.out_edges(n, data=True):
            for db, per_asm in d.get("connection", {}).items():
                ext_by_db[db].add(n)
                if g.nodes[w][NT] == "ensembl_gene":
                    for _asm, rels in per_asm.items():
                        base_genes[gid(w)].update(rels)
        if len(base_genes) >= 2:
            ambiguous.append(n)
            sets = list(base_genes.values())
            if any(not (a & b) for i, a in enumerate(sets) for b in sets[i + 1 :]):
                changed.append(n)
    for db in sorted(ext_by_db):
        for n in sample(ext_by_db[db], BUDGET["external_per_database"]):
            ids.append(("external_per_database", n, {"database": db}))
    syn = [n for n in by_type["external"] if n.startswith(prefix)]
    for n in sample(syn, BUDGET["synonym_stripped"]):
        ids.append(("synonym_stripped", n[len(prefix) :], {}))
    for n in sample(syn, BUDGET["synonym_prefixed"]):
        ids.append(("synonym_prefixed", n, {}))
    for n in sample(ambiguous, BUDGET["external_ambiguous"]):
        ids.append(("external_ambiguous", n, {}))
    for n in sample(changed, BUDGET["external_changed_gene"]):
        ids.append(("external_changed_gene", n, {}))

    # --- Identifiers that do not exist, and malformed input.
    gene_prefix = genes[0][: len(genes[0].split(".")[0]) - 11]  # e.g. ENSSSCG
    nonexist = [f"{gene_prefix}{rng.randrange(10**11):011d}" for _ in range(80)]
    nonexist += ["".join(rng.choice(string.ascii_uppercase + string.digits) for _ in range(rng.randint(3, 12))) for _ in range(80)]
    nonexist += [f"NM_{rng.randrange(10**9):09d}", f"XP_{rng.randrange(10**9):09d}.1", "P0DTC2-9", "ENSG", "LOC000000"]
    nonexist += [f"{gene_prefix}00000000000", f"{gene_prefix}99999999999.1", "0", "-1", "NaN", "None", "null"]
    for n in sorted(set(nonexist))[: BUDGET["nonexistent"]]:
        if n not in g.nodes:
            ids.append(("nonexistent", n, {}))

    valid_pool = sample(genes, 150) + sample(by_type["external"], 150) + sample(by_type["base_ensembl_gene"], 50)
    mal = ["", " ", "\t", "\n", ".", "_", "-", "..1", "ENS", "a" * 500, "_".join("ABCDEFGHIJKLM"), "-".join("ABCDEFGHIJKLMNOPQ")]
    mal += ["Ä", "ß", "ACTB​", "ＡＣＴＢ"]
    for v in valid_pool:
        k = rng.randrange(10)
        if k == 0:
            mal.append(f"  {v}  ")
        elif k == 1:
            mal.append(f"{v}\n")
        elif k == 2:
            mal.append(v.lower())
        elif k == 3:
            mal.append(v.upper())
        elif k == 4:
            mal.append(v.swapcase())
        elif k == 5:
            mal.append(f"{v.split('.')[0]}.999")
        elif k == 6:
            mal.append(f"{v.split('.')[0]}_{rng.randint(1, 30)}")
        elif k == 7:
            mal.append(f"{v.split('.')[0]}-{rng.randint(1, 30)}")
        elif k == 8:
            mal.append(f"{v}.1.1")
        else:
            mal.append(v.replace(".", "_").replace("-", "_"))
    for n in mal[: BUDGET["malformed_from_valid"]]:
        ids.append(("malformed", n, {}))
    # Non-string input.
    for n in (None, 0, 12345, 1.5, ["ACTB"]):
        ids.append(("malformed_non_string", n, {}))

    # --- Configurations per identifier.
    fd_weights = [(None, 6), ("ensembl_gene", 2), ("base_ensembl_gene", 1)] + [(d, 1) for d in dbs_main]
    fd_pool = [d for d, w in fd_weights for _ in range(w)]

    def cfg_random(extra):
        c = {
            "to_release": None if rng.random() < 0.5 else rng.choice(releases),
            "final_database": rng.choice(fd_pool),
            "strategy": rng.choice(["best", "all"]),
            "explain": rng.random() < 0.1,
        }
        if rng.random() < 0.25:
            c["from_release"] = extra.get("release") or rng.choice(releases)
        if rng.random() < 0.02 and dbs_not_main:
            c["final_database"] = rng.choice(dbs_not_main + ["not_a_database"])
        return c

    lines = []
    q = 0

    def add(entry, stratum, a):
        nonlocal q
        lines.append({"qid": f"{args.organism}-{q:06d}", "entry": entry, "stratum": stratum, "args": a})
        q += 1

    for stratum, ident, extra in ids:
        add("convert_identifier", stratum, {"identifier": ident})
        add("convert_identifier", stratum, {"identifier": ident, **cfg_random(extra)})
        if rng.random() < 0.4:
            c = cfg_random(extra)
            c.update(strategy="all", explain=True)
            add("convert_identifier", stratum, {"identifier": ident, **c})

    # --- Batches: convert_identifier_multiple + classify_multiple_conversion.
    pool = [i for s, i, _ in ids if isinstance(i, str)]
    for _ in range(50):
        batch = rng.sample(pool, 40)
        kw = {"to_release": rng.choice([None, None, rng.choice(releases)]), "final_database": rng.choice(fd_pool), "strategy": rng.choice(["best", "all"])}
        add("convert_identifier_multiple", "batch", {"identifier_list": batch, **kw})

    # --- infer_identifier_source on lists drawn from one release (and some noise).
    for _ in range(150):
        r = rng.choice(releases)
        present = [n for n in genes if g.is_active(n, r)]
        lst = rng.sample(sorted(present), min(60, len(present))) if present else []
        lst += rng.sample(pool, rng.choice([0, 0, 5, 20]))
        for mode in ("complete", "ensembl_release", "assembly", "assembly_ensembl_release"):
            for win in (True, False):
                add("infer_identifier_source", f"release_{r}", {"id_list": lst, "mode": mode, "report_only_winner": win})
    add("infer_identifier_source", "bad_mode", {"id_list": pool[:10], "mode": "nonsense", "report_only_winner": True})
    add("infer_identifier_source", "empty", {"id_list": [], "mode": "complete", "report_only_winner": True})

    # --- Track.convert with non-default parameters.
    for _ in range(3000):
        s, ident, extra = rng.choice(ids)
        if not isinstance(ident, str):
            continue
        a = {
            "from_id": ident,
            "to_release": rng.choice([None, rng.choice(releases)]),
            "final_database": rng.choice(fd_pool),
            "go_external": rng.random() < 0.7,
            "prioritize_to_one_filter": rng.random() < 0.5,
            "return_path": rng.random() < 0.3,
            # an argument of the api the published corpus was drawn with, which later trees do not have (rule
            # `lrg_copies`; run_corpus.py leaves it out); it stays so that the random stream is that of the corpus
            "deprioritize_lrg_genes": rng.random() < 0.7,
            "return_ensembl_alternative": rng.random() < 0.7,
        }
        if rng.random() < 0.25:
            a["from_release"] = extra.get("release") or rng.choice(releases)
        add("track_convert", s, a)

    for name in ("list_external_databases", "list_genome_assemblies", "list_external_databases_by_assembly",
                 "external_database_forms", "list_ensembl_releases"):
        add("list", name, {"method": name})
    for name in ("human", "Homo sapiens", "homo_sapiens", "9606", 9606, "mouse", "mus_musculus", 10090, "pig",
                 "sus_scrofa", "Sus scrofa", 9823, "rhesus", "macaca_mulatta", "Mmul_10", 9544, "macaca_fascicularis",
                 9541, "crab-eating macaque", "HUMAN", " human ", "nonexistent organism", ""):
        add("resolve_organism", "network", {"tentative_organism_name": name})

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        for ln in lines:
            fh.write(json.dumps(ln, sort_keys=True) + "\n")
    summary = {
        "organism": args.organism,
        "seed": args.seed,
        "releases": releases,
        "main_assembly": main_asm,
        "databases_main": dbs_main,
        "databases_other": dbs_not_main,
        "identifiers_by_stratum": dict(collections.Counter(s for s, _, _ in ids)),
        "calls_by_entry": dict(collections.Counter(x["entry"] for x in lines)),
        "candidates": {"split": len(split), "merge": len(merge), "retired": len(retired), "ambiguous": len(ambiguous),
                       "changed_gene": len(changed), "synonyms": len(syn)},
    }
    C.write_json(C.summary_file(args.out), summary)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
