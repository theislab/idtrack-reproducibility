"""The comparison sets of family G: seeded samples of the published sets of A1, A*, B, C and D1, frozen as G's own
published sets `g_<family>_v1` (docs/families/g.md, "The comparison sets").

Usage: python -m idtrack_reproducibility.families.g.sets --organism <o> [--families a1,astar,b,c,d1] [--draft <folder>]

Per organism and family, the lines are drawn from the latest line of every qid of the source family's sets (the base
set and every revision set of it): the cells below, equal allocation over the cells (a cell smaller than its share is
taken whole and the rest is shared by the others), one seeded permutation per cell (`random.Random(f"{SEED}:<organism>:
<family>:<cell>")` over the sorted qids); the first `N_WEB` lines of the allocation are the web services' lines (`tier
= web`), the first `N_LOCAL` (the same lines and more) the local tools' (`tier = local`).

- `a1` (family A, set `a1_v1`): cells = assembly era of the source release x distance in releases to 115 (0, 1-5, 6-15,
  16-30, 31+) x number of expected genes (0, 1, 2+).
- `astar` (family A*, set `v1`): every hard id (A*'s selection; the line's `id`, unversioned) from 5 origins (the first, the last and the quartiles of the
  releases it is asked from, to 115); the web lines are those of the first `N_WEB / 5` ids of a seeded permutation of
  the ids, the local lines every id's.
- `b` (family B: `b1_v1` with B2's strata, `b4_v1`, `b5_v1`, `b7_v1`): cells = the strata of B1, B2, B4, B5, B7.
- `c` (family C: `c1_v1`, `c2_v1`, `c3_v1` where published): cells = the strata.
- `d1` (family D, `d1_v1`): cells = the strata.

A G line is the source line with: `family = g`, `stratum = g_<family>_<source stratum>`, `source_family`,
`source_set`, `source_stratum`, `source_cmp`, `tier`, `task` (the comparison's task: `a` an old Ensembl id to its gene
at the target, `b` an Ensembl id at its own release, `c` a symbol or name, `d` a previous symbol, alias or synonym, `e`
an external identifier), and the genes compared without their versions (`cmp = base`; `expect`, `broad`,
`open_within`, `apart` and `support`'s keys without versions): most tools answer genes without versions, so every
subject (idtrack too) is judged at the gene. Nothing of the truth is computed here: the lines are the families'.
"""
import argparse
import gzip
import json
import os
import random
import zlib

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

SEED = 20261124
N_WEB = 10_000
N_LOCAL = 100_000
N_ORIGINS = 5
TARGET = 115
DIST_BINS = ((0, 0), (1, 5), (6, 15), (16, 30), (31, 10_000))
FAMILIES = {
    # G family: (source family, base sets, the source strata it takes: a prefix test, or None for every stratum)
    "a1": ("a", ["a1_v1"], None),
    "astar": ("astar", ["v1"], None),
    "b": ("b", ["b1_v1", "b4_v1", "b5_v1", "b7_v1"], ("b1_", "b2_", "b4_", "b5_", "b7_")),
    "c": ("c", ["c1_v1", "c2_v1", "c3_v1"], None),
    "d1": ("d", ["d1_v1"], None),
}


def decided(family, organism, name):
    """A revision set of the judge following decisions not built yet (`<set>_dec`, `<set>_r<n>_dec`): another judge,
    which the families score apart (`scores_decided`); G's lines are the judge users get (`scores_now`)."""
    return name.endswith("_dec")


def revision_closure(family, organism, bases):
    """The base sets that exist and every set that revises one of them (transitively), in publication order; the
    decided judge's revision sets left out (`decided`)."""
    names = [n for n in L.set_names(family, organism) if not decided(family, organism, n)]
    keep = [b for b in bases if b in names]
    changed = True
    while changed:
        changed = False
        for n in names:
            if n not in keep and set(L.manifest(family, organism, n).get("revises") or []) & set(keep):
                keep.append(n)
                changed = True
    return [n for n in names if n in keep]


def scan(family, organism, sets, columns, keep=None, value=None):
    """{qid: value} over the sets' files, the latest line of a qid winning: the value is the row of the given columns,
    or `value(row)`; `keep` (a function of the row) filters while reading (a set of tens of millions of lines holds only
    what is kept)."""
    import pyarrow.parquet as pq
    d = L.query_dir(family, organism)
    by = {}
    for s in sets:
        p = L._set_file(d, s)
        if p.endswith(".parquet"):
            pf = pq.ParquetFile(p)
            have = [c for c in columns if c in pf.schema_arrow.names]
            for b in pf.iter_batches(batch_size=500_000, columns=have):
                cols = b.to_pydict()
                for i, q in enumerate(cols["qid"]):
                    r = {c: cols[c][i] if c in cols else None for c in columns}
                    if keep is None or keep(r):
                        by[q] = r if value is None else value(r)
        else:
            with gzip.open(p, "rt") as fh:
                for raw in fh:
                    rec = json.loads(raw)
                    r = {c: rec.get(c) for c in columns}
                    if keep is None or keep(r):
                        by[rec["qid"]] = r if value is None else value(r)
    return by


def base_qids(family, organism, sets, bases):
    """The qids the base sets publish (a revision set of several sets also holds lines of the others)."""
    d = L.query_dir(family, organism)
    out = set()
    for b in bases:
        if b in sets:
            out.update(L._qid_column(L._set_file(d, b)).to_pylist())
    return out


def allocate(sizes, n):
    """Equal allocation of n over cells of the given sizes (water-filling: a cell smaller than its share is taken
    whole and the rest is shared by the others; a remainder smaller than the cells goes one each to the first cells in
    sorted order); {cell: k}."""
    out = {c: 0 for c in sizes}
    left = min(n, sum(sizes.values()))
    open_cells = sorted(c for c in sizes if sizes[c] > 0)
    while left > 0 and open_cells:
        share = left // len(open_cells)
        if share == 0:
            for c in open_cells[:left]:
                out[c] += 1
            break
        nxt = []
        for c in open_cells:
            t = min(share, sizes[c] - out[c])
            out[c] += t
            left -= t
            if out[c] < sizes[c]:
                nxt.append(c)
        open_cells = nxt
    return out


def permuted(organism, gfam, cell, qids):
    q = sorted(qids)
    random.Random(f"{SEED}:{organism}:{gfam}:{cell}").shuffle(q)
    return q


def dist_bin(d):
    for lo, hi in DIST_BINS:
        if lo <= d <= hi:
            return f"{lo}-{hi}" if hi < 10_000 else f"{lo}+"
    return "?"


def cell_function(gfam, organism):
    """The cell of a candidate row (qid, kind, stratum, from_release, expect)."""
    if gfam == "a1":
        from idtrack_reproducibility.harness import harness_lib as H
        era = H.gtf_assemblies(organism)

        def cell(r):
            fr, ne = r["from_release"], len(r["expect"] or [])
            return f"{era.get(fr, '?')}|{dist_bin(TARGET - fr)}|{min(ne, 2)}{'+' if ne >= 2 else ''}"
        return cell
    return lambda r: r["stratum"] or r["kind"]


def choose(gfam, organism):
    """(web qids, local qids, the allocation per cell, the source sets) of one G family."""
    src, bases, prefixes = FAMILIES[gfam]
    sets = revision_closure(src, organism, bases)
    if not sets:
        return [], [], {}, []
    if gfam == "astar":
        rows = scan(src, organism, sets, ["qid", "kind", "id", "from_release", "to_release", "stratum"],
                    keep=lambda r: r["kind"] == "astar_from" and r["to_release"] == TARGET and r["from_release"] is not None)
        base = base_qids(src, organism, sets, bases)
        per = {}
        for q, r in rows.items():
            if q in base:
                per.setdefault(r["id"], []).append((int(r["from_release"]), q))
        chosen = {}
        for u, xs in per.items():
            xs.sort()
            idx = sorted({round(i * (len(xs) - 1) / (N_ORIGINS - 1)) for i in range(N_ORIGINS)})
            chosen[u] = [xs[i][1] for i in idx]
        units = permuted(organism, gfam, "units", chosen)
        web, local, web_ids = [], [], 0
        for u in units:
            if len(web) + len(chosen[u]) <= N_WEB:
                web += chosen[u]
                web_ids += 1
            local += chosen[u]
        return web, local, {"ids": len(units), "lines": len(local), "web_ids": web_ids}, sets
    keep = (lambda r: (r["stratum"] or r["kind"]).startswith(prefixes)) if prefixes else None
    cell = scan(src, organism, sets, ["qid", "kind", "stratum", "from_release", "expect"], keep=keep,
                value=cell_function(gfam, organism))
    base = base_qids(src, organism, sets, bases)
    cell = {q: c for q, c in cell.items() if q in base}  # a revision set may also revise lines of other sets
    by_cell = {}
    for q, c in cell.items():
        by_cell.setdefault(c, []).append(q)
    sizes = {c: len(v) for c, v in by_cell.items()}
    kw, kl = allocate(sizes, N_WEB), allocate(sizes, N_LOCAL)
    web, local = [], []
    for c in sorted(by_cell):
        perm = permuted(organism, gfam, c, by_cell[c])
        web += perm[: kw.get(c, 0)]
        local += perm[: max(kl.get(c, 0), kw.get(c, 0))]
    return web, local, {c: [sizes[c], kw.get(c, 0), kl.get(c, 0)] for c in sorted(sizes)}, sets


def base_list(xs):
    return None if xs is None else sorted({C.base_id(x) for x in xs})


def task_of(line, src):
    st = L.stratum_of(line)
    if src in ("a", "astar") or (src == "d" and line["kind"] == "d1_id"):
        return "b" if line.get("from_release") == line.get("to_release") else "a"
    if src == "d":
        return "c"
    if src == "b":
        return "d" if st.startswith(("b4_", "b5_")) else "c"
    return "e"


def to_g(line, gfam, src, tier):
    g = dict(line)
    st = L.stratum_of(line)
    g.update(family="g", stratum=f"g_{gfam}_{st.lower()}", source_family=src, source_set=g.pop("set", None),
             source_stratum=st, source_cmp=line.get("cmp"), tier=tier, task=task_of(line, src), cmp="base")
    g["unit"] = L.unit_of(line)
    for k in ("expect", "broad", "open_within", "apart"):
        if g.get(k) is not None:
            g[k] = base_list(g[k])
    if isinstance(g.get("support"), dict):
        sup = {}
        for k, v in g["support"].items():
            sup[C.base_id(k)] = max(sup.get(C.base_id(k), v), v)
        g["support"] = sup
    return g


def draw(gfam, organism):
    """(G lines, manifest extras) of one G family and organism."""
    src = FAMILIES[gfam][0]
    web, local, alloc, sets = choose(gfam, organism)
    if not local:
        return [], {}
    lines = L.load_qids(src, organism, local, sets=sets)
    w = set(web)
    out = [to_g(lines[q], gfam, src, "web" if q in w else "local") for q in sorted(lines, key=lambda q: (zlib.crc32(q.encode()), q))]
    shas = {s: L.manifest(src, organism, s)["sha256"] for s in sets}
    return out, dict(seed=SEED, n_web=len(web), n_local=len(local), allocation=alloc, source_family=src, source_sets=shas)


def revise_judge(organism, dry=False):
    """The revision set `g_<family>_r1` of the lines whose source line was a decided judge's (published before
    `revision_closure` left those sets out): each such qid's line of the judge users get, as a G line of the same tier;
    returns {set: lines}."""
    out = {}
    for gfam, (src, bases, prefixes) in FAMILIES.items():
        name = f"g_{gfam}_v1"
        if name not in L.set_names("g", organism):
            continue
        bad = {}
        for q in L.iter_set(L._set_file(L.query_dir("g", organism), name)):
            if decided(src, organism, q.get("source_set") or ""):
                bad[q["qid"]] = q["tier"]
        if not bad:
            continue
        now = L.load_qids(src, organism, bad, sets=revision_closure(src, organism, bases))
        lines = [to_g(now[q], gfam, src, bad[q]) for q in sorted(bad)]
        rname = f"g_{gfam}_r1"
        out[rname] = len(lines)
        if dry or rname in L.set_names("g", organism):
            continue
        L.publish("g", organism, rname, lines,
                  f"family G's lines of {name} whose source line was a decided judge's revision set (`*_dec`): the "
                  f"line of the judge users get (the families' `scores_now`)", generator=__file__, revises=[name],
                  sources=[f"{src}/{s} sha256 {L.manifest(src, organism, s)['sha256']}"
                           for s in revision_closure(src, organism, bases)])
    return out


TRUTH_FIELDS = ("expect", "broad", "rule", "open", "open_within", "apart", "support", "spellings", "excluded",
                "expect_error", "or_error", "names_read_at", "locus", "locus_mapped")


def revise_from_sources(organism, dry=False):
    """The revision set `g_<family>_r<n>` of the lines whose source line has been revised since G drew it (a
    revision set of the source family published after G's sets): each such qid's current line of the source family
    (the judge users get), as a G line of the same tier; {set: lines}."""
    out = {}
    for gfam, (src, bases, prefixes) in FAMILIES.items():
        names = L.set_names("g", organism)
        if f"g_{gfam}_v1" not in names:
            continue
        mine = {}
        for n in [x for x in names if x == f"g_{gfam}_v1" or x.startswith(f"g_{gfam}_r")]:
            for q in L.iter_set(L._set_file(L.query_dir("g", organism), n)):
                mine[q["qid"]] = q
        now = L.load_qids(src, organism, mine, sets=revision_closure(src, organism, bases))
        changed = []
        for qid, g in mine.items():
            c = now.get(qid)
            if c is None:
                continue
            new = to_g(c, gfam, src, g["tier"])
            if any(new.get(k) != g.get(k) for k in TRUTH_FIELDS):
                changed.append(new)
        if not changed:
            continue
        k = 1 + max([int(x.rsplit("_r", 1)[1]) for x in names if x.startswith(f"g_{gfam}_r")] or [0])
        rname = f"g_{gfam}_r{k}"
        out[rname] = len(changed)
        if dry:
            continue
        L.publish("g", organism, rname, sorted(changed, key=lambda q: q["qid"]),
                  f"family G's lines whose source line the source family revised after G's sets were drawn: the "
                  f"source's current line (the judge users get)", generator=__file__,
                  revises=[x for x in names if x == f"g_{gfam}_v1" or x.startswith(f"g_{gfam}_r")],
                  sources=[f"{src}/{s} sha256 {L.manifest(src, organism, s)['sha256']}"
                           for s in revision_closure(src, organism, bases)])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.sets", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", required=True)
    ap.add_argument("--families", default=",".join(FAMILIES))
    ap.add_argument("--draft", default=None, help="write the drawn lines to this folder instead of publishing them")
    ap.add_argument("--revise-judge", action="store_true", help="publish g_<family>_r1 (revise_judge) and stop")
    ap.add_argument("--revise-sources", action="store_true", help="publish g_<family>_r<n> (revise_from_sources)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    if a.revise_judge:
        print(a.organism, revise_judge(a.organism, a.dry_run))
        return
    if a.revise_sources:
        print(a.organism, revise_from_sources(a.organism, a.dry_run))
        return
    for gfam in a.families.split(","):
        if not a.draft and f"g_{gfam}_v1" in L.set_names("g", a.organism):
            print(a.organism, f"g_{gfam}_v1", "published already")
            continue
        lines, extra = draw(gfam, a.organism)
        if not lines:
            print(f"{a.organism} {gfam}: no source set")
            continue
        name = f"g_{gfam}_v1"
        if a.draft:
            os.makedirs(a.draft, exist_ok=True)
            with gzip.open(os.path.join(a.draft, f"{a.organism}.{name}.jsonl.gz"), "wt") as fh:
                for q in lines:
                    fh.write(C.dumps(q) + "\n")
            json.dump(extra, open(os.path.join(a.draft, f"{a.organism}.{name}.extra.json"), "w"), indent=1, default=str)
            print(a.organism, name, len(lines), "web", extra["n_web"])
            continue
        desc = (f"family G's comparison lines from {extra['source_family']}'s sets {', '.join(extra['source_sets'])}: "
                f"{extra['n_web']} lines asked of the web services (tier web) and {extra['n_local']} of the local tools")
        L.publish("g", a.organism, name, lines, desc, generator=__file__,
                  sources=[f"{extra['source_family']}/{s} sha256 {h}" for s, h in extra["source_sets"].items()],
                  manifest_extra=extra)
        print(a.organism, name, len(lines), "published")


if __name__ == "__main__":
    main()
