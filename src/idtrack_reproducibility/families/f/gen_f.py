"""F1's datasets and query lines (docs/families/f.md, "The sets"): synthetic count matrices with a known truth.

Usage: gen_f.py --organism <o> [--set f1_v1] [--dry]

Per organism, five source releases R spread over the window (`choose_releases`: at least one per assembly era, the
rest by the highest average of releases per slot; the target release 115 itself is no source), and per release four
datasets, one per form of feature names: `ids` (Ensembl gene ids without version), `versioned` (the id with the
version it had at R), `symbols` (display names of R) and `mix` (all three, with duplicates): up to 20 datasets of
`N_FEATURES` features and `N_CELLS` cells per organism, named `r<RRR>_<form>`. A release whose GTF writes no versions
(`gtf_writes_versions`: most of its genes have none; pig and human before 77, mouse before 76) has no `versioned`
dataset, and its `mix` holds ids and names only (A1's few versioned lines there, genes of the core gene table that
the GTF lacks, are no form of the release).

Every feature is drawn (seeded) from a published line with its truth, never from idtrack:

- ids and versioned ids: family A's A1 lines from R to 115 (`a1_*` sets, the latest line of each qid), classed by the
  id's history (`<families>/a/patterns/<organism>.tsv.gz`: unchanged, born_late, renamed, merged, split, retired
  (`retired_none`)) and `resurrected` (an A1 line whose broad walk holds ids reached only through a return, rule
  `resurrection`);
- names: family B's B1 lines read at R to 115 (`b1_*` sets): `name_unique` (`b1_to115_unique`), `name_shared` (one
  name of several genes, `b1_to115_shared`), `name_moved` (`b2_moved_from`); a name holding `;` is left out (the
  harmonizer joins a merged feature's inputs with `;` in `Query ID`, so such a name could not be told apart there);
- made-unique copies (`made_unique`, rule `made_unique_copies`): `X.1` or `X-1` beside a unique name X of the same
  dataset, where no string of the truth store folds to the form (every GTF's gene ids and names, every BioMart table
  and VGNC gene set, every release's core cross-references and synonyms; `taken_forms`): the copy's truth is X's line;
- absent (`absent`): the accuracy harness's `nonexistent` lines (ids with the organism's prefix that never existed for
  the id forms, random letters and RefSeq-like accessions for the names), expect nothing;
- duplicates (`duplicate`, mix only): a gene asked twice in one dataset, as its id and as its symbol (an A1 and a B1
  line of R expecting the same one gene), or as its versioned and unversioned id; merge partners: for a merged id,
  another id of R whose line expects the same one gene, of any class (an n-to-1 group, extra `group`);
- the core (`core`, every dataset): `N_CORE` genes unchanged at every chosen release with a unique display name there,
  in each dataset's form, so that the intersection over the datasets is not empty.

The shares per class and form are `QUOTAS` (a class with fewer lines at R takes them all; the rest is filled with
unchanged ids or unique names); the counts drawn are in the manifest (`counts`, `drawn`).

A line of the set (`f1|<dataset>|<feature>`) is its source line's truth (expect, broad, rule, cmp, locus, apart, ...)
with the feature as `id`, `from_release` R, `to_release` 115, `strategy` `best` (the harmonizer reads the best view),
`family` f, `stratum` `f1_<class>`, `kind` `f1_<form>`, `unit` the source line's unit (lineage or name; a copy: X; an
absent string: itself) and the extra fields `dataset`, `form`, `source` (`<family>:<set>:<qid>`), `history_class`,
`group`, `copy_of`, `core`.

The datasets (`<families>/f/datasets/<organism>/<set>/<dataset>.h5ad`, written with h5py in AnnData's format,
`h5ad.py`): `X` counts drawn per feature and cell from a negative binomial (`counts`: mean m_f * s_c with m_f
log-normal (mu 0, sigma `SIGMA_MEAN`) and the cell's size factor s_c log-normal (0, `SIGMA_CELL`), dispersion
`DISPERSION`; float32, CSR), the layer `unspliced` drawn alike with a third of the mean (int32, CSR), `obs` (a
categorical, text, float, integer and boolean column, and one column of some forms only), `obsm` (`X_pca`, `X_umap`),
`obsp` (`connectivities`), `varm` (`PCs`), `var` (`feature_type`, `highly_variable`, `mean_input`) and `uns` (keys
equal in every dataset, keys that differ, a nested dict); seeds `seed_of(organism, dataset)`. The manifest records
each file's sha256, its seed and the count model.
"""
import argparse
import collections
import hashlib
import json
import os
import random
import sys

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

FORMS = ("ids", "versioned", "symbols", "mix")
N_RELEASES = 5
N_FEATURES = 2000
N_CELLS = 200
N_CORE = 200
TARGET = 115
SEED = 20261124
SIGMA_MEAN = 1.5
SIGMA_CELL = 0.3
DISPERSION = 2.0
ID_CLASSES = ("unchanged", "born_late", "renamed", "merged", "split", "retired", "resurrected")
NAME_CLASSES = ("name_unique", "name_shared", "name_moved")
# per form: (source, class) -> features (the core genes count in `unchanged` / `name_unique`); `fill` takes the rest
QUOTAS = {
    "ids": {("id", "unchanged"): 600, ("id", "born_late"): 300, ("id", "renamed"): 150, ("id", "merged"): 250,
            ("id", "split"): 250, ("id", "retired"): 250, ("id", "resurrected"): 100, ("absent", "absent"): 50,
            "partners": 50, "fill": ("id", "unchanged")},
    "versioned": {("ver", "unchanged"): 600, ("ver", "born_late"): 300, ("ver", "renamed"): 150, ("ver", "merged"): 250,
                  ("ver", "split"): 250, ("ver", "retired"): 250, ("ver", "resurrected"): 100,
                  ("absent", "absent"): 50, "partners": 50, "fill": ("ver", "unchanged")},
    "symbols": {("name", "name_unique"): 1150, ("name", "name_shared"): 250, ("name", "name_moved"): 300,
                ("copy", "made_unique"): 150, ("absent", "absent"): 50, "fill": ("name", "name_unique")},
    "mix": {("id", "unchanged"): 250, ("id", "born_late"): 100, ("id", "renamed"): 50, ("id", "merged"): 75,
            ("id", "split"): 50, ("id", "retired"): 50, ("id", "resurrected"): 25, ("ver", "unchanged"): 150,
            ("ver", "merged"): 50, ("ver", "split"): 50, ("ver", "retired"): 50, ("name", "name_unique"): 450,
            ("name", "name_shared"): 100, ("name", "name_moved"): 150, ("copy", "made_unique"): 50,
            ("absent", "absent"): 50, ("dup", "id_and_name"): 100, ("dup", "versioned_and_unversioned"): 50,
            "fill": ("name", "name_unique")},
}
TRUTH_KEYS_DROPPED = {"qid", "family", "stratum", "kind", "unit", "id", "id_codec", "strategy", "set", "from_release",
                      "to_release", "note", "era", "dist", "cross"}


def seed_of(*parts):
    return int(hashlib.sha256(":".join(str(p) for p in (SEED,) + parts).encode()).hexdigest()[:15], 16)


def choose_releases(releases, era_of, k=N_RELEASES):
    """`k` source releases spread over the window (releases before the target), at least one per assembly era: each era
    starts with one slot, every further slot goes to the era with the most releases per slot (D'Hondt), and an era's
    slots are its releases at evenly spaced positions."""
    rels = sorted(r for r in releases if r < TARGET)
    eras = collections.OrderedDict()
    for r in rels:
        eras.setdefault(era_of(r), []).append(r)
    if len(eras) > k:  # more eras than slots: the eras with the most releases
        keep = sorted(sorted(eras, key=lambda e: -len(eras[e]))[:k], key=lambda e: eras[e][0])
        eras = collections.OrderedDict((e, eras[e]) for e in keep)
    slots = {e: 1 for e in eras}
    for _ in range(min(k, len(rels)) - len(slots)):
        e = max(eras, key=lambda x: (len(eras[x]) / (slots[x] + 1), -eras[x][0]))
        slots[e] += 1
    out = []
    for e, rs in eras.items():
        n = min(slots[e], len(rs))
        out += [rs[min(len(rs) - 1, int((i + 0.5) * len(rs) / n))] for i in range(n)]
    return sorted(set(out))


def class_of_a1(line, history_class):
    """The class of an A1 line: `resurrected` where its broad walk holds ids reached only through a return (rule
    `resurrection`), else the id's history class (`retired_none` -> `retired`)."""
    b = line.get("broad")
    if b is not None and set(b) - set(line.get("expect") or []):
        return "resurrected"
    return {"retired_none": "retired"}.get(history_class, history_class)


def f_line(src, family, dataset, form, release, feature, cls, unit, **extra):
    """An F line from its source line: the source's truth, the feature as asked, the harmonizer's view."""
    q = {k: v for k, v in src.items() if k not in TRUTH_KEYS_DROPPED}
    q.update(qid=f"f1|{dataset}|{feature}", family="f", stratum=f"f1_{cls}", kind=f"f1_{form}", unit=unit, id=feature,
             from_release=release, to_release=TARGET, strategy="best", dataset=dataset, form=form,
             source=f"{family}:{src.get('set')}:{src['qid']}")
    q.update({k: v for k, v in extra.items() if v is not None})
    return q


def single(line):
    e = line.get("expect") or []
    return e[0] if len(e) == 1 else None


def draw(pool, n, rng, taken):
    """Up to n items of a pool (a list of (feature, line)) whose feature is not taken yet, in a seeded order."""
    items = [x for x in pool if x[0] not in taken]
    rng.shuffle(items)
    out = items[:n]
    taken.update(x[0] for x in out)
    return out


def dataset_lines(form, release, pools, core, rng, n=N_FEATURES, quotas=None):
    """The lines of one dataset: the core genes in the dataset's form, then each (source, class) quota drawn from the
    pools (`pools[(source, class)]`: lists of (feature, line, family, unit, history_class)), merge partners, made-unique
    copies, duplicates, then the fill. Returns the lines (in a seeded order) and the counts drawn per quota."""
    quotas = quotas or QUOTAS[form]
    dataset = f"r{release:03d}_{form}"
    taken, out, drawn = set(), [], collections.Counter()

    item_of = {}

    def add(item, cls, **extra):
        feature, line, family, unit, hist = item
        taken.add(feature)
        item_of[feature] = item
        out.append(f_line(line, family, dataset, form, release, feature, cls, unit, history_class=hist, **extra))

    # the core: each gene in the dataset's form (mix: by the gene's hash)
    for g, by_form in core:
        if form == "mix":
            keys = [k for k in ("id", "ver", "name") if k in by_form]
            key = keys[seed_of(g) % len(keys)]
        else:
            key = {"ids": "id", "versioned": "ver", "symbols": "name"}[form]
        item = by_form[key]
        if item[0] not in taken:
            add(item, "unchanged" if key != "name" else "name_unique", core=True)
            drawn[("core", key)] += 1
    for key, k in quotas.items():
        if not isinstance(key, tuple) or key[0] in ("copy", "dup"):
            continue
        src, cls = key
        have = drawn[("core", src)] if cls in ("unchanged", "name_unique") else 0
        pool = [(it[0], it) for it in pools.get(key, [])]
        for _, it in draw(pool, max(0, k - have), rng, set(taken)):
            add(it, cls)
            drawn[key] += 1
    # merge partners: another id of R expecting the same one gene as a merged id of the dataset
    if quotas.get("partners"):
        by_gene = pools.get(("by_gene", "ver" if form == "versioned" else "id"), {})
        merged = [x for x in out if x["stratum"] == "f1_merged"]
        rng.shuffle(merged)
        left = quotas["partners"]
        for x in merged:
            if left <= 0:
                break
            g = single(x)
            cand = [it for it in by_gene.get(g, []) if it[0] not in taken]
            if not cand:
                continue
            it = rng.choice(cand)
            group = f"partner:{g}"
            x["group"] = group
            cls = it[4] if it[4] in ID_CLASSES else "unchanged"
            add(it, {"retired_none": "retired"}.get(cls, cls), group=group)
            drawn[("partner", "id")] += 1
            left -= 1
    # made-unique copies of unique names of the dataset
    k = quotas.get(("copy", "made_unique"), 0)
    if k:
        allowed = pools.get(("copy_ok", None), set())
        names = [x for x in out if x["stratum"] == "f1_name_unique" and not x.get("core") and single(x)]
        rng.shuffle(names)
        for x in names:
            if drawn[("copy", "made_unique")] >= k:
                break
            sep = ".1" if drawn[("copy", "made_unique")] % 2 == 0 else "-1"
            form_ = x["id"] + sep
            if form_ in taken or form_ not in allowed:
                continue
            _f, line, fam, unit, hist = item_of[x["id"]]
            add((form_, line, fam, unit, hist), "made_unique", copy_of=x["id"], group=f"copy:{x['id']}")
            x["group"] = f"copy:{x['id']}"
            drawn[("copy", "made_unique")] += 1
    # duplicates: one gene as its id and its symbol; as its versioned and unversioned id
    k = quotas.get(("dup", "id_and_name"), 0)
    pairs = [p for p in pools.get(("pairs", "id_and_name"), []) if p[0][0] not in taken and p[1][0] not in taken]
    rng.shuffle(pairs)
    for a, b in pairs[: k // 2]:
        g = single(a[1])
        for it in (a, b):
            add(it, "duplicate", group=f"dup:{g}")
        drawn[("dup", "id_and_name")] += 2
    k = quotas.get(("dup", "versioned_and_unversioned"), 0)
    pairs = [p for p in pools.get(("pairs", "versioned_and_unversioned"), []) if p[0][0] not in taken and p[1][0] not in taken]
    rng.shuffle(pairs)
    for a, b in pairs[: k // 2]:
        g = single(a[1]) or a[0]
        for it in (a, b):
            add(it, "duplicate", group=f"dup:{g}")
        drawn[("dup", "versioned_and_unversioned")] += 2
    # the fill
    fill = quotas["fill"]
    pool = [(it[0], it) for it in pools.get(fill, [])]
    for _, it in draw(pool, max(0, n - len(out)), rng, set(taken)):
        add(it, fill[1])
        drawn[("fill",) + fill] += 1
    if len(out) > n:
        raise SystemExit(f"{dataset}: {len(out)} features drawn for {n}: the quotas exceed the dataset")
    rng.shuffle(out)
    return out, drawn


# ---------------- the counts and the metadata ----------------

def counts(seed, n_cells, n_features):
    """(X float32 CSR, unspliced int32 CSR, the per-feature means): negative binomial counts, mean m_f * s_c, dispersion
    `DISPERSION` (numpy's negative_binomial(r, r / (r + mean)))."""
    import numpy as np
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    m = rng.lognormal(0.0, SIGMA_MEAN, size=n_features)
    s = rng.lognormal(0.0, SIGMA_CELL, size=n_cells)
    mean = np.outer(s, m)
    r = DISPERSION
    x = rng.negative_binomial(r, r / (r + mean))
    u = rng.negative_binomial(r, r / (r + mean / 3.0))
    return sp.csr_matrix(x.astype(np.float32)), sp.csr_matrix(u.astype(np.int32)), m


def metadata(organism, dataset, form, release, seed, features, X):
    """obs, var, obsm, obsp, varm and uns of a dataset (module docstring), seeded."""
    import numpy as np
    import scipy.sparse as sp
    rng = np.random.default_rng(seed + 1)
    n = X.shape[0]
    types = ["acinar", "ductal", "alpha", "beta", "delta"]
    obs = {"index": [f"c{i:03d}" for i in range(n)],
           "columns": {"cell_type": {"categories": types, "codes": rng.integers(0, len(types), n).tolist()},
                       "donor": [f"{dataset}_d{i % 4}" for i in range(n)],
                       "n_counts": np.asarray(X.sum(axis=1)).ravel().astype(np.float64),
                       "batch_int": (np.arange(n) % 3).astype(np.int64),
                       "is_doublet": rng.random(n) < 0.05}}
    if form in ("symbols", "mix"):
        obs["columns"]["stage"] = [f"E{10 + i % 5}" for i in range(n)]
    if form == "ids":
        obs["columns"]["plate"] = (np.arange(n) // 50).astype(np.int64)
    var = {"index": list(features),
           "columns": {"feature_type": ["Gene Expression"] * len(features),
                       "highly_variable": rng.random(len(features)) < 0.2,
                       "mean_input": np.asarray(X.mean(axis=0)).ravel().astype(np.float64)}}
    k = 5
    rows = np.repeat(np.arange(n), k)
    cols = rng.integers(0, n, n * k)
    conn = sp.csr_matrix((rng.random(n * k).astype(np.float32), (rows, cols)), shape=(n, n))
    conn.sum_duplicates()
    obsm = {"X_pca": rng.normal(size=(n, 10)).astype(np.float32), "X_umap": rng.normal(size=(n, 2)).astype(np.float32)}
    varm = {"PCs": rng.normal(size=(len(features), 5))}
    uns = {"organism": organism, "f1_set": "f1_v1", "dataset": dataset, "source_release": int(release),
           "params": {"seed": int(seed % 1_000_000), "form": form}}
    return obs, var, obsm, {"connectivities": conn}, varm, uns


# ---------------- reading the sources ----------------

def a1_sets(org):
    return [s for s in L.set_names("a", org) if s.startswith("a1_")]


def b1_sets(org):
    return [s for s in L.set_names("b", org) if s.startswith("b1_")]


def read_parquet_rows(path, columns, filt):
    import pyarrow.parquet as pq
    return pq.read_table(path, columns=columns, filters=filt).to_pylist()


def patterns(org):
    import gzip
    out = {}
    with gzip.open(os.path.join(C.family_dir("a"), "patterns", f"{org}.tsv.gz"), "rt") as fh:
        next(fh)
        for x in fh:
            f = x.rstrip("\n").split("\t")
            out[f[0]] = f[2]
    return out


def gtf_writes_versions(org, release):
    """Whether the GTF of a release writes versions: most of its genes have a `gene_version`."""
    from idtrack_reproducibility.harness import harness_lib as H
    genes = H.gtf_genes(org, release) or {}
    return bool(genes) and sum(1 for r in genes.values() if r.get("gene_version")) > len(genes) / 2


def fold_index(forms):
    """{folded string: the forms that fold to it} (`X-1` and `X.1` fold alike: a string of the store takes both)."""
    from idtrack_reproducibility.families.e import truth_e as TE
    want = collections.defaultdict(set)
    for f in forms:
        want[TE.fold(f)].add(f)
    return want


def taken_forms(org, forms, log=print):
    """The candidate forms (made-unique copies) that a string of the truth store folds to (case and `-`, `_`, `.`
    taken as one): every GTF's gene ids and names, every BioMart table and VGNC gene set, every release's core
    cross-references (accessions and labels) and synonyms. A copy is asked only where its form is no such string."""
    import gzip
    from idtrack_reproducibility.families.c import truth_c as TC
    from idtrack_reproducibility.families.e import truth_e as TE
    from idtrack_reproducibility.harness import harness_lib as H
    want = fold_index(forms)
    hit = set()

    def see(s):
        if s:
            k = TE.fold(s)
            if k in want:
                hit.update(want[k])

    for rel in H.gtf_releases(org):
        for r in H.read_tsv(os.path.join(C.TRUTH, org, "gtf", f"genes_{rel}.tsv.gz")):
            see(r["gene_id"])
            see(r.get("gene_name"))
    for rel in H.biomart_releases(org):
        d = os.path.join(C.TRUTH, org, "biomart", str(rel))
        for f in sorted(os.listdir(d)):
            if f.endswith(".tsv"):
                t = H.biomart(org, rel, f[:-4])
                for row in (t[1] if t else []):
                    for s in row:
                        see(s)
        vg = H.vgnc_symbols(org, rel)
        for s in (vg[1] if vg else {}):
            see(s)
    xdir = os.path.join(C.TRUTH, org, "core_xref")
    for f in sorted(os.listdir(xdir)) if os.path.isdir(xdir) else []:
        if f.endswith(".tsv.gz"):
            for _ot, _oid, _g, _db, _disp, acc, lab, _v, _it in TC.read_xref(org, int(f.split(".")[0])):
                see(acc)
                see(lab)
    sdir = os.path.join(C.TRUTH, org, "core_synonyms")
    for f in sorted(os.listdir(sdir)) if os.path.isdir(sdir) else []:
        if f.endswith(".tsv.gz"):
            with gzip.open(os.path.join(sdir, f), "rt") as fh:
                next(fh)
                for x in fh:
                    see(x.split("\t", 1)[0])
    log(org, "made-unique candidates", len(forms), "taken by the truth store", len(hit))
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", default="f1_v1")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    org = a.organism
    log = lambda *x: print(*x, file=sys.stderr, flush=True)  # noqa: E731
    from idtrack_reproducibility.harness import harness_lib as H
    from idtrack_reproducibility.harness import score
    places = score.Places(org)
    qd_a = L.query_dir("a", org)
    rels_all = sorted({r["from_release"] for r in read_parquet_rows(os.path.join(qd_a, "a1_v1.parquet"), ["from_release"], None)})
    releases = choose_releases(rels_all, places.assembly)
    log(org, "releases", releases, {r: places.assembly(r) for r in releases})
    pat = patterns(org)
    # A1 lines of the chosen releases (published set: which qids exist; the latest line read below)
    a1 = read_parquet_rows(os.path.join(qd_a, "a1_v1.parquet"), ["qid", "id", "from_release", "to_release"],
                           [("from_release", "in", releases), ("to_release", "=", TARGET)])
    qb = L.query_dir("b", org)
    b1 = read_parquet_rows(os.path.join(qb, "b1_v1.parquet"), ["qid", "id", "from_release", "stratum"],
                           [("from_release", "in", releases), ("to_release", "=", TARGET),
                            ("stratum", "in", ["b1_to115_unique", "b1_to115_shared", "b2_moved_from"])])
    log(org, "A1 lines", len(a1), "B lines", len(b1))
    la = L.load_qids("a", org, {r["qid"] for r in a1}, sets=a1_sets(org))
    lb = L.load_qids("b", org, {r["qid"] for r in b1}, sets=b1_sets(org))
    absent = [q for q in L.load("harness", org, kinds=["nonexistent"])[0]]
    log(org, "latest lines A", len(la), "B", len(lb), "absent", len(absent))
    pools = {r: collections.defaultdict(list) for r in releases}
    by_gene = {r: {"id": collections.defaultdict(list), "ver": collections.defaultdict(list)} for r in releases}
    name_of = {r: {} for r in releases}  # gene -> unique-name item at R
    ids_of = {r: {} for r in releases}  # base id -> {"id": item, "ver": item}
    for q in la.values():
        if q.get("excluded") or q.get("open"):
            continue
        R = q["from_release"]
        base = C.base_id(q["id"])
        hist = pat.get(base, "unchanged")
        cls = class_of_a1(q, hist)
        src = "ver" if q["id"] != base else "id"
        item = (q["id"], q, "a", q.get("unit") or base, hist)
        pools[R][(src, cls)].append(item)
        ids_of[R].setdefault(base, {})[src] = item
        if single(q):
            by_gene[R][src][single(q)].append(item)
    cls_b = {"b1_to115_unique": "name_unique", "b1_to115_shared": "name_shared", "b2_moved_from": "name_moved"}
    for q in lb.values():
        if q.get("excluded") or q.get("open"):
            continue
        R = q["from_release"]
        cls = cls_b.get(q.get("stratum"))
        if not cls or ";" in str(q["id"]):
            continue  # (a name holding `;` cannot be told apart in the harmonizer's `;`-joined `Query ID`)
        g = single(q)
        item = (str(q["id"]), q, "b", q.get("unit") or str(q["id"]), None)
        pools[R][("name", cls)].append(item)
        if cls == "name_unique" and g:
            name_of[R].setdefault(g, item)
    for q in absent:
        note = q.get("note") or ""
        item = (q["id"], q, "harness", q["id"], None)
        for R in releases:
            pools[R][("absent_id" if "prefix" in note else "absent_name", "absent")].append(item)
    # releases whose GTF writes versions (most of its genes have one): only they have a `versioned` dataset; at the
    # others A1's few versioned lines (genes of the core gene table the GTF lacks) are no form of the release
    versioned_at = {R for R in releases if gtf_writes_versions(org, R)}
    for R in releases:
        if R not in versioned_at:
            for k in [k for k in pools[R] if k[0] == "ver"]:
                pools[R][k] = []
    # the core: genes unchanged at every chosen release, each with its id, its versioned id (where the GTF writes
    # versions) and its unique name at R
    core_by_release = {}
    common = None
    for R in releases:
        ok = set()
        for base, forms in ids_of[R].items():
            g = single(forms["id"][1]) if "id" in forms else None
            if g and forms["id"][4] == "unchanged" and g in name_of[R] and \
                    (R not in versioned_at or ("ver" in forms and single(forms["ver"][1]) == g)):
                ok.add(g)
        common = ok if common is None else common & ok
    rng0 = random.Random(seed_of(org, "core"))
    core_genes = sorted(common or [])
    rng0.shuffle(core_genes)
    core_genes = sorted(core_genes[:N_CORE])
    for R in releases:
        byg = {}
        for base, forms in ids_of[R].items():
            g = single(forms.get("id", (None, {}))[1]) if "id" in forms else None
            if g in core_genes:
                byg[g] = {"id": forms["id"], "name": name_of[R][g], **({"ver": forms["ver"]} if R in versioned_at else {})}
        core_by_release[R] = [(g, byg[g]) for g in core_genes]
    log(org, "core genes", len(core_genes))
    # duplicates' pools and the made-unique candidates
    for R in releases:
        p = pools[R]
        p[("by_gene", "id")] = by_gene[R]["id"]
        p[("by_gene", "ver")] = by_gene[R]["ver"]
        p[("pairs", "id_and_name")] = [(it, name_of[R][g]) for g, its in sorted(by_gene[R]["id"].items())
                                       for it in its[:1] if g in name_of[R] and g not in core_genes and len(its) == 1]
        p[("pairs", "versioned_and_unversioned")] = [(f["id"], f["ver"]) for b, f in sorted(ids_of[R].items())
                                                     if "id" in f and "ver" in f and f["id"][4] == "unchanged"
                                                     and R in versioned_at]
        p[("absent", "absent")] = []
    cands = set()
    for R in releases:
        for it in pools[R][("name", "name_unique")]:
            cands |= {it[0] + ".1", it[0] + "-1"}
    taken = taken_forms(org, cands, log)
    for R in releases:
        pools[R][("copy_ok", None)] = cands - taken
    # the datasets
    out_dir = os.path.join(C.family_dir("f"), "datasets", org, a.set)
    if not a.dry:
        os.makedirs(out_dir, exist_ok=True)
    all_lines, datasets, drawn_all = [], {}, {}
    for R in releases:
        for form in FORMS:
            if form == "versioned" and R not in versioned_at:
                continue  # (no versioned line of R: its GTF writes no versions)
            name = f"r{R:03d}_{form}"
            rng = random.Random(seed_of(org, name, "features"))
            p = dict(pools[R])
            p[("absent", "absent")] = p[("absent_id", "absent")] if form in ("ids", "versioned") else \
                p[("absent_name", "absent")] + (p[("absent_id", "absent")] if form == "mix" else [])
            lines, drawn = dataset_lines(form, R, p, core_by_release[R], rng)
            drawn_all[name] = {"/".join(str(x) for x in k if x is not None): v for k, v in sorted(drawn.items(), key=str)}
            feats = [q["id"] for q in lines]
            if len(set(feats)) != len(feats):
                raise SystemExit(f"{name}: a feature twice")
            seed = seed_of(org, name, "counts")
            X, U, _m = counts(seed, N_CELLS, len(feats))
            obs, var, obsm, obsp, varm, uns = metadata(org, name, form, R, seed, feats, X)
            path = os.path.join(out_dir, f"{name}.h5ad")
            if not a.dry:
                from idtrack_reproducibility.families.f import h5ad
                if os.path.exists(path):
                    raise SystemExit(f"{path} exists: a published set's datasets are never rewritten")
                h5ad.write(path, X, obs, var, layers={"unspliced": U}, obsm=obsm, obsp=obsp, varm=varm, uns=uns)
                os.chmod(path, 0o400)
            datasets[name] = dict(release=R, era=places.assembly(R), form=form, n_features=len(feats), n_cells=N_CELLS,
                                  seed=seed, file=os.path.relpath(path, C.family_dir("f")),
                                  sha256=None if a.dry else L.sha256(path),
                                  strata=dict(sorted(collections.Counter(q["stratum"] for q in lines).items())))
            all_lines += lines
            log(name, len(lines), datasets[name]["strata"])
    sources = [f"a:{s}:{L.manifest('a', org, s)['sha256']}" for s in a1_sets(org)] + \
              [f"b:{s}:{L.manifest('b', org, s)['sha256']}" for s in b1_sets(org)] + \
              [f"harness:{s}" for s in L.set_names("harness", org)] + \
              [f"a:patterns/{org}.tsv.gz:{L.sha256(os.path.join(C.family_dir('a'), 'patterns', f'{org}.tsv.gz'))}"]
    extra = dict(releases=releases, eras={str(r): places.assembly(r) for r in releases},
                 versioned_releases=sorted(versioned_at), datasets=datasets,
                 drawn=drawn_all, core_genes=len(core_genes), n_features=N_FEATURES, n_cells=N_CELLS,
                 seed=SEED, quotas={f: {"/".join(map(str, k)) if isinstance(k, tuple) else k: (list(v) if isinstance(v, tuple) else v)
                                        for k, v in q.items()} for f, q in QUOTAS.items()},
                 count_model=dict(kind="negative binomial", mean="m_f * s_c", m_f=f"lognormal(0, {SIGMA_MEAN})",
                                  s_c=f"lognormal(0, {SIGMA_CELL})", dispersion=DISPERSION, layer_unspliced="mean / 3",
                                  sampler="numpy.random.default_rng(seed).negative_binomial(r, r / (r + mean))"),
                 made_unique=dict(candidates=len(cands), taken_by_truth_store=len(taken)))
    if a.dry:
        print(json.dumps(dict(n=len(all_lines), counts=collections.Counter(q["stratum"] for q in all_lines), **extra),
                         default=str, indent=1)[:6000])
        return
    man = L.publish("f", org, a.set, all_lines, "F1: synthetic datasets with a known truth (gen_f.py docstring); each "
                    "line is a feature of a dataset with its source line's truth", __file__, sources=sources,
                    fmt="jsonl.gz", manifest_extra=extra)
    print(json.dumps({k: man[k] for k in ("name", "n", "sha256", "counts")}, indent=1))


if __name__ == "__main__":
    main()
