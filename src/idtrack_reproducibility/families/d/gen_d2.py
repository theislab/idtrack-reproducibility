"""Family D2: source inference (`infer_identifier_source` on d4c2dda, `infer_source` from 71690ab on), judged by Ensembl's
GTFs. Nothing of idtrack is read for the truth.

Usage: gen_d2.py --organism <o> [--publish <set>] [--dry] [--parts synthetic,matrices]

**Synthetic lists** (stratum `d2_synth_<size>_<noise>`): for every release R of the organism's window with a GTF, the
gene ids of R's main GTF (region primary) whose biotype Cell Ranger's `mkgtf` keeps (the reference builds' list:
protein_coding, the long non-coding biotypes of each era (lncRNA; before 94 lincRNA, antisense, sense_intronic,
sense_overlapping, 3prime_overlapping_ncRNA, bidirectional_promoter_lncRNA, macro_lncRNA), IG_* and TR_* genes and
their pseudogenes as the builds list them), unversioned (as a Cell Ranger feature table writes them); lists of 50, 500,
5,000 and all of them (seeded samples), each without noise and with 5 % and 20 % of noise of two kinds: ids of another
release R' that R does not hold (`other`), or strings that are no gene id of any release (`junk`: Ensembl-like ids
absent from every GTF of the window, and `NOT_A_GENE_<n>`), the noise replacing as many items; 5 seeds for the
samples, one draw for the lists of all ids (no sampling; one noise draw each, to keep the sets' size bounded).

**Rosetta's matrices** (stratum `d2_matrix_<ids|symbols|mixed>_<stated|measured>`): each source's feature list as rosetta asks
it (gen_d1.reduce; Ensembl ids unversioned, symbols as written), all features and only those with counts.

**Truth**: per release r of the window, the distinct items of the list that are genes of r, on the
assembly that serves r (mode `release`: a release is its own annotation): an id is a gene id of r's GTF (both regions),
or a gene of the core database that the GTF lacks (human 48 to 51, cynomolgus 112, as A*'s `TruthA` reads them); a
symbol is a name of a gene of r: the GTF's `gene_name` or the label of a cross-reference of a symbol database of r
(`SYMBOL_DBS`, Ensembl's core xrefs, `core_xref/<r>.tsv.gz`), without letter case for a mouse table written in
capitals (rosetta's rule: `MIXED_CASE`, `CAPITALS`); a list of both, as sterr2025's: either. `expect` is the exact tie
`T` (the releases where the most items are genes) and `band` the band `B` of the documented one-sided sign test at
`SIGNIFICANCE` (`band`: the best release is the lowest of `T`; another release is in `B` unless the items found at the
best and not there significantly outnumber those found there and not at the best); rule `within_band`: an answer is
right when it holds all of `T` and nothing outside `B`. A symbol list whose window holds a release without a core xref
extract is `open` (Ensembl's records read cannot say which names that release had). A line of a matrix with a stated
release (the authors' statement or the Cell Ranger reference named in rosetta's config) whose stated release is not in
`T` is `open` (the statement and the data disagree), with both recorded; `drawn` (the synthetic list's release),
`stated`, the shares at the best and the drawn release are kept.
"""
import argparse
import collections
import gzip
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

from idtrack_reproducibility.harness import harness_lib as H

FAMILY = "d"
SIZES = (50, 500, 5000, "all")
NOISE = (("none", 0.0), ("other5", 0.05), ("other20", 0.20), ("junk5", 0.05), ("junk20", 0.20))
SEEDS = (1, 2, 3, 4, 5)
LNC = {"lncRNA", "lincRNA", "antisense", "antisense_RNA", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA",
       "3prime_overlapping_ncrna", "bidirectional_promoter_lncRNA", "macro_lncRNA"}
IGTR = {"IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "IG_V_pseudogene", "IG_J_pseudogene", "IG_C_pseudogene",
        "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "TR_V_pseudogene", "TR_J_pseudogene"}
KEEP = {"protein_coding"} | LNC | IGTR
MIXED_CASE, CAPITALS = {"mus_musculus"}, 0.95  # rosetta's [gene_mapping] mixed_case_species, capitals_share: a table of
# symbols written in capitals for a species that writes them in mixed case is read without letter case (nowotschin2019)
PREFIX = {"homo_sapiens": "ENSG", "mus_musculus": "ENSMUSG", "sus_scrofa": "ENSSSCG", "macaca_mulatta": "ENSMMUG",
          "macaca_fascicularis": "ENSMFAG"}
# The release a matrix's reference names (rosetta's configs, [gene_ids] release and its comments): the authors'
# statement, or the Cell Ranger reference (refdata-gex-GRCh38-2020-A = Ensembl 98, GRCh38-3.0.0 = 93, mm10-2020-A =
# 98; hg19 references of Cell Ranger hold GENCODE 19 = Ensembl 75). Read by hand from each config, `stated_from` says where.
STATED = {
    "balboa2022_scisletmaturation": (98, "[gene_ids] release"),
    "goolam2016_4cellheterogeneity": (77, "[gene_ids] release"),
    "greywilson2025_hepatobiliarybud": (98, "[gene_ids] release; refdata-gex-GRCh38-2020-A"),
    "mima2025_ipscacinar": (101, "[gene_ids] release"),
    "pijuansala2019_gastrulationatlas": (92, "[gene_ids] release"),
    "suppinger2023_multimodalgastruloid": (99, "[gene_ids] release; GENCODE M24"),
    "sterr2025_pigneurog3chronoscrna": (101, "[gene_ids] release"),
    "yang2022_hhexscreen": (98, "[gene_ids] release; GENCODE v32"),
    "wiedenmann2021_humanpdlo": (96, "[gene_ids] release"),
    "yang2026_pigpostnatalpancreasmultiome": (101, "[gene_ids] release"),
    "ren2025_tlcembryomodel": (98, "config comment: refdata-gex-mm10-2020-A"),
    "sun2026_periembryoid": (98, "config comment: refdata-gex-GRCh38-2020-A"),
    "luo2023_ctcfcrispri": (93, "config comment: Cell Ranger GRCh38-3.0.0"),
    "ma2023_humanearlypancreas": (93, "config comment: Cell Ranger GRCh38-3.0.0"),
}


# The databases of Ensembl's core cross-references whose labels are gene symbols (a symbol counts at a
# release where Ensembl's records name a gene with it, not only the GTF's display name)
SYMBOL_DBS = {"HGNC", "MGI", "VGNC", "EntrezGene", "Uniprot_gn", "WikiGene"}
SIGNIFICANCE = 0.05  # the level of the one-sided sign test of the band (idtrack's documented `SOURCE_SIGNIFICANCE`)


def xref_symbols(org, rel):
    """The labels of the symbol databases' cross-references of a release (any object; core xref extract), or None
    where the release has no extract."""
    p = os.path.join(H.truth_dir(org), "core_xref", f"{rel}.tsv.gz")
    if not os.path.exists(p):
        return None
    out = set()
    with gzip.open(p, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[3] in SYMBOL_DBS and f[6].strip():
                out.add(f[6].strip())
    return out


def significantly_fewer(fewer, more):
    """The documented sign test: whether `fewer` (items found at the best and not at the other) exceeds `more` (found
    there and not at the best) more than chance would at SIGNIFICANCE (exact binomial up to 1,000, else normal)."""
    n = fewer + more
    if fewer <= more:
        return False
    if n <= 1000:
        p = sum(math.comb(n, k) for k in range(fewer, n + 1)) / 2 ** n
    else:
        z = (fewer - 0.5 - n / 2) / math.sqrt(n / 4)
        p = 0.5 * math.erfc(z / math.sqrt(2))
    return p < SIGNIFICANCE


def window(org):
    return [r for r in H.gtf_releases(org) if r <= 115]


class Presence:
    """Per release of the window: the gene ids (both regions), the kept ids of the main GTF, the gene names."""

    def __init__(self, org):
        self.org, self.rels = org, window(org)
        self.ids, self.ids_main, self.kept, self.names = {}, {}, {}, {}
        self.no_xref = []  # releases without a core xref extract: their names are the GTF's only
        from idtrack_reproducibility.families.astar.truth import CORE_GENES, core_gene_ids  # A*'s: genes of the core database a GTF lacks (human 48-51, cynomolgus 112)
        for r in self.rels:
            g = H.gtf_genes(org, r)
            self.ids_main[r] = set(g)  # the main GTF's ids: the pool the published lists' noise was drawn from
            core = set(core_gene_ids(org, r) or ()) if r in CORE_GENES.get(org, {}) else set()
            self.ids[r] = set(g) | core
            self.kept[r] = sorted(x for x, row in g.items() if row["region"] == "primary" and row["gene_biotype"] in KEEP)
            sym = xref_symbols(org, r)
            if sym is None:
                self.no_xref.append(r)
            self.names[r] = {row["gene_name"] for row in g.values() if row["gene_name"]} | (sym or set())
        self.every = set().union(*self.ids.values())
        import numpy as np
        self.np = np
        self.index, self.matrix = {}, {}
        folded = {r: {n.lower() for n in ns} for r, ns in self.names.items()}
        for kind, sets in (("ids", self.ids), ("symbols", self.names), ("symbols_fold", folded)):
            keys = sorted(set().union(*sets.values()))
            ix = {x: i for i, x in enumerate(keys)}
            m = np.zeros((len(keys) + 1, len(self.rels)), dtype=np.int32)  # the last row: an item of no release
            for j, r in enumerate(self.rels):
                m[[ix[x] for x in sets[r]], j] = 1
            self.index[kind], self.matrix[kind] = ix, m

    def presence(self, items, kind):
        """Boolean matrix (distinct items x releases): the item is a gene of the release (a mixed list: a gene id or a
        name of r)."""
        items = list(dict.fromkeys(items))
        m = 0
        for k in (("ids", "symbols") if kind == "mixed" else (kind,)):
            ix, mat = self.index[k], self.matrix[k]
            key = (lambda x: x.lower()) if k == "symbols_fold" else (lambda x: x)
            rows = self.np.fromiter((ix.get(key(x), len(mat) - 1) for x in items), dtype=self.np.int64, count=len(items))
            m = m + mat[rows]
        return self.np.asarray(m) > 0

    def truth(self, items, kind):
        """(the exact tie T, the band B, shares) of a list: T the releases at which the most distinct items are genes,
        B those the sign test cannot tell from the best (the lowest release of T); shares per release."""
        pres = self.presence(items, kind)
        n = max(pres.shape[0], 1)
        counts = pres.sum(axis=0)
        shares = {r: int(c) / n for r, c in zip(self.rels, counts)}
        top = counts.max()
        tie = [r for r, c in zip(self.rels, counts) if c == top]
        b = self.rels.index(tie[0])
        band = []
        for j, r in enumerate(self.rels):
            if counts[j] == 0:  # a release no item names is no origin of the list
                continue
            fewer = int((pres[:, b] & ~pres[:, j]).sum())
            more = int((pres[:, j] & ~pres[:, b]).sum())
            if j == b or not significantly_fewer(fewer, more):
                band.append(r)
        return tie, band, shares


def versionless_open(P, items, E, shares):
    """For a list of versioned ids: the releases whose GTF writes no versions (pig and human before 77) cannot say
    whether a version existed there; where the list's stable ids are genes of such a release with a share as large as
    the expected releases', the truth cannot decide whether that release belongs to the tie: the reason, else None."""
    if not hasattr(P, "versionless"):
        P.versionless = [r for r in P.rels if not any(row["gene_version"] for row in list(H.gtf_genes(P.org, r).values())[:50])]
    if not P.versionless:
        return None
    base = [x.rsplit(".", 1)[0] for x in items]
    _, _, sb = P.truth(base, "ids")
    best = shares[E[0]]
    hold = [r for r in P.versionless if sb[r] >= best]
    if hold:
        return (f"releases {hold[0]}..{hold[-1]} write no versions in their GTF and hold the list's stable ids as well as "
                f"{E[0]}..{E[-1]} hold its versions: whether they belong to the tie cannot be told from the GTFs")
    return None


def line(org, qid, stratum, items, P, kind, **extra):
    E, B, shares = P.truth(items, kind)
    q = dict(qid=qid, family=FAMILY, kind=f"d2_{extra.pop('label', kind)}", stratum=stratum, unit=extra.pop("unit"), id=f"<list of {len(items)}>",
             call="infer_source", batch=list(items), kwargs={"mode": "release"}, from_release=None, to_release=None,
             final_database=None, strategy="best", rule="within_band", cmp="versioned", expect=[str(r) for r in E],
             band=[str(r) for r in B], broad=None, best_share=round(shares[E[0]], 6), n_items=len(items), **extra)
    for k in ("drawn", "stated"):
        if q.get(k) is not None:
            q[f"{k}_share"] = round(shares.get(q[k], 0.0), 6)
            q[f"{k}_in_expected"] = q[k] in E
    if kind == "ids_versioned":
        why = versionless_open(P, items, E, shares)
        if why:
            q["open"] = why
    if kind in ("symbols", "symbols_fold", "mixed") and P.no_xref:
        q["open"] = (f"releases {','.join(map(str, P.no_xref))} have no core xref extract: which symbols named genes there "
                     "Ensembl's records read cannot say")
    if q.get("stated") is not None and q["stated"] not in E:
        q["open"] = (f"the stated release {q['stated']} fits the list worse ({q['stated_share']:.4f}) than "
                     f"{E[0]}..{E[-1]} ({q['best_share']:.4f}) by Ensembl's GTFs")
    return q


def synthetic(org, P):
    rels = P.rels
    pre = PREFIX[org]
    diff = {}  # (R, r) -> the ids of r that R does not hold, sorted (computed once)

    # the noise of another release is drawn from the main GTFs' ids (as the published sets were drawn; the truth reads
    # the GRCh37 site beside them)
    def other_ids(R, r):
        if (R, r) not in diff:
            diff[(R, r)] = sorted(P.ids_main[r] - P.ids_main[R])
        return diff[(R, r)]
    nd = {(R, r): len(P.ids_main[r] - P.ids_main[R]) for R in rels for r in rels if r != R}
    every = set().union(*P.ids_main.values())
    for R in rels:
        diff.clear()
        pool = P.kept[R]
        for size in SIZES:
            n = len(pool) if size == "all" else size
            if n > len(pool):
                continue
            for noise, frac in NOISE:
                for seed in (SEEDS if size != "all" else (1,)):
                    rnd = random.Random(f"{org}|{R}|{size}|{noise}|{seed}")
                    k = int(round(n * frac))
                    base = pool if size == "all" else rnd.sample(pool, n)
                    base = list(base)
                    other_rel = None
                    if noise.startswith("other"):
                        cands = [r for r in rels if r != R and nd[(R, r)] >= k]
                        if not cands:
                            continue
                        other_rel = rnd.choice(cands)
                        noise_items = rnd.sample(other_ids(R, other_rel), k)
                    elif noise.startswith("junk"):
                        noise_items = []
                        while len(noise_items) < k:
                            j = len(noise_items)
                            x = f"{pre}{rnd.randrange(10 ** 11):011d}" if j % 2 == 0 else f"NOT_A_GENE_{rnd.randrange(10 ** 6)}"
                            if x not in every and x not in noise_items:
                                noise_items.append(x)
                    else:
                        noise_items = []
                    if noise_items:  # the noise replaces that many items, so the list keeps its size
                        keep = base[: n - k] if size != "all" else rnd.sample(base, n - k)
                        items = keep + noise_items
                    else:
                        items = base
                    rnd.shuffle(items)
                    yield line(org, f"d2s|{R}|{size}|{noise}|{seed}", f"d2_synth_{size}_{noise}", items, P, "ids", unit=str(R),
                               drawn=R, size=str(size), noise=noise, seed=seed, other_release=other_rel)


def synthetic_more(org, P):
    """Set `d2_v2`: the same question on lists of Ensembl's display names of a release (`gene_name` of the kept genes;
    500 names, 5 seeds, and all) and of versioned ids (`<id>.<version>` of the kept genes, 500, 5 seeds; releases whose GTF
    writes versions). Truth: the releases whose GTF holds the largest share of the names (resp. id.version pairs)."""
    vers = {}
    for R in P.rels:
        g = H.gtf_genes(org, R)
        vers[R] = {x: row["gene_version"] for x, row in g.items() if row["gene_version"]}
    vsets = {r: {f"{x}.{v}" for x, v in vers[r].items()} for r in P.rels}
    P.ids_versioned = vsets
    keys = sorted(set().union(*vsets.values()))
    ix = {x: i for i, x in enumerate(keys)}
    m = P.np.zeros((len(keys) + 1, len(P.rels)), dtype=P.np.int32)
    for j, r in enumerate(P.rels):
        if vsets[r]:
            m[[ix[x] for x in vsets[r]], j] = 1
    P.index["ids_versioned"], P.matrix["ids_versioned"] = ix, m
    for R in P.rels:
        g = H.gtf_genes(org, R)
        names = sorted({g[x]["gene_name"] for x in P.kept[R] if g[x]["gene_name"]})
        for size in (500, "all"):
            for seed in (SEEDS if size != "all" else (1,)):
                if size != "all" and len(names) < size:
                    continue
                rnd = random.Random(f"{org}|{R}|names|{size}|{seed}")
                items = list(names) if size == "all" else rnd.sample(names, size)
                rnd.shuffle(items)
                yield line(org, f"d2n|{R}|{size}|{seed}", f"d2_synth_names_{size}", items, P, "symbols", unit=str(R), drawn=R,
                           size=str(size), noise="none", seed=seed)
        vs = sorted(f"{x}.{vers[R][x]}" for x in P.kept[R] if x in vers[R])
        if len(vs) >= 500:
            for seed in SEEDS:
                rnd = random.Random(f"{org}|{R}|versioned|{seed}")
                items = rnd.sample(vs, 500)
                yield line(org, f"d2v|{R}|500|{seed}", "d2_synth_versioned_500", items, P, "ids_versioned", unit=str(R), drawn=R,
                           size="500", noise="none", seed=seed)


def matrices(org, P):
    from idtrack_reproducibility.families.d import gen_d1 as G
    import pandas as pd
    for spec in G.dataset_specs(org):
        g = spec["gene_ids"]
        if g.get("mapping") != "idtrack":
            continue
        t = pd.read_csv(os.path.join(G.M02, f"{spec['dataset']}.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
        part = t[t["source"] == spec["file"]]
        qs, counted = {}, set()
        for a, cnt in zip(part["author_id"], part["counts"]):
            red = G.reduce(a, g)
            if red is None:
                continue
            q, kind, _ = red
            qs.setdefault(q, kind)
            if float(cnt or 0) > 0:
                counted.add(q)
        kinds = collections.Counter(qs.values())
        lkind = "mixed" if len(kinds) > 1 else "ids" if "id" in kinds else "symbols"
        cased = [x for x in qs if x.upper() != x.lower()]
        capitals = sum(x == x.upper() for x in cased) / max(len(cased), 1)
        truth_kind = "symbols_fold" if lkind == "symbols" and org in MIXED_CASE and capitals >= CAPITALS else lkind
        stated, how = STATED.get(spec["dataset"], (None, None))
        src = sorted(set(part["source_release"]))
        si = G.dataset_specs(org).index(spec)
        for which, items in (("all", list(qs)), ("counted", [q for q in qs if q in counted])):
            if not items:
                continue
            yield line(org, f"d2m|{spec['dataset']}|{si}|{which}", f"d2_matrix_{lkind}_{'stated' if stated else 'measured'}", items, P,
                       truth_kind, label=lkind, unit=spec["dataset"], capitals_share=round(capitals, 4), truth_kind=truth_kind, dataset=spec["dataset"], source=spec["file"], features=which, stated=stated,
                       stated_from=how, rosetta_source_release=",".join(src), mixed_kinds=dict(kinds))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--publish")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--parts", default="synthetic,matrices", help="synthetic, matrices (set d2_v1), more (set d2_v2)")
    a = ap.parse_args()
    P = Presence(a.organism)
    out = []
    if "synthetic" in a.parts:
        out += list(synthetic(a.organism, P))
    if "matrices" in a.parts:
        out += list(matrices(a.organism, P))
    if "more" in a.parts:
        out += list(synthetic_more(a.organism, P))
    c = collections.Counter(q["stratum"] for q in out)
    for k, v in sorted(c.items()):
        print(f"  {k:40s} {v}")
    m = [q for q in out if q["stratum"].startswith("d2_matrix")]
    for q in m:
        print(q["qid"], q["n_items"], "expect", q["expect"][:3], "..", q["expect"][-1:], "share", q["best_share"], "stated", q.get("stated"),
              q.get("stated_share"), "rosetta", q["rosetta_source_release"], "open" if q.get("open") else "")
    syn = [q for q in out if q["stratum"].startswith("d2_synth") and q.get("drawn_in_expected") is not None]
    print("synthetic: drawn release among the expected:", sum(q["drawn_in_expected"] for q in syn), "of", len(syn),
          "; expected a single release:", sum(len(q["expect"]) == 1 for q in syn))
    if a.publish and not a.dry:
        L.publish(FAMILY, a.organism, a.publish, out, description="D2: source inference, synthetic lists of every release and rosetta's "
                  "matrices, judged by the releases whose GTF holds the largest share of the list (gen_d2.py docstring)",
                  generator=os.path.abspath(__file__), sources=["gtf (genes tables of every release of the window)"],
                  manifest_extra=dict(sizes=[str(s) for s in SIZES], noise=[n for n, _ in NOISE], seeds=list(SEEDS)))
        print("published", a.publish)


if __name__ == "__main__":
    main()
