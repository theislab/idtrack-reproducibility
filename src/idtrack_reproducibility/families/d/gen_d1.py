"""Family D1: the feature lists of rosetta's 33 authors' matrices, each feature asked as rosetta asks, judged by place.

Usage: gen_d1.py --organism <o> [--publish <set>] [--dry] [--out <truth table .tsv.gz>]

Inputs (read only): rosetta's per-feature tables of its matrices (`r3_build_matrix_locus_dev/out/m02/<dataset>.tsv.gz`:
every author feature of a source mapped with idtrack, its counts in the source's cells, the source release rosetta used:
stated by the authors or measured) and the datasets' configs (`rosetta/stages/data_preparation/<dataset>/configs/
prepare.toml`: the declared type, `collapse_suffix`, `other_features`, the stated release). Nothing of idtrack is read
for the truth (the source release rosetta measured with idtrack's graph is part of the question, not of the truth).

**The question** (as rosetta asks, `rosetta.preprocessing.author_matrix.gene_ids` and `genes.map_features`): the
feature reduced as rosetta reduces it (the `collapse_suffix` cut, a version removed for `ensembl_gene_versioned`, an
Ensembl id or a symbol by the declared type), from the source release R to 115 (an Ensembl id that R's GTF does not
hold: without `from_release`, as rosetta asks it), `strategy` both. xiang2020's versioned ids are asked a second time
with their version (stratum `d1_id_versioned`).

**The truth** (locus; docs/methods.md, D, D1): the feature's gene model in Ensembl's GTF of R (merged exons; a symbol's
model: the gene(s) whose `gene_name` it is at R, exactly, else without letter case) against the models of every gene
of Ensembl's 115 GTF (main file: chromosomes and unplaced scaffolds) by exon bases shared on the strand. A model of
another assembly (human GRCh37, mouse GRCm38) is mapped exon by exon through Ensembl's chain file (`truth/<o>/chain/`):
every aligned block, the pieces on the sequence and strand that take most bases kept, unmapped below half of the
model's bases (the harness's `MIN_ALIGNED`). rosetta's thresholds: a model holds a gene where they share at least
`MIN_SHARED` (1) base; the one gene a model holds is its gene where it holds at least `MIN_SHARE` (0.5) of the model's
bases. Classes (`truth_class`):

- `one`: the model holds one gene, with at least half of the model's bases: expected (rule `set`, `cmp` base);
- `dominant`: the model holds several genes, and one of them is told by either of two rules that do not name two
  different genes: `bases` (the gene that shares most of the model's bases shares at least half of them, every other
  less than half of them and less than half of what that gene shares: `DOMINANCE`) or `mutual` (exactly one gene
  shares at least half of the model's bases and at least half of its own, and no fewer of the model's bases than any
  other: a readthrough beside the gene itself; revision `d1_r2`);
  expected, in its own stratum (a readthrough, a small RNA in an exon, a novel overlapping gene beside the gene
  itself, a gene grown since); the design's strict reading is the `one` stratum alone;
- `several`, `below_share`, `none`: no gene can be told; the line is `open` with `open_within` = the genes the model
  holds and every gene whose span overlaps the model's span on the strand: an answer of nothing or within it is
  undetermined, a gene elsewhere is wrong (scorer rule `open_within`, F0);
- `patch_at_source` (an id that R holds only on a patch or alternative locus: its place projected to the reference
  sequence through Ensembl's assembly exceptions, the harness's `assembly_map.Exceptions`): open within the genes of
  115 at that place, the id itself and its broad history walk;
- `unmapped` (no chain place), `not_at_source` (an id R's GTF does not hold), `name_not_at_source` (a symbol that no
  gene of R carries as its name), `made_unique` (`X.1`, `X-1` of a name `X` of R: open within `X`'s genes' places):
  open.

A symbol of several genes at R: every gene's own locus; when every one is `one` or `dominant`, the union is expected
(`set_exact`, the harness's rule for a name of several genes); otherwise open within all their places.

Beside the truth each line carries (extra fields): `dataset`, `sources`, `counts` (the counts of its author features
in their sources' cells, summed), `n_features`, `truth_class`, `locus_hits` (gene: shared bases), `model_bases`,
`assembly_mapped`, for ids Ensembl's history walk R -> 115 (`history`: the strict walk with living ids;
`history_broad`), for symbols `genes_at_source`, and rosetta's own result on its r3 code (`rosetta_target`,
`rosetta_class`), for the comparison of both with the place.
"""
import argparse
import bisect
import collections
import glob
import gzip
import json
import os
import re
import sys
import time
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

from idtrack_reproducibility.harness import assembly_map as AM
from idtrack_reproducibility.harness import harness_lib as H

from idtrack_reproducibility import config as _CFG
ROSETTA = _CFG.get().source("rosetta")  # the atlas project's checkout (its preparation configs)
M02 = _CFG.get().source("rosetta_m02")  # its per-feature locus tables
ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{11}$")
MADE_UNIQUE = re.compile(r"^(.+)[.-]([1-9]\d*)$")
MIN_SHARED = 1      # rosetta configs/preparation.toml [gene_mapping] locus_min_shared_bases
MIN_SHARE = 0.5     # rosetta configs/preparation.toml [gene_mapping] locus_min_share
DOMINANCE = 0.5     # `dominant`: every other gene shares less than this part of what the first one shares
TARGET = 115
FAMILY = "d"


# ---------------------------------------------------------------------------------------------- the datasets

def dataset_specs(organism=None):
    """One row per (dataset, source file) of rosetta's authors' matrices: species, the source's gene_ids config."""
    out = []
    for p in sorted(glob.glob(os.path.join(ROSETTA, "stages/data_preparation/*/configs/prepare.toml"))):
        with open(p, "rb") as fh:
            cfg = tomllib.load(fh)
        if cfg.get("path", {}).get("kind") != "author_matrix":
            continue
        u = p.split("/")[-3]
        species = cfg["inputs"]["annotation"].split("ensembl115_")[1]
        if organism and species != organism:
            continue
        srcs = cfg["path"]["source"]
        for s in (srcs if isinstance(srcs, list) else [srcs]):
            out.append(dict(dataset=u, species=species, file=s["file"], gene_ids=s.get("gene_ids", cfg.get("gene_ids", {}))))
    return out


def reduce(author_id, g):
    """(query, kind, versioned form or None) as rosetta reduces an author feature; None for an other_feature."""
    if author_id in set(g.get("other_features", [])):
        return None
    t = g["type"]
    q = author_id
    versioned = None
    if g.get("collapse_suffix") and t in ("ensembl_gene", "ensembl_gene_versioned"):
        q = re.sub(g["collapse_suffix"], "", q)
    if t == "ensembl_gene_versioned":
        versioned = q
        q = re.sub(r"\.\d+$", "", q)
    if t in ("ensembl_gene", "ensembl_gene_versioned"):
        kind = "id"
    elif t == "symbol_else_id":
        kind = "id" if ENSEMBL.match(q) else "symbol"
    else:
        kind = "symbol"
    return q, kind, versioned


def features(organism):
    """{(dataset, kind, query, R): record} over the organism's matrices, author features grouped."""
    import pandas as pd
    groups = {}
    for spec in dataset_specs(organism):
        t = pd.read_csv(os.path.join(M02, f"{spec['dataset']}.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
        part = t[t["source"] == spec["file"]]
        if not len(part):
            raise SystemExit(f"{spec['dataset']}: m02 has no rows of source {spec['file']}")
        g = spec["gene_ids"]
        if g.get("mapping") != "idtrack":
            continue
        for a, cnt, rel, rt, rc in zip(part["author_id"], part["counts"], part["source_release"], part["target_id"], part["class"]):
            red = reduce(a, g)
            if red is None:
                continue
            q, kind, versioned = red
            R = int(rel)
            for k, qq in ((kind, q), ("id_versioned", versioned)):
                if qq is None:
                    continue
                key = (spec["dataset"], k, qq, R)
                rec = groups.get(key)
                if rec is None:
                    rec = groups[key] = dict(dataset=spec["dataset"], kind=k, query=qq, R=R, counts=0, n_features=0, sources=set(),
                                             stated=g.get("release"), rosetta_target=set(), rosetta_class=set(), base=q)
                rec["counts"] += int(float(cnt or 0))
                rec["n_features"] += 1
                rec["sources"].add(spec["file"])
                if rt:
                    rec["rosetta_target"].add(rt)
                rec["rosetta_class"].add(rc)
    return groups


M16 = os.path.join(os.path.dirname(M02), "m16")


def rosetta_current(organism):
    """{(dataset, kind, query, R): (targets, classes, which)}: rosetta's own mapping of the features by its current rule:
    the r3 run of the place step (`out/m16/<dataset>.tsv.gz`, the 18 matrices whose source shares the target's
    assembly) where it exists, else the id route of `out/m02` (the other 15: no place step applies to them)."""
    import pandas as pd
    out = {}
    for spec in dataset_specs(organism):
        g = spec["gene_ids"]
        if g.get("mapping") != "idtrack":
            continue
        m02 = pd.read_csv(os.path.join(M02, f"{spec['dataset']}.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
        rel = dict(zip(m02["source"], m02["source_release"]))
        p16 = os.path.join(M16, f"{spec['dataset']}.tsv.gz")
        t, which = (pd.read_csv(p16, sep="\t", dtype=str, keep_default_na=False), "m16") if os.path.exists(p16) else (m02, "m02")
        part = t[t["source"] == spec["file"]]
        for a, rt, rc in zip(part["author_id"], part["target_id"], part["class"]):
            red = reduce(a, g)
            if red is None:
                continue
            q, kind, versioned = red
            for k, qq in ((kind, q), ("id_versioned", versioned)):
                if qq is None:
                    continue
                key = (spec["dataset"], k, qq, int(rel[spec["file"]]))
                tg, cl, _ = out.get(key, (set(), set(), which))
                if rt:
                    tg.add(rt)
                cl.add(rc)
                out[key] = (tg, cl, which)
    return out


# ---------------------------------------------------------------------------------------------- the models

class Models:
    """Gene models of one release (truth/<o>/exons/<release>.tsv.gz, fetch_exons.py)."""

    def __init__(self, org, rel):
        p = os.path.join(H.truth_dir(org), "exons", f"{rel}.tsv.gz")
        if not os.path.exists(p):
            raise SystemExit(f"{p} missing: fetch_exons.py {org}:{rel}")
        self.genes, self.names, self.lower = {}, collections.defaultdict(list), collections.defaultdict(list)
        with gzip.open(p, "rt") as fh:
            next(fh)
            for line in fh:
                gid, ver, name, bt, seq, strand, s, e, bases, exons, places = line.rstrip("\n").split("\t")
                ex = [tuple(map(int, x.split("-"))) for x in exons.split(",")]
                self.genes[gid] = dict(name=name, seq=seq, strand=strand, start=int(s), end=int(e), bases=int(bases), exons=ex,
                                       multi=";" in places, version=ver, biotype=bt)
                if name:
                    self.names[name].append(gid)
                    self.lower[name.lower()].append(gid)
        self._ix = self._sx = None

    def _index(self):
        ix, sx = collections.defaultdict(list), collections.defaultdict(list)
        for gid, g in self.genes.items():
            for s, e in g["exons"]:
                ix[(g["seq"], g["strand"])].append((s, e, gid))
            sx[(g["seq"], g["strand"])].append((g["start"], g["end"], gid))
        self._ix, self._sx = {}, {}
        for store, src in ((self._ix, ix), (self._sx, sx)):
            for k, bl in src.items():
                bl.sort()
                m, maxend = -1, []
                for b in bl:
                    m = max(m, b[1])
                    maxend.append(m)
                store[k] = (bl, maxend)

    def overlaps(self, seq, strand, exons):
        """{gene: exon bases shared on the strand} with a model given by its merged exons (1-based inclusive)."""
        if self._ix is None:
            self._index()
        out = collections.Counter()
        bl, maxend = self._ix.get((seq, strand), ([], []))
        for s, e in exons:
            i = bisect.bisect_left(maxend, s)
            while i < len(bl) and bl[i][0] <= e:
                a, b = max(s, bl[i][0]), min(e, bl[i][1])
                if a <= b:
                    out[bl[i][2]] += b - a + 1
                i += 1
        return out

    def spans(self, seq, strand, start, end):
        """Genes whose span overlaps [start, end] on the strand."""
        if self._sx is None:
            self._index()
        bl, maxend = self._sx.get((seq, strand), ([], []))
        i = bisect.bisect_left(maxend, start)
        out = set()
        while i < len(bl) and bl[i][0] <= end:
            if bl[i][1] >= start:
                out.add(bl[i][2])
            i += 1
        return out


def map_exons(chain, seq, strand, exons, bases):
    """A model's exons through a chain (Ensembl's chain files, 0-based blocks): (seq, strand, exons, aligned share) or
    (None, reason). Every aligned piece is kept on the sequence and strand that take most bases."""
    if seq not in chain.index:
        return None, f"sequence {seq} not in the chain"
    bl, maxend = chain.index[seq]
    groups = collections.defaultdict(list)
    for s, e in exons:
        s0, e0 = s - 1, e
        i = bisect.bisect_right(maxend, s0)
        while i < len(bl) and bl[i][0] < e0:
            t1, t2, qname, qstrand, qsize, q1 = bl[i]
            a, b = max(t1, s0), min(t2, e0)
            if a < b:
                qa, qb = q1 + (a - t1), q1 + (b - t1)
                if qstrand == "-":
                    qa, qb = qsize - qb, qsize - qa
                groups[(qname, qstrand)].append((qa + 1, qb))
            i += 1
    if not groups:
        return None, "no aligned block"
    (qname, qstrand), pieces = max(groups.items(), key=lambda kv: sum(b - a + 1 for a, b in kv[1]))
    n = sum(b - a + 1 for a, b in pieces)
    share = n / bases
    if share < AM.MIN_ALIGNED:
        return None, f"{share:.0%} of the model's bases aligned"
    new_strand = strand if qstrand == "+" else {"+": "-", "-": "+"}[strand]
    pieces.sort()
    merged = []
    for a, b in pieces:
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return (qname, new_strand, [tuple(x) for x in merged], round(share, 4)), ""


class Locus:
    """The place of a gene of a source release at 115, with rosetta's thresholds."""

    def __init__(self, org, target=TARGET):
        self.org, self.T = org, target
        self.asm = H.gtf_assemblies(org)
        self.target = Models(org, target)
        self.models = {}
        self.mapper = AM.Mapper(org)
        self.cache = {}

    def source(self, rel):
        if rel not in self.models:
            self.models[rel] = Models(self.org, rel)
        return self.models[rel]

    def of(self, gid, rel):
        """dict(cls, expect, open_within, hits, model_bases, mapped) of gene gid of release rel."""
        key = (gid, rel)
        if key in self.cache:
            return self.cache[key]
        g = self.source(rel).genes.get(gid)
        out = self._of(g, rel) if g else dict(cls="not_at_source", expect=[], open_within=None, hits={}, model_bases=0, mapped=None)
        self.cache[key] = out
        return out

    def _of(self, g, rel):
        seq, strand, exons, bases = g["seq"], g["strand"], g["exons"], g["bases"]
        mapped = None
        if self.asm.get(rel) != self.asm.get(self.T):
            p = self.mapper.path(self.asm[rel], self.asm[self.T])
            if p is None:
                return dict(cls="unmapped", expect=[], open_within=None, hits={}, model_bases=bases,
                            mapped=f"no chain from {self.asm[rel]} to {self.asm[self.T]}")
            share = 1.0
            for x, y in zip(p, p[1:]):
                res, why = map_exons(self.mapper.chain(x, y), seq, strand, exons, sum(b - a + 1 for a, b in exons))
                if res is None:
                    return dict(cls="unmapped", expect=[], open_within=None, hits={}, model_bases=bases, mapped=f"{x} -> {y}: {why}")
                seq, strand, exons, s = res
                share = min(share, s)
            mapped = dict(via=" -> ".join(p), seq=seq, strand=strand, start=exons[0][0], end=exons[-1][1], aligned=share)
        T = self.target
        hits = {k: v for k, v in T.overlaps(seq, strand, exons).items() if v >= MIN_SHARED}
        near = sorted(set(hits) | T.spans(seq, strand, exons[0][0], exons[-1][1]))
        res = dict(hits=dict(sorted(hits.items(), key=lambda kv: (-kv[1], kv[0]))), model_bases=bases, mapped=mapped,
                   multi=g["multi"])
        if not hits:
            return dict(res, cls="none", expect=[], open_within=near)
        if len(hits) == 1:
            (t, b), = hits.items()
            if b / bases >= MIN_SHARE:
                return dict(res, cls="one", expect=[t], open_within=None)
            return dict(res, cls="below_share", expect=[], open_within=near)
        order = sorted(hits, key=lambda t: (-hits[t], t))
        first, second = hits[order[0]], hits[order[1]]
        by_bases = {order[0]} if first / bases >= MIN_SHARE and second / bases < MIN_SHARE and second < DOMINANCE * first else set()
        mutual = [t for t, b in hits.items() if b / bases >= MIN_SHARE and b / T.genes[t]["bases"] >= MIN_SHARE]
        # the one gene that shares half of the model's bases and half of its own, and no fewer of the model's bases than
        # any other gene (a readthrough beside the gene; not a small new gene inside a gene that grew: CDC42 at 95)
        by_mutual = set(mutual) if len(mutual) == 1 and hits[mutual[0]] == first else set()
        chosen = by_bases | by_mutual
        if len(chosen) == 1:
            return dict(res, cls="dominant", expect=sorted(chosen), open_within=None, dominant_by=",".join(
                k for k, v in (("bases", by_bases), ("mutual", by_mutual)) if v))
        return dict(res, cls="several", expect=[], open_within=near)


# ---------------------------------------------------------------------------------------------- the lines

OPEN_TEXT = {
    "several": "the feature's model holds several genes of 115",
    "below_share": "the one gene the model holds has under half of the model's exon bases",
    "none": "the feature's model holds no gene of 115",
    "unmapped": "the model has no place on the assembly of 115",
    "not_at_source": "the id is no gene of the source release's GTF: no model",
    "patch_at_source": "the id is a gene of the source release only on a patch or alternative locus",
    "name_not_at_source": "the symbol is the name of no gene in the source release's GTF: no model",
    "made_unique": "a name made unique by the matrix's maker (X.1, X-1): no gene of the source release carries it",
    "name_mixed": "a name of several genes at the source release, not all of them with one gene at 115",
}


def symbol_truth(loc, sym, R):
    """(class, expect, open_within, genes at R, match, per-gene loci) of a symbol read at R."""
    M = loc.source(R)
    genes, match = M.names.get(sym, []), "exact"
    if not genes:
        genes, match = M.lower.get(sym.lower(), []), "case"
    if not genes:
        m = MADE_UNIQUE.match(sym)
        base = m and (M.names.get(m.group(1)) or M.lower.get(m.group(1).lower()))
        if base:
            per = {g: loc.of(g, R) for g in base}
            within = sorted({t for x in per.values() for t in list(x["hits"]) + list(x["expect"]) + list(x["open_within"] or [])})
            return "made_unique", [], within, sorted(base), "made_unique", per
        return "name_not_at_source", [], None, [], "none", {}
    per = {g: loc.of(g, R) for g in sorted(genes)}
    if len(genes) == 1:
        x = per[genes[0]]
        return x["cls"], x["expect"], x["open_within"], sorted(genes), match, per
    if all(x["cls"] in ("one", "dominant") for x in per.values()):
        return "name_several", sorted({t for x in per.values() for t in x["expect"]}), None, sorted(genes), match, per
    within = sorted({t for x in per.values() for t in list(x["hits"]) + list(x["expect"]) + list(x["open_within"] or [])})
    return "name_mixed", [], within, sorted(genes), match, per


def build(org, publish=None, dry=False, table=None):
    t0 = time.time()
    groups = features(org)
    print(f"{org}: {len(groups)} lines from {len({k[0] for k in groups})} datasets ({time.time() - t0:.0f} s)", flush=True)
    loc = Locus(org)
    tr = None
    if any(k[1] in ("id", "id_versioned") for k in groups):
        from idtrack_reproducibility.families.astar.truth import TruthA
        tr = TruthA(org, TARGET)
        print(f"{org}: history loaded ({time.time() - t0:.0f} s)", flush=True)
    asm_t = loc.asm[TARGET]
    exc = None
    walks = {}
    out, counts = [], collections.Counter()
    for (ds, kind, q, R), rec in sorted(groups.items()):
        same = loc.asm.get(R) == asm_t
        line = dict(family=FAMILY, kind=f"d1_{'symbol' if kind == 'symbol' else 'id'}", id=q, to_release=TARGET, final_database=None,
                    strategy="both", cmp="base", rule="set", broad=None, dataset=ds, sources=sorted(rec["sources"]),
                    counts=rec["counts"], n_features=rec["n_features"], source_release=R, stated_release=rec["stated"],
                    same_assembly=same, rosetta_target=sorted(rec["rosetta_target"]), rosetta_class=sorted(rec["rosetta_class"] - {""}))
        if kind in ("id", "id_versioned"):
            base = rec["base"]
            at_r = R in tr.present.get(base, ()) if tr else True
            line["from_release"] = R if at_r else None
            x = loc.of(base, R)
            cls, E, within = x["cls"], x["expect"], x["open_within"]
            line.update(locus_hits=x["hits"], model_bases=x["model_bases"], assembly_mapped=x.get("mapped"), dominant_by=x.get("dominant_by"))
            if (base, R) not in walks:
                if at_r:
                    hs = tr.walk_fwd(base, R, TARGET, True, living=True)
                    hb = tr.walk_fwd(base, R, TARGET, False, living=True, resurrect=True)
                else:
                    hs = hb = None
                walks[(base, R)] = (sorted(hs) if hs is not None else None, sorted(hb) if hb is not None else None)
            line["history"], line["history_broad"] = walks[(base, R)]
            line["unit"] = base
            row = tr.genes.get(R, {}).get(base) if (tr and at_r and cls == "not_at_source") else None
            if row is not None and row.get("region") not in ("primary", "core_only"):
                cls, place = "patch_at_source", None
                if exc is None:
                    exc = AM.Exceptions(org)
                p = exc.project(dict(seqname=row["seqname"], start=int(row["start"]), end=int(row["end"]), strand=row["strand"]), R)
                if p.get("unplaced"):
                    place = dict(unplaced=p["unplaced"])
                elif loc.asm[R] != asm_t:
                    place = loc.mapper.map(dict(p, assembly=loc.asm[R]), asm_t)
                else:
                    place = dict(p, assembly=asm_t)
                near = set() if "unmapped" in place or "unplaced" in place else \
                    loc.target.spans(place["seqname"], place["strand"], int(place["start"]), int(place["end"]))
                within = sorted(near | {base} | set(line["history_broad"] or []))
                line["patch_place"] = place
            if kind == "id_versioned":
                v = tr.ver(base, R) if tr else ""
                line["version_at_source"] = v
                line["version_asked"] = q.rsplit(".", 1)[1]
                if v and v != line["version_asked"]:
                    cls, E, within = "version_not_at_source", [], None
            stratum = "d1_id_versioned" if kind == "id_versioned" else f"d1_id_{'same' if same else 'cross'}"
            qid = f"d1|{ds}|{kind}|{q}|{R}"
        else:
            cls, E, within, genes, match, per = symbol_truth(loc, q, R)
            line["from_release"] = R
            line.update(genes_at_source=genes, name_match=match,
                        locus_hits={g: x["hits"] for g, x in per.items()}, assembly_mapped={g: x.get("mapped") for g, x in per.items()} or None)
            line["unit"] = q.upper()
            stratum = f"d1_symbol_{'same' if same else 'cross'}"
            qid = f"d1|{ds}|symbol|{q}|{R}"
            if cls == "name_several":
                line["rule"] = "set_exact"
        line["truth_class"] = cls
        line["stratum"] = f"{stratum}_{'dominant' if cls == 'dominant' else 'one' if cls in ('one', 'name_several') else 'open'}"
        line["qid"] = qid
        line["expect"] = list(E)
        if not E:
            line["open"] = OPEN_TEXT.get(cls, cls) if cls != "version_not_at_source" else \
                "the version asked is not the id's version at the source release"
            if within is not None:
                line["open_within"] = list(within)
        counts[(line["stratum"], cls)] += 1
        out.append(line)
    print(f"{org}: {len(out)} lines ({time.time() - t0:.0f} s)", flush=True)
    for (st, cls), n in sorted(counts.items()):
        print(f"  {st:32s} {cls:22s} {n}", flush=True)
    if table:
        with gzip.open(table, "wt") as fo:
            for line in out:
                fo.write(C.dumps(line) + "\n")
    if publish and not dry:
        srcs = sorted({f"exons/{r}" for r in list(loc.models) + [TARGET]} | {"chain", "core (stable_id_event, mapping_session 115)", "gtf (genes tables)"})
        L.publish(FAMILY, org, publish, out, description="D1: the features of rosetta's authors' matrices, as rosetta asks them (source "
                  "release to 115), judged by the place of their model in the source release's GTF against the 115 GTF (gen_d1.py docstring)",
                  generator=os.path.abspath(__file__), sources=srcs,
                  manifest_extra=dict(m02=M02, rosetta_commit=C.git_commit(ROSETTA)[0], min_shared=MIN_SHARED, min_share=MIN_SHARE,
                                      dominance=DOMINANCE))
        print(f"{org}: published {publish}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--publish")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--table")
    a = ap.parse_args()
    build(a.organism, a.publish, a.dry, a.table)
