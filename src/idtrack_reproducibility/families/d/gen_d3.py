"""Family D3: published marker lists, every symbol asked as a user of the list asks it. Nothing of idtrack is read.

Usage: gen_d3.py --organism <homo_sapiens|mus_musculus> [--publish <set>] [--dry]

**The lists** (`<families>/d/inputs/`, fetch_d3.py; human and mouse):

- PanglaoDB (`PanglaoDB_markers_27_Mar_2020.tsv.gz`, its `official gene symbol`; rows of `Hs` for human, `Mm` for
  mouse, written in capitals for both): dated 2020-03-27;
- CellMarker 2.0 (`all_cell_marker.txt`, the `marker` column: the string as the paper wrote it; `Human`, `Mouse`):
  each row dated by its paper's `year` (30 June of that year), so one string can be asked at several releases;
  CellMarker's own reading (`symbol`, NCBI `gene_id`) is kept beside it;
- MSigDB 2024.1 cell type signatures (`c8` human, `m8` mouse; dated by the files, 2024-08-09);
- rosetta's lineage table (`stages/reference_data/gene_sets/lineage_genes.tsv`, `gene_symbol_human`; human; read at 115,
  rosetta's release; rosetta's own gene id kept beside it).

**The question**: each distinct symbol to 115, from the release current at the list's date (`from`: the latest release
published before the date: Ensembl's archive index, and for the releases it does not list the date of the human GTF
plus one month) and without a release (`default`).

**Truth** (Ensembl's: principle 2's first kind): the genes whose GTF `gene_name` the symbol is at the release read
(exactly as written; else without letter case, `match: case`: PanglaoDB's mouse rows, papers' capitals), walked to
115 by Ensembl's history (the harness's strict walk with living ids); several genes: rule `set` (`best` any of them,
`all` all of them). `default` reads it as idtrack's rule `start_release` does: at 115 when the symbol is a name there, else at the
latest release where it is one. A symbol that is no name at the release read is open; one whose genes Ensembl retired
without successor is open within the genes that carry it at 115 (an answer elsewhere is wrong).

**The authority** (principle 2's second kind, beside the truth, never the judge): HGNC (approved symbol, else the one
entry whose previous symbol it is, else whose alias, else the entry a withdrawn symbol was merged into; its Ensembl gene
id) and MGI (MRK_List2 current symbols; MRK_List1 withdrawn symbols with their current marker; MRK_ENSEMBL's Ensembl
gene id), kept as `authority` and `authority_how`, measured by report_d3.py.
"""
import argparse
import collections
import datetime
import gzip
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

from idtrack_reproducibility.harness import harness_lib as H

FAMILY = "d"
INPUTS = os.path.join(C.BENCH, "d", "inputs")
from idtrack_reproducibility import config as _CFG
ROSETTA = _CFG.get().source("rosetta")  # the atlas project's checkout (its lineage gene table)
TARGET = 115


def release_dates():
    """{release: public date}: the archive index's month (its first day), else the human GTF's date plus one month."""
    out = {}
    for line in open(os.path.join(INPUTS, "ensembl_releases.tsv")).read().splitlines()[1:]:
        r, m = line.split("\t")
        out[int(r)] = datetime.date.fromisoformat(m + "-01")
    for line in open(os.path.join(INPUTS, "ensembl_gtf_dates.tsv")).read().splitlines()[1:]:
        r, d = line.split("\t")
        if d and int(r) not in out:
            out[int(r)] = datetime.date.fromisoformat(d) + datetime.timedelta(days=31)
    return out


def release_at(date, dates, window):
    ok = [r for r in window if r in dates and dates[r] <= date]
    return max(ok) if ok else min(window)


def lists(org):
    """[(source, symbol, date, extra)] of the organism's marker lists."""
    out = []
    sp = {"homo_sapiens": ("Hs", "Human"), "mus_musculus": ("Mm", "Mouse")}[org]
    with gzip.open(os.path.join(INPUTS, "panglaodb", "PanglaoDB_markers_27_Mar_2020.tsv.gz"), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_sp, i_sym = head.index("species"), head.index("official gene symbol")
        n = collections.Counter()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if sp[0] in f[i_sp].split() and f[i_sym].strip():
                n[f[i_sym].strip()] += 1
    out += [("panglaodb", s, datetime.date(2020, 3, 27), dict(n_rows=k)) for s, k in n.items()]
    with zipfile.ZipFile(os.path.join(INPUTS, "cellmarker", "all_cell_marker.zip")) as z, z.open("all_cell_marker.txt") as raw:
        import io
        fh = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        head = fh.readline().rstrip("\n").split("\t")
        ix = {k: head.index(k) for k in ("species", "marker", "symbol", "gene_id", "year")}
        n, cm = collections.Counter(), {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(head) or f[ix["species"]] != sp[1]:
                continue
            m = f[ix["marker"]].strip()
            y = f[ix["year"]].strip()
            if not m or not y or y.lower() == "nan":
                continue
            key = (m, int(float(y)))
            n[key] += 1
            cm.setdefault(key, (f[ix["symbol"]].strip(), f[ix["gene_id"]].strip().replace(".0", "")))
    out += [("cellmarker", m, datetime.date(y, 6, 30), dict(n_rows=k, year=y, cellmarker_symbol=cm[(m, y)][0], cellmarker_ncbi=cm[(m, y)][1]))
            for (m, y), k in n.items()]
    gmt = {"homo_sapiens": "c8.all.v2024.1.Hs.symbols.gmt", "mus_musculus": "m8.all.v2024.1.Mm.symbols.gmt"}[org]
    n = collections.Counter()
    for line in open(os.path.join(INPUTS, "msigdb", gmt)):
        for s in line.rstrip("\n").split("\t")[2:]:
            if s.strip():
                n[s.strip()] += 1
    out += [("msigdb", s, datetime.date(2024, 8, 9), dict(n_sets=k)) for s, k in n.items()]
    if org == "homo_sapiens":
        import csv
        rows = list(csv.DictReader(open(os.path.join(ROSETTA, "stages/reference_data/gene_sets/lineage_genes.tsv")), delimiter="\t"))
        ids = {}
        for r in rows:
            if r["species"] == "homo_sapiens":
                ids[r["gene_symbol_human"]] = r["gene_ids"]
        out += [("rosetta_lineage", s, None, dict(rosetta_gene=ids[s])) for s in sorted(ids)]
    return out


class Names:
    """Per release of the window: GTF gene_name -> genes (primary region), exactly and without letter case."""

    def __init__(self, org):
        self.rels = [r for r in H.gtf_releases(org) if r <= TARGET]
        self.exact, self.lower = {}, {}
        for r in self.rels:
            ex, lo = collections.defaultdict(set), collections.defaultdict(set)
            for gid, row in (H.gtf_genes(org, r, "primary") or {}).items():
                if row["gene_name"]:
                    ex[row["gene_name"]].add(gid)
                    lo[row["gene_name"].lower()].add(gid)
            self.exact[r], self.lower[r] = ex, lo

    def genes(self, s, r):
        g = self.exact[r].get(s)
        if g:
            return sorted(g), "exact"
        g = self.lower[r].get(s.lower())
        return (sorted(g), "case") if g else ([], "none")


class Authority:
    def __init__(self, org):
        self.org = org
        t = H.truth_dir(org)
        self.by = {}
        if org == "homo_sapiens":
            import csv
            rows = list(csv.DictReader(open(os.path.join(t, "hgnc", "hgnc_complete_set.txt")), delimiter="\t"))
            ens = {r["hgnc_id"]: r["ensembl_gene_id"] for r in rows}
            sym = {r["symbol"]: [r["hgnc_id"]] for r in rows}
            prev, alias = collections.defaultdict(set), collections.defaultdict(set)
            for r in rows:
                for x in r["prev_symbol"].strip('"').split("|"):
                    if x:
                        prev[x].add(r["hgnc_id"])
                for x in r["alias_symbol"].strip('"').split("|"):
                    if x:
                        alias[x].add(r["hgnc_id"])
            wd = collections.defaultdict(set)
            for r in csv.DictReader(open(os.path.join(t, "hgnc", "withdrawn.txt")), delimiter="\t"):
                for part in (r.get("MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)") or "").split(","):
                    p = part.strip().split("|")
                    if len(p) >= 1 and p[0].startswith("HGNC:"):
                        wd[r["WITHDRAWN_SYMBOL"]].add(p[0])
            self.tables = [("approved", sym), ("previous", prev), ("alias", alias), ("withdrawn", wd)]
            self.ens = ens
        else:
            import csv
            ens = {}
            for line in open(os.path.join(t, "mgi", "MRK_ENSEMBL.rpt")):
                f = line.rstrip("\n").split("\t")
                if len(f) > 5:
                    ens[f[0]] = f[5]
            cur, wd = {}, collections.defaultdict(set)
            for r in csv.DictReader(open(os.path.join(t, "mgi", "MRK_List1.rpt")), delimiter="\t"):
                if r["Status"] == "O" and r["Marker Type"] in ("Gene", "Pseudogene"):
                    cur[r["Marker Symbol"]] = [r["MGI Accession ID"]]
                elif r["Status"] == "W" and r.get("Current MGI Accession ID (if withdrawn)"):
                    for x in r["Current MGI Accession ID (if withdrawn)"].split("|"):
                        wd[r["Marker Symbol"]].add(x.strip())
            syn = collections.defaultdict(set)
            for r in csv.DictReader(open(os.path.join(t, "mgi", "MRK_List2.rpt")), delimiter="\t"):
                for x in (r.get("Marker Synonyms (pipe-separated)") or "").split("|"):
                    if x:
                        syn[x].add(r["MGI Accession ID"])
            self.tables = [("approved", cur), ("withdrawn", wd), ("synonym", syn)]
            self.ens = ens

    def of(self, s, fold):
        for how, tab in self.tables:
            hit = tab.get(s)
            if hit is None and fold:
                continue
            if hit:
                hit = sorted(hit)
                if len(hit) == 1 or how in ("approved", "withdrawn"):
                    return sorted({self.ens[h] for h in hit if self.ens.get(h)}), how
                return [], f"{how}_of_several"
        return [], "none"


def build(org, publish=None, dry=False):
    from idtrack_reproducibility.families.astar.truth import TruthA
    items = lists(org)
    dates = release_dates()
    N = Names(org)
    tr = TruthA(org, TARGET)
    A = Authority(org)
    walks = {}

    def walk(g, r):
        if (g, r) not in walks:
            walks[(g, r)] = (sorted(tr.walk_fwd(g, r, TARGET, True, living=True)) if r < TARGET else [g],
                             sorted(tr.walk_fwd(g, r, TARGET, False, living=True, resurrect=True)) if r < TARGET else [g])
        return walks[(g, r)]

    grouped = {}
    for src, s, d, extra in items:
        R = TARGET if d is None else release_at(d, dates, N.rels)
        key = (src, s, R)
        if key in grouped:
            for k in ("n_rows", "n_sets"):
                if k in extra:
                    grouped[key][k] = grouped[key].get(k, 0) + extra[k]
            continue
        grouped[key] = dict(extra, list_date=str(d) if d else "rosetta (115)")
    out, counts = [], collections.Counter()
    for (src, s, R), extra in sorted(grouped.items()):
        fold = s == s.upper() and s.upper() != s.lower() and org == "mus_musculus"
        auth, how = A.of(s, False)
        if not auth and fold:
            # a mouse symbol written in capitals: the MGI symbol it is without letter case
            cands = [k for k in A.tables[0][1] if k.upper() == s]
            if len(cands) == 1:
                auth, how = A.of(cands[0], False)
                how = "approved_case"
        for variant in ("from", "default"):
            if variant == "from":
                read = R
            else:
                read = max((r for r in N.rels if N.genes(s, r)[0]), default=None)
                read = TARGET if N.genes(s, TARGET)[0] else read
            genes, match = N.genes(s, read) if read else ([], "none")
            q = dict(qid=f"d3|{src}|{s}|{R}|{variant}", family=FAMILY, kind="d3_symbol", stratum=f"d3_{src}_{variant}", unit=s.upper(),
                     id=s, from_release=R if variant == "from" else None, to_release=TARGET, final_database=None, strategy="both",
                     rule="set", cmp="base", broad=None, source=src, release_of_list=R, read_at=read, genes_at_read=genes, name_match=match,
                     authority=auth, authority_how=how, **extra)
            if not genes:
                q["expect"] = []
                q["open"] = f"no gene carries the symbol as its Ensembl name at {read if read else 'any release of the window'}"
                cls = "not_a_name"
            else:
                E, B = set(), set()
                for g in genes:
                    e, b = walk(g, read)
                    E |= set(e)
                    B |= set(b)
                E = {C.base_id(x) for x in E}
                at115 = set(N.genes(s, TARGET)[0])
                if E:
                    q["expect"] = sorted(E)
                    if B - E:
                        q["broad"] = sorted({C.base_id(x) for x in B | E})
                    cls = "one" if len(E) == 1 else "several"
                else:
                    q["expect"] = []
                    q["open"] = f"Ensembl's history ends for the genes the symbol names at {read}"
                    q["open_within"] = sorted(at115 | {C.base_id(x) for x in B})
                    cls = "history_ends"
            q["truth_class"] = cls
            counts[(q["stratum"], cls, match)] += 1
            out.append(q)
    for k, v in sorted(counts.items()):
        print(f"  {k[0]:32s} {k[1]:14s} {k[2]:6s} {v}")
    print(f"{org}: {len(out)} lines")
    if publish and not dry:
        L.publish(FAMILY, org, publish, out, description="D3: marker lists (PanglaoDB, CellMarker 2.0, MSigDB c8/m8, rosetta's lineage "
                  "table), each symbol to 115 from the release of the list's date and by default; Ensembl's names and history as "
                  "truth, HGNC/MGI beside it (gen_d3.py docstring)", generator=os.path.abspath(__file__),
                  sources=["gtf (genes tables)", "core (stable_id_event 115)", "hgnc/hgnc_complete_set.txt", "hgnc/withdrawn.txt",
                           "mgi/MRK_List1.rpt", "mgi/MRK_List2.rpt", "mgi/MRK_ENSEMBL.rpt", "<families>/d/inputs (SOURCES.tsv)"])
        print("published", publish)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--publish")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    build(a.organism, a.publish, a.dry)
