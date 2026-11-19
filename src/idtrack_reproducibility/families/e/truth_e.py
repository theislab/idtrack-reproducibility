"""Family E's reading of the one truth store (nothing of idtrack): what a string can name, per organism.

- `Names(org)`: every string the truth store knows for the organism, as Ensembl's own tables write it: the GTF's
  `gene_id` and `gene_name` of every release, and every column but the gene id of every BioMart table of every
  release the store holds (display names, synonyms, HGNC/MGI/VGNC symbols and ids, UniProt gene names, WikiGene
  names, NCBI ids and labels, RefSeq, UniProt, RFAM accessions), plus VGNC's own symbols (`harness_lib.vgnc_symbols`).
  `fold(s)` is the reading idtrack's case and separator steps make equal (casefold; `-`, `_`, `.` as one character);
  `owners_folded(s)` the base gene ids of every string that folds to it; `exact115(s)` the genes that carry exactly
  `s` at 115 (GTF 115 names and ids, BioMart 115). Used to keep a form out of a judgement when it is itself another
  gene's identifier (`open`), and to judge cross-organism names.
- `Current(org)`: the genes of GTF 115 (id, version, name, region), transcripts and translations with their gene.

Cached per organism in `<cache>/e/<organism>.pkl`, keyed by the store's `SOURCES.tsv` size and time.
"""
import collections
import gzip
import os
import pickle
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H

SEP = re.compile(r"[-_.]")


def fold(s):
    return SEP.sub("-", s.casefold())


class Names:
    def __init__(self, org):
        self.org = org
        cache = os.path.join(C.CACHE, "e", f"{org}.pkl")
        src = os.path.join(C.TRUTH, org, "SOURCES.tsv")
        key = (os.path.getsize(src), int(os.path.getmtime(src)), 3)
        if os.path.exists(cache):
            d = pickle.load(open(cache, "rb"))
            if d.get("key") == key:
                self.__dict__.update(d)
                return
        self.key = key
        folded = collections.defaultdict(set)   # fold(string) -> base gene ids
        exact115 = collections.defaultdict(set)  # string at 115 -> versioned gene ids (GTF 115 version)
        display115 = collections.defaultdict(set)  # BioMart 115 display name -> versioned ids
        g115 = H.gtf_genes(org, 115) or {}
        ver = {g: f"{g}.{r['gene_version']}" if r.get("gene_version") else g for g, r in g115.items()}
        for rel in H.gtf_releases(org):
            for r in H.read_tsv(os.path.join(C.TRUTH, org, "gtf", f"genes_{rel}.tsv.gz")):
                for s in (r["gene_id"], r.get("gene_name") or ""):
                    if s:
                        folded[fold(s)].add(r["gene_id"])
                        if rel == 115:
                            exact115[s].add(ver.get(r["gene_id"], r["gene_id"]))
        for rel in H.biomart_releases(org):
            d = os.path.join(C.TRUTH, org, "biomart", str(rel))
            for f in sorted(os.listdir(d)):
                if not f.endswith(".tsv") or f == "transcripts.tsv":
                    continue
                t = H.biomart(org, rel, f[:-4])
                if t is None:
                    continue
                header, rows = t
                cols = [i for i, h in enumerate(header) if h not in ("ensembl_gene_id", "ensembl_gene_id_version",
                                                                        "external_gene_source", "gene_biotype")]
                for row in rows:
                    g = row[0]
                    for i in cols:
                        s = row[i] if i < len(row) else ""
                        if s:
                            folded[fold(s)].add(g)
                            if rel == 115:
                                exact115[s].add(ver.get(g, g))
                                if f == "genes.tsv" and header[i] in ("external_gene_name", "external_gene_id"):
                                    display115[s].add(ver.get(g, g))
            vg = H.vgnc_symbols(org, rel)
            if vg:
                for s, gs in vg[1].items():
                    folded[fold(s)] |= gs
                    if rel == 115:
                        exact115[s] |= {ver.get(g, g) for g in gs}
        self.folded = dict(folded)
        self.exact115 = dict(exact115)
        self.display115 = dict(display115)
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        tmp = f"{cache}.{os.getpid()}.tmp"
        pickle.dump(self.__dict__, open(tmp, "wb"), protocol=4)
        os.replace(tmp, cache)

    def owners_folded(self, s):
        return self.folded.get(fold(s), set())

    def known(self, s):
        """True when the string, in any case and separator form, is an identifier the truth knows."""
        return fold(s) in self.folded


class Current:
    def __init__(self, org):
        g = H.gtf_genes(org, 115) or {}
        self.genes = {k: dict(version=r.get("gene_version") or "", name=r.get("gene_name") or "", region=r["region"])
                      for k, r in g.items()}
        self.versioned = {k: (f"{k}.{v['version']}" if v["version"] else k) for k, v in self.genes.items()}
        self.tx, self.tl = {}, {}
        for r in H.gtf_tx(org, 115) or []:
            gv = self.versioned.get(r["gene_id"], r["gene_id"])
            self.tx[r["transcript_id"]] = (r.get("transcript_version") or "", gv)
            if r.get("protein_id"):
                self.tl[r["protein_id"]] = (r.get("protein_version") or "", gv)


_ = gzip
