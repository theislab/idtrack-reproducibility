"""A*'s truth: Ensembl's history walked as the harness walks it (`make_queries.Truth`), with two corrections that the
hardest lineages need, and a second, independent reading from Ensembl REST `archive/id`.

Nothing here reads idtrack. Sources: the truth store (`accuracy/truth/<organism>/`): GTFs of every release, the
release-115 core history tables, the core gene tables fetched by `fetch_core_genes.py`, REST records fetched by
`fetch_rest.py`.

Sessions are crossed where their effect appears (`make_queries.Truth.effective_release`; corrected in the harness
on 2026-11-18 after this family found rhesus 50 -> 52 and 52 -> 58 placed at their old release). One addition over
`make_queries.Truth`: **genes of the core database that the GTF lacks are present.** Human 48 to 51 (genes on unplaced
scaffolds) and cynomolgus 112 (47 genes): the core's gene stable ids (Ensembl's FTP MySQL dumps) are added to the
presence, with the core's version and no place (`region = core_only`).
"""
import gzip
import json
import os
import sys

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M

CORE_GENES = {  # organism -> {release: file in truth/<organism>/core/} (fetch_core_genes.py)
    "homo_sapiens": {48: "homo_sapiens_core_48_36j.gene_stable_id.txt.gz", 49: "homo_sapiens_core_49_36k.gene_stable_id.txt.gz",
                     50: "homo_sapiens_core_50_36l.gene_stable_id.txt.gz", 51: "homo_sapiens_core_51_36m.gene_stable_id.txt.gz"},
    "macaca_fascicularis": {112: "macaca_fascicularis_core_112_6.gene.txt.gz"},
}


def core_gene_ids(org, rel):
    """{stable id: version} of a core gene table dump (gene_stable_id: columns 1, 2; gene: columns 12, 13)."""
    f = CORE_GENES.get(org, {}).get(rel)
    if not f:
        return None
    out = {}
    with gzip.open(os.path.join(H.truth_dir(org), "core", f), "rt") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            sid, ver = (c[1], c[2]) if f.endswith("gene_stable_id.txt.gz") else (c[12], c[13])
            out[sid] = ver
    return out


class TruthA(M.Truth):
    """`make_queries.Truth` with the core genes the GTF lacks, and with rule `no_mito_nuclear_step` (decided): a step
    of Ensembl's history between the mitochondrial genome and a nuclear sequence is not followed (`mito_cut`; the sets
    published before the rule are revised by `revise_rules.py` and family A's `revise_a.py`)."""

    def __init__(self, org, target, mito_cut=True):
        super().__init__(org, target, mito_cut=mito_cut)
        self.core_only = {}
        for rel in CORE_GENES.get(org, {}):
            if rel not in self.genes:
                continue
            c = core_gene_ids(org, rel)
            extra = set(c) - set(self.genes[rel])
            # a GTF that writes no versions (human before 77) is compared without versions (`cmp` base): so are the
            # core's genes added to it
            versions = any(g["gene_version"] for g in self.genes[rel].values())
            for x in extra:
                self.genes[rel][x] = dict(gene_id=x, gene_version=c[x] if versions else "", gene_name="", gene_source="", gene_biotype="",
                                          seqname="", start="", end="", strand="", region="core_only")
                self.present[x].add(rel)
                self.ever.add(x)
            self.core_only[rel] = extra
        self._sides = {}

    def locus(self, gid, rel):
        g = self.genes.get(rel, {}).get(gid)
        if not g or g["region"] == "core_only":
            return None
        return super().locus(gid, rel)


class Rest:
    """Ensembl REST `archive/id` records (`truth/<organism>/rest/archive_astar.json`): REST's reading of where an id
    went, independent of the core tables and of the walk above."""

    def __init__(self, org, truth):
        p = os.path.join(H.truth_dir(org), "rest", "archive_astar.json")
        self.arch = json.load(open(p)) if os.path.exists(p) else {}
        self.meta = self.arch.pop("_meta", {})
        # where REST did not answer (rest.ensembl.org unreachable from 2026-11-18 ~15:30): the same records read from
        # the 116 core's tables (rest_from_core.py); `source` of a record says which
        p2 = os.path.join(H.truth_dir(org), "rest", "archive_astar_core116.json")
        if os.path.exists(p2):
            extra = json.load(open(p2))
            self.meta["core116"] = extra.pop("_meta", {})
            for k, v in extra.items():
                if k not in self.arch or "_error" in self.arch[k]:
                    self.arch[k] = v
        self.tr = truth
        self._memo = {}

    def last(self, x):
        rec = self.arch.get(x)
        if rec is None or "_error" in rec:
            return None
        return H.rel_int(rec.get("release"))

    def sources(self, x, t):
        """The sources of the records `at(x, t)` read: {'rest'}, {'core116'} or both."""
        out, seen, stack = set(), set(), [x]
        while stack:
            y = stack.pop()
            if y in seen or y not in self.arch or "_error" in self.arch[y]:
                continue
            seen.add(y)
            out.add(self.arch[y].get("source", "rest"))
            if (self.last(y) or 0) < t:
                stack.extend(p["stable_id"] for p in self.arch[y].get("possible_replacement") or [])
        return out

    def at(self, x, t):
        """(ids at release t that x is or became by REST, complete?): an id whose latest release (REST `release`)
        is at or after t stands (kept when the GTF/core has it at t); an earlier one is replaced by its
        `possible_replacement` ids, followed until they stand. For forward lines only (REST has no predecessors)."""
        key = (x, t)
        if key in self._memo:
            return self._memo[key]
        out, seen, stack, complete = set(), set(), [x], True
        while stack:
            y = stack.pop()
            if y in seen:
                continue
            seen.add(y)
            last = self.last(y)
            if last is None:
                complete = False
                continue
            if last >= t:
                if t in self.tr.present.get(y, ()):
                    out.add(y)
                continue
            for p in self.arch[y].get("possible_replacement") or []:
                stack.append(p["stable_id"])
        self._memo[key] = (out, complete)
        return out, complete
