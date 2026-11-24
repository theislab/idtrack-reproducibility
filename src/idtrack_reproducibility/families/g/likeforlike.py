"""The rule of equal truth (docs/families/g.md, "Like for like"): which of G's lines asked at release 115 have the
same truth at another release T, the release a tool answers at (116: Ensembl REST, MyGene.info, HGNC, UniProt; 114:
g:Profiler; 111: BridgeDb; 113: EnsDb), so that a tool answering at T is compared with idtrack's graphs of 115 on
them only.

Usage: python -m idtrack_reproducibility.families.g.likeforlike --organism <o> --releases 116,114,111 [--draft <folder>]

A line asked to 115 keeps its truth at T when, read from the truth store (Ensembl's GTFs and core tables of every
release, the core cross-references of 87 on; nothing of idtrack), with the harness's `make_queries.Truth`:

1. every gene its truth names (`expect`, `broad`, `open_within`) is in the GTFs of 115 and of T, and no mapping
   session between them maps it to another id or retires it;
2. its question resolves to the same genes at T as at 115:
   - an Ensembl id (tasks a, b): Ensembl's history walked from the line's release to T and to 115 (`walk_series`,
     both the strict and the wide walk) gives the same genes;
   - a symbol or name (task c): the genes whose GTF name or a core cross-reference's label is the string (compared
     without case) are the same at T and at 115; a previous symbol, alias or synonym (task d) also the genes of which
     it is a synonym (`core_synonyms`), which the truth store must hold at both releases;
   - an external identifier (task e): the genes a core cross-reference of that accession (without its version) names
     are the same at T and at 115.

A line asked at its own release (to_release not 115) is asked only of a tool's instance at that release and needs no
rule. Writes `<families>/g/likeforlike/<organism>/release_<T>.parquet` (qid, release, equal, reason) and
`summary.json` (counts by reason, the truth files read with their md5 from SOURCES.tsv).
"""
import argparse
import collections
import gzip
import json
import os

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.g import fetch as F
from idtrack_reproducibility.harness import harness_lib as H

TARGET = 115
EQUAL = "equal"


def out_dir(organism):
    return os.path.join(C.family_dir("g"), "likeforlike", organism)


class Release:
    """The name and cross-reference tables of one release, as the rule reads them."""

    def __init__(self, organism, rel):
        self.org, self.rel = organism, rel
        self.files = []
        g = H.gtf_genes(organism, rel) or {}
        self.genes = set(g)
        self.names = collections.defaultdict(set)
        for gid, row in g.items():
            if row.get("gene_name"):
                self.names[row["gene_name"].casefold()].add(gid)
        self.files.append(f"gtf/genes_{rel}.tsv.gz")
        self.acc = collections.defaultdict(set)
        p = os.path.join(H.truth_dir(organism), "core_xref", f"{rel}.tsv.gz")
        self.has_xref = os.path.exists(p)
        if self.has_xref:
            self.files.append(f"core_xref/{rel}.tsv.gz")
            for r in H.read_tsv(p):
                gid = C.base_id(r["gene_id"])
                if r.get("acc"):
                    self.acc[r["acc"].split(".")[0] if r["db_name"].startswith("RefSeq") else r["acc"]].add(gid)
                if r.get("label"):
                    self.names[r["label"].casefold()].add(gid)
        p = os.path.join(H.truth_dir(organism), "core_synonyms", f"{rel}.tsv.gz")
        self.has_synonyms = os.path.exists(p)
        self.synonyms = collections.defaultdict(set)
        if self.has_synonyms:
            self.files.append(f"core_synonyms/{rel}.tsv.gz")
            for r in H.read_tsv(p):
                self.synonyms[r["synonym"].casefold()].add(C.base_id(r["gene_id"]))


def stable(truth, gid, lo, hi):
    """The gene is in the GTFs of lo and hi and no session between them maps it to another id or retires it."""
    pr = truth.present.get(gid, ())
    if lo not in pr or hi not in pr:
        return False
    for s in truth.sessions_between(lo, hi):
        if gid in truth.succ[s] and truth.succ[s][gid] - {gid}:
            return False
        if gid in truth.succ[s] and not truth.succ[s][gid]:
            return False
    return True


def judge(line, truth, at115, atT, T):
    """(equal, reason) of one line asked to 115, at release T."""
    lo, hi = min(T, TARGET), max(T, TARGET)
    named = set()
    for k in ("expect", "broad", "open_within"):
        named |= {C.base_id(x) for x in (line.get(k) or [])}
    moved = sorted(g for g in named if not stable(truth, g, lo, hi))
    if moved:
        return False, f"a gene of the truth is not the same at {T}: {','.join(moved[:3])}"
    t = line.get("task")
    v = line["id"]
    if t in ("a", "b"):
        gid = C.base_id(v)
        r = line.get("from_release")
        if r is None or r not in truth.genes:
            return False, "no release of the question in the truth's window"
        for strict in (True, False):
            w = truth.walk_series(gid, r, [TARGET, T], strict)
            if w.get(TARGET, set()) != w.get(T, set()):
                return False, f"Ensembl's history from {r} reaches other genes at {T} than at 115"
        return True, EQUAL
    if t in ("c", "d"):
        key = v.casefold()
        if at115.names.get(key, set()) != atT.names.get(key, set()):
            return False, f"the string names other genes at {T} than at 115"
        if t == "d":
            if not (at115.has_synonyms and atT.has_synonyms):
                return False, f"the truth store has no synonyms of release {T if not atT.has_synonyms else TARGET}"
            if at115.synonyms.get(key, set()) != atT.synonyms.get(key, set()):
                return False, f"the string is a synonym of other genes at {T} than at 115"
        return True, EQUAL
    if t == "e":
        if not (at115.has_xref and atT.has_xref):
            return False, "the truth store has no core cross-references of both releases"
        a = v.split(".")[0] if v[:3] in ("NM_", "XM_", "NR_", "XR_", "NP_", "XP_") else v
        if at115.acc.get(a, set()) != atT.acc.get(a, set()):
            return False, f"the accession names other genes at {T} than at 115"
        return True, EQUAL
    return False, f"no rule for task {t}"


def run(organism, releases, draft=None):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from idtrack_reproducibility.harness import make_queries as M
    lines = [q for q in F.g_lines(organism, "local", None, draft) if q.get("to_release") == TARGET]
    truth = M.Truth(organism, max([TARGET] + list(releases)))
    at115 = Release(organism, TARGET)
    os.makedirs(out_dir(organism), exist_ok=True)
    srcs = {r["what"]: r for r in H.read_tsv(os.path.join(H.truth_dir(organism), "SOURCES.tsv"))}
    summary = {"organism": organism, "lines": len(lines), "releases": {}}
    for T in releases:
        atT = Release(organism, T)
        rows, counts = [], collections.Counter()
        for q in lines:
            eq, why = judge(q, truth, at115, atT, T)
            rows.append(dict(qid=q["qid"], release=T, equal=eq, reason=why, tier=q.get("tier"), task=q.get("task"),
                             stratum=q.get("stratum")))
            counts[(q.get("tier"), q.get("task"), why.split(":")[0])] += 1
        path = os.path.join(out_dir(organism), f"release_{T}.parquet")
        pq.write_table(pa.Table.from_pylist(rows), path + ".tmp")
        os.replace(path + ".tmp", path)
        files = sorted(set(at115.files + atT.files))
        summary["releases"][str(T)] = dict(
            counts=[dict(tier=k[0], task=k[1], reason=k[2], n=n) for k, n in sorted(counts.items())],
            truth_files=[dict(file=f, md5=(srcs.get(f.rsplit(".tsv.gz", 1)[0]) or {}).get("md5")) for f in files])
    json.dump(summary, open(os.path.join(out_dir(organism), "summary.json"), "w"), indent=1)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.likeforlike", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", required=True)
    ap.add_argument("--releases", default="116,114,111")
    ap.add_argument("--draft", default=None)
    a = ap.parse_args(argv)
    s = run(a.organism, [int(x) for x in a.releases.split(",")], a.draft)
    for T, d in s["releases"].items():
        for c in d["counts"]:
            print(a.organism, T, c["tier"], c["task"], c["n"], c["reason"])


if __name__ == "__main__":
    main()
