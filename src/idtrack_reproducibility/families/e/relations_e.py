"""E8's relations that need two calls in a row, and the explained path against Ensembl's tables.

Usage (the tree's Python):
  relations_e.py run --organism <o> --clone <clone> --graph <graph> --out <run>/<organism> [--count 3000]
  relations_e.py submit ... (one job)

Relations (no truth; each violation listed with the answers):

- composition: R1 -> R2 answers one gene h (`best`); then h R2 -> R3 must give what R1 -> R3 gives. A violation is
  read leg by leg against Ensembl's history (`make_queries.Truth`, strict walk with living ids): where the history
  itself composes for that triple (the walk R1 -> R2 is {h} and the walk of h R2 -> R3 equals the walk R1 -> R3) the
  violation is idtrack's (`refuted`); otherwise the history does not compose there either (splits and merges:
  `by_history`). The travel matrix's reading (accuracy README): a report, not a gate.
- return: an id that the history carries unchanged from R1 to R2 (same id and version in both GTFs, the strict walk
  {id}) answers itself R1 -> R2, and its answer back R2 -> R1 is itself.
- final database and back: a gene of 115 converted to a database X (`final_database`), and every id it gets read
  back (`from_database=X`, view `all`) reaches the gene.

The path against Ensembl's tables (trees of the convert generation: the path is data there, `PathStep`):
for history lines of the harness's published sets (kinds retired_*, reverse_*, unchanged_old_to_current,
cross_assembly, travel_*; a seeded sample), `explain=True`, every candidate of `all`:
- every `history` and `self_loop` step is an event of the release-115 core `stable_id_event` between its two ids
  (in either order) of a session crossed between its `old_release` and `new_release` (by the session's labels or by
  where its effect appears, `Truth.effective_release`); the versions of the event are reported apart;
- every `bridge_near`/`bridge_far` step is an identifier that the BioMart tables of its release attach to the gene
  (any table: names, synonyms, NCBI, RefSeq, UniProt...), where the store holds BioMart tables of that release, else
  undetermined; likewise the `final` step at the target release;
- the path is contiguous, starts at the query's node and ends at the candidate's gene.

Writes `<out>/e8/<relation>.jsonl.gz` (violations) and `<out>/e8_report.tsv`.
"""
import argparse
import collections
import gzip
import json
import os
import random
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L

SEED = 20261118
HISTORY_KINDS = ("retired_none", "retired_merged", "retired_renamed", "retired_split", "reverse_present", "reverse_absent",
                 "unchanged_old_to_current", "cross_assembly", "travel_born_late", "travel_handoff", "travel_merged",
                 "travel_renamed", "travel_retired_none", "travel_split", "travel_unchanged")


def base(x):
    return C.base_id(x)


class Ask:
    """One call through family E's adapter (best and all views)."""

    def __init__(self, a):
        from idtrack_reproducibility.families.e import adapter_e
        self.ad = adapter_e.AdapterE(a.clone, a.graph, a.organism, C.release_of_graph(a.graph), a.out, timeout=60)
        self.n = 0

    def __call__(self, ident, view="best", **kw):
        self.n += 1
        line = dict(qid=f"rel{self.n}", id=ident, strategy=view, **{k: v for k, v in kw.items() if v is not None})
        a = self.ad.ask(line)[0]
        if a["outcome"] != "answer":
            return None
        return set(a["got"])


def composition(ask, T, org, n, rng, out):
    rels = [r for r in T.rels]
    ids = sorted(T.present)
    stats = collections.Counter()
    with gzip.open(os.path.join(out, "e8", "composition.jsonl.gz"), "wt") as fo:
        tries = 0
        while stats["triples"] < n and tries < 50 * n:
            tries += 1
            g = rng.choice(ids)
            pr = sorted(T.present[g])
            r1 = rng.choice(pr)
            r2, r3 = rng.sample(rels, 2)
            if len({r1, r2, r3}) < 3:
                continue
            gv = T.vid(g, r1)
            a12 = ask(gv, from_release=r1, to_release=r2)
            if not a12 or len(a12) != 1:
                stats["not_one_gene"] += 1
                continue
            (h,) = a12
            stats["triples"] += 1
            b = ask(h, from_release=r2, to_release=r3)
            d = ask(gv, from_release=r1, to_release=r3)
            if b == d:
                stats["holds"] += 1
                continue
            w12 = T.walk_series(g, r1, [r2], True, living=True, resurrect=True)[r2]
            w13 = T.walk_series(g, r1, [r3], True, living=True, resurrect=True)[r3]
            w23 = T.walk_series(base(h), r2, [r3], True, living=True, resurrect=True)[r3]
            composes = w12 == {base(h)} and w23 == w13
            k = "refuted" if composes else "by_history"
            stats[k] += 1
            fo.write(json.dumps(dict(id=gv, r1=r1, r2=r2, r3=r3, h=h, h_to_r3=sorted(b or []), direct=sorted(d or []),
                                     walk_12=sorted(w12), walk_13=sorted(w13), walk_h3=sorted(w23), reading=k)) + "\n")
    return stats


def returns(ask, T, org, n, rng, out):
    stats = collections.Counter()
    ids = sorted(T.present)
    with gzip.open(os.path.join(out, "e8", "return.jsonl.gz"), "wt") as fo:
        tries = 0
        while stats["pairs"] < n and tries < 50 * n:
            tries += 1
            g = rng.choice(ids)
            pr = sorted(T.present[g])
            if len(pr) < 2:
                continue
            r1, r2 = sorted(rng.sample(pr, 2))
            if rng.random() < 0.5:
                r1, r2 = r2, r1
            if T.ver(g, r1) != T.ver(g, r2):
                continue
            w = T.walk_series(g, r1, [r2], True)[r2]
            if w != {g}:
                continue
            v1, v2 = T.vid(g, r1), T.vid(g, r2)
            fwd = ask(v1, from_release=r1, to_release=r2)
            stats["pairs"] += 1
            cmp_ = (lambda s: {base(x) for x in s}) if "." not in v1 else (lambda s: s)
            if fwd is None or cmp_(fwd) != cmp_({v2}):
                stats["forward_not_itself"] += 1
                fo.write(json.dumps(dict(check="forward", id=v1, r1=r1, r2=r2, got=sorted(fwd or []))) + "\n")
                continue
            back = ask(v2, from_release=r2, to_release=r1)
            if back is None or cmp_(back) != cmp_({v1}):
                stats["back_not_itself"] += 1
                fo.write(json.dumps(dict(check="back", id=v2, r1=r2, r2=r1, got=sorted(back or []))) + "\n")
            else:
                stats["holds"] += 1
    return stats


def roundtrip(ask, org, n, rng, out, dbs):
    stats = collections.Counter()
    g115 = H.gtf_genes(org, 115, "primary") or {}
    genes = sorted(f"{g}.{r['gene_version']}" for g, r in g115.items())
    with gzip.open(os.path.join(out, "e8", "roundtrip.jsonl.gz"), "wt") as fo:
        for gv in rng.sample(genes, min(n, len(genes))):
            for db in dbs:
                ids = ask(gv, view="all", final_database=db)
                if not ids:
                    stats[f"{db}|no_id"] += 1
                    continue
                for x in sorted(ids)[:3]:
                    stats[f"{db}|ids"] += 1
                    back = ask(x, view="all", from_database=db)
                    if back and gv in back:
                        stats[f"{db}|holds"] += 1
                    else:
                        stats[f"{db}|violation"] += 1
                        fo.write(json.dumps(dict(gene=gv, database=db, id=x, back=sorted(back or [])[:20])) + "\n")
    return stats


class Tables:
    """Ensembl's events (release-115 core) and BioMart's identifiers per release, for the path check."""

    def __init__(self, org, T):
        self.T = T
        _, ev, sessions = H.core_tables(org)
        self.ev = collections.defaultdict(list)  # frozenset of the two base ids (or one) -> [(old, ov, new, nv, session)]
        for o, ov, n, nv, s, sc in ev:
            if o and n:
                self.ev[frozenset((o, n))].append((o, ov, n, nv, s))
        self.sessions = sessions
        self.bm = {}

    def event(self, a, b, old_release, new_release):
        cands = self.ev.get(frozenset((base(a), base(b))), [])
        if not cands:
            return "no_event", None
        for o, ov, n, nv, s in cands:
            so, sn = H.rel_int(self.sessions[s]["old_release"]), H.rel_int(self.sessions[s]["new_release"])
            eff = self.T.eff.get(s)
            if (old_release, new_release) == (so, sn) or (old_release is not None and new_release is not None and eff is not None
                                                          and old_release < eff <= new_release):
                return "event", (o, ov, n, nv, s)
        return "event_of_another_session", cands[0]


    def tx_held(self, gene, ident, release):
        """A transcript or translation id shared through a bridge: the GTF of the release lists it under the gene."""
        key = ("tx", release)
        if key not in self.bm:
            rows = H.gtf_tx(self.T.org, release)
            d = collections.defaultdict(set) if rows is not None else None
            for r in rows or ():
                d[r["gene_id"]].add(r["transcript_id"])
                if r.get("protein_id"):
                    d[r["gene_id"]].add(r["protein_id"])
            self.bm[key] = d
        d = self.bm[key]
        if d is None:
            return None
        return base(ident) in d.get(base(gene), ())

    def tx_pair(self, a, b, release):
        key = ("pairs", release)
        if key not in self.bm:
            rows = H.gtf_tx(self.T.org, release)
            self.bm[key] = None if rows is None else {frozenset((r["transcript_id"], r["protein_id"])) for r in rows if r.get("protein_id")}
        d = self.bm[key]
        return None if d is None else frozenset((base(a), base(b))) in d

    def held(self, gene, ident, release, databases=()):
        """True / False when Ensembl's tables of the release attach the identifier to the gene; None when they cannot
        say: transcripts and translations through the GTF of the release; a synonym only when a table holds it (BioMart
        exports only the display name's synonyms: a missing synonym of another source is undetermined); other ids through
        the BioMart tables of the release, where the store has them."""
        if any(d in ("ensembl_transcript", "ensembl_translation") for d in databases) or C.ENS_ID.match(str(ident)):
            return self.tx_held(gene, ident, release)
        syn = ident.startswith("synonym_id::") or any(d.startswith("synonym_id::") for d in databases)
        h = self._bm_held(gene, ident, release)
        return None if (syn and h is False) else h

    def _bm_held(self, gene, ident, release):
        if release not in self.bm:
            if release not in set(H.biomart_releases(self.T.org)):
                self.bm[release] = None
            else:
                d = collections.defaultdict(set)
                root = os.path.join(C.TRUTH, self.T.org, "biomart", str(release))
                for f in os.listdir(root):
                    if not f.endswith(".tsv") or f == "transcripts.tsv":
                        continue
                    t = H.biomart(self.T.org, release, f[:-4])
                    if t is None:
                        continue
                    for row in t[1]:
                        for v in row[1:]:
                            if v:
                                d[row[0]].add(v)
                self.bm[release] = d
        d = self.bm[release]
        if d is None:
            return None
        ident = ident.split("::", 1)[1] if ident.startswith("synonym_id::") else ident
        return ident in d.get(base(gene), ())


def paths(ad, T, org, n, rng, out):
    stats = collections.Counter()
    tabs = Tables(org, T)
    hl, _, _ = L.load("harness", org)
    pool = [q for q in hl if q["kind"] in HISTORY_KINDS and C.ENS_ID.match(str(q["id"]))]
    pool.sort(key=lambda q: q["qid"])
    sample = rng.sample(pool, min(n, len(pool)))
    with gzip.open(os.path.join(out, "e8", "path.jsonl.gz"), "wt") as fo:
        for q in sample:
            kw = {k: q.get(k) for k in ("from_release", "to_release", "final_database") if q.get(k) is not None}
            try:
                r = ad.api.convert(q["id"], explain=True, budget_seconds=60, **kw)
            except Exception as e:  # noqa: BLE001
                stats["call_error"] += 1
                fo.write(json.dumps(dict(qid=q["qid"], problem="call_error", detail=str(e)[:300])) + "\n")
                continue
            stats["lines"] += 1
            steps = r.path.steps if r.path else ()
            for c in r.candidates:
                stats["candidates"] += 1
                chain = [steps[i] for i in c.path]
                probs = []
                if chain:
                    if chain[-1].target != c.gene:
                        probs.append(("end_not_gene", chain[-1].target))
                    if chain[0].kind != "identity" and chain[0].source != r.graph_id:
                        probs.append(("start_not_query", chain[0].source))
                    for s1, s2 in zip(chain, chain[1:]):
                        if s1.target != s2.source:
                            probs.append(("not_contiguous", f"{s1.target} -> {s2.source}"))
                            break
                for s in chain + [steps[i] for i in c.final_path]:
                    if s.kind in ("history", "self_loop"):
                        stats["history_steps"] += 1
                        how, ev = tabs.event(s.source, s.target, s.old_release, s.new_release)
                        stats[f"history_{how}"] += 1
                        if how != "event":
                            probs.append((f"history_{how}", [s.source, s.target, s.old_release, s.new_release, ev]))
                        else:
                            o, ov, nn, nv, sess = ev
                            got = {s.source, s.target}
                            want = {f"{o}.{ov}", f"{nn}.{nv}"}
                            if ov and nv and all(C.ENS_VERSIONED.match(x) for x in got):
                                k = "history_event_versions_equal" if got == want else "history_event_versions_differ"
                                stats[k] += 1
                                if got != want:
                                    fo.write(json.dumps(dict(qid=q["qid"], gene=c.gene, note="versions",
                                                             step=[s.source, s.target, s.old_release, s.new_release],
                                                             event=[o, ov, nn, nv, sess])) + "\n")
                    elif s.kind in ("bridge_near", "bridge_far", "final"):
                        stats[f"{s.kind}_steps"] += 1
                        is_gene = re.compile(r"^ENS[A-Z]*G\d{11}")
                        is_tx = re.compile(r"^ENS[A-Z]*[TP]\d{11}")
                        src, tgt = str(s.source), str(s.target)
                        if is_gene.match(src) or is_gene.match(tgt):
                            gene, ext = (src, tgt) if is_gene.match(src) else (tgt, src)
                            h = tabs.held(gene, ext, s.release, tuple(s.databases)) if s.release is not None else None
                        elif is_tx.match(src) and is_tx.match(tgt):
                            gene, ext = src, tgt  # a transcript and its translation: one row of the GTF of the release
                            h = tabs.tx_pair(src, tgt, s.release) if s.release is not None else None
                        else:
                            gene, ext = src, tgt  # an external id on a transcript or translation: no table of the store says
                            h = None
                        k = "undetermined" if h is None else ("held" if h else "not_held")
                        stats[f"{s.kind}_{k}"] += 1
                        if h is False:
                            probs.append((f"{s.kind}_not_held", [gene, ext, s.release, list(s.databases)]))
                if probs:
                    stats["candidates_with_problem"] += 1
                    fo.write(json.dumps(dict(qid=q["qid"], id=q["id"], kw=kw, gene=c.gene, problems=probs), default=str) + "\n")
    return stats


def run(a):
    from idtrack_reproducibility.harness import make_queries as MQ
    os.makedirs(os.path.join(a.out, "e8"), exist_ok=True)
    rng = random.Random(f"{SEED}:{a.organism}:e8")
    T = MQ.Truth(a.organism, 115)
    ask = Ask(a)
    t0 = time.time()
    res = {}
    res["composition"] = composition(ask, T, a.organism, a.n, rng, a.out)
    res["return"] = returns(ask, T, a.organism, a.n, rng, a.out)
    gen = C.generation(a.clone)
    fdb = ask.ad.api.final_databases() if gen == "convert" else ask.ad.api.valid_final_databases()
    want = []
    for pat in ("HGNC Symbol", "MGI Symbol", "VGNC Symbol", "UniProtKB Gene Name", "RefSeq_mRNA", "Uniprot/SWISSPROT", "NCBI gene",
                "NCBI gene (formerly Entrezgene)", "EntrezGene"):
        if pat in fdb:
            want.append(pat)
    res["roundtrip"] = roundtrip(ask, a.organism, max(1, a.n // 3), rng, a.out, want)
    if gen == "convert":
        res["path"] = paths(ask.ad, T, a.organism, a.n, rng, a.out)
    with open(os.path.join(a.out, "e8_report.tsv"), "w") as fo:
        fo.write("relation\torganism\tcount\tvalue\n")
        for k, st in res.items():
            for kk, v in sorted(st.items()):
                fo.write(f"{k}\t{a.organism}\t{kk}\t{v}\n")
    json.dump(dict(organism=a.organism, seconds=round(time.time() - t0, 1), calls=ask.n, clone=a.clone, commit=C.git_commit(a.clone)[0],
                   databases=want, generation=gen), open(os.path.join(a.out, "e8", "info.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["run", "submit"])
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--count", dest="n", type=int, default=3000)
    a = ap.parse_args()
    a.out, a.clone, a.graph = (os.path.abspath(x) for x in (a.out, a.clone, a.graph))
    if a.what == "run":
        if os.environ.get("PYTHONHASHSEED") != "0":
            os.environ["PYTHONHASHSEED"] = "0"
            os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
        return run(a)
    logs = os.path.join(a.out, "e8", "logs")
    os.makedirs(logs, exist_ok=True)
    from idtrack_reproducibility import jobs
    mem = jobs.memory_for(a.organism) + 16 if C.generation(a.clone) == "d4c2dda" else 24
    print(jobs.submit(f"e8_{a.organism}", jobs.module(C.python_for(a.clone), "idtrack_reproducibility.families.e.relations_e",
                                                      f"run --organism {a.organism} --clone {a.clone} --graph {a.graph} "
                                                      f"--out {a.out} --count {a.n}"),
                      logs, cpus=2, mem_gb=mem, time="12:00:00"))


if __name__ == "__main__":
    main()
