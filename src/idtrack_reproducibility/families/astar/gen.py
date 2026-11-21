"""A*: the query set of the hardest lineages of one organism (docs/methods.md, A*), published once, never rewritten.

Usage: gen.py <organism> --set <name> [--target 115] [--ids <file>] [--limit-ids N] [--dry-run]

For every selected id (`selection/<organism>.tsv.gz`, ids present in a GTF of the window; or the ids of `--ids`):
- `astar_from`: the stable id, from every release R where it exists to every release t of the window (forward,
  backward, in place);
- `astar_from_v`: the same with its version at R (where the GTF of R writes versions);
- `astar_default`, `astar_default_v`: without `from_release`, the stable id and each of its versions, to every t
  (rule `start_release`: read at t if it exists there (that version, for a versioned id), else at the latest release where it
  exists (with that version));
- `astar_db`: for every sixth id (crc32 % 6 == 0) the stable id from every R to every t with a BioMart table, to the
  name database of the travel matrix (HGNC Symbol, MGI Symbol, else UniProtKB Gene Name).
Every line has `strategy` `both` (scored under `best` and `all`).

Truth (truth.py; nothing of idtrack): the strict walk with living ids (`expect`), the broad walk with living and
returned ids (`broad`), both without a step between the mitochondrial genome and a nuclear sequence (rule
`no_mito_nuclear_step`; sets published before it are revised by revise_rules.py), copies on patches reported apart (rule `judge_copies_apart`), the locus rule with the place mapped across
assemblies (the corrected judge's `places`); BioMart's values of the expected genes for `astar_db`. The second reading: for
a line whose answer is read forward (t after the release read), REST `archive/id` followed through
`possible_replacement` (truth.Rest.at); where it is complete and gives other genes than the walk, the line carries
`open` (the reason) and `rest` (REST's genes), and is undetermined (docs/methods.md).

`stratum`: the direction (`same`, `forward`, `backward`, each `_handoff` when the two releases are on different
assemblies), `default`, `db`; `unit`: the id's lineage (the component's smallest id; the id itself when it has no
cross-id event). The patterns of an id are in the selection table (the report joins them).
"""
import argparse
import collections
import datetime
import json
import os
import sys
import time
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.astar import truth as TA
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M

FAMILY = "astar"
DATA = os.path.join(C.BENCH, FAMILY)
NAME_DB = {"homo_sapiens": "HGNC Symbol", "mus_musculus": "MGI Symbol"}  # others: UniProtKB Gene Name (travel matrix)
CHUNK = 400_000
GEN_SHA256 = L.sha256(os.path.abspath(__file__))  # of the code that runs (read before any later edit of the file)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


class Gen:
    def __init__(self, org, T, mito_cut=True):
        self.org, self.T = org, T
        self.tr = TA.TruthA(org, T, mito_cut=mito_cut)
        self.rest = TA.Rest(org, self.tr)
        self.ctx = M.Ctx78(org, T, truth=self.tr)
        self.asm = H.gtf_assemblies(org)
        self.db = NAME_DB.get(org, "UniProtKB Gene Name")
        self.grp, self.col = M.GENE_DB[self.db]
        self.bm_rels = [r for r in H.biomart_releases(org) if r in self.tr.rels]
        self.regions = {}
        self.mapper = None
        self.map_cache = {}
        from idtrack_reproducibility.harness import assembly_map
        self.exc = assembly_map.Exceptions(org)
        self.counts = collections.Counter()

    def apart_of(self, Ev, rel):
        """Rule `judge_copies_apart`'s copies apart, at every release: a gene on a haplotype, patch or alternative locus (GTF region
        `patch`, or before 82, where the GTF labels them `primary`, a sequence of the release's assembly exceptions:
        the regions hunt's J3, decided for idtrack's next version) is a copy."""
        if rel not in self.regions:
            exc = self.exc.rows(rel)
            where = {g: ("patch" if r["region"] == "patch" or (r["seqname"] in exc and any(e["type"] != "PAR" for e in
                                                                                               exc[r["seqname"]]))
                         else "primary") for g, r in (self.tr.genes.get(rel) or {}).items()}
            self.regions[rel] = where if "patch" in where.values() else {}
        where = self.regions[rel]
        if not where or len(Ev) < 2:
            return Ev, None
        sides = collections.defaultdict(list)
        for gene in Ev:
            sides[where.get(M.gene_base(gene))].append(gene)
        if None in sides or not sides["primary"] or not sides["patch"]:
            return Ev, None
        return sorted(sides["primary"]), sorted(sides["patch"])

    def mapped(self, g, R, loc, to_asm):
        key = (g, R, to_asm)
        if key not in self.map_cache:
            if self.mapper is None:
                from idtrack_reproducibility.harness import assembly_map
                self.mapper = assembly_map.Mapper(self.org)
            self.map_cache[key] = self.mapper.map(loc, to_asm)
        return self.map_cache[key]

    def stratum(self, R, t):
        if R == t:
            return "same"
        d = "forward" if t > R else "backward"
        return d + ("_handoff" if self.asm.get(R) != self.asm.get(t) else "")

    def truth_at(self, g, read, t, E, B):
        """(expect, broad, apart, locus, locus_mapped, open, rest) of gene g read at release `read`, at target t."""
        tr = self.tr
        e, b = E.get(t, set()), B.get(t, set()) | E.get(t, set())
        Ev, Bv = tr.vset(e, t), tr.vset(b, t)
        Ev2, apart = self.apart_of(Ev, t)
        out = dict(cmp=tr.cmp_for(t), expect=sorted(Ev2), broad=Bv if set(Bv) != set(Ev) else None)
        if apart:
            out["apart"] = apart
        if not e:
            loc = tr.locus(g, read)
            if loc:
                out["locus"] = loc
                to_asm = self.asm.get(t)
                proj = loc
                if self.exc.has(read) or self.exc.has(t):  # the corrected judge `exceptions`: places on the reference sequence
                    out["exceptions"] = True
                    proj = self.exc.project(loc, read)
                    if proj is not loc:
                        out["locus_projected"] = proj
                if to_asm and to_asm != loc["assembly"]:
                    out["locus_mapped"] = dict(unmapped=proj["unplaced"]) if proj.get("unplaced") else \
                        self.mapped(g, read, proj, to_asm)
            else:
                out["open"] = "Ensembl has no successor and the query's place is unknown (gene absent from the GTF)"
        if t > read and self.rest.arch and self.tr.mito_cut and g in self.tr.mito_lineage():
            # REST reads Ensembl's history as recorded, the steps rule `no_mito_nuclear_step` refuses among them: no
            # second reading of the walk the rule makes
            self.counts["rest not compared (rule no_mito_nuclear_step)"] += 1
        elif t > read and self.rest.arch:
            r, complete = self.rest.at(g, t)
            if not complete:
                self.counts["rest incomplete"] += 1
            elif r != e:
                out["open"] = "the walk of Ensembl's history and REST archive/id differ"
                out["rest"] = tr.vset(r, t)
                self.counts["rest differs"] += 1
            else:
                self.counts["rest agrees"] += 1
        return out

    def lines_of(self, g, unit):
        tr = self.tr
        rels = tr.rels
        rr = sorted(tr.present[g])
        walks = {}

        def walk(R):
            if R not in walks:
                walks[R] = (tr.walk_series(g, R, rels, True, living=True),
                            tr.walk_series(g, R, rels, False, living=True, resurrect=True))
            return walks[R]

        with_db = zlib.crc32(g.encode()) % 6 == 0
        pp = None
        for R in rr:
            E, B = walk(R)
            v = tr.ver(g, R)
            for t in rels:
                base = self.truth_at(g, R, t, E, B)
                st = self.stratum(R, t)
                for kind, q in (("astar_from", g), ("astar_from_v", f"{g}.{v}" if v else None)):
                    if q is None:
                        continue
                    yield dict(qid=f"{kind}|{q}|{R}|{t}", family=FAMILY, kind=kind, stratum=st, unit=unit, id=q,
                               from_release=R, to_release=t, strategy="both", rule="set", **base)
                if with_db and t in self.bm_rels and E.get(t):
                    if pp is None:
                        pp = {}
                    if t not in pp:
                        pp[t] = self.ctx.bm(t, self.grp, self.col)
                    if pp[t]:
                        e, b = E.get(t, set()), B.get(t, set()) | E.get(t, set())
                        vals = set().union(*(pp[t][0].get(x, set()) for x in e))
                        bvals = set().union(*(pp[t][0].get(x, set()) for x in b))
                        line = dict(qid=f"astar_db|{g}|{R}|{t}|{self.db}", family=FAMILY, kind="astar_db", stratum="db",
                                    unit=unit, id=g, from_release=R, to_release=t, final_database=self.db,
                                    strategy="both", rule="set", cmp="exact", expect=sorted(vals),
                                    broad=sorted(bvals) if bvals != vals else None)
                        sp = spellings_of(line["expect"])
                        if sp:
                            line["spellings"] = sp  # the corrected judge `spellings_once`
                        if base.get("open"):
                            line["open"] = base["open"]
                        yield line
        # without from_release (rule `start_release`)
        vers = sorted({tr.ver(g, r) for r in rr if tr.ver(g, r)}, key=int)
        for q, withv in [(g, rr)] + [(f"{g}.{v}", [r for r in rr if tr.ver(g, r) == v]) for v in vers]:
            kind = "astar_default_v" if q != g else "astar_default"
            for t in rels:
                read = t if t in withv else max(withv)
                E, B = walk(read)
                base = self.truth_at(g, read, t, E, B)
                yield dict(qid=f"{kind}|{q}|-|{t}", family=FAMILY, kind=kind, stratum="default", unit=unit, id=q,
                           from_release=None, to_release=t, strategy="both", rule="set", read_at=read, **base)


def spellings_of(values):
    """The corrected judge `spellings_once` (make_queries.rev_spellings_once): spellings of one value that differ only in
    letter case, grouped (`IFN-DELTA-8`, `IFN-delta-8`); an answer with either counts as that value."""
    groups = collections.defaultdict(list)
    for v in values:
        groups[v.casefold()].append(v)
    return sorted(sorted(g) for g in groups.values() if len(g) > 1)


def selected(org, ids_file=None, name=None):
    rows = H.read_tsv(os.path.join(DATA, "selection", f"{name or org}.tsv.gz"))
    out = [(r["id"], r["lineage"], r["patterns"]) for r in rows if r["releases"] != "0"]
    if ids_file:
        keep = {x.strip() for x in open(ids_file) if x.strip()}
        out = [x for x in out if x[0] in keep]
    return out


def publish_stream(org, name, gen_lines, description, sources, dry, generator=None, generator_sha256=None):
    """F0's `lines.publish`, streaming (the A* sets hold tens of millions of lines): the same file layout, schema and
    manifest fields; refuses an existing name."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    d = L.query_dir(FAMILY, org)
    os.makedirs(d, exist_ok=True)
    if L._set_file(d, name) or os.path.exists(os.path.join(d, f"{name}.manifest.json")):
        raise SystemExit(f"{FAMILY}/{org}/{name} is published; a published set is never rewritten")
    path = os.path.join(d, f"{name}.parquet")
    tmp = path + ".tmp"
    writer = None if dry else pq.ParquetWriter(tmp, L.schema(), compression="zstd")
    counts, n, buf = collections.Counter(), 0, []
    for q in gen_lines:
        counts[q["stratum"]] += 1
        n += 1
        if not dry:
            buf.append(L.to_row(q))
            if len(buf) >= CHUNK:
                writer.write_table(pa.Table.from_pylist(buf, schema=L.schema()))
                buf = []
        if n % 2_000_000 == 0:
            log("lines", n)
    if dry:
        return dict(n=n, counts=dict(counts))
    if buf:
        writer.write_table(pa.Table.from_pylist(buf, schema=L.schema()))
    writer.close()
    os.replace(tmp, path)
    seqs = [json.load(open(os.path.join(d, f)))["seq"] for f in os.listdir(d) if f.endswith(".manifest.json")]
    gen = generator or os.path.abspath(__file__)
    man = dict(name=name, family=FAMILY, organism=org, seq=max(seqs, default=0) + 1, format="parquet", n=n,
               counts=dict(sorted(counts.items())), sha256=L.sha256(path),
               published=datetime.datetime.now().isoformat(timespec="seconds"), description=description,
               generator=gen, generator_sha256=generator_sha256 or GEN_SHA256, revises=[], judge=None, sources=sources)
    json.dump(man, open(os.path.join(d, f"{name}.manifest.json"), "w"), indent=1)
    os.chmod(path, 0o400)
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--set", required=True)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--ids", default=None)
    ap.add_argument("--selection", default=None, help="selection file name (default: the organism)")
    ap.add_argument("--limit-ids", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org = a.organism
    log("truth", org)
    G = Gen(org, a.target)
    sel = selected(org, a.ids, a.selection)
    if a.limit_ids:
        sel = sel[: a.limit_ids]
    log("ids", len(sel), "rest records", len(G.rest.arch), G.rest.meta.get("releases"))

    def all_lines():
        for i, (g, unit, _) in enumerate(sel):
            yield from G.lines_of(g, unit)
            if i % 1000 == 0:
                log("ids done", i, dict(G.counts))

    sname = a.selection or org
    sources = [f"selection/{sname}.tsv.gz (sha256 {L.sha256(os.path.join(DATA, 'selection', f'{sname}.tsv.gz'))})"
               + (f"; ids of {a.ids}" if a.ids else ""),
               "truth: gtf/genes_<release>.tsv.gz of the window", f"core/{G.tr.core}.stable_id_event + mapping_session",
               "rest/archive_astar.json", "chain/*.chain.gz", f"biomart/<release>/{G.grp}.tsv"] + \
              [f"core/{f}" for f in TA.CORE_GENES.get(org, {}).values()]
    man = publish_stream(org, a.set, all_lines(),
                         f"A* hardest lineages: {len(sel)} ids; every release of the id to every release of the window "
                         f"(stable id and versioned), without from_release (id and each version), to {G.db} for every "
                         f"sixth id; strategy both; truth.py (effective session placement, core genes) and REST as a "
                         f"second reading", sources, a.dry_run)
    man["rest"] = dict(G.counts)
    log(json.dumps(man))
    if not a.dry_run:
        mp = os.path.join(L.query_dir(FAMILY, org), f"{a.set}.manifest.json")
        m = json.load(open(mp))
        m["rest_reading"] = dict(G.counts)
        m["rest_service_releases"] = G.rest.meta.get("releases")
        json.dump(m, open(mp, "w"), indent=1)


if __name__ == "__main__":
    main()
