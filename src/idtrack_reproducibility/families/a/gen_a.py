"""Family A, history travel exhaustive (docs/methods.md, A): the query sets of one organism and one stratum,
published once, never rewritten.

Usage: gen_a.py <organism> --stratum A1|A2|A3|A4|A5|A6|A8 --set <name> [--snapshot S (A8)] [--target 115] [--limit-ids N] [--dry-run]

Strata (every line `strategy: both`, scored under `best` and `all`; qid `<kind>|<query>|<from or ->|<to>`):
- A1 `a1`, `a1_v`: every gene id at every release R where the GTF has it, from R to the target T; the stable id and,
  where the GTF of R writes versions, the id with its version at R.
- A2 `a2`, `a2_v`: every gene id and every version it had, without `from_release`, to T (rule `start_release`: read at T if it
  exists there, with that version for a versioned id; else at the latest release where it exists in that form).
- A3 `a3`: every (id, R) to R-2, R-1, R+1, R+2 (releases of the window) and to the first release of the window, the
  stable id; `a3_v` the id with its version at R to R-1 and R+1.
- A4 `a4_fwd(_v)`, `a4_back(_v)`, `a4_ver`: every cross-id gene event of the release-115 `stable_id_event` in the
  window: the old id from the session's old release (the last GTF release before the session takes effect) to its
  new release, and the new id back; with the event's versions; every change of version of one id (`a4_ver`), the old
  version forward across its session.
- A5 `a5_fwd(_v)`, `a5_back(_v)`: every handoff between two assemblies (last release L of one, first F of the next):
  every id present at L to F, every id present at F to L; stratum `A5_<L>_<F>_<lives|dies|born>`.
- A6 `a6`: the full matrix (every release of the id to every release of the window), the stable id: every id of
  cynomolgus; for the other organisms a seeded sample of 1,000 ids per history class (the travel matrix's: unchanged,
  born_late, renamed, merged, split, retired_none, retired_branch) and assembly era of the id's first release.
- A8 `a8`, `a8_v`: A1 against the graph of an older snapshot S (issue #71): every (id, R <= S) to S.

Truth (nothing of idtrack): A*'s `truth.TruthA` (the harness's `make_queries.Truth` with sessions placed where their
effect appears, plus core genes the GTF lacks) read by A*'s `gen.Gen.truth_at` (strict walk with living ids =
`expect`; broad walk with returned ids = `broad`; copies on patches apart (rule `judge_copies_apart`); the locus rule with the place
projected through Ensembl's assembly exceptions and mapped across assemblies by the chain files). REST is not read
(it is A*'s second reading for its ids). A walk is computed directly for an id without a cross-id event and without
a gap in its releases (its walk is itself wherever it exists; checked against the full walk on a sample at start).

The walk follows rule `no_mito_nuclear_step` (TruthA): a step between the mitochondrial genome and a nuclear sequence
is not followed, also where the stable id stays the same (such an id has no direct walk).

`unit` (cluster of the bootstrap): the id's lineage (the smallest id of its component of cross-id events), else the
id. Extra fields: `era` (assembly of the release asked from), `dist` (target minus origin), `cross` (assembly differs).
"""
import argparse
import collections
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.astar import gen as AG  # A*'s generator: Gen.truth_at
from idtrack_reproducibility.families.astar import truth as TA
from idtrack_reproducibility.harness import harness_lib as H

FAMILY = "a"
ME = os.path.abspath(__file__)
ME_SHA256 = L.sha256(ME)
A6_PER = 1000
SEED = 20261118


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _int(x):
    try:
        int(x)
        return True
    except (TypeError, ValueError):
        return False


class ExceptionsA:
    """The harness's `assembly_map.Exceptions` without the rows that carry no coordinates (human 110 to 116 write
    `na` for GRC placements without them; `Exceptions.project` raises on them). Such a row cannot place anything; a
    sequence that had one is remembered (`na`), and a line whose place is on it is left open (GenA.truth_at)."""

    def __init__(self, exc):
        self.exc, self.na, self.filtered = exc, {}, {}

    def rows(self, rel):
        if rel not in self.filtered:
            keep, na = {}, set()
            for seq, rs in self.exc.rows(rel).items():
                good = [r for r in rs if all(_int(r[k]) for k in ("start", "end", "ref_start", "ref_end"))]
                if len(good) != len(rs):
                    na.add(seq)
                if good:
                    keep[seq] = good
            self.filtered[rel], self.na[rel] = keep, na
        return self.filtered[rel]

    def has(self, rel):
        return bool(self.rows(rel)) or bool(self.na.get(rel))

    def project(self, place, rel):
        self.rows(rel)
        saved = self.exc.cache.get(rel)
        self.exc.cache[rel] = self.filtered[rel]
        try:
            return self.exc.project(place, rel)
        finally:
            self.exc.cache[rel] = saved


class GenA(AG.Gen):
    def __init__(self, org, T):
        super().__init__(org, T)
        self.rest.arch = {}  # REST is A*'s second reading, not A's
        self.exc = ExceptionsA(self.exc)
        # A*'s TruthA gives the core genes the GTF lacks (human 48-51) the core table's version, while the GTF of
        # those releases writes none and their comparison is by stable id (`cmp: base`): an expected `<id>.<v>` then
        # never equals the answer. At a release whose GTF writes no versions, those genes carry none either.
        tr = self.tr
        self.core_versions_dropped = 0
        for r, extra in getattr(tr, "core_only", {}).items():
            gtf_versions = any((tr.genes[r][x].get("gene_version") or "") for x in tr.genes[r] if x not in extra)
            if not gtf_versions:
                for x in extra:
                    if tr.genes[r][x].get("gene_version"):
                        tr.genes[r][x]["gene_version"] = ""
                        self.core_versions_dropped += 1
        tr = self.tr
        self.cross = set()
        parent = {}

        def find(x):
            while parent.get(x, x) != x:
                parent[x] = parent.get(parent[x], parent[x])
                x = parent[x]
            return x

        for s in tr.succ:
            for o, ns in tr.succ[s].items():
                for n in ns:
                    if n and o and n != o:
                        self.cross.update((o, n))
                        a, b = find(o), find(n)
                        if a != b:
                            parent[max(a, b)] = min(a, b)
        self.unit = {x: find(x) for x in self.cross}
        self.idx = {r: i for i, r in enumerate(tr.rels)}
        self.contig = {}
        self.walk_cache = {}

    def truth_at(self, g, read, t, E, B):
        out = super().truth_at(g, read, t, E, B)
        loc = out.get("locus")
        if loc and not out.get("open") and any(loc["seqname"] in self.exc.na.get(r, ()) for r in (read, t)):
            out["open"] = "the place is on a sequence whose assembly exception has no coordinates ('na'): not projected"
            self.counts["open: exception without coordinates"] += 1
        return out

    def unit_of(self, g):
        return self.unit.get(g, g)

    def contiguous(self, g):
        if g not in self.contig:
            rr = sorted(self.tr.present[g])
            self.contig[g] = self.idx[rr[-1]] - self.idx[rr[0]] + 1 == len(rr)
        return self.contig[g]

    def walks(self, g, R, targets, fast=True):
        """(strict, broad) {target: ids} of g read at R."""
        tr = self.tr
        if fast and g not in self.cross and g not in tr.mito_steps() and self.contiguous(g):
            pr = tr.present[g]
            E = {t: ({g} if t in pr else set()) for t in targets}
            return E, E
        return (tr.walk_series(g, R, targets, True, living=True),
                tr.walk_series(g, R, targets, False, living=True, resurrect=True))

    def selfcheck(self, n=3000):
        """The direct walk equals the full walk for ids without cross-id events and gaps (a seeded sample)."""
        tr = self.tr
        pool = sorted(g for g in tr.present if g not in self.cross and g not in tr.mito_steps() and self.contiguous(g))
        rng = random.Random(f"{SEED}:{self.org}:selfcheck")
        bad = 0
        for g in rng.sample(pool, min(n, len(pool))):
            R = rng.choice(sorted(tr.present[g]))
            a = self.walks(g, R, tr.rels, fast=True)
            b = self.walks(g, R, tr.rels, fast=False)
            if a != b:
                bad += 1
        log("selfcheck direct walks", min(n, len(pool)), "differ", bad)
        if bad:
            raise SystemExit("direct walk differs from the full walk")

    def axes(self, R, t):
        a, b = self.asm.get(R, ""), self.asm.get(t, "")
        return dict(era=a, dist=(t - R) if R is not None else None, cross=a != b)

    def line(self, kind, stratum, g, q, R, t, base, **extra):
        rec = dict(qid=f"{kind}|{q}|{R if R is not None else '-'}|{t}", family=FAMILY, kind=kind, stratum=stratum,
                   unit=self.unit_of(g), id=q, from_release=R, to_release=t, strategy="both", rule="set", **base)
        rec.update(extra)
        return rec

    # ---------------------------------------------------------------- strata

    def a1(self, ids):
        tr, T = self.tr, self.T
        for g in ids:
            for R in sorted(tr.present[g]):
                E, B = self.walks(g, R, [T])
                base = self.truth_at(g, R, T, E, B)
                v = tr.ver(g, R)
                ax = self.axes(R, T)
                yield self.line("a1", "A1", g, g, R, T, base, **ax)
                if v:
                    yield self.line("a1_v", "A1", g, f"{g}.{v}", R, T, base, **ax)

    def a8(self, ids):
        """A8: A1 against a graph of an older snapshot S (`--snapshot`): every (id, R <= S) to S. The truth is read from
        the release-115 tables (Ensembl's history is cumulative; sessions after S are never crossed by a walk to S)."""
        tr, S = self.tr, self.snapshot
        for g in ids:
            for R in sorted(r for r in tr.present[g] if r <= S):
                E, B = self.walks(g, R, [S])
                base = self.truth_at(g, R, S, E, B)
                v = tr.ver(g, R)
                ax = self.axes(R, S)
                yield self.line("a8", "A8", g, g, R, S, base, snapshot=S, **ax)
                if v:
                    yield self.line("a8_v", "A8", g, f"{g}.{v}", R, S, base, snapshot=S, **ax)

    def a2(self, ids):
        tr, T = self.tr, self.T
        for g in ids:
            rr = sorted(tr.present[g])
            vers = sorted({tr.ver(g, r) for r in rr if tr.ver(g, r)}, key=int)
            for q, withv in [(g, rr)] + [(f"{g}.{v}", [r for r in rr if tr.ver(g, r) == v]) for v in vers]:
                read = T if T in withv else max(withv)
                E, B = self.walks(g, read, [T])
                base = self.truth_at(g, read, T, E, B)
                yield self.line("a2_v" if q != g else "a2", "A2", g, q, None, T, base, read_at=read,
                                **self.axes(read, T))

    def a3(self, ids):
        tr = self.tr
        rels, first = tr.rels, tr.rels[0]
        for g in ids:
            for R in sorted(tr.present[g]):
                i = self.idx[R]
                near = {rels[j] for j in (i - 2, i - 1, i + 1, i + 2) if 0 <= j < len(rels)}
                targets = sorted(near | ({first} if first != R else set()))
                if not targets:
                    continue
                E, B = self.walks(g, R, targets)
                v = tr.ver(g, R)
                for t in targets:
                    base = self.truth_at(g, R, t, E, B)
                    ax = self.axes(R, t)
                    st = "A3_first" if t == first and t not in near else f"A3_{'+' if t > R else '-'}{abs(self.idx[t] - i)}"
                    yield self.line("a3", st, g, g, R, t, base, **ax)
                    if v and abs(self.idx[t] - i) == 1 and t in near:
                        yield self.line("a3_v", st, g, f"{g}.{v}", R, t, base, **ax)

    def a4(self, ids=None):
        tr = self.tr
        core, ev, sessions = H.core_tables(tr.org, tr.T)
        seen = set()
        for o, ov, n, nv, s, sc in ev:
            if not (o and n):
                continue
            b, a = tr.sides(s)
            if b is None or a is None:
                continue
            if o != n:
                todo = [("a4_fwd", o, ov, b, a), ("a4_back", n, nv, a, b)]
            elif ov != nv:
                todo = [("a4_ver", o, ov, b, a)]
            else:
                continue
            for kind, g, v, R, t in todo:
                if R not in tr.present.get(g, ()):
                    # the event's id is not in the GTF of the release it is read at: rule `event_version`'s case, read there all
                    # the same (the truth is the walk; idtrack must not invent the id)
                    note = "the event's id is absent from the GTF of the release asked from"
                else:
                    note = None
                forms = [(kind, g)] if kind != "a4_ver" else []
                if v and tr.cmp_for(R) == "versioned":
                    forms.append((kind + ("_v" if kind != "a4_ver" else ""), f"{g}.{v}"))
                for k, q in forms:
                    key = (k, q, R, t)
                    if key in seen:
                        continue
                    seen.add(key)
                    E, B = self.walks(g, R, [t], fast=False)
                    base = self.truth_at(g, R, t, E, B)
                    ver_ok = (not v) or tr.ver(g, R) == v
                    if note:
                        # the id is in no gene table at R (human ASMPATCHG version-0 placeholders and the like): the
                        # builder leaves such events out (ids in no table are no nodes); Ensembl's history names a
                        # successor all the same, so the line is open (the reading of an id no release has is a rule)
                        base = dict(base, open="the event's id is not a gene of the release asked from (in no gene table there)")
                    if q != g and not ver_ok:
                        # the event names a version the id does not have at R in the GTF (mouse 102 -> 103: 3,179
                        # version jumps): rule `version_newer_than_graph`, a version the id never had is nothing; one it had at another
                        # release is open (the rule for a version read at the wrong release is not the walk's)
                        had = {tr.ver(g, r) for r in tr.present.get(g, ())}
                        if v not in had:
                            base = dict(base, expect=[], broad=None, apart=None, locus=None, locus_mapped=None,
                                        locus_projected=None, exceptions=None,
                                        note="rule `version_newer_than_graph`: a version the id never had in any GTF: nothing")
                            base = {k: x for k, x in base.items() if x is not None}
                        else:
                            base = dict(base, open="the event's version is the id's version at another release only")
                    st = "A4_" + kind.split("_")[1] + ("" if ver_ok else "_version_not_in_gtf")
                    yield self.line(k, st, g, q, R, t, base, session=s, score=sc, event_note=note, **self.axes(R, t))

    def a5(self, ids=None):
        tr = self.tr
        eras = tr.eras()
        for (a1, lo1, hi1), (a2, lo2, hi2) in zip(eras[:-1], eras[1:]):
            Lr, F = hi1, lo2
            for g in sorted(x for x, rr in tr.present.items() if Lr in rr or F in rr):
                rr = tr.present[g]
                fate = "lives" if (Lr in rr and F in rr) else ("dies" if Lr in rr else "born")
                st = f"A5_{Lr}_{F}_{fate}"
                for kind, R, t in (("a5_fwd", Lr, F), ("a5_back", F, Lr)):
                    if R not in rr:
                        continue
                    E, B = self.walks(g, R, [t])
                    base = self.truth_at(g, R, t, E, B)
                    v = tr.ver(g, R)
                    yield self.line(kind, st, g, g, R, t, base, **self.axes(R, t))
                    if v:
                        yield self.line(kind + "_v", st, g, f"{g}.{v}", R, t, base, **self.axes(R, t))

    def history_class(self, g):
        tr, T = self.tr, self.T
        rr = tr.present[g]
        if T in rr:
            return "unchanged" if len(rr) == len(tr.rels) else ("born_late" if min(rr) > tr.rels[0] else "changed_alive")
        Lr = max(rr)
        E = tr.walk_fwd(g, Lr, T, True, living=True)
        if not E:
            return "retired_none" if not tr.walk_fwd(g, Lr, T, False, living=True) else "retired_branch"
        if len(E) > 1:
            return "split"
        for s in tr.sessions_between(Lr, T):
            news = tr.succ[s].get(g, set()) - {None}
            if news and g not in news:
                return "merged" if any(len(tr.pred[s][n] - {None}) > 1 for n in news) else "renamed"
        return "renamed"

    def a6_ids(self):
        tr = self.tr
        ids = sorted(tr.present)
        if self.org == "macaca_fascicularis":
            return ids, {g: self.history_class(g) for g in ids}
        groups = collections.defaultdict(list)
        cls = {}
        for g in ids:
            cls[g] = self.history_class(g)
            groups[(cls[g], tr.assembly.get(min(tr.present[g]), ""))].append(g)
        chosen = []
        for k in sorted(groups):
            pop = sorted(groups[k])
            rng = random.Random(f"{SEED}:{self.org}:a6:{k[0]}:{k[1]}")
            chosen += pop if len(pop) <= A6_PER else sorted(rng.sample(pop, A6_PER))
        log("a6 groups", {f"{a}|{b}": len(v) for (a, b), v in sorted(groups.items())})
        return sorted(chosen), cls

    def a6(self, ids):
        tr = self.tr
        rels = tr.rels
        for g in ids:
            hc = self.a6_class[g]
            for R in sorted(tr.present[g]):
                E, B = self.walks(g, R, rels)
                for t in rels:
                    base = self.truth_at(g, R, t, E, B)
                    yield self.line("a6", f"A6_{self.stratum(R, t)}", g, g, R, t, base, history_class=hc,
                                    **self.axes(R, t))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--stratum", required=True, choices=["A1", "A2", "A3", "A4", "A5", "A6", "A8"])
    ap.add_argument("--snapshot", type=int, default=None, help="A8: the snapshot release of the older graph")
    ap.add_argument("--set", required=True)
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--limit-ids", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    org = a.organism
    log("truth", org)
    G = GenA(org, a.target)
    G.selfcheck()
    ids = sorted(G.tr.present)
    extra = {}
    if a.stratum == "A8":
        if not a.snapshot:
            raise SystemExit("A8 needs --snapshot")
        G.snapshot = a.snapshot
        ids = [g for g in ids if min(G.tr.present[g]) <= a.snapshot]
        extra["snapshot"] = a.snapshot
    if a.stratum == "A6":
        ids, G.a6_class = G.a6_ids()
        extra["a6_ids"] = len(ids)
        extra["a6_history_classes"] = dict(collections.Counter(G.a6_class[g] for g in ids))
    if a.limit_ids:
        ids = ids[: a.limit_ids]
    log("ids", len(ids), "with cross-id events", len(G.cross))
    gen = getattr(G, a.stratum.lower())

    def all_lines():
        k = 0
        for rec in gen(ids):
            k += 1
            if k % 1_000_000 == 0:
                log("lines", k, dict(G.counts))
            yield rec

    sources = ["truth: gtf/genes_<release>.tsv.gz of the window", f"core/{G.tr.core}.stable_id_event + mapping_session",
               "chain/*.chain.gz", "exceptions/<release>.tsv"] + [f"core/{f}" for f in TA.CORE_GENES.get(org, {}).values()]
    man = L.publish(FAMILY, org, a.set, all_lines(),
                    f"A {a.stratum} (docs/methods.md, A; gen_a.py docstring); strategy both; A*'s TruthA and "
                    f"Gen.truth_at (strict walk with living ids, broad, copies apart, locus with exceptions and chains)",
                    ME, sources=sources,
                    manifest_extra=dict(stratum=a.stratum, generator_sha256_at_start=ME_SHA256,
                                        astar_gen_sha256=L.sha256(AG.__file__), astar_truth_sha256=L.sha256(TA.__file__),
                                        make_queries_sha256=L.sha256(AG.M.__file__),
                                        assembly_map_sha256=L.sha256(os.path.join(C.ACCURACY, "assembly_map.py")),
                                        target=a.target, limit_ids=a.limit_ids, **extra),
                    dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k not in ("counts",)})
    log("counts", (man or {}).get("counts"))


if __name__ == "__main__":
    main()
