"""Family E's relations of one question asked two ways: a batch against single calls (`e8_batch`) and one release
deciding a string's reading (`e8_onerelease`), from Ensembl's tables only.

Usage: batch_e.py --organism <o> [--sets e8_batch,e8_onerelease] [--dry]

No truth decides these lines: each is a relation (rule `relation`, `op: equal`, under `best` and `all`) between two
lines of the set, the reference line being `open` ("the reference of the relation").

**`e8_batch`** (rule [`batch_agrees_with_single`]): `convert_many` must answer each item as `convert` answers it
alone with the same arguments, but for the copies rule `made_unique_copies` defines.

- *whole-release lists*: for every release R of the window whose core database is read (`core_display/<R>.tsv.gz`),
  the list of every display name of R's genes (`truth_b.clean`, sorted), converted as one batch at R
  (`from_release=R`, `to_release=R`: stratum `e8_batch_release_atR`) and without a release (`e8_batch_release_default`);
  every item is Ensembl's display name at R, so no item is a copy and every item must agree. The lines ask a seeded
  sample of `SAMPLE` positions per list and mode, and every position whose item is another item of the list followed
  by `-<digits>`, `.<digits>` or `_<digits>` (the form the copy rule reads; strata `..._suffix`). The lists are
  published beside the set (`e8_batch.lists.json.gz`, key `release|<R>`; lines name them by `batch_ref`).
- *mixed lists* (`e8_batch_mixed`, without a release): `MIXED` lists of `MIXED_N` items per organism, seeded, as a
  user's list mixes them: display names of the last release and of older releases, synonyms of the last release,
  gene ids with and without version, repeated items, and made-unique copies `X-1`, `X.1`, `X_1` of items X of the list
  that no table of the truth store holds at any release (no display name, label of any cross-reference, synonym or
  GTF name). A copy's batch answer must be its item X's single answer (`e8_batch_mixed_copy`); every other item's,
  its own.
- *the attack hunt's lists* (`e8_batch_hunt`): the pairs of finding AT-2 (`mml-mir-129`, `mml-mir-129-2` at 58;
  `Krtap13`, `Krtap13-1` at 48; `OCLN`, `OCLN_2` at 61; `SLA`, `SLA-10` without a release), both items, where the
  organism's display records hold them.

The reference lines (`e8_batch_single`): `convert(item)` with the same arguments, one per distinct question.

**`e8_onerelease`** (rule [`one_release_reading`]): a string read with `to_release=r` alone must be answered as with
`from_release=r, to_release=r` where it exists at r. The strings: a name at r (a display name or a label of a name
database of r's core cross-references) that is no name at the last release and is a synonym there
(`e8_onerelease_name_then_synonym`), or a synonym at r that is no name at r and is a name at the last release
(`e8_onerelease_synonym_then_name`); one release r per string (seeded among the releases where it holds), up to
`ONE_RELEASE_N` strings per kind (seeded). Reference: `e8_onerelease_reference`.
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import random
import re
import sys
import time

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.b import truth_b as TB
from idtrack_reproducibility.harness import harness_lib as H

FAMILY = "e"
SEED = 20261121
T = 115
SAMPLE = 1000          # positions asked per whole-release list and mode (beside every suffixed position)
MIXED, MIXED_N = 10, 400
ONE_RELEASE_N = 3000
SUFFIX = re.compile(r"^(.+)[-._](\d+)$")
ENS_ID = re.compile(rf"^{H.STABLE_ID_FORM}$")
REF_OPEN = "the reference of the relation (batch against single)"
REF_OPEN_1R = "the reference of the relation (one release decides the reading)"
HUNT_PAIRS = [("macaca_mulatta", ["mml-mir-129", "mml-mir-129-2"], 58),
              ("mus_musculus", ["Krtap13", "Krtap13-1"], 48),
              ("homo_sapiens", ["OCLN", "OCLN_2"], 61),
              ("sus_scrofa", ["SLA", "SLA-10"], None)]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def rank(*parts):
    """A seeded order key (the smallest come first): independent of the order of reading."""
    return hashlib.sha256("|".join(str(p) for p in (SEED,) + parts).encode()).hexdigest()


def copy_form(item, given):
    """The item X of `given` that `item` extends by a numbered suffix (the copy rule's form), or None."""
    m = SUFFIX.match(item)
    if m and m.group(1) in given and not ENS_ID.match(m.group(1)):
        return m.group(1)
    return None


def release_list(org, R):
    """The display names of R's genes (`core_display`, genes of R's GTF), sorted, or None."""
    disp = TB.core_display(org, R)
    G = H.gtf_genes(org, R)
    if disp is None or not G:
        return None
    return sorted({lab for g, (_, lab) in disp.items() if g in G and TB.clean(lab) and lab != g})


def rel_tag(r):
    return "-" if r is None else str(r)


class Batch:
    def __init__(self, org):
        self.org = org
        self.rels = [r for r in TB.core_display_releases(org) if r <= T and r in set(H.gtf_releases(org))]
        self.lists, self.singles, self.counts = {}, {}, collections.Counter()
        self.extra = {}

    def single(self, item, fr, to):
        qid = f"e8b|single|{rel_tag(fr)}|{rel_tag(to)}|{item}"
        if qid not in self.singles:
            self.singles[qid] = dict(qid=qid, family=FAMILY, kind="e8_batch", stratum="e8_batch_single", unit=item,
                                     id=item, from_release=fr, to_release=to, final_database=None, strategy="both",
                                     rule="set", cmp="exact", expect=[], open=REF_OPEN)
        return qid

    def batch_line(self, qid, stratum, key, pos, item, fr, to, with_qid, unit=None):
        self.counts[stratum] += 1
        return dict(qid=qid, family=FAMILY, kind="e8_batch", stratum=stratum, unit=unit or item, id=item,
                    call="convert_many", batch_ref={"set": "e8_batch", "key": key}, pos=pos, from_release=fr,
                    to_release=to, final_database=None, strategy="both", rule="relation", cmp="exact", expect=[],
                    relation={"op": "equal", "with": with_qid})

    def whole_release(self):
        out = []
        for R in self.rels:
            items = release_list(self.org, R)
            if not items:
                continue
            key = f"release|{R}"
            self.lists[key] = items
            given = set(items)
            for mode, fr, to in (("atR", R, R), ("default", None, None)):
                take = set(sorted(range(len(items)), key=lambda i: rank(self.org, key, mode, items[i]))[:SAMPLE])
                for i, v in enumerate(items):
                    suf = copy_form(v, given) is not None
                    if i not in take and not suf:
                        continue
                    st = f"e8_batch_release_{mode}" + ("_suffix" if suf else "")
                    out.append(self.batch_line(f"e8b|batch|{mode}|{R}|{v}", st, key, i, v, fr, to, self.single(v, fr, to)))
            log(self.org, "release list", R, len(items), "items")
        return out

    def known_strings(self, cands):
        """The candidates that any table of the truth store holds at any release (names, labels, synonyms, GTF names)."""
        seen = set()
        d = H.truth_dir(self.org)
        for sub, col in (("core_display", 4), ("core_xref", 6), ("core_synonyms", 0)):
            p = os.path.join(d, sub)
            for f in sorted(os.listdir(p)) if os.path.isdir(p) else []:
                if not f.endswith(".tsv.gz"):
                    continue
                with gzip.open(os.path.join(p, f), "rt") as fh:
                    next(fh)
                    for line in fh:
                        c = line.rstrip("\n").split("\t")
                        if len(c) > col and c[col] in cands:
                            seen.add(c[col])
        for R in H.gtf_releases(self.org):
            for r in (H.gtf_genes(self.org, R) or {}).values():
                if r.get("gene_name") in cands:
                    seen.add(r["gene_name"])
        return seen

    def mixed(self):
        rng = random.Random(f"{SEED}|{self.org}|mixed")
        last = self.rels[-1]
        cur = release_list(self.org, last) or []
        curset = set(cur)
        old = sorted({v for R in self.rels[:-1] for v in (release_list(self.org, R) or []) if v not in curset})
        syn = sorted(v for v in (TB.core_synonyms(self.org, last) or {}) if TB.clean(v) and v not in curset)
        G = H.gtf_genes(self.org, last) or {}
        gids = sorted(G)
        plan = []
        for k in range(MIXED):
            items = rng.sample(cur, min(len(cur), int(MIXED_N * 0.4))) + rng.sample(old, min(len(old), int(MIXED_N * 0.2))) \
                + rng.sample(syn, min(len(syn), int(MIXED_N * 0.1)))
            for g in rng.sample(gids, min(len(gids), int(MIXED_N * 0.15))):
                items.append(f"{g}.{G[g]['gene_version']}" if rng.random() < 0.5 and G[g]["gene_version"] else g)
            names = [x for x in items if not ENS_ID.match(x.split(".")[0])]
            copies = [(f"{x}{rng.choice('-._')}{rng.choice((1, 1, 2, 3))}", x) for x in rng.sample(names, min(len(names), int(MIXED_N * 0.1)))]
            dups = rng.sample(items, min(len(items), int(MIXED_N * 0.05)))
            plan.append((items, copies, dups))
        known = self.known_strings({c for _, cs, _ in plan for c, _ in cs})
        out = []
        for k, (items, copies, dups) in enumerate(plan):
            lst = items + [c for c, _ in copies if c not in known] + dups
            rng.shuffle(lst)
            key = f"mixed|{k}"
            self.lists[key] = lst
            base = {c: x for c, x in copies if c not in known}
            for i, v in enumerate(lst):
                if v in base:
                    out.append(self.batch_line(f"e8b|batch|mixed|{k}|{i}", "e8_batch_mixed_copy", key, i, v, None, None,
                                               self.single(base[v], None, None)))
                else:
                    out.append(self.batch_line(f"e8b|batch|mixed|{k}|{i}", "e8_batch_mixed", key, i, v, None, None,
                                               self.single(v, None, None)))
        self.extra["mixed_copies_dropped_as_known"] = sum(1 for _, cs, _ in plan for c, _ in cs if c in known)
        return out

    def hunt(self):
        out = []
        for k, (org, items, R) in enumerate(HUNT_PAIRS):
            if org != self.org:
                continue
            key = f"hunt|{k}"
            self.lists[key] = items
            for i, v in enumerate(items):
                out.append(self.batch_line(f"e8b|batch|hunt|{k}|{i}", "e8_batch_hunt", key, i, v, R, R,
                                           self.single(v, R, R)))
        return out

    def lines(self):
        out = self.whole_release() + self.mixed() + self.hunt()
        self.counts["e8_batch_single"] = len(self.singles)
        # the lines of one list together (a batch is converted once per list and arguments: adapter_e)
        out.sort(key=lambda q: (q["batch_ref"]["key"], rel_tag(q["from_release"]), q["pos"]))
        return out + sorted(self.singles.values(), key=lambda q: q["qid"])


def names_at(org, R):
    """The names of R: display names and labels of the name databases of R's core cross-references (clean strings)."""
    out = set()
    disp = TB.core_display(org, R)
    if disp:
        out |= {lab for _, lab in disp.values() if TB.clean(lab)}
    cx = TB.core_names(org, R)
    if cx:
        out |= {v for v in cx[0] if TB.clean(v)}
    return out


def one_release(org):
    rels = [r for r in TB.core_display_releases(org) if r <= T]
    last = rels[-1]
    n_last = names_at(org, last)
    syn_last = {v for v in (TB.core_synonyms(org, last) or {}) if TB.clean(v)}
    a_cand = syn_last - n_last
    where = {"name_then_synonym": collections.defaultdict(list), "synonym_then_name": collections.defaultdict(list)}
    syn_rels = []
    for r in rels[:-1]:
        nr = names_at(org, r)
        for s in nr & a_cand:
            where["name_then_synonym"][s].append(r)
        sr = TB.core_synonyms(org, r)
        if sr is None:
            continue
        syn_rels.append(r)
        for s in (set(sr) & n_last) - nr:
            if TB.clean(s):
                where["synonym_then_name"][s].append(r)
    out, counts = [], collections.Counter()
    for kind, d in where.items():
        pick = sorted(d, key=lambda s: rank(org, kind, s))[:ONE_RELEASE_N]
        for s in sorted(pick):
            r = random.Random(rank(org, kind, s, "release")).choice(sorted(d[s]))
            ref = f"e8r|from|{s}|{r}"
            out.append(dict(qid=ref, family=FAMILY, kind="e8_onerelease", stratum="e8_onerelease_reference", unit=s, id=s,
                            from_release=r, to_release=r, final_database=None, strategy="both", rule="set", cmp="exact",
                            expect=[], open=REF_OPEN_1R, holds=kind))
            out.append(dict(qid=f"e8r|to|{s}|{r}", family=FAMILY, kind="e8_onerelease", stratum=f"e8_onerelease_{kind}",
                            unit=s, id=s, from_release=None, to_release=r, final_database=None, strategy="both",
                            rule="relation", cmp="exact", expect=[], relation={"op": "equal", "with": ref}, holds=kind))
            counts[kind] += 1
    extra = dict(candidates={k: len(d) for k, d in where.items()}, asked=dict(counts), last=last,
                 synonym_releases=syn_rels)
    return out, extra


def sources(org):
    d = H.truth_dir(org)
    out = []
    for sub in ("core_display", "core_xref", "core_synonyms"):
        p = os.path.join(d, sub)
        out += [f"{org}/{sub}/{f}" for f in sorted(os.listdir(p)) if f.endswith(".tsv.gz")] if os.path.isdir(p) else []
    return out + [f"{org}/gtf/genes_<every release>.tsv.gz"]


def publish_lists(org, lists, dry):
    """`e8_batch.lists.json.gz` beside the set (written once, read-only); its sha256."""
    p = os.path.join(L.query_dir(FAMILY, org), "e8_batch.lists.json.gz")
    if dry:
        return None
    if os.path.exists(p):
        raise SystemExit(f"{p} exists: a published set's lists are never rewritten")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with gzip.open(p + ".tmp", "wt") as fo:
        json.dump(lists, fo, sort_keys=True)
    os.replace(p + ".tmp", p)
    os.chmod(p, 0o400)
    return L.sha256(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="e8_batch,e8_onerelease")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    for s in a.sets.split(","):
        if s == "e8_batch":
            b = Batch(a.organism)
            lines = b.lines()
            log("e8_batch", len(lines), dict(sorted(b.counts.items())))
            if a.dry:
                print(dict(set=s, n=len(lines), counts=dict(sorted(b.counts.items())), lists=len(b.lists),
                           list_items=sum(len(x) for x in b.lists.values()), **b.extra))
                continue
            sha = publish_lists(a.organism, b.lists, a.dry)
            man = L.publish(FAMILY, a.organism, s, lines, "E8: a batch against single calls (whole-release lists of "
                            "display names, mixed lists, the attack hunt's lists)", __file__, judge="now",
                            sources=sources(a.organism), fmt="parquet",
                            manifest_extra=dict(b.extra, lists="e8_batch.lists.json.gz", lists_sha256=sha,
                                                lists_n=len(b.lists), seed=SEED, sample=SAMPLE))
            log("published", s, man["n"])
        elif s == "e8_onerelease":
            lines, extra = one_release(a.organism)
            if a.dry:
                print(dict(set=s, n=len(lines), **extra))
                continue
            man = L.publish(FAMILY, a.organism, s, lines, "E8: one release decides a string's reading (to_release alone "
                            "against from_release and to_release)", __file__, judge="now", sources=sources(a.organism),
                            manifest_extra=dict(extra, seed=SEED, cap=ONE_RELEASE_N))
            log("published", s, man["n"])
        else:
            raise SystemExit(f"unknown set {s}")


if __name__ == "__main__":
    main()
