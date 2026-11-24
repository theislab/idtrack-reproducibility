"""Revision sets of family B: every set regenerated with the current truth and compared, qid by qid, with what is
published; only the lines whose truth changed are published, under their qids (F0's rule: a published set is never
rewritten).

Usage: revise_b.py --organism <o> --name r1 --why "<what changed in the truth>" [--sets b1,b4,b5,b6,b7,b8,b9] [--dry]

For each set family S (`b1`, ...):
- `S_<name>`: the lines of the `now` judge (every set of S whose name holds no `_dec`) whose regenerated line differs
  from the line the published sets give (`revised_by` names the revision); a published qid that the generator no
  longer makes is published again with `excluded` (left out of the scoring, counted in run.json). The regenerated
  line is first given the corrections of the judge published for S since its base sets, in their publication order
  (`replay` of `chain`: the harness's revisions of decided rules of each `--rules` set, each `--correction` set's
  correction), as they were given to the published lines; otherwise a line a correction revised would be published
  again without it. A line whose regenerated base is its published base but whose replayed line is not its published
  line is a correction the replay does not reproduce: then nothing is published (`replay mismatches` in the log).
  With `--dry --lines-to <file>` every changed line is written beside its published line, with the fields that
  differ (a change of a field that asks, `QUESTION`, is a new question: a stored answer does not answer it);
- `S_<name>_dec` (B6, B7): the lines the `decided` judge must read after `S_<name>` (sets in order: `S_v1`, `S_dec`,
  `S_<name>`, `S_<name>_dec`): every qid whose decided line (its regenerated `_dec` line, else its regenerated line)
  differs from what that order gives;
- `S_v2`: lines with qids not published before (a new set, never a revision): B5's synonyms that Ensembl's records
  of every cross-reference add at a BioMart release to BioMart's synonym table (`name_src: core_synonyms`).

`--rules <revision>,...` (with `--name`): instead of regenerating, the latest lines of each set family (judge `now`)
are revised by the accuracy harness's revisions of decided rules (`make_queries.REVISIONS`, judge `rules`:
`par_shared_names`, `capitals_mgi` and `lower_case_symbols` (rule `single_case_symbols`), `missing_values`,
`stray_punctuation`), read with B's truth context (`truth_b.NamesB.C`), and the lines they change are published as
`S_<name>`, judge `rules`. A family with a decided judge (`S_dec`: B6, B7) also gets `S_<name>_dec`: the latest lines
of its `_dec` sets with the same revisions, and every one of them whose qid `S_<name>` holds, so that the decided
truth and the rule both stand whichever order the sets are read in (the decided judge's, or publication order).
The `now` judge reads the sets without `_dec` in their name, the `decided` judge every set (the family's methods, docs/families/b.md).

`--correction places_not_comparable` (with `--name`; the corrected judge, judge `78`): a B1 line of a name read at a
release R and asked at a target t (`b1_name_to_current`) whose genes Ensembl's history carries to nothing at t expects
nothing; where one gene carries the name it holds that gene's place, and an answer there is undetermined by the locus
rule. A name that several genes carry (`b1_*_shared`) holds no place, so any answer was wrong, also where R and t lie
on two assemblies that no chain file of Ensembl joins (Sscrofa9 or Sscrofa10.2 and Sscrofa11.1, the macaques'
assemblies), where the locus rule calls the answer of a single carrier undetermined (places not comparable). Such a
line (no expected gene, no place, R and t on two assemblies without a chain) is `open` with that reason, as the line of
a single carrier is undetermined: published as `S_<name>`. A line whose name a decided rule says is not held at the
release read (`not_held`) expects nothing because nothing is found, and is left alone; `--correction places_undo
--undo <name>` publishes again, as `S_<name>`, the line before `S_<undo>` of every such line that set opened.

`--correction support_from_xrefs` (with `--name`; the corrected judge, judge `78`): rule `support_ranks_shared_name`
(a `best` answer that is part of a name's expected genes, whose genes carry the name in strictly more name sources
than every gene left out, is undetermined) needs each expected gene's support; the harness reads it from BioMart's
tables, which most early releases lack. A B1 line of a name several genes carry (`set_exact`) gets `support`: per
expected gene, the number of Ensembl's name sources that give the gene the name at the release read (each database
of the core cross-references whose label it is, the display name; for a line asked at a later target, the largest
support of the genes whose strict walk reaches it), where the genes differ; published as `S_<name>` (`_support`).

`--correction support_names` (with `--name`; judge `78`): the support of `support_from_xrefs` counted Ensembl's
synonyms as a name source; idtrack reads a string as the names that carry it (rule `names_before_synonyms`) and
counts the databases attaching those (`Candidate.entry_support`), so a synonym is no source (BioMart's synonyms in
the generator's support, Ensembl's core synonyms in `support_from_xrefs`'s): every B1 line of a name several genes
carry, with or without `support`, gets the support the published sets give it in their order, of the name sources
alone (`_support`), and the lines where that differs from what they hold are published as `S_<name>` (`support`
set, changed, or dropped where the genes no longer differ).

`--correction held_in_core` (with `--name`; the corrected judge, judge `78`): B4's and B5's lines of a string that the
truth found in no Ensembl table of any release (BioMart's tables and the core extracts it read when it was generated)
expect nothing; the core extracts of every release (`core_xref` labels, `core_display`, `core_synonyms`, fetched
since) hold some of them. Such a line is read as rule `start_release` reads a string (at the target where it is held
there, else at the latest release where it is; names before synonyms), its genes walked to the target; published as
`S_<name>`. A line rule `mirbase_genes` decided is left alone (`decided_by_mirbase`); `--correction places_undo --undo
held<k>` puts back, as `S_<name>`, the lines of `S_held<k>` it had revised before this.
"""
import argparse
import collections
import copy
import functools
import gzip
import hashlib
import json
import os
import sys
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.b import gen_b as GB
from idtrack_reproducibility.families.b import gen_b2 as G2
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.b import truth_b as TB

FAMILY = "b"
log = GB.log


# The fields that explain a line in words or say what the truth read (`note`, `name_src`): no part of its truth, so a
# change of their wording (a decision's number become a rule's name) or of what was read with the same result is no
# revision. `open` and `excluded` count by whether they are set, not by their wording.
WORDING = ("note", "name_src")
PRESENCE = ("open", "excluded")


def canon(line):
    """A line as compared between its published and its regenerated form: the truth and the question, without the
    bookkeeping of the revisions that published it (`revised_by`, `revision_note`) and without its wording (`WORDING`,
    `PRESENCE`), so that a later revision does not publish again the lines an earlier one revised, nor lines whose truth
    is the same in other words."""
    r = L.to_row({k: v for k, v in line.items() if k not in ("set", "revised_by", "revision_note")})
    r.pop("h", None)
    r.pop("revised_by", None)
    if r.get("extra"):
        x = json.loads(r["extra"])
        x.pop("set", None)
        for k in WORDING:
            x.pop(k, None)
        for k in PRESENCE:
            if k in x:
                x[k] = bool(x[k])
        r["extra"] = json.dumps(x, sort_keys=True) if x else None
    for k in WORDING:
        r.pop(k, None)
    for k in PRESENCE:
        if k in r:
            r[k] = bool(r[k])
    return json.dumps(r, sort_keys=True, default=str)


def latest(org, sets):
    """{qid: line} over the sets in their order (a later set's line replaces an earlier one)."""
    out = {}
    d = L.query_dir(FAMILY, org)
    for s in sets:
        for ln in L.iter_set(L._set_file(d, s)):
            ln["set"] = s
            out[ln["qid"]] = ln
    return out


def family_sets(org, S):
    names = L.set_names(FAMILY, org)
    return [s for s in names if s == S or s.startswith(S + "_")]


def regenerate(cx, S):
    """(base lines, dec lines) of set family S with the current truth."""
    if S == "b1":
        return GB.GenB1(cx.nb).lines(), []
    if S == "b8":
        return GB.b8_lines(cx.nb), []
    if S == "b9":
        from idtrack_reproducibility.families.b import gen_b9 as B9
        return list(B9.GenB9(cx.org).lines()), []
    got = {}

    def capture(cx_, name, lines_, desc, dry, revises=None):
        got[name] = list(lines_)
        return {"n": len(got[name])}
    real = G2.publish
    G2.publish = capture
    try:
        {"b4": G2.b4, "b5": G2.b5, "b6": G2.b6, "b7": G2.b7}[S](cx, True)
    finally:
        G2.publish = real
    return got.get(f"{S}_v1", []), got.get(f"{S}_dec", [])


def tag(line, name, why):
    new = dict(line)
    new.pop("set", None)
    new["revised_by"] = list(line.get("revised_by") or []) + [f"b_{name}"]
    new["revision_note"] = why
    return new


def digest(line):
    """A short digest of a line's truth and question (`canon`)."""
    return hashlib.blake2b(canon(line).encode(), digest_size=16).digest()


def digests(org, sets):
    """{qid: (digest, excluded)} of the latest lines over the sets in their order (kept small: the large sets hold
    millions of lines)."""
    out = {}
    d = L.query_dir(FAMILY, org)
    for s in sets:
        for ln in L.iter_set(L._set_file(d, s)):
            out[ln["qid"]] = (digest(ln), bool(ln.get("excluded")))
    return out


def lines_of(org, sets, qids):
    """{qid: line} over the sets in their order, for the given qids only."""
    out = {}
    d = L.query_dir(FAMILY, org)
    for s in sets:
        for ln in L.iter_set(L._set_file(d, s)):
            if ln["qid"] in qids:
                ln["set"] = s
                out[ln["qid"]] = ln
    return out


def is_correction(org, s):
    """A set that corrects the judge (the harness's revisions of decided rules, or a `--correction`), not a base set
    (a set of version 1, `_dec`, `_v2`, or a regeneration revision)."""
    m = L.manifest(FAMILY, org, s)
    return bool(m.get("rules") or m.get("correction"))


def chain(org, S, decided=False):
    """The corrections of the judge published for set family S, in their publication order: [(set, kind, rules)],
    kind `rules` (the harness's revisions of decided rules, `--rules`) or the name of a `--correction`. The `now` judge
    reads the sets without `_dec`; the decided judge every set, where the `_dec` set of a `--rules` revision applies
    the same rules as its `now` set (counted once)."""
    out, seen = [], set()
    for s in sorted(family_sets(org, S), key=lambda x: L.manifest(FAMILY, org, x)["seq"]):
        if "_dec" in s and not decided:
            continue
        m = L.manifest(FAMILY, org, s)
        if m.get("rules"):
            if m.get("revision") in seen:
                continue
            seen.add(m.get("revision"))
            out.append((s, "rules", list(m["rules"])))
        elif m.get("correction"):
            out.append((s, m["correction"], None))
    return out


class Readers:
    """The truth readers the corrections of the judge read, each made once when first needed."""

    def __init__(self, org, nb=None):
        self.org, self._nb = org, nb

    @functools.cached_property
    def nb(self):
        return self._nb or TB.NamesB(self.org, log=lambda *a: None)

    @functools.cached_property
    def tr(self):
        from idtrack_reproducibility.harness import harness_lib as H
        from idtrack_reproducibility.harness import make_queries as M
        return M.Truth(self.org, H.ORGANISMS[self.org][2])

    @functools.cached_property
    def asm(self):
        from idtrack_reproducibility.harness import harness_lib as H
        return H.gtf_assemblies(self.org)

    @functools.cached_property
    def mapper(self):
        from idtrack_reproducibility.harness import assembly_map
        return assembly_map.Mapper(self.org)


def replay(org, lines, ch, rd):
    """The lines with the corrections `ch` (`chain`) applied in order, as each correction set applied its correction
    to the latest lines when it was published. A correction set that put back lines an earlier correction had revised
    against a decided rule (`places_not_comparable_undone`) needs no step: the corrections leave those lines alone
    (`not_held`, `decided_by_mirbase`)."""
    from idtrack_reproducibility.harness import make_queries as M
    by = {ln["qid"]: ln for ln in lines}
    for _s, kind, rules in ch:
        if kind == "rules":
            by = {ln["qid"]: ln for ln in M.apply_revisions(org, list(by.values()), rd.nb.C, names=rules)}
        elif kind != "places_not_comparable_undone":
            for ln in CORRECTIONS[kind](org, list(by.values()), rd):
                by[ln["qid"]] = ln
    return list(by.values())


def chunks(it, n):
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) >= n:
            yield buf
            buf = []
    if buf:
        yield buf


def differing(old, new):
    """The fields in which two lines differ (their truth and question, `canon`)."""
    a, b = json.loads(canon(old)), json.loads(canon(new))
    out = sorted(k for k in set(a) | set(b) if k != "extra" and a.get(k) != b.get(k))
    xa, xb = json.loads(a.get("extra") or "{}"), json.loads(b.get("extra") or "{}")
    return out + sorted(f"extra.{k}" for k in set(xa) | set(xb) if xa.get(k) != xb.get(k))


# the fields that ask (a change of one of them changes the question, which a stored answer no longer answers)
QUESTION = ("id", "from_release", "to_release", "final_database", "strategy", "from_database", "extra.from_database_rows")


def revise(cx, S, name, why, dry, rd):
    """A regeneration revision of set family S (the module's docstring): every regenerated line, with the corrections
    of the judge published for S replayed on it (`replay`), against the line the published sets give."""
    org = cx.org
    fam = family_sets(org, S)
    if not fam:
        log(org, S, "no published set")
        return
    if f"{S}_{name}" in fam:
        log(org, S, f"{S}_{name} published before: skipped")
        return
    now_sets = [s for s in fam if "_dec" not in s]
    base_sets = [s for s in now_sets if not is_correction(org, s)]
    ch = chain(org, S)
    log(org, S, "base sets", base_sets, "corrections replayed", [(s, k, r) for s, k, r in ch])
    pub, base = digests(org, now_sets), digests(org, base_sets)
    base_lines, dec = regenerate(cx, S)
    changed, new, seen, regen_base, mismatch = [], [], set(), {}, []
    keep_base = S in ("b6", "b7")
    for chunk in chunks(base_lines, 100_000):
        bdig = {ln["qid"]: digest(ln) for ln in chunk}
        if keep_base:
            regen_base.update((ln["qid"], copy.deepcopy(ln)) for ln in chunk)
        for ln in replay(org, chunk, ch, rd):
            q = ln["qid"]
            seen.add(q)
            old = pub.get(q)
            if old is None:
                new.append(ln)
                continue
            if digest(ln) != old[0] or old[1]:
                if not old[1] and base.get(q, (None,))[0] == bdig[q]:
                    mismatch.append(q)  # the regenerated base is the published one, the replay is not the published line
                changed.append(tag(ln, name, why))
    gone = {q for q, (_, exc) in pub.items() if q not in seen and not exc}
    old_lines = lines_of(org, now_sets, {ln["qid"] for ln in changed} | gone)
    for q in sorted(gone):
        changed.append(tag(dict(old_lines[q], excluded=f"no longer asked: {why}"), name, why))
    summarise(org, S, changed, old_lines, set(mismatch), gone, new)
    if mismatch:
        raise SystemExit(f"{org} {S}: {len(mismatch):,} lines whose regenerated base is the published one but whose "
                         "replayed corrections are not the published line: a correction the replay does not reproduce; "
                         "nothing published")
    if changed:
        L.publish(FAMILY, org, f"{S}_{name}", changed, f"{S}: revision {name}: {why}", __file__, revises=now_sets,
                  judge=f"b_{name}", dry=dry,
                  manifest_extra=dict(revision=name, why=why, changed=len(changed) - len(gone), excluded=len(gone),
                                      replayed=[s for s, _, _ in ch]))
    if new:
        L.publish(FAMILY, org, f"{S}_v2", new, f"{S}: lines asked for the first time ({name}: {why})", __file__,
                  judge=f"b_{name}", dry=dry)
    if S in ("b6", "b7"):
        # what the decided judge reads after this revision: the sets in publication order, then `S_<name>`
        dec_sets = [s for s in fam if s not in now_sets]
        order = sorted(fam, key=lambda s: L.manifest(FAMILY, org, s)["seq"])
        pub_dec0 = latest(org, order)
        dec_base = latest(org, [s for s in order if not is_correction(org, s)])
        pub_dec = dict(pub_dec0)
        for ln in changed:
            pub_dec[ln["qid"]] = ln
        chd = chain(org, S, decided=True)
        dec_qids = {ln["qid"] for ln in dec}
        want = sorted(dec_qids | {q for q in pub_dec if pub_dec[q].get("set") in dec_sets})
        start = [copy.deepcopy(ln) for ln in dec] + [copy.deepcopy(regen_base[q]) for q in want
                                                      if q not in dec_qids and q in regen_base]
        sdig = {ln["qid"]: digest(ln) for ln in start}
        target = {ln["qid"]: ln for ln in replay(org, start, chd, rd)}
        out, mism = [], []
        for q in want:
            t = target.get(q)
            if t is None:  # no longer asked
                t = dict(pub_dec[q], excluded=f"no longer asked: {why}")
            if canon(t) != canon(pub_dec[q]) or (t.get("excluded") and not pub_dec[q].get("excluded")):
                if q in dec_base and sdig.get(q) == digest(dec_base[q]) and canon(t) != canon(pub_dec0[q]):
                    mism.append(q)
                out.append(tag(t, f"{name}_dec", why))
        log(org, S, f"decided judge: {len(out):,} lines; replay mismatches {len(mism):,}")
        if mism:
            raise SystemExit(f"{org} {S}: {len(mism):,} decided lines the replay does not reproduce; nothing published")
        if LINES_TO and dry and out:
            with gzip.open(LINES_TO, "at") as fo:
                for t in out:
                    o = pub_dec[t["qid"]]
                    fo.write(json.dumps(dict(set=f"{S}_{name}_dec", qid=t["qid"], kind=t.get("kind"),
                                             stratum=t.get("stratum"), fields=differing(o, t), old=o, new=t),
                                        sort_keys=True, default=str) + "\n")
        if out:
            L.publish(FAMILY, org, f"{S}_{name}_dec", out, f"{S}: revision {name} under the decided judge: {why}", __file__,
                      revises=fam + ([f"{S}_{name}"] if changed else []), judge=f"b_{name}_dec", dry=dry)


def summarise(org, S, changed, old_lines, mismatch, gone, new):
    """Log the revision's lines by kind and stratum and by the fields that changed; with --dry and --lines-to, write
    every changed line beside its published line (for reading against the truth files before publishing)."""
    by_kind, by_field, asks = collections.Counter(), collections.Counter(), collections.Counter()
    fo = gzip.open(LINES_TO, "at") if LINES_TO else None
    for ln in changed:
        q = ln["qid"]
        o = old_lines.get(q)
        f = ["excluded"] if q in gone else differing(o, ln)
        by_kind[(ln.get("kind"), ln.get("stratum"))] += 1
        by_field[",".join(f)] += 1
        if any(x in QUESTION for x in f):
            asks[(ln.get("kind"), ln.get("stratum"))] += 1
        if fo:
            fo.write(json.dumps(dict(set=f"{S}", qid=q, kind=ln.get("kind"), stratum=ln.get("stratum"), fields=f,
                                     mismatch=q in mismatch, old=o, new=ln), sort_keys=True, default=str) + "\n")
    if fo:
        fo.close()
    log(org, S, f"changed {len(changed) - len(gone):,}; no longer asked {len(gone):,}; new {len(new):,}; "
        f"replay mismatches {len(mismatch):,}")
    log(org, S, "by kind and stratum", {f"{k}|{s}": n for (k, s), n in sorted(by_kind.items(), key=str)})
    log(org, S, "by the fields that changed", dict(by_field.most_common()))
    log(org, S, "question changed (ask again)", {f"{k}|{s}": n for (k, s), n in sorted(asks.items(), key=str)})


def revise_rules(org, nb, S, name, rules, dry):
    """`--rules`: the harness's revisions of decided rules applied to the latest `now` lines of set family S."""
    from idtrack_reproducibility.harness import make_queries as M
    fam = family_sets(org, S)
    if not fam:
        log(org, S, "no published set")
        return
    now_sets = [s for s in fam if "_dec" not in s]
    if f"{S}_{name}" in fam:
        log(org, S, f"{S}_{name} published before")
    else:
        _rules_set(org, S, f"{S}_{name}", latest(org, now_sets), now_sets, name, rules, dry, nb)
    dec_sets = [s for s in fam if "_dec" in s]
    if dec_sets and f"{S}_{name}_dec" not in fam:
        # the decided judge reads the `now` sets, then the `_dec` sets: its lines get the same revisions, published last,
        # with every decided line of a qid the `now` revision holds (read in publication order, the `now` revision would
        # otherwise replace the decided truth)
        pub_dec = latest(org, sorted(dec_sets, key=lambda x: L.manifest(FAMILY, org, x)["seq"]))
        now_rev = family_sets(org, S)
        must = {ln["qid"] for ln in L.iter_set(L._set_file(L.query_dir(FAMILY, org), f"{S}_{name}"))} \
            if f"{S}_{name}" in now_rev else set()
        _rules_set(org, S, f"{S}_{name}_dec", pub_dec, now_rev, name, rules, dry, nb, must=must)


def _rules_set(org, S, set_name, pub, revises, name, rules, dry, nb, judge="rules", must=()):
    """Publish `set_name`: the lines of `pub` the revisions change, and every line of `pub` whose qid is in `must`
    (a decided line read after the `now` revision of its qid, revised or not)."""
    from idtrack_reproducibility.harness import make_queries as M
    before = {q: canon(ln) for q, ln in pub.items()}
    lines = [dict(ln) for ln in pub.values() if not ln.get("excluded")]
    for ln in lines:
        ln.pop("set", None)
    out = []
    for ln in M.apply_revisions(org, lines, nb.C, names=rules):
        if canon(ln) != before[ln["qid"]] or ln["qid"] in must:
            out.append(ln)
    by = {}
    for ln in out:
        for r in ln.get("revised_by") or []:
            if r in rules:
                by[r] = by.get(r, 0) + 1
    log(org, S, f"{set_name}: lines {len(lines):,}; revised {len(out):,} {by}")
    if out:
        L.publish(FAMILY, org, set_name, out, f"{S}: the accuracy harness's revisions of decided rules ({', '.join(rules)})"
                  + (" on the decided judge's lines" if set_name.endswith("_dec") else ""),
                  __file__, revises=revises, judge=judge, dry=dry, manifest_extra=dict(revision=name, rules=rules, revised=by))


def not_held(ln):
    """A line whose name a decided rule says is no name at the release read (`mirbase_genes`: a miRBase entry of two
    loci carried by no other database): it expects nothing because nothing is found, not for want of a successor."""
    return "mirbase_genes" in (ln.get("revised_by") or []) and "not found" in (ln.get("note") or "")


def decided_by_mirbase(ln):
    """A line whose reading rule `mirbase_genes` decided (a revision of that decided rule published before)."""
    return "mirbase_genes" in (ln.get("revised_by") or [])


def places_undo(org, S, name, wrong, dry, keep=None):
    """The lines that a correction published as `S_<wrong>` though a decided rule had decided them before it
    (`places_not_comparable`: `not_held`; `held_in_core`: `decided_by_mirbase`): their line before `S_<wrong>`,
    published again as `S_<name>`."""
    keep = keep or not_held
    fam = family_sets(org, S)
    if f"{S}_{name}" in fam or f"{S}_{wrong}" not in fam:
        log(org, S, f"{S}_{name} published before" if f"{S}_{name}" in fam else f"no {S}_{wrong}")
        return
    before = latest(org, [s for s in fam if "_dec" not in s and s != f"{S}_{wrong}"
                          and L.manifest(FAMILY, org, s)["seq"] < L.manifest(FAMILY, org, f"{S}_{wrong}")["seq"]])
    out = []
    for ln in L.iter_set(L._set_file(L.query_dir(FAMILY, org), f"{S}_{wrong}")):
        old = before.get(ln["qid"])
        if old is not None and keep(old):
            new = dict(old, revised_by=list(old.get("revised_by") or []) + [f"{wrong}_undone"])
            new.pop("set", None)
            out.append(new)
    log(org, S, f"{S}_{name}: lines of {S}_{wrong} put back {len(out):,}")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: lines that {S}_{wrong} revised though rule `mirbase_genes` had "
                  "decided them: their truth before it", __file__,
                  revises=[s for s in fam if "_dec" not in s], judge="rules", dry=dry,
                  manifest_extra=dict(revision=name, correction="places_not_comparable_undone"))


def _core_support(org, todo, tr):
    """{qid: support or None} from Ensembl's core records at the release read: per expected gene, the number of name
    sources that give the gene the string there (each database of the core cross-references whose label it is, and the
    display name, `core_display`); for a line asked at a later target, the largest support of the genes whose strict
    walk reaches the expected gene. None where the expected genes do not differ."""
    import gzip
    from idtrack_reproducibility.harness import harness_lib as H
    from idtrack_reproducibility.harness import make_queries as M
    want = collections.defaultdict(set)
    for ln in todo:
        want[ln["from_release"]].add(ln["id"])
    count = collections.defaultdict(lambda: collections.defaultdict(set))  # (release, string) -> gene -> sources
    for R, vs in want.items():
        for sub, col, src in (("core_xref", "label", "db_name"), ("core_display", "label", None)):
            p = os.path.join(H.truth_dir(org), sub, f"{R}.tsv.gz")
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt") as fh:
                head = fh.readline().rstrip("\n").split("\t")
                iv, ig = head.index(col), head.index("gene_id")
                isrc = head.index(src) if src else None
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) > iv and f[iv] in vs:
                        kind = f[isrc] if isrc is not None else sub
                        count[(R, f[iv])][M.gene_base(f[ig])].add(kind)
    out = {}
    for ln in todo:
        R, t, v = ln["from_release"], ln["to_release"], ln["id"]
        by = count.get((R, v), {})
        by_base = {M.gene_base(e): e for e in ln["expect"]}
        sup = {e: 0 for e in ln["expect"]}
        for g, srcs in by.items():
            reach = {g} if t == R else tr.walk_fwd(g, R, t, True, living=True)
            for x in reach:
                if x in by_base:
                    sup[by_base[x]] = max(sup[by_base[x]], len(srcs))
        out[ln["qid"]] = sup if len(set(sup.values())) > 1 else None
    return out


def _generated_support(nb, ln):
    """The support B1's generator gives a line (`gen_b.GenB1.string_lines`: `truth_b.NamesB.support` of the genes whose
    display name the string is at the release read, walked to the target for a later one), or None where the
    release's name sources are not read or the expected genes do not differ."""
    from idtrack_reproducibility.harness import make_queries as M
    R, t, v = ln["from_release"], ln["to_release"], ln["id"]
    by_base = {M.gene_base(e): e for e in ln["expect"]}
    if t == R:
        sup = {e: nb.support(R, b, v) for b, e in by_base.items()}
        if None in sup.values():
            return None
    else:
        allg = set(nb.display(R).get(v, ())) & nb.gset(R)
        if not allg or nb.support(R, next(iter(sorted(allg))), v) is None:
            return None
        sup = {e: 0 for e in ln["expect"]}
        for g in sorted(allg):
            n = nb.support(R, g, v) or 0
            for x in nb.walk(g, R, t)[0]:
                if x in by_base:
                    sup[by_base[x]] = max(sup[by_base[x]], n)
    return sup if len(set(sup.values())) > 1 else None


def _support(org, todo, nb=None):
    """{qid: support or None} of B1 lines of a name several genes carry (rule `support_ranks_shared_name`), each by the
    procedure that gave its published support, of name sources alone: a synonym is none (rule `names_before_synonyms`:
    idtrack reads a string as the names that carry it and counts the databases attaching those,
    `Candidate.entry_support`). At a release with BioMart's tables, the generator's support (`_generated_support`:
    BioMart's name sources, which it counted with BioMart's synonyms); where that leaves the genes equal, and at the
    other releases, the core records' (`_core_support`, the correction `support_from_xrefs`, which counted Ensembl's
    core synonyms). The generator's support at the other releases never counted a synonym: such a line keeps it
    (`_KEEP`)."""
    from idtrack_reproducibility.families.b import truth_b as TBm
    nb = nb or TBm.NamesB(org, log=lambda *a: None)
    out = {}
    for ln in sorted(todo, key=lambda r: r["from_release"]):  # the release's tables are read in order
        by_xrefs = "support_from_xrefs" in (ln.get("revised_by") or [])
        if ln["from_release"] in nb.bm_rels and not by_xrefs:
            out[ln["qid"]] = _generated_support(nb, ln)
        elif ln.get("support") and not by_xrefs:
            out[ln["qid"]] = _KEEP
        else:
            out[ln["qid"]] = None
    rest = [ln for ln in todo if out[ln["qid"]] is None]
    out.update(_core_support(org, rest, nb.tr))
    return out


_KEEP = object()  # a line whose support no synonym entered: kept as it is


def _support_todo(lines, with_support):
    return [ln for ln in lines if ln.get("kind") in ("b1_name_at_release", "b1_name_to_current")
            and ln.get("rule") == "set_exact" and len(ln.get("expect") or []) > 1 and not ln.get("excluded")
            and (with_support or not ln.get("support")) and not ln.get("open") and ln.get("from_release") is not None]


def support_xrefs_lines(org, lines, rd):
    """The lines `support_from_xrefs` revises (the module's docstring), revised."""
    todo = _support_todo(lines, False)
    sup = _core_support(org, todo, rd.tr) if todo else {}
    out = []
    for ln in todo:
        if sup[ln["qid"]]:
            new = dict(ln, support=sup[ln["qid"]], revised_by=list(ln.get("revised_by") or []) + ["support_from_xrefs"])
            new.pop("set", None)
            out.append(new)
    return out


def support_names_lines(org, lines, rd, moved=None):
    """The lines `support_names` revises (the module's docstring), revised; `moved` counts how their support moved."""
    todo = _support_todo(lines, True)
    sup = _support(org, todo, rd.nb) if todo else {}
    out = []
    for ln in todo:
        new_sup = sup[ln["qid"]]
        if new_sup is _KEEP or (new_sup or None) == (ln.get("support") or None):
            continue
        new = dict(ln, revised_by=list(ln.get("revised_by") or []) + ["support_names"])
        new.pop("set", None)
        if new_sup:
            new["support"] = new_sup
        else:
            new.pop("support", None)
        if moved is not None:
            moved["set" if new_sup and not ln.get("support") else ("dropped" if not new_sup else "changed")] += 1
        out.append(new)
    return out


NOT_HELD_NOTE = "no Ensembl table of any release holds the string"
LINES_TO = None  # with --dry: the revision lines of a correction are written here, for reading before publishing


def held_lines(org, lines, rd):
    """The lines `held_in_core` revises (the module's docstring), revised."""
    from idtrack_reproducibility.harness import make_queries as M
    todo = [ln for ln in lines if (ln.get("note") or "").startswith(NOT_HELD_NOTE) and not ln.get("excluded")
            and not ln.get("expect") and not decided_by_mirbase(ln)]
    if not todo:
        return []
    tr = rd.tr
    found = M.core_readings(org, tr.rels, {ln["id"] for ln in todo}, casefold=False)
    out = []
    for ln in todo:
        at = found.get(ln["id"], {})
        t = ln["to_release"]
        held = sorted({r for r, _ in at if r <= t})
        if not held:
            continue
        R = ln["from_release"] if ln.get("from_release") in held else (t if t in held else held[-1])
        genes = at.get((R, "name")) or at.get((R, "synonym"))  # rule `names_before_synonyms`
        E = set().union(*(tr.walk_series(g, R, [t], True, living=True)[t] for g in genes))
        B = E | set().union(*(tr.walk_series(g, R, [t], False, living=True, resurrect=True)[t] for g in genes))
        Ev, Bv = M._vers(types.SimpleNamespace(tr=tr), E, t, ln["cmp"]), M._vers(types.SimpleNamespace(tr=tr), B, t, ln["cmp"])
        new = dict(ln, expect=Ev, broad=Bv if set(Bv) != set(Ev) else None, read_at=R,
                   rule="set_exact" if len(genes) > 1 and len(Ev) > 1 else "set",
                   note=f"{ln['note']}; held at {R} in Ensembl's core records (the corrected judge, held_in_core)",
                   revised_by=list(ln.get("revised_by") or []) + ["held_in_core"])
        if not E and len(genes) == 1:
            loc = tr.locus(next(iter(genes)), R)
            if loc:
                new["locus"] = loc
        new.pop("set", None)
        out.append(new)
    return out


def places_lines(org, lines, rd):
    """The lines `places_not_comparable` revises (the module's docstring), revised."""
    out = []
    for ln in lines:
        R, t = ln.get("from_release"), ln.get("to_release")
        if ln.get("kind") != "b1_name_to_current" or ln.get("excluded") or ln.get("expect") or ln.get("broad") or \
                ln.get("locus") or ln.get("open") or R is None or rd.asm.get(R) == rd.asm.get(t) or \
                rd.mapper.path(rd.asm.get(R), rd.asm.get(t)) is not None or not_held(ln):
            continue
        new = dict(ln, open=f"Ensembl's history gives the genes that carry the name at {R} no successor at {t}, and "
                            f"{rd.asm.get(R)} and {rd.asm.get(t)} have no assembly chain: an answer cannot be compared by place",
                   revised_by=list(ln.get("revised_by") or []) + ["places_not_comparable"])
        new.pop("set", None)
        out.append(new)
    return out


# the corrections of the judge (`--correction`), each a function of the latest lines: the lines it revises, revised
CORRECTIONS = {"places_not_comparable": places_lines, "support_from_xrefs": support_xrefs_lines,
               "support_names": support_names_lines, "held_in_core": held_lines}
DESCRIPTIONS = {
    "places_not_comparable": ("names of several genes without successor, read on an assembly without a chain to the "
                              "target's: open (the corrected judge, places_not_comparable)"),
    "support_from_xrefs": ("support per gene from Ensembl's records at the release read, for rule "
                           "`support_ranks_shared_name` (the corrected judge, support_from_xrefs)"),
    "support_names": ("support per gene from Ensembl's name sources at the release read, synonyms not counted, for rule "
                      "`support_ranks_shared_name` (the corrected judge, support_names)"),
    "held_in_core": ("strings the truth found in no Ensembl table, held by Ensembl's core records of a release (the "
                     "corrected judge, held_in_core)"),
}


def correct(org, S, name, correction, dry):
    """`--correction <correction>` on set family S: the lines it revises among the latest `now` lines, published as
    `S_<name>` (judge `78`)."""
    fam = family_sets(org, S)
    if not fam or f"{S}_{name}" in fam:
        log(org, S, "no set" if not fam else f"{S}_{name} published before")
        return
    pub = latest(org, [s for s in fam if "_dec" not in s])
    moved = collections.Counter()
    fn = CORRECTIONS[correction]
    rd = Readers(org)
    out = fn(org, list(pub.values()), rd, moved) if correction == "support_names" else fn(org, list(pub.values()), rd)
    log(org, S, f"{S}_{name}: revised {len(out):,}", dict(moved) if moved else "")
    if dry and LINES_TO and out:
        with gzip.open(LINES_TO, "at") as fo:
            for r in out:
                fo.write(json.dumps(r, sort_keys=True, default=str) + "\n")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: {DESCRIPTIONS[correction]}", __file__,
                  revises=[s for s in fam if "_dec" not in s], judge="78", dry=dry,
                  manifest_extra=dict(revision=name, correction=correction))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--why", default=None, help="what changed in the truth (required without --rules)")
    ap.add_argument("--sets", default="b4,b5,b6,b7,b8,b1")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--rules", default=None, help="the harness's revisions of decided rules to apply (no regeneration)")
    ap.add_argument("--correction", default=None, choices=["places_not_comparable", "places_undo", "support_from_xrefs",
                                                           "support_names", "held_in_core"],
                    help="a correction of the judge (`places_undo` with --undo <name of the places set>)")
    ap.add_argument("--undo", default=None)
    ap.add_argument("--lines-to", default=None, help="with --dry: write the revision lines here (a regeneration "
                                                      "revision: each beside its published line)")
    a = ap.parse_args()
    global LINES_TO
    LINES_TO = a.lines_to
    if a.correction:
        for S in a.sets.split(","):
            if a.correction == "places_undo":
                places_undo(a.organism, S, a.name, a.undo, a.dry,
                            keep=decided_by_mirbase if a.undo.startswith("held") else not_held)
            else:
                correct(a.organism, S, a.name, a.correction, a.dry)
        return
    if a.rules:
        nb = TB.NamesB(a.organism, log=log)
        for S in a.sets.split(","):
            revise_rules(a.organism, nb, S, a.name, a.rules.split(","), a.dry)
        return
    if not a.why:
        ap.error("--why is required without --rules")
    sets = a.sets.split(",")
    # B9 reads its own truth (gen_b9): the context of the other sets (the history, BioMart) only where they are asked
    cx = G2.Ctx(a.organism) if set(sets) - {"b9"} else types.SimpleNamespace(org=a.organism, nb=types.SimpleNamespace(_walk={}))
    rd = Readers(a.organism, nb=getattr(cx, "nb", None) if set(sets) - {"b9"} else None)
    for S in sets:
        t0 = time.time()
        revise(cx, S, a.name, a.why, a.dry, rd)
        cx.nb._walk.clear()
        log(a.organism, S, f"{time.time() - t0:.0f}s")
    _ = TB


if __name__ == "__main__":
    main()
