"""Revision sets of family B: every set regenerated with the current truth and compared, qid by qid, with what is
published; only the lines whose truth changed are published, under their qids (F0's rule: a published set is never
rewritten).

Usage: revise_b.py --organism <o> --name r1 --why "<what changed in the truth>" [--sets b1,b4,b5,b6,b7,b8,b9] [--dry]

For each set family S (`b1`, ...):
- `S_<name>`: the lines of the `now` judge (every set of S whose name holds no `_dec`) whose regenerated line differs
  (`revised_by` names the revision); a published qid that the generator no longer makes is published again with
  `excluded` (left out of the scoring, counted in run.json);
- `S_<name>_dec` (B6, B7): the lines the `decided` judge must read after `S_<name>` (sets in order: `S_v1`, `S_dec`,
  `S_<name>`, `S_<name>_dec`): every qid whose decided line (its regenerated `_dec` line, else its regenerated line)
  differs from what that order gives;
- `S_v2`: lines with qids not published before (none expected; a new set, never a revision).

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


def canon(line):
    """A line as compared between its published and its regenerated form: the truth and the question, without the
    bookkeeping of the revisions that published it (`revised_by`, `revision_note`), so that a later revision does not
    publish again the lines an earlier one revised."""
    r = L.to_row({k: v for k, v in line.items() if k not in ("set", "revised_by", "revision_note")})
    r.pop("h", None)
    r.pop("revised_by", None)
    if r.get("extra"):
        x = json.loads(r["extra"])
        x.pop("set", None)
        r["extra"] = json.dumps(x, sort_keys=True) if x else None
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


def revise(cx, S, name, why, dry):
    org = cx.org
    fam = family_sets(org, S)
    if not fam:
        log(org, S, "no published set")
        return
    if f"{S}_{name}" in fam:
        log(org, S, f"{S}_{name} published before: skipped")
        return
    now_sets = [s for s in fam if "_dec" not in s]
    pub_now = latest(org, now_sets)
    base, dec = regenerate(cx, S)
    changed, new, seen, regen_base = [], [], set(), {}
    keep_base = bool(dec) or S in ("b6", "b7")
    for ln in base:
        q = ln["qid"]
        seen.add(q)
        if keep_base:
            regen_base[q] = ln
        old = pub_now.get(q)
        if old is None:
            new.append(ln)
        elif canon(old) != canon(ln) or old.get("excluded"):
            changed.append(tag(ln, name, why))
    gone = [q for q in pub_now if q not in seen and not pub_now[q].get("excluded")]
    for q in gone:
        changed.append(tag(dict(pub_now[q], excluded=f"no longer asked: {why}"), name, why))
    log(org, S, f"regenerated {len(seen):,}; changed {len(changed) - len(gone):,}; no longer asked {len(gone):,}; new {len(new):,}")
    if changed:
        L.publish(FAMILY, org, f"{S}_{name}", changed, f"{S}: revision {name}: {why}", __file__, revises=now_sets,
                  judge=f"b_{name}", dry=dry,
                  manifest_extra=dict(revision=name, why=why, changed=len(changed) - len(gone), excluded=len(gone)))
    if new:
        L.publish(FAMILY, org, f"{S}_v2", new, f"{S}: lines asked for the first time ({name}: {why})", __file__,
                  judge=f"b_{name}", dry=dry)
    if S in ("b6", "b7"):
        # what the decided judge reads after this revision: v1 < dec < r < r_dec
        dec_sets = [s for s in fam if s not in now_sets]
        pub_dec = latest(org, sorted(set(now_sets) | set(dec_sets), key=lambda s: L.manifest(FAMILY, org, s)["seq"]))
        for ln in changed:
            pub_dec[ln["qid"]] = ln
        regen_dec = {ln["qid"]: ln for ln in dec}
        out = []
        for q in sorted(set(regen_dec) | {q for q in pub_dec if any(pub_dec[q].get("set") == s for s in dec_sets)}):
            target = regen_dec.get(q) or regen_base.get(q)
            if target is None:  # no longer asked
                target = dict(pub_dec[q], excluded=f"no longer asked: {why}")
            if canon(target) != canon(pub_dec[q]) or (target.get("excluded") and not pub_dec[q].get("excluded")):
                out.append(tag(target, f"{name}_dec", why))
        log(org, S, f"decided judge: {len(out):,} lines")
        if out:
            L.publish(FAMILY, org, f"{S}_{name}_dec", out, f"{S}: revision {name} under the decided judge: {why}", __file__,
                      revises=fam + ([f"{S}_{name}"] if changed else []), judge=f"b_{name}_dec", dry=dry)


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


def _support_todo(pub, with_support):
    return [ln for ln in pub.values() if ln.get("kind") in ("b1_name_at_release", "b1_name_to_current")
            and ln.get("rule") == "set_exact" and len(ln.get("expect") or []) > 1 and not ln.get("excluded")
            and (with_support or not ln.get("support")) and not ln.get("open") and ln.get("from_release") is not None]


def support_from_xrefs(org, S, name, dry):
    """`--correction support_from_xrefs` on set family S (the module's docstring)."""
    fam = family_sets(org, S)
    if not fam or f"{S}_{name}" in fam:
        log(org, S, "no set" if not fam else f"{S}_{name} published before")
        return
    pub = latest(org, [s for s in fam if "_dec" not in s])
    todo = _support_todo(pub, False)
    from idtrack_reproducibility.harness import harness_lib as H
    from idtrack_reproducibility.harness import make_queries as M
    sup = _core_support(org, todo, M.Truth(org, H.ORGANISMS[org][2]))
    out = []
    for ln in todo:
        if sup[ln["qid"]]:
            new = dict(ln, support=sup[ln["qid"]], revised_by=list(ln.get("revised_by") or []) + ["support_from_xrefs"])
            new.pop("set", None)
            out.append(new)
    log(org, S, f"{S}_{name}: revised {len(out):,} of {len(todo):,}")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: support per gene from Ensembl's records at the release read, for "
                  "rule `support_ranks_shared_name` (the corrected judge, support_from_xrefs)", __file__,
                  revises=[s for s in fam if "_dec" not in s], judge="78", dry=dry,
                  manifest_extra=dict(revision=name, correction="support_from_xrefs"))


def support_names(org, S, name, dry):
    """`--correction support_names` on set family S (the module's docstring)."""
    fam = family_sets(org, S)
    if not fam or f"{S}_{name}" in fam:
        log(org, S, "no set" if not fam else f"{S}_{name} published before")
        return
    pub = latest(org, [s for s in fam if "_dec" not in s])
    todo = _support_todo(pub, True)
    sup = _support(org, todo)
    out, moved = [], collections.Counter()
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
        moved["set" if new_sup and not ln.get("support") else ("dropped" if not new_sup else "changed")] += 1
        out.append(new)
    log(org, S, f"{S}_{name}: revised {len(out):,} of {len(todo):,}", dict(moved))
    if LINES_TO and dry and out:
        import gzip
        with gzip.open(LINES_TO, "at") as fo:
            for r in out:
                fo.write(json.dumps(r, default=str) + "\n")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: support per gene from Ensembl's name sources at the release "
                  "read, synonyms not counted, for rule `support_ranks_shared_name` (the corrected judge, support_names)",
                  __file__, revises=[s for s in fam if "_dec" not in s], judge="78", dry=dry,
                  manifest_extra=dict(revision=name, correction="support_names"))


NOT_HELD_NOTE = "no Ensembl table of any release holds the string"
LINES_TO = None  # with --dry: the revision lines of a correction are written here, for reading before publishing


def held_in_core(org, S, name, dry):
    """`--correction held_in_core` on set family S (the module's docstring)."""
    from idtrack_reproducibility.harness import harness_lib as H
    from idtrack_reproducibility.harness import make_queries as M
    fam = family_sets(org, S)
    if not fam or f"{S}_{name}" in fam:
        log(org, S, "no set" if not fam else f"{S}_{name} published before")
        return
    pub = latest(org, [s for s in fam if "_dec" not in s])
    todo = [ln for ln in pub.values() if (ln.get("note") or "").startswith(NOT_HELD_NOTE) and not ln.get("excluded")
            and not ln.get("expect") and not decided_by_mirbase(ln)]
    tr = M.Truth(org, H.ORGANISMS[org][2])
    found = M.core_readings(org, tr.rels, {ln["id"] for ln in todo}, casefold=False) if todo else {}
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
    log(org, S, f"{S}_{name}: revised {len(out):,} of {len(todo):,}")
    if dry and LINES_TO and out:
        import gzip
        with gzip.open(LINES_TO, "at") as fo:
            for r in out:
                fo.write(json.dumps(r, sort_keys=True, default=str) + "\n")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: strings the truth found in no Ensembl table, held by Ensembl's "
                  "core records of a release (the corrected judge, held_in_core)", __file__,
                  revises=[s for s in fam if "_dec" not in s], judge="78", dry=dry,
                  manifest_extra=dict(revision=name, correction="held_in_core"))


def places_not_comparable(org, S, name, dry):
    """`--correction places_not_comparable` on set family S (the module's docstring)."""
    from idtrack_reproducibility.harness import assembly_map
    from idtrack_reproducibility.harness import harness_lib as H
    fam = family_sets(org, S)
    if not fam or f"{S}_{name}" in fam:
        log(org, S, "no set" if not fam else f"{S}_{name} published before")
        return
    asm, mapper = H.gtf_assemblies(org), assembly_map.Mapper(org)
    pub = latest(org, [s for s in fam if "_dec" not in s])
    out = []
    for ln in pub.values():
        R, t = ln.get("from_release"), ln.get("to_release")
        if ln.get("kind") != "b1_name_to_current" or ln.get("excluded") or ln.get("expect") or ln.get("broad") or \
                ln.get("locus") or ln.get("open") or R is None or asm.get(R) == asm.get(t) or \
                mapper.path(asm.get(R), asm.get(t)) is not None or not_held(ln):
            continue
        new = dict(ln, open=f"Ensembl's history gives the genes that carry the name at {R} no successor at {t}, and "
                            f"{asm.get(R)} and {asm.get(t)} have no assembly chain: an answer cannot be compared by place",
                   revised_by=list(ln.get("revised_by") or []) + ["places_not_comparable"])
        new.pop("set", None)
        out.append(new)
    log(org, S, f"{S}_{name}: revised {len(out):,}")
    if out:
        L.publish(FAMILY, org, f"{S}_{name}", out, f"{S}: names of several genes without successor, read on an assembly "
                  "without a chain to the target's: open (the corrected judge, places_not_comparable)", __file__,
                  revises=[s for s in fam if "_dec" not in s], judge="78", dry=dry,
                  manifest_extra=dict(revision=name, correction="places_not_comparable"))


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
    ap.add_argument("--lines-to", default=None, help="with --dry: write the revision lines of a correction here")
    a = ap.parse_args()
    global LINES_TO
    LINES_TO = a.lines_to
    if a.correction:
        for S in a.sets.split(","):
            if a.correction == "places_undo":
                places_undo(a.organism, S, a.name, a.undo, a.dry,
                            keep=decided_by_mirbase if a.undo.startswith("held") else not_held)
            elif a.correction == "support_from_xrefs":
                support_from_xrefs(a.organism, S, a.name, a.dry)
            elif a.correction == "support_names":
                support_names(a.organism, S, a.name, a.dry)
            elif a.correction == "held_in_core":
                held_in_core(a.organism, S, a.name, a.dry)
            else:
                places_not_comparable(a.organism, S, a.name, a.dry)
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
    for S in sets:
        t0 = time.time()
        revise(cx, S, a.name, a.why, a.dry)
        cx.nb._walk.clear()
        log(a.organism, S, f"{time.time() - t0:.0f}s")
    _ = TB


if __name__ == "__main__":
    main()
