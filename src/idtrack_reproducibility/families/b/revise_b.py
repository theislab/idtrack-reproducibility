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
"""
import argparse
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--why", default=None, help="what changed in the truth (required without --rules)")
    ap.add_argument("--sets", default="b4,b5,b6,b7,b8,b1")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--rules", default=None, help="the harness's revisions of decided rules to apply (no regeneration)")
    a = ap.parse_args()
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
