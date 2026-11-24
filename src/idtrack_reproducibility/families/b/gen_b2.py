"""Family B's sets B4 (previous symbols), B5 (synonyms and aliases), B6 (forms) and B7 (spreadsheet and typography).

Usage: gen_b2.py --organism <o> [--sets b4,b5,b6,b7] [--dry]

Two truths (DESIGN principle 2): every line's `expect` is **Ensembl's** (the GTF's display names, BioMart's and the
core xrefs' name sources and synonyms at the release, the history's walks: as B1 judges, `gen_b.GenB1`); the
**authority's** gene (HGNC, MGI, VGNC, NCBI: the record's own Ensembl id, at 115) is the line's `authority` field,
never a judge: the report gives the share where idtrack's answer is the authority's gene and the share where Ensembl
and the authority differ (a limitation of Ensembl's data, not an answer of idtrack).

- **B4** (`b4_v1`): every previous symbol (HGNC `prev_symbol` and the withdrawn list's symbols merged into an approved
  one; MGI's withdrawn gene symbols with their current marker; VGNC `prev_symbol`), asked without a release
  (`b4_prev_default`: rule `start_release`, read at 115: names before synonyms, rules `names_before_synonyms` and `name_before_synonym`; a string that is no exact
  name or synonym of 115 is read at the latest release whose Ensembl records hold it, its genes walked to 115
  (`Ctx.read_where_held`, where the records of every release are read; else `open` with `open_within` = the walks to
  115 of every gene that ever carried it, so an answer outside every gene linked to it is wrong) and from the last release whose GTF has it as a display name
  (`b4_prev_from`: B1's judge of that string at that release, to 115).
- **B5** (`b5_v1`; `b5_v2`: the synonyms Ensembl's records of every cross-reference add at those releases,
  `name_src: core_synonyms`): every BioMart synonym of every BioMart release R, read at R, to R (`b5_syn_atR`: the harness's
  `synonym_to_gene` rule: names before synonyms, rule `name_before_synonym`) and to 115 (`b5_syn_to115`: their walks); every
  authority alias (HGNC `alias_symbol`, MGI synonyms of current gene markers, VGNC `alias_symbol`, NCBI `Synonyms`
  where gene_info is in the store) without a release (`b5_alias_<authority>`, judged as B4's default lines).
- **B6** (`b6_v1`, with the revision set `b6_dec` for decisions not built yet): every display name of 115 in lower,
  upper and title case and with white space around it (`b6_case`, `b6_space`: exhaustive); a seeded sample of 20,000
  (name, form) pairs of separators changed, non-breaking space, zero-width space, byte-order mark, CR/LF and
  full-width letters (`b6_sep`, `b6_invisible`, `b6_fullwidth`). Truth: the original's genes (`set`) for what the docs
  of the reference trees promise (rule `spelling_steps`: case and the dash/underscore swap; white space trimmed); `safe` (the
  gene or nothing) for forms no built rule reads; a form that is itself an exact name or synonym of another gene at
  115 is that string's question (left out); a form that equals another gene's string in another letter case is `open`
  within the original's genes and those. `b6_dec`: rule `unicode_normalised` (NFKC, zero-width and BOM removed: the gene) and,
  in mouse, rule `single_case_symbols` (a string of one letter case that is the MGI symbol of exactly one gene without case:
  that gene).
- **B7** (`b7_v1`, with `b7_dec` for Q1 of DESIGN section 7): every symbol of a date-like family (a month name or
  abbreviation and a day: `SEPT7`, `MARCH1`, `DEC1`, `OCT4`, ...) of any release or authority, in the forms Excel
  writes (`7-Sep`, `07-Sep`, `Sep-07`, `2007-09-07`, `9/7/2007`, the serial `39332`); Greek letters and their
  spellings of a curated list (`TNFα`, `IL-1β`, `IFN-γ`, `TGF-β1`, ...); suffixes (`-AS1` on a gene without one,
  `_at`, `_s_at`). `safe` (Q1 not built: the symbol's genes or nothing right, another gene wrong); `b7_dec` (Q1
  decided: read where it resolves to exactly one gene, else nothing): `set` with that gene, or nothing (the suffixes
  stay `safe`: Q1 does not cover them). A string that no table of any release holds expects nothing.
"""
import argparse
import collections
import datetime
import os
import random
import re
import sys
import time
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families.b import gen_b as GB
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.b import truth_b as TB

FAMILY, T, SEED = "b", TB.T, 20261119
log = GB.log


class Ctx:
    def __init__(self, org):
        self.org = org
        self.skip_published = False  # the generator's own reruns skip a set published before (main sets it)
        self.nb = TB.NamesB(org, log=log)
        self.g1 = GB.GenB1(self.nb)
        self.GT = set(self.nb.genes(T))
        self.disp_T = self.nb.display(T)
        self.src_T = self.nb.sources(T)[0]
        self.syn_T = self.nb.synonyms(T) or {}
        lower = collections.defaultdict(set)
        for d in (self.disp_T, self.src_T, self.syn_T):
            for v, gs in d.items():
                lower[v.lower()] |= gs
        self.lower_T = lower

    def any_T(self, v):
        """Genes that carry v at 115 as display name, in another name source or as a synonym (the harness's `broad`
        of a current name)."""
        return (set(self.disp_T.get(v, ())) | set(self.src_T.get(v, ())) | set(self.syn_T.get(v, ()))) & self.GT

    def exact_T(self, v):
        """(name genes, synonym genes) of string v at 115, exactly as written."""
        names = (set(self.disp_T.get(v, ())) | set(self.src_T.get(v, ()))) & self.GT
        return names, set(self.syn_T.get(v, ())) & self.GT

    def default_line(self, qid, stratum, v, unit, authority=None, auth_db=None):
        """A string asked without a release, read at 115 (rule `start_release`): names before synonyms (rules `names_before_synonyms` and `name_before_synonym`)."""
        nb, tr, GT = self.nb, self.nb.tr, self.GT
        disp = set(self.disp_T.get(v, ())) & GT
        names, syn = self.exact_T(v)
        E = disp or names or syn
        rec = dict(qid=qid, family=FAMILY, kind="b_string_default", stratum=stratum, unit=unit, id=v, from_release=None,
                   to_release=T, final_database=None, strategy="both", cmp="versioned", locus=None, broad=None)
        if authority is not None:
            rec["authority"] = tr.vset({g for g in authority if g in GT}, T)
            rec["authority_db"] = auth_db
        if E:
            prim = {g for g in E if nb.region(g, T) != "patch"}
            Ex, apart = (prim, E - prim) if prim else (E, set())
            B = names | syn | E
            rec.update(rule=GB.rule_of(len(Ex)), expect=tr.vset(Ex, T))
            if set(B) != set(Ex):
                rec["broad"] = tr.vset(B, T)
            if apart:
                rec["apart"] = tr.vset(apart, T)
            if len(Ex) > 1:
                s = {tr.vid(g, T): nb.support(T, g, v) for g in Ex}
                if len(set(s.values())) > 1:
                    rec["support"] = s
            return rec
        # no exact name or synonym of 115: read at the latest release where it is active, which the tables do not fix
        ever = set(nb.ever().get(v, ()))
        within = set()
        for R in nb.rels:
            G = nb.genes(R)
            here = {g for g in ever if g in G}
            if here:
                within |= nb.walks(here, R)[1] if R < T else here
        within |= set(self.lower_T.get(v.lower(), ())) & GT
        if not within:  # no table of any release holds the string, in any letter case: nothing is the answer
            rec.update(rule="set", expect=[], note="no Ensembl table of any release holds the string")
            return rec
        rec.update(rule="set", expect=[], open="no exact name or synonym of release 115: read at the latest release where "
                   "it is active (rule `start_release`), which the truth's tables do not fix", open_within=tr.vset(within & GT, T))
        return rec

    def records_of_every_release(self):
        """Whether Ensembl's records of every release of the window are read (`core_xref`, `core_display`,
        `core_synonyms` of each release in the truth store)."""
        d = H.truth_dir(self.org)
        return all(os.path.exists(os.path.join(d, sub, f"{R}.tsv.gz")) for R in self.nb.rels
                   for sub in ("core_xref", "core_display", "core_synonyms"))

    def read_where_held(self, lines):
        """Rule `start_release` on the lines asked without a release whose string no name or synonym of 115 holds
        (`default_line`'s last branches), where Ensembl's records of every release are read
        (`records_of_every_release`): the string is read at the latest release whose records hold it (a label of any
        database of the core cross-references, the display name, or a synonym of any cross-reference;
        `make_queries.core_readings`), names before synonyms (rule `names_before_synonyms`), and its genes walked to
        115 (the strict walk with living ids `expect`, the broad walk `broad`; one gene without a successor holds its
        place, `locus`), as the corrected judge's `held_in_core` reads the strings its tables had not found. The
        line names the release read (`read_at`). A string that a gene of 115 carries in another letter case (idtrack's
        spelling steps may read it there first), one that no record holds as written, and every such line without
        the records of every release, stays as `default_line` made it: open within every gene that ever carried the
        string (the release idtrack reads it at cannot be fixed), or nothing where no table holds it. Lines are
        revised in place; returns how many."""
        todo = [ln for ln in lines if ln.get("kind") == "b_string_default" and not ln.get("expect")
                and (ln.get("open") or ln.get("note") == "no Ensembl table of any release holds the string")
                and not set(self.lower_T.get(ln["id"].lower(), ())) & self.GT]
        if not todo or not self.records_of_every_release():
            return 0
        tr, rels = self.nb.tr, self.nb.rels
        found = GB.M.core_readings(self.org, rels, {ln["id"] for ln in todo}, casefold=False)
        n = 0
        for ln in todo:
            at = found.get(ln["id"], {})
            held = sorted({r for r, _ in at if r <= T})
            if not held:  # no record holds the string as written: the line stays as `default_line` made it
                continue
            R = held[-1]
            genes = at.get((R, "name")) or at.get((R, "synonym"))  # rule `names_before_synonyms`
            E = set().union(*(tr.walk_series(g, R, [T], True, living=True)[T] for g in genes))
            Bw = E | set().union(*(tr.walk_series(g, R, [T], False, living=True, resurrect=True)[T] for g in genes))
            Ev, Bv = tr.vset(E, T), tr.vset(Bw, T)
            ln.pop("open", None)
            ln.pop("open_within", None)
            ln.update(expect=Ev, broad=Bv if set(Bv) != set(Ev) else None, read_at=R,
                      rule="set_exact" if len(genes) > 1 and len(Ev) > 1 else "set",
                      note=f"no exact name or synonym of release 115: read at {R}, the latest release whose records hold it "
                           "(rule `start_release`)")
            if not E and len(genes) == 1:
                loc = tr.locus(next(iter(genes)), R)
                if loc:
                    ln["locus"] = loc
            n += 1
        return n


# ---------------- B4 ----------------

def authority_prev(org, nb):
    """{previous symbol: (authority genes (base ids), authority db)} of the organism's authority."""
    out = collections.defaultdict(lambda: [set(), None])
    h = TB.hgnc(org)
    if h:
        by_id = {r["id"]: r for r in h}
        for r in h:
            for p in r["prev"]:
                out[p][0] |= {r["ensembl"]} if r["ensembl"] else set()
                out[p][1] = "HGNC"
        for r in TB.read_tsv(os.path.join(H.truth_dir(org), "hgnc", "withdrawn.txt")):
            if r["STATUS"] != "Merged/Split":
                continue
            for rep in r["MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)"].split(","):
                hid = rep.strip().split("|")[0]
                t = by_id.get(hid)
                if t:
                    out[r["WITHDRAWN_SYMBOL"]][0] |= {t["ensembl"]} if t["ensembl"] else set()
                    out[r["WITHDRAWN_SYMBOL"]][1] = "HGNC"
    m = TB.mgi(org)
    if m:
        for r in m["withdrawn"]:
            out[r["symbol"]][0] |= set(r["ensembl"])
            out[r["symbol"]][1] = "MGI"
    v = TB.vgnc(org)
    if v:
        for r in v:
            for p in r["prev"]:
                out[p][0] |= {r["ensembl"]} if r["ensembl"] else set()
                out[p][1] = "VGNC"
    return {k: (a, b) for k, (a, b) in out.items() if TB.clean(k)}


def b4(cx, dry):
    nb, g1 = cx.nb, cx.g1
    prev = authority_prev(cx.org, nb)
    log(cx.org, "B4 previous symbols", len(prev))
    last = {}
    for R in nb.rels:
        for v in nb.display(R):
            if v in prev:
                last[v] = R
    out = []
    for p in sorted(prev):
        a, db = prev[p]
        out.append(cx.default_line(f"b4|{p}|default|{T}", "b4_prev_default", p, p, authority=a, auth_db=db))
    by_rel = collections.defaultdict(list)
    for p, R in last.items():
        if R < T:
            by_rel[R].append(p)
    for R in sorted(by_rel):
        for p in sorted(by_rel[R]):
            for ln in g1.string_lines(p, R, prefix="b4", at=False, to=True, strata=lambda st, moved: "b4_prev_from"):
                a, db = prev[p]
                ln["authority"] = nb.tr.vset({g for g in a if g in cx.GT}, T)
                ln["authority_db"] = db
                out.append(ln)
        nb._walk.clear()
    log(cx.org, "B4 lines without a release read where Ensembl's records hold them", cx.read_where_held(out))
    return publish(cx, "b4_v1", out, "B4 previous symbols of the authorities: without a release (read at 115) and from "
                   "the last release whose GTF displays them", dry)


# ---------------- B5 ----------------

def b5(cx, dry):
    nb, tr, org = cx.nb, cx.nb.tr, cx.org
    out = []
    for R in sorted(nb.bm_rels):
        if R not in nb.tr.genes:
            continue
        syn = nb.synonyms(R)
        if not syn:
            continue
        disp, src, G = nb.display(R), nb.sources(R)[0], nb.genes(R)
        prim = {g for g, r in G.items() if r["region"] != "patch"}
        bm_syn = (nb.C.bm(R, "synonym", 1) or (None, {}))[1]  # BioMart's synonyms (of display xrefs): the rest is the core's
        n0 = len(out)
        for s in sorted(syn):
            if not TB.clean(s) or s in disp:
                continue
            sg = syn[s] & prim
            if not sg:
                continue
            names = (set(disp.get(s, ())) | set(src.get(s, ()))) & prim
            E = names if names and names != sg else sg  # rule `name_before_synonym` (make_queries.name_before_synonym)
            B = (set(src.get(s, ())) | set(syn[s]) | set(disp.get(s, ()))) & nb.gset(R)
            cmpR = tr.cmp_for(R)
            vs = (lambda x, r=R, c=cmpR: tr.vset(x, r) if c == "versioned" else sorted(x))
            rec = dict(qid=f"b5|syn|{s}|{R}|{R}", family=FAMILY, kind="b5_synonym", stratum="b5_syn_at115" if R == T else "b5_syn_atR",
                       unit=s, id=s, from_release=R, to_release=R, final_database=None, strategy="both", rule=GB.rule_of(len(E)),
                       cmp=cmpR, expect=vs(E), broad=vs(B) if B != E else None, locus=None,
                       name_src="biomart" if s in bm_syn else "core_synonyms")
            if len(E) > 1:
                sup = {tr.vid(g, R): nb.support(R, g, s, synonyms=not names) for g in E}  # read as a synonym: its sources
                if len(set(sup.values())) > 1:
                    rec["support"] = sup
            out.append(rec)
            if R < T:
                Es, Bs = nb.walks(E, R)
                Bs = set(Bs) | nb.walks(B, R)[1] | Es
                loc = tr.locus(next(iter(E)), R) if not Es and len(E) == 1 else None
                if not Es:
                    Bs |= cx.exact_T(s)[0] | cx.exact_T(s)[1]
                Et, At = GB.apart_split(nb, Es, T)
                ln = dict(qid=f"b5|syn|{s}|{R}|{T}", family=FAMILY, kind="b5_synonym_to_current", stratum="b5_syn_to115", unit=s,
                          id=s, from_release=R, to_release=T, final_database=None, strategy="both",
                          rule="set" if len(E) == 1 else GB.rule_of(len(Et)), cmp="versioned", expect=tr.vset(Et, T),
                          broad=tr.vset(Bs | At, T) if (Bs | At) != Et else None, locus=loc,
                          name_src="biomart" if s in bm_syn else "core_synonyms")
                if At:
                    ln["apart"] = tr.vset(At, T)
                if loc is not None:
                    ln = GB.M.apply_revisions(org, [ln], nb.C, names=GB.LOCUS_REVISIONS)[0]
                out.append(ln)
        nb._walk.clear()
        log(org, "B5 synonyms of BioMart", R, len(out) - n0)
    # the authorities' aliases, asked without a release
    auth = []
    h = TB.hgnc(org)
    if h:
        auth += [("hgnc", a, {r["ensembl"]} if r["ensembl"] else set()) for r in h for a in r["alias"]]
    m = TB.mgi(org)
    if m:
        auth += [("mgi", a, set(r["ensembl"])) for r in m["current"] if r["type"] == "Gene" for a in r["alias"]]
    v = TB.vgnc(org)
    if v:
        auth += [("vgnc", a, {r["ensembl"]} if r["ensembl"] else set()) for r in v for a in r["alias"]]
    n = TB.ncbi(org)
    if n:
        auth += [("ncbi", a, set(r["ensembl"])) for r in n for a in r["alias"]]
    per = collections.defaultdict(set)
    for who, a, ens in auth:
        if TB.clean(a):
            per[(who, a)] |= ens
    for (who, a), ens in sorted(per.items()):
        out.append(cx.default_line(f"b5|{who}|{a}|default|{T}", f"b5_alias_{who}", a, a, authority=ens, auth_db=who.upper()))
    log(org, "B5 authority aliases", len(per), "sources", sorted({w for w, _ in per}))
    log(org, "B5 lines without a release read where Ensembl's records hold them", cx.read_where_held(out))
    return publish(cx, "b5_v1", out, "B5 BioMart synonyms of every BioMart release (to R, to 115) and the authorities' "
                   "aliases without a release", dry)


# ---------------- B6 ----------------

FULLWIDTH = {c: chr(ord(c) + 0xFEE0) for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-"}


def b6(cx, dry):
    nb, tr, org, GT = cx.nb, cx.nb.tr, cx.org, cx.GT
    rng = random.Random(f"{SEED}|{org}|b6")
    base, dec = [], []
    mgi_lower = collections.defaultdict(set)
    if org == "mus_musculus":
        pp = nb.C.bm(T, "mgi", 1)
        for v, gs in (pp[1] if pp else {}).items():
            mgi_lower[v.lower()] |= gs & GT

    def add(v, E, apart, form, how, stratum, rule_now, dec_rule=None):
        names, syn = cx.exact_T(form)
        if names:
            return  # the form is itself a name of 115: its own question (B1, B5)
        Ev = tr.vset(E, T)
        rec = dict(qid=f"b6|{how}|{v}", family=FAMILY, kind="b6_form", stratum=stratum, unit=v, id=form, from_release=None,
                   to_release=T, final_database=None, strategy="both", cmp="versioned", broad=None, locus=None, note=f"{how} of {v}")
        if apart:
            rec["apart"] = tr.vset(apart, T)
        if syn - E:  # an exact synonym of another gene: rule `names_before_synonyms` reads it before a case or separator match
            rec.update(rule="set_exact" if len(syn) > 1 else "set", expect=tr.vset(syn, T), broad=tr.vset(syn | E, T),
                       note=rec["note"] + "; an exact synonym of another gene at 115 (rule `names_before_synonyms`)")
        else:
            other = (set(cx.lower_T.get(form.lower(), ())) & GT) - E
            Bv = cx.any_T(v) | E | apart  # the original string's genes in every source: what it may name
            if other - Bv and rule_now == "set":
                rec.update(rule="set", expect=[], open="another gene carries the form in another letter case",
                           open_within=tr.vset(E | other | Bv, T))
            elif rule_now == "set":
                rec.update(rule="set", expect=Ev, broad=tr.vset(Bv, T) if Bv != E else None)
            else:  # safe: nothing, or genes the original string names in some source
                rec.update(rule="safe", expect=tr.vset(Bv, T))
        base.append(rec)
        if dec_rule:
            r2 = dict(rec, revised_by=["b6_dec"], note=rec["note"] + f"; decided: {dec_rule[1]}")
            for k in ("open", "open_within", "broad"):
                r2.pop(k, None)
            Bd = cx.any_T(v) | set(dec_rule[0])
            r2.update(rule="set_exact" if len(dec_rule[0]) > 1 else "set", expect=tr.vset(dec_rule[0], T),
                      broad=tr.vset(Bd, T) if Bd != set(dec_rule[0]) else None)
            if (r2["rule"], r2["expect"]) != (rec["rule"], rec["expect"]) or rec.get("open"):
                dec.append(r2)

    pool = []
    for v in sorted(cx.disp_T):
        if not TB.clean(v):
            continue
        allg = set(cx.disp_T[v]) & GT
        prim = {g for g in allg if nb.region(g, T) != "patch"}
        E, apart = (prim, allg - prim) if prim else (allg, set())
        if not E:
            continue
        # exhaustive: letter case and white space
        for how, form in (("lower", v.lower()), ("upper", v.upper()), ("title", v[:1].upper() + v[1:].lower())):
            if form == v:
                continue
            d = None
            if org == "mus_musculus" and how in ("lower", "upper"):
                m = mgi_lower.get(v.lower(), set())
                if len(m) == 1:
                    d = (m, "rule `single_case_symbols`: the MGI symbol of exactly one gene without letter case")
            add(v, E, apart, form, how, "b6_case", "set", d)
        add(v, E, apart, f" {v} ", "space", "b6_space", "set")
        add(v, E, apart, f"\t{v}\r\n", "tab_crlf", "b6_space", "set")
        pool.append((v, E, apart))
    # sampled: separators, invisible characters, full-width letters
    forms = []
    for v, E, apart in pool:
        if re.search(r"[-_]", v):
            forms += [(v, E, apart, "sep_swap", v.translate(str.maketrans("-_", "_-")), "b6_sep", "set", None)]
            for how, rep in (("sep_dot", "."), ("sep_space", " "), ("sep_none", "")):
                forms.append((v, E, apart, how, re.sub(r"[-_]", rep, v), "b6_sep", "safe", None))
        d84 = (E, "rule `unicode_normalised`: NFKC, zero-width characters and byte-order marks removed")
        forms += [(v, E, apart, "nbsp", f" {v} ", "b6_invisible", "safe", d84),
                  (v, E, apart, "zwsp", f"{v[:len(v) // 2]}​{v[len(v) // 2:]}", "b6_invisible", "safe", d84),
                  (v, E, apart, "bom", f"﻿{v}", "b6_invisible", "safe", d84),
                  (v, E, apart, "fullwidth", "".join(FULLWIDTH.get(c, c) for c in v), "b6_fullwidth", "safe", d84)]
    forms = [f for f in forms if f[4] != f[0] and unicodedata.normalize("NFKC", f[4]).strip() != "" ]
    for v, E, apart, how, form, st, rule, d in sorted(rng.sample(forms, min(20000, len(forms))), key=lambda x: (x[0], x[3])):
        add(v, E, apart, form, how, st, rule, d)
    have = L.set_names(FAMILY, org) if cx.skip_published else []  # a rerun of the generator after an interruption
    m1 = None if "b6_v1" in have else publish(cx, "b6_v1", base, "B6 forms of every display name of 115: letter case and white space (all), separators, "
                 "invisible characters and full-width letters (20,000 sampled); the docs of the reference trees", dry)
    m2 = publish(cx, "b6_dec", dec, "B6 under decisions not built (84; mouse 55, 80)", dry, revises=["b6_v1"]) if dec else None
    return m1, m2


# ---------------- B7 ----------------

MONTHS = {"JAN": 1, "JANUARY": 1, "FEB": 2, "FEBRUARY": 2, "MAR": 3, "MARCH": 3, "APR": 4, "APRIL": 4, "MAY": 5, "JUN": 6,
          "JUNE": 6, "JUL": 7, "JULY": 7, "AUG": 8, "AUGUST": 8, "SEP": 9, "SEPT": 9, "SEPTEMBER": 9, "OCT": 10,
          "OCTOBER": 10, "NOV": 11, "NOVEMBER": 11, "DEC": 12, "DECEMBER": 12}
DATE_RE = re.compile(r"^(" + "|".join(sorted(MONTHS, key=len, reverse=True)) + r")[-_]?(\d{1,2})$", re.I)
ABBR = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
GREEK = [("TNFα", "TNF"), ("TNF-α", "TNF"), ("TNF-alpha", "TNF"), ("IL-1β", "IL1B"), ("IL1β", "IL1B"), ("IL-1beta", "IL1B"),
         ("IL-1α", "IL1A"), ("IFN-γ", "IFNG"), ("IFNγ", "IFNG"), ("IFN-gamma", "IFNG"), ("IFN-β", "IFNB1"), ("TGF-β1", "TGFB1"),
         ("TGFβ1", "TGFB1"), ("TGF-beta1", "TGFB1"), ("TGF-β2", "TGFB2"), ("TGF-β3", "TGFB3"), ("TGF-α", "TGFA"), ("HIF-1α", "HIF1A"),
         ("HIF1α", "HIF1A"), ("PGC-1α", "PPARGC1A"), ("PPARγ", "PPARG"), ("PPAR-γ", "PPARG"), ("PPARα", "PPARA"), ("PPARδ", "PPARD"),
         ("RXRα", "RXRA"), ("RARα", "RARA"), ("RARβ", "RARB"), ("ERα", "ESR1"), ("ERβ", "ESR2"), ("C/EBPα", "CEBPA"),
         ("C/EBPβ", "CEBPB"), ("GSK-3β", "GSK3B"), ("GSK3β", "GSK3B"), ("PKCα", "PRKCA"), ("β-catenin", "CTNNB1"),
         ("α-SMA", "ACTA2"), ("IL-1Ra", "IL1RN"), ("TNF-β", "LTA"), ("HNF-1α", "HNF1A"), ("HNF1β", "HNF1B"), ("HNF-4α", "HNF4A"),
         ("HNF4α", "HNF4A"), ("ROR-γt", "RORC"), ("RORγt", "RORC"), ("IκBα", "NFKBIA"), ("NF-κB1", "NFKB1"), ("Ngn3", "NEUROG3"),
         ("Pdx-1", "PDX1"), ("Nkx6.1", "NKX6-1"), ("Nkx2.2", "NKX2-2"), ("MafA", "MAFA")]


def excel_forms(month, day):
    d = datetime.date(2007, month, day) if day <= 28 or month != 2 else None
    out = [(f"{day}-{ABBR[month]}", "d-Mon"), (f"{day:02d}-{ABBR[month]}", "dd-Mon"), (f"{ABBR[month]}-{day:02d}", "Mon-dd")]
    if d:
        out += [(d.isoformat(), "iso"), (f"{month}/{day}/2007", "us"), (str((d - datetime.date(1899, 12, 30)).days), "serial")]
    return out


def b7(cx, dry):
    nb, tr, org, GT = cx.nb, cx.nb.tr, cx.org, cx.GT
    # every string of a date-like symbol family, with the genes (at 115) that carry it or carried it
    fam = collections.defaultdict(set)  # (month, day) -> genes at 115
    strings = collections.defaultdict(set)
    auth = []
    for r in (TB.hgnc(org) or []) + (TB.vgnc(org) or []):
        auth += [(s, {r["ensembl"]} if r["ensembl"] else set()) for s in [r["symbol"]] + r["prev"] + r["alias"]]
    m = TB.mgi(org)
    if m:
        auth += [(s, set(r["ensembl"])) for r in m["current"] for s in [r["symbol"]] + r["alias"]]
        auth += [(r["symbol"], set(r["ensembl"])) for r in m["withdrawn"]]
    for s, ens in auth:
        mm = DATE_RE.match(s)
        if mm and 1 <= int(mm.group(2)) <= 31:
            k = (MONTHS[mm.group(1).upper()], int(mm.group(2)))
            fam[k] |= ens & GT
            strings[k].add(s)
    for R in nb.rels:
        for v, gs in nb.display(R).items():
            mm = DATE_RE.match(v)
            if mm and 1 <= int(mm.group(2)) <= 31:
                k = (MONTHS[mm.group(1).upper()], int(mm.group(2)))
                strings[k].add(v)
                here = set(gs) & set(nb.genes(R))
                fam[k] |= (nb.walks(here, R)[0] if R < T else here) & GT
        nb._walk.clear()
    base, dec = [], []

    def add(qid, stratum, form, genes, unit, note, decided=True, reading=frozenset()):
        names, syn = cx.exact_T(form)
        if names or syn:
            return  # the form is a name or synonym of 115 in its own right
        rec = dict(qid=qid, family=FAMILY, kind="b7_form", stratum=stratum, unit=unit, id=form, from_release=None, to_release=T,
                   final_database=None, strategy="both", rule="safe", cmp="versioned", expect=tr.vset(set(genes) | reading, T),
                   broad=None, locus=None, note=note)
        base.append(rec)
        if decided:  # Q1 covers the spreadsheet forms and Greek letters, not the suffixes
            dec.append(dict(rec, revised_by=["b7_dec"], rule="set", expect=tr.vset(genes, T) if len(genes) == 1 else [],
                            note=note + "; decided (Q1): read where it resolves to exactly one gene, else nothing"))

    for (mo, day), genes in sorted(fam.items()):
        if not genes:
            continue
        for form, how in excel_forms(mo, day):
            add(f"b7|date|{how}|{mo}|{day}", "b7_date", form, genes, f"{mo}-{day}",
                f"Excel's {how} of {','.join(sorted(strings[(mo, day)])[:6])}")
    for form, sym in GREEK:
        g = set()
        for v, gs in cx.disp_T.items():
            if v.upper() == sym.upper():
                g |= set(gs) & GT
        if g:
            add(f"b7|greek|{form}", "b7_greek", form, g, sym, f"typography of {sym}",
                reading=set().union(*(cx.any_T(v) for v in cx.disp_T if v.upper() == sym.upper())))
    # suffixes: -AS1 on a gene without one; _at, _s_at (probe-set suffixes), on a seeded sample of current names
    rng = random.Random(f"{SEED}|{org}|b7")
    ever_lower = {k.lower() for k in nb.ever()}
    pool = sorted(v for v in cx.disp_T if TB.clean(v) and len(set(cx.disp_T[v]) & GT) == 1)
    for v in rng.sample(pool, min(2000, len(pool))):
        g = set(cx.disp_T[v]) & GT
        for how, form in (("as1", f"{v}-AS1"), ("at", f"{v}_at"), ("s_at", f"{v}_s_at")):
            if form.lower() in ever_lower:  # a string of some gene in some letter case: not a suffix on v
                continue
            add(f"b7|suffix|{how}|{v}", "b7_suffix", form, g, v, f"{how} suffix on {v}", decided=False, reading=cx.any_T(v))
    have = L.set_names(FAMILY, org) if cx.skip_published else []  # a rerun of the generator after an interruption
    m1 = None if "b7_v1" in have else publish(cx, "b7_v1", base, "B7 spreadsheet dates, Greek letters and suffixes: the symbol's genes or nothing (safe)", dry)
    m2 = publish(cx, "b7_dec", dec, "B7 under Q1 decided (exactly one gene, else nothing)", dry, revises=["b7_v1"])
    return m1, m2


def publish(cx, name, lines_, desc, dry, revises=None):
    if not lines_:
        log(cx.org, name, "no lines")
        return None
    man = L.publish(FAMILY, cx.org, name, lines_, desc, __file__, revises=revises, judge="b_v1",
                    sources=["gtf/genes_<R>.tsv.gz", "biomart/<R>/", "core_xref/<R>.tsv.gz", "hgnc/", "mgi/", "vgnc/", "ncbi/"],
                    dry=dry)
    log(cx.org, name, "dry" if dry else "published", man["n"], man["counts"])
    return man


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--sets", default="b4,b5,b6,b7")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    cx = Ctx(a.organism)
    cx.skip_published = not a.dry
    for s in a.sets.split(","):
        have = L.set_names(FAMILY, a.organism)
        if not a.dry and f"{s}_v1" in have and (s not in ("b6", "b7") or f"{s}_dec" in have):
            log(a.organism, s, "published before: skipped (a rerun after an interruption)")
            continue
        t0 = time.time()
        {"b4": b4, "b5": b5, "b6": b6, "b7": b7}[s](cx, a.dry)
        log(a.organism, s, f"{time.time() - t0:.0f}s")
    _ = C


if __name__ == "__main__":
    main()
