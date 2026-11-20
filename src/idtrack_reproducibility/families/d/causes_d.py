"""Family D's causes (a plugin of F0's causes.py: `--plugin families/d/causes_d.py`).

D1 judges by place (gen_d1.py); Ensembl's history is the other reading of Ensembl's own data. A wrong answer of D1
that follows Ensembl's history where the place says otherwise is idtrack following Ensembl's history (class b: a
decision the truth contradicts), never a fault of the reading; a wrong answer that leaves a successor of the history
out is idtrack's (class a). The tests read the D1 line (its `history`, `open`, `open_within`, `genes_at_source`) by
qid, and run before the generic tests for D1's lines only.

- `d1_history_to_another_place` (b): the answer is Ensembl's strict walk R -> 115, which names genes other than the
  place's (rule `place_moves_marked`, finding hh-2: Ensembl's history carries the id to another place, idtrack follows it).
- `d1_history_split_place_one` (b): the answer is the strict walk, several genes of which the place holds one
  (`all` returns every successor of a split, the place shows one gene; or `best` ties between them).
- `d1_history_part_not_place` (b): the answer is a part of the strict walk's genes, not (only) the place's (`best`
  chose another successor of the split).
- `d1_no_successor_nothing` (b, rule `bridge_same_place`): Ensembl's history retires the id; idtrack answers nothing; the place holds
  a gene of 115 (Ensembl gave the place a new id without an event).
- `d1_successor_dropped` (a): Ensembl's history has a successor, idtrack answers nothing.
- `d1_open_follows_history` (b): an open line (the place names no one gene) answered by Ensembl's history (strict or
  broad walk) to a gene outside the place.
- `d1_name_history_to_another_place` (b): a symbol line whose answer is the strict walk of the gene(s) the GTF of R
  names so, which leads elsewhere than the place.
- `d1_name_not_a_name_at_source` (e): a symbol line answered with a gene whose name the symbol is not at R by Ensembl's
  GTF, while the GTF names another gene so: the generic name tests (other name sources, synonyms) are tried first;
  this one only names what they leave (the answer's reading, from the explained path: `match`, the release read).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import causes as K
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

_LINES = {}


def d1_line(c):
    """The latest line of the row's qid over every published D1 set (or $D1_SETS)."""
    if not str(c.r.get("qid", "")).startswith("d1|"):
        return None
    org = c.w.org
    if org not in _LINES:
        sets = os.environ.get("D1_SETS", "").split(",") if os.environ.get("D1_SETS") else \
            [s for s in L.set_names("d", org) if s.startswith("d1")]
        lines, _, _ = L.load("d", org, sets)
        _LINES[org] = {q["qid"]: q for q in lines if q["qid"].startswith("d1|")}
    return _LINES[org].get(c.r["qid"])


def _sets(q):
    return set(q.get("expect") or []), set(q.get("history") or []), set(q.get("history_broad") or [])


def d1_history_to_another_place(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_id" or q.get("open") or not q.get("history"):
        return None
    E, Hs, _ = _sets(q)
    if c.gotb == Hs and not (Hs & E):
        return f"{q['id']} from {q['source_release']}: history -> {','.join(sorted(Hs))[:80]}; place -> {','.join(sorted(E))}"
    return None


def d1_history_split_place_one(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_id" or q.get("open") or not q.get("history"):
        return None
    E, Hs, _ = _sets(q)
    if c.gotb == Hs and E < Hs:
        return f"{q['id']} from {q['source_release']}: history -> {len(Hs)} genes; place -> {','.join(sorted(E))}"
    return None


def d1_history_part_not_place(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_id" or q.get("open") or not q.get("history"):
        return None
    E, Hs, _ = _sets(q)
    if c.gotb and c.gotb < Hs:
        return f"{q['id']}: answered {','.join(sorted(c.gotb))[:80]} of the history's {len(Hs)}; place -> {','.join(sorted(E))}"
    return None


def d1_no_successor_nothing(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_id" or q.get("open") or q.get("history") is None:
        return None
    if not q["history"] and not c.got:
        return f"{q['id']} retired by Ensembl after {q['source_release']}; place -> {','.join(q['expect'])}"
    return None


def d1_successor_dropped(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_id" or not q.get("history"):
        return None
    if not c.got:
        return f"{q['id']}: history -> {','.join(q['history'])[:80]}; nothing answered"
    return None


def d1_open_follows_history(c):
    q = d1_line(c)
    if not q or not q.get("open") or q["kind"] != "d1_id":
        return None
    _, Hs, Hb = _sets(q)
    if c.gotb and c.gotb <= (Hs | Hb):
        return f"{q['id']} ({q['truth_class']}): answered {','.join(sorted(c.gotb))[:80]}, Ensembl's walk, outside the place"
    return None


def d1_name_history_to_another_place(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or not q.get("genes_at_source") or not c.gotb:
        return None
    tr, R = c.w.tr, q["source_release"]
    walk = set()
    for g in q["genes_at_source"]:
        walk |= tr.walk_fwd(g, R, 115, True, living=True)
    walk = {C.base_id(x) for x in walk}
    if c.gotb == walk and not (walk & set(q.get("expect") or [])):
        return f"{q['id']} at {R}: GTF gene(s) {','.join(q['genes_at_source'])[:60]} -> {','.join(sorted(walk))[:60]}; place elsewhere"
    return None


def d1_name_not_a_name_at_source(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or not c.gotb:
        return None
    if c.gotb & set(q.get("genes_at_source") or []):
        return None
    m = c.res.get("match") if c.res else None
    return (f"{q['id']} at {q['source_release']}: GTF names {','.join(q.get('genes_at_source') or [])[:60] or 'no gene'}; "
            f"answered {','.join(sorted(c.gotb))[:60]} (match {m}, read at {c.read_at})")


def _name_walk(c, q):
    tr, R = c.w.tr, q["source_release"]
    walk = set()
    for g in q.get("genes_at_source") or []:
        walk |= tr.walk_fwd(g, R, 115, True, living=True)
    return {C.base_id(x) for x in walk}


def d1_name_in_other_source(c):
    """Every gene the answer adds was reached through a name edge whose database is not Ensembl's display name of that
    gene at R (the GTF does not name it so): a name in another source, or a synonym (rule `names_before_synonyms`)."""
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or not c.res:
        return None
    named = set(q.get("genes_at_source") or [])
    extra = c.gotb - set(q.get("expect") or []) - set(q.get("open_within") or []) - named - _name_walk(c, q)
    if not extra:
        return None
    parts = []
    for g in sorted(extra):
        ps = c.paths_to(g)
        if not ps:
            return None
        e0, d0 = ps[0]["edges"][0], ps[0]["edge_dbs"][0]
        if not d0:
            return None
        parts.append(f"{g} via {','.join(sorted({x.split('@')[0] for x in d0}))[:60]} at {e0[3] if len(e0) > 3 else '?'}")
    return f"{q['id']} at {q['source_release']}: " + "; ".join(parts)[:300]


def d1_name_several_part(c):
    """A symbol that several genes of R carry as their Ensembl name: `all` returns the successors of only some of them
    (the graph does not give the string to every gene so at R: miRBase names, rule `mirbase_genes`)."""
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or len(q.get("genes_at_source") or []) < 2 or not c.gotb:
        return None
    E = set(q.get("expect") or [])
    if E and c.gotb < E:
        return f"{q['id']} at {q['source_release']}: {len(q['genes_at_source'])} genes named so; answered {len(c.gotb)} of {len(E)} place genes"
    return None


def d1_name_no_successor_nothing(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or c.got or q.get("open") or not q.get("genes_at_source"):
        return None
    if not _name_walk(c, q):
        return f"{q['id']} at {q['source_release']}: Ensembl retires {','.join(q['genes_at_source'])[:60]}; place -> {','.join(q['expect'])[:60]}"
    return None


def d1_name_not_found(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or c.got or q.get("open"):
        return None
    if c.res and c.res.get("match") is None:
        return f"{q['id']} at {q['source_release']}: not found (no match); GTF names {','.join(q.get('genes_at_source') or [])[:60]}"
    return None


def d1_name_successor_dropped(c):
    q = d1_line(c)
    if not q or q["kind"] != "d1_symbol" or c.got or q.get("open"):
        return None
    w = _name_walk(c, q)
    if w:
        return f"{q['id']} at {q['source_release']}: read ({c.res.get('match') if c.res else '?'}), walk -> {','.join(sorted(w))[:60]}; nothing"
    return None


_D3 = {}


def d3_line(c):
    """The latest D3 line of the row's qid."""
    if not str(c.r.get("qid", "")).startswith("d3|"):
        return None
    org = c.w.org
    if org not in _D3:
        sets = [s for s in L.set_names("d", org) if s.startswith("d3")]
        lines, _, _ = L.load("d", org, sets)
        _D3[org] = {q["qid"]: q for q in lines if q["qid"].startswith("d3|")}
    return _D3[org].get(c.r["qid"])


def d3_read_elsewhere(c):
    """D3 `default`: idtrack read the symbol at another release than Ensembl's names say (the truth reads it at 115
    where it is a display name there, else at the latest release where it is one; idtrack reads any name database)."""
    q = d3_line(c)
    if not q or not c.res or q.get("from_release") is not None:
        return None
    ra = c.res.get("from_release")
    if ra is not None and q.get("read_at") is not None and int(ra) != int(q["read_at"]):
        return f"{q['id']}: read at {ra} (match {c.res.get('match')}), Ensembl's names at {q['read_at']}"
    return None


def d3_name_in_other_source(c):
    """D3: the genes the answer adds carry the string in another name source or as a synonym at the release read, not as
    Ensembl's display name (rule `names_before_synonyms`; the truth reads display names only)."""
    q = d3_line(c)
    if not q or not c.res:
        return None
    extra = c.gotb - set(q.get("expect") or []) - set(q.get("open_within") or []) - set(q.get("broad") or [])
    if not extra:
        return None
    parts = []
    for g in sorted(extra):
        ps = c.paths_to(g)
        if not ps or not ps[0]["edge_dbs"][0]:
            return None
        e0, d0 = ps[0]["edges"][0], ps[0]["edge_dbs"][0]
        parts.append(f"{g} via {','.join(sorted({x.split('@')[0] for x in d0}))[:60]} at {e0[3] if len(e0) > 3 else '?'}")
    return f"{q['id']} at {q.get('read_at')}: " + "; ".join(parts)[:300]


def d3_nothing(c):
    q = d3_line(c)
    if not q or c.got or not q.get("expect"):
        return None
    m = c.res.get("match") if c.res else None
    return f"{q['id']} at {q.get('read_at')}: nothing answered (match {m}); Ensembl -> {','.join(q['expect'])[:60]}"


def d3_part_of_several(c):
    q = d3_line(c)
    if not q or not c.gotb or len(q.get("expect") or []) < 2:
        return None
    if c.gotb < set(q["expect"]):
        return f"{q['id']} at {q.get('read_at')}: {len(c.gotb)} of the {len(q['expect'])} genes Ensembl names so"
    return None


MISSING = {"nan", "na", "none", "null", "n/a", ""}


def d4_missing_value_answered(c):
    """D4: a missing-value string in a list (`nan`, `NA`, `None`, `null`) answered (nothing, or a gene) where 71690ab's
    specification, rule `missing_values`, raises MissingValueError for that item."""
    if not str(c.r.get("qid", "")).startswith("d4|"):
        return None
    if str(c.id).strip().lower() in MISSING:
        return f"{c.id!r} in a list: answered {','.join(sorted(c.gotb)) or 'nothing'}; rule `missing_values` raises MissingValueError for it"
    return None


MIRBASE_FORM = __import__("re").compile(r"^[a-z]{3}[-_.]?(mir|let|lin)[-_.]", __import__("re").I)


def d4_mirbase_form(c):
    """D4: a miRBase name, in any of its separator forms (the harness's separator variants), not found: miRBase names
    are not in the graph (rule `mirbase_genes`), as the item's single call."""
    if not str(c.r.get("qid", "")).startswith("d4|") or c.got:
        return None
    if MIRBASE_FORM.match(str(c.id)):
        return f"{c.id}: a miRBase name (form of the harness's line); not found"
    return None


def d2_names_read_at_last_release(c):
    """D2: a list of symbols (or of symbols and ids) whose inferred releases are not the releases its Ensembl display
    names fit best: infer_source resolves every string at the graph's last release (its docstring), so a name that
    changed is not counted and a name new since counts at the newest releases only, and it counts a name in every name
    database, so neighbouring releases where the names live in another source tie (rule `inference_per_release`, decided,
    reads them at each candidate release)."""
    if not str(c.r.get("qid", "")).startswith("d2") or c.kind not in ("d2_symbols", "d2_mixed"):
        return None
    return f"inferred {','.join(sorted(c.got))[:40]}; Ensembl's names fit {','.join(sorted(c.exp))[:40]}"


def d4_separator_synonym(c):
    """D4: a case or separator variant of a name (the harness's `case_variant`, `separator_variant` lines) that idtrack reads
    as a synonym of other genes, at the latest release where that spelling is one (match `synonym:separator`)."""
    if not str(c.r.get("qid", "")).startswith("d4|") or not c.res:
        return None
    if str(c.res.get("match") or "").startswith(("synonym:separator", "synonym:case")):
        return f"{c.id}: read as a synonym at {c.read_at} by the separator search; got {','.join(sorted(c.gotb))[:80]}"
    return None


import __main__  # noqa: E402

for mod in {id(m): m for m in (K, __main__) if hasattr(m, "register") and hasattr(m, "TESTS")}.values():
    first = mod.TESTS[0].__name__[2:]
    for name, cls, dec, text, fn in (
            ("d1_history_to_another_place", "b", "rule `place_moves_marked`; finding hh-2",
             "D1: the answer is Ensembl's history walk, which leads to genes other than the place's", d1_history_to_another_place),
            ("d1_history_split_place_one", "b", "Ensembl's split; the place shows one gene",
             "D1: the answer is every successor of a split whose place holds one of them", d1_history_split_place_one),
            ("d1_history_part_not_place", "b", "ranking among Ensembl's successors",
             "D1: `best` answers a successor of the split that is not the place's gene", d1_history_part_not_place),
            ("d1_no_successor_nothing", "b", "rule `bridge_same_place`",
             "D1: Ensembl retires the id without successor, idtrack answers nothing; the place holds a 115 gene", d1_no_successor_nothing),
            ("d1_successor_dropped", "a", "-", "D1: Ensembl's history has a successor; idtrack answers nothing", d1_successor_dropped),
            ("d1_open_follows_history", "b", "rule `place_moves_marked`",
             "D1: the place names no one gene; the answer is Ensembl's walk to a gene outside the place", d1_open_follows_history),
            ("d1_name_history_to_another_place", "b", "rule `place_moves_marked`",
             "D1: a symbol read as the GTF's gene of R, whose history leads elsewhere than the place", d1_name_history_to_another_place)):
        mod.register(name, cls, dec, text, fn, before=first)
    for name, cls, dec, text, fn in (
            ("d1_name_no_successor_nothing", "b", "rule `bridge_same_place`",
             "D1: a symbol read as the GTF's gene of R, which Ensembl retires without successor; nothing answered, the place "
             "holds a 115 gene", d1_name_no_successor_nothing),
            ("d1_name_not_found", "c", "-", "D1: a symbol that R's GTF names is not found by idtrack at R", d1_name_not_found),
            ("d1_name_several_part", "c", "rule `mirbase_genes`",
             "D1: a name of several genes at R: `all` returns the successors of only some of them", d1_name_several_part),
            ("d1_name_successor_dropped", "a", "-", "D1: a symbol read at R whose gene's history has a successor; nothing answered",
             d1_name_successor_dropped),
            ("d1_name_in_other_source", "e", "rule `names_before_synonyms`",
             "D1: the answer adds genes that carry the string in another name source or as a synonym, not as Ensembl's "
             "display name at R (the truth reads only the display name)", d1_name_in_other_source)):
        mod.register(name, cls, dec, text, fn, before="external_bridge")
    mod.register("d1_name_not_a_name_at_source", "e", "rules `names_before_synonyms` and `name_before_synonym`",
                 "D1: the symbol answered with a gene it does not name in R's GTF (read in another name source or as a "
                 "synonym, or at another release)", d1_name_not_a_name_at_source, before="external_bridge")
    for name, cls, dec, text, fn in (
            ("d3_read_elsewhere", "e", "rule `start_release`",
             "D3 default: idtrack reads the symbol at another release than Ensembl's display names place it", d3_read_elsewhere),
            ("d3_name_in_other_source", "e", "rule `names_before_synonyms`",
             "D3: the answer adds genes that carry the string in another name source or as a synonym", d3_name_in_other_source),
            ("d3_part_of_several", "c", "rule `mirbase_genes`",
             "D3: a name of several genes; `all` returns only some of them", d3_part_of_several),
            ("d3_nothing", "c", "-", "D3: Ensembl names a gene so at the release read; idtrack answers nothing", d3_nothing)):
        mod.register(name, cls, dec, text, fn, before="external_bridge")
    mod.register("d2_names_read_at_last_release", "a", "rule `inference_per_release`",
                 "D2: a symbol list's releases inferred from names read at the graph's last release, in every name database", d2_names_read_at_last_release,
                 before=first)
    mod.register("d4_separator_synonym", "e", "rule `names_before_synonyms`; the separator search",
                 "D4: a separator variant read as another gene's synonym at an older release by the separator search",
                 d4_separator_synonym, before="external_bridge")
    mod.register("d4_mirbase_form", "c", "rule `mirbase_genes`", "D4: a miRBase name in a separator form, not found (miRBase not in the graph)",
                 d4_mirbase_form, before="external_bridge")
    mod.register("d4_missing_value_answered", "a", "rule `missing_values` (71690ab)",
                 "D4: a missing-value string in a list is answered instead of raising MissingValueError for that item",
                 d4_missing_value_answered, before=first)
    mod.register_kinds(history=("d1_id",), names=("d1_symbol", "d3_symbol"))
