"""Family C's causes (a plugin of F0's causes.py: `--plugin families/c/causes_c.py`).

Each test reads the scored row (query, expected, got, reason), the explain record where there is one, and Ensembl's
xref extracts of family C (`truth_c.read_xref`), never idtrack's code:

- `accession_version_never_had` (a; rule `accession_versions`, decided: the answer to the open choice Q2, docs/methods.md): a RefSeq or UniProt accession written with a
  version that Ensembl never gave it in the graph's window (below every version, in a gap, or above every version of
  an accession the graph no longer holds) is answered with the accession's gene: the version is dropped
  (`match: version`), where rule `version_newer_than_graph`'s rule answers nothing.
- `ncbi_label_name_by_release` (a; rule `one_kind_per_database`): the NCBI symbols of a gene asked under their one name `NCBI gene
  (formerly Entrezgene)` at a release where the graph files them under `EntrezGene` (to 87) or `NCBI gene` (88 to 99).
- `xref_on_other_form` (c; the externals configuration): a value Ensembl attaches to the gene at the release (on its
  gene, transcript or translation) that the graph does not hold because the configuration reads that database from
  another object form (the value is on an object of another form only).
- `db_not_held_at_release` (c): the database's values at the release are outside the configuration's releases for the
  release's assembly (rule `primary_assembly_only`: the graph holds the primary assembly only).
- `accession_old_version_read_at_target` (a; rule `accession_versions`, decided): an accession written
  with a version it had at an older release, not at the target, is read at the target without its version (the genes
  the accession names there), not at the release where it had that version.
- `value_of_another_database` (e; the judge reads one database): the string is also a value (accession or label) of
  another database of Ensembl, or the GTF's gene name, of the returned gene at the release read: without
  `from_database` the string names both genes, and the truth of one database cannot decide.
- `rfam_label_is_accession` (d; the truth): an RFAM label line whose label is written as an RFAM accession (Ensembl
  writes the accession as the label at some releases): the same string is the accession of the entry, which the
  label's line does not expect.
- `no_successor_answer_elsewhere` (b; rule `bridge_same_place`): a line without an explained path whose truth is a history walk with
  no successor and whose answer lies elsewhere (the judge's reason): the class of `external_bridge`, not confirmed
  by a path because only a sample of each unit is explained.
"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import causes as K
from idtrack_reproducibility.families.c import truth_c as TC

VER_NOTHING = {"c2_ver_lower", "c2_ver_gap", "c2_ver_higher_retired", "c3_ver_lower", "c3_ver_gap",
               "c3_ver_higher_retired"}
_XREF = collections.OrderedDict()  # what the truth's readers read per release, least recently used first
XREF_CACHE = 6  # releases' tables held at once (each a few hundred MB for the largest organisms)


def _cached(key, build):
    """`build()` cached under `key`, keeping the `XREF_CACHE` most recently used entries."""
    if key in _XREF:
        _XREF.move_to_end(key)
        return _XREF[key]
    while len(_XREF) >= XREF_CACHE:
        _XREF.popitem(last=False)
    _XREF[key] = build()
    return _XREF[key]


def accession_version_never_had(c):
    if c.kind not in VER_NOTHING or not c.got or c.exp:
        return None
    m = c.res.get("match") if c.res else None
    return f"{c.id} -> {','.join(sorted(c.got))[:80]} (match {m or 'version, by the answer'}): a version never had"


def ncbi_label_name_by_release(c):
    if c.db != "NCBI gene (formerly Entrezgene)" or not c.missing and not (set(c.exp) - set(c.got)):
        return None
    rel = c.to or 0
    if rel >= 100:
        return None
    name = "EntrezGene" if rel < 88 else "NCBI gene"
    return f"{c.id} at {rel}: the NCBI symbols are under `{name}` in the graph; asked under the name of 100 on"


def _xref(org, rel):
    """gene base -> {(db_name, db_display_name, acc, label, object type)} of a release (cached)."""
    def build():
        m = collections.defaultdict(set)
        for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, rel):
            m[gid].add((db, disp, acc, lab, ot.lower()))
        return m
    return _cached(("xref", org, rel), build)


SCOPE = {}


def _scope(org):
    if org not in SCOPE:
        from idtrack_reproducibility.harness import harness_lib as H
        SCOPE[org] = TC.Scope(org, CFG.get().tree("d4c2dda")["source"], H.gtf_assemblies(org))
    return SCOPE[org]


def xref_on_other_form(c):
    if not c.db or not str(c.r.get("stratum", "")).startswith("c6:"):
        return None
    missing = set(c.exp) - set(c.got)
    if not missing or c.to is None:
        return None
    org = c.w.org
    sc = _scope(org)
    form = sc.form.get(c.db)
    rows = _xref(org, c.to).get(TC.base(c.id), set())
    forms = collections.defaultdict(set)
    for db, disp, acc, lab, ot in rows:
        for v in (acc, lab):
            if v in missing:
                forms[v].add(ot)
    other = [v for v in missing if forms.get(v) and form not in forms[v]]
    if other and len(other) == len(missing):
        return f"{c.id} at {c.to}: {','.join(sorted(other))[:80]} on {','.join(sorted(set().union(*(forms[v] for v in other))))} " \
               f"only; the configuration reads {c.db} from {form}"
    return None


def db_not_held_at_release(c):
    if not c.db or not str(c.r.get("stratum", "")).startswith("c6:") or c.got:
        return None
    sc = _scope(c.w.org)
    if c.to is not None and not sc.included(c.db, c.to) and c.db != "NCBI gene (formerly Entrezgene)":
        return f"{c.db} is not held at {c.to} by the configuration"
    return None


def _strings(org, rel):
    """value -> genes (base) of every database of the extract at a release, and the GTF's gene names (cached)."""
    def build():
        m = collections.defaultdict(set)
        for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, rel):
            for v in (acc, lab):
                if v:
                    m[v].add(gid)
        from idtrack_reproducibility.harness import harness_lib as H
        for g, r in (H.gtf_genes(org, rel) or {}).items():
            if r.get("gene_name"):
                m[r["gene_name"]].add(g)
        return m
    return _cached(("strings", org, rel), build)


def value_of_another_database(c):
    if not c.extra or c.db or C_ENS.match(str(c.id)):
        return None
    rels = {r for r in (c.fr, c.read_at, c.to) if r}
    carriers = set()
    for r in rels:
        carriers |= _strings(c.w.org, r).get(c.id, set())
    hit = sorted(g for g in c.extra if g in carriers)
    if hit and len(hit) == len(c.extra):
        return f"{c.id}: also a value or name of {','.join(hit)[:100]} at {','.join(str(r) for r in sorted(rels))}"
    return None


def accession_old_version_read_at_target(c):
    if c.kind not in ("c2_oldver", "c3_oldver") or c.to is None:
        return None
    acc = str(c.id).rsplit(".", 1)[0]
    at_t = _strings(c.w.org, c.to).get(acc, set())
    if c.gotb and c.gotb <= at_t and c.gotb != c.expb:
        return f"{c.id}: answered with the accession's genes at {c.to} ({','.join(sorted(c.gotb))[:80]}), not read where " \
               f"it had that version"
    return None


def rfam_label_is_accession(c):
    if not (str(c.kind).startswith("c4_rfam_label") or c.kind == "c7_rfam_label"):
        return None
    if re.match(r"^RF\d{5}$", str(c.id)):
        return f"{c.id}: the label line's string is an RFAM accession"
    return None


def no_successor_answer_elsewhere(c):
    if c.paths or not str(c.r.get("reason", "")).startswith("Ensembl has no successor; not at the query's place"):
        return None
    return f"{c.id}: no successor; answer {','.join(sorted(c.gotb))[:80]} elsewhere (no explained path in the sample)"


import re  # noqa: E402
from idtrack_reproducibility.families import common as C0
C_ENS = C0.ENS_ID
import __main__  # noqa: E402

K.register_kinds(databases={"c6"})
for mod in {id(m): m for m in (K, __main__) if hasattr(m, "register") and hasattr(m, "TESTS")}.values():
    mod.register_kinds(databases={"c6"})
    mod.register("accession_version_never_had", "a", "rule `accession_versions` (decided)", "an accession version Ensembl never gave it in the window is answered with the "
                 "accession's gene (the version is dropped)", accession_version_never_had, before="entrezgene_label_rows")
    mod.register("ncbi_label_name_by_release", "a", "rule `one_kind_per_database` (decided)", "the NCBI symbols carry their name of 100 on "
                 "only from 100; before, they are filed under EntrezGene (to 87) or NCBI gene (88-99)",
                 ncbi_label_name_by_release, before="entrezgene_label_rows")
    mod.register("xref_on_other_form", "c", "externals configuration (one object form per database)",
                 "Ensembl attaches the value to the gene through an object of another form than the one the "
                 "configuration reads the database from", xref_on_other_form, before="name_in_other_source")
    mod.register("rfam_label_is_accession", "d", "-", "the label line's string is the RFAM accession (Ensembl writes "
                 "the accession as the label at some releases): the line should expect the entry's genes",
                 rfam_label_is_accession, before="entrezgene_label_rows")
    mod.register("accession_old_version_read_at_target", "a", "rule `accession_versions` (decided)",
                 "an accession version had at an older release is read at the target without its version",
                 accession_old_version_read_at_target, before="entrezgene_label_rows")
    mod.register("value_of_another_database", "e", "-", "the string is also a value or name of the returned gene in "
                 "another database at the release read; without from_database it names both",
                 value_of_another_database, before="name_in_other_source")
    mod.register("no_successor_answer_elsewhere", "b", "rule `external_bridge`; rule `bridge_same_place` (decided)",
                 "Ensembl's history has no successor and the answer lies elsewhere; no explained path in the sample "
                 "(external_bridge's class)", no_successor_answer_elsewhere)
    mod.register("db_not_held_at_release", "c", "externals configuration; rule `primary_assembly_only`",
                 "the configuration does not hold the database at the release's assembly", db_not_held_at_release,
                 before="name_in_other_source")
