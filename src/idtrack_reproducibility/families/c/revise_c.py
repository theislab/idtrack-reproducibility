"""Revision sets of family C: corrections of the judge where its truth is wrong, published under the qids of the lines
they change (a published set is never rewritten; `lines.publish(revises=...)`), with the evidence in the docstring.

Usage: revise_c.py <organism> --correction <name> --set <new set name> [--dry-run] [--redo <published revision set>]

Corrections:
- `rfam_label_accession`: an RFAM label line (`c4_rfam_label_*`, `c7_rfam_label`) whose label is written as an RFAM
  accession (`RF00003`): Ensembl writes the entry's accession as its label at some releases (core xrefs: `RFAM` rows
  with `display_label` = `dbprimary_acc`), so the string is also the entry's accession and names every gene of the
  entry; the label line expected only the genes of the label rows: the entry's other genes go to `broad` (not
  refuted). Evidence: the core xref extracts (`core_xref/<release>.tsv.gz`).
- `rfam_label_rows`: RFAM label lines under the decided rules `one_kind_per_database` and `rfam_display_label`: the
  line asks the label rows (`from_database_rows: labels`) and expects only the genes whose display name the label is
  (the function's docstring).
- `par_shared_names`: rule `par_shared_names` on C's lines (the harness's revision of that name; the function's
  docstring).
- `other_database_broad`: a line without `from_database` whose string is also a value (accession or label) of another
  Ensembl database, or the GTF's gene name, of other genes at the release read, or, for a line asked without a
  release and without a version, at the target where the string is then read (rule `start_release`; an accession
  with a version is read where it had it, rule `accession_versions`; the core xref extracts and the GTF of that
  release): without `from_database` the string names those genes too, so an answer within them is not refuted;
  they go to `broad` (their walks to the target for a line read at an older release), and to `expect` of a `safe`
  line (the string names them); the copies the line expects apart are not others. The corrected judge's
  `names_broad` for family C.
- `display_names_at_release`, `no_mito_nuclear_step`, `mirbase_genes`: the decided rules of those names on C's lines
  (each function's docstring).
- `version_not_judged`, `versions_in_window`: rule `accession_versions` as built after `next`'s first round
  (row I140) on the version strata (`VERSION_STRATA`), which the other corrections leave out (each
  function's docstring).
"""
import argparse
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.c import truth_c as TC

FAMILY = "c"
ME = os.path.abspath(__file__)


def log(*a):
    import time
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rfam_label_accession(org, lines, ctx):
    """RFAM label lines whose label is written as an RFAM accession (`RF00003`): Ensembl writes an entry's accession as
    its label at some releases, so the string is also the entry's accession and names every gene of the entry at the
    release read and at the target, where the string is read when no release is given (core xrefs of database `RFAM`
    with that accession). Those genes (walked to the target from an older release) are not refuted by the truth: they
    go to `broad`; the label rows' genes stay `expect`."""
    from idtrack_reproducibility.families.a import gen_a as GA
    car = ctx["car"] if "car" in ctx else ctx.setdefault("car", Carriers(org))
    for rec in sorted(lines, key=lambda r: (r.get("from_release") or r.get("read_at") or r["to_release"])):
        k = rec["kind"]
        if "rfam_label_accession" in (rec.get("revised_by") or []) or rec.get("excluded") or \
                not (k.startswith("c4_rfam_label") or k == "c7_rfam_label") or not re.match(r"^RF\d{5}$", str(rec["id"])):
            continue
        G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))  # built once, on the first such line
        R = rec.get("from_release") or rec.get("read_at") or rec["to_release"]
        t = rec["to_release"]
        genes = {g for g in car.entry(R).get(str(rec["id"]), set()) if g in (G.tr.genes.get(R) or {})}
        at_t = {g for g in car.entry(t).get(str(rec["id"]), set()) if g in (G.tr.genes.get(t) or {})}
        ext = set(G.tr.vset(at_t, t))  # the entry's genes at the target, where the string is read without a release
        if R == t:
            ext |= set(G.tr.vset(genes, t))
        else:
            for g in sorted(genes):
                E, B = G.walks(g, R, [t])
                ext |= set(G.tr.vset(E.get(t, set()) | B.get(t, set()), t))
        ext -= set(rec.get("expect") or []) | set(rec.get("broad") or []) | set(rec.get("apart") or [])
        if not ext:
            continue
        new = dict(rec)
        new["broad"] = sorted(set(rec.get("broad") or []) | set(rec.get("expect") or []) | ext)
        new["revised_by"] = list(rec.get("revised_by") or []) + ["rfam_label_accession"]
        yield new


class Carriers:
    """value -> genes (base) of every database of the core xref extract at a release and the GTF's gene names."""

    def __init__(self, org):
        self.org, self.cache = org, {}

    def at(self, rel):
        if rel not in self.cache:
            if len(self.cache) > 3:
                self.cache.clear()
            from idtrack_reproducibility.harness import harness_lib as H
            m = collections.defaultdict(set)
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(self.org, rel):
                for v in (acc, lab):
                    if v:
                        m[v].add(gid)
            for g, r in (H.gtf_genes(self.org, rel) or {}).items():
                if r.get("gene_name"):
                    m[r["gene_name"]].add(g)
            self.cache[rel] = m
        return self.cache[rel]

    def entry(self, rel):
        """accession -> genes (base) of database `RFAM` at a release (core xref extract)."""
        key = ("rfam", rel)
        if key not in self.cache:
            m = collections.defaultdict(set)
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(self.org, rel):
                if db == "RFAM" and acc:
                    m[acc].add(gid)
            self.cache[key] = m
        return self.cache[key]


def other_database_broad(org, lines, ctx):
    """A string that other databases or the GTF give to other genes at the release read: those genes (walked to the
    target) in `broad` (in `expect` of a `safe` line)."""
    from idtrack_reproducibility.families.a import gen_a as GA
    G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))
    tr = G.tr
    car = ctx["car"] if "car" in ctx else ctx.setdefault("car", Carriers(org))
    # lines sorted by the release read, so that one release's carriers are read once
    lines = sorted(lines, key=lambda r: (r.get("from_release") or r.get("read_at") or r["to_release"]))
    for rec in lines:
        if rec.get("from_database") or rec.get("final_database") or rec.get("excluded") or rec.get("open") or \
                "other_database_broad" in (rec.get("revised_by") or []):
            continue
        q = str(rec["id"])
        R = rec.get("from_release") or rec.get("read_at") or rec["to_release"]
        t = rec["to_release"]
        exp_b = {C.base_id(x) for x in rec.get("expect") or []}
        gs = car.at(R).get(q, set())
        if not gs and "." in q:
            gs = car.at(R).get(q.rsplit(".", 1)[0], set())  # an accession with its version: its unversioned value
        others = {g for g in gs if g in (tr.genes.get(R) or {})}
        at_t = set()
        # without a release an unversioned string is read at the target where it is held (rule `start_release`); a
        # version of an accession is read where it had that version (rule `accession_versions`), not here
        if rec.get("from_release") is None and R != t and not re.search(r"\.\d+$", q):
            gt = car.at(t).get(q, set()) or (car.at(t).get(q.rsplit(".", 1)[0], set()) if "." in q else set())
            at_t = set(tr.vset({g for g in gt if g in (tr.genes.get(t) or {})}, t))
        if R == t:
            extra = others - exp_b
            ext = set(tr.vset(extra, t))
        else:
            ext = set(at_t)
            for g in sorted(others):
                E, B = G.walks(g, R, [t])
                ext |= set(tr.vset(E.get(t, set()) | B.get(t, set()), t))
        # the line's own genes and the copies it expects reported apart (rule `judge_copies_apart`: a copy on a patch
        # or haplotype that carries the value is the line's, answered apart, never inside the answer) are no others
        ext -= set(rec.get("expect") or []) | set(rec.get("apart") or [])
        if not ext:
            continue
        new = dict(rec)
        if rec.get("rule") == "safe":
            new["expect"] = sorted(set(rec.get("expect") or []) | ext)
        else:
            new["broad"] = sorted(set(rec.get("broad") or []) | set(rec.get("expect") or []) | ext)
        new["other_carriers"] = sorted(ext)[:50]
        new["revised_by"] = list(rec.get("revised_by") or []) + ["other_database_broad"]
        yield new


class LabelNames:
    """Per release: the GTF's gene names, and label -> genes of the gene-level cross-references of every database but
    `RFAM` (core xref extract): what tells which RFAM label is a gene's display name."""

    def __init__(self, org):
        self.org, self.cache = org, {}

    def at(self, rel):
        if rel not in self.cache:
            if len(self.cache) > 2:
                self.cache.clear()
            from idtrack_reproducibility.harness import harness_lib as H
            other = collections.defaultdict(set)
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(self.org, rel):
                if ot == "Gene" and db != "RFAM" and lab:
                    other[lab].add(gid)
            names = {g: r.get("gene_name") for g, r in (H.gtf_genes(self.org, rel) or {}).items()}
            self.cache[rel] = (names, other)
        return self.cache[rel]


def rfam_label_rows(org, lines, ctx):
    """Rules `one_kind_per_database` and `rfam_display_label` (decided, built from idtrack's next version on): an RFAM
    label line (`c4_rfam_label_*`, `c7_rfam_label`) asks the label rows of `RFAM`, which a graph that records a name
    of their own (`record.label_names`) holds under that name (`RFAM label`): the line says so (`from_database_rows:
    labels`; the adapter asks the graph's name, a graph without the record `RFAM`). A label written as the entry's
    accession (`RF00003`) stays with the accession rows (rule `one_kind_per_database`: a label that repeats its
    accession is no row of its own where Ensembl names the labels like the accessions) and is not touched. And a
    label names only the gene whose display name it is (Ensembl's `gene.display_xref_id`): of the genes whose RFAM
    row carries the label at the release read, a gene whose display name in the GTF of that release is another
    string is no gene of the label (it leaves `expect`; recorded in `label_not_display`); a gene whose display name
    is the label while a cross-reference of another database of the gene carries the same label (an HGNC symbol
    projected to a snoRNA) has a display cross-reference the extract cannot tell (`broad`, recorded in
    `label_display_unknown`); the truth is computed again by C's generator (`gen_c.GenC.at_release`, `walked`) on
    the remaining genes. Where no gene is left, the label is no row of that release: the line is `excluded` (C's
    generator makes no line of a value not held at the release read), or `open` where only genes of the unknown kind
    carry it. A label the graph never holds (`c4_rfam_label_not_held`) keeps its truth. Evidence: the core xref
    extracts (`core_xref/<release>.tsv.gz`) and the GTFs (`gtf/genes_<release>.tsv.gz`)."""
    from idtrack_reproducibility.families.c import gen_c as GC
    names_of = ctx["names"] if "names" in ctx else ctx.setdefault("names", LabelNames(org))
    for rec in sorted(lines, key=lambda r: (r.get("from_release") or r.get("read_at") or r["to_release"])):
        k = rec["kind"]
        if "rfam_label_rows" in (rec.get("revised_by") or []) or rec.get("excluded") or \
                not (k.startswith("c4_rfam_label") or k == "c7_rfam_label") or re.match(r"^RF\d{5}$", str(rec["id"])):
            continue
        gc = ctx["gc"] if "gc" in ctx else ctx.setdefault("gc", GC.GenC(org, 115, ["rfam_label"], authority=False))
        lab = str(rec["id"])
        R = rec.get("from_release") or rec.get("read_at") or rec["to_release"]
        t = rec["to_release"]
        new = dict(rec)
        new["from_database_rows"] = "labels"
        new["revised_by"] = list(rec.get("revised_by") or []) + ["rfam_label_rows"]
        h = gc.X.at("rfam_label", lab, R)
        genes = set(h[0]) if h else set()
        names, other = names_of.at(R)
        display = {g for g in genes if names.get(g) == lab and g not in other.get(lab, ())}
        unknown = {g for g in genes if names.get(g) == lab and g in other.get(lab, ())}
        not_display = genes - display - unknown
        if k.endswith("_not_held") or not genes or not (unknown or not_display):
            # a label the graph never holds (`_not_held`: nothing, or its gene, is right) or one whose genes all carry it
            # as their display name: the truth stands
            yield new
            continue
        tr = gc.tr
        new["label_not_display"] = sorted(not_display)[:50]
        if unknown:
            new["label_display_unknown"] = sorted(unknown)[:50]
        if not display and not unknown:
            # the label is no gene's display name at the release read: the graph holds no row of it there, and C's
            # generator makes no line of a value not held at that release (`gen_c.GenC.held_strata`)
            new["excluded"] = (f"rule `rfam_display_label`: {lab} is the display name of no gene at {R} (the genes whose "
                               f"RFAM row carries it are named otherwise), so no label row of it is held there")
            yield new
            continue
        if not display:
            new["open"] = (f"rule `rfam_display_label`: at {R} {lab} is the display name only of genes where another "
                           "database's cross-reference carries the same label; which is the display cross-reference, and so "
                           "whether a label row is held there, the extract cannot tell")
            for f in ("open_within",):
                new.pop(f, None)
            yield new
            continue
        tru = gc.at_release(display, R) if R == t else gc.walked(display, R, t)
        if R == t:
            unk = set(tr.vset({g for g in unknown if g in (tr.genes.get(R) or {})}, t))
        else:
            unk = set()
            for g in sorted(unknown):
                if g in (tr.genes.get(R) or {}):
                    E, B = gc.G.walks(g, R, [t])
                    unk |= set(tr.vset(E.get(t, set()) | B.get(t, set()), t))
        for f in ("expect", "broad", "apart", "open", "open_within", "locus", "locus_mapped", "locus_projected", "exceptions"):
            new.pop(f, None)
        new.update({f: x for f, x in tru.items() if x is not None})
        unk -= set(new.get("expect") or []) | set(new.get("apart") or [])
        if unk:
            new["broad"] = sorted(set(new.get("broad") or []) | set(new.get("expect") or []) | unk)
        yield new


def par_shared_names(org, lines, ctx):
    """Rule `par_shared_names` (decided, built from idtrack's next version on): a value that the X and Y copies of a
    pseudoautosomal gene share names the X gene, the Y copy reported apart: the accuracy harness's revision of the
    same name (`make_queries.rev_par_shared_names`, Ensembl's assembly exceptions and GTF places) applied to C's lines."""
    from idtrack_reproducibility.harness import make_queries as M
    from idtrack_reproducibility.harness import harness_lib as H
    C78 = ctx["C78"] if "C78" in ctx else ctx.setdefault("C78", M.Ctx78(org, H.ORGANISMS[org][2]))
    todo = [r for r in lines if "par_shared_names" not in (r.get("revised_by") or [])]
    for v in M.rev_par_shared_names(org, todo, C78).values():
        yield from v


class HeldEarly:
    """Rule `display_names_at_release` (decided) as Ensembl's records give it: at a release, the names Ensembl gives a
    gene as its display name (`core_display/<release>`: the display cross-reference's label; the name, not the
    accession of the record it is filed under, which is held only where a cross-reference of the gene carries it), and the
    NCBI gene ids of any of its forms (`core_xref`: `EntrezGene` on the gene, a transcript or a translation), read
    for the values asked only; {value: {release: genes}}."""

    def __init__(self, org, values):
        self.org, self.want = org, set()
        self.at = collections.defaultdict(lambda: collections.defaultdict(set))
        self.add(values)

    def add(self, values):
        """Read the records of more values (a correction reads the published sets one after another)."""
        from idtrack_reproducibility.harness import harness_lib as H
        org = self.org
        want = set(values) - self.want
        if not want:
            return
        self.want |= want
        for rel in TC.xref_releases(org):
            for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, rel):
                if db == "EntrezGene" and acc in want:
                    self.at[acc][rel].add(gid)
            p = os.path.join(H.truth_dir(org), "core_display", f"{rel}.tsv.gz")
            if os.path.exists(p):
                for r in H.read_tsv(p):
                    if r.get("label") in want:
                        self.at[r["label"]][rel].add(TC.base(r["gene_id"]))


def display_names_at_release(org, lines, ctx):
    """Rule `display_names_at_release` (decided, built from idtrack's next version on): every display name of every
    release is a name of its gene there, whatever database Ensembl filed it under, and NCBI gene ids are read on every
    form (the gene, its transcripts, its translations). A value the reference tree's configuration never holds (rule
    `safe`: `c1_not_held`, `c3_entry`, `c4_*_not_held`, ...) that Ensembl gives a gene in one of these ways at some
    release of the window (`HeldEarly`) is held there now: the line is read as a held value is (rule `start_release`,
    C's generator `gen_c.default_truth`): at the target where it is held there, else at the latest release where it
    is, its genes walked to the target (`gen_c.walked_truth`: strict walk with living ids, broad, copies apart, the
    locus rule), rule `set`; `read_at` says where. Values held nowhere so keep their truth."""
    from idtrack_reproducibility.families.a import gen_a as GA
    from idtrack_reproducibility.families.c import gen_c as GC
    todo = [r for r in lines if r.get("rule") == "safe" and not r.get("excluded") and not r.get("final_database")
            and "display_names_at_release" not in (r.get("revised_by") or [])]
    if not todo:
        return
    values = {str(r["id"]) for r in todo}
    held = ctx["held"] if "held" in ctx else ctx.setdefault("held", HeldEarly(org, values))
    held.add(values)
    G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))
    for rec in todo:
        t = rec["to_release"]
        at = {r: g for r, g in held.at.get(str(rec["id"]), {}).items() if r <= t}
        if not at:
            continue
        R = rec["from_release"] if rec["from_release"] in at else (t if t in at else max(at))
        tru = GC.walked_truth(G, at[R], R, t)
        new = dict(rec)
        for f in ("expect", "broad", "apart", "open", "open_within", "locus", "locus_mapped", "locus_projected",
                  "exceptions", "partial"):
            new.pop(f, None)
        new.update({f: x for f, x in tru.items() if x is not None})
        new.update(rule="set", read_at=R, revised_by=list(rec.get("revised_by") or []) + ["display_names_at_release"])
        new.pop("set", None)
        yield new


def no_mito_nuclear_step(org, lines, ctx):
    """Rule `no_mito_nuclear_step` (decided, built from idtrack's next version on): a history step between the
    mitochondrial genome and a nuclear sequence is not followed, and the versions of a stable id on its two sides are
    two genes (family A's truth, `truth.TruthA`, follows it). A line read at one release and asked at another whose
    value names, at the release read (`core_xref`, any database's accession or label), a gene of the lineages the rule
    changes (`Truth.mito_lineage`) is computed again (`gen_c.walked_truth`, as C's generator) with and without the
    rule; where they differ the rule's truth replaces the line's."""
    from idtrack_reproducibility.families.a import gen_a as GA
    from idtrack_reproducibility.families.c import gen_c as GC
    G = ctx["G"] if "G" in ctx else ctx.setdefault("G", GA.GenA(org, 115))
    lineage = G.tr.mito_lineage()
    if not lineage:
        return
    G0 = ctx["G0"] if "G0" in ctx else ctx.setdefault("G0", _gen_a_without_rule(org))
    carriers = ctx["mito_carriers"] if "mito_carriers" in ctx else ctx.setdefault("mito_carriers", _carriers_of(org, lineage))
    for rec in lines:
        if rec.get("excluded") or rec.get("final_database") or "no_mito_nuclear_step" in (rec.get("revised_by") or []):
            continue
        R = rec.get("from_release") or rec.get("read_at")
        t = rec["to_release"]
        if R is None or R == t:
            continue
        v = str(rec["id"])
        genes = carriers.get((R, v)) or carriers.get((R, re.sub(r"\.\d+$", "", v)))
        if not genes:
            continue
        a, b = GC.walked_truth(G, genes, R, t), GC.walked_truth(G0, genes, R, t)
        if a == b:
            continue
        new = dict(rec)
        for f in ("expect", "broad", "apart", "open", "open_within", "locus", "locus_mapped", "locus_projected",
                  "exceptions", "partial"):
            new.pop(f, None)
        new.update({f: x for f, x in a.items() if x is not None})
        new["revised_by"] = list(rec.get("revised_by") or []) + ["no_mito_nuclear_step"]
        new.pop("set", None)
        yield new


def _gen_a_without_rule(org):
    """Family A's generator with the history read as Ensembl recorded it (no rule `no_mito_nuclear_step`)."""
    from idtrack_reproducibility.families.a import gen_a as GA
    from idtrack_reproducibility.families.astar import truth as TA
    G0 = GA.GenA(org, 115)
    G0.tr = TA.TruthA(org, 115, mito_cut=False)
    return G0


def _carriers_of(org, genes):
    """{(release, value): genes} of the cross-references (accession and label, any database) of the given genes."""
    out = collections.defaultdict(set)
    for rel in TC.xref_releases(org):
        for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, rel):
            if gid in genes:
                for v in (acc, lab):
                    if v:
                        out[(rel, v)].add(gid)
    return out


def mirbase_genes(org, lines, ctx):
    """Rule `mirbase_genes` (decided, built from idtrack's next version on): a miRBase entry that Ensembl gives to genes
    at two loci is no identifier of them: the accuracy harness's revision of the same name
    (`make_queries.rev_mirbase_genes`: the loci from Ensembl's cross-references and the GTFs) applied to C's lines."""
    from idtrack_reproducibility.harness import make_queries as M
    from idtrack_reproducibility.harness import harness_lib as H
    C78 = ctx["C78"] if "C78" in ctx else ctx.setdefault("C78", M.Ctx78(org, H.ORGANISMS[org][2]))
    todo = [r for r in lines if "mirbase_genes" not in (r.get("revised_by") or [])]
    for v in M.rev_mirbase_genes(org, todo, C78).values():
        yield from v


# The strata of accessions written with a version (c2 RefSeq, c3 UniProt): the lines the version corrections read.
VERSION_STRATA = re.compile(r"^c[23]_(at_v|default_v|oldver|ver_[a-z_]+)$")
VERSION_FAMILY = {"c2": "refseq", "c3": "uniprot"}
TRUTH_FIELDS = ("expect", "broad", "apart", "open", "open_within", "locus", "locus_mapped", "locus_projected",
                "exceptions", "partial", "cmp")


def _version_ctx(org, ctx):
    """C's generator over the versioned databases (Ensembl's cross-references of every release, family A's walks) and
    the edit counters of rule `accession_versions` measured from the same records (`truth_c.version_change_shares`), built once."""
    from idtrack_reproducibility.families.c import gen_c as GC
    if "GV" not in ctx:
        ctx["GV"] = GC.GenC(org, 115, ["refseq", "uniprot"], authority=False, continuity=ctx.get("continuity", False))
        dbs = set(TC.FAMILIES["refseq"]["dbs"]) | set(TC.FAMILIES["uniprot"]["dbs"])
        ctx["shares"] = TC.version_change_shares(org, dbs, upto=115)
        ctx["counters"] = TC.edit_counter_dbs(ctx["shares"])
        log("version shares", {db: f"{c}/{n} = {c / n:.3f}" for db, (n, c) in sorted(ctx["shares"].items()) if n},
            "edit counters", sorted(ctx["counters"]))
    return ctx["GV"], ctx["counters"]


# The fields of a line that hold its truth: a versioned line read as another line of the sets takes these from it.
REF_FIELDS = TRUTH_FIELDS + ("rule", "other_carriers", "primaries", "spellings", "core_places")
OWN_SETS = ("version_not_judged", "version_above_latest", "versions_in_window")  # the version corrections' own sets


def _with_truth(rec, tru, name, **extra):
    """The line with a new truth (`tru`'s fields of REF_FIELDS), the correction named in `revised_by`; None when the
    truth is the line's own."""
    if all((tru.get(f) or None) == (rec.get(f) or None) for f in REF_FIELDS):
        return None
    new = dict(rec)
    for f in REF_FIELDS + ("q2", "read_as", "latest_version", "version", "read_as_name"):
        new.pop(f, None)
    new.update({f: x for f, x in tru.items() if f in REF_FIELDS and x is not None})
    new.setdefault("rule", "set")
    new["revised_by"] = list(rec.get("revised_by") or []) + [name]
    new.update(**extra)
    new.pop("set", None)
    return new


def _references(org, ctx, qids):
    """{qid: latest line} of the lines a versioned line is read as, from every published set but the version
    corrections' own (their truth carries every correction published for them: par_shared_names, other_database_broad,
    no_mito_nuclear_step, ...)."""
    want = set(qids) - set(ctx.setdefault("refs", {}))
    if want:
        sets = [n for n in L.resolve_sets(FAMILY, org, None, "revised") if not n.startswith(OWN_SETS)]
        ctx["refs"].update(L.load_qids(FAMILY, org, want, sets=sets))
    return ctx["refs"]


def _nothing(g, t):
    return dict(cmp=g.tr.cmp_for(t), expect=[], rule="set")


def _read_as(org, ctx, lines, name, reading):
    """The lines of the version strata asked without a release, each read as `reading(rec, acc, ver, fam, g)` says:
    (qid of the line it is read as, read_at, extra fields) or None for nothing (or to leave the line: `"leave"`). A
    line read as a line of the sets takes that line's latest truth (`_references`; a line read as itself, its truth
    before the version corrections' sets); where the sets hold no such line, the truth C's generator gives the reading
    (`gen_c.GenC.at_release`, `walked`). And rule `display_names_at_release` comes first where it reads the string
    later: a versioned string that Ensembl gives a gene as its display name (`HeldEarly`: the label of RefSeq's
    `RefSeq_dna_predicted` before 64, `XM_001083595.1`) is read as that name at the latest release where it is one,
    when that is later than the version's own reading (`read_as_name`), its genes walked. Only lines whose truth that changes are given: running a correction again
    (`<name>_v2`) publishes what its earlier sets made wrong, and puts back what they changed without cause."""
    todo = []
    for rec in lines:
        st = rec["stratum"]
        if not VERSION_STRATA.match(st) or rec.get("from_release") is not None or rec.get("excluded"):
            continue
        m = TC.VERSIONED_ACC.match(str(rec["id"]))
        if m:
            todo.append((rec, m.group(1), int(m.group(2))))
    if not todo:
        return
    g, counters = _version_ctx(org, ctx)
    held = ctx["held"] if "held" in ctx else ctx.setdefault("held", HeldEarly(org, ()))
    held.add({str(rec["id"]) for rec, _, _ in todo})
    plan = []
    for rec, acc, ver in todo:
        fam = VERSION_FAMILY[rec["stratum"][:2]]
        how = reading(rec, acc, ver, fam, g, counters)
        if how == "leave":
            continue
        # rule `display_names_at_release`: the whole string is a gene's display name at a release; read there where
        # that is later than the version's own reading (idtrack's steps read the string as a name before its version)
        t = rec["to_release"]
        at = {r: gs for r, gs in held.at.get(str(rec["id"]), {}).items() if r <= t}
        if at and max(at) > (how[1] if how else -1):
            Rd = max(at)
            how = (None, Rd, dict(read_as_name=True), at[Rd])
        plan.append((rec, acc, fam, how))
    refs = _references(org, ctx, {h[0] for _, _, _, h in plan if h and h[0]})
    missing = collections.Counter()
    for rec, acc, fam, how in plan:
        t = rec["to_release"]
        if how is None:
            tru, extra = _nothing(g, t), dict(read_at=None)
        else:
            q, R, extra = how[:3]
            ref = refs.get(q) if q else None  # a line read as itself: its truth before the version corrections' sets
            if ref is not None and not ref.get("excluded"):
                tru = ref
            else:
                missing[rec["stratum"] if q else "read as a display name"] += 1
                genes = how[3] if len(how) > 3 else g.X.at(fam, acc, R)[0]
                tru = dict(g.at_release(genes, t) if R == t else g.walked(genes, R, t), rule="set")
            extra = dict(extra, read_at=R) if q == rec["qid"] or not q else dict(extra, read_at=R, read_as=q)
        if how is not None and q == rec["qid"] and ref is not None:  # put back: the line as before the corrections' sets
            extra["revised_by"] = list(ref.get("revised_by") or []) + [name]
        new = _with_truth(rec, tru, name, **extra)
        if new:
            yield new
    if missing:
        log(name, "read as lines the sets do not hold (truth from the generator):", dict(missing))


def version_not_judged(org, lines, ctx):
    """Rule `accession_versions`, its edit counters (built in idtrack's next version after its first round): a version dates an accession only where it
    is a property of the sequence record that holds over time. Of the pairs of consecutive releases at which one link
    (the accession, the Ensembl object with its version, the database) holds, the share in which Ensembl's stored
    `xref.version` differs, measured from Ensembl's cross-references of every release (`truth_c.version_change_shares`);
    above one half (`truth_c.EDIT_COUNTER_SHARE`) the database's number counts its records' edits (UniProt's entry
    version) and a typed version is not judged: the accession is read without it, never nothing. A line of a version
    stratum (`VERSION_STRATA`) asked without a release, whose accession only databases of that kind carry, is read as
    the accession's own line without a version (`c3_default|<acc>|-|<t>` where the graph holds it at the target, else
    `c3_retired|...`; its latest truth, every correction of it included; `read_as` names it; `version: not judged`);
    an accession the graph never holds: nothing. The lines read at a release (`c3_at_v`) already expect the
    accession's reading there. On the five organisms' records Swiss-Prot and TrEMBL are such databases, RefSeq's are
    not (the shares are logged)."""

    def reading(rec, acc, ver, fam, g, counters):
        kinds = set(TC.FAMILIES[fam]["dbs"]) & counters
        if not kinds:
            return "leave"  # its numbers date sequences: judged (`versions_in_window`)
        if kinds != set(TC.FAMILIES[fam]["dbs"]):  # a family of both kinds would need each accession's own databases
            raise NotImplementedError(f"{fam}: some but not all databases are edit counters ({sorted(counters)})")
        t = rec["to_release"]
        R, _ = g.default_truth(fam, acc)
        if R is None:
            return None
        part = rec["stratum"][:2]
        return (f"{part}_default|{acc}|-|{t}" if R == t else f"{part}_retired|{acc}|-|{t}"), R, \
            dict(version="not judged (an edit counter)")

    yield from _read_as(org, ctx, lines, "version_not_judged", reading)


def _judged_reading(rec, acc, ver, fam, g, counters):
    """Rule `accession_versions` for a database whose number dates a sequence, on the versions the graph's window had
    (up to the line's target): a version it had, read at the latest release with it; above every version it had, the
    latest, where the accession is attached at the line's target (the graph's last release); any other version (above
    the versions of an accession no longer attached there, or one the window never held), nothing. The line
    read as: `<part>_default_v|<acc>.<v>|-|<t>` where the version is held at the target, else `<part>_oldver|...`."""
    if set(TC.FAMILIES[fam]["dbs"]) & counters:
        return "leave"  # not judged (`version_not_judged`)
    t = rec["to_release"]
    had = {int(x) for x in g.X.versions_in_graph(fam, acc, t) if x.isdigit()}
    h = g.X.at(fam, acc, t)
    held_t = bool(h and h[2])  # the accession is attached at the graph's last release
    read = ver if ver in had else (max(had) if had and ver > max(had) and held_t else None)
    R = g.X.last_with_version(fam, acc, str(read), t) if read is not None else None
    if R is None:
        return None
    part = rec["stratum"][:2]
    q = f"{part}_default_v|{acc}.{read}|-|{t}" if R == t else f"{part}_oldver|{acc}.{read}|-|{t}"
    return q, R, ({"latest_version": read} if read != ver else {})


def versions_in_window(org, lines, ctx):
    """Rule `accession_versions` as built in idtrack's next version after its first round (row I140: an accession's
    versions are spellings of its links), on the versions the graph's window had. A version it had is read at the
    latest release with it; a version above every version the accession had is read as the latest where the accession
    is attached at the graph's last release (row I150: elsewhere it is an unknown version, nothing, as rule
    `version_newer_than_graph` reads an Ensembl id); any other version, nothing. And the window ends at the line's
    target: the generator read the versions over every release of the extracts, which reach one release beyond the
    graphs' window (116, where the reference tree's configuration lists it), so a version written only at 116 made a
    version below it a gap (`c2_ver_gap`: nothing). Every line of a version stratum asked without a release, of a
    database whose number dates a sequence (not the edit counters of `version_not_judged`), is read so
    (`Xrefs.versions_in_graph(..., upto)`; `_judged_reading`), as the line of the version read (its latest truth,
    every correction of it included; `read_as` names it), or expects nothing. The sets `version_above_latest_v1`
    (2026-11-21) and the first `versions_in_window` sets read the latest also for accessions no longer attached at the
    target; this correction run again puts those lines back to nothing."""
    yield from _read_as(org, ctx, lines, "versions_in_window", _judged_reading)


# The strata of held_strata (gen_c) whose reading rule `start_release` gives, by the prefix of their kind: the family read.
HELD_PREFIX = {"c1": "ncbi", "c2": "refseq", "c3": "uniprot", "c4_rfam_acc": "rfam", "c4_rfam_label": "rfam_label",
               "c5_ccds": "ccds", "c5_ottg": "ottg", "c5_ottt": "ottt", "c5_ottp": "ottp", "c5_vega_gene": "vega_gene",
               "c5_clone": "clone"}
HELD_KIND = re.compile(r"^(c1|c2|c3|c4_rfam_acc|c4_rfam_label|c5_ccds|c5_ottg|c5_ottt|c5_ottp|c5_vega_gene|c5_clone)_"
                       r"(default|retired|not_held)$")


def database_continuity(org, lines, ctx):
    """Rule `database_filing` (rows I149 and I157, built in idtrack after its second round): a database of
    Ensembl that continues one the graph includes (the links that move between them from one release to the next:
    `truth_c.continuations`, read from Ensembl's cross-references) is that database, so its values are held at its
    releases: RefSeq's accessions that Ensembl filed as `RefSeq_dna` and `RefSeq_dna_predicted` before release 64 are
    RefSeq_mRNA's, RefSeq_ncRNA's and their predicted databases' there, Rfam's accessions on transcripts before 63
    Rfam's. Two kinds of line read again with that reading (`gen_c.GenC(..., continuity=True)`): the versioned lines
    asked without a release (rule `accession_versions` as `versions_in_window` reads it, the versions those databases
    wrote counted as the accession's versions: a version `c2_ver_lower` said the accession never had, which
    `RefSeq_dna_predicted` wrote, is read at the latest release with it); and the lines of a value read by rule
    `start_release` (`<part>_default`, `_retired`, `_not_held` of rule `safe`): at the latest release where the graph
    holds it now, its genes walked to the target (`gen_c.GenC.default_truth`). Only lines whose truth changes are
    given."""
    from idtrack_reproducibility.families.c import gen_c as GC
    ctx["continuity"] = True
    lines = list(lines)
    yield from _read_as(org, ctx, [r for r in lines if VERSION_STRATA.match(r["stratum"])], "database_continuity",
                        _judged_reading)
    todo = []
    for rec in lines:
        m = HELD_KIND.match(rec.get("kind") or "")
        if not m or rec.get("excluded") or rec.get("final_database") or rec.get("from_release") is not None:
            continue
        if m.group(2) == "not_held" and rec.get("rule") != "safe":
            continue  # a not-held value another correction already reads where held (display_names_at_release)
        if rec["to_release"] != 115:
            continue
        todo.append((rec, HELD_PREFIX[m.group(1)]))
    if not todo:
        return
    if "GK" not in ctx:
        ctx["GK"] = GC.GenC(org, 115, sorted(set(HELD_PREFIX.values())), authority=False, continuity=True)
    g = ctx["GK"]
    for rec, fam in todo:
        Rd, tru = g.default_truth(fam, str(rec["id"]))
        if Rd is None:
            continue
        if Rd == rec.get("read_at") and rec.get("rule") != "safe":
            continue  # read where it was: its truth (and every correction of it) stands
        new = _with_truth(rec, dict(tru, rule="set"), "database_continuity", read_at=Rd)
        if new:
            yield new


VERSION_CORRECTIONS = {"version_not_judged", "versions_in_window"}

CORRECTIONS = {"rfam_label_accession": rfam_label_accession, "other_database_broad": other_database_broad,
               "rfam_label_rows": rfam_label_rows, "par_shared_names": par_shared_names,
               "display_names_at_release": display_names_at_release, "no_mito_nuclear_step": no_mito_nuclear_step,
               "mirbase_genes": mirbase_genes, "version_not_judged": version_not_judged,
               "versions_in_window": versions_in_window, "database_continuity": database_continuity}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--correction", required=True, choices=sorted(CORRECTIONS))
    ap.add_argument("--set", required=True)
    ap.add_argument("--sets", default=None, help="the sets read (default: every published set)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--lines-to", default=None, help="with --dry-run: the revision lines, for reading before publishing")
    ap.add_argument("--redo", default=None, help="a published revision set of this correction, made again (the correction "
                    "fixed since): every line of it is corrected anew from its line before the set, and the lines that "
                    "differ from the published ones are published as --set")
    a = ap.parse_args()
    org = a.organism
    if a.redo:
        return redo(org, a)
    sets = a.sets.split(",") if a.sets else None
    names = L.resolve_sets(FAMILY, org, sets, "revised")
    # the latest line of a qid: a revision set's line replaces the base set's (sets in the order of their seq)
    latest, base_sets, revision_sets = {}, [], []
    for n in L.resolve_sets(FAMILY, org, None, "revised"):
        if L.manifest(FAMILY, org, n).get("revises"):
            revision_sets.append(n)
            for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), n)):
                latest[r["qid"]] = r
        elif n in names:
            base_sets.append(n)
    names = base_sets
    # the sets this revision revises: the base sets read and the revision sets published on them (every qid of the
    # revision is published in one of them; naming only these keeps the check's qids to the sets read)
    revises = [n for n in L.resolve_sets(FAMILY, org, None, "revised") if n in set(base_sets) | set(revision_sets)]
    ctx, out = {}, []
    d = L.query_dir(FAMILY, org)
    for n in names:
        p = L._set_file(d, n)
        if a.correction in VERSION_CORRECTIONS:  # the version strata only
            recs = (latest.get(r["qid"], r) for r in L.iter_set(p) if VERSION_STRATA.match(r["stratum"]))
        elif a.correction == "database_continuity":  # the version strata and the lines read by start_release
            recs = (latest.get(r["qid"], r) for r in L.iter_set(p)
                    if VERSION_STRATA.match(r["stratum"]) or HELD_KIND.match(r.get("kind") or ""))
        else:
            recs = (latest.get(r["qid"], r) for r in L.iter_set(p)
                    if not r["stratum"].startswith("c6") and "_ver_" not in r["stratum"])
        k = len(out)
        out.extend(CORRECTIONS[a.correction](org, recs, ctx))
        log(n, "revised", len(out) - k)
    by = collections.Counter(r["stratum"] for r in out)
    log("revised", len(out), dict(by))
    if not out:
        log("nothing to revise")
        return
    if a.dry_run and a.lines_to:
        import gzip
        with gzip.open(a.lines_to, "wt") as fo:
            for r in out:
                fo.write(C.dumps(r) + "\n")
    doc = CORRECTIONS[a.correction].__doc__ or ""
    man = L.publish(FAMILY, org, a.set, out, f"revision `{a.correction}` (revise_c.py docstring)", ME,
                    revises=revises, judge="c-" + a.correction, sources=["core_xref/<release>.tsv.gz", "gtf/genes_<release>.tsv.gz"],
                    manifest_extra=dict(correction=a.correction, revision=doc.strip() or a.correction), dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k != "counts"})


def redo(org, a):
    """`--redo`: a published revision set of the correction made again from the lines before it (main's docstring)."""
    names = L.resolve_sets(FAMILY, org, None, "revised")
    before = names[:names.index(a.redo)]
    published = {r["qid"]: r for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), a.redo))}
    prior = L.load_qids(FAMILY, org, set(published), sets=before)
    for r in prior.values():
        r.pop("set", None)
    again = {r["qid"]: r for r in CORRECTIONS[a.correction](org, list(prior.values()), {})}
    keys = ("expect", "broad", "apart", "rule", "open", "revised_by", "other_carriers")
    out = []
    for q, old in published.items():
        new = again.get(q, prior[q])  # where the correction no longer applies, the line before it is the truth again
        if any((new.get(k) or None) != (old.get(k) or None) for k in keys):
            out.append(new)
    log("redo", a.redo, "lines", len(published), "differ", len(out),
        dict(collections.Counter((r["stratum"], q in again) for q, r in ((r["qid"], r) for r in out))))
    if not out:
        return
    doc = CORRECTIONS[a.correction].__doc__ or ""
    man = L.publish(FAMILY, org, a.set, out, f"revision `{a.correction}` made again ({a.redo}; revise_c.py docstring)", ME,
                    revises=names, judge="c-" + a.correction, sources=["core_xref/<release>.tsv.gz", "gtf/genes_<release>.tsv.gz"],
                    manifest_extra=dict(correction=a.correction, redo=a.redo, revision=doc.strip() or a.correction),
                    dry=a.dry_run, log=log)
    log("published", {k: v for k, v in (man or {}).items() if k != "counts"})


if __name__ == "__main__":
    main()
