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


CORRECTIONS = {"rfam_label_accession": rfam_label_accession, "other_database_broad": other_database_broad,
               "rfam_label_rows": rfam_label_rows, "par_shared_names": par_shared_names,
               "display_names_at_release": display_names_at_release, "no_mito_nuclear_step": no_mito_nuclear_step,
               "mirbase_genes": mirbase_genes}


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
    latest, base_sets = {}, []
    for n in L.resolve_sets(FAMILY, org, None, "revised"):
        if L.manifest(FAMILY, org, n).get("revises"):
            for r in L.iter_set(L._set_file(L.query_dir(FAMILY, org), n)):
                latest[r["qid"]] = r
        elif n in names:
            base_sets.append(n)
    names = base_sets
    ctx, out = {}, []
    d = L.query_dir(FAMILY, org)
    for n in names:
        p = L._set_file(d, n)
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
                    revises=L.resolve_sets(FAMILY, org, None, "revised"), judge="c-" + a.correction, sources=["core_xref/<release>.tsv.gz", "gtf/genes_<release>.tsv.gz"],
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
