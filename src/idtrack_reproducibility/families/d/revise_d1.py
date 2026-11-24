"""Revision sets of D1 (a published set is never rewritten: the changed lines under their published qids).

Usage: revise_d1.py --organism <o> --revision <name> --set <new set name> [--dry]

Revisions:

- `name_several_set`: a symbol that several genes carry as their name in R's GTF (`truth_class` name_several) was
  published with `set_exact` (the harness's rule for a name of several genes: every view must return them all). The
  matrix's feature is one of those genes (a counts matrix holds one row per name; R's and Seurat's make.unique name
  the others `X.1`, which are the `made_unique` lines) and which one cannot be told from the name: `best` answering any
  non-empty part of the expected genes is right, `all` must still return them all: rule `set`. Found on d4c2dda's
  first run (MALAT1 at 77: the gene and a 95-base gene inside it both named MALAT1).
- `names_broad` (the corrected judge's `names_broad`, rule `judge_corrections`; the accuracy harness's revision of that
  name and family B's truth read names so): a symbol that other genes carry as a name in another name source at the
  release read cannot tell those genes from the genes whose GTF name it is, so the genes those carriers become at 115
  are `broad` (an answer among them is undetermined), or, on an open line with `open_within`, within the place (an
  answer among them is undetermined, not wrong). The name sources at R (`other_names`): BioMart's name tables (HGNC,
  MGI, UniProt gene name, NCBI label, WikiGene), the VGNC symbol and the Vega gene name where an archive answers (the
  harness's `Ctx78.name_genes`), and Ensembl's core cross-references of R in the name databases (`truth_b.NAME_DB`:
  HGNC, MGI, VGNC, UniProt gene name, NCBI gene, WikiGene, Vega gene, clone-based names), at every release; synonyms
  are no names (rule `names_before_synonyms`), the string is read exactly as written, and only genes of R's primary
  assembly count (copies on patches and haplotypes are rule `judge_copies_apart`'s). A carrier becomes at 115 what
  Ensembl's history walks it to (the broad walk with living ids, as the harness's `names_broad` walks it). Only the
  symbol lines change (`d1_symbol`; D3's lines by `revise_d3 --revision names_broad`).
"""
import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = 115  # the release every D line is asked to (gen_d1.TARGET)
from idtrack_reproducibility.families import common as C_
from idtrack_reproducibility.families import lines as L

HERE = os.path.dirname(os.path.abspath(__file__))


def name_several_set(q):
    if q.get("truth_class") == "name_several" and q.get("rule") == "set_exact":
        return dict(q, rule="set")
    return None


def other_names(org, R, C):
    """{string: {gene ids}} of the strings genes carry as a name in a name source other than the GTF at R (the module's
    docstring, `names_broad`)."""
    from idtrack_reproducibility.families.b import truth_b as TB
    from idtrack_reproducibility.harness import make_queries as M
    out = collections.defaultdict(set)
    for grp, col in M.NAME_SOURCES:
        pp = C.bm(R, grp, col)
        if pp:
            for v, gs in pp[1].items():
                out[v] |= gs
    for pp in (C.vgnc(R), C.vega(R)):
        if pp:
            for v, gs in pp[1].items():
                out[v] |= gs
    cx = TB.core_names(org, R)
    if cx is not None:
        for v, gs in cx[0].items():
            out[v] |= gs
    return {v: {C_.base_id(g) for g in gs} for v, gs in out.items()}


def names_broad(org, lines, kind, genes_field, read_field, log=print, tr=None, names_of=None):
    """`names_broad` (the module's docstring) on the latest lines of one kind (`d1_symbol`, `d3_symbol`): the revised
    lines. `genes_field` holds the genes the truth read the symbol as (the GTF's names at the release read),
    `read_field` the release it read. `tr` (the history: `genes`, `walk_fwd`) and `names_of` (R -> `other_names`) are
    read from the truth store unless given."""
    todo = collections.defaultdict(list)
    for q in lines:
        if q.get("kind") == kind and not q.get("excluded") and q.get(read_field) is not None and \
                (q.get("expect") or q.get("open_within") is not None):
            todo[q[read_field]].append(q)
    if not todo:
        return []
    if tr is None:
        from idtrack_reproducibility.families.astar.truth import TruthA
        tr = TruthA(org, TARGET)
    if names_of is None:
        from idtrack_reproducibility.harness import make_queries as M
        C = M.Ctx78(org, TARGET, truth=tr)
        names_of = lambda R: other_names(org, R, C)  # noqa: E731
    out, walks = [], {}
    for R in sorted(todo):
        names = names_of(R)
        # genes of R's primary assembly: copies on patches and haplotypes are the judge's `judge_copies_apart`, not names
        genes_r = {g for g, row in tr.genes.get(R, {}).items() if row.get("region") in (None, "primary", "core_only")}
        for q in todo[R]:
            others = (names.get(q["id"], set()) & genes_r) - {C_.base_id(g) for g in q.get(genes_field) or []}
            add = set()
            for g in sorted(others):
                if (g, R) not in walks:
                    walks[(g, R)] = {C_.base_id(x) for x in (tr.walk_fwd(g, R, TARGET, False, living=True)
                                                              if R < TARGET else {g})}
                add |= walks[(g, R)]
            new = None
            if q.get("expect"):
                base = set(q.get("broad") or q["expect"])
                if add - base:
                    new = dict(q, broad=sorted(base | add))
            elif add - set(q["open_within"]):
                new = dict(q, open_within=sorted(set(q["open_within"]) | add))
            if new is not None:
                new.pop("set", None)
                new["revised_by"] = list(q.get("revised_by") or []) + ["names_broad"]
                new["names_broad_genes"] = sorted(others)
                out.append(new)
        log(f"{org}: release {R}: {len(todo[R])} lines read, {sum(1 for q in out if q[read_field] == R)} revised")
    return out


REVISIONS = {"name_several_set": name_several_set}
# Revisions over every latest line at once (a function of the lines and the truth, not of one line):
TRUTH_REVISIONS = {"names_broad": lambda org, lines: names_broad(org, lines, "d1_symbol", "genes_at_source", "from_release")}
# Revisions that regenerate the lines with gen_d1.py's current truth and publish every line whose truth changed:
# - `mutual_first` (2026-11-19): the `dominant` class's `mutual` rule chose a small new gene lying inside a gene that
#   grew (CDC42 at 95: the model shares 3,598 bases with CDC42 at 115 and 2,914 with ENSG00000289694, which it holds
#   whole, while CDC42 holds under half of its own bases inside the model); the rule now also asks that the gene share
#   no fewer of the model's bases than any other gene. Lines whose truth changes are published again (dominant ->
#   open, mostly; a name of several genes whose genes' classes change). `name_several_set` is applied to them too.
REGENERATE = {"mutual_first"}
TRUTH_FIELDS = ("expect", "open", "open_within", "truth_class", "stratum", "rule")


def report(org, out, show=3):
    """The revised lines per stratum and what changed (broad or open_within), with a few examples."""
    by = collections.Counter((q["stratum"], "broad" if q.get("expect") else "open_within") for q in out)
    for (st, what), n in sorted(by.items()):
        print(f"  {st:32s} {what:12s} {n}")
    shown = collections.Counter()
    for q in out:
        if shown[q["stratum"]] < show:
            shown[q["stratum"]] += 1
            print("  e.g.", q["qid"], "genes", q.get("genes_at_source") or q.get("genes_at_read"), "expect", q.get("expect"),
                  "others", q["names_broad_genes"], "broad", q.get("broad"), "within", q.get("open_within"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--revision", required=True, choices=sorted(REVISIONS) + sorted(REGENERATE) + sorted(TRUTH_REVISIONS))
    ap.add_argument("--set", required=True)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--lines-to", default=None, help="with --dry: the revised lines, for reading before publishing")
    a = ap.parse_args()
    d1_sets = [s for s in L.set_names("d", a.organism) if s.startswith("d1")]
    lines, _, _ = L.load("d", a.organism, d1_sets)
    out = []
    if a.revision in REGENERATE:
        from idtrack_reproducibility.families.d import gen_d1
        latest = {q["qid"]: q for q in lines if q["qid"].startswith("d1|")}
        for q in gen_d1.build(a.organism, dry=True):
            q = name_several_set(q) or q
            p = latest.get(q["qid"])
            if p is None:
                raise SystemExit(f"{q['qid']}: not published; a revision only revises")
            if any(p.get(k) != q.get(k) for k in TRUTH_FIELDS):
                q["revised_by"] = list(p.get("revised_by") or []) + [a.revision]
                out.append(q)
        lines = list(latest.values())
    if a.revision in TRUTH_REVISIONS:
        out = TRUTH_REVISIONS[a.revision](a.organism, [q for q in lines if q["qid"].startswith("d1|")])
        report(a.organism, out)
    for q in (lines if a.revision in REVISIONS else []):
        if not q["qid"].startswith("d1|"):
            continue
        new = REVISIONS[a.revision](q)
        if new is not None:
            new.pop("set", None)
            new["revised_by"] = list(q.get("revised_by") or []) + [a.revision]
            out.append(new)
    if a.dry and a.lines_to:
        import gzip
        import json
        with gzip.open(a.lines_to, "wt") as fo:
            for q in out:
                fo.write(json.dumps(q, sort_keys=True, default=str) + "\n")
    print(f"{a.organism}: {len(out)} lines revised by {a.revision}")
    if out and not a.dry:
        revises = sorted({q["set"] for q in lines if q["qid"].startswith("d1|")})
        L.publish("d", a.organism, a.set, out, description=f"D1 revision {a.revision} (revise_d1.py docstring)",
                  generator=os.path.abspath(__file__), revises=revises,
                  judge="78" if a.revision in TRUTH_REVISIONS else None)
        print("published", a.set)


if __name__ == "__main__":
    main()
