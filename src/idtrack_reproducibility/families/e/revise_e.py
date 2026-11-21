"""Revision sets of family E (a published set is never rewritten; a change of truth is a revision set).

Usage: revise_e.py --organism <o> --revision <name> [--dry]

- `e5_r1`: the line `e|e5_database_spelling|from_database|name_db` assumed that the name chosen for E5 is a symbol of
  the organism's symbol database (`gen_e.NAME_DB`); the truth says whether it is (BioMart 115: `hgnc`, `mgi`; VGNC's
  own sets through BioMart's VGNC ids). Where it is not (cynomolgus: USP24 is a UniProt gene name, no VGNC symbol), the
  relation has no premise: the line is `open` with the reason.
- `e5_r2`: the same line as a plain line (rule `set`, open): F0's scorer applies `open` to every rule but `relation`.
- `e1_dec2` (judge `decided`): rule `missing_values` as built: the strings pandas reads as missing by default (its own
  list) are refused too, so E1's lines of such a spelling outside the earlier members (`<NA>`, `#N/A`, `#NA`: strata
  `e1_single_other`, `e1_batch_other`) expect `MissingValueError` under the decided judge, as the members do; the
  lines with the source database named stay open (looked up there).
- `e2_names_every_release` (both judges: the truth's reach): a made-unique form (`X-1`, `X.1`, `X_1`: strata
  `e2_makeunique_*`, `e2_all_makeunique`, `hunt_IN-3`) that Ensembl gives exactly, as a value of any database or as the
  GTF's gene name, to other genes at a release of the window names those genes too (idtrack reads an identifier at
  the latest release where it exists, rule `start_release`). The truth read names at BioMart's releases only; it now
  reads every release's core cross-references (`core_xref/<release>.tsv.gz`) and GTF: those genes (walked to the
  target) join `expect` (rule `safe`: the string names them). Synonyms are read only at the releases that have
  BioMart's synonym table or a core synonym extract (`core_synonyms/<release>.tsv.gz`): every such line records them
  (`names_read_at`), and the scorer calls an answer that read the form itself as a synonym (`synonym:exact`) at another release undetermined.
  Published as `<set>_r1` (judge `now`) and `<set>_dec2` (the decided judge's line, revised alike).
- `e2_case_r1` (judge `now`): the documented rule of mouse's letter case (an exact name or synonym before a match
  without letter case; `e2_case`) names every gene that carries the form exactly at the release read (115: lines
  without a release): the generator read BioMart's synonyms (those of a gene's display cross-reference only) and left
  the MGI symbol's own gene out of the exact carriers. The revision reads Ensembl's synonyms of every cross-reference
  of a name database at 115 (`core_synonyms/115.tsv.gz`, fetch_b.py synonyms) beside them: a line whose form some
  gene carries exactly expects every such gene (the symbol's own gene too); `open` is left only where no gene
  carries it exactly. Published as `e2_case_r1`.
- `stray_decided` (judge `decided`; published as `e1_dec3`, `e2_dec3`): rule `stray_punctuation` (decided, built in
  idtrack's next version) composed with the others: punctuation that ends a sentence and quotes at the ends of a string
  are no part of an identifier, so a string that is a missing value once they are cut (`None.`; strata
  `e1_single_other`, `e1_batch_other`, no source database named) expects `MissingValueError`, as the missing value
  itself does; and E2's junk strings whose cut form is a name (`CLEC1?`: stratum `e2_junk_wild`, the wildcard `?` is
  terminal punctuation) expect the genes the harness's revision `stray_punctuation` gives them (`make_queries`).
- `e8_identity_dec` (judge `decided`): rule `lrg_answers_gene` (decided): an LRG id answers the Ensembl gene
  the LRG record describes. A line asked by an LRG gene id (`LRG_<n>_g<k>`) at a release expects, under the decided
  judge, the gene that Ensembl's core database of that release links to the record (the gene's xref of database
  `ENS_LRG_gene`, "LRG display in Ensembl", whose accession is the LRG gene id: `core_xref/<release>.tsv.gz`); a line
  whose release has no such link is left as it is.
"""
import argparse
import collections
import csv
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.e import gen_e as G


REL_AS_SET = False


def e5_r1(org):
    qid = "e|e5_database_spelling|from_database|name_db"
    lines, _, _ = L.load("e", org, ["e5"])
    q = next((x for x in lines if x["qid"] == qid), None)
    if q is None:
        return []
    db = G.NAME_DB[org]
    if db == "VGNC Symbol":
        t = H.vgnc_symbols(org, 115)
        holds = bool(t and t[1].get(q["id"]))
    else:
        t = H.biomart_pairs(org, 115, "hgnc" if db == "HGNC Symbol" else "mgi", 1)
        holds = bool(t and t[1].get(q["id"]))
    if holds:
        return []
    r = dict(q, open=f"{q['id']} is not a {db} of this organism at 115 (the line's premise)", revised_by=["e5_r1"])
    if REL_AS_SET:
        # e5_r2: F0's scorer reads `open` before every rule but `relation`: the line becomes a plain line, open
        r.update(rule="set", expect=[], revised_by=["e5_r1", "e5_r2"])
        r.pop("relation", None)
    r.pop("set", None)
    return [r]


LRG_GENE = re.compile(r"^(LRG_\d+_g\d+)(\.\d+)?$")


def lrg_genes(org, rel):
    """{LRG gene id: {Ensembl gene id (versioned)}} of a release, from Ensembl's core xrefs (`ENS_LRG_gene`); None when
    the truth store has no core xrefs of that release."""
    p = os.path.join(H.truth_dir(org), "core_xref", f"{rel}.tsv.gz")
    if not os.path.exists(p):
        return None
    out = {}
    with gzip.open(p, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if r["db_name"] == "ENS_LRG_gene" and r["object_type"] == "Gene":
                out.setdefault(r["acc"], set()).add(r["gene_id"])
    return out


def e8_identity_dec(org):
    """The decided judge's lines of rule `lrg_answers_gene` (see the module's docstring) and the truth files read."""
    lines, _, _ = L.load("e", org, ["e8_identity"])
    out, files, cache = [], set(), {}
    for q in lines:
        m = LRG_GENE.match(q["id"] or "")
        if not m:
            continue
        rel = q.get("from_release") or q.get("to_release")
        if rel not in cache:
            cache[rel] = lrg_genes(org, rel)
        genes = (cache[rel] or {}).get(m.group(1))
        if not genes:
            continue
        files.add(os.path.join(H.truth_dir(org), "core_xref", f"{rel}.tsv.gz"))
        expect = sorted(C.base_id(g) for g in genes) if q.get("cmp") == "base" else sorted(genes)
        r = dict(q, expect=expect, revised_by=["decided"],
                 note=((q.get("note") or "") + "; " if q.get("note") else "") + "the gene the LRG record describes "
                      "(rule `lrg_answers_gene`; Ensembl's xref ENS_LRG_gene)")
        r.pop("set", None)
        out.append(r)
    return out, sorted(files)


def e1_dec2(org):
    """The decided judge's E1 lines of rule `missing_values` for pandas' spellings (see the module's docstring)."""
    from idtrack_reproducibility.harness import make_queries as MQ
    lines, _, _ = L.load("e", org, ["e1", "e1_dec"])
    out = []
    for q in lines:
        if q.get("stratum") not in ("e1_single_other", "e1_batch_other") or not MQ.is_missing_value(q.get("id") or ""):
            continue
        r = dict(q, open=None, expect_error=G.E_MISSING, revised_by=["decided"],
                 note=((q.get("note") or "") + "; " if q.get("note") else "") + "a missing value of pandas' list (rule "
                      "`missing_values`)")
        r.pop("set", None)
        out.append(r)
    return out


def stray_decided(org):
    """{set name: revised lines} of `stray_decided` (the module's docstring), from the decided judge's latest lines."""
    from idtrack_reproducibility.families.e import score_e as SE
    from idtrack_reproducibility.harness import make_queries as MQ
    lines, _, _ = L.load("e", org, SE.sets_for(org, "decided", only=["e1", "e2"]))
    out = collections.defaultdict(list)
    junk = []
    for q in lines:
        st, s = q.get("stratum"), q.get("id")
        if not isinstance(s, str) or MQ.strip_stray_punctuation(s) == s.strip():
            continue
        if st in ("e1_single_other", "e1_batch_other") and not q.get("from_database") and \
                MQ.is_missing_value(MQ.strip_stray_punctuation(s)) and q.get("expect_error") != G.E_MISSING:
            r = dict(q, open=None, expect_error=G.E_MISSING, revised_by=list(q.get("revised_by") or []) + ["stray_decided"],
                     note=((q.get("note") or "") + "; " if q.get("note") else "") + "a missing value once stray "
                          "punctuation is cut (rules `stray_punctuation`, `missing_values`)")
            r.pop("set", None)
            out["e1_dec3"].append(r)
        elif st == "e2_junk_wild":
            junk.append(q)
    if junk:
        C = MQ.Ctx78(org, H.ORGANISMS[org][2])
        for v in MQ.rev_stray_punctuation(org, junk, C).values():
            for r in v:
                r["revised_by"] = [x for x in r["revised_by"] if x != "stray_punctuation"] + ["stray_decided"]
                out["e2_dec3"].append(r)
    return out


def e2_case_r1(org):
    """The `now` judge's `e2_case` lines that Ensembl's synonyms of every cross-reference at 115 change (module
    docstring), and the truth files read."""
    from idtrack_reproducibility.families.b import truth_b as TB
    from idtrack_reproducibility.families.e import truth_e as T
    lines, _, _ = L.load("e", org, ["e2_case"])
    if not lines:
        return [], []
    dec = {q["qid"]: q for q in L.load("e", org, ["e2_case", "e2_case_dec"])[0]}
    names = T.Names(org)
    core = TB.core_synonyms(org, 115) or {}
    cur = T.Current(org)
    out = []
    for q in lines:
        f = str(q["id"])
        exact = set(names.exact115.get(f, set())) | {cur.versioned.get(g, g) for g in core.get(f, set())}
        if not exact:
            continue
        own = (dec.get(q["qid"]) or {}).get("expect") or []
        want = sorted(exact)
        if want == sorted(q.get("expect") or []) and not q.get("open"):
            continue
        r = dict(q, expect=want, open=None, revised_by=list(q.get("revised_by") or []) + ["e2_case_r1"],
                 note=((q.get("note") or "") + "; " if q.get("note") else "") +
                 f"genes carrying the form exactly at 115 (BioMart and Ensembl's synonyms of every xref)"
                 + ("; the MGI symbol's own gene among them" if set(own) & exact else ""))
        r.pop("set", None)
        out.append(r)
    return out, [os.path.join(H.truth_dir(org), "core_synonyms", "115.tsv.gz")]


MAKE_UNIQUE = ("e2_makeunique_", "e2_all_makeunique", "hunt_IN-3")


def e2_names_every_release(org, log=print):
    """{set to publish: (revised lines, the set it revises)} of `e2_names_every_release` (module docstring)."""
    from idtrack_reproducibility.families.a import gen_a as GA
    from idtrack_reproducibility.families.c import truth_c as TC
    from idtrack_reproducibility.families.e import score_e as SE
    now = {q["qid"]: q for q in L.load("e", org, SE.sets_for(org, "now"))[0] if q["stratum"].startswith(MAKE_UNIQUE)}
    dec = {q["qid"]: q for q in L.load("e", org, SE.sets_for(org, "decided"))[0] if q["qid"] in now}
    forms = {str(q["id"]) for q in now.values()}
    hits = collections.defaultdict(lambda: collections.defaultdict(set))  # form -> release -> genes
    xdir = os.path.join(H.truth_dir(org), "core_xref")
    rels = sorted(int(f.split(".")[0]) for f in os.listdir(xdir) if f.endswith(".tsv.gz")) if os.path.isdir(xdir) else []
    for rel in rels:
        for _ot, _oid, gid, _db, _disp, acc, lab, _v, _it in TC.read_xref(org, rel):
            for v in (acc, lab):
                if v in forms:
                    hits[v][rel].add(gid)
        for g, r in (H.gtf_genes(org, rel) or {}).items():
            if r.get("gene_name") in forms:
                hits[r["gene_name"]][rel].add(g)
    log(org, "core xref releases", len(rels), "forms", len(forms), "named at some release", len(hits))
    syn_rels = sorted({r for r in H.biomart_releases(org) if os.path.exists(os.path.join(H.truth_dir(org), "biomart",
                                                                                            str(r), "synonym.tsv"))} |
                      {int(f.split(".")[0]) for f in os.listdir(os.path.join(H.truth_dir(org), "core_synonyms"))
                       if f.endswith(".tsv.gz")} if os.path.isdir(os.path.join(H.truth_dir(org), "core_synonyms")) else set())
    log(org, "releases whose synonyms the truth reads", syn_rels)
    G = GA.GenA(org, 115) if hits else None
    out = collections.defaultdict(list)
    for qid, q in sorted(now.items()):
        h = hits.get(str(q["id"])) or {}
        t = q["to_release"]
        ext = set()
        for rel, genes in h.items():
            for g in genes:
                if g not in (G.tr.genes.get(rel) or {}):
                    continue
                E, B = G.walks(g, rel, [t])
                ext |= set(G.tr.vset(E.get(t, set()) | B.get(t, set()), t))
        for which, line in (("now", q), ("decided", dec.get(qid, q))):
            add = ext - set(line.get("expect") or [])
            r = dict(line, names_read_at=syn_rels, revised_by=list(line.get("revised_by") or []) + ["e2_names_every_release"],
                     note=((line.get("note") or "") + "; " if line.get("note") else "") +
                     "the form names other genes at a release of the window (core xrefs, GTF)")
            if line.get("rule") == "safe":
                r["expect"] = sorted(set(line.get("expect") or []) | add)
            else:
                r["broad"] = sorted(set(line.get("broad") or []) | set(line.get("expect") or []) | add)
            base = q["set"]
            r.pop("set", None)
            name = f"{base}_r1" if which == "now" else f"{base}_dec2"
            out[(name, base, which)].append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--revision", required=True, choices=["e5_r1", "e5_r2", "e1_dec2", "e8_identity_dec",
                                                         "e2_names_every_release", "e2_case_r1", "stray_decided"])
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    if a.revision == "stray_decided":
        for name, rev in sorted(stray_decided(a.organism).items()):
            if name in L.set_names("e", a.organism):
                print(a.organism, name, "published")
                continue
            base = name.rsplit("_dec", 1)[0]
            revises = [x for x in L.set_names("e", a.organism) if x == base or x.startswith(base + "_")]
            m = L.publish("e", a.organism, name, rev, "E: stray punctuation cut before the other decided rules (judge: "
                          "the decided rule `stray_punctuation` as built)", __file__, revises=revises, judge="decided",
                          sources=["make_queries.strip_stray_punctuation", "the harness's truth (name genes)"], dry=a.dry)
            print(a.organism, name, m["n"] if m else len(rev), collections.Counter(r["stratum"] for r in rev))
        return
    if a.revision == "e2_names_every_release":
        from idtrack_reproducibility.families.e import score_e as SE
        for (name, base, which), rev in sorted(e2_names_every_release(a.organism).items()):
            if name in L.set_names("e", a.organism):
                print(a.organism, name, "published")
                continue
            revises = [base] if which == "now" else \
                [base] + [s for s in L.set_names("e", a.organism) if SE.DEC.match(s) and SE.DEC.match(s).group("base") == base]
            m = L.publish("e", a.organism, name, rev, "E2 made-unique forms that name other genes at a release of the "
                          "window (the truth reads every release's core xrefs and GTF)", __file__, revises=revises,
                          judge=which, sources=[os.path.join(H.truth_dir(a.organism), "core_xref"),
                                                os.path.join(H.truth_dir(a.organism), "gtf")], dry=a.dry)
            print(a.organism, name, m["n"])
        return
    if a.revision == "e2_case_r1":
        rev, files = e2_case_r1(a.organism)
        if not rev:
            print(a.organism, a.revision, "no line changes")
            return
        import collections as _c
        print(a.organism, a.revision, len(rev), dict(_c.Counter(r["stratum"] for r in rev)))
        m = L.publish("e", a.organism, a.revision, rev, "E2 mouse case: every gene carrying the form exactly at 115 "
                      "(Ensembl's synonyms of every cross-reference beside BioMart's)", __file__, revises=["e2_case"],
                      judge="now", sources=files, dry=a.dry)
        print(a.organism, a.revision, m["n"])
        return
    if a.revision == "e1_dec2":
        rev = e1_dec2(a.organism)
        if not rev:
            print(a.organism, a.revision, "no line changes")
            return
        m = L.publish("e", a.organism, a.revision, rev, "E1: pandas' spellings of a missing value refused (judge: the "
                      "decided rule `missing_values` as built)", __file__, revises=["e1", "e1_dec"], judge="decided",
                      sources=["pandas._libs.parsers.STR_NA_VALUES"], dry=a.dry)
        print(a.organism, a.revision, m["n"])
        return
    if a.revision == "e8_identity_dec":
        rev, files = e8_identity_dec(a.organism)
        if not rev:
            print(a.organism, a.revision, "no line changes")
            return
        m = L.publish("e", a.organism, a.revision, rev, "E8 identity: an LRG id answers the gene its record describes "
                      "(judge: the decided rule `lrg_answers_gene`)", __file__, revises=["e8_identity"],
                      judge="decided", sources=files, dry=a.dry)
        print(a.organism, a.revision, m["n"])
        return
    global REL_AS_SET
    REL_AS_SET = a.revision == "e5_r2"
    rev = e5_r1(a.organism)
    if not rev:
        print(a.organism, a.revision, "no line changes")
        return
    m = L.publish("e", a.organism, a.revision, rev, "E5: a relation whose premise the truth does not hold, open", __file__,
                  revises=["e5"] + (["e5_r1"] if a.revision == "e5_r2" else []), judge="now", sources=[os.path.join(H.truth_dir(a.organism), "SOURCES.tsv")], dry=a.dry)
    print(a.organism, a.revision, m["n"])


if __name__ == "__main__":
    main()
