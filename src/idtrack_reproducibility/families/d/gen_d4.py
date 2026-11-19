"""Family D4: mixed lists, as users paste them: Ensembl ids of several releases (current, retired, versioned, with a
version newer than the graph or unknown), transcript and translation ids, symbols (unique, of several genes, synonyms,
case and separator variants, a name that moved to another gene), NCBI gene ids, RefSeq and UniProt accessions, strings
that are nothing, missing values and duplicates, in one `convert_many` call to 115 without a release (a list of
several origins has none). Nothing of idtrack is read.

Usage: gen_d4.py --organism <o> [--publish <set>] [--dry] [--lists 2000]

**Items and their truth**: each item is a published line of the accuracy harness (`family="harness"`, the latest line
of its qid, the revised judge) asked without `from_release` to 115 and without a final database; the item's line in
D4 keeps that line's truth whole (`expect`, `broad`, `rule`, `cmp`, `apart`, `locus`, `locus_mapped`, `support`,
`spellings` and the rest), so F0's scorer judges the item exactly as the harness judges the single call (`src_qid`
names it). Missing values (`nan`, `NA`, `None`, `null`, ``, ` `) expect rule `missing_values`'s error for that item
(`expect_error` MissingValueError, d4c2dda's ValueError), as family E's `e1_batch_d52`. The harness's kind
`nonexistent_version` (superseded by rule `version_newer_than_graph`) is not drawn. Strata: `d4_<category>` of the item; `d4_issue85` is
issue #85's list for human (`TP53`, `ENSG00000141510`, `NOT_A_REAL_GENE`).

**Lists**: 2,000 per organism (seeded, `random.Random("<organism>|<list>")`): a size of 5, 10, 20, 50 or 100; the
categories' weights drawn per list (so some lists are almost all of one kind, others even mixtures); with probability
0.3 one to three items repeated at other positions, and with probability 0.3 one or two missing values inserted.
Each item of a list is one line (`call: convert_many`, `batch`, `pos`; family E's adapter), so a list's answer is read
item by item.
"""
import argparse
import collections
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

FAMILY = "d"
CATEGORIES = {
    "id_current": ("gene_unversioned_to_self", "gene_versioned_to_self"),
    "id_old": ("retired_renamed", "retired_merged", "retired_split", "retired_none", "retired_none_branch",
               "version_newer_than_graph", "version_unknown"),
    "transcript": ("transcript_to_gene", "translation_to_gene"),
    "symbol": ("name_unique_to_gene", "name_ambiguous_to_set", "synonym_to_gene", "case_variant", "separator_variant",
               "name_changed_gene"),
    "ncbi": ("ncbi_to_gene",),
    "refseq": ("refseq_to_gene",),
    "uniprot": ("uniprot_to_gene",),
    "nothing": ("nonexistent",),
}
TRUTH_KEYS = ("expect", "broad", "rule", "cmp", "apart", "locus", "locus_mapped", "locus_projected", "exceptions", "support",
              "spellings", "revised_by", "open", "open_within")
MISSING = ["nan", "NA", "None", "null", "", " "]
SIZES = (5, 10, 20, 50, 100)
EXPECT_MISSING = "MissingValueError|d4c2dda:ValueError"


def pools(org):
    lines, _, _ = L.load("harness", org)
    kind_cat = {k: c for c, ks in CATEGORIES.items() for k in ks}
    out = collections.defaultdict(list)
    for q in lines:
        c = kind_cat.get(q["kind"])
        if c is None or q.get("from_release") is not None or int(q.get("to_release") or 0) != 115 or q.get("final_database"):
            continue
        out[c].append(q)
    for c in out:
        out[c].sort(key=lambda q: q["qid"])
    return out


def item(q, cat):
    d = {k: q.get(k) for k in TRUTH_KEYS if q.get(k) is not None}
    d.setdefault("expect", [])
    return dict(d, id=q["id"], category=cat, src_qid=q["qid"], src_kind=q["kind"], unit=C.base_id(q["id"]))


def make_lists(org, P, n_lists):
    cats = sorted(P)
    for k in range(n_lists):
        rnd = random.Random(f"{org}|{k}")
        size = rnd.choice(SIZES)
        w = [rnd.random() ** 2 for _ in cats]
        items = []
        for _ in range(size):
            c = rnd.choices(cats, weights=w)[0]
            items.append(item(rnd.choice(P[c]), c))
        if rnd.random() < 0.3:
            for _ in range(rnd.randint(1, 3)):
                x = dict(rnd.choice(items), duplicate=True)
                items.insert(rnd.randrange(len(items) + 1), x)
        if rnd.random() < 0.3:
            for _ in range(rnd.randint(1, 2)):
                items.insert(rnd.randrange(len(items) + 1), dict(id=rnd.choice(MISSING), category="missing", expect=[], rule="set",
                                                                 cmp="versioned", expect_error=EXPECT_MISSING, unit="missing"))
        yield k, items


def lines_of(org, k, items, stratum_override=None):
    batch = [x["id"] for x in items]
    comp = dict(collections.Counter(x["category"] for x in items))
    for pos, x in enumerate(items):
        q = dict(x)
        cat = q.pop("category")
        q.update(qid=f"d4|{k}|{pos}", family=FAMILY, kind=f"d4_{cat}", stratum=stratum_override or f"d4_{cat}",
                 call="convert_many", batch=batch, pos=pos, from_release=None, to_release=115, final_database=None,
                 strategy="both", list=k, list_size=len(batch), composition=comp)
        q["id"] = str(x["id"])
        q.setdefault("broad", None)
        yield q


def build(org, publish=None, dry=False, n_lists=2000):
    P = pools(org)
    print(org, {c: len(v) for c, v in sorted(P.items())})
    out = []
    for k, items in make_lists(org, P, n_lists):
        out += list(lines_of(org, k, items))
    if org == "homo_sapiens":
        names = {q["id"]: q for q in P.get("symbol", []) if q["kind"] == "name_unique_to_gene"}
        ids = {q["id"]: q for q in P.get("id_current", []) if q["kind"] == "gene_unversioned_to_self"}
        junk = dict(id="NOT_A_REAL_GENE", category="nothing", expect=[], rule="set", cmp="versioned", unit="NOT_A_REAL_GENE",
                    note="issue #85's third item: a string that is nothing (no harness line; no gene, name or id is it)")
        if "TP53" in names and "ENSG00000141510" in ids:
            out += list(lines_of(org, "issue85", [item(names["TP53"], "symbol"), item(ids["ENSG00000141510"], "id_current"), junk],
                                 stratum_override="d4_issue85"))
        else:
            print("issue #85's items are not both harness lines: TP53", "TP53" in names, "ENSG00000141510", "ENSG00000141510" in ids)
    c = collections.Counter(q["stratum"] for q in out)
    print(org, len(out), "lines;", dict(sorted(c.items())))
    if publish and not dry:
        L.publish(FAMILY, org, publish, out, description="D4: mixed lists in one convert_many call to 115; each item judged by its "
                  "harness line's truth (gen_d4.py docstring)", generator=os.path.abspath(__file__),
                  sources=["the accuracy harness's published sets (revised judge)"], manifest_extra=dict(lists=n_lists))
        print("published", publish)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--publish")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--lists", type=int, default=2000)
    a = ap.parse_args()
    build(a.organism, a.publish, a.dry, a.lists)
