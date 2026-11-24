"""F1's truth: what each feature is, and the harmonized output it implies (docs/families/f.md, "Truth").

Usage: truth_f.py --organism <o> --check [--set f1_v1]        the F lines whose source line has a newer truth
       truth_f.py --organism <o> --symbol-judge --name <set>   publish the `symbol` judge's lines of a set (below)
       truth_f.py --organism <o> --copies --name <set>         publish the revision `copies_taken` (below)

A feature's expected gene(s) at 115 are its line's `expect`, copied by the generator from the A or B line it was drawn
from (the latest line of that qid, so the revisions published before the set), and judged by the common scorer with
the source's rule. A revision of a source line published later is a revision set of F (`--check` lists the lines it
concerns). Nothing of idtrack is read.

`decided`: a feature's truth names one destination when it expects one gene (and its broad walk allows no other, and
it is not open), or nothing (and no place rule could accept a gene: no `locus`, no `broad`); every other feature (several
genes, a wider set, a place rule, open) is undetermined for the matrices and kept apart, never forced.

The expected outputs follow from the decided features by the rules `HarmonizeFeatures` documents (its class
docstring and `output_plan`; written here again, not imported, so that the package is checked against its own words):
a gene reached by one input feature is one output feature named by the gene; several input features on one gene: with
`merge="sum"` one feature whose values are their sum; with `merge="flag"` the input that names the gene itself (its
stable id; or, under the default final database, one of the gene's names there) is the gene's feature and the others
keep their input names as merge partners, and where none or several name it all keep their names and the gene has no
feature; every other input keeps its name (`keep_unconverted`); an unconverted name equal to a converted feature's name
(in the dataset, or in any dataset when merging) is renamed `unconverted:<name>`. Union: every feature of any dataset;
intersection: the features of every dataset; `uns`: the entries equal in every dataset.

Whether an input names the gene itself is decided by the truth for Ensembl ids (the id's stable id is the gene's) and
in the configuration `gene` (whose only name of a gene is its id); for a name under the default final database the
truth does not read the gene's names there, and a `flag` group holding a name and other inputs is undetermined.

The revision `copies_taken` (`f1_r1`): `f1_v1`'s generator mapped each folded string of the truth store to one of
the made-unique forms that fold to it, so a form folding like another (`X-1` and `X.1`) was not seen where only the
other's spelling is a string of the store, and some forms that are names of other genes (mouse `Ighv5-9-1`, a GTF name
beside `Ighv5-9`) were drawn as copies. The revision reads the store again for every `f1_made_unique` line (the
generator's `taken_forms`, corrected) and opens each line whose form is a string of the store: such a form is a name
of its own, read as itself by rule `made_unique_copies`, and its copy's truth does not hold.

The `symbol` judge: under the default final database the harmonizer restricts a tie to the genes with an id there; the
truth does not read which genes have one, so a feature expecting exactly several genes (rule `set_exact`, or the
support of rule `support_ranks_shared_name`) is open within them under that configuration (`--symbol-judge` publishes
those lines as a revision set; an answer outside the expected genes stays wrong).
"""
import argparse
import collections
import json
import sys

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

CONVERTED = ("converted", "converted_no_final_id")


def decided(line):
    """("gene", g) or ("nothing", None) where the truth names one destination; (None, why) otherwise."""
    if line.get("open"):
        return None, "open"
    E = list(line.get("expect") or [])
    B = line.get("broad")
    if B is not None and set(B) - set(E):
        return None, "a wider set is also allowed"
    if len(E) == 1:
        return "gene", E[0]
    if len(E) > 1:
        return None, "several genes"
    if line.get("locus") or line.get("locus_mapped") or line.get("locus_projected"):
        return None, "no successor: a gene at the place may be answered"
    return "nothing", None


def own_by_truth(feature, gene, config):
    """Whether the input names the gene itself (see the module docstring): True, False, or None (not read)."""
    if C.ENS_ID.match(feature):
        return C.base_id(feature) == C.base_id(gene)
    return False if config == "gene" else None


def plan(inputs, merge):
    """The output features of one dataset by the documented rules. `inputs`: [(name, gene or None, own, status)] in file
    order, `status` the input's own status (`converted`, `converted_no_final_id` for a gene; anything else: not
    converted, kept as `unconverted_<status>`), `own` True/False/None (None: not known, the group is `undetermined`
    under flag). Returns [{name, members, status, gene}] in the order of each output's first input."""
    groups = collections.OrderedDict()
    for i, (name, gene, own, st) in enumerate(inputs):
        if gene is not None and st in CONVERTED:
            groups.setdefault(gene, []).append(i)
    out, done = [], set()
    for i, (name, gene, own, st) in enumerate(inputs):
        if gene is not None and st in CONVERTED:
            members = groups[gene]
            suffix = "" if st == "converted" else "_no_final_id"
            if len(members) == 1:
                out.append(dict(name=gene, members=[name], status=st, gene=gene))
                continue
            if merge == "sum":
                if gene not in done:
                    done.add(gene)
                    out.append(dict(name=gene, members=[inputs[j][0] for j in members], status=f"merged{suffix}", gene=gene))
                continue
            owns = [inputs[j][2] for j in members]
            if any(o is None for o in owns):
                out.append(dict(name=name, members=[name], status="undetermined", gene=gene))
                continue
            mine = [j for j, o in zip(members, owns) if o]
            if len(mine) == 1 and mine[0] == i:
                out.append(dict(name=gene, members=[name], status=st, gene=gene))
            else:
                out.append(dict(name=name, members=[name], status="unconverted_merge_partner", gene=gene))
        else:
            out.append(dict(name=name, members=[name], status=f"unconverted_{st}", gene=None))
    converted = {o["name"] for o in out if not o["status"].startswith("unconverted") and o["status"] != "undetermined"}
    for o in out:
        if o["status"].startswith("unconverted") and o["name"] in converted:
            o["name"] = f"unconverted:{o['name']}"
    return out


def merged_names(per_dataset, features):
    """The merged object's feature names from each dataset's outputs ({dataset: [(name, status)]}): unconverted names
    equal to a converted name of any dataset renamed `unconverted:<name>`, then the union or the intersection.
    Returns ({dataset: {name as merged: name in the dataset}}, the merged names)."""
    converted = {n for outs in per_dataset.values() for n, st in outs if not str(st).startswith("unconverted")}
    renamed = {}
    for d, outs in per_dataset.items():
        renamed[d] = {(f"unconverted:{n}" if str(st).startswith("unconverted") and n in converted else n): n for n, st in outs}
    sets = [set(r) for r in renamed.values()]
    names = set().union(*sets) if features == "union" else (set.intersection(*sets) if sets else set())
    return renamed, names


def uns_kept(unses):
    """The `uns` of a merge whose entries are kept where they are the same in every dataset (recursively)."""
    if not unses:
        return {}
    out = {}
    for k in unses[0]:
        vals = [u.get(k, KeyError) for u in unses]
        if any(v is KeyError for v in vals):
            continue
        if all(isinstance(v, dict) for v in vals):
            sub = uns_kept(vals)
            if sub:
                out[k] = sub
        elif all(json.dumps(v, sort_keys=True) == json.dumps(vals[0], sort_keys=True) for v in vals):
            out[k] = vals[0]
    return out


def source_of(line):
    family, set_, qid = line["source"].split(":", 2)
    return family, set_, qid


TRUTH_FIELDS = ("expect", "broad", "rule", "cmp", "apart", "locus", "locus_mapped", "locus_projected", "open",
                "open_within", "support", "spellings", "names_read_at", "excluded")


def stale(org, set_name="f1_v1"):
    """The F lines whose source line has a newer truth than the one copied: [(F qid, source, fields that differ)]."""
    from idtrack_reproducibility.families.f import gen_f as G
    lines, _, _ = L.load("f", org, [set_name])
    want = collections.defaultdict(set)
    for q in lines:
        fam, _s, qid = source_of(q)
        want[fam].add(qid)
    latest = {}
    for fam, qids in want.items():
        sets = {"a": G.a1_sets(org), "b": G.b1_sets(org)}.get(fam)
        latest[fam] = L.load_qids(fam, org, qids, sets=sets) if fam != "harness" else \
            {q["qid"]: q for q in L.load("harness", org, kinds=["nonexistent"])[0] if q["qid"] in qids}
    out = []
    for q in lines:
        fam, _s, qid = source_of(q)
        src = latest[fam].get(qid)
        if src is None:
            out.append((q["qid"], q["source"], ["source line not found"]))
            continue
        diff = [k for k in TRUTH_FIELDS if json.dumps(src.get(k), sort_keys=True) != json.dumps(q.get(k), sort_keys=True)]
        if diff:
            out.append((q["qid"], q["source"], diff))
    return out


def symbol_judge_lines(lines):
    """The lines the `symbol` judge opens (module docstring)."""
    out = []
    for q in lines:
        E = list(q.get("expect") or [])
        if len(E) > 1 and (q.get("rule") == "set_exact" or q.get("support")) and not q.get("open"):
            r = {k: v for k, v in q.items() if k != "set"}
            r.update(open="under the default final database the harmonizer restricts a tie to the genes with an id "
                          "there, which the truth does not read", open_within=E,
                     revised_by=list(q.get("revised_by") or []) + ["symbol_tie"])
            out.append(r)
    return out


def copies_taken_lines(org, lines, log=print):
    """The `f1_made_unique` lines whose form is a string of the truth store, opened (the revision `copies_taken`)."""
    from idtrack_reproducibility.families.f import gen_f as G
    copies = [q for q in lines if q.get("stratum") == "f1_made_unique"]
    taken = G.taken_forms(org, {str(q["id"]) for q in copies}, log)
    out = []
    for q in copies:
        if str(q["id"]) in taken:
            r = {k: v for k, v in q.items() if k != "set"}
            r.update(open="the form is a string of the truth store (a name of its own, read as itself by rule "
                          "`made_unique_copies`): not a copy of the feature it extends",
                     revised_by=list(q.get("revised_by") or []) + ["copies_taken"])
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--set", default="f1_v1")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--symbol-judge", action="store_true")
    ap.add_argument("--copies", action="store_true")
    ap.add_argument("--name", default=None)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    if a.check:
        rows = stale(a.organism, a.set)
        for r in rows[:50]:
            print("\t".join([r[0], r[1], ",".join(r[2])]))
        print(f"{a.organism}: {len(rows)} lines of {a.set} whose source has a newer truth", file=sys.stderr)
        return 1 if rows else 0
    if a.symbol_judge:
        lines, _, _ = L.load("f", a.organism, [a.set])
        rev = symbol_judge_lines(lines)
        if a.dry:
            print(len(rev), collections.Counter(q["stratum"] for q in rev))
            return 0
        man = L.publish("f", a.organism, a.name, rev, f"F1's `symbol` judge: the lines of {a.set} expecting exactly several "
                        "genes, open within them under the default final database (truth_f.py docstring)", __file__,
                        revises=[a.set], judge="symbol", fmt="jsonl.gz")
        print(json.dumps({k: man[k] for k in ("name", "n", "counts")}))
        return 0
    if a.copies:
        lines, _, _ = L.load("f", a.organism, [a.set])
        rev = copies_taken_lines(a.organism, lines, lambda *x: print(*x, file=sys.stderr))
        print(len(rev), [q["qid"] for q in rev][:10], file=sys.stderr)
        if a.dry or not rev:
            return 0
        man = L.publish("f", a.organism, a.name, rev, f"F1's revision `copies_taken` of {a.set}: made-unique forms that are "
                        "strings of the truth store, opened (truth_f.py docstring)", __file__, revises=[a.set],
                        judge=None, fmt="jsonl.gz")
        print(json.dumps({k: man[k] for k in ("name", "n", "counts")}))
        return 0
    ap.error("--check, --symbol-judge or --copies")


if __name__ == "__main__":
    sys.exit(main())
