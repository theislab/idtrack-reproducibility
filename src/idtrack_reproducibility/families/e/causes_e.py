"""A cause for every wrong answer and error of family E (the output format of F0's `../causes.py`).

Usage: causes_e.py --organism <o> --run <run>/<organism> [--judges now,decided]

F0's classifier reads the harness's kinds through explained paths; family E's wrong answers are of other kinds:
inputs the docs refuse or a decision reads otherwise. Its causes are read from the line (stratum, judge, the form),
the answer (status, match kind, the release the identifier was read at, the error class) and the truth store
(which releases have BioMart tables: a name read exactly at a release without them is one the truth cannot see).
Every row gets a cause; `unexplained` rows are read by hand and given one here.

Writes `<run>/scores_<judge>/causes.tsv.gz` and `causes_summary.tsv` (the common columns: stratum, cause, class,
rule, answers, units). Classes as F0's: a = idtrack wrong under its own decisions (here mostly a decision not
built yet: `decided` judge), b = idtrack follows its rule and the truth (the user's meaning) contradicts it,
c = the graph lacks what the truth has, d = the truth or the query is wrong, e = the truth cannot decide.
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import causes as K
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import lines as L

csv.field_size_limit(1 << 30)
FAMILY = "e"

# cause -> (class, decision or finding, one line)
CAUSES = {
    "missing_value_read_as_name": ("a", "rule `missing_values` (finding IN-1; decided)",
                                   "a spelling of a missing value is read as a name and answers a gene"),
    "missing_value_not_refused": ("a", "rule `missing_values` (finding IN-1; decided)",
                                  "a spelling of a missing value is looked up (not found) instead of refused"),
    "unicode_not_normalised": ("a", "rule `unicode_normalised` (finding IN-4; decided)",
                               "zero-width characters, a BOM, full-width letters or a Unicode dash: not normalised, not found"),
    "stray_punctuation_artifact_name": ("a", "rule `stray_punctuation` (finding IN-11; decided)",
                                        "a symbol typed with a stray `;`, `.`, `,` or quotes reaches another gene through a "
                                        "name of Ensembl's tables that carries the same stray character"),
    "stray_punctuation_not_read": ("a", "rule `stray_punctuation` (finding IN-11; decided)", "stray punctuation and quotes are not read away"),
    "separator_joins_suffix": ("b", "rule `made_unique_copies` (finding IN-3; decided)",
                               "a make-unique suffix (`-1`, `.1`, `_1`) is joined by the case or separator step into "
                               "another gene's name"),
    "suffix_form_is_name_at_unread_release": ("e", "-", "the suffixed form is itself a name (exactly) of another gene at a "
                                              "release whose name tables the truth does not read: it may be what the user meant"),
    "internal_name_accepted": ("a", "rule `no_synonym_target` (finding IN-6; decided)",
                               "the graph's bookkeeping name (`synonym_id::`) is read as an identifier"),
    "long_version_value_error": ("a", "rule `bounded_arguments` (finding IN-8; decided)",
                                 "a version of more than 4,300 digits raises Python's ValueError, not an input error"),
    "long_input_not_refused": ("a", "rule `bounded_arguments` (finding IN-8; decided)",
                               "an identifier or version longer than any real one is not refused with an input error"),
    "unversioned_tx_case_not_found": ("a", "rule `one_case_schemes` (finding IN-5; decided)",
                                      "an unversioned transcript or translation id in another letter case is not found"),
    "case_exact_synonym_first": ("a", "rule `single_case_symbols` (finding IN-2; decided)",
                                 "a mouse symbol in capitals or lower case is an exact synonym of another gene, read "
                                 "before the symbol in another case"),
    "lowercase_exact_synonym_other_organism": ("b", "IN-2's residual (rule `single_case_symbols` covers organisms with mixed-case symbols only)",
                                               "a symbol typed in another case is exactly a synonym of another gene, read first"),
    "source_inference_reads_names_at_snapshot": ("a", "rule `inference_per_release` (finding IN-7; decided)",
                                                 "source inference reads every name at the snapshot: a list of an old "
                                                 "release is placed at a later one"),
    "synonym_database_offered": ("a", "rule `no_synonym_target` (decided)", "a synonym database (`synonym_id::`) is accepted as a target"),
    "separator_group_one_taken": ("a", "rule `one_case_schemes` (finding IN-10; decided)",
                                  "a separator reading that reaches several names takes one silently"),
    "nonalnum_name_in_graph": ("a", "rule `missing_values`, builder part (decided)",
                               "a string without a letter or digit (`.`, `-`, `_`) reaches a gene through a placeholder "
                               "synonym of Ensembl's tables (EntrezGene's `.`), by the separator step"),
    "mouse_hgnc_projection": ("d", "harness_lib.EXCLUDED (mouse HGNC Symbol); the corrected judge",
                              "the mouse gene's display name is a human HGNC symbol projected onto it; the graph does not "
                              "hold HGNC symbols for mouse, the truth reads BioMart's display name"),
    "infer_duplicates_counted": ("a", "docs/methods.md E8 (the relation); docs: 'each identifier counts once'",
                                 "a repeated identifier counts once per occurrence in source inference; a list with "
                                 "duplicates can be placed at another release than the same list without them"),
    "line_premise_wrong": ("d", "revision e5_r1", "the line assumed a name of the organism's symbol database that the truth "
                           "does not hold (revised)"),
    "lrg_gene_version_form": ("d", "-", "a GTF of NCBI36 lists LRG genes as `LRG_<n>_g1`; the graph writes them `LRG_<n>_g1.1`; "
                              "the judge's unversioned comparison does not take the `.1` off such ids: the same gene"),
    "docs_case_rule_judge_short": ("e", "-", "judge `now` reads the docs' rule (an exact synonym before the case step) from the "
                                   "tables of 115 only; idtrack reads every release's synonyms, so its tie holds other genes "
                                   "too (the MGI symbol's gene is in neither): the judge cannot reproduce the docs' answer"),
    "unexplained": ("?", "-", "read by hand"),
}


def answers_of(run, qids):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    out = {}
    for p in glob.glob(os.path.join(run, "answers", "*.parquet")):
        t = pq.read_table(p, columns=["qid", "view", "status", "match", "got", "from_release", "error_type", "outcome"])
        t = t.filter(pc.is_in(t.column("qid"), value_set=pa.array(sorted(qids))))
        for r in t.to_pylist():
            out[(r["qid"], r["view"])] = r
    return out


def classify(line, row, ans, judge, bm_rels):
    st = row["stratum"]
    got = [x for x in (row.get("got") or "").split(",") if x]
    a = ans or {}
    m = a.get("match") or ""
    o = row["outcome"]
    if st.startswith("e1_") or st == "hunt_IN-1":
        return ("missing_value_read_as_name", f"answered {','.join(got[:3])} ({m})") if got else \
            ("missing_value_not_refused", a.get("status") or o)
    if st.startswith(("e2_zw_", "e2_fullwidth_", "e2_dash_")) or (st == "hunt_IN-4" and not got):
        if not got:
            return "unicode_not_normalised", line.get("note") or ""
    if st.startswith(("e2_punct_", "e2_all_punct")) or st == "hunt_IN-11" or (st == "hunt_IN-4" and got):
        if got:
            return "stray_punctuation_artifact_name", f"{line['id']!r} -> {','.join(got[:3])} ({m} at {a.get('from_release')})"
        return "stray_punctuation_not_read", line.get("note") or ""
    if st.startswith(("e2_makeunique_", "e2_all_makeunique")) or st == "hunt_IN-3":
        rel = a.get("from_release")
        if m in ("exact", "synonym:exact") and rel not in bm_rels:
            return "suffix_form_is_name_at_unread_release", f"{line['id']} {m} at {rel}"
        return "separator_joins_suffix", f"{line['id']} -> {','.join(got[:3])} ({m} at {rel})"
    if st.startswith("e2_internal_") or st == "hunt_IN-6":
        return "internal_name_accepted", f"{line['id'][:40]} -> {','.join(got[:3])} ({m})"
    if st == "e2_longversion" or st == "hunt_IN-8":
        if o == "error" and "ValueError" in (row.get("reason") or "") or a.get("error_type") == "ValueError":
            return "long_version_value_error", (row.get("reason") or "")[:100]
        return "long_input_not_refused", f"version of {len(line['id'].rsplit('.', 1)[-1])} digits"
    if st.startswith("e2_long"):
        return "long_input_not_refused", f"{len(line['id'])} characters"
    if st.startswith("e2_txcase_") or st == "hunt_IN-5":
        return "unversioned_tx_case_not_found", line["id"]
    if st.startswith("e2_case_") and judge == "now":
        return "docs_case_rule_judge_short", f"{line['id']} -> {','.join(got[:3])} ({m}); judge expected {','.join(line.get('expect', [])[:3])}"
    if st.startswith("e2_case_") or (st == "hunt_IN-2" and C.ENS_ID.match(got[0] if got else "") and "ENSMUSG" in (got[0] if got else "")):
        return "case_exact_synonym_first", f"{line['id']} -> {','.join(got[:3])} ({m})"
    if st == "hunt_IN-2":
        return "lowercase_exact_synonym_other_organism", f"{line['id']} -> {','.join(got[:3])} ({m})"
    if st == "hunt_IN-7":
        return "source_inference_reads_names_at_snapshot", f"release {line['expect']} -> {','.join(got[:5])}"
    if st == "hunt_IN-10":
        return "separator_group_one_taken", f"{line['id']} -> {','.join(got[:3])}"
    if st.startswith("e2_junk_") and got and m.endswith("separator") and not any(ch.isalnum() for ch in line["id"]):
        return "nonalnum_name_in_graph", f"{line['id']!r} -> {','.join(got[:3])} ({m})"
    if st.startswith("e4_symbol_") and not got and line.get("expect") and C.ENS_ID.match(line["expect"][0]) \
            and line["expect"][0].startswith("ENSMUSG") and line["id"].isupper():
        disp = [r for r in (H.biomart_genes("mus_musculus", 115) or []) if r["name"] == line["id"]]
        if disp and all(r["source"] == "HGNC Symbol" for r in disp):
            return "mouse_hgnc_projection", f"{line['id']} is the HGNC display name of {line['expect'][0]}"
    if st.startswith("e8_identity") and got == [f"{line['id']}.1"] and line["id"].startswith("LRG_"):
        return "lrg_gene_version_form", f"{line['id']} -> {got[0]}"
    if st == "e8_infer_duplicates":
        return "infer_duplicates_counted", (row.get("reason") or "")[:150]
    if st == "e5_database_spelling" and "from_database|name_db" in row["qid"]:
        return "line_premise_wrong", "the name is no symbol of the organism's symbol database (revision e5_r1)"
    if st == "e5_database_synonym":
        return "synonym_database_offered", line.get("kwargs", {}).get("final_database", "")
    return "unexplained", f"{st}: {(row.get('reason') or '')[:120]} ({m})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--judges", default="now,decided")
    a = ap.parse_args()
    bm_rels = set(H.biomart_releases(a.organism))
    for judge in a.judges.split(","):
        folder = os.path.join(a.run, f"scores_{judge}")
        rows = []
        for f in ("wrong.tsv.gz",):
            p = os.path.join(folder, f)
            if os.path.exists(p):
                rows += [dict(r, outcome="wrong") for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t", quoting=csv.QUOTE_NONE)]
        p = os.path.join(folder, "errors.tsv")
        if os.path.exists(p):
            rows += [dict(r, outcome="error", unit=r.get("unit", "")) for r in csv.DictReader(open(p), delimiter="\t", quoting=csv.QUOTE_NONE)]
        rj = json.load(open(os.path.join(folder, "run.json")))
        lines, _, _ = L.load(FAMILY, a.organism, list(rj["sets"]))
        by = {q["qid"]: q for q in lines}
        ans = answers_of(a.run, {r["qid"] for r in rows})
        summ = collections.defaultdict(lambda: [0, set()])
        with gzip.open(os.path.join(folder, "causes.tsv.gz"), "wt") as fo:
            fo.write("organism\tqid\tview\tkind\tstratum\tunit\toutcome\tcause\tclass\tdetail\n")
            for r in rows:
                q = by.get(r["qid"], {"id": "", "expect": []})
                view = (r.get("variant") or "best,").split(",")[0]
                cause, detail = classify(q, r, ans.get((r["qid"], view)), judge, bm_rels)
                unit = r.get("unit") or L.unit_of(q) if q.get("id") else ""
                fo.write("\t".join(str(x).replace("\t", " ").replace("\n", " ")[:400] for x in
                                   (a.organism, r["qid"], view, r["kind"], r["stratum"], unit, r["outcome"], cause,
                                    CAUSES[cause][0], detail)) + "\n")
                s = summ[(r["stratum"], cause)]
                s[0] += 1
                s[1].add(unit)
        K.write_summary(os.path.join(folder, "causes_summary.tsv"), summ, CAUSES)
        un = sum(n for (st, c), (n, _) in summ.items() if c == "unexplained")
        print(f"{a.organism} {judge}: {len(rows)} wrong or error rows; {un} unexplained")


if __name__ == "__main__":
    main()
