"""Family E's rows of the known-limitations table (docs/methods.md, "Known limitations"), generated from its scored runs.

Usage: limits_e.py --runs <run>/<organism> [...] --out <tsv>

Three sources, every number counted from the runs (never typed):
1. the causes of the wrong answers and errors (`causes_e.py`, judge `now`: what users get from the tree), mapped
   by `MAP` (cause -> limitation); a cause that a decided rule removes is written with that rule as `fixed_by`, so
   that the next tree measures it again;
2. what idtrack never reads (rule `safe`, answered nothing): forms of real identifiers that are not found (the
   share of such lines not answered), per group of strata;
3. what the docs leave open but users meet (`open` lines that named a gene: `undetermined_answered` of
   `e_summary.tsv`): missing-value spellings read as names, numbers of another organism read in every numeric
   database, symbols of another organism.

Columns: F0's (`../limitations.py`).
"""
import argparse
import csv
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import limitations as LIM
from idtrack_reproducibility.families.e import causes_e as K

HDR = ["limitation", "group", "title", "what", "organism", "answers", "units", "of_answers", "share", "family", "strata", "run",
       "commit", "graph_digest", "class", "cause", "rule", "user_action", "fixed_by", "measured_by", "since"]

# cause -> (limitation, group, title, what, user action, fixed by)
MAP = {
    "missing_value_read_as_name": ("L-missing-values", "inputs", "missing values read as names",
                                   "`nan`, `NA`, `N/A`, `None`, `null`, `NaT`, `unknown` (what pandas and R write for a missing "
                                   "value) are looked up as names and can answer a gene", "drop missing values before converting; "
                                   "pass `from_database`", "decided: rule `missing_values`"),
    "missing_value_not_refused": ("L-missing-values", "inputs", "missing values looked up", "as L-missing-values, not found",
                                  "drop missing values before converting", "decided: rule `missing_values`"),
    "unicode_not_normalised": ("L-unicode", "inputs", "invisible characters and Unicode look-alikes are not read",
                               "zero-width characters, a byte-order mark, full-width letters and Unicode dashes make an identifier "
                               "not found", "clean the strings (NFKC, strip format characters)", "decided: rule `unicode_normalised`"),
    "stray_punctuation_artifact_name": ("L-stray-punctuation", "inputs", "stray punctuation reaches another gene",
                                        "`X;`, `X.`, `'X'` match names of Ensembl's tables that carry the same stray character "
                                        "(UniProt gene names of release 78, HGNC symbols of release 50, MGI synonyms of 65-67)",
                                        "strip quotes and trailing punctuation", "decided: rule `stray_punctuation`"),
    "stray_punctuation_not_read": ("L-stray-punctuation", "inputs", "stray punctuation is not read away",
                                   "an identifier with quotes or a trailing `;`, `.`, `,` is not found", "strip them", "decided: rule `stray_punctuation`"),
    "separator_joins_suffix": ("L-make-unique", "inputs", "make-unique suffixes reach other genes",
                               "`X-1`, `X.1`, `X_1` (scanpy's and R's make-unique) are read through the case and separator "
                               "steps as another gene's name", "convert the names before making them unique", "decided: rule `made_unique_copies`"),
    "suffix_form_is_name_at_unread_release": ("L-truth-names-unread", "truth", "names of releases without BioMart tables",
                                              "the suffixed form is exactly a name of another gene at a release whose name tables "
                                              "the truth does not read; the judge cannot say which the user meant",
                                              "none needed", "the revision e2_names_every_release (names of every release "
                                              "from the core xrefs; synonyms of releases without a synonym table "
                                              "undetermined); the stored scorings predate it"),
    "internal_name_accepted": ("L-internal-names", "inputs", "the graph's internal names are read",
                               "`synonym_id::X` and `<id>.Void` are read as identifiers", "do not type them", "decided: rule `no_synonym_target`"),
    "long_version_value_error": ("L-long-inputs", "inputs", "very long versions raise a bare ValueError",
                                 "a version of more than 4,300 digits raises Python's ValueError", "none needed (no real input)",
                                 "decided: rule `bounded_arguments`"),
    "long_input_not_refused": ("L-long-inputs", "inputs", "inputs longer than any real identifier are looked up",
                               "no bound on the length of an identifier or a version", "none needed", "decided: rule `bounded_arguments`"),
    "unversioned_tx_case_not_found": ("L-transcript-case", "inputs", "unversioned transcript ids in another case",
                                      "`enst00000016171` is not found while `ENST00000016171` is", "upper-case Ensembl ids",
                                      "decided: rule `one_case_schemes`"),
    "case_exact_synonym_first": ("L-mouse-case", "names", "mouse symbols in capitals or lower case reach other genes",
                                 "an exact synonym of another gene is read before the MGI symbol in another case",
                                 "pass `from_database='MGI Symbol'`", "decided: rule `single_case_symbols`"),
    "lowercase_exact_synonym_other_organism": ("L-case-synonym", "names", "a symbol in another case that is another gene's synonym",
                                               "in organisms whose symbols are written in one case (pig, macaques, human), a "
                                               "symbol typed in another case that is exactly another gene's synonym answers that gene",
                                               "type symbols as the organism writes them, or pass `from_database`", "none (no rule decided)"),
    "source_inference_reads_names_at_snapshot": ("L-infer-names", "source inference", "lists of names placed at a later release",
                                                 "names renamed after the list's release vote for later releases",
                                                 "infer the release from Ensembl ids", "decided: rule `inference_per_release`"),
    "synonym_database_offered": ("L-internal-names", "inputs", "synonym databases accepted as targets",
                                 "`synonym_id::<database>` is accepted as `final_database`", "use the database itself", "decided: rule `no_synonym_target`"),
    "separator_group_one_taken": ("L-separator-groups", "names", "names equal up to separators",
                                  "a third spelling takes one of several names silently", "type the exact spelling", "decided: rule `one_case_schemes`"),
    # the limitation F0 measures on the harness's names (families/limitations.py): the same row
    "synonym_spelling_merged": ("L-spellings", "names", "two genes with one string in two letter cases",
                                "the exact step cannot tell `Tc` from `tc`", "type the exact spelling (decided: rule `spelling_per_edge`)",
                                "decided: rule `spelling_per_edge`"),
    "infer_duplicates_counted": ("L-infer-duplicates", "source inference", "repeated identifiers count once each time",
                                 "a list with repeated identifiers can be placed at another release than the same list without "
                                 "repeats", "remove repeated identifiers before inferring", "decided: rule `inference_distinct_ids`"),
    "lrg_record_not_its_gene": ("L-lrg-records", "identifiers", "LRG ids answer the LRG record",
                                "an LRG gene id is answered with the LRG record, not the Ensembl gene it describes",
                                "convert the Ensembl gene id", "decided: rule `lrg_answers_gene`"),
    "mouse_hgnc_projection": ("L-mouse-hgnc-names", "names", "human symbols projected onto mouse genes",
                              "a mouse gene whose display name is a human HGNC symbol (no MGI symbol) is not found by that name",
                              "use the Ensembl id", "none (the graph holds no HGNC symbols for mouse)"),
    "nonalnum_name_in_graph": ("L-missing-values", "inputs", "strings without a letter or digit reach genes",
                               "`_`, `-`, `.` are read through EntrezGene's placeholder synonym `.` by the separator step",
                               "drop such strings", "decided: rule `missing_values`, builder part"),
}

# strata whose `safe` lines answered nothing: forms never read (coverage, not wrong)
NOT_READ = {
    "L-format-characters": ("e2_bidi_", "inputs", "bidi marks and soft hyphens are not read",
                            "an identifier with a left-to-right or right-to-left mark, an isolate or a soft hyphen is not found",
                            "strip format characters", "none (no rule decided)"),
    "L-control-characters": ("e2_ctrl_", "inputs", "control characters are not read",
                             "an identifier with NUL, an escape sequence or another control character is not found",
                             "strip control characters", "none (no rule decided)"),
    "L-lookalikes": ("e2_lookalike_", "inputs", "letters of other scripts are not read",
                     "a Cyrillic or Greek look-alike letter makes an identifier not found", "none", "none (by design)"),
    "L-emoji-brackets": ("e2_emoji_|e2_brackets_|e2_innerspace_", "inputs", "emoji, brackets and inner spaces are not read",
                         "an identifier with an emoji, in brackets, or with a space inside is not found", "clean the strings",
                         "none (by design)"),
    "L-prefixes": ("e2_prefix", "inputs", "database prefixes are not read",
                   "`ensembl:ENSG...`, `UniProtKB:P04637`, `RefSeq:NM_...` are not found", "strip the prefix", "none (no rule decided)"),
    "L-number-forms": ("e2_numbers_ncbi", "external identifiers", "NCBI gene ids with zeros, as floats or in scientific notation",
                       "`07157`, `7157.0`, `7.157e+03`, `GeneID:7157` are not found", "write NCBI ids as integers", "none (no rule decided)"),
    "L-version-forms": ("e2_versionform_", "inputs", "versions written with `-`, `_` or `v`", "`ENSG...-5`, `ENSG..._5`, `ENSG....v5` are not found",
                        "write `<id>.<version>`", "none (no rule decided)"),
}

# strata whose `open` lines named a gene: what users meet where the docs leave it open
OPEN_ANSWERED = {
    "L-missing-values": ("e1_single_d52|e1_single_d81", "inputs", "missing values read as names (single calls)",
                         "the share of rule `missing_values`'s and 81's spellings that answer a gene on this tree",
                         "drop missing values before converting", "decided: rule `missing_values`"),
    "L-bare-numbers": ("e4_ncbi_foreign", "external identifiers", "numbers of another organism read in this organism's databases",
                       "an NCBI gene id of another organism answers a gene here through numeric ids of other databases (old Havana "
                       "numbers, RefSeq rows with NCBI numbers)", "pass `from_database` for NCBI gene ids",
                       "none (rule `bare_numbers`: a number is read in every numeric database, each named)"),
    "L-foreign-symbols": ("e4_symbol_absent", "names", "symbols of another organism",
                          "a symbol of another organism that no table of this organism names (that the truth reads) answers a gene",
                          "convert symbols in their own organism's graph", "none (the truth's reach)"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--judge", default="now")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for run in a.runs:
        folder = os.path.join(run, f"scores_{a.judge}")
        info = json.load(open(os.path.join(folder, "run.json")))
        org = info["organism"]
        common = dict(organism=org, family="e", run=LIM.where(info.get("run") or run), commit=(info.get("commit") or "")[:12],
                      graph_digest=",".join(info.get("graph_digests") or []), since=(info.get("commit") or "")[:12])
        rep = {}
        for r in csv.DictReader(open(os.path.join(folder, "report.tsv")), delimiter="\t"):
            if r["view"] == "*":
                rep[r["stratum"]] = int(r["n"])
        agg = {}
        p = os.path.join(folder, "causes.tsv.gz")
        if os.path.exists(p):
            for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t", quoting=csv.QUOTE_NONE):
                x = agg.setdefault(r["cause"], [0, set(), set()])
                x[0] += 1
                x[1].add(r["unit"])
                x[2].add(r["stratum"])
        for cause, (n, units, strata) in sorted(agg.items()):
            lim = MAP.get(cause, (f"UNMAPPED-{cause}", "", "UNMAPPED", "", "", ""))
            of = sum(rep.get(s, 0) for s in strata)
            cls, dec, _ = K.CAUSES.get(cause, ("?", "-", ""))
            rows.append(dict(common, limitation=lim[0], group=lim[1], title=lim[2], what=lim[3], answers=n, units=len(units),
                             of_answers=of, share=f"{n / of:.5f}" if of else "", strata=",".join(sorted(strata)),
                             **{"class": cls}, cause=cause, rule=dec, user_action=lim[4], fixed_by=lim[5],
                             measured_by=f"e: causes_e on {LIM.where(folder)}"))
        summ = list(csv.DictReader(open(os.path.join(run, "e_summary.tsv")), delimiter="\t"))
        summ = [r for r in summ if r["judge"] == a.judge]
        for lid, (pat, group, title, what, act, fixed) in NOT_READ.items():
            sel = [r for r in summ if any(r["stratum"].startswith(p) for p in pat.split("|"))]
            n = sum(int(r["n"]) for r in sel)
            if not n:
                continue
            not_read = sum(int(r["n"]) - int(r["answered"]) for r in sel)
            rows.append(dict(common, limitation=lid, group=group, title=title, what=what, answers=not_read, units="",
                             of_answers=n, share=f"{not_read / n:.5f}", strata=",".join(r["stratum"] for r in sel),
                             **{"class": "coverage"}, cause="not_read", rule="-", user_action=act, fixed_by=fixed,
                             measured_by=f"e: e_summary.tsv of {LIM.where(run)} (lines answered nothing)"))
        for lid, (pat, group, title, what, act, fixed) in OPEN_ANSWERED.items():
            sel = [r for r in summ if r["stratum"] in pat.split("|")]
            n = sum(int(r["n"]) for r in sel)
            if not n:
                continue
            ans = sum(int(r["undetermined_answered"]) for r in sel)
            rows.append(dict(common, limitation=lid, group=group, title=title, what=what, answers=ans, units="",
                             of_answers=n, share=f"{ans / n:.5f}", strata=",".join(r["stratum"] for r in sel),
                             **{"class": "open"}, cause="open_answered", rule="-", user_action=act, fixed_by=fixed,
                             measured_by=f"e: e_summary.tsv of {LIM.where(run)} (open lines that named a gene)"))
    with open(a.out, "w") as fo:
        fo.write("\t".join(HDR) + "\n")
        for r in rows:
            fo.write("\t".join(str(r[h]).replace("\t", " ") for h in HDR) + "\n")
    print(f"{len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
