"""The reading of the sampled disagreements of family G (docs/families/g.md, "The reading of the disagreements").

Usage: python -m idtrack_reproducibility.families.g.reading [--organisms ...] [--tree main] [--out <folder>]

Reads `<families>/g/report/disagreements.tsv` (the seeded sample, `sampled`), G's lines, the answers of each tool
(form best: the release it answered at), the truth store (Ensembl's GTF, names, synonyms and cross-references of
that release, Ensembl's history) and, for idtrack's wrong lines, the causes of the families' own runs of the tree
(`<families>/<source family>/runs/<tree>/<organism>/causes.tsv.gz`, or a scoring folder's) with their rows of the
known-limitations table (`results/limitations/limitations.tsv`). Nothing is fetched. Writes into
`<families>/g/report/`:

  reading.tsv          every sampled line: subject, organism, task, qid, the two verdicts, the answer of the side read
                       (`got`: the tool's, or idtrack's), the class, its group and its evidence (for idtrack the judge's
                       reason), and its weight (the lines of its sampling cell, subject x organism x task x outcome
                       pair, per sampled line)
  reading_classes.tsv  per side (tool, idtrack), subject, task and class: its group, the sampled lines, their share of
                       the side's sampled lines, the weighted share (the estimated share of all the subject's
                       disagreements of that side and task, every organism together), three examples
  reading_summary.tsv  per subject and task: the sampled lines where the tool is wrong and the weighted shares of the
                       groups below; the sampled lines where idtrack is wrong and the weighted share with a known
                       limitation

A sampled line where the tool is wrong is classed by what the tool did, read in this order:

  the service could not read the input   the tool answered an error for the input (e.g. Ensembl REST's 400 on a
                                         Greek letter)
  another annotation's genes only        the tool named only genes of another Ensembl annotation of the species (pig's
                                         breed assemblies), which its adapter leaves out (`other_annotation_only`)
  not a gene of the release              an id the tool returned is no gene of the release it answers at in Ensembl's
                                         GTF (a retired id, another assembly's or annotation's id)
  a gene on a patch                      a returned gene lies on a patch or alternative locus of that release
  a part of the truth's genes            the tool returned some of the truth's genes, not all (the view `all`)
  the authority's gene                   the returned gene is the one the authority (HGNC, MGI, VGNC, NCBI, UniProt)
                                         names for the input, which the line records; the truth reads Ensembl's own
                                         tables (principle 2: their disagreement is Ensembl's, not the tool's fault)
  a reading of another release           the line asks a name, id or accession as an earlier release read it; at the
                                         tool's release Ensembl's own tables read it as the tool does (another gene
                                         or none): the tool answers what it is built for, today's reading
  another gene the string names          at the tool's release Ensembl's tables also name the returned gene by the
                                         string (a shared name, a synonym, a cross-reference): a valid reading the
                                         judge does not allow
  a successor Ensembl's history allows   (an Ensembl id) the returned gene is reached by Ensembl's history from the
                                         line's release by the wide walk, not by the truth's
  no answer                              the tool named nothing where the truth names a gene, and Ensembl's tables of
                                         the tool's release name it by the input (or its history reaches it): a miss
  a synonym, which the tool does not read  (a tool that reads display names only: Ensembl REST's `lookup/symbol`) the
                                         truth's gene has the input as a synonym only
  no answer, not in Ensembl's tables     the tool named nothing, and Ensembl's tables of its release do not name the
                                         truth's gene by the input either (the truth reads another source)
  truly wrong                            another current gene, by none of the readings above

The classes fall in the groups of the comparison's reading (`GROUPS`): `truly wrong` (truly wrong; no answer),
`right under another reading` (another annotation's genes, an id of another release, a gene on a patch, a part of the
truth's genes, the authority's gene, another gene the string names, a successor Ensembl's history allows), `not a
question the tool claims` (a reading of another release; a synonym a display-name lookup does not read; no answer where
Ensembl's tables of the tool's release do not name the gene), `the service failed` (the service could not read the
input).

A sampled line where idtrack is wrong is classed by its cause in the families' runs and that cause's row of the
known-limitations table (`L-...`; the family's and organism's row, else the row of any family that states the cause;
the cause `unexplained` names none), else by the source family's map of its wrong answers' classes
(`limitations_map_outcomes.tsv`, `limits_outcomes.classify`) on idtrack's answer. A line no known limitation states is
read here by what idtrack did (`idtrack_class`, against Ensembl's tables of the release the line asks to): a string
that stands for a missing value, refused (`L-missing-values`); a symbol in another case that is another gene's synonym
(the answer's match `synonym:case`) or a string no Ensembl table of any release holds, answered (the graph holds
Ensembl's strings only: one held in another case) (`L-case-synonym`); else a new cause (`new: ...`): a gene where the truth names none, the truth's copies on patches or
haplotypes in the answer, a name on patches answered with another patch's gene, another gene.
"""
import argparse
import collections
import csv
import gzip
import json
import os

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families.g import fetch as F
from idtrack_reproducibility.families.g import likeforlike as LL
from idtrack_reproducibility.harness import harness_lib as H

TARGET = 115

GROUPS = {
    "truly wrong": "truly wrong", "no answer": "truly wrong",
    "not a gene of the release": "right under another reading", "a gene on a patch": "right under another reading",
    "another annotation's genes only": "right under another reading",
    "a part of the truth's genes": "right under another reading", "the authority's gene": "right under another reading",
    "another gene the string names": "right under another reading",
    "a successor Ensembl's history allows": "right under another reading",
    "a reading of another release": "not a question the tool claims",
    "a synonym, which the tool does not read": "not a question the tool claims",
    "no answer, not in Ensembl's tables": "not a question the tool claims",
    "the service could not read the input": "the service failed",
}


def releases_of(subject, organism):
    """{qid: (release, the tool's status)} of a subject's answers (form best)."""
    import pyarrow.parquet as pq
    d = os.path.join(C.family_dir("g"), "runs", subject, "best", organism, "answers")
    out = {}
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        if f.endswith(".parquet"):
            t = pq.read_table(os.path.join(d, f), columns=["qid", "flags", "status"]).to_pylist()
            for r in t:
                if r["flags"]:
                    out[r["qid"]] = (json.loads(r["flags"]).get("release"), r["status"])
    return out


class At:
    """Ensembl's tables of one release, as the reading uses them."""

    def __init__(self, organism, rel):
        self.rel = rel
        self.gtf = H.gtf_genes(organism, rel) or {}
        self.names = LL.Release(organism, rel)

    def genes_named(self, line, synonyms=True):
        v = line["id"]
        if line["task"] == "e":
            return self.names.acc.get(v.split(".")[0] if v[:3] in ("NM_", "XM_", "NR_", "XR_", "NP_", "XP_") else v, set())
        k = v.casefold()
        return self.names.names.get(k, set()) | (self.names.synonyms.get(k, set()) if synonyms else set()) | \
            self.names.acc.get(v, set())


DISPLAY_NAMES_ONLY = {"ensembl_rest_lookup"}  # `lookup/symbol` reads display names, never synonyms (its documentation)


def tool_class(line, got, tool_outcome, at, truth, status=None, synonyms=True):
    exp = set(line.get("expect") or [])
    g = {x for x in got if x}
    if tool_outcome == "error":
        return "the service could not read the input", ""
    if not g and status == "other_annotation_only":
        return "another annotation's genes only", "the tool named genes of another annotation of the species only"
    if g and any(x not in at.gtf for x in g):
        return "not a gene of the release", ",".join(sorted(x for x in g if x not in at.gtf)[:3])
    if g and any((at.gtf.get(x) or {}).get("region") == "patch" for x in g):
        return "a gene on a patch", ",".join(sorted(x for x in g if (at.gtf.get(x) or {}).get("region") == "patch")[:3])
    if g and g < exp:
        return "a part of the truth's genes", f"{len(g)} of {len(exp)}"
    auth = {C.base_id(x) for x in (line.get("authority") or line.get("auth") or []) if x and x != "none"}
    if g and g <= auth and not g <= exp:
        return "the authority's gene", f"{line.get('authority_db') or 'the authority'} names {','.join(sorted(g)[:3])}"
    named = at.genes_named(line) if line["task"] in ("c", "d", "e") else set()
    # the line reads its input at a release other than the target and the tool's: its own (from_release) or where the
    # truth read it (read_at); a line read at the target and kept by the rule of equal truth reads alike at the tool's
    old = line.get("to_release") == TARGET and (line.get("from_release") not in (None, TARGET, at.rel)
                                                 or line.get("read_at") not in (None, TARGET, at.rel))
    if line["task"] in ("c", "d", "e"):
        if g and g <= named:
            if old and not (exp & named):
                return "a reading of another release", f"at {at.rel} the input names {','.join(sorted(named)[:3])}"
            return "another gene the string names", f"at {at.rel} the input also names {','.join(sorted(g)[:3])}"
        if not g:
            if exp & named and not synonyms and not exp & at.genes_named(line, synonyms=False):
                return "a synonym, which the tool does not read", f"at {at.rel} {','.join(sorted(exp & named)[:3])} has the input as a synonym"
            if exp & named:
                return "no answer", f"at {at.rel} Ensembl names {','.join(sorted(exp & named)[:3])} by the input"
            if old:
                return "a reading of another release", f"at {at.rel} the input names {','.join(sorted(named)[:3]) or 'no gene'}"
            return "no answer, not in Ensembl's tables", "Ensembl's tables of the release do not name the truth's gene by the input either"
    if line["task"] in ("a", "b") and truth is not None and line.get("from_release") in truth.genes:
        w = truth.walk_series(C.base_id(line["id"]), line["from_release"], [at.rel], False).get(at.rel, set())
        if g and g <= w:
            return "a successor Ensembl's history allows", f"the wide walk to {at.rel} reaches {','.join(sorted(g)[:3])}"
        if not g:
            return "no answer", f"Ensembl's history reaches {','.join(sorted(w)[:3]) or 'nothing'} at {at.rel}"
    if not g:
        return "no answer", ""
    return "truly wrong", f"returned {','.join(sorted(g)[:3])}; the truth {','.join(sorted(exp)[:3]) or 'nothing'}"


def idtrack_class(line, got, at, match=None, reason=None):
    """What idtrack did on a wrong line that no known limitation states (Ensembl's tables of the release the line asks
    to, `at`):
    (class, the limitation it falls under, or None for a cause the table does not state)."""
    exp = set(line.get("expect") or [])
    g = {C.base_id(x) for x in got if x}
    if "MissingValueError" in (reason or ""):
        return "a string that stands for a missing value, refused (rule `missing_values`)", "L-missing-values"
    if (match or "").startswith("synonym:case"):
        return "a symbol in another case that is another gene's synonym", "L-case-synonym"
    region = lambda xs: {(at.gtf.get(x) or {}).get("region", "absent") for x in xs}  # noqa: E731
    if line["task"] in ("c", "d") and g and not exp:
        # the graph holds Ensembl's strings only: an answer to a string no Ensembl table of any release holds is a
        # string held in another case (idtrack reads names without case where the organism writes them in one case)
        if (line.get("note") or "").startswith("no Ensembl table of any release holds the string"):
            return "a string Ensembl's tables hold only in another case, read as that name", "L-case-synonym"
        return "a gene where the truth names none", None
    missing, extra = exp - g, g - exp
    if g and missing and not extra and region(missing) <= {"patch"}:
        # idtrack answers the primary assembly's gene and reports copies apart (`alternative_locus_genes`)
        return "the truth's copies on patches or haplotypes in the answer (idtrack answers the primary gene)", None
    if extra and missing and region(extra | missing) <= {"patch"}:
        return "a name on patches answered with another patch's gene", None
    return "another gene", None


def causes_of(source_family, organism, tree):
    """{qid: (cause, class)} of the families' run of the tree (view all; a family that scores in folders: `scores_now`)."""
    out = {}
    base = os.path.join(C.family_dir(source_family), "runs", tree, organism)
    for p in (os.path.join(base, "causes.tsv.gz"), os.path.join(base, "scores_now", "causes.tsv.gz")):
        if os.path.exists(p):
            with gzip.open(p, "rt", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                    if r["view"] == "all":
                        out[r["qid"]] = (r["cause"], r["class"])
    return out


UNEXPLAINED = "unexplained"


def limitations():
    """{(family, organism, cause): limitation} and {cause: limitation} of the known-limitations table: a cause is one
    word of the families' classifier, and the row that states it in one family states it for every family (an
    `unexplained` cause names no limitation)."""
    from idtrack_reproducibility.check import ROOT
    p = os.path.join(ROOT, "results", "limitations", "limitations.tsv")
    out, by_cause = {}, {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t", quoting=csv.QUOTE_NONE):
            if not r["cause"] or r["cause"] == UNEXPLAINED or UNEXPLAINED in r["limitation"]:
                continue
            out.setdefault((r["family"], r["organism"], r["cause"]), r["limitation"])
            by_cause.setdefault(r["cause"], r["limitation"])
    return out, by_cause


def outcome_maps():
    """{family: the rows of its map of wrong answers' classes} (`families/<family>/limitations_map_outcomes.tsv`, read by
    `limits_outcomes.classify`: a class the answer's own fields name, by the line's source stratum)."""
    from idtrack_reproducibility.families import limits_outcomes as LO
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = {}
    for fam in os.listdir(here):
        p = os.path.join(here, fam, "limitations_map_outcomes.tsv")
        if os.path.exists(p):
            out[fam] = LO.load_map(p)
    return out


def idtrack_outcomes(organism):
    """{qid: (status, match, reason, got)} of idtrack's scoring on G's lines (view all)."""
    import pyarrow.parquet as pq
    d = os.path.join(C.family_dir("g"), "runs", "idtrack", "stated", organism, "outcomes")
    out = {}
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.endswith(".parquet"):
                cols = ["qid", "view", "status", "match", "reason", "got"]
                for r in pq.read_table(os.path.join(d, f), columns=cols).to_pylist():
                    if r["view"] == "all":
                        got = r["got"] if isinstance(r["got"], str) else ",".join(r["got"] or [])
                        out[r["qid"]] = (r["status"], r["match"], r["reason"], got)
    return out


def run(organisms, tree="main", out_dir=None):
    import pandas as pd
    from idtrack_reproducibility.harness import make_queries as M
    rep = os.path.join(C.family_dir("g"), "report")
    out_dir = out_dir or rep
    d = pd.read_csv(os.path.join(rep, "disagreements.tsv"), sep="\t", keep_default_na=False)
    cell = ["subject", "organism", "task", "tool", "idtrack"]
    d["sampled"] = d.sampled.astype(str) == "True"
    w = d.groupby(cell).agg(n=("qid", "size"), k=("sampled", "sum"))
    weight = {key: n / k for key, n, k in zip(w.index, w.n, w.k) if k}
    d = d[d.sampled]
    lim, lim_by_cause = limitations()
    from idtrack_reproducibility.families import limits_outcomes as LO
    omaps = outcome_maps()
    rows = []
    for org in organisms:
        sub = d[d.organism == org]
        if not len(sub):
            continue
        lines = {q["qid"]: q for q in F.g_lines(org, "local")}
        truth = M.Truth(org, 116)
        ats, rel_by = {}, {}
        causes, idto = {}, None
        for subj in sorted(set(sub.subject)):
            rel_by[subj] = releases_of(subj, org)
        for r in sub.itertuples():
            ln = lines.get(r.qid)
            if ln is None:
                continue
            got = [x for x in str(r.got).split(",") if x]
            base = dict(subject=r.subject, organism=org, task=r.task, qid=r.qid, tool=r.tool, idtrack=r.idtrack,
                        disagreement=r.cls, got=",".join(got[:5]), expected=",".join((ln.get("expect") or [])[:5]),
                        stratum=ln.get("stratum"), weight=weight[(r.subject, org, r.task, r.tool, r.idtrack)])
            if r.cls in ("the tool wrong", "both wrong"):
                rel, status = rel_by[r.subject].get(r.qid) or (TARGET, None)
                rel = rel or TARGET
                if rel not in ats:
                    ats[rel] = At(org, rel)
                cls, ev = tool_class(ln, got, r.tool, ats[rel], truth, status, r.subject not in DISPLAY_NAMES_ONLY)
                rows.append(dict(base, side="tool", release=rel, cls=cls, group=GROUPS[cls], evidence=ev))
            if r.cls in ("idtrack wrong", "both wrong"):
                src = ln["source_family"]
                if src not in causes:
                    causes[src] = causes_of(src, org, tree)
                c = causes[src].get(r.qid)
                if idto is None:
                    idto = idtrack_outcomes(org)
                st, mt, rs, ig = idto.get(r.qid, (None, None, None, ""))
                L_ = None
                if c:
                    L_ = lim.get((src, org, c[0])) or lim_by_cause.get(c[0])
                if not L_ and src in omaps:  # the family's map of its wrong answers' classes, on idtrack's G answer
                    m = LO.classify(omaps[src], ln.get("source_stratum") or "", dict(status=st, match=mt, reason=rs))
                    if m is not None:
                        L_ = m["limitation"]
                        c = (m["when"], m["class"]) if not c or c[0] == UNEXPLAINED else c
                if c and L_:
                    cls = f"{L_}: {c[0]}"
                else:  # read here: what idtrack did, under a known limitation where one states it, else a new cause
                    rel_i = ln.get("to_release") or TARGET  # the release the line asks to (115, or its own)
                    if rel_i not in ats:
                        ats[rel_i] = At(org, rel_i)
                    k, L_ = idtrack_class(ln, [x for x in ig.split(",") if x], ats[rel_i], mt, rs)
                    why = f"{c[0]} in the families' runs" if c else "no cause in the families' runs"
                    cls = f"{L_}: {k}" if L_ else f"new: {k} ({why})"
                rows.append(dict(base, side="idtrack", release=TARGET, cls=cls, got=",".join(ig.split(",")[:5]),
                                 group="a known limitation" if L_ else "no known limitation",
                                 evidence=rs or (c[1] if c else "")))  # idtrack's answer and the judge's reason
        print(org, "read", sum(1 for x in rows if x["organism"] == org), flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(out_dir, "reading.tsv"), sep="\t", index=False)
    out, summ = [], {}
    for (side, subj, task), g in t.groupby(["side", "subject", "task"]):
        n, wn = len(g), g.weight.sum()
        for cls, h in g.groupby("cls"):
            ex = h.head(3)
            out.append(dict(side=side, subject=subj, task=task, cls=cls, group=h.group.iloc[0], lines=len(h),
                            share=round(len(h) / n, 4), weighted_share=round(h.weight.sum() / wn, 4),
                            examples=" | ".join(f"{x.organism}:{x.qid} got {x.got or 'nothing'} expected {x.expected or 'nothing'}"
                                                f"{(' (' + x.evidence + ')') if x.evidence else ''}" for x in ex.itertuples())))
        s = summ.setdefault((subj, task), dict(subject=subj, task=task))
        if side == "tool":
            s["tool_wrong_sampled"] = n
            for grp in dict.fromkeys(GROUPS.values()):
                s[grp] = round(g.weight[g.group == grp].sum() / wn, 4)
        else:
            s["idtrack_wrong_sampled"] = n
            s["idtrack_known_limitation"] = round(g.weight[g.group == "a known limitation"].sum() / wn, 4)
    pd.DataFrame(out).sort_values(["side", "subject", "task", "lines"], ascending=[True, True, True, False]) \
        .to_csv(os.path.join(out_dir, "reading_classes.tsv"), sep="\t", index=False)
    sm = pd.DataFrame(list(summ.values())).sort_values(["subject", "task"])
    for col in ("tool_wrong_sampled", "idtrack_wrong_sampled"):
        if col in sm:
            sm[col] = sm[col].fillna(0).astype(int)
    sm.to_csv(os.path.join(out_dir, "reading_summary.tsv"), sep="\t", index=False)
    return len(t)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.reading", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--tree", default="main")
    ap.add_argument("--out", default=None, help="the folder of the tables (default `<families>/g/report/`)")
    a = ap.parse_args(argv)
    print(run(CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(","), a.tree, a.out))


if __name__ == "__main__":
    main()
