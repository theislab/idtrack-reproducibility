"""Scoring of family G: each tool's answers from its stored responses, scored by the common scorer on G's sets
(docs/families/g.md, "Scoring").

Usage: python -m idtrack_reproducibility.families.g.score --subjects <s,...|all> --organisms <o,...|all>
                                                         [--forms best,stated] [--tree <tree>]

Per subject, form of the input and organism: the answers (`<families>/g/runs/<subject>/<form>/<organism>/answers/
shard_0_of_1.0.parquet`, the runner's schema): for a line the tool cannot pose, outcome `not_supported` with the
reason; for a line it was asked, its adapter's answer (adapters.py) under both views (a tool gives one set: `best`
judges it as a part of the truth, `all` as the whole); a line asked of no stored response has no answer (`missing`).
The scorer (`families/scorer.py`) writes the record's rows beside them, and `run.json` names the subject (`kind =
tool`, `name`, `version`, `asked_at`, `endpoint`). The subject `idtrack` is the tree's answers to the same lines from
the families' own runs (`<families>/<source family>/runs/<tree>/<organism>/answers/`, read only; the lines' qids are
the families'), scored on G's lines like every tool.
"""
import argparse
import glob
import json
import os
import types

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families.g import adapters as A
from idtrack_reproducibility.families.g import fetch as F
from idtrack_reproducibility.families.g import local as LO
from idtrack_reproducibility.families.g import tools as T

WEB = ("ensembl_rest", "ensembl_rest_lookup", "biomart", "gprofiler", "mygene", "hgnc", "uniprot", "bridgedb", "gget",
       "seurat")
LOCAL = ("pyensembl", "orgdb", "limma", "hgnchelper", "ensdb", "ncbi_files")
SUBJECTS = WEB + LOCAL
ENDPOINT = {"ensembl_rest": "https://rest.ensembl.org (archive/id, xrefs/symbol)",
            "ensembl_rest_lookup": "https://rest.ensembl.org (lookup/symbol)",
            "biomart": "https://<archive>.archive.ensembl.org/biomart/martservice",
            "gprofiler": "https://biit.cs.ut.ee/gprofiler/api/convert/convert/", "mygene": "https://mygene.info/v3/query",
            "hgnc": "https://rest.genenames.org/fetch", "uniprot": "https://rest.uniprot.org/idmapping",
            "bridgedb": "https://webservice.bridgedb.org/<organism>/xrefsBatch", "gget": "gget info (Ensembl REST lookup)",
            "seurat": "Seurat UpdateSymbolList (HGNC REST)", "pyensembl": "pyensembl on Ensembl's GTF",
            "orgdb": "org.*.eg.db (AnnotationDbi::mapIds)", "limma": "limma alias2Symbol",
            "hgnchelper": "HGNChelper checkGeneSymbols", "ensdb": "AnnotationHub EnsDb 113 (ensembldb)",
            "ncbi_files": "NCBI gene files of the truth store"}


def raw_tool(subject):
    return "ensembl_rest" if subject == "ensembl_rest_lookup" else subject


def versions_by_release(subject, organism):
    """{release: [version folders]} of a subject's stored responses for an organism (from each folder's passes, or the
    folder's name for pyensembl)."""
    out = {}
    tool = raw_tool(subject)
    for v in A.raw_versions(tool, organism):
        d = os.path.join(C.family_dir("g"), "raw", tool, v, organism)
        rels = set()
        p = os.path.join(d, "passes.jsonl")
        if os.path.exists(p):
            rels = {json.loads(x)["release"] for x in open(p) if x.strip()}
        elif tool == "pyensembl":
            rels = {int(v)}
        else:
            rels = {None}  # a local tool: one snapshot
        for r in rels:
            out.setdefault(r, []).append(v)
    return out


def posable(subject, line, organism, gprofiler_archives=None):
    if subject in T.TOOLS or subject == "ensembl_rest_lookup":
        tool = T.TOOLS[raw_tool(subject)]
        if subject == "gprofiler" and gprofiler_archives is not None:
            tool.instances = gprofiler_archives
        if subject == "ensembl_rest_lookup" and line.get("task") not in ("c", "d"):
            return None, "lookup/symbol looks up symbols only"
        return tool.posable(line, organism)
    if subject == "gget":
        return LO.gget_posable(line, organism)
    if subject in LO.R_TOOLS:
        return LO.r_posable(subject, line, organism)
    if subject == "pyensembl":
        return LO.posable(line, organism)
    if subject == "ncbi_files":
        return LO.ncbi_posable(line, organism)
    raise ValueError(subject)


def tier_of(subject):
    if subject in LO.R_TOOLS:
        return LO.R_TOOLS[subject][3]
    return "local" if subject in LOCAL else "web"


def answer_rows(subject, organism, form):
    """(rows of the runner's answer schema, the versions read, the asking dates) of one subject, form and organism."""
    lines = F.g_lines(organism, tier_of(subject))
    vb = versions_by_release(subject, organism)
    archives = None
    if subject == "gprofiler":
        archives = {int(v[1: v.index("_")]): v for v in A.raw_versions("gprofiler", organism)
                    if v.startswith("e") and v.count("_") == 2}  # the archives (the current one has a fourth part)
    by_rel, rows = {}, []
    for q in lines:
        rel, why = posable(subject, q, organism, archives)
        if rel is None:
            for view in ("best", "all"):
                rows.append(dict(qid=q["qid"], view=view, outcome="not_supported", got=[], error=why, seconds=None, calls=0))
        else:
            by_rel.setdefault(rel, []).append(q)
    used, dates = set(), set()
    for rel, qs in by_rel.items():
        vs = vb.get(rel) or vb.get(None) or []
        got = {}
        for v in vs:
            got.update(A.answers(subject, organism, v, qs, release=rel if rel in (116, 114, 113, 111) else None, form=form))
            used.add(v)
            st = A.K.RawStore.of(raw_tool(subject), v, organism)
            dates |= {r["asked_at"][:10] for r in st.rows.values() if r.get("asked_at")}
        for q in qs:
            a = got.get(q["qid"])
            if a is None:
                continue
            for view in ("best", "all"):
                rows.append(dict(qid=q["qid"], view=view, outcome=a.outcome, got=list(a.got), status=a.status,
                                 error=a.error, seconds=a.seconds, calls=1, from_release=None,
                                 flags=json.dumps({"release": rel, "form": form})))
    return rows, sorted(used), sorted(dates)


def idtrack_rows(organism, tree):
    """The tree's answers to G's lines from the families' runs (read only)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    lines = F.g_lines(organism, "local")
    by_src = {}
    for q in lines:
        by_src.setdefault(q["source_family"], set()).add(q["qid"])
    tabs, found = [], []
    for src, qids in by_src.items():
        run = os.path.join(C.family_dir(src), "runs", tree, organism)
        files = C.answer_files(run)
        if not files:
            continue
        found.append(CFG.in_root(run))
        vs = pa.array(sorted(qids))
        for p in files:
            t = pq.read_table(p)
            tabs.append(t.filter(pc.is_in(t.column("qid"), value_set=vs)))
    return tabs, found


def score_folder(run, organism):
    from idtrack_reproducibility.families import scorer as S
    from idtrack_reproducibility.harness import score as HS
    a = types.SimpleNamespace(family="g", organism=organism, run=run, sets=None, within=None, strata=None, kinds=None,
                              limit=None, judge="revised", scores_to=run, parts=1, part=None, merge=False)
    for old in glob.glob(os.path.join(run, "outcomes", "part_*")):
        os.remove(old)
    S.run_part(a, 0, HS.Places(organism))
    S.merge(a)


def write_answers(run, rows):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from idtrack_reproducibility.families import runner as R
    d = os.path.join(run, "answers")
    os.makedirs(d, exist_ok=True)
    for old in glob.glob(os.path.join(d, "*.parquet")):
        os.remove(old)  # the answers are made again from the stored responses
    for r in rows:
        r["h"] = C.crc(r["qid"])
    path = os.path.join(d, "shard_0_of_1.0.parquet")
    pq.write_table(pa.Table.from_pylist(rows, schema=R.answer_schema()), path + ".tmp", compression="zstd")
    os.replace(path + ".tmp", path)


def score_subject(subject, organism, form, tree="main"):
    run = os.path.join(C.family_dir("g"), "runs", subject, form, organism)
    if subject == "idtrack":
        import pyarrow as pa
        import pyarrow.parquet as pq
        tabs, found = idtrack_rows(organism, tree)
        d = os.path.join(run, "answers")
        if tabs:
            os.makedirs(d, exist_ok=True)
            for old in glob.glob(os.path.join(d, "families_*.parquet")):
                os.remove(old)
            pq.write_table(pa.concat_tables(tabs, promote_options="default"),
                           os.path.join(d, "families_0_of_1.0.parquet"))
        elif not glob.glob(os.path.join(d, "load_*.json")):
            return None  # neither the families' runs of the tree nor the runner's answers on G's sets
        # the runner's own answers on G's sets (`runner --family g`, for lines no family run holds) stay beside them
        subj = dict(kind="tree", name=tree, version=tree, asked_at=None, endpoint=";".join(found) or "runner on G's sets")
    else:
        rows, used, dates = answer_rows(subject, organism, form)
        if not used:
            return None
        write_answers(run, rows)
        subj = dict(kind="tool", name=subject, version=",".join(used), asked_at=",".join(dates) if dates else None,
                    endpoint=ENDPOINT[subject])
    score_folder(run, organism)
    p = os.path.join(run, "run.json")
    info = json.load(open(p))
    if subject != "idtrack":
        info["subject"] = subj
    info["g"] = dict(form=form, tier=tier_of(subject) if subject != "idtrack" else "local", tree=tree if subject == "idtrack" else None,
                     sources=subj["endpoint"] if subject == "idtrack" else None)
    json.dump(info, open(p, "w"), indent=1)
    return run


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.score", description=__doc__.split("\n\n")[0])
    ap.add_argument("--subjects", default="all")
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--forms", default="best,stated")
    ap.add_argument("--tree", default="main", help="the idtrack tree whose families' answers are idtrack's rows")
    a = ap.parse_args(argv)
    subs = list(SUBJECTS) + ["idtrack"] if a.subjects == "all" else a.subjects.split(",")
    orgs = CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")
    for s in subs:
        for o in orgs:
            for form in (["stated"] if s == "idtrack" else a.forms.split(",")):
                r = score_subject(s, o, form, a.tree)
                print(s, o, form, CFG.in_root(r) if r else "nothing stored", flush=True)


if __name__ == "__main__":
    main()
