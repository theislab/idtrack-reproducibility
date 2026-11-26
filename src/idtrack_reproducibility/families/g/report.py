"""The tables of family G from its scored folders (docs/families/g.md, "Tables").

Usage: python -m idtrack_reproducibility.families.g.report [--organisms ...] [--out <folder>]
       python -m idtrack_reproducibility.families.g.report --publish   (after `report` and `reading`: PUBLISHED below)

Reads every scored folder `<families>/g/runs/<subject>/<form>/<organism>/` (the record's rows), G's lines (tier,
task, the source family's fields), the rule of equal truth (`<families>/g/likeforlike/<organism>/release_<T>.parquet`)
and the raw index (`<families>/g/raw/`), and writes into `--out` (default `<families>/g/report/`):

  tool_task.tsv        subject x form x organism x task: the tier's lines (without the lines their source family
                       withdrew, `excluded`, counted apart), `not_supported`, posed (answered),
                       `not_answered` (posable, no usable response stored: not answered by the service on the dates
                       asked), left out
                       by the rule of equal truth, kept, right / wrong / undetermined / error under `all` and `best`,
                       the rate with its Wilson interval, coverage (kept lines with a gene), and idtrack's rate on the
                       same kept lines (the tree's answers, `all` and `best`)
  not_supported.tsv    subject x organism x reason: the lines a subject cannot pose, counted by reason
  by_distance.tsv      A1's and A*'s lines by distance in releases to 115: n, right, rate per subject (form best)
  forms.tsv            per subject and organism: the posed lines whose two readings differ, and how many of them
                       change their outcome between the forms
  speed.tsv            per subject and organism: requests, seconds of the requests, lines asked, seconds per 1,000 lines
  failures.tsv         per organism and task: the sets of subjects wrong (form best, view all) on the lines kept for
                       every subject of the set, counted (an UpSet's table)
  disagreements.tsv    every line where a subject and idtrack disagree (form best, view all), its class, and a seeded
                       sample of at most 200 per subject, task and organism stratified by the outcome pair (`sampled`)
  disagreement_classes.tsv  the classes counted: `idtrack wrong` (the subject right, idtrack wrong), `the tool wrong`
                       (idtrack right), `both wrong`, `the truth cannot decide` (either undetermined)
"""
import argparse
import collections
import glob
import json
import os
import random

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import report as R
from idtrack_reproducibility.families.g import fetch as F

COUNTED = ("right", "wrong", "error")
DIST_BINS = ((0, 0), (1, 5), (6, 15), (16, 30), (31, 10_000))


def dist_bin(d):
    for lo, hi in DIST_BINS:
        if lo <= d <= hi:
            return f"{lo}-{hi}" if hi < 10_000 else f"{lo}+"
    return "?"


def scored(organism):
    """{(subject, form): scored folder} of an organism."""
    out = {}
    for p in glob.glob(os.path.join(C.family_dir("g"), "runs", "*", "*", organism, "run.json")):
        d = os.path.dirname(p)
        subject, form = d.split(os.sep)[-3], d.split(os.sep)[-2]
        out[(subject, form)] = d
    return out


def frame(folder):
    import pandas as pd
    import pyarrow.parquet as pq
    o = R.read_outcomes(folder, columns=["qid", "view", "outcome", "got", "reason"])
    a = pd.concat([pq.read_table(p, columns=["qid", "view", "flags"]).to_pandas() for p in C.answer_files(folder)])
    rel = {q: json.loads(f).get("release") for q, f in zip(a["qid"], a["flags"]) if isinstance(f, str) and f}
    o["release"] = o["qid"].map(rel)
    return o


def likeforlike(organism):
    import pandas as pd
    out = {}
    for p in glob.glob(os.path.join(C.family_dir("g"), "likeforlike", organism, "release_*.parquet")):
        t = pd.read_parquet(p)
        out[int(p.rsplit("_", 1)[1].split(".")[0])] = dict(zip(t["qid"], zip(t["equal"], t["reason"])))
    return out


def kept_mask(o, lines, lf):
    """(kept, left out reason) per row: a line asked at release 115 of a tool answering at T != 115 is kept when the
    rule of equal truth keeps it at T; every other posed line is kept."""
    kept, why = [], []
    for q, rel in zip(o["qid"], o["release"]):
        ln = lines.get(q)
        if ln is None or rel is None or rel != rel or ln.get("to_release") != 115 or int(rel) == 115:
            kept.append(True)
            why.append("")
            continue
        eq, r = lf.get(int(rel), {}).get(q, (False, f"no table of release {int(rel)}"))
        kept.append(bool(eq))
        why.append("" if eq else r.split(":")[0])
    return kept, why


def run(organisms, out):
    import pandas as pd
    os.makedirs(out, exist_ok=True)
    tt, ns, dist, forms, speed, fail, dis = [], [], [], [], [], [], []
    for org in organisms:
        sc = scored(org)
        if not sc:
            continue
        lines = {q["qid"]: q for q in F.g_lines(org, "local")}
        lf = likeforlike(org)
        idt = None
        if ("idtrack", "stated") in sc:
            idt = frame(sc[("idtrack", "stated")])
            idt = {(q, v): o for q, v, o in zip(idt["qid"], idt["view"], idt["outcome"])}
        per_subject = {}
        for (subj, form), folder in sorted(sc.items()):
            if subj == "idtrack":
                continue
            info = json.load(open(os.path.join(folder, "run.json")))
            tier = (info.get("g") or {}).get("tier", "web")
            o = frame(folder)
            o["task"] = o["qid"].map(lambda q: lines.get(q, {}).get("task"))
            o["kept"], o["left_out"] = kept_mask(o, lines, lf)
            if form == "best":
                per_subject[subj] = o
            tier_all = [q for q in lines.values() if tier == "local" or q.get("tier") == "web"]
            tier_lines = [q for q in tier_all if not q.get("excluded")]  # withdrawn by the source family: never scored
            for task in sorted({q["task"] for q in tier_all}):
                n_lines = sum(1 for q in tier_lines if q["task"] == task)
                sub = o[(o.task == task)]
                row = dict(subject=subj, form=form, organism=org, task=task, tier=tier, lines=n_lines,
                           excluded=sum(1 for q in tier_all if q["task"] == task and q.get("excluded")))
                for view in ("all", "best"):
                    v = sub[sub.view == view]
                    row["not_supported"] = int((v.outcome == "not_supported").sum())
                    posed = v[v.outcome != "not_supported"]
                    row["posed"] = len(posed)
                    # lines the tool can pose that no usable stored response answers: not answered by the service on
                    # the dates asked (never scored, never dropped: the next pass asks them and they join the tables)
                    row["not_answered"] = n_lines - row["not_supported"] - len(posed)
                    row["left_out"] = int((~posed.kept).sum())
                    k = posed[posed.kept]
                    row["kept"] = len(k)
                    c = k.outcome.value_counts().to_dict()
                    for x in ("right", "wrong", "undetermined", "error"):
                        row[f"{x}_{view}"] = int(c.get(x, 0))
                    den = sum(int(c.get(x, 0)) for x in COUNTED)
                    row[f"rate_{view}"] = round(c.get("right", 0) / den, 4) if den else None
                    lo, hi = R.wilson(int(c.get("right", 0)), den) if den else (None, None)
                    row[f"wilson_low_{view}"], row[f"wilson_high_{view}"] = lo, hi
                    if view == "all":
                        row["coverage"] = round(float((k.got.fillna("") != "").mean()), 4) if len(k) else None
                    if idt is not None:
                        io = [idt.get((q, view)) for q in k.qid]
                        ic = collections.Counter(x for x in io if x)
                        iden = sum(ic[x] for x in COUNTED)
                        row[f"idtrack_n_{view}"] = sum(ic.values())
                        row[f"idtrack_rate_{view}"] = round(ic["right"] / iden, 4) if iden else None
                tt.append(row)
            v = o[o.view == "all"]
            for reason, n in v[v.outcome == "not_supported"].reason.value_counts().items():
                ns.append(dict(subject=subj, form=form, organism=org, reason=reason, lines=int(n)))
            # speed: the raw index of the subject's stored responses
        for subj in per_subject:
            rows = []
            for f in glob.glob(os.path.join(C.family_dir("g"), "raw", "ensembl_rest" if subj == "ensembl_rest_lookup" else subj,
                                            "*", org, "index.jsonl")):
                if "/_versions/" in f:
                    continue
                rows += [r for r in (json.loads(x) for x in open(f) if x.strip()) if "@" in r["endpoint"]]  # the passes, not the pilot
            if subj.startswith("ensembl_rest"):
                rows = [r for r in rows if r["endpoint"].startswith("lookup/symbol") == (subj == "ensembl_rest_lookup")]
            secs = sum(r.get("seconds") or 0 for r in rows)
            nl = len({q for r in rows for q in r.get("qids") or []})
            speed.append(dict(subject=subj, organism=org, requests=len(rows), seconds=round(secs, 1), lines=nl,
                              seconds_per_1000_lines=round(1000 * secs / nl, 2) if nl else None))
        # forms
        for (subj, form), folder in sc.items():
            if form != "stated" or (subj, "best") not in sc:
                continue
            from idtrack_reproducibility.families.g import tools as T
            b, s = per_subject.get(subj), frame(folder)
            changed = {q["qid"] for q in T.in_form(list(lines.values()), "best")}
            bm = {(q, v): x for q, v, x in zip(b.qid, b.view, b.outcome)}
            sm = s[(s.view == "all") & (s.outcome != "not_supported")]
            diff = [q for q in sm.qid if q in changed]
            moved = sum(1 for q, x in zip(sm.qid, sm.outcome) if q in changed and bm.get((q, "all")) != x)
            forms.append(dict(subject=subj, organism=org, posed=len(sm), readings_differ=len(diff), outcome_moved=moved))
        # by distance (A1 and A*: the source release's distance to 115)
        for subj, o in list(per_subject.items()) + ([("idtrack", None)] if idt is not None else []):
            src = o if o is not None else frame(sc[("idtrack", "stated")]).assign(kept=True)
            v = src[(src.view == "all") & (src.outcome != "not_supported")]
            if "kept" in v:
                v = v[v.kept]
            c = collections.defaultdict(collections.Counter)
            for q, x in zip(v.qid, v.outcome):
                ln = lines.get(q) or {}
                if ln.get("source_family") in ("a", "astar") and ln.get("from_release") is not None:
                    c[dist_bin(115 - int(ln["from_release"]))][x] += 1
            for b_, cc in sorted(c.items()):
                den = sum(cc[x] for x in COUNTED)
                dist.append(dict(subject=subj, organism=org, distance=b_, n=sum(cc.values()), right=cc["right"],
                                 rate=round(cc["right"] / den, 4) if den else None))
        # failures' overlap and disagreements (form best, view all, kept lines)
        wrong_by = collections.defaultdict(set)
        kept_by = collections.defaultdict(set)
        for subj, o in per_subject.items():
            v = o[(o.view == "all") & (o.outcome != "not_supported") & o.kept]
            for q, x, t in zip(v.qid, v.outcome, v.task):
                kept_by[(t, q)].add(subj)
                if x in ("wrong", "error"):
                    wrong_by[(t, q)].add(subj)
            if idt is None:
                continue
            pairs = collections.defaultdict(list)
            for q, x, t, got in zip(v.qid, v.outcome, v.task, v.got):
                y = idt.get((q, "all"))
                if y is None or x == y:
                    continue
                cls = ("idtrack wrong" if x == "right" and y in ("wrong", "error") else
                       "the tool wrong" if y == "right" and x in ("wrong", "error") else
                       "the truth cannot decide" if "undetermined" in (x, y) else "both wrong")
                pairs[(t, x, y)].append(dict(subject=subj, organism=org, task=t, qid=q, tool=x, idtrack=y, cls=cls,
                                             got=got, expected=",".join((lines.get(q) or {}).get("expect") or [])))
            by_task = collections.defaultdict(list)
            for (t, x, y), rs in pairs.items():
                by_task[t].append(rs)
            for t, groups in by_task.items():
                rng = random.Random(f"{C.SEED}:{org}:{subj}:{t}")
                quota, left = {}, 200
                groups = sorted(groups, key=len)
                for i, g in enumerate(groups):
                    quota[id(g)] = min(len(g), left // (len(groups) - i))
                    left -= quota[id(g)]
                for g in groups:
                    pick = {r["qid"] for r in rng.sample(sorted(g, key=lambda r: r["qid"]), quota[id(g)])}
                    for r in g:
                        dis.append(dict(r, sampled=r["qid"] in pick))
        combos = collections.Counter()
        for (t, q), subs in wrong_by.items():
            combos[(t, ",".join(sorted(subs)), len(kept_by[(t, q)]))] += 1
        for (t, subs, nk), n in sorted(combos.items()):
            fail.append(dict(organism=org, task=t, wrong=subs, subjects_kept=nk, lines=n))
    tabs = dict(tool_task=tt, not_supported=ns, by_distance=dist, forms=forms, speed=speed, failures=fail,
                disagreements=dis)
    for name, rows in tabs.items():
        pd.DataFrame(rows).to_csv(os.path.join(out, f"{name}.tsv"), sep="\t", index=False)
    d = pd.DataFrame(dis)
    if len(d):
        d.groupby(["subject", "organism", "task", "cls"]).agg(lines=("qid", "size"), sampled=("sampled", "sum")) \
            .reset_index().to_csv(os.path.join(out, "disagreement_classes.tsv"), sep="\t", index=False)
    return {k: len(v) for k, v in tabs.items()}


# The tables published in the repository (results/g/main/all/: the tables a manuscript cites, one row per subject,
# task, class or organism) and the per-line tables that stay in the data root (`<families>/g/report/`, listed in the
# family's data manifest): every disagreement (`disagreements.tsv`, about 1 million lines) and every sampled line read
# (`reading.tsv`), as the other families keep their per-line tables (wrong.tsv.gz, causes.tsv.gz) in their runs.
PUBLISHED = ("tool_task", "not_supported", "by_distance", "forms", "speed", "failures", "disagreement_classes",
             "reading_classes", "reading_summary")
PER_LINE = ("disagreements", "reading")
PUBLISH_DIR = os.path.join("results", "g", "main", "all")
GENERATOR = "python -m idtrack_reproducibility.families.g.report --publish"


def publish(root):
    """The published tables (PUBLISHED) copied from `<families>/g/report/` into `results/g/main/all/`, each with its
    `.provenance` (sources: the table in the data root, every registered scored folder's run.json, the store's summary),
    and the document's tables written from them (`write_doc`)."""
    import datetime
    import shutil
    from idtrack_reproducibility import results as RES
    data = CFG.get().root
    rep = os.path.join(C.family_dir("g"), "report")
    out = os.path.join(root, PUBLISH_DIR)
    os.makedirs(out, exist_ok=True)
    runs = [os.path.join(f, "run.json") for e, _, f in RES.registered_runs(root, data) if e["family"] == "g"]
    runs = [p for p in runs if os.path.exists(p)]
    store = os.path.join(C.family_dir("g"), "raw", "STORE_SUMMARY.tsv")
    now = datetime.datetime.now().isoformat(timespec="seconds")
    common = [(os.path.relpath(p, data), RES.sha256(p)) for p in runs + [store]]
    for name in PUBLISHED:
        src = os.path.join(rep, f"{name}.tsv")
        dst = os.path.join(out, f"{name}.tsv")
        shutil.copyfile(src, dst + ".tmp")
        os.replace(dst + ".tmp", dst)
        rel = os.path.relpath(dst, root)
        prov = [dict(file=rel, sha256="", generated_by=GENERATOR, source=s, source_sha256=h, generated_at=now)
                for s, h in [(os.path.relpath(src, data), RES.sha256(src))] + common]
        RES.write_tsv(dst + ".provenance", RES.PROV_COLS, prov)
    write_doc(root)
    return len(PUBLISHED)


ORGANISMS = ("sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis")
ORGANISM_NAMES = {"sus_scrofa": "pig", "mus_musculus": "mouse", "homo_sapiens": "human", "macaca_mulatta": "rhesus",
                  "macaca_fascicularis": "cynomolgus"}
SUBJECT_NAMES = {"ensembl_rest": "Ensembl REST", "ensembl_rest_lookup": "Ensembl REST `lookup/symbol`",
                 "biomart": "BioMart", "gprofiler": "g:Profiler", "mygene": "MyGene.info", "hgnc": "HGNC",
                 "uniprot": "UniProt", "bridgedb": "BridgeDb", "pyensembl": "pyensembl", "orgdb": "org.\\*.eg.db",
                 "ensdb": "EnsDb 113", "seurat": "Seurat", "gget": "gget", "limma": "limma", "hgnchelper": "HGNChelper",
                 "ncbi_files": "NCBI's files"}
GROUPS = ("truly wrong", "right under another reading", "not a question the tool claims", "the service failed")


def _pct(x):
    return "–" if x in (None, "", "nan") or x != x else f"{100 * float(x):.0f}"


def doc_tables(root):
    """{block: Markdown table} of the document's tables, from the published tables (`results/g/main/all/`): `rates`
    (each tool's rate and idtrack's on the same kept lines, form best, view all, with the kept lines) and `reading`
    (the reading's groups per tool and task, weighted shares)."""
    import pandas as pd
    d = os.path.join(root, PUBLISH_DIR)
    tt = pd.read_csv(os.path.join(d, "tool_task.tsv"), sep="\t")
    tt = tt[tt.form == "best"]
    lines = ["| tool | task | " + " | ".join(ORGANISM_NAMES[o] for o in ORGANISMS) + " |",
             "| --- | --- | " + " | ".join("---" for _ in ORGANISMS) + " |"]
    for (s, t), g in tt.groupby(["subject", "task"]):
        if g.rate_all.isna().all():
            continue
        cells = []
        for o in ORGANISMS:
            r = g[g.organism == o]
            if not len(r) or r.rate_all.isna().iloc[0] or r.kept.iloc[0] == 0:
                cells.append("–")
                continue
            r = r.iloc[0]
            cells.append(f"{_pct(r.rate_all)} / {_pct(r.idtrack_rate_all)} ({int(r.kept):,})")
        lines.append(f"| {SUBJECT_NAMES.get(s, s)} | {t} | " + " | ".join(cells) + " |")
    out = {"rates": "\n".join(lines)}
    rs = pd.read_csv(os.path.join(d, "reading_summary.tsv"), sep="\t")
    lines = ["| tool | task | tool wrong, read | " + " | ".join(GROUPS) + " | idtrack wrong, read | of them a known limitation |",
             "| --- | --- | --- | " + " | ".join("---" for _ in GROUPS) + " | --- | --- |"]
    for _, r in rs.iterrows():
        lines.append(f"| {SUBJECT_NAMES.get(r.subject, r.subject)} | {r.task} | {int(r.tool_wrong_sampled):,} | "
                     + " | ".join(_pct(r[g]) for g in GROUPS)
                     + f" | {int(r.idtrack_wrong_sampled):,} | {_pct(r.idtrack_known_limitation)} |")
    out["reading"] = "\n".join(lines)
    return out


DOC = os.path.join("docs", "families", "g.md")


def doc_block(text, name):
    """The text between the document's markers of a table (`<!-- table g:<name> -->` ... `<!-- /table -->`)."""
    start = f"<!-- table g:{name} -->\n"
    i = text.index(start) + len(start)
    return i, text.index("\n<!-- /table -->", i)


def write_doc(root):
    """The document's tables (between their markers) written from the published tables."""
    p = os.path.join(root, DOC)
    text = open(p).read()
    for name, table in doc_tables(root).items():
        i, j = doc_block(text, name)
        text = text[:i] + table + text[j:]
    open(p + ".tmp", "w").write(text)
    os.replace(p + ".tmp", p)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.report", description=__doc__.split("\n\n")[0])
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--out", default=None)
    ap.add_argument("--publish", action="store_true", help="publish the tables in results/g/ and write the document's "
                    "tables from them (after `report` and `reading`)")
    a = ap.parse_args(argv)
    if a.publish:
        from idtrack_reproducibility.check import ROOT
        print(publish(ROOT), "tables published")
        return
    orgs = CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")
    print(run(orgs, a.out or os.path.join(C.family_dir("g"), "report")))


if __name__ == "__main__":
    main()
