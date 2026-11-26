"""`idtrack-repro results`: the repository's results, generated from the scored runs of the data root.

The registry `results.toml` names, per benchmark, tree and organism, the scored folder in the data root and the tables
to take from it. For every table this command writes `results/<family>/<tree>/<organism>/<name>` (a copy of the table
the pipeline wrote into the scored folder), a row of `results/provenance.tsv` (the file's sha256, the source in the
data root with its sha256, when the source was scored, the idtrack commit and the graph's content digest), and, for
tables of rates, the rows of `results/index.tsv` (one per family, stratum, organism and tree: n, the outcomes, the
rate and its interval). Regenerating replaces: files of the registry's benchmarks that the registry no longer names
are removed, so `results/` holds only the current results.

Index columns: family, stratum, organism, tree, scope (`exhaustive`: every line of the family's sets; `panel`: the
family's panel, docs/methods.md "Measurement design"; `sets`: every line of the sets the entry names in `sets`, the
others not asked on that tree; `modest`: the family's modest selection, up to 500 sampled lines per stratum,
docs/modest.md; a registry entry says `scope`, default `exhaustive`), commit,
graph_digest, n, right, wrong, undetermined, timeout, error, rate, low, high, interval (`cluster` for the cluster
bootstrap over the stratum's units, `wilson` where the table has no clustered interval), file.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import json
import os
import shutil
import sys
import tomllib

from . import config as CFG

INDEX_COLS = ["family", "stratum", "organism", "tree", "scope", "commit", "graph_digest", "n", "right", "wrong",
              "undetermined", "timeout", "error", "rate", "low", "high", "interval", "file"]
# every line of the family's sets, its panel (docs/methods.md, "Measurement design"), or every line of the sets the
# entry names (`sets`; the family's other sets not asked on that tree)
SCOPES = ("exhaustive", "panel", "sets", "modest")
PROV_COLS = ["file", "sha256", "generated_by", "source", "source_sha256", "scored_at", "commit", "graph_digest",
             "generated_at"]
OWN = {"index.tsv", "provenance.tsv"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_tsv(path, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    os.replace(path + ".tmp", path)


def registry(root):
    with open(os.path.join(root, "results.toml"), "rb") as fh:
        return tomllib.load(fh)


def scored_info(folder, rel_jsons):
    """When the source was scored (the latest of its run.json files) and the graph digests they record."""
    at, digests = [], set()
    for rj in rel_jsons:
        p = os.path.join(folder, rj)
        if os.path.exists(p):
            d = json.load(open(p))
            if d.get("scored_at"):
                at.append(d["scored_at"])
            digests |= {x for x in (d.get("graph_digests") or []) if x and x != "None"}
    return (max(at) if at else ""), digests


def index_rows(kind, path, fam, org, tree, commit, digest, rel, scope="exhaustive"):
    if scope not in SCOPES:
        raise SystemExit(f"{rel}: scope must be one of {SCOPES}, not {scope!r}")
    rows = []
    t = read_tsv(path)
    if kind == "report":  # the common report.tsv: one row per stratum and view; `*` = both views
        for r in t:
            if r.get("view") != "*":
                continue
            boot = r.get("boot_low") not in (None, "", "nan")
            rows.append(dict(family=fam, stratum=r["stratum"], organism=org, tree=tree, scope=scope, commit=commit,
                             graph_digest=digest, n=r["n"], right=int(r["right"]) + int(r.get("expected_error") or 0),
                             wrong=r["wrong"], undetermined=r["undetermined"], timeout=r["timeout"], error=r["error"],
                             rate=r["rate"], low=r["boot_low"] if boot else r["wilson_low"],
                             high=r["boot_high"] if boot else r["wilson_high"], interval="cluster" if boot else "wilson",
                             file=rel))
    elif kind == "scores":  # the harness's scores.tsv: one row per kind and variant; `*` = the whole kind
        for r in t:
            if r.get("variant") != "*":
                continue
            rows.append(dict(family=fam, stratum=r["kind"], organism=org, tree=tree, scope=scope, commit=commit, graph_digest=digest,
                             n=r["n"], right=r["right"], wrong=r["wrong"], undetermined=r["undetermined"],
                             timeout=r["timeout"], error=r["error"], rate=r["rate"], low=r["ci_low"], high=r["ci_high"],
                             interval="wilson", file=rel))
    else:
        raise SystemExit(f"unknown index kind {kind!r}")
    return rows


def collect(root, only=None, dry=False):
    reg = registry(root)
    cfg = CFG.get()
    data = cfg.root
    res = os.path.join(root, "results")
    now = datetime.datetime.now().isoformat(timespec="seconds")
    index, prov, written, problems = [], [], set(), []
    graphs = reg.get("graphs", {})
    for e in reg.get("run", []):
        fam, tree = e["family"], e["tree"]
        if only and fam not in only:
            continue
        for org in e["organisms"]:
            folder = os.path.join(data, e["folder"].format(organism=org))
            at, digests = scored_info(folder, e.get("scored", ["run.json"]))
            digest = graphs.get(e.get("graphs", tree), {}).get(org, "")
            if digests and digest and digest not in digests:
                problems.append(f"{fam} {tree} {org}: the run's graph digest {sorted(digests)} is not the registry's {digest}")
            for name, src in e["files"].items():
                sp = os.path.join(folder, src)
                rel = os.path.join("results", fam, tree, org, name)
                if not os.path.exists(sp):
                    if name not in e.get("optional", []):  # a table the step writes for some organisms only
                        problems.append(f"{rel}: no source {sp}")
                    continue
                dst = os.path.join(root, rel)
                if not dry:
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copyfile(sp, dst + ".tmp")
                    os.replace(dst + ".tmp", dst)
                written.add(rel)
                src_rel = os.path.relpath(sp, data)
                prov.append(dict(file=rel, sha256=sha256(sp), generated_by=e.get("generated_by", {}).get(name, e.get("generator", "")),
                                 source=src_rel, source_sha256=sha256(sp), scored_at=at, commit=e["commit"],
                                 graph_digest=digest, generated_at=now))
                kind = e.get("index", {}).get(name)
                if kind:
                    index += index_rows(kind, sp, fam, org, tree, e["commit"], digest, rel, e.get("scope", "exhaustive"))
    for e in reg.get("generated", []):  # tables a command of the package writes straight into results/
        if only and e["family"] not in only:
            continue
        rels = []
        for rel in e["files"]:  # a pattern names every file the command writes (`*`)
            rels += sorted(os.path.relpath(x, root) for x in glob.glob(os.path.join(root, rel))) if "*" in rel else [rel]
        for rel in rels:
            p = os.path.join(root, rel)
            if not os.path.exists(p):
                problems.append(f"{rel}: not generated yet ({e['generated_by']})")
                continue
            written.add(rel)
            rows = read_tsv(p + ".provenance") if os.path.exists(p + ".provenance") else []
            for r in rows:
                r["file"], r["sha256"] = rel, sha256(p)
            prov += rows or [dict(file=rel, sha256=sha256(p), generated_by=e["generated_by"], generated_at=now)]
            kind = e.get("index")
            if kind:
                index += index_rows(kind, p, e["family"], e.get("organism", ""), e.get("tree", ""), e.get("commit", ""),
                                    "", rel, e.get("scope", "exhaustive"))
    if only:  # keep the rows of the other benchmarks
        keep = lambda rows: [r for r in rows if r["file"].split("/")[1] not in only]  # noqa: E731
        index = keep(read_tsv(os.path.join(res, "index.tsv"))) + index if os.path.exists(os.path.join(res, "index.tsv")) else index
        prov = keep(read_tsv(os.path.join(res, "provenance.tsv"))) + prov if os.path.exists(os.path.join(res, "provenance.tsv")) else prov
    if not dry:
        for p in glob.glob(os.path.join(res, "**", "*"), recursive=True):
            rel = os.path.relpath(p, root)
            fam = rel.split("/")[1] if rel.count("/") >= 2 else None
            if os.path.isfile(p) and rel not in written and os.path.basename(p) not in OWN and not p.endswith(".provenance") \
                    and (not only or fam in only) and fam in {e["family"] for e in reg.get("run", []) + reg.get("generated", [])}:
                os.remove(p)  # a result the registry no longer names: regenerating replaces
        index.sort(key=lambda r: (r["family"], r["tree"], r["organism"], r["stratum"], r["file"]))
        prov.sort(key=lambda r: r["file"])
        write_tsv(os.path.join(res, "index.tsv"), INDEX_COLS, index)
        write_tsv(os.path.join(res, "provenance.tsv"), PROV_COLS, prov)
    return index, prov, problems


def rows_behind(root, data, problems):
    """Every registered run's tables rest on rows in the data root: each scored folder holds its rows of the record
    (docs/records.md: `outcomes/*.parquet`, as many as its run.json counts, of the record's version), unless the
    registry says what it keeps instead (`rows = "<what>"`)."""
    import pyarrow.parquet as pq
    from .records import RECORD_VERSION
    for e in registry(root).get("run", []):
        if e.get("rows"):
            continue
        for org in e["organisms"]:
            base = os.path.join(data, e["folder"].format(organism=org))
            for rj in e.get("scored") or ["run.json"]:
                folder = os.path.dirname(os.path.join(base, rj))
                rel = os.path.relpath(folder, data)
                parts = glob.glob(os.path.join(folder, "outcomes", "*.parquet"))
                info = json.load(open(os.path.join(folder, "run.json"))) if os.path.exists(os.path.join(folder, "run.json")) else {}
                if not parts:
                    problems.append(("results", rel, "no rows behind its tables (no outcomes in the data root)"))
                    continue
                n = sum(pq.ParquetFile(p).metadata.num_rows for p in parts)
                want = info.get("n_rows")
                if not n or (want is not None and n != want):
                    problems.append(("results", rel, f"{n} outcome rows where run.json counts {want}"))
                if int(info.get("record_version", 0)) != RECORD_VERSION:
                    problems.append(("results", rel, f"rows of record version {info.get('record_version', 0)}, not "
                                     f"{RECORD_VERSION}: score it again (idtrack-repro rescore)"))


def registered_runs(root, data):
    """[(entry, organism, scored folder)] of every run the registry names (each scored folder of an entry's `scored`)."""
    out = []
    for e in registry(root).get("run", []):
        for org in e["organisms"]:
            base = os.path.join(data, e["folder"].format(organism=org))
            out += [(e, org, os.path.dirname(os.path.join(base, rj))) for rj in e.get("scored") or ["run.json"]]
    return out


def subjects_behind(root, data, problems):
    """A registered run's scored folders are of the entry's subject and reading: run.json names the entry's `subject`
    (an entry without one is a tree's own run, whose run.json names no subject; docs/records.md, "Another tool") and,
    where the entry names a `form`, the family's block of run.json records that form."""
    for e, org, folder in registered_runs(root, data):
        p = os.path.join(folder, "run.json")
        if not os.path.exists(p):
            continue  # reported by rows_behind
        info = json.load(open(p))
        rel = os.path.relpath(folder, data)
        named = (info.get("subject") or {}).get("name")
        if named != e.get("subject"):
            problems.append(("results", rel, f"the folder's subject is {named or 'the tree'}, the registry's "
                             f"{e.get('subject') or 'the tree'}"))
        if e.get("form") and (info.get(e["family"]) or {}).get("form") != e["form"]:
            problems.append(("results", rel, f"the folder's reading of the input is "
                             f"{(info.get(e['family']) or {}).get('form')}, the registry's {e['form']}"))


def judged_today(info):
    """The scoring key the common scorer gives a scored folder's selection today (`scorer.scoring_key`: its sets with
    the revisions published since, the record version and the code of the judge), or None where run.json records no
    key (a scoring made before keys were recorded)."""
    import types
    from .families import scorer as SC
    if not info.get("key"):
        return None
    a = types.SimpleNamespace(family=info["family"], organism=info["organism"], judge=info.get("judge"),
                              sets=list((info.get("sets") or {}).keys()), within=info.get("within"),
                              strata=info.get("strata"), kinds=info.get("kinds"), limit=info.get("limit"),
                              parts=info.get("parts"))
    return SC.scoring_key(a)


def judges_behind(root, data, problems):
    """A family whose registered runs compare subjects (entries with `subject`) is a comparison under one judge:
    every scored folder of the family holds the scoring key the scorer gives it today, so that every subject is
    judged by today's judge on today's truth (a family of one subject records the judge it was scored with, and
    `rescore` brings it to today's)."""
    runs = registered_runs(root, data)
    compared = {e["family"] for e, _, _ in runs if e.get("subject")}
    for e, org, folder in runs:
        p = os.path.join(folder, "run.json")
        if e["family"] not in compared or not os.path.exists(p):
            continue
        info = json.load(open(p))
        now = judged_today(info)
        if now != info.get("key"):
            problems.append(("results", os.path.relpath(folder, data), "scored under another judge or truth than "
                             "today's (its scoring key is not the scorer's today): score it again"))


def tables_behind(root, data, problems):
    """A family's generated tables are made from its registered runs: each generated file of a family with runs names
    every registered scored folder (its run.json) among the sources of its provenance."""
    reg = registry(root)
    runs = registered_runs(root, data)
    for e in reg.get("generated", []):
        want = {os.path.relpath(os.path.join(f, "run.json"), data) for g, _, f in runs if g["family"] == e["family"]}
        if not want:
            continue
        for rel in e["files"]:
            for p in sorted(glob.glob(os.path.join(root, rel))) if "*" in rel else [os.path.join(root, rel)]:
                if not os.path.exists(p + ".provenance"):
                    problems.append(("results", os.path.relpath(p, root), "no provenance of its sources"))
                    continue
                have = {r["source"] for r in read_tsv(p + ".provenance")}
                for w in sorted(want - have)[:3]:
                    problems.append(("results", os.path.relpath(p, root), f"not made from the registered run {w}"))


def stores_behind(root, data, problems):
    """A family's store of fetched responses (a `store` entry) matches its manifest: every file the store's manifest
    lists (`manifest`, written by the family's store command) is in the family's data manifest with the same size and
    sha256 (the data manifest is checked against the data root by the manifests' check), or, where the data manifest's
    area leaves such a file out, is in the data root with its size."""
    import re
    from . import manifests as M
    stores = registry(root).get("store", [])
    areas = {a["name"]: a for a in M.areas(root)} if stores else {}
    for s in stores:
        mp = os.path.join(data, s["manifest"])
        name = f"families/{s['family']}"
        dm = [p for p in (os.path.join(root, "manifests", name + x) for x in (".tsv.gz", ".tsv")) if os.path.exists(p)]
        if not os.path.exists(mp) or not dm:
            problems.append(("results", s["manifest"], "the store's manifest or the family's data manifest is missing"))
            continue
        have = {r["path"]: r for r in M.read_manifest(dm[0])}
        exclude = [re.compile(x) for x in (areas.get(name) or {}).get("exclude", [])]
        for r in read_tsv(mp):
            sub = os.path.join(r["tool"], r["version"], r["organism"], r["file"]) if r["area"] == "store" else r["file"]
            path = os.path.join(s["folder"], sub)
            m = have.get(path)
            if m is None:
                if any(x.search(path) for x in exclude) and os.path.exists(os.path.join(data, path)) and \
                        os.path.getsize(os.path.join(data, path)) == int(r["bytes"]):
                    continue
                problems.append(("results", path, "in the store's manifest, not in the family's data manifest"))
            elif m["bytes"] != r["bytes"] or m["sha256"] != r["sha256"]:
                problems.append(("results", path, "the store's manifest and the data manifest differ in size or sha256"))


def check(root, problems, skipped):
    """The results' checks of `idtrack-repro check` (see check.py)."""
    res = os.path.join(root, "results")
    ip, pp = os.path.join(res, "index.tsv"), os.path.join(res, "provenance.tsv")
    if not os.path.exists(ip) or not os.path.exists(pp):
        problems.append(("results", "results/", "no index.tsv or provenance.tsv: run `idtrack-repro results`"))
        return
    prov = {r["file"]: r for r in read_tsv(pp)}
    for r in read_tsv(ip):
        if not os.path.exists(os.path.join(root, r["file"])):
            problems.append(("results", "results/index.tsv", f"row of a missing file: {r['file']}"))
        if r["file"] not in prov:
            problems.append(("results", "results/index.tsv", f"{r['file']} has no provenance row"))
    for p in glob.glob(os.path.join(res, "**", "*"), recursive=True):
        rel = os.path.relpath(p, root)
        if not os.path.isfile(p) or os.path.basename(p) in OWN or p.endswith(".provenance"):
            continue
        if rel not in prov:
            problems.append(("results", rel, "not generated by the pipeline (no provenance row)"))
        elif prov[rel]["sha256"] != sha256(p):
            problems.append(("results", rel, "edited after it was generated (sha256 differs from its provenance)"))
    data = CFG.get().root
    if not os.path.isdir(data):
        skipped.append(f"results staleness: no data root at {data}")
        return
    rows_behind(root, data, problems)
    subjects_behind(root, data, problems)
    judges_behind(root, data, problems)
    tables_behind(root, data, problems)
    stores_behind(root, data, problems)
    for rel, r in prov.items():
        if not r.get("source"):
            continue
        sp = os.path.join(data, r["source"])
        if not os.path.exists(sp):
            problems.append(("results", rel, f"its source {r['source']} is not in the data root"))
            continue
        if r.get("source_sha256") and sha256(sp) != r["source_sha256"]:
            problems.append(("results", rel, f"stale: {r['source']} changed since the result was generated"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro results", description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", default=None, help="comma-separated benchmarks (families) to regenerate")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args(argv)
    from .check import ROOT
    index, prov, problems = collect(ROOT, set(a.only.split(",")) if a.only else None, a.dry)
    print(f"{len(prov)} result files, {len(index)} index rows")
    for p in problems:
        print("problem:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
