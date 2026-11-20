"""`idtrack-repro results`: the repository's results, generated from the scored runs of the data root.

The registry `results.toml` names, per benchmark, tree and organism, the scored folder in the data root and the tables
to take from it. For every table this command writes `results/<family>/<tree>/<organism>/<name>` (a copy of the table
the pipeline wrote into the scored folder), a row of `results/provenance.tsv` (the file's sha256, the source in the
data root with its sha256, when the source was scored, the idtrack commit and the graph's content digest), and, for
tables of rates, the rows of `results/index.tsv` (one per family, stratum, organism and tree: n, the outcomes, the
rate and its interval). Regenerating replaces: files of the registry's benchmarks that the registry no longer names
are removed, so `results/` holds only the current results.

Index columns: family, stratum, organism, tree, scope (`exhaustive`: every line of the family's sets; `panel`: the
family's panel, docs/methods.md "Measurement design"; a registry entry says `scope`, default `exhaustive`), commit,
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
SCOPES = ("exhaustive", "panel")  # every line of the family's sets, or its panel (docs/methods.md, "Measurement design")
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
