"""A pass of one tool over G's lines: its version calls, then every line it can pose, asked through the client and
stored as fetched (docs/families/g.md, "Fetching").

Usage: python -m idtrack_reproducibility.families.g.fetch --tool <tool> --organisms <o,...|all> [--tier web|local]
                                                         [--per-family N] [--draft <folder>]

- The lines are G's published sets of the organism (`g_<family>_v1`; `--draft`: the drawn lines in a folder, before
  they are published), of the tier asked (`web`: the web services' lines; `local`: every line); `--per-family N`: the
  first N lines of each family in the set's order (a pilot).
- The version calls of the tool are asked first (salted with the pass, so every pass asks them again) into
  `<families>/g/raw/<tool>/_versions/<organism>/`; the version they report names the folder of the pass's
  responses, `<families>/g/raw/<tool>/<version>/<organism>/`, and the pass is recorded in its `passes.jsonl`.
- A line already held by a stored response of the same endpoint is not asked again (client.RawStore.answered), so a
  pass resumes where an earlier one stopped, whatever its batches were.
"""
import argparse
import gzip
import json
import os
import time
import urllib.parse

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.g import client as K
from idtrack_reproducibility.families.g import tools as T

G_FAMILIES = ("a1", "astar", "b", "c", "d1")
# Hosts that failed a request after its retries in this pass (an archive of BioMart that answers 500): the rest of
# the pass leaves their lines for the next pass instead of waiting out every request again.
FAILED_HOSTS = set()
GROUP = {"ensembl_rest": "ensembl_rest", "biomart": "ensembl_biomart", "gprofiler": "gprofiler", "mygene": "mygene",
         "hgnc": "hgnc", "uniprot": "uniprot", "bridgedb": "bridgedb"}


def g_lines(organism, tier="web", per_family=None, draft=None):
    """G's lines of one organism and tier, in each set's order (the first `per_family` of each family)."""
    out = []
    for gfam in G_FAMILIES:
        name = f"g_{gfam}_v1"
        if draft:
            p = os.path.join(draft, f"{organism}.{name}.jsonl.gz")
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt") as fh:
                xs = [json.loads(x) for x in fh]
        else:
            names = L.set_names("g", organism)
            if name not in names:
                continue
            xs = list(L.iter_set(L._set_file(L.query_dir("g", organism), name)))
            # a revision set of the family (g_<family>_r<n>) replaces its lines: the line of the latest set is the truth
            for r in [n for n in names if n.startswith(f"g_{gfam}_r")]:
                rev = {q["qid"]: q for q in L.iter_set(L._set_file(L.query_dir("g", organism), r))}
                xs = [rev.get(q["qid"], q) for q in xs]
        xs = [q for q in xs if tier == "local" or q.get("tier") == "web"]
        out += xs[:per_family] if per_family else xs
    return out


def version_of(tool, rows, store):
    """The version a tool's version calls report (the folder of its responses)."""
    body = {r["endpoint"]: store.body(r) for r in rows if r.get("file") and not r.get("error")}
    if tool.name == "ensembl_rest":
        return str(json.loads(body["info/data"])["releases"][0])
    if tool.name == "gprofiler":
        return json.loads(body["data_versions"])["gprofiler_version"]
    if tool.name == "mygene":
        return str(json.loads(body["metadata"])["build_version"])
    if tool.name == "hgnc":
        return json.loads(body["info"])["lastModified"][:10]
    if tool.name == "bridgedb":
        # the block of the Ensembl gene-product database among the service's data sources
        blocks, cur = [], {}
        for x in body["properties"].decode().splitlines():
            if "\t" not in x:
                continue
            k, v = x.split("\t", 1)
            if k == "DATASOURCENAME" and cur:
                blocks.append(cur)
                cur = {}
            cur[k] = v
        blocks.append(cur)
        ens = [b for b in blocks if b.get("DATASOURCENAME") == "Ensembl"]
        return f"{ens[0].get('DATASOURCEVERSION', '?')}_{ens[0].get('BUILDDATE', '?')}" if ens else "unknown"
    if tool.name == "uniprot":
        hdr = headers_of(store, [r for r in rows if r["endpoint"] == "idmapping/fields"][0]["request_id"])
        return {k.lower(): v for k, v in hdr.items()}.get("x-uniprot-release", "unknown")
    raise ValueError(tool.name)


def headers_of(store, rid):
    out = {}
    p = os.path.join(store.folder, "headers.jsonl")
    with open(p) as fh:
        for raw in fh:
            r = json.loads(raw)
            if r["request_id"] == rid:
                out = r["headers"]
    return out


def ask_versions(tool, organism, client, pass_id, calls):
    vs = K.RawStore.of(tool.name, "_versions", organism)
    rows = [client.fetch(vs, salt=pass_id, **c) for c in calls]
    vs.compile()
    return rows, vs


def record_pass(store, **kw):
    with open(os.path.join(store.folder, "passes.jsonl"), "a") as fh:
        fh.write(json.dumps(kw, sort_keys=True, default=str) + "\n")


def ask(store, client, reqs, log):
    done_by = {}
    n = failed = 0
    for r in reqs:
        ep = r.get("endpoint")
        if ep not in done_by:
            done_by[ep] = store.answered(ep)
        if r.get("qids") and set(r["qids"]) <= done_by[ep]:
            continue
        row = client.fetch(store, **r)
        n += 1
        if not row.get("error"):
            done_by[ep].update(row.get("qids") or [])
            failed = 0
        else:
            failed += 1
            host = urllib.parse.urlparse(r["url"]).hostname
            if failed >= 1:  # a server that failed a request after its retries: the rest waits for the next pass
                FAILED_HOSTS.add(host)
                log(f"{store.tool} {store.organism}: a request to {host} failed after its retries ({row['error']}); "
                    f"the host's remaining requests of this pass are left for the next pass")
                break
        if n % 200 == 0:
            store.compile()
            log(f"{store.tool} {store.organism}: {n} requests")
    store.compile()
    return n


def pass_tool(tool_name, organism, lines, client, log=print, releases=None):
    """One pass of a tool over the lines of one organism; {release: (version, requests asked, lines posed)};
    `releases`: only the lines the tool answers at these releases."""
    tool = T.TOOLS[tool_name]
    pass_id = K.now()
    by_rel, out = {}, {}
    for q in lines:
        rel, _ = tool.posable(q, organism)
        if rel is not None:
            by_rel.setdefault(rel, []).append(q)
    if tool_name == "gprofiler":
        rows, vs = ask_versions(tool, organism, client, pass_id, tool.version_calls(organism))
        archives = T.gprofiler_archives(vs.body([r for r in rows if r["endpoint"] == "list_archives"][0]))
        cur = version_of(tool, rows, vs)
        tool.current = int(cur[1: cur.index("_")])
        tool.instances = {r: u for r, u in archives.items()}
        by_rel = {}
        for q in lines:
            rel, _ = tool.posable(q, organism)
            if rel is not None:
                by_rel.setdefault(rel, []).append(q)
    for rel, qs in sorted(by_rel.items()):
        if releases and rel not in releases:
            continue
        t0 = time.time()
        if tool_name == "biomart":
            host = urllib.parse.urlparse(tool.instances[rel]).hostname
            if host in FAILED_HOSTS:
                log(f"biomart {organism} release {rel}: {host} failed earlier in this pass; its {len(qs)} lines are "
                    f"left for the next pass")
                out[rel] = (None, 0, len(qs))
                continue
            rows, vs = ask_versions(tool, organism, client, pass_id, tool.version_calls(organism, rel)[:1])
            if rows[0].get("error") or not rows[0].get("file"):
                FAILED_HOSTS.add(host)
                log(f"biomart {organism} release {rel}: the archive did not answer its registry ({rows[0].get('error')});"
                    f" its {len(qs)} lines are left for the next pass")
                out[rel] = (None, 0, len(qs))
                continue
            rows += ask_versions(tool, organism, client, pass_id, tool.version_calls(organism, rel)[1:])[0]
            frow = [r for r in rows if r["endpoint"] == "filters"][0]
            if frow.get("error") or not frow.get("file"):
                log(f"biomart {organism} release {rel}: the archive did not answer its version calls "
                    f"({frow.get('error')}); its {len(qs)} lines are left for the next pass")
                out[rel] = (None, 0, len(qs))
                continue
            fl = vs.body(frow).decode(errors="replace")
            filters = {x.split("\t")[0] for x in fl.splitlines() if x.strip()}
            version = str(rel)
            build = lambda xs: tool.requests(xs, organism, rel, filters)  # noqa: E731
        elif tool_name == "gprofiler":
            base = T.GPROFILER if rel == tool.current else tool.instances[rel]
            # an archived instance has no version call of its own (its data_versions answers 404): its version is
            # the name of its address in the archive list stored with the pass
            version = cur if rel == tool.current else base.rstrip("/").rsplit("/", 1)[-1]
            build = lambda xs: tool.requests(xs, organism, base)  # noqa: E731
        else:
            rows, vs = ask_versions(tool, organism, client, pass_id, tool.version_calls(organism))
            version = version_of(tool, rows, vs)
            build = (lambda xs: tool.requests(xs, organism)) if tool_name not in ("hgnc", "uniprot") else None
        store = K.RawStore.of(tool.name, version, organism)
        n = 0
        for form in T.FORMS:  # the line as it stands, and the careful user's reading (docs/families/g.md)
            fl = T.in_form(qs, form)
            if not fl:
                continue
            if tool_name == "hgnc":
                n += ask_hgnc(store, client, tool, fl, log, form)
            elif tool_name == "uniprot":
                n += ask_uniprot(store, client, fl, organism, log, form)
            else:
                reqs = build(fl)
                for r in reqs:
                    r["endpoint"] = f"{r['endpoint']}@{form}"
                n += ask(store, client, reqs, log)
        record_pass(store, pass_id=pass_id, release=rel, version=version, lines=len(qs), requests_asked=n,
                    lines_best_differs=len(T.in_form(qs, "best")),
                    started=pass_id, finished=K.now(), seconds=round(time.time() - t0, 1),
                    host_group=GROUP[tool_name], host_interval=K.HOSTS[GROUP[tool_name]].interval)
        out[rel] = (version, n, len(qs))
        log(f"{tool_name} {organism} release {rel} ({version}): {len(qs)} lines, {n} requests")
    return out


def ask_hgnc(store, client, tool, qs, log, form):
    """HGNC: per line, its fields in order (an approved symbol, then a previous, then an alias symbol, as HGNC's
    multi-symbol checker reads a term), until one finds a record."""
    n = 0
    for q in qs:
        for field in tool.chain(q):
            r = tool.request(q, field)
            r["endpoint"] = f"{r['endpoint']}@{form}"
            row = client.fetch(store, **r)
            n += 1
            if row.get("error"):
                break
            try:
                found = json.loads(store.body(row))["response"]["numFound"]
            except (KeyError, ValueError):
                found = 0
            if found:
                break
        if n % 200 == 0:
            store.compile()
    store.compile()
    return n


def ask_uniprot(store, client, qs, organism, log, form, poll=5, max_polls=720):
    """UniProt's ID mapping: one job of every accession (at most 100,000), polled until it is finished or failed,
    then its results streamed; every response stored. A job that failed (`jobStatus` ERROR) or whose results did not
    come is asked again by a new job at the next pass (its request salted by the attempt)."""
    ids = sorted({q["id"] for q in qs})
    qids = [q["qid"] for q in qs]
    if store.answered(f"idmapping/stream@{form}") >= set(qids):
        return 0
    attempt = sum(1 for r in store.rows.values() if r["endpoint"] == f"idmapping/run@{form}")
    run = client.fetch(store, "POST", f"{T.UNIPROT}/idmapping/run", data={"from": "UniProtKB_AC-ID", "to": "Ensembl",
                       "ids": ",".join(ids)}, qids=qids, endpoint=f"idmapping/run@{form}",
                       salt=str(attempt) if attempt else None)
    if not run.get("file") or run.get("error"):
        return 1
    job = json.loads(store.body(run))["jobId"]
    n, status = 1, None
    for i in range(max_polls):
        st = client.fetch(store, "GET", f"{T.UNIPROT}/idmapping/status/{job}", allow_redirects=False,
                          salt=f"{attempt}:{i}", endpoint=f"idmapping/status@{form}")
        n += 1
        body = json.loads(store.body(st)) if st.get("file") else {}
        status = body.get("jobStatus")
        if status == "FINISHED" or "results" in body or st["http_status"] == 303:
            status = "FINISHED"
            break
        if status == "ERROR":
            log(f"uniprot {organism}: job {job} failed ({str(body.get('errors'))[:200]}); asked again at the next pass")
            return n
        time.sleep(poll)
    if status == "FINISHED":
        client.fetch(store, "GET", f"{T.UNIPROT}/idmapping/stream/{job}", params={"format": "json"}, qids=qids,
                     endpoint=f"idmapping/stream@{form}")
        n += 1
    store.compile()
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.fetch", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tool", required=True, choices=sorted(T.TOOLS))
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--tier", default="web", choices=["web", "local"])
    ap.add_argument("--per-family", type=int, default=None)
    ap.add_argument("--draft", default=None)
    ap.add_argument("--releases", default=None, help="only the lines the tool answers at these releases")
    a = ap.parse_args(argv)
    rels = {int(x) for x in a.releases.split(",")} if a.releases else None
    orgs = CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")
    client = K.Client(log=lambda *x: print(*x, flush=True))
    for o in orgs:
        lines = g_lines(o, a.tier, a.per_family, a.draft)
        if lines:
            pass_tool(a.tool, o, lines, client, log=lambda *x: print(*x, flush=True), releases=rels)


if __name__ == "__main__":
    main()
