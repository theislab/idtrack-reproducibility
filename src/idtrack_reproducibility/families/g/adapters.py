"""The adapters of family G: a tool's answer to each line, read from its stored responses only (docs/families/g.md,
"From responses to answers").

`answers(tool, organism, version, lines)` -> {qid: Answer}: the genes the tool named for the line (`got`, Ensembl gene
ids without versions), the outcome (`answer`, `error` where the tool's response for the line failed or is missing),
the tool's own status where it gives one, and the line's share of its request's seconds. Each adapter reads the
response as the tool's documentation describes it; nothing is corrected by hand:

- Ensembl REST `archive/id`: a current id is its own answer (`is_current`), a retired one its `possible_replacement`
  ids; an id the response does not hold: nothing. `xrefs/symbol`: the gene objects it lists. `lookup/symbol`: the
  gene of the symbol's entry.
- BioMart: the `ensembl_gene_id` of the rows whose filter value is the line's value (compared without case, as the
  mart's filter matches); a response without its `[success]` line is an error.
- g:Profiler: the `converted` ENSG ids of the line's `incoming` value (not `N/A` or `None`, its marks of no conversion).
- MyGene.info: the `ensembl.gene` ids of every hit of the line's term; the ids of Ensembl's other annotations of the
  same species (pig's breed assemblies: ids absent from Ensembl's GTF of the tool's release) are not the organism's
  reference genes and are left out, read from that GTF.
- HGNC: the `ensembl_gene_id` of the documents of the first field of the line's chain that found one.
- UniProt: the Ensembl genes (`to`) of the line's accession (`from`), other annotations left out as for MyGene.
- BridgeDb: the `En:` ids of the line's id.
"""
import collections
import json
import os

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families.g import client as K
from idtrack_reproducibility.families.g import tools as T
from idtrack_reproducibility.harness import harness_lib as H


class Answer:
    __slots__ = ("outcome", "got", "status", "seconds", "error", "endpoint")

    def __init__(self, outcome, got=(), status=None, seconds=None, error=None, endpoint=None):
        self.outcome, self.got, self.status = outcome, sorted({C.base_id(x) for x in got}), status
        self.seconds, self.error, self.endpoint = seconds, error, endpoint


def _rows(store, endpoint=None):
    """The usable stored responses of an endpoint (client.usable): a failed or refused request answers nothing (its
    lines have no answer from it, and the next pass asks them again)."""
    return [r for r in store.rows.values() if K.usable(r) and (endpoint is None or r["endpoint"] == endpoint
            or (endpoint.endswith("*") and r["endpoint"].startswith(endpoint[:-1])))]


def _share(r):
    n = max(1, len(r.get("qids") or []))
    return (r.get("seconds") or 0.0) / n


def _by_qid(rows):
    out = {}
    for r in sorted(rows, key=lambda r: r["asked_at"] or ""):
        for q in r.get("qids") or []:
            out[q] = r
    return out


_REFERENCE = {}


def reference_genes(organism, release):
    """The organism's reference genes at a release (Ensembl's GTF of it), to tell them from other annotations."""
    key = (organism, release)
    if key not in _REFERENCE:
        _REFERENCE[key] = set(H.gtf_genes(organism, release) or ())
    return _REFERENCE[key]


def rest(store, lines, endpoint):
    out = {}
    by = _by_qid(_rows(store, endpoint))
    cache = {}
    for q in lines:
        r = by.get(q["qid"])
        if r is None:
            continue
        if r.get("error"):
            out[q["qid"]] = Answer("error", error=r["error"], seconds=_share(r), endpoint=endpoint)
            continue
        if r["request_id"] not in cache:
            cache[r["request_id"]] = json.loads(store.body(r)) if r.get("file") else None
        doc = cache[r["request_id"]]
        if endpoint == "archive/id":
            hit = {C.base_id(x["id"]): x for x in (doc or []) if isinstance(x, dict) and x.get("id")}.get(C.base_id(q["id"]))
            if hit is None:
                out[q["qid"]] = Answer("answer", [], "not_found", _share(r), endpoint=endpoint)
            elif str(hit.get("is_current")) == "1":
                out[q["qid"]] = Answer("answer", [hit["latest"]], "current", _share(r), endpoint=endpoint)
            else:
                rep = [x["stable_id"] for x in hit.get("possible_replacement") or []]
                out[q["qid"]] = Answer("answer", rep, "replaced" if rep else "retired", _share(r), endpoint=endpoint)
        elif endpoint == "xrefs/symbol":
            if r["http_status"] != 200 or not isinstance(doc, list):
                # the service could not answer the input (a 400 with its error, e.g. a database error on a Greek letter)
                err = (doc or {}).get("error") if isinstance(doc, dict) else None
                out[q["qid"]] = Answer("error", error=f"HTTP {r['http_status']}: {str(err)[:200]}", seconds=_share(r),
                                       endpoint=endpoint)
            else:
                got = [x["id"] for x in doc if x.get("type") == "gene"]
                out[q["qid"]] = Answer("answer", got, "found" if got else "not_found", _share(r), endpoint=endpoint)
        elif endpoint == "lookup/symbol":
            hit = (doc or {}).get(q["id"]) if isinstance(doc, dict) else None
            out[q["qid"]] = Answer("answer", [hit["id"]] if hit and hit.get("id") else [], "found" if hit else "not_found",
                                   _share(r), endpoint=endpoint)
    return out


def biomart(store, lines, organism):
    out = {}
    rows = _rows(store, "martservice:*")
    parsed = {}
    for r in rows:
        if r.get("error") or not r.get("file"):
            parsed[r["request_id"]] = None
            continue
        txt = store.body(r).decode(errors="replace")
        ok = txt.rstrip().endswith("[success]")
        table = collections.defaultdict(set)
        # the rows are (ensembl_gene_id, the filter's attribute); a query on ensembl_gene_id holds the id itself
        # first (its second column is the id's version)
        key_col = 0 if r["endpoint"] == "martservice:ensembl_gene_id" else 1
        if ok:
            ls = txt.rstrip().splitlines()[:-1]
            for x in ls[1:]:
                f = x.split("\t")
                if len(f) >= 2 and f[0]:
                    table[f[key_col].casefold()].add(f[0])
        parsed[r["request_id"]] = (ok, table, r)
    by_value = collections.defaultdict(list)
    for r in rows:
        p = parsed[r["request_id"]]
        for q in r.get("qids") or []:
            by_value[q].append((r, p))
    for q in lines:
        hits = by_value.get(q["qid"])
        if not hits:
            continue
        f, v = T.biomart_filter(q)
        got, secs, err = set(), 0.0, None
        for r, p in hits:
            secs += _share(r)
            if p is None or not p[0]:
                err = r.get("error") or "the response has no [success] line"
                continue
            got |= p[1].get(v.casefold(), set())
        if err and not got:
            out[q["qid"]] = Answer("error", error=err, seconds=secs, endpoint="martservice")
        else:
            out[q["qid"]] = Answer("answer", got, "found" if got else "not_found", secs, endpoint="martservice")
    return out


def gprofiler(store, lines):
    out = {}
    by = _by_qid(_rows(store, "convert"))
    cache = {}
    for q in lines:
        r = by.get(q["qid"])
        if r is None:
            continue
        if r.get("error"):
            out[q["qid"]] = Answer("error", error=r["error"], seconds=_share(r), endpoint="convert")
            continue
        if r["request_id"] not in cache:
            m = collections.defaultdict(set)
            for x in json.loads(store.body(r)).get("result") or []:
                if x.get("converted") not in (None, "", "N/A", "None"):  # g:Profiler's marks of no conversion
                    m[x["incoming"]].add(x["converted"])
                else:
                    m.setdefault(x["incoming"], set())
            cache[r["request_id"]] = m
        got = cache[r["request_id"]].get(q["id"], set())
        out[q["qid"]] = Answer("answer", got, "found" if got else "not_found", _share(r), endpoint="convert")
    return out


def mygene(store, lines, organism, release):
    out = {}
    ref = reference_genes(organism, release)
    by = _by_qid(_rows(store, "query:*"))
    cache = {}
    for q in lines:
        r = by.get(q["qid"])
        if r is None:
            continue
        if r.get("error"):
            out[q["qid"]] = Answer("error", error=r["error"], seconds=_share(r), endpoint=r["endpoint"])
            continue
        if r["request_id"] not in cache:
            m = collections.defaultdict(set)
            for h in json.loads(store.body(r)):
                g = h.get("ensembl.gene")
                g = [g] if isinstance(g, str) else (g or [])
                m[h["query"]].update(g)
            cache[r["request_id"]] = m
        allg = cache[r["request_id"]].get(q["id"], set())
        got = {g for g in allg if g in ref}
        st = "found" if got else ("other_annotation_only" if allg else "not_found")
        out[q["qid"]] = Answer("answer", got, st, _share(r), endpoint=r["endpoint"])
    return out


def hgnc(store, lines):
    out = {}
    per = collections.defaultdict(list)
    for r in _rows(store):
        for x in r.get("qids") or []:
            per[x].append(r)
    for q in lines:
        rs = sorted(per.get(q["qid"], []), key=lambda r: r["asked_at"] or "")
        if not rs:
            continue
        got, st, secs, err = set(), "not_found", sum(_share(r) for r in rs), None
        for r in rs:
            if r.get("error") or r.get("http_status") != 200:
                err = r.get("error") or f"HTTP {r.get('http_status')}: the service could not read the input"
                break
            docs = json.loads(store.body(r)).get("response", {}).get("docs") or []
            ens = {d["ensembl_gene_id"] for d in docs if d.get("ensembl_gene_id")}
            if docs:
                got, st = ens, r["endpoint"].split("/", 1)[1]
                break
        out[q["qid"]] = Answer("error", error=err, seconds=secs) if err and not got else \
            Answer("answer", got, st, secs, endpoint="fetch")
    return out


def uniprot(store, lines, organism, release):
    out = {}
    ref = reference_genes(organism, release)
    rows = sorted(_rows(store, "idmapping/stream"), key=lambda r: r["asked_at"] or "")
    m, secs, err = collections.defaultdict(set), {}, {}
    for r in rows:
        if r.get("error"):
            for q in r.get("qids") or []:
                err[q] = r["error"]
            continue
        for x in json.loads(store.body(r)).get("results") or []:
            m[x["from"]].add(x["to"])
        for q in r.get("qids") or []:
            secs[q] = _share(r)
            err.pop(q, None)
    for q in lines:
        if q["qid"] in err:
            out[q["qid"]] = Answer("error", error=err[q["qid"]])
        elif q["qid"] in secs:
            allg = {C.base_id(x) for x in m.get(q["id"], set())}
            got = {g for g in allg if g in ref}
            st = "found" if got else ("other_annotation_only" if allg else "not_found")
            out[q["qid"]] = Answer("answer", got, st, secs[q["qid"]], endpoint="idmapping")
    return out


def bridgedb(store, lines, organism):
    out = {}
    by = _by_qid(_rows(store, "xrefsBatch/*"))
    cache = {}
    for q in lines:
        r = by.get(q["qid"])
        if r is None:
            continue
        if r.get("error"):
            out[q["qid"]] = Answer("error", error=r["error"], seconds=_share(r), endpoint=r["endpoint"])
            continue
        if r["request_id"] not in cache:
            m = {}
            for x in store.body(r).decode(errors="replace").splitlines():
                f = x.split("\t")
                if len(f) >= 3:
                    m[f[0]] = [y[3:] for y in f[2].split(",") if y.startswith("En:")]
            cache[r["request_id"]] = m
        sys_ = T.TOOLS["bridgedb"].system(q, organism)
        key = q["id"]
        got = cache[r["request_id"]].get(key, [])
        out[q["qid"]] = Answer("answer", got, "found" if got else "not_found", _share(r), endpoint=r["endpoint"])
    return out


class FormView:
    """A store's rows of one form of the input (endpoints `<endpoint>@<form>`), with the form's suffix taken off, so
    that the adapters read them by their endpoints."""

    def __init__(self, store, form):
        self.folder, self.tool, self.version, self.organism = store.folder, store.tool, store.version, store.organism
        self._store = store
        self.rows = {}
        for k, r in store.rows.items():
            if r["endpoint"].endswith("@" + form):
                self.rows[k] = dict(r, endpoint=r["endpoint"][: -len(form) - 1])

    def body(self, row):
        return self._store.body(row)


def gget(store, lines):
    """gget info: an id that its result holds as a gene is its own answer (the gene at Ensembl's current release)."""
    out = {}
    by = _by_qid(_rows(store, "info"))
    cache = {}
    for q in lines:
        r = by.get(q["qid"])
        if r is None:
            continue
        if r.get("error"):
            out[q["qid"]] = Answer("error", error=r["error"], seconds=_share(r), endpoint="info")
            continue
        if r["request_id"] not in cache:
            doc = json.loads(store.body(r)) or {}
            cache[r["request_id"]] = {C.base_id(k): v for k, v in doc.items()} if isinstance(doc, dict) else {}
        hit = cache[r["request_id"]].get(C.base_id(q["id"]))
        ok = isinstance(hit, dict) and str(hit.get("object_type", "Gene")).lower() == "gene"
        out[q["qid"]] = Answer("answer", [C.base_id(q["id"])] if ok else [], "found" if ok else "not_found",
                               _share(r), endpoint="info")
    return out


def calls(store, lines, organism=None, release=None):
    """A local tool's stored calls (`calls` responses: one JSON line per line asked, `result` its Ensembl genes);
    with `release`, the ids of Ensembl's other annotations of the species left out (as for MyGene.info)."""
    import gzip as _gz
    out = {}
    ref = reference_genes(organism, release) if release else None
    for r in sorted(_rows(store, "calls"), key=lambda r: r["asked_at"] or ""):
        body = store.body(r)
        if body[:2] == b"\x1f\x8b":
            body = _gz.decompress(body)
        per = (r.get("seconds") or 0.0) / max(1, len(r.get("qids") or []))
        for x in body.decode().splitlines():
            c = json.loads(x)
            allg = {C.base_id(g) for g in c.get("result") or []}
            got = {g for g in allg if ref is None or g in ref}
            st = "found" if got else ("other_annotation_only" if allg else "not_found")
            out[c["qid"]] = Answer("answer", got, st, c.get("seconds", per), endpoint="calls")
    want = {q["qid"] for q in lines}
    return {k: v for k, v in out.items() if k in want}


def answers(tool, organism, version, lines, release=None, form="best"):
    """{qid: Answer} of one tool's stored responses (one version and organism) for the lines, in one form of the
    input (`best`: the lines whose value changes are read from the `best` requests, the others from the `stated`)."""
    if form == "best":
        changed = T.in_form(lines, "best")
        ids = {q["qid"] for q in changed}
        out = _answers(tool, organism, version, T.in_form([q for q in lines if q["qid"] not in ids], "stated"),
                       release, "stated")
        out.update(_answers(tool, organism, version, changed, release, "best"))
        return out
    return _answers(tool, organism, version, T.in_form(lines, "stated"), release, "stated")


def _answers(tool, organism, version, lines, release, form):
    store = FormView(K.RawStore.of("ensembl_rest" if tool == "ensembl_rest_lookup" else tool, version, organism), form)
    if tool == "ensembl_rest":
        out = rest(store, [q for q in lines if q["task"] in ("a", "b")], "archive/id")
        out.update(rest(store, [q for q in lines if q["task"] in ("c", "d", "e")], "xrefs/symbol"))
        return out
    if tool == "ensembl_rest_lookup":
        return rest(store, lines, "lookup/symbol")
    if tool == "biomart":
        return biomart(store, lines, organism)
    if tool == "gprofiler":
        return gprofiler(store, lines)
    if tool == "mygene":
        return mygene(store, lines, organism, release or 116)
    if tool == "hgnc":
        return hgnc(store, lines)
    if tool == "uniprot":
        return uniprot(store, lines, organism, release or 116)
    if tool == "bridgedb":
        return bridgedb(store, lines, organism)
    if tool == "gget":
        return gget(store, lines)
    if tool in ("pyensembl", "orgdb", "limma", "hgnchelper", "seurat", "ensdb"):
        return calls(store, lines)
    if tool == "ncbi_files":
        return calls(store, lines, organism, release or 116)
    raise ValueError(tool)


def raw_versions(tool, organism):
    """The version folders of a tool's stored responses for an organism."""
    d = os.path.join(C.family_dir("g"), "raw", tool)
    if not os.path.isdir(d):
        return []
    return sorted(v for v in os.listdir(d) if v != "_versions" and os.path.isdir(os.path.join(d, v, organism)))
