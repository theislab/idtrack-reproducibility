"""The local tools of family G asked from Python: pyensembl (docs/families/g.md, "The tools").

Usage: python -m idtrack_reproducibility.families.g.local --tool pyensembl --organisms <o,...|all> [--tier local]
                                                           [--per-family N] [--draft <folder>] [--releases 115,...]

pyensembl reads one Ensembl release's GTF. The GTF of each release a line needs is downloaded by the client from the
address pyensembl itself builds (`ensembl_url_templates.make_gtf_url`, Ensembl's FTP site over HTTPS) into the raw
store (`<families>/g/raw/pyensembl/<release>/<organism>/`, with its URL, size, sha256, md5 and Last-Modified),
then given to `pyensembl.Genome(..., gtf_path_or_url=<that file>)`, which indexes it into its cache
(`PYENSEMBL_CACHE_DIR` = `<families>/g/raw/pyensembl/_cache`, rebuilt from the GTF when missing). The transcript and
protein FASTA files pyensembl's installer also fetches are not needed for genes and names and are not fetched.

The calls, one per line (task b: `gene_by_id` of the id in each form of the input (tools.in_form): the id itself when the release
holds it; task c: `gene_ids_of_gene_name`, the genes named so in the GTF), are stored as the tool returned them
(`calls_<pass>.jsonl.gz`: qid, function, argument, result or the error's class and message, seconds) and indexed in
`index.parquet` like a web service's responses, so the adapter reads stored answers only.
"""
import argparse
import gzip
import json
import os
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families.g import client as K
from idtrack_reproducibility.families.g import fetch as F
from idtrack_reproducibility.families.g import tools as T

SPECIES = {"sus_scrofa": "sus_scrofa", "mus_musculus": "mus_musculus", "homo_sapiens": "homo_sapiens",
           "macaca_mulatta": "macaca_mulatta", "macaca_fascicularis": "macaca_fascicularis"}
TASKS = ("b", "c")


def posable(line, organism):
    if line.get("task") not in TASKS:
        return None, f"pyensembl does not claim task {line.get('task')}"
    r = line.get("to_release")
    if r is None or r < 40 or r > 116:
        return None, f"no Ensembl release {r} in pyensembl's range (40 to 116)"
    if line["task"] == "b" and line.get("from_release") != r:
        return None, "an id at another release than the line's target (history) is not pyensembl's task"
    from pyensembl.species import find_species_by_name
    try:
        find_species_by_name(SPECIES[organism]).which_reference(r)
    except ValueError as e:  # pyensembl's own table of the releases of each species' assemblies
        return None, f"pyensembl: {e}"
    return r, ""


def genome(organism, release, client):
    os.environ.setdefault("PYENSEMBL_CACHE_DIR", os.path.join(C.family_dir("g"), "raw", "pyensembl", "_cache"))
    import pyensembl
    from pyensembl import ensembl_url_templates as U
    from pyensembl.species import find_species_by_name
    sp = find_species_by_name(SPECIES[organism])
    url = U.make_gtf_url(release, sp)
    store = K.RawStore.of("pyensembl", release, organism)
    row = K.download(store, url, endpoint="gtf", version=f"{pyensembl.__version__}:{release}", client=client)
    store.compile()
    if row.get("error"):
        raise SystemExit(f"{url}: {row['error']}")
    g = pyensembl.Genome(reference_name=sp.which_reference(release), annotation_name="ensembl",
                         annotation_version=release, gtf_path_or_url=os.path.join(store.folder, row["file"]))
    g.index()
    return g, store, url


def ask(g, store, lines, pass_id, form):
    """The calls of the lines in one form of the input (tools.in_form), stored as one response."""
    name = f"calls_{form}_{pass_id.replace(':', '')}"
    done = store.answered(f"calls@{form}")
    todo = [q for q in lines if q["qid"] not in done]
    if not todo:
        return 0
    t_all = time.monotonic()
    out = []
    for q in todo:
        t0 = time.monotonic()
        fn = "gene_by_id" if q["task"] == "b" else "gene_ids_of_gene_name"
        arg = q["id"]
        try:
            res, err = ([g.gene_by_id(arg).gene_id] if fn == "gene_by_id" else sorted(g.gene_ids_of_gene_name(arg))), None
        except ValueError as e:  # pyensembl's documented answer for an id or name the release does not hold
            res, err = [], f"ValueError: {e}"
        out.append(json.dumps(dict(qid=q["qid"], function=fn, argument=arg, result=res, error=err,
                                   seconds=round(time.monotonic() - t0, 6))) + "\n")
    content = gzip.compress("".join(out).encode(), mtime=0)  # the calls as one response, kept under its sha256
    row = dict(request_id=K.request_id("CALL", name, {"pass": pass_id}), tool="pyensembl", tool_version=store.version,
               organism=store.organism, endpoint=f"calls@{form}", method="CALL", url=name,
               params={"functions": ["gene_by_id", "gene_ids_of_gene_name"]}, body=None,
               qids=[q["qid"] for q in todo], asked_at=pass_id, http_status=None,
               seconds=round(time.monotonic() - t_all, 3), bytes=None, sha256=None, file=None, error=None)
    store.put(row, content, "jsonl.gz", {})
    store.compile()
    return len(todo)


class NcbiFiles:
    """NCBI's gene files of the truth store (read only), resolved as NCBI's README describes: `gene_info` (Symbol,
    Synonyms), `gene_history` (a discontinued GeneID to its current one), `gene2refseq` (an accession to its GeneID),
    `gene2ensembl` (a GeneID to Ensembl's genes)."""

    def __init__(self, organism):
        from idtrack_reproducibility.harness import harness_lib as H
        import glob
        d = os.path.join(H.truth_dir(organism), "ncbi")
        self.files = sorted(os.path.relpath(p, H.truth_dir(organism)) for p in glob.glob(os.path.join(d, "*")))
        srcs = {r["what"]: r for r in H.read_tsv(os.path.join(H.truth_dir(organism), "SOURCES.tsv"))}
        dates = sorted({(srcs.get(f) or {}).get("fetched", "")[:10] for f in self.files} - {""})
        self.version = dates[-1] if dates else "unknown"
        info = glob.glob(os.path.join(d, "*.gene_info.gz"))[0]
        self.current, self.symbol, self.synonym = set(), {}, {}
        for r in H.read_tsv(info):
            g = r["GeneID"]
            self.current.add(g)
            self.symbol.setdefault(r["Symbol"], set()).add(g)
            for x in (r.get("Synonyms") or "-").split("|"):
                if x != "-":
                    self.synonym.setdefault(x, set()).add(g)
        self.history, self.disc_symbol = {}, {}
        for r in H.read_tsv(os.path.join(d, "gene_history.tsv.gz")):
            self.history[r["Discontinued_GeneID"]] = r["GeneID"]
            self.disc_symbol.setdefault(r["Discontinued_Symbol"], set()).add(r["GeneID"])
        self.ens = {}
        for r in H.read_tsv(os.path.join(d, "gene2ensembl.tsv.gz")):
            self.ens.setdefault(r["GeneID"], set()).add(r["Ensembl_gene_identifier"])
        self.refseq = {}
        for r in H.read_tsv(os.path.join(d, "gene2refseq.tsv.gz")):
            for k in ("RNA_nucleotide_accession.version", "protein_accession.version"):
                a = r.get(k)
                if a and a != "-":
                    self.refseq.setdefault(a, set()).add(r["GeneID"])
                    self.refseq.setdefault(a.split(".")[0], set()).add(r["GeneID"])

    def call(self, line):
        """(function, GeneIDs, Ensembl genes) of one line as asked (its `id` in the form asked)."""
        v, t = line["id"], line["task"]
        kind = T.source_kind(line)
        if t in ("c", "d"):
            g = set(self.symbol.get(v, ()))
            fn = "gene_info Symbol"
            if not g and t == "d":
                g = set(self.synonym.get(v, ())) or {x for x in self.disc_symbol.get(v, ()) if x != "-"}
                fn = "gene_info Synonyms, then gene_history Discontinued_Symbol"
        elif kind == "ncbi":
            fn, g = "gene_info GeneID, then gene_history", set()
            if v in self.current:
                g = {v}
            elif self.history.get(v) not in (None, "-"):
                g = {self.history[v]}
        else:  # RefSeq: the accession as asked (with its version where it has one)
            fn, g = "gene2refseq", set(self.refseq.get(v, ()))
        return fn, sorted(g), sorted({e for x in g for e in self.ens.get(x, ())})


def ncbi_posable(line, organism):
    if line.get("task") not in ("c", "d", "e"):
        return None, f"NCBI's files do not claim task {line.get('task')}"
    if line.get("task") == "e" and T.source_kind(line) not in ("ncbi", "refseq"):
        return None, "NCBI's files map NCBI gene ids and RefSeq accessions, not this database"
    if line.get("to_release") != 115:
        return None, f"NCBI's files of the day answer for today, not at release {line.get('to_release')}"
    return 116, ""


def ask_ncbi(organism, lines, pass_id):
    nf = NcbiFiles(organism)
    store = K.RawStore.of("ncbi_files", nf.version, organism)
    n = 0
    for form in T.FORMS:
        fl = [q for q in T.in_form(lines, form) if q["qid"] not in store.answered(f"calls@{form}")]
        if not fl:
            continue
        rows = []
        for q in fl:
            fn, gid, ens = nf.call(q)
            rows.append(dict(qid=q["qid"], function=fn, argument=q["id"], geneids=gid, result=ens, error=None))
        content = gzip.compress("".join(json.dumps(r) + "\n" for r in rows).encode(), mtime=0)
        row = dict(request_id=K.request_id("CALL", "ncbi_files", {"form": form, "pass": pass_id}), tool="ncbi_files",
                   tool_version=nf.version, organism=organism, endpoint=f"calls@{form}", method="CALL",
                   url="truth store: " + ",".join(nf.files), params={"files": nf.files}, body=None,
                   qids=[q["qid"] for q in fl], asked_at=pass_id, http_status=None, seconds=None, bytes=None,
                   sha256=None, file=None, error=None)
        store.put(row, content, "jsonl.gz", {})
        n += len(fl)
    store.compile()
    return n


R_TOOLS = {
    # tool: (organisms, tasks, the release it is judged at, tier)
    "orgdb": (("homo_sapiens", "mus_musculus", "macaca_mulatta"), ("c", "d", "e"), 115, "local"),
    "limma": (("homo_sapiens", "mus_musculus"), ("c", "d"), 115, "local"),
    "hgnchelper": (("homo_sapiens", "mus_musculus"), ("c", "d"), 115, "local"),
    "seurat": (("homo_sapiens",), ("c", "d"), 116, "web"),
    "ensdb": (("sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis"), ("b", "c", "e"), 113, "local"),
}


def r_posable(tool, line, organism):
    orgs, tasks, rel, _ = R_TOOLS[tool]
    if organism not in orgs:
        return None, {"orgdb": "no org package with Ensembl genes for this organism (org.Ss.eg.db has no ENSEMBL column; "
                               "no cynomolgus package)"}.get(tool, f"{tool} has no data for {organism}")
    if line.get("task") not in tasks:
        return None, f"{tool} does not claim task {line.get('task')}"
    if tool == "orgdb" and organism == "macaca_mulatta" and line["task"] == "d":
        return None, "org.Mmu.eg.db has no ALIAS column"
    k = T.source_kind(line)
    if line["task"] == "e" and (k not in ("ncbi", "refseq", "uniprot") or (tool == "ensdb" and k != "ncbi")):
        return None, f"{tool} maps no identifiers of this database"
    if line.get("to_release") != 115 and not (tool == "ensdb" and line.get("to_release") == rel):
        return None, f"{tool} answers at its own snapshot, not at release {line.get('to_release')}"
    if line["task"] == "b" and line.get("from_release") != line.get("to_release"):
        return None, "an id at another release (history) is not its task"
    return rel, ""


def ask_r(tool, organism, lines, rscript, pass_id, workdir):
    """One R tool for one organism: per form, the lines written as JSON, `r_tools.R` run, its output stored as the
    tool's response (one JSON line per line asked)."""
    import subprocess
    here = os.path.dirname(os.path.abspath(__file__))
    n, store = 0, None
    for form in T.FORMS:
        fl = T.in_form(lines, form)
        if not fl:
            continue
        stamp = pass_id.replace(":", "")  # the files of a pass are its own: no pass writes over another's
        inp = os.path.join(workdir, f"{tool}_{organism}_{form}_{stamp}.lines.json")
        outp = os.path.join(workdir, f"{tool}_{organism}_{form}_{stamp}.out.json")
        json.dump([dict(qid=q["qid"], id=q["id"], task=q["task"], kind=T.source_kind(q)) for q in fl], open(inp, "w"))
        t0 = time.monotonic()
        env = dict(os.environ)
        env.setdefault("ANNOTATION_HUB_CACHE", os.path.join(C.family_dir("g"), "raw", "ensdb", "_cache"))
        subprocess.run([rscript, os.path.join(here, "r_tools.R"), tool, organism, inp, outp], check=True, env=env)
        res = json.load(open(outp))
        version = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in res["version"])
        store = K.RawStore.of(tool, version, organism)
        content = gzip.compress("".join(json.dumps(dict(c, result=list(c.get("result") or []))) + "\n"
                                        for c in res["calls"]).encode(), mtime=0)
        row = dict(request_id=K.request_id("CALL", f"r_tools.R {tool}", {"form": form, "pass": pass_id}), tool=tool,
                   tool_version=res["version"], organism=organism, endpoint=f"calls@{form}", method="CALL",
                   url=f"r_tools.R {tool}", params={"package_metadata": res.get("package_metadata")}, body=None,
                   qids=[q["qid"] for q in fl], asked_at=pass_id, http_status=None,
                   seconds=round(time.monotonic() - t0, 3), bytes=None, sha256=None, file=None, error=None)
        if store.answered(f"calls@{form}") >= {q["qid"] for q in fl}:
            continue
        store.put(row, content, "jsonl.gz", {})
        store.compile()
        n += len(fl)
    return n


GGET_BATCH = 200  # batches of 1,000 were answered 500 and 504 by Ensembl REST (2026-11-24)


def gget_posable(line, organism):
    if line.get("task") != "b":
        return None, f"gget info does not claim task {line.get('task')} (an id's gene at the current release only)"
    if line.get("to_release") != 115:
        return None, f"gget info answers at Ensembl's current release only, not at {line.get('to_release')}"
    return 116, ""


def ask_gget(organism, lines, client, pass_id):
    """`gget info` (Ensembl REST's lookup through gget's own code, without its NCBI, UniProt and PDB parts) for the
    lines' ids, 200 per call at Ensembl REST's interval (the same host as the REST pass: one process asks it); each
    call's result (`json=True`) stored as the tool returned it."""
    import importlib.metadata as md
    import gget
    version = md.version("gget")
    store = K.RawStore.of("gget", version, organism)
    n = 0
    for form in T.FORMS:
        fl = T.in_form(lines, form)
        done = store.answered(f"info@{form}")
        todo = [q for q in fl if q["qid"] not in done]
        for i in range(0, len(todo), GGET_BATCH):
            b = todo[i: i + GGET_BATCH]
            ids = sorted({q["id"] for q in b})
            client._wait_turn("ensembl_rest")
            t0, err, res = time.monotonic(), None, None
            try:
                res = gget.info(ids, ncbi=False, uniprot=False, pdb=False, json=True, verbose=False)
            except Exception as e:  # noqa: BLE001  (the tool's failure is its answer, recorded)
                err = f"{type(e).__name__}: {e}"
            content = json.dumps(res, sort_keys=True, default=str).encode()
            row = dict(request_id=K.request_id("CALL", "gget.info", {"form": form}, ids), tool="gget", tool_version=version,
                       organism=organism, endpoint=f"info@{form}", method="CALL", url="gget.info",
                       params={"ncbi": False, "uniprot": False, "pdb": False, "json": True}, body=ids,
                       qids=[q["qid"] for q in b], asked_at=K.now(), http_status=None,
                       seconds=round(time.monotonic() - t0, 3), bytes=None, sha256=None, file=None, error=err)
            store.put(row, content, "json", {})
            n += 1
    store.compile()
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.local", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tool", default="pyensembl", choices=["pyensembl", "gget", "ncbi_files"] + list(R_TOOLS))
    ap.add_argument("--rscript", default="Rscript", help="the Rscript of the R environment (environments/g_r.yml)")
    ap.add_argument("--workdir", default=None, help="where the R tools' input and output files are written (default: "
                    "<families>/g/raw/_r_work, kept with the store)")
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--tier", default="local", choices=["web", "local"])
    ap.add_argument("--per-family", type=int, default=None)
    ap.add_argument("--draft", default=None)
    ap.add_argument("--releases", default=None)
    a = ap.parse_args(argv)
    rels = {int(x) for x in a.releases.split(",")} if a.releases else None
    client = K.Client()
    for o in (CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")):
        if a.tool in R_TOOLS:
            tier = R_TOOLS[a.tool][3]
            qs = [q for q in F.g_lines(o, tier, a.per_family, a.draft) if r_posable(a.tool, q, o)[0] is not None]
            if qs:
                wd = a.workdir or os.path.join(C.family_dir("g"), "raw", "_r_work")
                os.makedirs(wd, exist_ok=True)
                print(f"{a.tool} {o}: {len(qs)} lines, {ask_r(a.tool, o, qs, a.rscript, K.now(), wd)} answered", flush=True)
            continue
        if a.tool == "ncbi_files":
            qs = [q for q in F.g_lines(o, a.tier, a.per_family, a.draft) if ncbi_posable(q, o)[0] is not None]
            print(f"ncbi_files {o}: {len(qs)} lines, {ask_ncbi(o, qs, K.now())} answered", flush=True)
            continue
        if a.tool == "gget":
            qs = [q for q in F.g_lines(o, a.tier, a.per_family, a.draft) if gget_posable(q, o)[0] is not None]
            print(f"gget {o}: {len(qs)} lines, {ask_gget(o, qs, client, K.now())} calls", flush=True)
            continue
        by = {}
        for q in F.g_lines(o, a.tier, a.per_family, a.draft):
            r, _ = posable(q, o)
            if r is not None and (not rels or r in rels):
                by.setdefault(r, []).append(q)
        for r, qs in sorted(by.items()):
            g, store, url = genome(o, r, client)
            n = sum(ask(g, store, T.in_form(qs, form), K.now(), form) for form in T.FORMS)
            print(f"pyensembl {o} {r}: {len(qs)} lines, {n} calls", flush=True)


if __name__ == "__main__":
    main()
