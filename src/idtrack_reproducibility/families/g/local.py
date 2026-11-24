"""The local tools of family G asked from Python: pyensembl (docs/families/g.md, "The tools").

Usage: python -m idtrack_reproducibility.families.g.local --tool pyensembl --organisms <o,...|all> [--tier local]
                                                           [--per-family N] [--draft <folder>] [--releases 115,...]

pyensembl reads one Ensembl release's GTF. The GTF of each release a line needs is downloaded by the client from the
address pyensembl itself builds (`ensembl_url_templates.make_gtf_url`, Ensembl's FTP site over HTTPS) into the raw
store (`<families>/g/raw/pyensembl/<release>/<organism>/`, with its URL, size, sha256, md5 and Last-Modified),
then given to `pyensembl.Genome(..., gtf_path_or_url=<that file>)`, which indexes it into its cache
(`PYENSEMBL_CACHE_DIR` = `<families>/g/raw/pyensembl/_cache`, rebuilt from the GTF when missing). The transcript and
protein FASTA files pyensembl's installer also fetches are not needed for genes and names and are not fetched.

The calls, one per line (task b: `gene_by_id` of the id without its version: the id itself when the release holds
it; task c: `gene_ids_of_gene_name`, the genes named so in the GTF), are stored as the tool returned them
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


def ask(g, store, lines, pass_id):
    out = os.path.join(store.folder, f"calls_{pass_id.replace(':', '')}.jsonl.gz")
    done = store.answered("calls")
    todo = [q for q in lines if q["qid"] not in done]
    if not todo:
        return 0
    t_all = time.monotonic()
    with gzip.open(out + ".tmp", "wt") as fh:
        for q in todo:
            t0 = time.monotonic()
            if q["task"] == "b":
                fn, arg = "gene_by_id", C.base_id(q["id"])
            else:
                fn, arg = "gene_ids_of_gene_name", q["id"]
            try:
                if fn == "gene_by_id":
                    res, err = [g.gene_by_id(arg).gene_id], None
                else:
                    res, err = sorted(g.gene_ids_of_gene_name(arg)), None
            except ValueError as e:  # pyensembl's documented answer for an id or name the release does not hold
                res, err = [], f"ValueError: {e}"
            fh.write(json.dumps(dict(qid=q["qid"], function=fn, argument=arg, result=res, error=err,
                                     seconds=round(time.monotonic() - t0, 6))) + "\n")
    os.replace(out + ".tmp", out)
    with open(out, "rb") as fh:
        content = fh.read()
    row = dict(request_id=K.request_id("CALL", out, {"pass": pass_id}), tool="pyensembl", tool_version=store.version,
               organism=store.organism, endpoint="calls", method="CALL", url=os.path.basename(out),
               params={"functions": ["gene_by_id", "gene_ids_of_gene_name"]}, body=None,
               qids=[q["qid"] for q in todo], asked_at=pass_id, http_status=None,
               seconds=round(time.monotonic() - t_all, 3), bytes=None, sha256=None, file=None, error=None)
    store.put(row, content, "jsonl.gz", {})
    os.remove(out)  # its content is in the store under its sha256
    store.compile()
    return len(todo)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.families.g.local", description=__doc__.split("\n\n")[0])
    ap.add_argument("--tool", default="pyensembl", choices=["pyensembl"])
    ap.add_argument("--organisms", default="all")
    ap.add_argument("--tier", default="local", choices=["web", "local"])
    ap.add_argument("--per-family", type=int, default=None)
    ap.add_argument("--draft", default=None)
    ap.add_argument("--releases", default=None)
    a = ap.parse_args(argv)
    rels = {int(x) for x in a.releases.split(",")} if a.releases else None
    client = K.Client()
    for o in (CFG.ORGANISMS if a.organisms == "all" else a.organisms.split(",")):
        by = {}
        for q in F.g_lines(o, a.tier, a.per_family, a.draft):
            r, _ = posable(q, o)
            if r is not None and (not rels or r in rels):
                by.setdefault(r, []).append(q)
        for r, qs in sorted(by.items()):
            g, store, url = genome(o, r, client)
            n = ask(g, store, qs, K.now())
            print(f"pyensembl {o} {r}: {len(qs)} lines, {n} asked", flush=True)


if __name__ == "__main__":
    main()
