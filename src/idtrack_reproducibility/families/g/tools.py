"""The tools of family G: how each is asked, at which release, which lines it can be asked (docs/families/g.md,
"The tools").

Each tool says its organisms, the tasks it claims (sets.task_of: `a` an old Ensembl id to its gene, `b` an Ensembl id
at a release, `c` a symbol, `d` a previous symbol, alias or synonym, `e` an external identifier), the releases it
answers at (`instances`: {release: where}), its version call, and how a list of lines becomes its requests
(`requests`: one dict per request: method, url, params, data or json_body, qids, endpoint, ext, headers). A line is
asked of a tool only when the tool can pose it (`posable`): the organism and the task are the tool's, and the line's
release is one the tool answers at: a line asked to release 115 (the graphs' release) is asked at the tool's own
release (its instance at 115 where it has one; else its current one, compared by the rule of equal truth,
likeforlike.py), a line asked at another release R only of an instance at R. Everything else is `not_supported`,
with the reason.
"""
import json
import re
import shlex
import urllib.parse
from xml.sax.saxutils import quoteattr

TARGET = 115
REST = "https://rest.ensembl.org"
REST_SPECIES = {o: o for o in ("sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis")}
BIOMART_DATASET = {"sus_scrofa": "sscrofa_gene_ensembl", "mus_musculus": "mmusculus_gene_ensembl",
                   "homo_sapiens": "hsapiens_gene_ensembl", "macaca_mulatta": "mmulatta_gene_ensembl",
                   "macaca_fascicularis": "mfascicularis_gene_ensembl"}
# Ensembl's archive hosts whose BioMart answers (probed 2026-11-24; docs/families/g.md): release -> host.
BIOMART_HOSTS = {105: "dec2021", 106: "apr2022", 107: "jul2022", 108: "oct2022", 109: "feb2023", 110: "jul2023",
                 111: "jan2024", 112: "may2024", 113: "oct2024", 114: "may2025", 115: "sep2025", 116: "jun2026"}
GPROFILER = "https://biit.cs.ut.ee/gprofiler"
GPROFILER_ORGANISM = {"sus_scrofa": "sscrofa", "mus_musculus": "mmusculus", "homo_sapiens": "hsapiens",
                      "macaca_mulatta": "mmulatta", "macaca_fascicularis": "mfascicularis"}
MYGENE = "https://mygene.info/v3"
TAXID = {"sus_scrofa": 9823, "mus_musculus": 10090, "homo_sapiens": 9606, "macaca_mulatta": 9544,
         "macaca_fascicularis": 9541}
HGNC = "https://rest.genenames.org"
UNIPROT = "https://rest.uniprot.org"
BRIDGEDB = "https://webservice.bridgedb.org"
BRIDGEDB_ORGANISM = {"sus_scrofa": "Sus scrofa", "mus_musculus": "Mus musculus", "homo_sapiens": "Homo sapiens",
                     "macaca_mulatta": "Macaca mulatta"}
JSON = {"Content-Type": "application/json", "Accept": "application/json"}


def source_kind(line):
    """What the line's input is, for a tool that needs to be told (an external identifier's database)."""
    src, st = line.get("source_family"), line.get("source_stratum") or ""
    if src == "c":
        return {"c1": "ncbi", "c2": "refseq", "c3": "uniprot"}.get(st[:2], "external")
    return {"a": "ensembl", "b": "symbol", "astar": "ensembl"}.get(src) or ("ensembl" if line.get("kind") == "d1_id" else "symbol")


# An identifier with its version: an Ensembl stable id, a RefSeq accession or a UniProt accession, then `.<n>`.
VERSIONED = re.compile(r"^((?:ENS[A-Z]*[GTPE]\d{11})|(?:[A-Z]{2}_\d+)|(?:[OPQ]\d[A-Z\d]{3}\d)|"
                       r"(?:[A-NR-Z]\d(?:[A-Z][A-Z\d]{2}\d){1,2}))\.\d+$")
FORMS = ("best", "stated")


def value(line, form):
    """What a tool is asked for a line in a form of the input (docs/families/g.md, "Two readings of the input"):
    `stated`, the line's input exactly as it stands; `best`, the reading a careful user gives a tool that takes no
    versions (none of the tools asked documents versioned input): an identifier's version removed (tasks a, b, e;
    names are asked as they stand)."""
    v = line["id"]
    if form == "best" and line.get("task") in ("a", "b", "e"):
        m = VERSIONED.match(v)
        if m:
            return m.group(1)
    return v


def in_form(lines, form):
    """The lines as a tool is asked them in a form: copies whose `id` is the value asked; for `best`, only the lines
    whose value differs from the stated one (the others share the stated requests)."""
    out = []
    for q in lines:
        v = value(q, form)
        if form == "stated" or v != q["id"]:
            out.append(dict(q, id=v, stated_id=q["id"]))
    return out


def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i: i + n]


class Tool:
    name = ""
    organisms = ()
    tasks = ()
    instances = {}      # release -> where (host, url)
    current = None      # the release a line to 115 is asked at when there is no instance at 115

    def release_for(self, line):
        """The release the tool answers this line at, or None (and why)."""
        r = line.get("to_release")
        if r == TARGET:
            if TARGET in self.instances:
                return TARGET, ""
            return self.current, ""
        if r in self.instances:
            return r, ""
        return None, f"no instance of {self.name} at release {r}"

    def posable(self, line, organism):
        if organism not in self.organisms:
            return None, f"{self.name} has no data for {organism}"
        if line.get("task") not in self.tasks:
            return None, f"{self.name} does not claim task {line.get('task')}"
        rel, why = self.release_for(line)
        if rel is None:
            return None, why
        why = self.input_problem(line, organism)
        if why:
            return None, why
        return rel, ""

    def input_problem(self, line, organism):
        return ""


class EnsemblRest(Tool):
    name = "ensembl_rest"
    organisms = tuple(REST_SPECIES)
    tasks = ("a", "b", "c", "d", "e")
    instances = {}
    current = 116  # read again from info/data at every pass (version)

    def version_calls(self, organism):
        return [dict(method="GET", url=f"{REST}/info/data/", headers=JSON, endpoint="info/data"),
                dict(method="GET", url=f"{REST}/info/software", headers=JSON, endpoint="info/software")]

    def requests(self, lines, organism):
        sp = REST_SPECIES[organism]
        ids = [q for q in lines if q["task"] in ("a", "b")]
        names = [q for q in lines if q["task"] in ("c", "d", "e")]
        out = []
        for b in chunks(ids, 400):  # the wiki's advice for archive lookups ("up to 400 IDs per batch")
            out.append(dict(method="POST", url=f"{REST}/archive/id", json_body={"id": [q["id"] for q in b]},
                            headers=JSON, qids=[q["qid"] for q in b], endpoint="archive/id"))
        for b in chunks([q for q in names if q["task"] in ("c", "d")], 200):  # 20 symbols took 6 s in the pilot
            out.append(dict(method="POST", url=f"{REST}/lookup/symbol/{sp}", json_body={"symbols": [q["id"] for q in b]},
                            headers=JSON, qids=[q["qid"] for q in b], endpoint="lookup/symbol"))
        for q in names:
            out.append(dict(method="GET", url=f"{REST}/xrefs/symbol/{sp}/{urllib.parse.quote(q['id'], safe='')}",
                            params={"object_type": "gene"}, headers=JSON, qids=[q["qid"]], endpoint="xrefs/symbol"))
        return out

    def input_problem(self, line, organism):
        if "/" in line["id"] and line["task"] in ("c", "d", "e"):
            return "a string with a slash cannot be a path segment of xrefs/symbol"
        return ""


BIOMART_REFSEQ = (("NM_", "refseq_mrna"), ("XM_", "refseq_mrna_predicted"), ("NR_", "refseq_ncrna"),
                  ("XR_", "refseq_ncrna_predicted"), ("NP_", "refseq_peptide"), ("XP_", "refseq_peptide_predicted"))


def biomart_filter(line):
    """(filter, value) BioMart is asked with for the line, or (None, why)."""
    k = source_kind(line)
    t = line["task"]
    v = line["id"]
    if t in ("a", "b"):
        return "ensembl_gene_id", v
    if t == "c":
        return "external_gene_name", v
    if t == "d":
        return "external_synonym", v
    if k == "ncbi":
        return ("entrezgene_id", v) if v.isdigit() else (None, "not an NCBI gene id")
    if k == "refseq":
        for pre, f in BIOMART_REFSEQ:
            if v.startswith(pre):
                return f, v
        return None, "a RefSeq identifier of no prefix BioMart filters on"
    if k == "uniprot":
        return "uniprot", v
    return None, "no filter"


class BioMart(Tool):
    name = "biomart"
    organisms = tuple(BIOMART_DATASET)
    tasks = ("b", "c", "d", "e")
    instances = {r: f"https://{h}.archive.ensembl.org/biomart/martservice" for r, h in BIOMART_HOSTS.items()}

    def version_calls(self, organism, release=TARGET):
        u = self.instances[release]
        ds = BIOMART_DATASET[organism]
        return [dict(method="GET", url=u, params={"type": "registry"}, ext="xml", endpoint="registry"),
                dict(method="GET", url=u, params={"type": "version", "mart": "ENSEMBL_MART_ENSEMBL"}, ext="txt",
                     endpoint="version"),
                dict(method="GET", url=u, params={"type": "filters", "dataset": ds}, ext="tsv", endpoint="filters"),
                dict(method="GET", url=u, params={"type": "attributes", "dataset": ds}, ext="tsv", endpoint="attributes")]

    def input_problem(self, line, organism):
        f, v = biomart_filter(line)
        if f is None:
            return v
        if "," in v:
            return "a value with a comma cannot be one of BioMart's comma-separated filter values"
        return ""

    def query(self, organism, filt, values, attrs):
        ds = BIOMART_DATASET[organism]
        a = "".join(f"<Attribute name={quoteattr(x)} />" for x in attrs)
        return ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" '
                'formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6" completionStamp="1">'
                f'<Dataset name={quoteattr(ds)} interface="default"><Filter name={quoteattr(filt)} '
                f'value={quoteattr(",".join(values))} />{a}</Dataset></Query>')

    def requests(self, lines, organism, release=TARGET, filters=None):
        """One XML query per filter and 500 values (about 500 values per query, the advice of biomaRt's maintainer);
        UniProt accessions are asked of Swiss-Prot's and TrEMBL's filters both. `filters`: the dataset's filters at
        that release (a line whose filter the dataset lacks is left for `not_supported`)."""
        u = self.instances[release]
        groups = {}
        for q in lines:
            f, v = biomart_filter(q)
            for ff in (("uniprotswissprot", "uniprotsptrembl") if f == "uniprot" else (f,)):
                if filters is None or ff in filters:
                    groups.setdefault(ff, {}).setdefault(v, []).append(q["qid"])
        out = []
        for f, vals in sorted(groups.items()):
            keys = sorted(vals)
            for b in chunks(keys, 500):
                attrs = ["ensembl_gene_id"] + ([f] if f != "ensembl_gene_id" else ["ensembl_gene_id_version"])
                out.append(dict(method="POST", url=u, data={"query": self.query(organism, f, b, attrs)}, ext="tsv",
                                qids=sorted({x for v in b for x in vals[v]}), endpoint=f"martservice:{f}"))
        return out


class GProfiler(Tool):
    name = "gprofiler"
    organisms = tuple(GPROFILER_ORGANISM)
    tasks = ("b", "c", "d", "e")
    current = 114
    instances = {}  # the archives (api/util/list_archives) at their releases, added at a pass

    def version_calls(self, organism, base=GPROFILER):
        return [dict(method="GET", url=f"{base}/api/util/data_versions", headers=JSON, endpoint="data_versions"),
                dict(method="GET", url=f"{GPROFILER}/api/util/list_archives", headers=JSON, endpoint="list_archives")]

    def requests(self, lines, organism, base=GPROFILER):
        out = []
        for b in chunks(lines, 1000):
            body = {"organism": GPROFILER_ORGANISM[organism], "query": [q["id"] for q in b], "target": "ENSG",
                    "numeric_ns": "ENTREZGENE_ACC"}
            out.append(dict(method="POST", url=f"{base}/api/convert/convert/", json_body=body, headers=JSON,
                            qids=[q["qid"] for q in b], endpoint="convert"))
        return out


MYGENE_SCOPES = {"b": "ensembl.gene", "c": "symbol", "d": "symbol,alias", "ncbi": "entrezgene,retired",
                 "refseq": "refseq", "uniprot": "uniprot"}


class MyGene(Tool):
    name = "mygene"
    organisms = tuple(TAXID)
    tasks = ("b", "c", "d", "e")
    current = 116

    def scopes(self, line):
        return MYGENE_SCOPES.get(line["task"]) or MYGENE_SCOPES.get(source_kind(line))

    def version_calls(self, organism):
        return [dict(method="GET", url=f"{MYGENE}/metadata", endpoint="metadata")]

    def input_problem(self, line, organism):
        if "," in line["id"]:
            return "a term with a comma cannot be one of the comma-separated terms of a batch"
        if not self.scopes(line):
            return "no scope for this input"
        return ""

    def requests(self, lines, organism):
        groups = {}
        for q in lines:
            groups.setdefault(self.scopes(q), []).append(q)
        out = []
        for sc, qs in sorted(groups.items()):
            for b in chunks(qs, 1000):
                terms = sorted({q["id"] for q in b})
                # the service splits `q` at commas and white space as a shell would: a term with a space or a quote
                # is quoted (shlex), so that it stays one term; hits per term: at most 50 (terms x size <= 100,000)
                out.append(dict(method="POST", url=f"{MYGENE}/query",
                                data={"q": ",".join(shlex.quote(t) if re.search(r"[\s'\"\\]", t) else t for t in terms),
                                      "scopes": sc, "species": str(TAXID[organism]),
                                      "fields": "ensembl.gene,symbol,taxid", "size": "50", "dotfield": "true"},
                                qids=[q["qid"] for q in b], endpoint=f"query:{sc}"))
        return out


HGNC_CHAIN = {"c": ("symbol", "prev_symbol", "alias_symbol"), "d": ("symbol", "prev_symbol", "alias_symbol"),
              "ncbi": ("entrez_id",), "refseq": ("refseq_accession",), "uniprot": ("uniprot_ids",)}


class Hgnc(Tool):
    name = "hgnc"
    organisms = ("homo_sapiens",)
    tasks = ("c", "d", "e")
    current = 116

    def chain(self, line):
        return HGNC_CHAIN.get(line["task"]) or HGNC_CHAIN.get(source_kind(line))

    def version_calls(self, organism):
        return [dict(method="GET", url=f"{HGNC}/info", headers={"Accept": "application/json"}, endpoint="info")]

    def input_problem(self, line, organism):
        if not self.chain(line):
            return "no HGNC field for this input"
        if "/" in line["id"]:
            return "a string with a slash cannot be a path segment of fetch"
        return ""

    def request(self, line, field):
        v = line["id"]
        return dict(method="GET", url=f"{HGNC}/fetch/{field}/{urllib.parse.quote(v, safe='')}",
                    headers={"Accept": "application/json"}, qids=[line["qid"]], endpoint=f"fetch/{field}")


class UniProt(Tool):
    name = "uniprot"
    organisms = ("sus_scrofa", "mus_musculus", "homo_sapiens", "macaca_mulatta", "macaca_fascicularis")
    tasks = ("e",)
    current = 116

    def input_problem(self, line, organism):
        return "" if source_kind(line) == "uniprot" else "UniProt's ID mapping is asked with UniProt accessions only"

    def version_calls(self, organism):
        return [dict(method="GET", url=f"{UNIPROT}/configure/idmapping/fields", headers=JSON, endpoint="idmapping/fields")]


BRIDGEDB_SYSTEM = {"ensembl": "En", "ncbi": "L", "refseq": "Q", "uniprot": "S"}


class BridgeDb(Tool):
    name = "bridgedb"
    organisms = tuple(BRIDGEDB_ORGANISM)
    tasks = ("b", "c", "e")
    current = 111

    def version_calls(self, organism):
        o = urllib.parse.quote(BRIDGEDB_ORGANISM[organism])
        return [dict(method="GET", url=f"{BRIDGEDB}/{o}/properties", ext="txt", endpoint="properties"),
                dict(method="GET", url=f"{BRIDGEDB}/contents", ext="txt", endpoint="contents")]

    def system(self, line, organism):
        if line["task"] == "c":
            return "H" if organism == "homo_sapiens" else None
        return BRIDGEDB_SYSTEM.get(source_kind(line))

    def input_problem(self, line, organism):
        return "" if self.system(line, organism) else "no BridgeDb data source for this input in this organism"

    def requests(self, lines, organism):
        o = urllib.parse.quote(BRIDGEDB_ORGANISM[organism])
        groups = {}
        for q in lines:
            groups.setdefault(self.system(q, organism), []).append(q)
        out = []
        for sys_, qs in sorted(groups.items()):
            for b in chunks(qs, 100):
                ids = [q["id"] for q in b]
                out.append(dict(method="POST", url=f"{BRIDGEDB}/{o}/xrefsBatch/{sys_}", params={"dataSource": "En"},
                                data="\n".join(ids), headers={"Content-Type": "text/plain"}, ext="txt",
                                qids=[q["qid"] for q in b], endpoint=f"xrefsBatch/{sys_}"))
        return out


TOOLS = {t.name: t for t in (EnsemblRest(), BioMart(), GProfiler(), MyGene(), Hgnc(), UniProt(), BridgeDb())}


def gprofiler_archives(listing):
    """{release: base URL} of g:Profiler's archived instances with the current API (`gprofiler_archive3/e<r>_...`),
    from api/util/list_archives (the latest revision of each release)."""
    out = {}
    for a in sorted(json.loads(listing) if isinstance(listing, (str, bytes)) else listing, key=lambda x: x["date"]):
        m = re.search(r"/gprofiler_archive3/e(\d+)_eg\d+_p\d+$", a["url"].rstrip("/"))
        if m:
            out[int(m.group(1))] = a["url"].rstrip("/")
    return out
