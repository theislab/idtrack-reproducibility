"""A0, graph fidelity (docs/methods.md, A): the graph against Ensembl's tables, no search.

Usage: a0_fidelity.py <organism> --clone <idtrack clone> --graph <graph .npz> --out <folder> [--checks P,H,T,N]

Checks (each over every release of the graph's window that the truth has):

- P presence: every GTF gene (id.version; the id alone where the GTF of the release writes no versions) of every
  release is a backbone node active at that release, and no other backbone node is. Read two ways: the history view
  (`active_releases`: the ranges of the history edges) and the cross-reference view (`combined`: the releases at
  which the base id and the transcripts link the gene). The truth's genes of the core database that the GTF lacks
  (human 48-51, cynomolgus 112; A*'s `truth.TruthA`) are tagged `core_only`; the GTF's patch genes `patch`.
- H history: every gene event of the release-115 `stable_id_event` whose session lies in the window (placed where
  its effect appears, `make_queries.Truth.effective_release`) is the edge it should be: between two ids the edge
  old.v -> new.v at the session's two releases; a version change the edge id.v -> id.w; a birth `<id>.Void` (or
  `.Retired`) -> id.v entering at the session's new release; a retirement id.v -> `<id>.Retired` leaving at its old
  release; an id mapped onto itself with the same version an edge only by the builder's rule (score below 1 or the
  same session maps the id to or from another id). And the converse: every history edge of the graph is an event of
  Ensembl's (between ids), a change of version between two consecutive GTFs, or a birth, retirement or return
  that the GTFs show.
- T transcripts and translations: every GTF transcript (id.version) of a release links its gene (id.version) at that
  release, every translation its transcript; and the converse.
- N names and cross-references: every value of a BioMart table of a release in the window (display name, synonym,
  NCBI gene id and label, HGNC, MGI, VGNC symbol, UniProt gene name, Swiss-Prot, TrEMBL, RefSeq, RFAM, WikiGene) is
  linked at that release to the gene (directly, through a transcript or through a translation) under the graph's
  database of that table; and the converse for the tables the graph includes. A table whose database the graph does
  not include is counted as `database_not_in_graph` (the externals file of the organism chooses it).

Outputs in --out: `summary.tsv` (check, class, by_design, count, distinct ids, example), `summary.json`, and one
`<check>.tsv.gz` with every row that is not `ok`. Nothing here reads idtrack's answers; the graph is read through its
accessors (`history_edges`, `xref_edges`, `active_releases`, `combined`), which both api generations have.
"""
import argparse
import collections
import gzip
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M
from idtrack_reproducibility.families.astar import truth as TA  # A*'s truth: core genes the GTF lacks

VOID, RETIRED = "Void", "Retired"

# BioMart table family -> (tables as (group, column), graph databases that carry it; None = any name database)
NAME_FAMILIES = {
    "display": ([("genes", None)], None),
    "synonym": ([("synonym", 1)], "synonym"),
    "ncbi": ([("entrez", 1), ("entrez", 2)], {"EntrezGene", "NCBI gene", "NCBI gene (formerly Entrezgene)", "NCBI gene ID"}),
    "hgnc": ([("hgnc", 1)], {"HGNC Symbol"}),
    "mgi": ([("mgi", 1)], {"MGI Symbol"}),
    "vgnc": ([("vgnc", None)], {"VGNC Symbol", "VGNC"}),
    "uniprot_gn": ([("uniprot_gn", 1)], {"UniProtKB Gene Name"}),
    "swissprot": ([("swissprot", 1)], {"Uniprot/SWISSPROT"}),
    "sptrembl": ([("sptrembl", 1)], {"Uniprot/SPTREMBL"}),
    "refseq_mrna": ([("refseq_mrna", 1)], {"RefSeq_mRNA"}),
    "refseq_mrna_predicted": ([("refseq_mrna_predicted", 1)], {"RefSeq_mRNA_predicted"}),
    "refseq_ncrna": ([("refseq_ncrna", 1)], {"RefSeq_ncRNA"}),
    "refseq_ncrna_predicted": ([("refseq_ncrna_predicted", 1)], {"RefSeq_ncRNA_predicted"}),
    "refseq_peptide": ([("refseq_peptide", 1)], {"RefSeq_peptide"}),
    "refseq_peptide_predicted": ([("refseq_peptide_predicted", 1)], {"RefSeq_peptide_predicted"}),
    "rfam": ([("rfam", 1)], {"RFAM"}),
    "wikigene": ([("wikigene", 1)], {"WikiGene"}),
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def split_node(name):
    i = name.rfind(".")
    return (name[:i], name[i + 1:]) if i > 0 else (name, "")


class Tally:
    """Counts per (check, class) with distinct ids and a few examples; rows that are not `ok` to a file."""

    def __init__(self, out, check, by_design):
        self.check, self.by_design = check, by_design
        self.n = collections.Counter()
        self.ids = collections.defaultdict(set)
        self.ex = collections.defaultdict(list)
        self.fh = gzip.open(os.path.join(out, f"{check}.tsv.gz"), "wt")
        self.header = None

    def add(self, cls, uid, row=None):
        self.n[cls] += 1
        self.ids[cls].add(uid)
        if row is not None and not cls.split(":")[-1].startswith("ok"):
            if self.header is None:
                self.header = ["class"] + list(row)
                self.fh.write("\t".join(self.header) + "\n")
            self.fh.write("\t".join([cls] + [str(row[k]) for k in self.header[1:]]) + "\n")
            if len(self.ex[cls]) < 5:
                self.ex[cls].append(row)

    def close(self):
        self.fh.close()

    def design(self, cls):
        """The reason a class is by design: the class, or its tail after a ':' or a '|' (`node_missing:...`)."""
        base = cls.split("|", 1)[0]
        parts = [cls] + [cls.split(":", i)[-1] for i in range(1, cls.count(":") + 1)] + [cls.rsplit("|", 1)[-1]]
        if base != cls:  # a region suffix (`|patch`, `|not_in_gtf`) does not change why the class is by design
            parts += [base] + [base.split(":", i)[-1] for i in range(1, base.count(":") + 1)]
        for c in parts:
            if c in self.by_design:
                return self.by_design[c]
        return ""

    def rows(self):
        for cls in sorted(self.n):
            yield dict(check=self.check, cls=cls, by_design=self.design(cls), count=self.n[cls],
                       ids=len(self.ids[cls]), example=json.dumps(self.ex[cls][:2], default=str) if self.ex[cls] else "")


class GraphReader:
    def __init__(self, clone, graph, org, workdir):
        from idtrack_reproducibility.families import adapter as AD
        self.ad = AD.Adapter(clone, graph, org, C.release_of_graph(graph), workdir)
        self.G = self.ad.graph
        st = self.G._frozen()
        self.names = st.names
        self.rels = sorted(int(r) for r in self.G.releases) if hasattr(self.G, "releases") else None
        self.backbone, self.tx, self.tl, self.ext = [], [], [], set()
        types = collections.Counter()
        for i, n in enumerate(st.names):
            t = st.node_type(i)
            types[t] += 1
            if t == "ensembl_gene":
                self.backbone.append(n)
            elif t == "ensembl_transcript":
                self.tx.append(n)
            elif t == "ensembl_translation":
                self.tl.append(n)
        self.types = dict(types)
        self.nonref = getattr(self.G, "non_reference_releases", None)
        try:
            self.dbs = sorted(self.G.available_external_databases)
        except Exception:  # noqa: BLE001
            self.dbs = []
        self.manifest = (getattr(self.G, "graph", None) or {}).get("manifest") or {}
        serve = getattr(self.G, "serving_assembly", None) or getattr(self.G, "backbone_assembly", None) or {}
        self.serve = {int(r): int(a) for r, a in dict(serve).items()}
        # per assembly, the mask of the releases it serves: a link counts only on the assembly that serves its release
        self.serve_mask = collections.defaultdict(int)
        for r, asm in self.serve.items():
            self.serve_mask[asm] |= 1 << r

    def served(self, asm, rs):
        """(mask of releases of `rs` on the assembly that serves them, mask of the others)."""
        m = self.mask(rs)
        s = m & self.serve_mask.get(int(asm), 0)
        return s, m & ~s

    def mask(self, rs):
        return rs.mask if hasattr(rs, "mask") else sum(1 << int(r) for r in rs)


def bits(mask):
    r, out = 0, []
    while mask:
        if mask & 1:
            out.append(r)
        mask >>= 1
        r += 1
    return out


# ---------------------------------------------------------------- P presence

def check_presence(gr, tr, rels, out):
    G = gr.G
    by = {"core_only_absent": "truth: core gene the GTF lacks (A*'s truth adds it)",
          "graph_extra|lrg": "truth: the GTF lacks LRG records (hh-5); the core has them",
          "nonref_mark_extra|gtf_without_patch_file": "truth: the GTF of the release has no patch file, so every gene reads primary (J3)"}
    t = Tally(out, "presence", by)
    hist = collections.defaultdict(lambda: collections.defaultdict(set))  # R -> id -> {versions} (history view)
    xref = collections.defaultdict(lambda: collections.defaultdict(set))  # R -> id -> {versions} (xref view)
    nonref = collections.defaultdict(set)  # R -> ids marked non-reference
    relmask = sum(1 << r for r in rels)
    errors = collections.Counter()
    for k, node in enumerate(gr.backbone):
        gid, v = split_node(node)
        if v in (VOID, RETIRED):
            continue
        try:
            m = G.active_releases(node).mask & relmask
        except Exception as e:  # noqa: BLE001
            errors[type(e).__name__] += 1
            t.add("graph_ranges_error", gid, dict(release="", id=gid, version=v, region="", note=str(e)[:200]))
            m = 0
        for r in bits(m):
            hist[r][gid].add(v)
        try:
            xm = 0
            for per in G.combined(node).values():
                for asm, rs in per.items():
                    xm |= gr.served(asm, rs)[0]
            xm &= relmask
        except KeyError:
            xm = 0
        for r in bits(xm):
            xref[r][gid].add(v)
        if gr.nonref is not None:
            for r in bits(gr.mask(gr.nonref(node)) & relmask):
                nonref[r].add(gid)
        if k % 200000 == 0:
            log("presence: graph nodes", k)
    patch_file = {r for r in rels if any(g["region"] == "patch" for g in (tr.genes.get(r) or {}).values())}
    for r in rels:
        gtf = tr.genes.get(r) or {}
        versioned = tr.cmp_for(r) == "versioned"
        for gid, g in gtf.items():
            v = g["gene_version"] if versioned else ""
            region = g["region"]
            hv, xv = hist[r].get(gid, set()), xref[r].get(gid, set())
            row = dict(release=r, id=gid, version=v, region=region, graph_history=",".join(sorted(hv)),
                       graph_xref=",".join(sorted(xv)))
            inh = (v in hv) if versioned else bool(hv)
            inx = (v in xv) if versioned else bool(xv)
            if inh and inx:
                cls = "ok"
            elif inh:
                cls = "ok_history_not_xref"
            elif inx:
                cls = "xref_not_history"
            elif hv or xv:
                cls = "other_version"
            elif gid in tr.present and region == "core_only":
                cls = "core_only_absent"
            else:
                cls = "absent"
            if region != "primary":
                cls += f"|{region}"
            t.add(cls, gid, row)
            if gr.nonref is not None and region in ("primary", "patch"):
                marked = gid in nonref[r]
                if (region == "patch") != marked:
                    t.add(f"nonref_mark_{'missing' if region == 'patch' else 'extra'}" +
                          ("" if r in patch_file or region == "patch" else "|gtf_without_patch_file"), gid,
                          dict(release=r, id=gid, version=v, region=region, graph_history="", graph_xref=""))
        for gid, hv in hist[r].items():
            g = gtf.get(gid)
            if g is None:
                t.add("graph_extra" + ("|lrg" if gid.startswith("LRG_") else ""), gid, dict(release=r, id=gid, version=",".join(sorted(hv)), region="",
                                               graph_history=",".join(sorted(hv)), graph_xref=""))
            elif versioned and len(hv) > 1:
                t.add("graph_two_versions", gid, dict(release=r, id=gid, version=g["gene_version"], region=g["region"],
                                                      graph_history=",".join(sorted(hv)), graph_xref=""))
        log("presence release", r, "gtf", len(gtf), "graph", len(hist[r]))
    t.close()
    return t, hist


# ---------------------------------------------------------------- H history

def builder_self_rule(ev):
    """{(session, id, version)} of ids that a session also maps to or from another id (the builder's rule)."""
    beside = set()
    for o, ov, n, nv, s, sc in ev:
        if o and n and o != n:
            beside.add((s, o, ov))
            beside.add((s, n, nv))
    return beside


def check_history(gr, tr, rels, out, hist=None):
    G = gr.G
    hist = hist or {}

    def tables_have(node, r):
        """The id of `node` is in the graph's tables at release r with another version (Ensembl's event names a
        version its gene tables do not have there: rule `event_version`'s kind)."""
        i, v = split_node(node)
        vs = (hist.get(r) or {}).get(i, set())
        return bool(vs) and v not in vs
    by = {
        "self_same_version_no_edge_by_rule": "builder's rule: an id mapped onto itself with score 1 and no other event "
                                              "of the id in the session is no edge",
        "outside_window": "session outside the graph's window",
        "version_zero": "Ensembl's version 0 (builder leaves it out)",
        "birth:node_missing:id_in_no_gtf": "builder: an event whose id is in no id table is left out",
        "retirement:node_missing:id_in_no_gtf": "builder: an event whose id is in no id table is left out",
        "between_ids:node_missing:id_in_no_gtf": "builder: an event whose id is in no id table is left out",
        "edge_version_change_unversioned_gtf": "the GTFs of the two releases write no versions: not checkable",
        "birth:other_release": "Ensembl's event contradicts its gene tables (the id exists earlier); the builder follows the tables",
        "version_change:other_pair": "the builder's documented rule: an event that does not fit its session's pair is placed at the nearest pair it fits",
        "between_ids:other_pair": "the builder's documented rule: an event that does not fit its session's pair is placed at the nearest pair it fits",
        "event_version_not_in_tables": "Ensembl's event names a version its gene tables do not have there (rule `event_version`'s kind); the builder follows the tables",
        "lrg": "truth: the GTF lacks LRG records (hh-5); the core has them",
    }
    t = Tally(out, "history", by)
    first, last = rels[0], rels[-1]
    pos = {r: i for i, r in enumerate(rels)}
    edge_cache = {}

    def out_edges(node):
        if node not in edge_cache:
            try:
                edge_cache[node] = [(e.node, e.old_release, e.new_release) for e in G.history_edges(node)]
            except KeyError:
                edge_cache[node] = None
        return edge_cache[node]

    index = G._frozen().index
    exists = lambda x: x in index  # noqa: E731
    core, ev, sessions = H.core_tables(tr.org, tr.T)
    beside = builder_self_rule(ev)
    event_keys = collections.defaultdict(set)  # (o.ov, n.nv) -> {(before, after)}
    for o, ov, n, nv, s, sc in ev:
        before, after = tr.sides(s)
        e_rel = tr.eff[s]
        row = dict(session=s, old=f"{o}.{ov}" if o else "", new=f"{n}.{nv}" if n else "", score=sc,
                   recorded=f"{sessions[s]['old_release']}->{sessions[s]['new_release']}", sides=f"{before}->{after}",
                   graph="")
        uid = o or n
        if before is None or after is None or e_rel <= first or e_rel > last or before < first:
            t.add("outside_window", uid)
            continue
        if (o and str(ov) == "0") or (n and str(nv) == "0"):
            t.add("version_zero", uid, row)
            continue
        if o and n:
            on, nn = f"{o}.{ov}", f"{n}.{nv}"
            event_keys[(on, nn)].add((before, after))
            kind = "between_ids" if o != n else ("version_change" if ov != nv else "self_same_version")
            if kind == "self_same_version" and not (sc < 1.0 or (s, o, ov) in beside):
                # the builder's rule: no edge expected
                es = out_edges(on) or []
                has = any(x == nn and a == before and b == after for x, a, b in es)
                t.add("self_same_version_edge_without_rule" if has else "self_same_version_no_edge_by_rule", uid,
                      row if has else None)
                continue
            if not exists(on) or not exists(nn):
                miss = [x for x in (on, nn) if not exists(x)]
                g_has = [x for x in miss if tr.present.get(split_node(x)[0])]
                other_v = (not exists(on) and tables_have(on, before)) or (not exists(nn) and tables_have(nn, after))
                cls = f"{kind}:node_missing" + ("" if g_has else ":id_in_no_gtf") + (":event_version_not_in_tables" if other_v and g_has else "")
                row["graph"] = "missing " + ",".join(miss)
                t.add(cls, uid, row)
                continue
            es = out_edges(on) or []
            pairs = [(a, b) for x, a, b in es if x == nn]
            if (before, after) in pairs:
                t.add(f"ok_{kind}", uid)
            elif pairs:
                row["graph"] = ";".join(f"{a}->{b}" for a, b in pairs)
                t.add(f"{kind}:other_pair", uid, row)
            else:
                t.add(f"{kind}:no_edge", uid, row)
        elif n:  # birth
            nn = f"{n}.{nv}"
            if not exists(nn):
                t.add("birth:node_missing" + ("" if tr.present.get(n) else ":id_in_no_gtf") +
                      (":event_version_not_in_tables" if tr.present.get(n) and tables_have(nn, after) else ""), uid, row)
                continue
            try:
                ins = [(e.node, e.old_release, e.new_release) for e in G.history_edges(nn, reverse=True)]
            except KeyError:
                ins = []
            origin = [(x, a, b) for x, a, b in ins if split_node(x)[1] in (VOID, RETIRED) and split_node(x)[0] == n]
            if any(b == after for x, a, b in origin):
                t.add("ok_birth", uid)
            elif origin:
                row["graph"] = ";".join(f"{x}:{a}->{b}" for x, a, b in origin)
                t.add("birth:other_release", uid, row)
            else:
                row["graph"] = ";".join(f"{x}:{a}->{b}" for x, a, b in ins)
                t.add("birth:no_origin_edge", uid, row)
        elif o:  # retirement
            on = f"{o}.{ov}"
            if not exists(on):
                t.add("retirement:node_missing" + ("" if tr.present.get(o) else ":id_in_no_gtf") +
                      (":event_version_not_in_tables" if tr.present.get(o) and tables_have(on, before) else ""), uid, row)
                continue
            es = out_edges(on) or []
            ret = [(x, a, b) for x, a, b in es if x == f"{o}.{RETIRED}"]
            if any(a == before for x, a, b in ret):
                t.add("ok_retirement", uid)
            elif ret:
                row["graph"] = ";".join(f"{a}->{b}" for x, a, b in ret)
                t.add("retirement:other_release", uid, row)
            else:
                row["graph"] = ";".join(f"{x}:{a}->{b}" for x, a, b in es)
                t.add("retirement:no_retired_edge", uid, row)
    log("history: events read", sum(t.n.values()))
    # the converse: every history edge of the graph
    for k, node in enumerate(gr.backbone):
        es = out_edges(node) or []
        gid, v = split_node(node)
        for x, a, b in es:
            xid, xv = split_node(x)
            row = dict(session="", old=node, new=x, score="", recorded="", sides=f"{a}->{b}", graph="edge")
            if b == math.inf or (isinstance(b, float) and math.isinf(b)):
                cls = "edge_ok_last_loop" if tr.ver(gid, last) == v or (tr.cmp_for(last) == "base" and last in tr.present.get(gid, ())) else "edge_last_loop_not_in_gtf"
            elif v in (VOID, RETIRED):
                pres = tr.present.get(gid, ())
                vok = tr.cmp_for(b) != "versioned" or tr.ver(gid, b) == xv
                if b in pres and (a not in pres or a == b) and vok:
                    cls = "edge_ok_birth_or_return"
                elif b in pres and (a not in pres or a == b):
                    cls = "edge_birth_other_version_in_gtf"
                else:
                    cls = "edge_birth_not_in_gtf"
            elif xv == RETIRED:
                pres = tr.present.get(gid, ())
                vok = tr.cmp_for(a) != "versioned" or tr.ver(gid, a) == v
                cls = ("edge_ok_retirement" if vok else "edge_retirement_other_version_in_gtf") if a in pres and b not in pres else "edge_retirement_not_in_gtf"
            elif xid != gid:
                pairs = event_keys.get((node, x), set())
                cls = "edge_ok_event" if (a, b) in pairs else ("edge_event_other_pair" if pairs else "edge_without_event")
            elif xv != v:
                vf = tr.cmp_for(a) == "versioned" and tr.cmp_for(b) == "versioned"
                ok = (tr.ver(gid, a) == v and tr.ver(gid, b) == xv) if vf else (a in tr.present.get(gid, ()) and b in tr.present.get(gid, ()))
                cls = "edge_ok_version_change" if ok else ("edge_version_change_not_in_gtf" if vf else "edge_version_change_unversioned_gtf")
            else:
                cls = "edge_ok_self_event" if (node, x) in event_keys else "edge_self_without_event"
            if gid.startswith("LRG_") and "not_in_gtf" in cls:
                cls += "|lrg"
            t.add(cls, gid, row)
        if k % 200000 == 0:
            log("history converse: nodes", k)
    t.close()
    return t


# ---------------------------------------------------------------- T transcripts

def check_transcripts(gr, tr, rels, out):
    G = gr.G
    t = Tally(out, "transcripts", {"tx:link_only_on_another_assembly": "a link held only on an assembly that does not serve the release (counted per id)",
                                   "tl:link_only_on_another_assembly": "a link held only on an assembly that does not serve the release (counted per id)","tx:graph_extra|lrg": "truth: the GTF lacks LRG records (hh-5)",
                                   "tl:graph_extra|lrg": "truth: the GTF lacks LRG records (hh-5)",
                                   "tx:graph_extra|parent_core_only": "truth: transcripts of core genes the GTF lacks (hh-5, J2)",
                                   "tl:graph_extra|parent_core_only": "truth: translations of core genes the GTF lacks (hh-5, J2)"})
    index = G._frozen().index
    relmask = sum(1 << r for r in rels)
    # graph: transcript base -> [(transcript node, gene node, mask)], translation base -> [(node, transcript node, mask)]
    links = {}
    for kind, nodes, target_type in (("tx", gr.tx, "ensembl_gene"), ("tl", gr.tl, "ensembl_transcript")):
        d = collections.defaultdict(list)
        for node in nodes:
            for tgt, con in G.xref_edges(node):
                m = other = 0
                for db, asm, rs in con.items:
                    a_, b_ = gr.served(asm, rs)
                    m |= a_
                    other |= b_
                d[split_node(node)[0]].append((node, tgt, m & relmask))
                if other & relmask & ~m:
                    t.add(f"{kind}:link_only_on_another_assembly", split_node(node)[0])
        links[kind] = d
        log("transcripts: graph", kind, len(d))
    seen = collections.defaultdict(set)  # kind -> {(base, R)} found in the GTF
    for r in rels:
        rows = H.gtf_tx(tr.org, r) or []
        versioned = tr.cmp_for(r) == "versioned"
        for x in rows:
            gene = x["gene_id"]
            gv = tr.ver(gene, r) if versioned else ""
            for kind, sid, sv, parent, pv in (("tx", x["transcript_id"], x.get("transcript_version", ""), gene, gv),
                                              ("tl", x.get("protein_id", ""), x.get("protein_version", ""),
                                               x["transcript_id"], x.get("transcript_version", ""))):
                if not sid:
                    continue
                seen[kind].add((sid, r))
                cand = links[kind].get(sid, [])
                at = [(n, tg) for n, tg, m in cand if m >> r & 1]
                row = dict(kind=kind, release=r, id=f"{sid}.{sv}" if sv else sid, parent=f"{parent}.{pv}" if pv else parent,
                           region=x.get("region", ""), graph=";".join(f"{n}->{tg}" for n, tg in at)[:300])
                if not cand:
                    cls = "absent"
                elif not at:
                    cls = "not_linked_at_release"
                else:
                    def okv(node, want):
                        return not want or node == f"{split_node(node)[0]}.{want}"
                    good = [(n, tg) for n, tg in at if split_node(tg)[0] == parent and (not sv or okv(n, sv)) and (not pv or okv(tg, pv))]
                    if good:
                        cls = "ok"
                    elif any(split_node(tg)[0] == parent for n, tg in at):
                        cls = "other_version"
                    else:
                        cls = "other_parent"
                if x.get("region", "primary") != "primary":
                    cls += f"|{x['region']}"
                t.add(f"{kind}:{cls}", sid, row)
        log("transcripts release", r, len(rows))
    for kind in ("tx", "tl"):
        for sid, cand in links[kind].items():
            for n, tg, m in cand:
                for r in bits(m):
                    if (sid, r) not in seen[kind]:
                        pg = split_node(tg)[0]
                        core_only = ((tr.genes.get(r) or {}).get(pg) or {}).get("region") == "core_only" if kind == "tx" else \
                            any(((tr.genes.get(r) or {}).get(split_node(t2)[0]) or {}).get("region") == "core_only"
                                for n2, t2, m2 in links["tx"].get(pg, []) if m2 >> r & 1)
                        t.add(f"{kind}:graph_extra" + ("|parent_core_only" if core_only else "") + ("|lrg" if sid.startswith("LRG_") else ""), sid, dict(kind=kind, release=r, id=n, parent=tg, region="", graph="edge"))
    t.close()
    return t


# ---------------------------------------------------------------- N names and cross-references

def truth_pairs(org, r, fam, meta=None):
    """gene base id -> {values} of one table family at BioMart release r; None when no table. `meta` (a dict) gets
    the display name's source per (gene, value) and the table columns that are absent at the release."""
    tables, _ = NAME_FAMILIES[fam]
    out, found = collections.defaultdict(set), False
    meta = meta if meta is not None else {}
    meta.setdefault("source", {})
    meta.setdefault("absent", set())
    for grp, col in tables:
        if grp == "genes":
            rows = H.biomart_genes(org, r)
            if rows is None:
                continue
            found = True
            for x in rows:
                if x["name"]:
                    out[x["gene_id"]].add(x["name"])
                    meta["source"][(x["gene_id"], x["name"])] = x["source"]
        elif grp == "vgnc":
            pp = H.vgnc_symbols(org, r)
            if pp is None:
                continue
            found = True
            for g, vs in pp[0].items():
                out[g] |= vs
        else:
            g2 = grp
            if grp == "swissprot" and os.path.exists(os.path.join(H.truth_dir(org), "biomart", str(r), "swissprot_acc.tsv")):
                g2 = "swissprot_acc"  # the corrected judge: BioMart 75 writes entry names; the accessions are this table
            t = H.biomart(org, r, g2)
            if t is not None and len(t[0]) <= col:
                meta["absent"].add((grp, col))
                continue
            pp = H.biomart_pairs(org, r, g2, col)
            if pp is None:
                continue
            found = True
            for g, vs in pp[0].items():
                out[g] |= vs
    return out if found else None


def check_names(gr, tr, rels, hist, out):
    G = gr.G
    by = {"database_not_in_graph": "the organism's externals file does not include the database",
          "source_not_in_graph": "the display name's source database is not in the graph (miRBase: rule `mirbase_genes` adds it)",
          "graph_extra|label_column_absent": "truth: BioMart's NCBI table of the release has no label column",
          "vgnc:missing": "authority: VGNC's own current and previous symbols of the gene's VGNC id, not Ensembl's label then",
          "vgnc:graph_extra": "authority: VGNC's own current and previous symbols of the gene's VGNC id, not Ensembl's label then",
          "vgnc:under_another_database": "authority: VGNC's own symbols; Ensembl wrote the gene's label under another database",
          "vgnc:graph_extra_gene_without_table_row": "authority: VGNC's own symbols (no VGNC id in BioMart for the gene)",
          "links_on_another_assembly_not_read": "links of a release on an assembly that does not serve it (human GRCh37 at 76 on): not part of the graph's answers at that release; counted, not compared",
          "graph_extra_gene_without_table_row|lrg": "truth: BioMart has no rows for LRG records (hh-5)",
          "other_spelling": "rule `spelling_steps`: a node is named as the value was written at its first release",
          "rfam:graph_extra": "rule `rfam_display_label`: the graph holds an RFAM family's accession and its label as two ids"}
    t = Tally(out, "names", by)
    bm = [r for r in H.biomart_releases(tr.org) if r in rels]
    bmask = sum(1 << r for r in bm)
    types = {}
    store = G._frozen()

    def ntype(x):
        if x not in types:
            types[x] = store.node_type(store.index[x])
        return types[x]

    # graph: (release, gene base) -> {(database, value)}
    graph_vals = collections.defaultdict(set)
    other_assembly = [0]  # (value, gene, release) links on an assembly that does not serve the release (not read)

    def collect(node, mask, gbase, depth):
        for src, con in G.xref_edges(node, reverse=True):
            st = ntype(src)
            for db, asm, rs in con.items:
                m, other = gr.served(asm, rs)
                m &= mask
                if other & mask and st == "external":
                    other_assembly[0] += bin(other & mask).count("1")
                if not m:
                    continue
                if st in ("ensembl_transcript", "ensembl_translation") and depth < 2:
                    collect(src, m, gbase, depth + 1)
                elif st == "external":
                    for r in bits(m):
                        graph_vals[(r, gbase)].add((db, src))

    for k, node in enumerate(gr.backbone):
        gid, v = split_node(node)
        if v in (VOID, RETIRED):
            continue
        try:
            m = G.active_releases(node).mask & bmask
        except Exception:  # noqa: BLE001
            continue
        if m:
            collect(node, m, gid, 0)
        if k % 200000 == 0:
            log("names: graph nodes", k)
    dbs = set(gr.dbs)
    log("names: graph databases", sorted(dbs))
    t.n["links_on_another_assembly_not_read"] += other_assembly[0]
    t.ids["links_on_another_assembly_not_read"]
    for r in bm:
        gtf = tr.genes.get(r) or {}
        for fam, (tables, gdbs) in NAME_FAMILIES.items():
            meta = {}
            tp = truth_pairs(tr.org, r, fam, meta)
            if tp is None:
                continue
            if gdbs is None:
                match_db = lambda db: not db.startswith("synonym_id::")  # noqa: E731
                in_graph = True
            elif gdbs == "synonym":
                match_db = lambda db: db.startswith("synonym_id::")  # noqa: E731
                in_graph = any(d.startswith("synonym_id::") for d in dbs)
            else:
                match_db = lambda db, s=gdbs: db in s  # noqa: E731
                in_graph = bool(gdbs & dbs)
            genes_seen = set()
            for g, vals in tp.items():
                genes_seen.add(g)
                region = (gtf.get(g) or {}).get("region", "not_in_gtf")
                gv = {(db, val) for db, val in graph_vals.get((r, g), ())}
                gvals = {val[len("synonym_id::"):] if val.startswith("synonym_id::") else val for db, val in gv if match_db(db)}
                gany = {val for db, val in gv}
                glower = {x.lower() for x in gvals}
                for val in vals:
                    row = dict(release=r, family=fam, gene=g, value=val, region=region, graph_dbs="")
                    src = meta["source"].get((g, val))
                    if src is not None:
                        row["graph_dbs"] = f"source={src}"
                    if not in_graph:
                        cls = "database_not_in_graph"
                    elif val in gvals:
                        cls = "ok"
                    elif src and src not in dbs and fam == "display":
                        cls = "source_not_in_graph"
                    elif g not in hist[r]:
                        cls = "gene_not_in_graph_at_release"
                    elif val.lower() in glower:
                        cls = "other_spelling"
                    elif val in gany or f"synonym_id::{val}" in gany:
                        cls = "under_another_database"
                        row["graph_dbs"] = ",".join(sorted(db for db, x in gv if x in (val, f"synonym_id::{val}")))
                    else:
                        cls = "missing"
                    if region not in ("primary",) and cls != "ok":
                        cls += f"|{region}"
                    t.add(f"{fam}:{cls}", f"{g}|{val}", row)
                if in_graph and gdbs not in (None,):
                    tvals = set(vals)
                    tl = {x.lower() for x in tvals}
                    for val in gvals - tvals:
                        if val.lower() in tl:
                            continue
                        lab = "|label_column_absent" if fam == "ncbi" and ("entrez", 2) in meta["absent"] and not val.isdigit() else ""
                        t.add(f"{fam}:graph_extra{lab}", f"{g}|{val}", dict(release=r, family=fam, gene=g, value=val,
                                                                       region=region, graph_dbs=""))
            if in_graph and gdbs not in (None,):  # genes with graph values for the family but no truth row
                for (rr, g), gv in graph_vals.items():
                    if rr != r or g in genes_seen:
                        continue
                    for db, val in gv:
                        if match_db(db):
                            t.add(f"{fam}:graph_extra_gene_without_table_row" + ("|lrg" if g.startswith("LRG_") else ""), f"{g}|{val}",
                                  dict(release=r, family=fam, gene=g, value=val, region=(gtf.get(g) or {}).get("region", "not_in_gtf"), graph_dbs=db))
            log("names release", r, fam, len(tp))
    t.close()
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--checks", default="P,H,T,N")
    a = ap.parse_args()
    org = a.organism
    os.makedirs(a.out, exist_ok=True)
    checks = a.checks.split(",")
    log("truth", org)
    tr = TA.TruthA(org, C.release_of_graph(a.graph))
    log("graph", a.graph)
    gr = GraphReader(a.clone, a.graph, org, a.out)
    grels = gr.rels or tr.rels
    rels = [r for r in grels if r in tr.genes]
    log("graph releases", grels[0], grels[-1], len(grels), "with truth", len(rels), "types", gr.types)
    tallies = []
    hist = None
    if "P" in checks or "N" in checks or "H" in checks:
        tp, hist = check_presence(gr, tr, rels, a.out)
        if "P" in checks:
            tallies.append(tp)
    if "H" in checks:
        tallies.append(check_history(gr, tr, rels, a.out, hist))
    if "T" in checks:
        tallies.append(check_transcripts(gr, tr, rels, a.out))
    if "N" in checks:
        tallies.append(check_names(gr, tr, rels, hist, a.out))
    rows = [r for t in tallies for r in t.rows()]
    with open(os.path.join(a.out, "summary.tsv"), "w") as fh:
        fh.write("check\tclass\tby_design\tcount\tids\texample\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in ("check", "cls", "by_design", "count", "ids", "example")) + "\n")
    commit, dirty = C.git_commit(a.clone)
    hist_report = gr.manifest.get("history") if isinstance(gr.manifest, dict) else None
    json.dump(dict(organism=org, clone=a.clone, commit=commit, dirty=dirty, graph=a.graph, load=gr.ad.info,
                   graph_releases=[grels[0], grels[-1]], releases_with_truth=len(rels), node_types=gr.types,
                   databases=gr.dbs, sources=["gtf/genes_<r>.tsv.gz, gtf/tx_<r>.tsv.gz", f"core/{tr.core}",
                                              "biomart/<r>/*.tsv", "vgnc/", "core/*gene_stable_id* (TruthA)"],
                   a0_sha256=L_sha(__file__), astar_truth_sha256=L_sha(TA.__file__),
                   make_queries_sha256=L_sha(M.__file__),
                   manifest_history=str(hist_report)[:20000] if hist_report else None,
                   counts={f"{r['check']}:{r['cls']}": r["count"] for r in rows}),
              open(os.path.join(a.out, "summary.json"), "w"), indent=1, default=str)
    log("done")


def L_sha(p):
    import hashlib
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


if __name__ == "__main__":
    main()
