"""A cause for every wrong answer (F0): the residual reading's classifier as a library every family extends.

Usage: causes.py --organism <o> --scored <scored folder> [--explain <explain .jsonl.gz>] [--und]

The classifier of the reading of every wrong answer of the reference tree's harness run (docs/families/f0.md,
"Causes"): its tests are functions over any family's scored folder
(`wrong.tsv.gz` or, with --und, `undetermined.tsv.gz`) and the explain records of `explain.py` (default
`<scored>/explain_wrong.jsonl.gz`, or `explain_undetermined.jsonl.gz`). Every cause is a test
`t_<cause>(ctx) -> detail or None`; the first test that holds names the cause; a row no test explains is
`unexplained` and is read by hand. Classes (the reading's): a = idtrack is wrong under its own decisions; b = idtrack
follows a decision that the truth contradicts; c = the graph lacks what the truth has; d = the truth is wrong or the
query invalid; e = the truth cannot decide.

Families extend it: `register_kinds(history=..., names=...)` for their kinds (a kind not registered is read by the
shape of its line: a `final_database` -> a database query, an Ensembl id -> history, else a name), and
`register(cause, cls, decision, text, test, before=None)` for a new cause (its test placed before the cause named,
else before the generic ones at the end).

Writes `<scored>/causes.tsv.gz` (organism, qid, view, kind, stratum, unit, cause, class, detail) and
`<scored>/causes_summary.tsv` (stratum, cause, class, rule, answers, units); `--und` writes `causes_und.*`.
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import pickle
import re
import sys

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as MQ

csv.field_size_limit(1 << 30)

# cause -> (class, decision or finding, one line)
CAUSE_CLASS = {
    "entrezgene_label_rows": ("a", "rule `one_kind_per_database`",
                              "EntrezGene up to release 87 holds NCBI labels beside the numbers: the label comes back "
                              "as an NCBI gene id and ties with the number"),
    "biomart75_swissprot_entry_name": ("d", "-", "BioMart 75's attribute `uniprot_swissprot` is the entry name "
                                       "(display label); the accession is `uniprot_swissprot_accession`"),
    "one_spelling_per_db_release": ("b", "rules `spelling_steps` and `spelling_of_release` (9d point 3)",
                                    "one database writes two spellings of one value at one release; the graph keeps one"),
    "mirbase_not_in_graph": ("c", "rule `mirbase_genes`",
                             "miRBase names are not in the graph (externals YAML: miRBase `Include: false`)"),
    "strict_walk_drops_living_id": ("d", "rules `id_lives_on` and `own_id_first`",
                                    "the strict walk leaves an id that exists on both sides of a session when Ensembl "
                                    "recorded only events of other ids into or out of it"),
    "old_name_split_judged_as_set": ("d", "-", "one gene carried the name at R and Ensembl split it; the line wants "
                                     "every successor (`set_exact`)"),
    "support_ranks_shared_name": ("b", "rules `ties`, `from_database`, `support_counts_all`, `all_order` and `support_ranks_shared_name`",
                                  "a display name of several genes: `best` returns the gene(s) that more databases "
                                  "name so; the truth wants the whole set (`set_exact`)"),
    "rfam_family_label": ("b", "rule `one_kind_per_database` renames, does not remove; rule `rfam_display_label`",
                          "an RFAM family label attaches the string to genes whose name it is not"),
    "vgnc_symbol_not_in_truth": ("d", "-", "the gene carries the string as its VGNC symbol; BioMart's VGNC column is "
                                 "the VGNC id, so the truth never sees it"),
    "synonym_spelling_merged": ("c", "rules `spelling_steps` and `spelling_of_release` (9d point 3)",
                                "two genes carry the string in two letter cases; the graph has one node for every case "
                                "of a string and one spelling per database and release"),
    "name_in_other_source": ("e", "rule `names_before_synonyms`",
                             "the returned gene carries the string at the release read in a name source other than "
                             "the display name; the truth counts only the display name"),
    "variant_lacks_broad": ("e", "rule `names_before_synonyms`",
                            "a case or separator variant of a name that another gene carries in another name source"),
    "copies_apart_before_82": ("d", "rule `judge_copies_apart` (its limit)",
                               "the missing genes are copies on patches or haplotypes that the answer reports apart; "
                               "the truth cannot mark copies at target releases before 82"),
    "gtf_lacks_core_gene": ("d", "-", "the returned gene is in the core `gene` table of the target release but not in "
                            "that release's GTF, which the truth reads as the list of genes"),
    "path_through_resurrection": ("e", "graph builder: `<id>.Retired -> <id>.<v>` when an id comes back",
                                  "the answer follows an id that Ensembl retired and later brought back under the "
                                  "same stable id; Ensembl's events end the line at the retirement, the graph continues"),
    "event_version_contradicts_tables": ("c", "rule `event_version`",
                                         "Ensembl's event names an old version the id did not have at the session; "
                                         "the builder leaves the event out, the truth reads it by stable id"),
    "mouse_hgnc_projection": ("d", "harness_lib.EXCLUDED (mouse HGNC Symbol)",
                              "the mouse gene's display name is a human HGNC symbol projected onto it"),
    "name_in_db_not_in_truth": ("d", "-", "the returned gene carries the string in a database whose names the truth "
                                "does not read (Vega gene, Havana gene, VGNC symbol)"),
    "external_bridge": ("b", "rule `external_bridge`",
                        "Ensembl's history has no successor; idtrack bridges through an external id to a gene that "
                        "does not overlap the query's gene"),
    "newer_version_kind_reversed": ("d", "rule `version_newer_than_graph`", "the harness's kind `nonexistent_version` (a current id with its "
                                 "largest version + 37, expecting nothing): rule `version_newer_than_graph` answers a version above every "
                                 "version of a living id with its gene, flagged; the kind is kept and left out of the acceptance"),
    "history_branch": ("e", "-", "(undetermined only) the answer follows Ensembl's history, a branch outside the strict walk"),
    "rule_decided_after_tree": ("a", "rule `judge_follows_rules` (a rule idtrack decided after this tree)",
                                "the line follows a rule decided for idtrack after this tree (the harness's judge `rules`); "
                                "the answer is the one the judge before the rule expected"),
    "unexplained": ("?", "-", "no test holds"),
}

MIR = re.compile(r"^[a-z]{3}[-_](mir|let|lin)[-_]", re.I)
HIST_KINDS = {"retired_renamed", "retired_merged", "retired_split", "retired_none", "retired_none_branch",
              "unchanged_old_to_current", "reverse_present", "reverse_absent", "cross_assembly", "travel_born_late",
              "travel_handoff", "travel_merged", "travel_renamed", "travel_retired_none", "travel_split", "travel_unchanged"}
NAME_KINDS = {"name_unique_to_gene", "name_ambiguous_to_set", "synonym_to_gene", "case_variant", "separator_variant",
              "name_old_to_current", "name_changed_gene", "ncbi_to_gene", "travel_name"}
DB_KINDS = {"gene_to_db", "old_gene_to_db"}
NAME_SRC = [("hgnc", 1), ("mgi", 1), ("uniprot_gn", 1), ("entrez", 2), ("wikigene", 1)]


def register_kinds(history=(), names=(), databases=()):
    HIST_KINDS.update(history)
    NAME_KINDS.update(names)
    DB_KINDS.update(databases)


def base(x):
    return C.base_id(x)


def split(s):
    return [x for x in (s or "").split(",") if x]


# ---------------- the world of one organism: the truth's readers, cached ----------------

class World:
    def __init__(self, org):
        self.org = org
        self._tr = None
        self._bm, self._disp, self._core, self._events = {}, {}, {}, None

    @property
    def tr(self):
        if self._tr is None:
            self._tr = truth(self.org)
        return self._tr

    def bm_pairs(self, rel, group, col):
        k = (rel, group, col)
        if k not in self._bm:
            self._bm[k] = H.biomart_pairs(self.org, rel, group, col)
        return self._bm[k]

    def display(self, rel):
        """gene -> (display name, source) at a BioMart release."""
        if rel not in self._disp:
            self._disp[rel] = {r["gene_id"]: (r["name"], r["source"]) for r in (H.biomart_genes(self.org, rel) or [])}
        return self._disp[rel]

    def name_sources(self, rel, gene, string):
        out = []
        d = self.display(rel).get(gene)
        if d and d[0] == string:
            out.append("display:" + d[1])
        for grp, col in NAME_SRC + [("synonym", 1)]:
            pp = self.bm_pairs(rel, grp, col)
            if pp and string in pp[0].get(gene, ()):
                out.append(grp)
        return out

    def bm_rel(self, rel):
        return rel if rel in H.biomart_releases(self.org) else None

    def core_ids(self, rel):
        """Gene stable ids of the core `gene` table of a release (the tables idtrack downloaded to build the reference
        graph, `<graphs>/<organism>/idtrack_tables/`, read only)."""
        if rel in self._core:
            return self._core[rel]
        import h5py
        ids = None
        for f in glob.glob(os.path.join(CFG.get().graphs, self.org, "idtrack_tables", "*",
                                        f"ens{rel}_mysql_gene_stable_id__cols=gene_id+stable_id+version.h5")):
            with h5py.File(f, "r", locking=False) as h:
                blob = h["c1.blob"][()].tobytes()
                off = h["c1.offsets"][()]
                ids = (ids or set()) | {blob[off[i]:off[i + 1]].decode() for i in range(len(off) - 1)}
        self._core[rel] = ids
        return ids

    def events(self):
        if self._events is None:
            self._events = collections.defaultdict(list)
            for o, ov, n, nv, s, sc in H.core_tables(self.org, 115)[1]:
                self._events[o].append((ov, n, nv, s))
        return self._events


def truth(org):
    """make_queries.Truth of the organism, cached in `<cache>/truth/` under the sha256 of make_queries.py (a change of
    the harness's reading of the history makes a new cache)."""
    import hashlib
    key = hashlib.sha256(open(MQ.__file__, "rb").read()).hexdigest()[:12]
    mine = os.path.join(C.CACHE, "truth", f"truth_{org}_{key}.pickle")
    if os.path.exists(mine):
        with open(mine, "rb") as fh:
            tr = pickle.load(fh)
        tr.succ = collections.defaultdict(dict, tr.succ)
        tr.pred = collections.defaultdict(dict, tr.pred)
        return tr
    tr = MQ.Truth(org, 115)
    os.makedirs(os.path.dirname(mine), exist_ok=True)
    tr.succ = {s: dict(v) for s, v in tr.succ.items()}
    tr.pred = {s: dict(v) for s, v in tr.pred.items()}
    with open(mine + ".tmp", "wb") as fh:
        pickle.dump(tr, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(mine + ".tmp", mine)
    tr.succ = collections.defaultdict(dict, tr.succ)
    tr.pred = collections.defaultdict(dict, tr.pred)
    return tr


# ---------------- context of one row ----------------

class Ctx:
    def __init__(self, r, d, world):
        self.r, self.w = r, world
        self.kind = r["kind"]
        self.q = json.loads(r["query"])
        self.rec = d["rec"] if d else None
        self.res = (d or {}).get("res") or {}
        self.got = split(r["got"])
        self.exp = split(r["expected"])
        self.broad = split(r.get("broad", ""))
        self.gotb = {base(x) for x in self.got}
        self.expb = {base(x) for x in self.exp}
        self.extra = self.gotb - self.expb
        self.missing = self.expb - self.gotb
        self.id = self.q["id"]
        self.fr = self.q.get("from_release")
        self.to = self.q.get("to_release")
        self.db = self.q.get("final_database")
        self.paths = self.res.get("paths") or []
        self.read_at = self.res.get("from_release") or self.fr or self.to
        tid = [x for x in (self.res.get("target_id") or []) if not str(x).startswith("...")]
        self.same_answer = {base(x) for x in tid} == self.gotb if self.res else None
        self.domain = domain(self.kind, self.q)

    def paths_to(self, gene_base):
        return [p for p in self.paths if base(p.get("ensembl_gene")) == gene_base]

    def entry_dbs(self, gene_base):
        out = set()
        for p in self.paths_to(gene_base):
            for e, dbs in zip(p["edges"], p["edge_dbs"]):
                if dbs:
                    out |= {x.split("@")[0] for x in dbs}
                    break
        return out

    def jumps(self, gene_base):
        js = [p.get("external_jump") for p in self.paths_to(gene_base)]
        return min(js) if js else None


def domain(kind, q):
    if kind in HIST_KINDS:
        return "history"
    if kind in NAME_KINDS:
        return "name"
    if kind in DB_KINDS:
        return "db"
    if q.get("final_database"):
        return "db"
    if C.ENS_ID.match(str(q.get("id"))):
        return "history"
    return "name"


# ---------------- the tests (as in the residual reading; kinds read through `domain`) ----------------

def t_entrezgene_label_rows(c):
    if c.domain != "db" or c.db != "EntrezGene":
        return None
    extra = [x for x in c.got if x not in c.exp]
    if not extra or any(x.isdigit() for x in extra) or not set(c.exp) <= set(c.got):
        return None
    lab = c.w.bm_pairs(c.to, "entrez", 2)
    labels = lab[0].get(base(c.id), set()) if lab else set()
    return f"extra {','.join(extra)}; BioMart {c.to} NCBI label of the gene: {','.join(sorted(labels)) or '-'}; " \
           f"label equal: {set(extra) <= labels}"


def t_biomart75_swissprot_entry_name(c):
    if c.domain != "db" or c.db != "Uniprot/SWISSPROT" or c.to != 75:
        return None
    if c.exp and all(re.match(r"^[A-Z0-9]+_[A-Z0-9]+$", x) for x in c.exp) and c.got and \
            all(re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", x) for x in c.got):
        return f"expected entry name {','.join(c.exp)}, got accession {','.join(c.got)}"
    return None


def t_one_spelling_per_db_release(c):
    if c.domain != "db":
        return None
    miss = set(c.exp) - set(c.got)
    if miss and all(any(m.lower() == g.lower() for g in c.got) for m in miss) and set(c.got) <= set(c.exp):
        return f"missing {','.join(sorted(miss))}: another case of {','.join(c.got)}"
    return None


def _orig_name(c):
    note = (c.rec or {}).get("note") or ""
    m = re.match(r"^(lower|upper|separators) of (.*)$", note)
    return m.group(2) if m else c.id


def t_mirbase_not_in_graph(c):
    if c.domain != "name":
        return None
    name = _orig_name(c)
    if not MIR.match(name):
        return None
    rel = c.w.bm_rel(c.fr) if c.fr else 115
    srcs = {s for g, (n, s) in (c.w.display(rel) if rel else {}).items() if n == name}
    if "miRBase" not in srcs:
        return None
    other = c.extra and "found through " + ",".join(sorted(set().union(*(c.entry_dbs(g) for g in c.extra))))
    return f"{name}: display source at {rel} {','.join(sorted(srcs))}; answer: {'nothing' if not c.got else other}"


def t_strict_walk_drops_living_id(c):
    if c.domain != "history":
        return None
    tr = c.w.tr
    qb = base(c.id)
    if qb not in c.gotb or qb in c.expb or c.to not in tr.present.get(qb, ()):
        return None
    start = c.fr if c.fr is not None else max(r for r in tr.present.get(qb, (0,)) if r <= 115)
    a, b = sorted((start, c.to))
    gap = [r for r in tr.rels if a <= r <= b and r not in tr.present[qb]]
    evs = []
    for s in tr.sessions_between(a, b):
        so = tr.succ[s].get(qb, set()) - {None}
        sp = tr.pred[s].get(qb, set()) - {None}
        if (so and qb not in so) or (sp and qb not in sp):
            evs.append(f"{tr.sessions[s]['old_release']}->{tr.sessions[s]['new_release']}: "
                       f"out {','.join(sorted(so))[:60] or '-'} in {','.join(sorted(sp))[:60] or '-'}")
    if gap:
        return None
    return f"{qb} in the GTF of every release {a}..{b}; events: {' | '.join(evs)[:300]}"


def t_vgnc_symbol_not_in_truth(c):
    if c.domain != "name" or not c.extra:
        return None
    if all(c.entry_dbs(g) and c.entry_dbs(g) <= {"VGNC Symbol"} for g in c.extra):
        return f"{','.join(sorted(c.extra))} entered through VGNC Symbol only"
    return None


def t_rfam_family_label(c):
    if c.domain != "name" or not c.extra:
        return None
    name = _orig_name(c)
    rel = c.w.bm_rel(c.read_at) or 115
    rfam_label = any(n == name and s == "RFAM" for n, s in c.w.display(rel).values())
    seen = [g for g in c.extra if c.paths_to(g)]
    if seen and all(c.entry_dbs(g) <= {"RFAM"} and c.entry_dbs(g) for g in seen) and (rfam_label or len(seen) == len(c.extra)):
        disp = c.w.display(115)
        names = collections.Counter(disp.get(g, ("?", "?"))[0] for g in c.extra)
        return f"{len(c.extra)} extra genes entered through RFAM only; their display names {dict(names.most_common(4))}"
    return None


def t_copies_apart_before_82(c):
    if c.domain != "history" or not c.missing or c.extra or (c.to or 999) >= 82:
        return None
    apart = {base(x) for x in (c.res.get("alternative_locus_genes") or [])}
    if c.missing <= apart:
        return f"missing {','.join(sorted(c.missing))} reported apart; target release {c.to}"
    return None


def t_gtf_lacks_core_gene(c):
    if c.domain != "history" or not c.extra:
        return None
    tr = c.w.tr
    outside = c.extra - {base(x) for x in c.broad}  # the genes the line's wider set does not allow
    absent = {g for g in outside if c.to not in tr.present.get(g, ())}
    if not outside or absent != outside:
        return None
    ids = c.w.core_ids(c.to)
    if ids and absent <= ids:
        return f"{','.join(sorted(absent))} in the core gene table of {c.to}, not in its GTF " \
               f"(GTF genes {len(tr.genes.get(c.to, {}))}, core genes {len(ids)})"
    return None


def t_path_through_resurrection(c):
    if not c.extra:
        return None
    hits = []
    for g in c.extra:
        for p in c.paths_to(g):
            if any(str(x).endswith(".Retired") for e in p["edges"] for x in e[:2]):
                hits.append(g)
                break
    if hits and set(hits) == c.extra:
        return f"{','.join(sorted(hits))} reached through " + ",".join(sorted(
            {x for g in hits for p in c.paths_to(g) for e in p["edges"] for x in map(str, e[:2]) if x.endswith(".Retired")}))
    return None


def t_event_version_contradicts_tables(c):
    if c.domain != "history" or not c.missing:
        return None
    tr = c.w.tr
    qb = base(c.id)
    out = []
    for ov, n, nv, s in c.w.events().get(qb, []):
        if n in c.missing:
            old = tr.sessions[s]["old_release"]
            have = tr.ver(qb, H.rel_int(old))
            if have and str(ov) != str(have):
                out.append(f"{qb}.{ov}->{n}.{nv} in session {s} ({old}->{tr.sessions[s]['new_release']}); "
                           f"{qb} has version {have} at {old}")
    if out and len(out) >= len(c.missing):
        return "; ".join(out)
    return None


def t_mouse_hgnc_projection(c):
    if c.w.org != "mus_musculus" or c.domain != "name" or c.got:
        return None
    name = _orig_name(c)
    srcs = {s for n, s in c.w.display(115).values() if n == name}
    if srcs == {"HGNC Symbol"}:
        return f"{name}: display source HGNC Symbol only; no MGI symbol"
    return None


def t_name_in_db_not_in_truth(c):
    if c.domain != "name" or not c.extra:
        return None
    foreign = {"VGNC Symbol", "Vega gene", "Havana gene", "Vega_gene"}
    if all(c.entry_dbs(g) and c.entry_dbs(g) <= foreign and not (c.jumps(g) or 0) for g in c.extra):
        return "; ".join(f"{g} entered through {'+'.join(sorted(c.entry_dbs(g)))}" for g in sorted(c.extra))
    return None


def t_old_name_split_judged_as_set(c):
    if c.kind not in ("name_old_to_current", "travel_name") or c.q["strategy"] != "best" or c.extra or not c.missing or not c.gotb:
        return None
    rel = c.w.bm_rel(c.fr)
    if rel is None:
        return None
    src = {g for g, (n, s) in c.w.display(rel).items() if n == c.id and rel in c.w.tr.present.get(g, ())}
    if len(src) == 1:
        (g,) = src
        return f"{c.id} at {rel}: one gene {g}, split by Ensembl into {len(c.expb)}; best returned {len(c.gotb)}"
    return None


def t_support_ranks_shared_name(c):
    if c.domain != "name" or c.q["strategy"] != "best" or c.extra or not c.missing or not c.gotb:
        return None
    tr = c.w.tr
    rel = c.w.bm_rel(c.read_at) or 115
    if c.fr is not None and c.fr != c.to:
        src = {g for g, (n, s) in c.w.display(rel).items() if n == c.id and rel in tr.present.get(g, ())}
        lead = {g: tr.walk_fwd(g, rel, c.to, True) for g in src}
        sup = {g: c.w.name_sources(rel, g, c.id) for g in src}
        to_got = [g for g in src if lead[g] & c.gotb]
        to_miss = [g for g in src if lead[g] & c.missing and not lead[g] & c.gotb]
        if not to_got or not to_miss:
            return None
        got_n = min(len(sup[g]) for g in to_got)
        miss_n = max(len(sup[g]) for g in to_miss)
    else:
        sup = {g: c.w.name_sources(rel, g, c.id) for g in c.expb}
        got_n = min(len(sup[g]) for g in c.gotb if g in sup) if c.gotb & set(sup) else -1
        miss_n = max(len(sup[g]) for g in c.missing)
    verdict = "more" if got_n > miss_n else ("equal" if got_n == miss_n else "fewer")
    s = "; ".join(f"{g}{'*' if g in c.gotb else ''}:{'+'.join(v) or '-'}" for g, v in sorted(sup.items()))[:400]
    return f"returned genes have {verdict} BioMart {rel} name sources than the missing ones; {s}"


def t_synonym_spelling_merged(c):
    if c.domain != "name" or not c.extra:
        return None
    syn = c.w.bm_pairs(115, "synonym", 1)
    spell = {}
    for g in c.extra:
        vs = {v for v in (syn[0].get(g, set()) if syn else set()) if v.lower() == c.id.lower()}
        for grp, col in NAME_SRC:
            pp = c.w.bm_pairs(115, grp, col)
            if pp:
                vs |= {v for v in pp[0].get(g, set()) if v.lower() == c.id.lower()}
        spell[g] = vs
    if all(spell[g] and c.id not in spell[g] for g in c.extra):
        return f"read as {c.res.get('graph_id')} ({c.res.get('match')}); the extra genes carry " + \
               "; ".join(f"{g}:{','.join(sorted(v))}" for g, v in sorted(spell.items()))
    return None


def t_name_in_other_source(c):
    if c.domain != "name" or not c.extra:
        return None
    if any((c.jumps(g) or 0) > 0 for g in c.extra):
        return None
    rel = c.w.bm_rel(c.read_at)
    if rel is None:
        return None
    name = _orig_name(c)
    out = []
    for g in sorted(c.extra):
        entries = {base(p["entry_gene"]) for p in c.paths_to(g)}
        srcs = set()
        for e in entries:
            srcs |= {s for s in c.w.name_sources(rel, e, name) if not s.startswith("display")}
        if not srcs:
            return None
        out.append(f"{g}<-{','.join(sorted(entries))}:{'+'.join(sorted(srcs))}")
    return f"BioMart {rel}: " + "; ".join(out)[:400]


def t_variant_lacks_broad(c):
    if c.kind not in ("case_variant", "separator_variant") or not c.extra:
        return None
    return t_name_in_other_source(c)


def t_external_bridge(c):
    if not c.extra:
        return None
    js = {g: c.jumps(g) for g in c.extra}
    if not all(j and j > 0 for j in js.values()):
        return None
    loc = (c.rec or {}).get("locus")
    parts = []
    for g in sorted(c.extra):
        p = c.paths_to(g)[0]
        steps = list(zip(p["edges"], p["edge_dbs"]))
        bridge = None
        for (e1, d1), (e2, d2) in zip(steps, steps[1:]):
            if d1 and d2 and e1[1] == e2[0] and not C.ENS_ID.match(str(e1[1])) and len(e1) > 3 and len(e2) > 3 \
                    and e1[3] is not None and e2[3] is not None and (e1[3] != c.read_at or e2[3] != c.read_at
                                                                      or base(e1[0]) != base(e2[1])):
                bridge = (e1, d1, e2, d2)
                break
        if bridge is None:
            parts.append(f"{g} jump {js[g]} (bridge not located)")
            continue
        e1, d1, e2, d2 = bridge
        x, y = int(e1[3]), int(e2[3])
        dbs = sorted({d.split("@")[0] for d in d1 + d2})
        kind = "same release" if x == y else ("next release" if abs(y - x) == 1 else f"gap {abs(y - x)}")
        lg = c.w.tr.locus(g, c.to)
        dist = "n/a"
        if loc and lg:
            if loc["assembly"] != lg["assembly"]:
                dist = "other assembly"
            elif loc["seqname"] != lg["seqname"]:
                dist = "other sequence"
            else:
                d = max(lg["start"] - loc["end"], loc["start"] - lg["end"], 0)
                dist = f"{d}" + ("" if loc["strand"] == lg["strand"] else " other strand")
        parts.append(f"{g} via {e1[1]} [{'+'.join(dbs)}] {kind} ({x}->{y}); distance {dist}")
    return "; ".join(parts)[:500]


def t_newer_version_kind_reversed(c):
    if c.kind != "nonexistent_version":
        return None
    return f"{c.id}: answered {','.join(c.got) or 'nothing'} (rule `version_newer_than_graph`: a version newer than the graph)"


def t_rule_decided_after_tree(c):
    rules = [x for x in str(c.r.get("rule") or "").replace("revised:", "").split("+") if x in MQ.RULES_JUDGE]
    before = c.r.get("_before")
    if not rules or before is None:
        return None
    allowed = {base(x) for x in (before.get("expect") or []) + (before.get("broad") or [])}
    if c.gotb <= allowed and (c.gotb or not before.get("expect")):
        return f"rule {'+'.join(rules)}: the answer is the line's expectation before the rule"
    return None


TESTS = [t_rule_decided_after_tree, t_newer_version_kind_reversed, t_entrezgene_label_rows, t_biomart75_swissprot_entry_name, t_one_spelling_per_db_release,
         t_mirbase_not_in_graph, t_mouse_hgnc_projection, t_strict_walk_drops_living_id, t_copies_apart_before_82,
         t_gtf_lacks_core_gene, t_path_through_resurrection, t_event_version_contradicts_tables,
         t_vgnc_symbol_not_in_truth, t_name_in_db_not_in_truth, t_rfam_family_label,
         t_old_name_split_judged_as_set, t_support_ranks_shared_name, t_synonym_spelling_merged, t_variant_lacks_broad,
         t_name_in_other_source, t_external_bridge]


def register(cause, cls, decision, text, test, before=None):
    """Add a family's cause: its class, the decision or finding behind it, one line, and its test (placed before the
    test of `before`, else before the generic `t_name_in_other_source` and `t_external_bridge`)."""
    CAUSE_CLASS[cause] = (cls, decision, text)
    test.__name__ = "t_" + cause
    names = [t.__name__[2:] for t in TESTS]
    i = names.index(before) if before in names else names.index("name_in_other_source")
    TESTS.insert(i, test)


def release_order(rows):
    """The rows' indices in the order of the releases their queries read (the start release, then the target): the
    truth's readers cache what they read per release (a release's tables are read once while its rows are classified),
    so rows met in the order of the scored file (by the hash of the qid) would read a release's tables again for
    nearly every row."""
    def key(i):
        q = json.loads(rows[i]["query"]) if rows[i].get("query") else {}
        fr, to = q.get("from_release"), q.get("to_release")
        return (fr is not None, fr or 0, to or 0)
    return sorted(range(len(rows)), key=key)


def classify_rows(org, rows, explain, und=False):
    """(cause, detail) for each row (a wrong.tsv.gz or undetermined.tsv.gz row), in the rows' order; explain: (qid,
    view) -> record. The rows are classified in the order of the releases they read (`release_order`); every row is
    classified on its own, so the order changes no cause."""
    w = World(org)
    out = [None] * len(rows)
    # lines that follow a rule decided after the tree (the harness's judge `rules`): their line before the rule
    ruled = [r for r in rows if any(x in MQ.RULES_JUDGE for x in str(r.get("rule") or "").replace("revised:", "").split("+"))]
    if ruled:
        from idtrack_reproducibility.families import lines as L
        sets = [n for n in L.set_names("harness", org) if L.manifest("harness", org, n).get("judge") != "rules"]
        before = L.load_qids("harness", org, {r["qid"] for r in ruled}, sets=sets)
        for r in ruled:
            r["_before"] = before.get(r["qid"])
    for i in release_order(rows):
        r = rows[i]
        view = r["variant"].split(",")[0]
        c = Ctx(r, explain.get((r["qid"], view)), w)
        cause, detail = "unexplained", ""
        for t in TESTS:
            try:
                dt = t(c)
            except Exception as e:  # noqa: BLE001  a test that cannot be applied explains nothing
                dt = None
                print("test error", t.__name__, r["qid"], type(e).__name__, e, file=sys.stderr)
            if dt:
                cause, detail = t.__name__[2:], dt
                break
        if cause == "unexplained" and und and c.extra and all(c.jumps(g) == 0 for g in c.extra):
            cause, detail = "history_branch", f"{','.join(sorted(c.extra))[:100]} by history only; {r['reason'][:80]}"
        if cause == "unexplained":
            detail = (f"{'no explain record; ' if c.rec is None else ''}match {c.res.get('match')} read at {c.read_at}; "
                      f"extra {','.join(sorted(c.extra))[:120]}; missing {','.join(sorted(c.missing))[:120]}; "
                      f"reason {r['reason'][:80]}")
        out[i] = (view, cause, detail, c.same_answer)
    return out


def write_summary(path, summ, classes):
    """`causes_summary.tsv`: stratum, cause, class, rule (the rule or finding behind the cause, docs/rules.md), answers,
    units; `summ` is {(stratum, cause): [answers, {units}]}."""
    with open(path + ".tmp", "w") as fo:
        fo.write("stratum\tcause\tclass\trule\tanswers\tunits\n")
        for (st, cause), (n, us) in sorted(summ.items()):
            cls, rule = (classes.get(cause) or ("?", "-"))[:2]
            fo.write(f"{st}\t{cause}\t{cls}\t{rule}\t{n}\t{len(us)}\n")
    os.replace(path + ".tmp", path)


def summary_from(causes_path, out_path, classes=None):
    """A causes table's summary written again (e.g. after the rules' names of a cause changed): one pass over it."""
    summ = collections.defaultdict(lambda: [0, set()])
    with gzip.open(causes_path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            s = summ[(r["stratum"], r["cause"])]
            s[0] += 1
            s[1].add(r["unit"])
    write_summary(out_path, summ, CAUSE_CLASS if classes is None else classes)


def summary_only(scored, suffix=""):
    """`causes<suffix>.tsv.gz` with its class column, and `causes<suffix>_summary.tsv`, written again from the stored
    rows under `CAUSE_CLASS` as it is now (the causes themselves are not classified again)."""
    path = os.path.join(scored, f"causes{suffix}.tsv.gz")
    summ = collections.defaultdict(lambda: [0, set()])
    with gzip.open(path, "rt") as fh, gzip.open(path + ".tmp", "wt") as fo:
        rd = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        fo.write("\t".join(rd.fieldnames) + "\n")
        for r in rd:
            r["class"] = CAUSE_CLASS.get(r["cause"], (r["class"],))[0]
            fo.write("\t".join(r[k] for k in rd.fieldnames) + "\n")
            x = summ[(r["stratum"], r["cause"])]
            x[0] += 1
            x[1].add(r["unit"])
    os.replace(path + ".tmp", path)
    write_summary(os.path.join(scored, f"causes{suffix}_summary.tsv"), summ, CAUSE_CLASS)
    print(open(os.path.join(scored, f"causes{suffix}_summary.tsv")).read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--scored", required=True)
    ap.add_argument("--explain", default=None)
    ap.add_argument("--und", action="store_true")
    ap.add_argument("--plugin", action="append", default=[], help="a family's module that registers its kinds and causes")
    ap.add_argument("--summary-only", action="store_true", help="classify nothing: write the stored causes' class column and "
                    "causes_summary.tsv again under the current classes and rules of the causes (renamed rules, a cause "
                    "moved to another class)")
    a = ap.parse_args()
    import importlib.util
    for p in a.plugin:
        spec = importlib.util.spec_from_file_location(os.path.basename(p)[:-3], p)
        spec.loader.exec_module(importlib.util.module_from_spec(spec))
    if a.summary_only:
        return summary_only(a.scored, "_und" if a.und else "")
    src = "undetermined" if a.und else "wrong"
    with gzip.open(os.path.join(a.scored, f"{src}.tsv.gz"), "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    ep = a.explain or os.path.join(a.scored, f"explain_{src}.jsonl.gz")
    explain = {}
    if os.path.exists(ep):
        with gzip.open(ep, "rt") as fh:
            for x in fh:
                d = json.loads(x)
                explain[(d["qid"], d["view"])] = d
    res = classify_rows(a.organism, rows, explain, und=a.und) if rows else []  # nothing to classify: empty tables
    suffix = "_und" if a.und else ""
    summ = collections.defaultdict(lambda: [0, set()])
    out = os.path.join(a.scored, f"causes{suffix}.tsv.gz")  # written beside, then renamed: never in place
    with gzip.open(out + ".tmp", "wt") as fo:
        fo.write("organism\tqid\tview\tkind\tstratum\tunit\tcause\tclass\tdetail\n")
        for r, (view, cause, detail, same) in zip(rows, res):
            fo.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in
                               (a.organism, r["qid"], view, r["kind"], r["stratum"], r["unit"], cause, CAUSE_CLASS[cause][0], detail)) + "\n")
            s = summ[(r["stratum"], cause)]
            s[0] += 1
            s[1].add(r["unit"])
    os.replace(out + ".tmp", out)
    write_summary(os.path.join(a.scored, f"causes{suffix}_summary.tsv"), summ, CAUSE_CLASS)
    diff = sum(1 for x in res if x[3] is False)
    print(f"{a.organism}: {len(rows)} {src} answers; {sum(1 for x in res if x[1] == 'unexplained')} unexplained; "
          f"{diff} explained calls gave another answer than the scored one")
    print(open(os.path.join(a.scored, f"causes{suffix}_summary.tsv")).read())


if __name__ == "__main__":
    main()
