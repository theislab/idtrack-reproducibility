"""A* revision sets for a rule decided after the sets were published, whose truth the generator now writes (judge
`rules`): rule `no_mito_nuclear_step` (truth.TruthA, `mito_cut`).

The lines of the ids whose lineage holds a step the rule refuses (`Truth.mito_steps`; every line of their units, gene
lines and the transcript and translation lines of their genes) are made again twice, by the generator as it is now
(`gen.Gen.lines_of`, `gen_tx.lines`) with and without the rule, and compared with the latest published line of each
qid: a field the rule changes takes its value under the rule; every other field stays as published (a line already as
the rule, published by an earlier run, is left: a rerun publishes only what is missing). A line that
changes goes into the revision set of the set that published it (`<set>r<k>`, the next free name; `revises` that set
and its revisions), under its qid. Nothing else is touched.

Usage: revise_rules.py <organism> [--dry-run [--lines-to <jsonl.gz>]]   (--lines-to: the revision lines a dry run would
       publish, for reading them before they are published)
"""
import argparse
import collections
import os
import re

from idtrack_reproducibility.families import lines as L
from idtrack_reproducibility.families.astar import gen as GN
from idtrack_reproducibility.families.astar import gen_tx as GT
from idtrack_reproducibility.harness import make_queries as M

RULE = "no_mito_nuclear_step"
ME = os.path.abspath(__file__)


def made(G, org, genes, unit_of):
    """{qid: line} of the gene lines and the transcript and translation lines of `genes`, by generator G."""
    out = {}
    for g in sorted(genes):
        for q in G.lines_of(g, unit_of[g]):
            out[q["qid"]] = q
    txs = GT.collect(org, G.tr, set(genes))
    for q in GT.lines(G, txs, unit_of, G.T):
        out[q["qid"]] = q
    return out


def db_without_successor(G, line):
    """The truth of an `astar_db` line whose walk the rule leaves without a successor at the target (the generator
    makes no db line of such a walk, so the published line holds the values of the walk without the rule), or None
    where the walk reaches a gene. On the nuclear side of the step the id names a nuclear gene that has no ancestor or
    successor across it: nothing is expected. On the mitochondrial side the decided answer is a bridge to the
    mitochondrial gene of the other side (Ensembl gave it a new id), whose values in the name database the judge
    cannot tell from Ensembl's history: the line is open."""
    g, R, t = M.gene_base(line["id"]), line["from_release"], line["to_release"]
    if G.tr.walk_series(g, R, [t], True, living=True)[t]:
        return None
    on_mito = (G.tr.genes.get(R, {}).get(g) or {}).get("seqname") in G.tr.mito_seqs(R)
    if on_mito:
        return dict(expect=[], broad=None, open=f"rule `{RULE}`: Ensembl's history gives the mitochondrial gene no "
                    f"successor at {t}; a bridge to the mitochondrial gene of the other side is not judged by the "
                    "database's values")
    return dict(expect=[], broad=None)


def gene_line_without_successor(G, line):
    """The truth of a gene line (not `astar_db`) whose walk the rule leaves without a successor at the target: the
    generator makes no line of a walk that reaches nothing (the versioned lines, `astar_*_v`), so the published line
    holds the successors of the walk without the rule. Read where the line reads (`from_release`, or `read_at` for a
    line asked without a release), as `db_without_successor` judges it: on the nuclear side of the step nothing is
    expected; on the mitochondrial side the line is open. None where the walk under the rule reaches a gene."""
    g, t = M.gene_base(line["id"]), line["to_release"]
    R = line.get("from_release") or line.get("read_at")
    if R is None or t <= R or G.tr.walk_series(g, R, [t], True, living=True)[t]:
        return None
    on_mito = (G.tr.genes.get(R, {}).get(g) or {}).get("seqname") in G.tr.mito_seqs(R)
    if on_mito:
        return dict(expect=[], broad=None, open=f"rule `{RULE}`: Ensembl's history gives the mitochondrial gene no "
                    f"successor at {t}")
    return dict(expect=[], broad=None)


def base_of(name):
    """The published set a revision set revises: `v1r3` -> `v1`, `tx1` -> `tx1`."""
    m = re.match(r"^(.*?\d+)(r\d+)?$", name)
    return m.group(1) if m else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--target", type=int, default=115)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--lines-to", default=None)
    a = ap.parse_args()
    org = a.organism
    G = GN.Gen(org, a.target)
    steps = G.tr.mito_steps()
    GN.log("ids of refused steps", sorted(steps))
    if not steps:
        GN.log("no step the rule refuses: no revision")
        return
    names = L.set_names("astar", org)
    # the units (lineages) of the refused ids, then every published line of those units (the latest line of each qid)
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    d = L.query_dir("astar", org)
    units, latest = set(), {}
    for s in names:
        p = L._set_file(d, s)
        if p.endswith(".parquet"):  # in batches: a set holds up to tens of millions of lines
            want = pa.array(sorted(steps))
            for b in pq.ParquetFile(p).iter_batches(batch_size=1_000_000, columns=["id", "unit"]):
                base = pc.replace_substring_regex(b.column("id"), pattern=r"\.\d+$", replacement="")
                units |= set(b.filter(pc.is_in(base, value_set=want)).column("unit").to_pylist())
        else:
            units |= {q["unit"] for q in L.iter_set(p) if M.gene_base(q["id"]) in steps}
    GN.log("units", sorted(units))
    for s in names:
        p = L._set_file(d, s)
        rows = []
        if p.endswith(".parquet"):
            for batch in pq.ParquetFile(p).iter_batches(batch_size=500_000):
                rows += [L.from_row(r) for r in
                         batch.filter(pc.is_in(batch.column("unit"), value_set=pa.array(sorted(units)))).to_pylist()]
        else:
            rows = [q for q in L.iter_set(p) if q.get("unit") in units]
        for q in rows:
            q["set"] = s
            latest[q["qid"]] = q
    unit_of = {}
    for q in latest.values():
        unit_of[M.gene_base(q.get("gene_at_origin") or q["id"])] = q["unit"]
    genes = set(unit_of)
    GN.log("published lines of the units", len(latest), "genes", len(genes))
    with_rule = made(G, org, genes, unit_of)
    without = made(GN.Gen(org, a.target, mito_cut=False), org, genes, unit_of)
    out = collections.defaultdict(list)
    c = collections.Counter()
    for qid, old in latest.items():
        x, y = with_rule.get(qid), without.get(qid)
        if x is None and old["kind"] == "astar_db":  # the generator makes no db line of a walk that reaches nothing
            tru = db_without_successor(G, old)
            if tru is not None:
                y = y or {k: old.get(k) for k in tru}
                x = dict(y, **tru)
                x.pop("spellings", None)
        if x is None and old["kind"] != "astar_db":  # no line of a walk that reaches nothing
            tru = gene_line_without_successor(G, old)
            if tru is not None:
                y = y or {k: old.get(k) for k in tru}
                x = dict(y, **tru)
        if x is None or y is None:
            c[("not made again", old["kind"])] += 1
            continue
        keys = {k for k in set(x) | set(y) if x.get(k) != y.get(k)}
        if not keys:
            c["unchanged by the rule"] += 1
            continue
        new = dict(old)
        for k in keys:
            if x.get(k) is None:
                new.pop(k, None)
            else:
                new[k] = x[k]
        new.pop("set", None)
        if all(new.get(k) == old.get(k) for k in keys):
            c["already as the rule"] += 1
            continue
        new["revised_by"] = list(old.get("revised_by") or []) + [RULE]
        out[base_of(old["set"])].append(new)
        c[("changed", old["stratum"])] += 1
    GN.log(org, dict(c))
    for base, rows in sorted(out.items()):
        chain = [s for s in names if base_of(s) == base]
        k = max([int(m.group(1)) for s in chain for m in [re.match(rf"^{re.escape(base)}r(\d+)$", s)] if m], default=0) + 1
        name = f"{base}r{k}"
        if a.dry_run:
            GN.log("dry run", name, len(rows), "revises", chain)
            if a.lines_to:
                import gzip
                from idtrack_reproducibility.families import common as C
                with gzip.open(a.lines_to, "at") as fo:
                    for q in rows:
                        fo.write(C.dumps(q) + "\n")
            continue
        man = L.publish("astar", org, name, rows, f"rule `{RULE}` (decided after {base}): the lines of the lineages "
                        f"it changes (revise_rules.py)", ME, revises=chain, judge="rules",
                        sources=[f"as {base}"])
        GN.log("published", name, man["n"], man["sha256"])


if __name__ == "__main__":
    main()
