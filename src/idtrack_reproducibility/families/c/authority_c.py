"""Family C's authority measure (docs/methods.md principle 2): where each database's own file (`auth` of a line: the
genes the database itself attaches to its identifier now) and Ensembl's truth (`expect`) differ, the difference is
Ensembl's and is reported beside idtrack's rate, never as a wrong answer of idtrack; and what a user who trusts the
database gets from idtrack.

Usage: authority_c.py --organism <o> --run <run folder>/<organism> [--sets ...]   writes <run>/authority.tsv

Per stratum and view, over the lines with an `auth` (`["none"]`: the database attaches no gene or has retired the
identifier): `lines`; `truths_agree` (the authority's genes equal Ensembl's expected genes, by stable id);
`truths_differ` (and its share: the size of the limitation "Ensembl against the authority"); `auth_only_none` (the
authority has nothing, Ensembl a gene), `ens_only_none` (Ensembl nothing, the authority a gene); and for idtrack's
answers: `idtrack_as_auth` (best: nothing when the authority has nothing, else a non-empty part of its genes; all:
exactly its genes) and its share.
"""
import argparse
import collections
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = os.path.dirname(HERE)
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L


def main():
    import pyarrow.parquet as pq
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--sets", default=None)
    a = ap.parse_args()
    d = L.query_dir("c", a.organism)
    lines = {}
    for s in L.resolve_sets("c", a.organism, a.sets.split(",") if a.sets else None, "revised"):
        p = L._set_file(d, s)
        for b in pq.ParquetFile(p).iter_batches(batch_size=500_000, columns=["qid", "stratum", "expect", "extra"]) \
                if p.endswith(".parquet") else []:
            for r in b.to_pylist():
                ex = r["extra"]
                if not ex or '"auth"' not in ex:
                    continue
                auth = json.loads(ex).get("auth")
                if auth is None:
                    continue
                lines[r["qid"]] = (r["stratum"], {C.base_id(x) for x in (r["expect"] or [])},
                                   set() if auth == ["none"] else set(auth))
    stats = collections.defaultdict(collections.Counter)
    seen = set()
    for f in sorted(glob.glob(os.path.join(a.run, "answers", "*.parquet"))):
        t = pq.read_table(f, columns=["qid", "view", "outcome", "got"]).to_pylist()
        for r in t:
            ln = lines.get(r["qid"])
            if ln is None or (r["qid"], r["view"]) in seen:
                continue
            seen.add((r["qid"], r["view"]))
            st, E, A = ln
            k = stats[(st, r["view"])]
            k["lines"] += 1
            k["truths_agree" if E == A else "truths_differ"] += 1
            if E and not A:
                k["auth_only_none"] += 1
            if A and not E:
                k["ens_only_none"] += 1
            if r["outcome"] != "answer":
                k["idtrack_no_answer"] += 1
                continue
            g = {C.base_id(x) for x in (r["got"] or [])}
            ok = (g == A) if r["view"] == "all" else ((not g and not A) or (g and A and g <= A))
            k["idtrack_as_auth"] += bool(ok)
    out = os.path.join(a.run, "authority.tsv")
    cols = ["lines", "truths_agree", "truths_differ", "auth_only_none", "ens_only_none", "idtrack_as_auth",
            "idtrack_no_answer"]
    with open(out, "w") as fo:
        fo.write("organism\tstratum\tview\t" + "\t".join(cols) + "\tshare_truths_differ\tshare_idtrack_as_auth\n")
        for (st, v), k in sorted(stats.items()):
            n = k["lines"] or 1
            fo.write(f"{a.organism}\t{st}\t{v}\t" + "\t".join(str(k[c]) for c in cols) +
                     f"\t{k['truths_differ'] / n:.4f}\t{k['idtrack_as_auth'] / n:.4f}\n")
    print(open(out).read())


if __name__ == "__main__":
    main()
