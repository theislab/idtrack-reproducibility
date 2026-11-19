"""Check family C's reading of Ensembl's xrefs (the core dumps) against Ensembl's two other publications of them: the
FTP `tsv/` tables (87 on) and BioMart (the archive releases fetched by the harness). Nothing of idtrack.

Usage: check_truth_c.py <organism> [<organism> ...]   writes <families>/c/truth_check/<organism>.tsv

One row per (source, release, database): values in both with the same genes, in both with other genes, only in the
other source, only in the core reading (with up to three examples each).
"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families.c import truth_c as TC
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.families import common as C

BIOMART = {  # BioMart group -> (column, Ensembl db names)
    "entrez": (1, {"EntrezGene"}),
    "refseq_mrna": (1, {"RefSeq_mRNA"}), "refseq_ncrna": (1, {"RefSeq_ncRNA"}), "refseq_peptide": (1, {"RefSeq_peptide"}),
    "refseq_mrna_predicted": (1, {"RefSeq_mRNA_predicted"}), "refseq_peptide_predicted": (1, {"RefSeq_peptide_predicted"}),
    "swissprot": (1, {"Uniprot/SWISSPROT"}), "sptrembl": (1, {"Uniprot/SPTREMBL"}), "rfam": (1, {"RFAM"}),
}
TSV_DBS = {"entrez": {"EntrezGene"}, "refseq": None, "uniprot": {"Uniprot/SWISSPROT", "Uniprot/SPTREMBL", "Uniprot_isoform"}}


def core_pairs(org, rel, dbs):
    out = collections.defaultdict(set)
    for ot, oid, gid, db, disp, acc, lab, ver, it in TC.read_xref(org, rel):
        if db in dbs and acc:
            out[(db, acc)].add(gid)
    return out


def compare(a, b):
    """a: the other source, b: the core reading; per db: (same, differ, only_a, only_b, examples)."""
    out = collections.defaultdict(lambda: [0, 0, 0, 0, []])
    for k in set(a) | set(b):
        db = k[0]
        r = out[db]
        if k in a and k in b:
            if a[k] == b[k]:
                r[0] += 1
            else:
                r[1] += 1
                if len(r[4]) < 3:
                    r[4].append(f"{k[1]}: other {','.join(sorted(a[k]))[:60]} core {','.join(sorted(b[k]))[:60]}")
        elif k in a:
            r[2] += 1
            if len(r[4]) < 3:
                r[4].append(f"{k[1]} only other")
        else:
            r[3] += 1
    return out


def main():
    os.makedirs(os.path.join(C.BENCH, "c", "truth_check"), exist_ok=True)
    for org in sys.argv[1:]:
        rows = []
        xr = set(TC.xref_releases(org))
        for rel in sorted(xr):
            for kind, dbs in TSV_DBS.items():
                t = TC.tsv_pairs(org, rel, kind)
                if t is None:
                    continue
                dbs = dbs or {k[0] for k in t}
                t = {k: {TC.base(g) for g in v} for k, v in t.items() if k[0] in dbs}
                for db, r in compare(t, core_pairs(org, rel, dbs)).items():
                    rows.append(("tsv", rel, db, *r[:4], "; ".join(r[4])))
            for grp, (col, dbs) in BIOMART.items():
                pp = H.biomart_pairs(org, rel, grp, col) if rel in H.biomart_releases(org) else None
                if not pp:
                    continue
                db = next(iter(dbs))
                bm = {(db, v): set(gs) for v, gs in pp[1].items()}
                for d, r in compare(bm, core_pairs(org, rel, dbs)).items():
                    rows.append((f"biomart:{grp}", rel, d, *r[:4], "; ".join(r[4])))
            print(org, rel, flush=True)
        p = os.path.join(C.BENCH, "c", "truth_check", f"{org}.tsv")
        with open(p, "w") as fo:
            fo.write("source\trelease\tdb\tsame_genes\tother_genes\tonly_in_source\tonly_in_core\texamples\n")
            for r in rows:
                fo.write("\t".join(str(x) for x in r) + "\n")
        tot = collections.Counter()
        for r in rows:
            tot[(r[0].split(":")[0], "same")] += r[3]
            tot[(r[0].split(":")[0], "other")] += r[4]
            tot[(r[0].split(":")[0], "only_source")] += r[5]
            tot[(r[0].split(":")[0], "only_core")] += r[6]
        print(org, dict(tot), flush=True)


if __name__ == "__main__":
    main()
