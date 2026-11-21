"""The judge following rules decided in idtrack's first round of fixes (`no_mito_nuclear_step`, `mirbase_genes`), on an
invented truth store: the mitochondrial cut of the history walk, the miRBase entries of two loci, and the revisions."""
import pytest

from conftest import GENE_COLS, ORG, write_gz
from idtrack_reproducibility.harness import make_queries as M

X, Y, N, P, Q, R1, R2, S = (f"ENSSSCG0000000010{i}" for i in range(8))
XREF_COLS = ["object_type", "object_id", "gene_id", "db_name", "db_display_name", "acc", "label", "version", "info_type"]
SYN_COLS = ["synonym", "gene_id", "db_name", "label", "object_type"]


@pytest.fixture
def store(data_root):
    """Releases 1 and 2 on assembly Sscrofa1, 3 on Sscrofa2 (no chain between them).

    X: a gene of the mitochondrial genome at 1 and 2 (COX1); at 3 Ensembl gives its id to a nuclear lncRNA of
    chromosome 4, with no event (the id lives on); Y: the mitochondrial COX1 at 3, at X's old place. N: a nuclear gene
    with an Mt_tRNA_pseudogene beside it on chromosome 1 (no mitochondrial sequence).
    miRBase entries (core xrefs): `ssc-mir-a` on P and Q at 1, two places of chromosome 1 (two loci); `ssc-mir-b` on R1
    at 1-2 and R2 at 3, both on chromosome 2 (one gene at a time across the assembly change: kept); `ssc-mir-c` on R1
    at 1-2 (chromosome 2) and on S at 3 (chromosome 5: whole sequences of different names, two loci). `ssc-mir-a` is
    also an NCBI synonym of P at 2.
    """
    t = data_root / "truth" / ORG
    genes = {
        1: [(X, 1, "COX1", "protein_coding", "MT", 5000, 6500, "+", "primary"),
            ("ENSSSCG00000000201", 1, "MT-TF", "Mt_tRNA", "MT", 1, 70, "+", "primary"),
            (N, 1, "GN", "protein_coding", "1", 1, 1000, "+", "primary"),
            ("ENSSSCG00000000202", 1, "", "Mt_tRNA_pseudogene", "1", 2000, 2070, "+", "primary"),
            (P, 1, "ssc-mir-a", "miRNA", "1", 10000, 10080, "+", "primary"),
            (Q, 1, "ssc-mir-a", "miRNA", "1", 50000, 50080, "-", "primary"),
            (R1, 1, "ssc-mir-b", "miRNA", "2", 700, 780, "+", "primary")],
        2: [(X, 1, "COX1", "protein_coding", "MT", 5000, 6500, "+", "primary"),
            ("ENSSSCG00000000201", 1, "MT-TF", "Mt_tRNA", "MT", 1, 70, "+", "primary"),
            (N, 1, "GN", "protein_coding", "1", 1, 1000, "+", "primary"),
            (P, 1, "ssc-mir-a", "miRNA", "1", 10000, 10080, "+", "primary"),
            (R1, 1, "ssc-mir-b", "miRNA", "2", 700, 780, "+", "primary")],
        3: [(X, 2, "", "lncRNA", "4", 900, 1900, "-", "primary"),
            (Y, 1, "COX1", "protein_coding", "MT", 5000, 6500, "+", "primary"),
            ("ENSSSCG00000000201", 1, "MT-TF", "Mt_tRNA", "MT", 1, 70, "+", "primary"),
            (N, 1, "GN", "protein_coding", "1", 1, 1000, "+", "primary"),
            (R2, 1, "ssc-mir-b", "miRNA", "2", 650, 730, "+", "primary"),
            (S, 1, "ssc-mir-c", "miRNA", "5", 100, 180, "+", "primary")],
    }
    for r, rows in genes.items():
        write_gz(str(t / "gtf" / f"genes_{r}.tsv.gz"), GENE_COLS, rows)
    core = "sus_scrofa_core_3_1"
    write_gz(str(t / "core" / f"{core}.mapping_session.txt.gz"), None, [
        (1, "sus_scrofa_core_1_1", "sus_scrofa_core_2_1", "1", "2", "Sscrofa1", "Sscrofa1", "2026-01-01 00:00:00"),
        (2, "sus_scrofa_core_2_1", "sus_scrofa_core_3_1", "2", "3", "Sscrofa1", "Sscrofa2", "2026-02-01 00:00:00")])
    write_gz(str(t / "core" / f"{core}.stable_id_event.txt.gz"), None, [(N, 1, N, 1, 2, "gene", 1.0)])
    xr = {1: [("Gene", f"{P}.1", f"{P}.1", "miRBase", "miRBase", "MI1", "ssc-mir-a", 0, "MISC"),
              ("Gene", f"{Q}.1", f"{Q}.1", "miRBase", "miRBase", "MI1", "ssc-mir-a", 0, "MISC"),
              ("Gene", f"{R1}.1", f"{R1}.1", "miRBase", "miRBase", "MI2", "ssc-mir-b", 0, "MISC"),
              ("Gene", f"{R1}.1", f"{R1}.1", "miRBase", "miRBase", "MI3", "ssc-mir-c", 0, "MISC")],
          2: [("Gene", f"{P}.1", f"{P}.1", "miRBase", "miRBase", "MI1", "ssc-mir-a", 0, "MISC"),
              ("Gene", f"{R1}.1", f"{R1}.1", "miRBase", "miRBase", "MI2", "ssc-mir-b", 0, "MISC"),
              ("Gene", f"{R1}.1", f"{R1}.1", "miRBase", "miRBase", "MI3", "ssc-mir-c", 0, "MISC")],
          3: [("Gene", f"{R2}.1", f"{R2}.1", "miRBase", "miRBase", "MI2", "ssc-mir-b", 0, "MISC"),
              ("Gene", f"{S}.1", f"{S}.1", "miRBase", "miRBase", "MI3", "ssc-mir-c", 0, "MISC")]}
    for r, rows in xr.items():
        write_gz(str(t / "core_xref" / f"{r}.tsv.gz"), XREF_COLS, rows)
    write_gz(str(t / "core_synonyms" / "2.tsv.gz"), SYN_COLS, [("ssc-mir-a", P, "EntrezGene", "MIRA", "Gene")])
    (t / "sequences").mkdir(parents=True, exist_ok=True)
    for r, names in ((2, ["1", "2", "4", "5", "MT"]), (3, ["1", "2", "4", "5", "MT"])):
        with open(t / "sequences" / f"{r}.tsv", "w") as fo:
            fo.write("seqname\tcoord_system\tcoord_system_rank\tlength\tkaryotype_rank\tcodon_table\n")
            for i, n in enumerate(names):
                fo.write(f"{n}\tchromosome\t1\t1000000\t{i + 1}\t{'2' if n == 'MT' else ''}\n")
    with open(t / "SOURCES.tsv", "w") as fo:
        fo.write("what\turl\tbytes\tmd5\tfetched\n")
        for r, a in ((1, "Sscrofa1"), (2, "Sscrofa1"), (3, "Sscrofa2")):
            fo.write(f"gtf/{r}\thttps://example.org/release-{r}/gtf/sus_scrofa/Sus_scrofa.{a}.{r}.gtf.gz\t1\tx\t"
                     f"2026-11-11T00:00:00\n")
    return t


def test_the_mitochondrial_genome_is_read_from_its_own_biotypes(store):
    t = M.Truth(ORG, 3)
    assert t.mito_seqs(1) == {"MT"}  # an Mt_tRNA_pseudogene of chromosome 1 makes no mitochondrial sequence
    assert t.mito_steps() == {X}


def test_a_reused_mitochondrial_id_is_two_genes(store):
    plain, cut = M.Truth(ORG, 3), M.Truth(ORG, 3, mito_cut=True)
    assert plain.walk_fwd(X, 2, 3, True, living=True) == {X}
    assert cut.walk_fwd(X, 2, 3, True, living=True) == set()
    assert cut.walk_back(X, 3, 1, True, living=True) == set()
    assert cut.walk_fwd(N, 1, 3, True, living=True) == {N}  # a nuclear gene walks as before


def test_the_mitochondrial_revision_judges_by_place(store):
    C = M.Ctx78(ORG, 3)
    lines = [dict(qid="x", kind="cross_assembly", id=X, from_release=2, to_release=3, final_database=None,
                  strategy="best", rule="set", cmp="versioned", expect=[f"{X}.2"], broad=None, note=None),
             dict(qid="n", kind="cross_assembly", id=N, from_release=1, to_release=3, final_database=None,
                  strategy="best", rule="set", cmp="versioned", expect=[f"{N}.1"], broad=None, note=None)]
    out = {r["qid"]: r for v in M.rev_no_mito_nuclear_step(ORG, lines, C).values() for r in v}
    assert set(out) == {"x"}
    assert out["x"]["expect"] == [] and out["x"]["locus"]["seqname"] == "MT"
    assert out["x"]["revised_by"] == ["no_mito_nuclear_step"]


def test_mirbase_entries_of_two_loci(store):
    ml = M.MirbaseLoci(ORG, M.Ctx78(ORG, 3))
    assert ml.left_out("ssc-mir-a")      # two genes at one release at two places
    assert ml.left_out("SSC-MIR-A")      # without letter case
    assert not ml.left_out("ssc-mir-b")  # one gene at a time, chromosome 2 on both assemblies
    assert ml.left_out("ssc-mir-c")      # chromosome 2, then chromosome 5 of the next assembly
    assert not ml.left_out("ssc-mir-z")  # no entry


def test_a_left_out_entry_is_read_through_the_other_databases(store):
    C = M.Ctx78(ORG, 3)
    lines = [dict(qid="a", kind="name_ambiguous_to_set", id="ssc-mir-a", from_release=2, to_release=2, final_database=None,
                  strategy="all", rule="set_exact", cmp="versioned", expect=[f"{P}.1", f"{Q}.1"], broad=None, note=None),
             dict(qid="c", kind="name_unique_to_gene", id="ssc-mir-c", from_release=3, to_release=3, final_database=None,
                  strategy="best", rule="set", cmp="versioned", expect=[f"{S}.1"], broad=None, note=None),
             dict(qid="b", kind="name_unique_to_gene", id="ssc-mir-b", from_release=3, to_release=3, final_database=None,
                  strategy="best", rule="set", cmp="versioned", expect=[f"{R2}.1"], broad=None, note=None),
             dict(qid="v", kind="separator_variant", id="ssc_mir_c", from_release=None, to_release=3, final_database=None,
                  strategy="best", rule="set", cmp="versioned", expect=[f"{S}.1"], broad=None, note="separators of ssc-mir-c")]
    out = {r["qid"]: r for v in M.rev_mirbase_genes(ORG, lines, C).values() for r in v}
    assert set(out) == {"a", "c", "v"}            # ssc-mir-b is kept: its line stands
    assert out["a"]["expect"] == [f"{P}.1"]       # the NCBI synonym of P at 2
    assert out["c"]["expect"] == [] and out["v"]["expect"] == []  # no other database carries ssc-mir-c: nothing


def test_a_batch_item_without_a_value():
    import math

    import pandas as pd

    from idtrack_reproducibility.families.e import revise_e as RE
    assert all(RE.is_missing_item(v) for v in (None, math.nan, float("nan"), pd.NA, pd.NaT))
    assert not any(RE.is_missing_item(v) for v in ("", "nan", 0, b"TP53", object()))


def test_a_name_asked_by_its_variants():
    rec = dict(id="ｓｓｃ－ｍｉｒ－１２８５", note="fullwidth of ssc-mir-1285; decided: NFKC")
    assert "ssc-mir-1285" in M._names_asked(rec)
    assert "mml-mir-922" in M._names_asked(dict(id="mml-mir-922;", note=None))
    assert "ssc-mir-9788-1" in M._names_asked(dict(id="ssc_mir_9788_1", note="sep_swap of ssc-mir-9788-1"))
