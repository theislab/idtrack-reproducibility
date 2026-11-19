"""The place test across assemblies: chain files and assembly exceptions, on invented coordinates."""
import gzip
import os

from idtrack_reproducibility.harness import assembly_map as AM

from conftest import ORG


def write_chain(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fo:
        fo.write(text)


# chr1 of asmA aligns to chrB of asmB shifted by 1000: blocks 0-500 and 600-1000 (a gap of 100 on both sides)
FORWARD = "chain 1000 chr1 5000 + 0 1000 chrB 9000 + 1000 2000 1\n500 100 100\n400\n\n"
# chr2 of asmA aligns to the reverse strand of chrC (size 3000): asmA 0-1000 = chrC reverse 0-1000 = forward 2000-3000
REVERSE = "chain 1000 chr2 5000 + 0 1000 chrC 3000 - 0 1000 2\n1000\n\n"


def test_forward_chain_maps_by_offset(tmp_path):
    p = str(tmp_path / "c.chain.gz")
    write_chain(p, FORWARD)
    ch = AM.Chain(p)
    place, share, why = ch.map("chr1", 101, 200, "+")
    assert place == dict(seqname="chrB", start=1101, end=1200, strand="+") and share == 1.0 and why == ""
    assert ch.map("chrX", 1, 10, "+")[0] is None


def test_reverse_chain_flips_strand_and_coordinates(tmp_path):
    p = str(tmp_path / "c.chain.gz")
    write_chain(p, REVERSE)
    place, share, _ = AM.Chain(p).map("chr2", 1, 100, "+")
    assert place == dict(seqname="chrC", start=2901, end=3000, strand="-")


def test_a_place_mostly_in_a_gap_is_unmapped(tmp_path):
    p = str(tmp_path / "c.chain.gz")
    write_chain(p, FORWARD)
    place, share, why = AM.Chain(p).map("chr1", 471, 610, "+")  # 30 of 140 bases aligned before the gap, 10 after
    assert place is None and share < AM.MIN_ALIGNED and "aligned" in why


def test_mapper_through_a_third_assembly(truth_store):
    write_chain(str(truth_store / "chain" / "asmA_to_asmB.chain.gz"), FORWARD)
    write_chain(str(truth_store / "chain" / "asmB_to_asmC.chain.gz"),
                "chain 1000 chrB 9000 + 1000 2000 chrZ 9000 + 0 1000 3\n1000\n\n")
    m = AM.Mapper(ORG)
    out = m.map(dict(assembly="asmA", seqname="chr1", start=101, end=200, strand="+"), "asmC")
    assert (out["seqname"], out["start"], out["end"], out["via"]) == ("chrZ", 101, 200, "asmA -> asmB -> asmC")
    assert "unmapped" in m.map(dict(assembly="asmC", seqname="chrZ", start=1, end=2, strand="+"), "asmA")


def test_overlap_needs_sequence_strand_and_coordinates():
    place = dict(seqname="1", start=100, end=200, strand="+")
    assert AM.overlaps(place, dict(seqname="1", start=200, end=300, strand="+"))
    assert not AM.overlaps(place, dict(seqname="1", start=201, end=300, strand="+"))
    assert not AM.overlaps(place, dict(seqname="1", start=150, end=160, strand="-"))
    assert not AM.overlaps(place, None)


def test_a_y_par_place_is_projected_onto_x(truth_store):
    os.makedirs(truth_store / "exceptions")
    (truth_store / "exceptions" / "3.tsv").write_text(
        "seqname\tstart\tend\ttype\tref_seqname\tref_start\tref_end\tori\tseq_length\tref_length\tfrom_release\n"
        "Y\t10001\t2781479\tPAR\tX\t10001\t2781479\t1\t57227415\t156040895\t3\n"
        "Y\t56887903\t57217415\tPAR\tX\t155701383\t156030895\t1\t57227415\t156040895\t3\n")
    ex = AM.Exceptions(ORG)
    p = ex.project(dict(seqname="Y", start=56887903 + 10, end=56887903 + 20, strand="+"), 3)
    assert (p["seqname"], p["start"], p["end"], p["projected_from"]) == ("X", 155701393, 155701403, "Y")
    outside = dict(seqname="Y", start=3000000, end=3000100, strand="+")
    assert ex.project(outside, 3) == outside
    assert ex.project(dict(seqname="1", start=1, end=2, strand="+"), 3)["seqname"] == "1"
