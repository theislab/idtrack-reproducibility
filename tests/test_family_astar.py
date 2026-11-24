"""Family A* (the hardest lineages): its selection's building blocks on invented inputs."""
from idtrack_reproducibility.families.astar import gen as G
from idtrack_reproducibility.families.astar import selection as SEL


def test_lineage_components_are_named_by_their_smallest_id():
    name, size = SEL.components([("B", "C"), ("C", "D"), ("X", "Y"), ("A", "D")])
    assert name["D"] == "A" and name["B"] == "A" and name["Y"] == "X"
    assert size == {"A": 4, "X": 2}


def test_overlapping_genes_on_the_same_strand_only():
    genes = {"G1": dict(seqname="1", strand="+", start="100", end="200", region="primary"),
             "G2": dict(seqname="1", strand="+", start="150", end="300", region="primary"),
             "G3": dict(seqname="1", strand="-", start="150", end="300", region="primary"),
             "G4": dict(seqname="1", strand="+", start="301", end="400", region="primary"),
             "G5": dict(seqname="1", strand="+", start="120", end="130", region="core_only")}
    ov = SEL.Overlaps(genes)
    assert ov.others("1", "+", 100, 200, "G1") == ["G2"]
    assert ov.others("2", "+", 1, 10, "G1") == []


def test_spellings_of_one_value_are_grouped_by_letter_case():
    assert G.spellings_of(["IFN-DELTA-8", "IFN-delta-8", "TP53", "Tp53x"]) == [["IFN-DELTA-8", "IFN-delta-8"]]
    assert G.spellings_of(["a", "b"]) == []


def test_recorded_ids_are_read_for_a_published_selection_only():
    rec = SEL._recorded("homo_sapiens.r3")
    assert {r["pattern"] for r in rec} == {"hunt_history", "travel_weak_cell"}
    assert SEL._recorded("sus_scrofa") == []  # published, took no external id
    assert SEL._recorded("sus_scrofa.r99") is None  # a new selection reads the current sources
    assert SEL.weak_cell_ids("sus_scrofa", "sus_scrofa.r99") == set()


def test_a_versioned_line_the_mito_rule_leaves_without_successor_expects_nothing_on_the_nuclear_side():
    from idtrack_reproducibility.families.astar import revise_rules as RR

    class Tr:  # an id on a nuclear contig at 48-54 whose stable id Ensembl gives a mitochondrial gene at 55
        genes = {48: {"ENSG00000000001": dict(seqname="NT_1")}, 55: {"ENSG00000000001": dict(seqname="MT")}}

        def walk_series(self, g, R, ts, forward, living):
            return {t: set() for t in ts}  # the rule refuses the step: no successor

        def mito_seqs(self, r):
            return {"MT"}

    class Gen:
        tr = Tr()

    nuclear = dict(id="ENSG00000000001.1", kind="astar_from_v", from_release=48, to_release=100, expect=["ENSG00000000001.2"])
    assert RR.gene_line_without_successor(Gen, nuclear) == dict(expect=[], broad=None)
    default = dict(id="ENSG00000000001.1", kind="astar_default_v", to_release=100, read_at=48, expect=["ENSG00000000001.2"])
    assert RR.gene_line_without_successor(Gen, default) == dict(expect=[], broad=None)
    mito = dict(id="ENSG00000000001.2", kind="astar_from_v", from_release=55, to_release=100)
    assert RR.gene_line_without_successor(Gen, mito)["open"].startswith("rule `no_mito_nuclear_step`")
    assert RR.gene_line_without_successor(Gen, dict(nuclear, to_release=48)) is None  # not across the step
