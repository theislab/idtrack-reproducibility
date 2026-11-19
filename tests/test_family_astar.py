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
