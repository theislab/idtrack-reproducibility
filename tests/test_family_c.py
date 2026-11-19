"""Family C (external identifiers): its readers of Ensembl's cross-references and of the graph's scope, the classes of
accession versions, its samples, a correction of the judge and its causes, on invented values."""
import collections
import gzip
import importlib
import sys
import types

import pytest

from idtrack_reproducibility.families.c import gen_c as GC
from idtrack_reproducibility.families.c import revise_c as RC
from idtrack_reproducibility.families.c import truth_c as TC

XREF_COLS = ["object_type", "object_id", "gene", "db_name", "db_display_name", "acc", "label", "version", "info_type"]


def write_xref(truth, org, rel, rows):
    d = truth / org / "core_xref"
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / f"{rel}.tsv.gz", "wt") as fo:
        fo.write("\t".join(XREF_COLS) + "\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")


class Scope:
    def __init__(self, inc, form):
        self.inc, self.form = inc, form

    def included(self, name, rel):
        return rel in self.inc.get(name, ())


def test_versions_are_cut_from_ensembl_ids_only():
    assert TC.base("ENSG00000141510.20") == "ENSG00000141510" and TC.base("LRG_321.1") == "LRG_321"
    assert TC.base("NM_000546.6") == "NM_000546.6" and TC.base("ENSG00000141510") == "ENSG00000141510"


def test_mirbase_files_served_as_html_are_unwrapped(tmp_path):
    p = tmp_path / "hairpin.fa.gz"
    with gzip.open(p, "wt") as fo:
        fo.write("<p>&gt;hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop<br>UGGGAUGAGG<br></p>")
    assert TC.mirbase_text(str(p)).splitlines() == [">hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop", "UGGGAUGAGG"]
    q = tmp_path / "plain.gff3.gz"
    with gzip.open(q, "wt") as fo:
        fo.write("# plain\nchr1\t.\tmiRNA\n")
    assert TC.mirbase_text(str(q)) == "# plain\nchr1\t.\tmiRNA\n"


def test_the_scope_holds_a_database_only_at_releases_served_by_the_assembly_it_lists(tmp_path, monkeypatch):
    cfg = tmp_path / "tree" / "idtrack" / "default_config"
    cfg.mkdir(parents=True)
    (cfg / "homo_sapiens_externals_modified.yml").write_text("(read by the stand-in below)\n")
    config = {"homo_sapiens": {
        "transcript": {"RefSeq_mRNA": {"Assembly": {"38": {"Include": True, "Ensembl release": "100,101"},
                                                    "37": {"Include": True, "Ensembl release": "75"}}}},
        "gene": {"miRBase": {"Assembly": {"38": {"Include": False, "Ensembl release": "100,101"}}}}}}
    monkeypatch.setitem(sys.modules, "yaml", types.SimpleNamespace(safe_load=lambda fh: config))
    sc = TC.Scope("homo_sapiens", str(tmp_path / "tree"), {75: "GRCh37", 100: "GRCh38", 101: "GRCh37"})
    assert sc.included("RefSeq_mRNA", 100) and sc.included("RefSeq_mRNA", 75)
    assert not sc.included("RefSeq_mRNA", 101)  # listed under 38, but the release's GTF is of GRCh37
    assert not sc.included("miRBase", 100) and sc.form["RefSeq_mRNA"] == "transcript"


def test_cross_references_are_held_by_the_graph_only_under_its_name_and_object_form(tmp_path, monkeypatch):
    monkeypatch.setattr(TC, "TRUTH", str(tmp_path))
    write_xref(tmp_path, "sus_scrofa", 100, [
        ("Translation", 1, "ENSSSCG00000000001.3", "RefSeq_peptide", "RefSeq peptide", "NP_000001", "NP_000001", 2, "SEQUENCE_MATCH"),
        ("Transcript", 5, "ENSSSCG00000000002", "RefSeq_peptide", "RefSeq peptide", "NP_000002", "NP_000002", 1, "SEQUENCE_MATCH"),
        ("Gene", 1, "ENSSSCG00000000001", "EntrezGene", "NCBI gene (formerly Entrezgene)", "7157", "TP53", 0, "DEPENDENT"),
        ("Gene", 9, "ENSMUSG00000059552", "EntrezGene", "NCBI gene (formerly Entrezgene)", "22059", "Trp53", 0, "DEPENDENT"),
        ("Gene", 3, "ENSSSCG00000000003", "RFAM", "RFAM", "RF00005", "tRNA", 0, "DIRECT")])
    write_xref(tmp_path, "sus_scrofa", 101, [
        ("Translation", 1, "ENSSSCG00000000001", "RefSeq_peptide", "RefSeq peptide", "NP_000001", "NP_000001", 3, "SEQUENCE_MATCH")])
    scope = Scope({"RefSeq_peptide": {100, 101}, "EntrezGene": {100}}, {"RefSeq_peptide": "translation", "EntrezGene": "gene"})
    X = TC.Xrefs("sus_scrofa", scope, families=["refseq", "ncbi", "ncbi_label", "rfam"], log=lambda *a: None)
    assert X.at("refseq", "NP_000001", 100) == (frozenset({"ENSSSCG00000000001"}), frozenset({"2"}), True)
    assert X.at("refseq", "NP_000002", 100)[2] is False  # on a transcript: the configuration reads it from translations
    assert X.at("ncbi", "7157", 100) == (frozenset({"ENSSSCG00000000001"}), frozenset(), True)  # accession: db_name
    assert X.at("ncbi_label", "TP53", 100)[2] is False  # a label goes under db_display_name, which this scope leaves out
    assert "22059" not in X.H["ncbi"]  # another organism's gene
    assert X.labels["rfam"]["RF00005"] == {"tRNA"} and X.at("rfam", "RF00005", 100)[2] is False
    assert X.last_in_graph("refseq", "NP_000001", 115) == 101 and X.last_in_graph("refseq", "NP_000002", 115) is None
    assert X.last_with_version("refseq", "NP_000001", "2", 115) == 100
    assert X.versions_in_graph("refseq", "NP_000001") == {"2", "3"}


def test_accession_versions_are_classed_by_what_the_window_had():
    held = lambda gs, vs: (frozenset(gs), frozenset(vs), True)  # noqa: E731
    X = object.__new__(TC.Xrefs)
    X.H = {"refseq": {"NM_1": {100: held({"G1"}, {"2"}), 105: held({"G1"}, {"3"}), 115: held({"G1"}, {"5"})},
                      "NM_9": {100: held({"G9"}, {"1", "2"})}}}
    gc = object.__new__(GC.GenC)
    gc.T, gc.X, gc.org, gc.counts = 115, X, "sus_scrofa", collections.Counter()
    gc.tr = types.SimpleNamespace(cmp_for=lambda R: "versioned")
    gc.walked = lambda gs, R, t: dict(cmp="versioned", expect=sorted(gs), walked_from=R)
    gc.at_release = lambda gs, R: dict(cmp="versioned", expect=sorted(gs))
    lines = list(gc.version_strata("refseq", "c2", catalog={"NM_1": {4, 7}}))
    got = {(r["kind"], r["id"]): r for r in lines}
    assert set(got) == {("c2_oldver", "NM_1.2"), ("c2_oldver", "NM_1.3"), ("c2_oldver", "NM_9.1"), ("c2_oldver", "NM_9.2"),
                        ("c2_ver_lower", "NM_1.1"), ("c2_ver_gap", "NM_1.4"), ("c2_ver_higher", "NM_1.6"),
                        ("c2_ver_higher", "NM_1.10"), ("c2_ver_catalog", "NM_1.7"), ("c2_ver_higher_retired", "NM_9.3")}
    assert got[("c2_oldver", "NM_1.2")]["walked_from"] == 100 and got[("c2_oldver", "NM_1.3")]["walked_from"] == 105
    for k in (("c2_ver_lower", "NM_1.1"), ("c2_ver_gap", "NM_1.4"), ("c2_ver_higher_retired", "NM_9.3")):
        assert got[k]["expect"] == []  # a version never had: nothing
    for k in (("c2_ver_higher", "NM_1.6"), ("c2_ver_catalog", "NM_1.7")):
        assert got[k]["expect"] == ["G1"]  # above every version of an accession held at the target: its gene
    assert got[("c2_ver_gap", "NM_1.4")]["qid"] == "c2_ver_gap|NM_1.4|-|115" and got[("c2_ver_gap", "NM_1.4")]["cls"] == "NM"


def test_samples_are_seeded_sorted_and_whole_when_small():
    pop = [f"v{i:03d}" for i in range(300)]
    a = GC.sample("sus_scrofa", "c7_other", pop, 10)
    assert len(a) == 10 and a == sorted(a) and a == GC.sample("sus_scrofa", "c7_other", list(reversed(pop)), 10)
    assert a != GC.sample("macaca_mulatta", "c7_other", pop, 10)
    assert GC.sample("sus_scrofa", "x", ["b", "a"], 10) == ["a", "b"]


def test_rfam_labels_written_as_accessions_are_left_out_of_the_label_lines():
    lines = [dict(qid="1", kind="c4_rfam_label_at", id="RF00005"), dict(qid="2", kind="c4_rfam_label_at", id="tRNA"),
             dict(qid="3", kind="c4_rfam_acc_at", id="RF00005"), dict(qid="4", kind="c7_rfam_label", id="RF00001"),
             dict(qid="5", kind="c4_rfam_label_retired", id="RF00002", revised_by=["rfam_label_accession"])]
    out = list(RC.rfam_label_accession("sus_scrofa", lines, {}))
    assert [r["qid"] for r in out] == ["1", "4"]
    assert all(r["excluded"] and r["revised_by"] == ["rfam_label_accession"] for r in out)
    assert "excluded" not in lines[0]  # the published line is not changed


@pytest.fixture
def causes_c(monkeypatch):
    """C's plugin registered into a copy of the classifier's registry (the plugin registers itself on import)."""
    from idtrack_reproducibility.families import causes as K
    monkeypatch.setattr(K, "TESTS", list(K.TESTS))
    monkeypatch.setattr(K, "CAUSE_CLASS", dict(K.CAUSE_CLASS))
    for name in ("HIST_KINDS", "NAME_KINDS", "DB_KINDS"):
        monkeypatch.setattr(K, name, set(getattr(K, name)))
    sys.modules.pop("idtrack_reproducibility.families.c.causes_c", None)
    yield importlib.import_module("idtrack_reproducibility.families.c.causes_c"), K
    sys.modules.pop("idtrack_reproducibility.families.c.causes_c", None)


def test_causes_of_versions_never_had_and_of_ncbi_symbols_by_release(causes_c):
    CC, K = causes_c
    c = dict(kind="c2_ver_gap", got=["ENSSSCG00000000001.3"], exp=[], res={"match": "version"}, id="NM_1.4")
    assert "never had" in CC.accession_version_never_had(types.SimpleNamespace(**c))
    assert CC.accession_version_never_had(types.SimpleNamespace(**dict(c, exp=["G1"]))) is None
    assert CC.accession_version_never_had(types.SimpleNamespace(**dict(c, kind="c2_ver_higher"))) is None
    n = dict(db="NCBI gene (formerly Entrezgene)", missing={"TP53"}, exp=["TP53"], got=[], id="ENSSSCG00000000001.3")
    assert "`EntrezGene`" in CC.ncbi_label_name_by_release(types.SimpleNamespace(**n, to=87))
    assert "`NCBI gene`" in CC.ncbi_label_name_by_release(types.SimpleNamespace(**n, to=88))
    assert CC.ncbi_label_name_by_release(types.SimpleNamespace(**n, to=100)) is None
    assert CC.rfam_label_is_accession(types.SimpleNamespace(kind="c4_rfam_label_at", id="RF00005"))
    assert CC.rfam_label_is_accession(types.SimpleNamespace(kind="c4_rfam_label_at", id="tRNA")) is None
    assert CC.rfam_label_is_accession(types.SimpleNamespace(kind="c4_rfam_acc_at", id="RF00005")) is None
    names = [t.__name__ for t in K.TESTS]
    assert names.index("t_accession_version_never_had") < names.index("t_entrezgene_label_rows")
    assert K.CAUSE_CLASS["accession_version_never_had"][0] == "a" and K.CAUSE_CLASS["rfam_label_is_accession"][0] == "d"
    assert "c6" in K.DB_KINDS
