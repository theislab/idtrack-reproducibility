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


def test_an_rfam_label_written_as_the_accession_widens_its_line_by_the_entry_s_genes():
    class Tr:
        genes = {100: {"G1": 1, "G2": 1, "G3": 1}}

        @staticmethod
        def vset(genes, rel):
            return [f"{g}.1" for g in sorted(genes)]

    class G:
        tr = Tr()

    class Car:
        @staticmethod
        def entry(rel):
            return {"RF00005": {"G1", "G2", "G3"}, "RF00001": {"G1"}}

    def line(qid, kind, ident, expect, **kw):
        return dict(qid=qid, kind=kind, id=ident, from_release=100, to_release=100, expect=expect, broad=None, **kw)
    lines = [line("1", "c4_rfam_label_at", "RF00005", ["G1.1"]), line("2", "c4_rfam_label_at", "tRNA", ["G1.1"]),
             line("3", "c4_rfam_acc_at", "RF00005", ["G1.1"]), line("4", "c7_rfam_label", "RF00001", ["G1.1"]),
             line("5", "c4_rfam_label_retired", "RF00005", ["G1.1"], revised_by=["rfam_label_accession"])]
    out = list(RC.rfam_label_accession("sus_scrofa", lines, {"car": Car(), "G": G()}))
    assert [r["qid"] for r in out] == ["1"]  # line 4: the entry has no other gene; 2, 3: no label written as accession
    assert out[0]["expect"] == ["G1.1"] and out[0]["broad"] == ["G1.1", "G2.1", "G3.1"]
    assert out[0]["revised_by"] == ["rfam_label_accession"] and lines[0]["broad"] is None


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


def test_a_value_another_database_gives_at_the_target_widens_a_line_asked_without_a_release():
    class Tr:
        genes = {90: {"G1": 1}, 115: {"G2": 1, "G3": 1}}

        @staticmethod
        def vset(genes, rel):
            return [f"{g}.1" for g in sorted(genes)]

    class G:
        tr = Tr()

        @staticmethod
        def walks(g, R, targets):
            return {t: {"G2"} for t in targets}, {t: set() for t in targets}

    class Car:
        @staticmethod
        def at(rel):
            return {"2846624": {"G1"}} if rel == 90 else {"2846624": {"G3"}}

    lines = [dict(qid="1", kind="c2_retired", stratum="c2_retired", id="2846624", from_release=None, read_at=90,
                  to_release=115, expect=["G2.1"], broad=None, rule="set"),
             dict(qid="2", kind="c2_at", stratum="c2_at", id="2846624", from_release=90, to_release=90,
                  expect=["G1.1"], broad=None, rule="set")]
    out = {r["qid"]: r for r in RC.other_database_broad("macaca_mulatta", lines, {"G": G(), "car": Car()})}
    assert set(out) == {"1"} and out["1"]["broad"] == ["G2.1", "G3.1"] and out["1"]["expect"] == ["G2.1"]


def test_an_accession_with_a_version_is_not_read_at_the_target():
    class Tr:
        genes = {90: {"G1": 1}, 115: {"G2": 1, "G3": 1}}

        @staticmethod
        def vset(genes, rel):
            return [f"{g}.1" for g in sorted(genes)]

    class G:
        tr = Tr()

        @staticmethod
        def walks(g, R, targets):
            return {t: {"G2"} for t in targets}, {t: set() for t in targets}

    class Car:
        @staticmethod
        def at(rel):
            return {} if rel == 90 else {"NM_1": {"G3"}}

    line = dict(qid="1", kind="c2_oldver", stratum="c2_oldver", id="NM_1.2", from_release=None, read_at=90,
                to_release=115, expect=["G2.1"], broad=None, rule="set")
    assert list(RC.other_database_broad("x", [line], {"G": G(), "car": Car()})) == []


def test_a_set_revised_from_a_published_one_compares_the_truth_not_its_wording():
    from idtrack_reproducibility.families import lines as L
    from idtrack_reproducibility.families.c import gen_c as G
    line = dict(qid="c2_ver_gap|NM_1.3|-|115", kind="c2_ver_gap", id="NM_1.3", from_release=None, to_release=115,
                final_database=None, strategy="both", rule="set", cmp="versioned", expect=[], broad=None,
                note="a version never had", q2="nothing (an older wording of the rule)", read_at=115)
    same = dict(line, note="another wording", q2="nothing (rule `version_newer_than_graph`'s rule)")
    other = dict(line, read_at=89)
    assert G.truth_hash(L.to_row(line)) == G.truth_hash(L.to_row(same))
    assert G.truth_hash(L.to_row(line)) != G.truth_hash(L.to_row(other))
    assert G.truth_hash(L.to_row(line)) != G.truth_hash(L.to_row(dict(line, expect=["ENSMUSG00000000001.1"])))


def test_a_copy_the_line_expects_apart_is_no_other_carrier():
    class Tr:
        genes = {75: {"G1": 1, "G9": 1}}

        @staticmethod
        def vset(genes, rel):
            return [f"{g}.1" for g in sorted(genes)]

    class G:
        tr = Tr()

    class Car:
        @staticmethod
        def at(rel):
            return {"100039209": {"G1", "G9"}}  # G9: the copy on a patch, carrying the same NCBI id

    line = dict(qid="1", kind="c1_at", stratum="c1_at", id="100039209", from_release=75, to_release=75,
                expect=["G1.1"], apart=["G9.1"], broad=None, rule="set")
    assert list(RC.other_database_broad("x", [line], {"G": G(), "car": Car()})) == []
    other = dict(line, apart=None)
    assert [r["broad"] for r in RC.other_database_broad("x", [other], {"G": G(), "car": Car()})] == [["G1.1", "G9.1"]]


def test_a_database_whose_version_changes_with_most_releases_counts_edits(tmp_path, monkeypatch):
    """Rule `accession_versions`' edit counters measured from Ensembl's records: the share of consecutive-release pairs of one link (accession,
    object with version, database) in which the version differs; a version 0 is none, a gap breaks the pair, and the
    bound is strict (one half is not above one half)."""
    monkeypatch.setattr(TC, "TRUTH", str(tmp_path))
    sp, tr, rs = "Uniprot/SWISSPROT", "Uniprot/SPTREMBL", "RefSeq_peptide"
    rows = {
        100: [("Translation", "P1.1", "G1", sp, "", "P00001", "P00001", 10, "X"), ("Translation", "P1.1", "G1", rs, "", "NP_1", "NP_1.1", 1, "X"),
              ("Translation", "P2.1", "G2", tr, "", "A00002", "A00002", 5, "X"), ("Translation", "P3.1", "G3", tr, "", "A00003", "A00003", 0, "X")],
        101: [("Translation", "P1.1", "G1", sp, "", "P00001", "P00001.11", 11, "X"), ("Translation", "P1.1", "G1", rs, "", "NP_1", "NP_1.1", 1, "X"),
              ("Translation", "P3.1", "G3", tr, "", "A00003", "A00003", 0, "X")],
        102: [("Translation", "P1.1", "G1", sp, "", "P00001", "P00001.12", 12, "X"), ("Translation", "P1.1", "G1", rs, "", "NP_1", "NP_1.2", 2, "X"),
              ("Translation", "P2.1", "G2", tr, "", "A00002", "A00002.7", 7, "X")],
        103: [("Translation", "P1.1", "G1", sp, "", "P00001", "P00001.12", 12, "X"), ("Translation", "P1.1", "G1", rs, "", "NP_1", "NP_1.2", 2, "X"),
              ("Translation", "P2.1", "G2", tr, "", "A00002", "A00002.7", 7, "X")],
        116: [("Translation", "P1.1", "G1", rs, "", "NP_1", "NP_1.9", 9, "X")],
    }
    for rel, rr in rows.items():
        write_xref(tmp_path, "sus_scrofa", rel, rr)
    shares = TC.version_change_shares("sus_scrofa", {sp, tr, rs}, upto=115)
    # Swiss-Prot 10 -> 11 -> 12 -> 12: 2 of 3; TrEMBL A00002: no pair across the gap at 101, 102 -> 103 the same: 0 of 1
    # (A00003's version 0 is none); RefSeq 1 -> 1 -> 2 -> 2: 1 of 3 (release 116 lies beyond the window)
    assert shares == {sp: (3, 2), tr: (1, 0), rs: (3, 1)}
    assert TC.edit_counter_dbs(shares) == {sp}
    assert TC.edit_counter_dbs({sp: (4, 2)}) == set()  # exactly one half: not above it


class _GV:
    """C's generator over the versioned databases, invented: values, their genes and versions per release."""

    def __init__(self, H):
        X = object.__new__(TC.Xrefs)
        X.H = H
        self.X, self.T = X, 115
        self.tr = types.SimpleNamespace(cmp_for=lambda R: "versioned")

    def at_release(self, gs, R):
        return dict(cmp="versioned", expect=sorted(gs))

    def walked(self, gs, R, t):
        return dict(cmp="versioned", expect=sorted(f"{g}w" for g in gs))

    def default_truth(self, fam, v):
        return GC.GenC.default_truth(self, fam, v)


def _ver_line(stratum, ident, expect, fr=None, to=115, **kw):
    return dict(qid=f"{stratum}|{ident}|{fr or '-'}|{to}", stratum=stratum, kind=stratum, id=ident, from_release=fr,
                to_release=to, expect=expect, rule="set", cmp="versioned", **kw)


def _held(gs, vs):
    return frozenset(gs), frozenset(vs), True


def test_a_typed_version_of_an_edit_counter_is_not_judged_and_reads_the_accession_s_own_line():
    gv = _GV({"uniprot": {"P00001": {100: _held({"G1"}, {"3"}), 115: _held({"G2"}, {"9"})},
                          "A00002": {100: _held({"G3"}, {"4"}), 105: _held({"G3"}, {"6"})}},
              "refseq": {"NM_1": {100: _held({"G4"}, {"1"})}}})
    # the accessions' own lines without a version, as a correction of another rule left them (a wider set)
    refs = {"c3_default|P00001|-|115": _ver_line("c3_default", "P00001", ["G2"], broad=["G2", "G7"]),
            "c3_retired|A00002|-|115": _ver_line("c3_retired", "A00002", ["G3w"])}
    ctx = {"GV": gv, "counters": set(TC.FAMILIES["uniprot"]["dbs"]), "refs": refs, "held": _Held({})}
    lines = [_ver_line("c3_ver_lower", "P00001.2", []), _ver_line("c3_oldver", "P00001.3", ["G1w"]),
             _ver_line("c3_ver_higher_retired", "A00002.7", []), _ver_line("c3_at_v", "P00001.3", ["G1"], fr=100, to=100),
             _ver_line("c3_ver_higher", "P00001.10", ["G2"], broad=["G2", "G7"]), _ver_line("c2_ver_lower", "NM_1.0", []),
             _ver_line("c3_default", "P00001", ["G2"])]
    out = {r["qid"]: r for r in RC.version_not_judged("sus_scrofa", lines, ctx)}
    assert set(out) == {"c3_ver_lower|P00001.2|-|115", "c3_oldver|P00001.3|-|115", "c3_ver_higher_retired|A00002.7|-|115"}
    r = out["c3_ver_lower|P00001.2|-|115"]
    assert r["expect"] == ["G2"] and r["broad"] == ["G2", "G7"]  # the own line's truth, its correction included
    assert r["read_at"] == 115 and r["read_as"] == "c3_default|P00001|-|115" and r["version"].startswith("not judged")
    r = out["c3_ver_higher_retired|A00002.7|-|115"]
    assert r["expect"] == ["G3w"] and r["read_at"] == 105 and r["revised_by"] == ["version_not_judged"] and "q2" not in r
    assert not list(RC.version_not_judged("sus_scrofa", [dict(r)], ctx))  # a line revised once is left
    # RefSeq's numbers date sequences: not this correction's
    assert not list(RC.version_not_judged("sus_scrofa", [_ver_line("c2_ver_lower", "NM_1.0", [])], dict(ctx, counters=set())))


def test_a_version_above_every_version_had_is_the_latest_only_where_the_accession_is_attached_at_the_target():
    gv = _GV({"refseq": {"NM_1": {100: _held({"G1"}, {"2"}), 104: _held({"G1"}, {"3"}), 106: (frozenset({"G1"}), frozenset({"4"}), False)},
                         "NM_2": {110: _held({"G2"}, {"1"}), 113: _held({"G2"}, {"2"}), 115: _held({"G2"}, {"2"})},
                         "NM_3": {100: _held({"G1"}, {"2"}), 113: _held({"G1"}, {"3"}), 116: _held({"G1"}, {"5"})}},
              "uniprot": {"P00001": {100: _held({"G3"}, {"5"})}}})
    refs = {"c2_oldver|NM_1.3|-|115": _ver_line("c2_oldver", "NM_1.3", ["G1w"], read_at=104),
            "c2_default_v|NM_2.2|-|115": _ver_line("c2_default_v", "NM_2.2", ["G2"], other_carriers=["G9"])}
    ctx = {"GV": gv, "counters": set(TC.FAMILIES["uniprot"]["dbs"]), "refs": refs, "held": _Held({})}
    # NM_1 (versions 2 and 3; 4 at 106 is outside the graph's scope) is not attached at the target: 4 above its
    # versions is an unknown version, nothing (an earlier set's reading of the latest is put back); NM_2 is attached
    # at the target: 3 above its versions is the latest, read as NM_2.2; NM_3.4: version 5 lies at 116, beyond the
    # window, and NM_3 is not attached at 115: nothing
    lines = [_ver_line("c2_ver_higher_retired", "NM_1.4", ["G1w"], read_at=104, latest_version=3, revised_by=["version_above_latest"]),
             _ver_line("c2_ver_higher", "NM_2.3", ["G2"]), _ver_line("c3_ver_higher_retired", "P00001.6", []),
             _ver_line("c2_ver_gap", "NM_3.4", [])]
    out = {r["qid"]: r for r in RC.versions_in_window("sus_scrofa", lines, ctx)}
    assert set(out) == {"c2_ver_higher_retired|NM_1.4|-|115", "c2_ver_higher|NM_2.3|-|115"}
    back = out["c2_ver_higher_retired|NM_1.4|-|115"]
    assert back["expect"] == [] and back["read_at"] is None and "latest_version" not in back
    assert back["revised_by"] == ["version_above_latest", "versions_in_window"]
    r = out["c2_ver_higher|NM_2.3|-|115"]
    assert r["expect"] == ["G2"] and r["other_carriers"] == ["G9"] and r["read_at"] == 115 and r["latest_version"] == 2
    assert not list(RC.versions_in_window("sus_scrofa", list(out.values()), ctx))  # run again on its own output: nothing


def test_versions_beyond_the_graph_s_window_are_no_versions_the_accession_had():
    gv = _GV({"refseq": {"NM_3": {100: _held({"G1"}, {"2"}), 113: _held({"G1"}, {"3"}), 115: _held({"G1"}, {"3"}), 116: _held({"G1"}, {"5"})},
                         "NM_4": {100: _held({"G2"}, {"1"}), 115: _held({"G2"}, {"3"})},
                         "NM_5": {116: _held({"G3"}, {"1"})}}})
    refs = {"c2_default_v|NM_3.3|-|115": _ver_line("c2_default_v", "NM_3.3", ["G1"]),
            "c2_default_v|NM_4.3|-|115": _ver_line("c2_default_v", "NM_4.3", ["G2"], broad=["G2", "G8"])}
    ctx = {"GV": gv, "counters": set(TC.FAMILIES["uniprot"]["dbs"]), "refs": refs, "held": _Held({})}
    assert gv.X.versions_in_graph("refseq", "NM_3") == {"2", "3", "5"}
    assert gv.X.versions_in_graph("refseq", "NM_3", 115) == {"2", "3"}
    lines = [_ver_line("c2_ver_gap", "NM_3.4", []), _ver_line("c2_ver_gap", "NM_4.2", []),
             _ver_line("c2_default_v", "NM_3.3", ["G1"]), _ver_line("c2_default_v", "NM_4.3", ["G2"], broad=["G2", "G8"]),
             _ver_line("c2_ver_lower", "NM_5.0", []), _ver_line("c2_at_v", "NM_3.5", [], fr=116, to=116)]
    out = list(RC.versions_in_window("sus_scrofa", lines, ctx))
    # NM_3.4: version 5 lies at 116, beyond the window: 4 is above the versions 2 and 3 the window had, and NM_3 is
    # attached at the target: read as NM_3.3; a line read as itself keeps its truth, whatever corrections made it
    assert [r["qid"] for r in out] == ["c2_ver_gap|NM_3.4|-|115"]
    assert out[0]["expect"] == ["G1"] and out[0]["read_at"] == 115 and out[0]["latest_version"] == 3
    assert out[0]["revised_by"] == ["versions_in_window"]
    # run again on its own output: nothing; a line an earlier set changed without cause is put back
    assert not list(RC.versions_in_window("sus_scrofa", [out[0]], ctx))
    wrong = dict(_ver_line("c2_default_v", "NM_4.3", ["G2", "G8"]), revised_by=["versions_in_window"], latest_version=9)
    back = list(RC.versions_in_window("sus_scrofa", [wrong], ctx))
    assert len(back) == 1 and back[0]["expect"] == ["G2"] and back[0]["broad"] == ["G2", "G8"]
    assert "latest_version" not in back[0] and "read_as" not in back[0]


class _Held:
    """`revise_c.HeldEarly`, invented: {string: {release: genes whose display name it is there}}."""

    def __init__(self, at):
        self.at = at

    def add(self, values):
        pass


def test_a_versioned_string_that_is_a_display_name_is_read_as_that_name_where_that_is_later():
    gv = _GV({"refseq": {"XM_9": {70: _held({"G1"}, {"2"})}, "XM_8": {70: _held({"G2"}, {"2"}), 115: _held({"G2"}, {"3"})}}})
    refs = {"c2_oldver|XM_9.2|-|115": _ver_line("c2_oldver", "XM_9.2", ["G1w"], read_at=70),
            "c2_default_v|XM_8.3|-|115": _ver_line("c2_default_v", "XM_8.3", ["G2"])}
    held = _Held({"XM_9.1": {60: {"G5"}, 62: {"G5"}}, "XM_8.3": {60: {"G7"}}})
    ctx = {"GV": gv, "counters": set(TC.FAMILIES["uniprot"]["dbs"]), "refs": refs, "held": held}
    lines = [_ver_line("c2_ver_lower", "XM_9.1", []), _ver_line("c2_default_v", "XM_8.3", ["G2"])]
    out = list(RC.versions_in_window("sus_scrofa", lines, ctx))
    # XM_9.1: no version reading (below version 2), Ensembl's display name of G5 at 60 and 62: read at 62; XM_8.3: its
    # version reading at the target is later than the display name at 60
    assert [r["qid"] for r in out] == ["c2_ver_lower|XM_9.1|-|115"]
    assert out[0]["expect"] == ["G5w"] and out[0]["read_at"] == 62 and out[0]["read_as_name"] is True


def _filing_releases(tmp_path):
    """RefSeq's accessions as `RefSeq_dna` (on transcripts) at 62 and 63, as RefSeq_mRNA from 64; Rfam's accessions on
    transcripts at 62, on genes from 63; and a database (GO) whose links never move into an included one."""
    g1, g2, g3 = "ENSMUSG00000000001", "ENSMUSG00000000002", "ENSMUSG00000000003"
    dna = [("Transcript", 1, g1, "RefSeq_dna", "RefSeq DNA", "NM_000001", "NM_000001", 1, "DIRECT"),
           ("Transcript", 2, g2, "RefSeq_dna", "RefSeq DNA", "NM_000002", "NM_000002", 1, "DIRECT"),
           ("Transcript", 3, g3, "GO", "GO", "GO:0000001", "GO:0000001", 0, "DIRECT")]
    write_xref(tmp_path, "mus_musculus", 62, dna + [("Transcript", 1, g1, "RFAM", "RFAM", "RF00001", "5S_rRNA", 0, "DIRECT")])
    write_xref(tmp_path, "mus_musculus", 63, [r[:7] + (2,) + r[8:] if r[5] == "NM_000001" else r for r in dna] +
               [("Gene", 1, g1, "RFAM", "RFAM", "RF00001", "5S_rRNA", 0, "DIRECT")])
    write_xref(tmp_path, "mus_musculus", 64, [
        ("Transcript", 1, g1, "RefSeq_mRNA", "RefSeq mRNA", "NM_000001", "NM_000001", 2, "DIRECT"),
        ("Transcript", 2, g2, "RefSeq_mRNA", "RefSeq mRNA", "NM_000002", "NM_000002", 1, "DIRECT"),
        ("Transcript", 3, g3, "GO", "GO", "GO:0000001", "GO:0000001", 0, "DIRECT"),
        ("Gene", 1, g1, "RFAM", "RFAM", "RF00001", "5S_rRNA", 0, "DIRECT")])
    return Scope({"RefSeq_mRNA": {64}, "RFAM": {63, 64}}, {"RefSeq_mRNA": "transcript", "RFAM": "gene"})


def test_a_database_whose_links_move_into_an_included_one_continues_it(tmp_path, monkeypatch):
    monkeypatch.setattr(TC, "TRUTH", str(tmp_path))
    scope = _filing_releases(tmp_path)
    out, ev = TC.continuations("mus_musculus", scope, log=lambda *a: None)
    assert out[("RefSeq_dna", None)] == ("RefSeq_mRNA", frozenset({62, 63}))  # every release it has links
    assert out[("RFAM", "transcript")] == ("RFAM", frozenset({62}))  # an included name on another form: its run
    assert ("GO", None) not in out  # its links stay where they are
    X = TC.Xrefs("mus_musculus", scope, families=["refseq", "rfam"], log=lambda *a: None, continued=out)
    assert X.at("refseq", "NM_000001", 62) == (frozenset({"ENSMUSG00000000001"}), frozenset({"1"}), True)
    assert X.versions_in_graph("refseq", "NM_000001") == {"1", "2"}  # the version RefSeq_dna wrote is the accession's
    assert X.at("rfam", "RF00001", 62)[2] is True
    plain = TC.Xrefs("mus_musculus", scope, families=["refseq", "rfam"], log=lambda *a: None)
    assert plain.at("refseq", "NM_000001", 62) is None and plain.at("rfam", "RF00001", 62)[2] is False


def test_a_database_that_ends_hands_over_what_it_still_had(tmp_path, monkeypatch):
    """Row I157: few of the ending database's links move, but they are most of the receiver's links there."""
    monkeypatch.setattr(TC, "TRUTH", str(tmp_path))
    old = [("Transcript", i, f"ENSMUSG{i:011d}", "RefSeq_dna_predicted", "x", f"XM_{i:06d}", f"XM_{i:06d}", 1, "DIRECT")
           for i in range(1, 11)]
    write_xref(tmp_path, "mus_musculus", 63, old)
    write_xref(tmp_path, "mus_musculus", 64, [("Transcript", i, f"ENSMUSG{i:011d}", "RefSeq_mRNA_predicted", "x",
                                               f"XM_{i:06d}", f"XM_{i:06d}", 1, "DIRECT") for i in (1, 2, 3)])
    scope = Scope({"RefSeq_mRNA_predicted": {64}}, {"RefSeq_mRNA_predicted": "transcript"})
    out, _ = TC.continuations("mus_musculus", scope, log=lambda *a: None)
    assert out[("RefSeq_dna_predicted", None)][0] == "RefSeq_mRNA_predicted"  # 3 of 10 move: 3 of the receiver's 3


def test_a_value_held_now_through_a_continued_database_is_read_where_it_is_held():
    """Correction `database_continuity`: a not-held line of rule `safe` (and a retired line) read by rule
    `start_release` at the latest release where the continued reading holds the value; a line read where it was keeps
    its truth."""
    class G:
        def default_truth(self, fam, v):
            return {"RF00605": (57, dict(cmp="versioned", expect=["ENSSSCG00000028053.2"])),
                    "NM_000009": (70, dict(cmp="versioned", expect=["ENSSSCG00000000009.3"]))}.get(v, (None, None))
    lines = [dict(qid="c4_rfam_acc_not_held|RF00605|-|115", kind="c4_rfam_acc_not_held", stratum="c4_rfam_acc_not_held",
                  id="RF00605", from_release=None, to_release=115, rule="safe", cmp="versioned", expect=[], read_at=57),
             dict(qid="c2_retired|NM_000009|-|115", kind="c2_retired", stratum="c2_retired", id="NM_000009",
                  from_release=None, to_release=115, rule="set", cmp="versioned", expect=["ENSSSCG00000000009.3"], read_at=70),
             dict(qid="c4_rfam_acc_not_held|RF00001|-|115", kind="c4_rfam_acc_not_held", stratum="c4_rfam_acc_not_held",
                  id="RF00001", from_release=None, to_release=115, rule="safe", cmp="versioned", expect=[])]
    out = list(RC.database_continuity("sus_scrofa", lines, {"GK": G()}))
    assert [r["qid"] for r in out] == ["c4_rfam_acc_not_held|RF00605|-|115"]
    assert out[0]["rule"] == "set" and out[0]["read_at"] == 57 and out[0]["expect"] == ["ENSSSCG00000028053.2"]
    assert out[0]["revised_by"] == ["database_continuity"]


def test_a_ccds_accession_without_transcripts_has_no_authority():
    """CCDS's sequence table, keyed by accession (`truth_c.Authority.ccds_tx`, a defaultdict of sets): an accession
    it does not hold, with or without its version, gives no transcript (the c5 generation once joined two empty
    tuples with `|`)."""
    tx = collections.defaultdict(set, {"CCDS11118": {"ENST1"}, "CCDS11118.1": {"ENST2"}})
    assert GC.ccds_transcripts(tx, "CCDS11118.1") == {"ENST1", "ENST2"}
    assert GC.ccds_transcripts(tx, "CCDS11118") == {"ENST1"}
    assert GC.ccds_transcripts({}, "CCDS9.2") == set()
