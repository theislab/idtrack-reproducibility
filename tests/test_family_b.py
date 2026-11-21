"""Family B (names over time): the names it asks, the judges of B1/B2, of a string read without a release and of B8,
the spreadsheet forms, on invented names and genes."""
import unicodedata

from idtrack_reproducibility.families.b import gen_b as GB
from idtrack_reproducibility.families.b import gen_b2 as G2
from idtrack_reproducibility.families.b import report_b as RB
from idtrack_reproducibility.families.b import truth_b as TB

T = TB.T


class History:
    """Versions and comparisons of an invented history: every gene has version 1."""

    def vid(self, g, R):
        return f"{g}.1"

    def vset(self, genes, R):
        return sorted(f"{g}.1" for g in genes)

    def cmp_for(self, R):
        return "versioned"

    def locus(self, g, R):
        return None


class Names:
    """The parts of `truth_b.NamesB` the generators read: display names per release, walks per (gene, release); no
    other name source, no synonym, no support."""

    org = "sus_scrofa"

    def __init__(self, display, walks=None, regions=None, marks=False, rels=None, ever=None, by_gene=None):
        self.tr = History()
        self._display, self._walks, self._regions = display, walks or {}, regions or {}
        self.marks, self.rels, self._ever, self._by_gene = marks, rels or sorted(display), ever or {}, by_gene or {}
        self.bm_rels = set()
        self._walk = {}

    def display(self, R):
        return self._display[R]

    def genes(self, R):
        return {g: dict(region=self._regions.get(g, "primary")) for gs in self._display[R].values() for g in gs}

    def gset(self, R):
        return frozenset(self.genes(R))

    def sources(self, R):
        by_gene = self._by_gene.get(R, {})
        return {}, "core_xref", {}, by_gene

    def synonyms(self, R):
        return {}

    def marks_patches(self, R):
        return self.marks

    def region(self, g, R):
        return self._regions.get(g, "primary")

    def ever(self, kind="all"):
        return self._ever

    def support(self, R, g, v):
        return None

    def walk(self, g, R, b=T):
        return self._walks[(g, R)]

    def walks(self, genes, R, b=T):
        s, br = set(), set()
        for g in genes:
            a, c = self.walk(g, R)
            s |= a
            br |= c
        return s, br


def test_only_name_strings_are_asked_as_names():
    assert TB.clean("TP53") and TB.clean("Hla-a") and TB.clean("MT-CO1")
    for v in ("", " TP53", "TP53 ", "ENSG00000141510", "ENSSSCG00000017950.5", "7157"):
        assert not TB.clean(v)


def test_copies_on_patches_are_apart_only_where_the_release_marks_them():
    nb = Names({100: {}}, regions={"P1": "primary", "H1": "patch", "H2": "patch"}, marks=True)
    assert GB.apart_split(nb, {"P1", "H1"}, 100) == ({"P1"}, {"H1"})
    assert GB.apart_split(nb, {"H1", "H2"}, 100) == ({"H1", "H2"}, set())  # only copies carry it: they are the answer
    nb.marks = False
    assert GB.apart_split(nb, {"P1", "H1"}, 100) == ({"P1", "H1"}, set())
    assert GB.rule_of(1) == "set" and GB.rule_of(2) == "set_exact"


def test_a_name_another_gene_carries_now_is_a_moved_name_judged_by_the_old_genes_history():
    # release 100: ABC is gene A's name, XYZ gene C's; at 115 gene B (no descendant of A) carries ABC, C2 (C's) XYZ
    nb = Names({100: {"ABC": {"A"}, "XYZ": {"C"}, "SHR": {"D", "E"}},
                T: {"ABC": {"B"}, "XYZ": {"C2"}, "SHR": {"D2", "E2"}, "OTHER": {"A2"}}},
               walks={("A", 100): ({"A2"}, {"A2"}), ("C", 100): ({"C2"}, {"C2"}), ("D", 100): ({"D2"}, {"D2"}),
                      ("E", 100): ({"E2"}, {"E2"})})
    g = GB.GenB1(nb)
    at, to = g.string_lines("ABC", 100)
    assert (at["qid"], at["stratum"], at["expect"], at["rule"]) == ("b1|ABC|100|100", "b1_atR_unique", ["A.1"], "set")
    assert (to["qid"], to["stratum"], to["expect"]) == ("b1|ABC|100|115", "b2_moved_from", ["A2.1"])
    assert g.moved == {"ABC"}
    _, to = g.string_lines("XYZ", 100)
    assert to["stratum"] == "b1_to115_unique" and to["expect"] == ["C2.1"] and g.moved == {"ABC"}
    at, to = g.string_lines("SHR", 100)  # a name of two genes (B3): every view must hold both
    assert at["stratum"] == "b1_atR_shared" and at["rule"] == "set_exact" and at["expect"] == ["D.1", "E.1"]
    assert to["rule"] == "set_exact" and to["expect"] == ["D2.1", "E2.1"]
    assert g.string_lines("NONE", 100) == []


def test_a_string_without_a_release_is_read_names_before_synonyms_else_open_within_its_genes_walks():
    nb = Names({100: {"OLD1": {"G9"}}, T: {"NAME1": {"A"}, "G9B": {"G9b"}}}, rels=[100, T],
               walks={("G9", 100): ({"G9b"}, {"G9b", "G9c"})}, ever={"OLD1": {"G9"}})
    cx = object.__new__(G2.Ctx)
    cx.nb = nb
    cx.GT = {"A", "B", "C", "G9b"}
    cx.disp_T, cx.src_T, cx.syn_T = {"NAME1": {"A"}}, {}, {"NAME1": {"B"}, "SYN1": {"C"}}
    cx.lower_T = {"name1": {"A", "B"}, "syn1": {"C"}}
    r = cx.default_line("q1", "b4_prev_default", "NAME1", "NAME1", authority={"A", "Z"}, auth_db="HGNC")
    assert r["expect"] == ["A.1"] and r["broad"] == ["A.1", "B.1"] and r["rule"] == "set"  # the name, not the synonym
    assert r["authority"] == ["A.1"] and r["authority_db"] == "HGNC"  # an authority's gene absent at 115 is left out
    assert cx.default_line("q2", "s", "SYN1", "SYN1")["expect"] == ["C.1"]
    r = cx.default_line("q3", "s", "OLD1", "OLD1")  # no name or synonym of 115: open within its genes' walks
    assert r["expect"] == [] and r["open"] and r["open_within"] == ["G9b.1"]
    r = cx.default_line("q4", "s", "NOWHERE", "NOWHERE")
    assert r["expect"] == [] and "open" not in r and r["note"]


def test_names_to_their_database_keep_both_spellings_of_one_value_and_skip_releases_without_it():
    nb = Names({100: {"Hsp90ab1": {"G1"}, "ABC": {"G2"}, "ENSSSCG00000000003": {"G3"}}, 101: {"ABC": {"G2"}}},
               by_gene={100: {"G1": {"Uniprot_gn": {"HSP90AB1", "Hsp90ab1"}}, "G2": {"Uniprot_gn": {"ABC"}}}, 101: {}})
    lines = list(GB.b8_lines(nb))
    assert [x["qid"] for x in lines] == ["b8|G1|100|UniProtKB Gene Name", "b8|G2|100|UniProtKB Gene Name"]
    g1, g2 = lines
    assert g1["expect"] == ["HSP90AB1", "Hsp90ab1"] and g1["spellings"] == [["HSP90AB1", "Hsp90ab1"]]
    assert "spellings" not in g2 and g2["cmp"] == "exact" and g2["stratum"] == "b8_core_xref" and g2["id"] == "G2.1"


def test_spreadsheet_forms_of_a_date_like_symbol():
    forms = {how: f for f, how in G2.excel_forms(9, 7)}
    assert forms == {"d-Mon": "7-Sep", "dd-Mon": "07-Sep", "Mon-dd": "Sep-07", "iso": "2007-09-07", "us": "9/7/2007",
                     "serial": "39332"}
    assert {how for _, how in G2.excel_forms(2, 29)} == {"d-Mon", "dd-Mon", "Mon-dd"}  # 2007 has no 29 February
    for s in ("SEPT7", "Sept7", "MARCH1", "DEC1", "OCT4", "SEP-15"):
        assert G2.DATE_RE.match(s)
    for s in ("SEPTIN7", "MARCHF1", "DECR1", "OCT"):
        assert not G2.DATE_RE.match(s)


def test_full_width_forms_fold_back_to_the_name():
    v = "HLA-A1"
    form = "".join(G2.FULLWIDTH.get(c, c) for c in v)
    assert form != v and all(0xFF00 < ord(c) < 0xFF5F for c in form)
    assert unicodedata.normalize("NFKC", form) == v  # what rule `unicode_normalised` reads (b6_dec)


def test_names_of_several_genes_are_reported_as_b3():
    assert RB.group_of("b1_atR_shared") == RB.group_of("b1_default_shared") == "b3"
    assert RB.group_of("b1_atR_unique") == "b1" and RB.group_of("b2_moved_from") == "b2"
    assert RB.group_of("b8_core_xref") == "b8" and RB.group_of("b5_alias_ncbi") == "b5"


# ---------------- B9: display names at their own release (gen_b9) ----------------

def _core(truth_store, sub, R, header, rows):
    from tests.conftest import write_gz
    write_gz(str(truth_store / sub / f"{R}.tsv.gz"), header, rows)


DISPLAY_COLS = ["gene_id", "db_name", "db_display_name", "acc", "label", "info_type"]
XREF_COLS = ["object_type", "object_id", "gene_id", "db_name", "db_display_name", "acc", "label", "version", "info_type"]


def test_the_era_starts_where_every_display_database_is_one_of_today():
    from idtrack_reproducibility.families.b import gen_b9 as B9
    assert B9.era_boundary({1: {"HGNC_curated_gene"}, 2: {"HGNC"}, 3: {"HGNC", "RFAM"}}, 3) == 2
    # a database of the past in the middle of the window keeps every release before it early
    assert B9.era_boundary({1: {"HGNC"}, 2: {"Uniprot/SPTREMBL"}, 3: {"HGNC"}}, 3) == 3
    assert B9.db_code("Uniprot/SPTREMBL") == "uniprot_sptrembl" and B9.db_code("HGNC_curated_gene") == "hgnc_curated_gene"


def test_b9_asks_each_display_name_at_its_release_by_its_filing_database(truth_store):
    from idtrack_reproducibility.families.b import gen_b9 as B9
    G1, G2, G3, G9 = "ENSSSCG00000000001", "ENSSSCG00000000002", "ENSSSCG00000000003", "ENSSSCG00000000009"
    _core(truth_store, "core_display", 1, DISPLAY_COLS, [(f"{G1}.1", "HGNC_curated_gene", "HGNC", "1", "GA", "PROJ"),
                                                         (f"{G3}.1", "HGNC", "HGNC", "3", "GC", "DIRECT"),
                                                         ("ENSSSCG00000000077.1", "HGNC", "HGNC", "7", "GX", "DIRECT")])
    _core(truth_store, "core_display", 2, DISPLAY_COLS, [(f"{G1}.1", "HGNC", "HGNC", "1", "GA", "DIRECT"),
                                                         (f"{G9}.1", "HGNC", "HGNC", "1", "GA", "DIRECT"),
                                                         (f"{G2}.1", "HGNC", "HGNC", "2", "ENSSSCG00000000002", "DIRECT")])
    _core(truth_store, "core_display", 3, DISPLAY_COLS, [(f"{G2}.2", "HGNC", "HGNC", "2", "GB", "DIRECT")])
    # at 1 another gene carries GA as an NCBI label and GC as a synonym: both go to broad, never to expect
    _core(truth_store, "core_xref", 1, XREF_COLS, [("Gene", f"{G9}.1", f"{G9}.1", "EntrezGene", "NCBI gene", "5", "GA", "0", "DEPENDENT")])
    _core(truth_store, "core_synonyms", 1, ["synonym", "gene_id", "db_name", "label", "object_type"], [("GC", G9, "HGNC", "GB0", "Gene")])
    g = B9.GenB9("sus_scrofa")
    lines = {q["qid"]: q for q in g.lines()}
    assert g.era_from == 2 and g.era(1) == "early" and g.era(3) == "late"
    ga1 = lines["b9|GA|1"]
    assert ga1["stratum"] == "b9_early_hgnc_curated_gene" and ga1["expect"] == [G1] and ga1["broad"] == [G1, G9]
    assert ga1["from_release"] == ga1["to_release"] == 1 and ga1["rule"] == "set" and ga1["cmp"] == "base"
    assert ga1["strategy"] == "both" and ga1["db"] == "HGNC_curated_gene"
    assert lines["b9|GC|1"]["broad"] == [G3, G9]
    # a name two genes carry at once: one line, every carrier expected (rule `set`: any part under best, all under all)
    ga2 = lines["b9|GA|2"]
    assert ga2["stratum"] == "b9_shared_late" and ga2["expect"] == [G1, G9]
    assert "b9|GB|3" in lines and lines["b9|GB|3"]["stratum"] == "b9_late_hgnc"
    # a gene's own id as its label is no name; a carrier the GTF lacks is no question
    assert not any(q.startswith("b9|ENSSSCG") for q in lines) and "b9|GX|1" not in lines
    assert g.excluded["display carriers not in the GTF of the release"] == 1


def test_b9_puts_patch_copies_apart_and_keeps_a_name_only_copies_carry(truth_store):
    from idtrack_reproducibility.families.b import gen_b9 as B9
    from tests.conftest import GENE_COLS, write_gz
    write_gz(str(truth_store / "gtf" / "genes_3.tsv.gz"), GENE_COLS, [
        ("ENSSSCG00000000002", 2, "GB", "protein_coding", "1", 1000, 1250, "+", "primary"),
        ("ENSSSCG00000000004", 1, "GB", "protein_coding", "CHR_P1", 1000, 1250, "+", "patch"),
        ("ENSSSCG00000000006", 1, "GF", "protein_coding", "CHR_P2", 5, 50, "+", "patch")])
    _core(truth_store, "core_display", 3, DISPLAY_COLS, [("ENSSSCG00000000002.2", "HGNC", "HGNC", "2", "GB", "DIRECT"),
                                                         ("ENSSSCG00000000004.1", "HGNC", "HGNC", "2", "GB", "DIRECT"),
                                                         ("ENSSSCG00000000006.1", "HGNC", "HGNC", "6", "GF", "DIRECT")])
    lines = {q["qid"]: q for q in B9.GenB9("sus_scrofa").lines()}
    gb = lines["b9|GB|3"]
    assert gb["expect"] == ["ENSSSCG00000000002"] and gb["apart"] == ["ENSSSCG00000000004"]
    assert lines["b9|GF|3"]["stratum"] == "b9_late_patchonly" and lines["b9|GF|3"]["expect"] == ["ENSSSCG00000000006"]


def test_b9_puts_copies_apart_before_82_by_the_assembly_exceptions(truth_store):
    """Before release 82 the GTF marks no copy: a gene on a haplotype (Ensembl's assembly exceptions) that carries a
    name beside the reference gene goes apart (the corrected judge's `copies_apart_early`); the line is then a name of
    one gene, in its database's stratum."""
    from idtrack_reproducibility.families.b import gen_b9 as B9
    from tests.conftest import GENE_COLS, write_gz
    write_gz(str(truth_store / "gtf" / "genes_2.tsv.gz"), GENE_COLS, [
        ("ENSSSCG00000000001", 1, "GA", "protein_coding", "1", 100, 200, "+", "primary"),
        ("ENSSSCG00000000005", 1, "GA", "protein_coding", "HAP1", 100, 200, "+", "primary")])
    (truth_store / "exceptions").mkdir()
    with open(truth_store / "exceptions" / "2.tsv", "w") as fo:
        fo.write("seqname\tstart\tend\ttype\tref_seqname\tref_start\tref_end\tori\tseq_length\tref_length\tfrom_release\n")
        fo.write("HAP1\t50\t300\tHAP\t1\t50\t300\t1\t1000\t1000\t2\n")
    _core(truth_store, "core_display", 2, DISPLAY_COLS, [("ENSSSCG00000000001.1", "HGNC", "HGNC", "1", "GA", "DIRECT"),
                                                         ("ENSSSCG00000000005.1", "HGNC", "HGNC", "1", "GA", "DIRECT")])
    g = B9.GenB9("sus_scrofa", releases=[2])
    ga = {q["qid"]: q for q in g.lines()}["b9|GA|2"]
    assert ga["expect"] == ["ENSSSCG00000000001"] and ga["apart"] == ["ENSSSCG00000000005"]
    assert ga["stratum"].startswith("b9_") and "shared" not in ga["stratum"] and ga["revised_by"] == ["copies_apart_early"]


def test_b9_expects_the_gene_an_lrg_record_describes(truth_store):
    """Rule `lrg_answers_gene`: a display name that is an LRG record's expects the Ensembl gene Ensembl's core database
    links to the LRG gene (`ENS_LRG_gene`), not the LRG gene; without such a link the line keeps the LRG gene."""
    from idtrack_reproducibility.families.b import gen_b9 as B9
    from tests.conftest import GENE_COLS, write_gz
    write_gz(str(truth_store / "gtf" / "genes_2.tsv.gz"), GENE_COLS, [
        ("ENSSSCG00000000001", 1, "GA", "protein_coding", "1", 100, 200, "+", "primary"),
        ("LRG_1_g1", 1, "LRG_1", "LRG_gene", "LRG_1", 1, 100, "+", "primary"),
        ("LRG_2_g1", 1, "LRG_2", "LRG_gene", "LRG_2", 1, 100, "+", "primary")])
    _core(truth_store, "core_display", 2, DISPLAY_COLS, [("LRG_1_g1.1", "LRG", "LRG", "LRG_1", "LRG_1", "DIRECT"),
                                                         ("LRG_2_g1.1", "LRG", "LRG", "LRG_2", "LRG_2", "DIRECT")])
    _core(truth_store, "core_xref", 2, XREF_COLS, [("Gene", "ENSSSCG00000000001.1", "ENSSSCG00000000001.1", "ENS_LRG_gene",
                                                    "LRG display in Ensembl", "LRG_1_g1", "LRG_1", "0", "DIRECT")])
    g = B9.GenB9("sus_scrofa", releases=[2])
    by = {q["qid"]: q for q in g.lines()}
    assert by["b9|LRG_1|2"]["expect"] == ["ENSSSCG00000000001"] and "lrg_answers_gene" in by["b9|LRG_1|2"]["note"]
    assert by["b9|LRG_2|2"]["expect"] == ["LRG_2_g1"] and g.lrg_lines == 1


def test_a_revision_compares_the_truth_not_the_revisions_bookkeeping():
    from idtrack_reproducibility.families.b import revise_b as RB
    line = dict(qid="b9|GA|2", kind="b9_display_at_release", id="GA", from_release=2, to_release=2, expect=["G1"],
                rule="set", cmp="base", strategy="both")
    published = RB.tag(dict(line, set="b9_r1"), "r1", "copies apart before 82")
    assert RB.canon(published) == RB.canon(line)
    assert RB.canon(dict(line, expect=["G2"])) != RB.canon(published)
