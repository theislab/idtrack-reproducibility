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
