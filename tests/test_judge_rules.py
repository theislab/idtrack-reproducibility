"""The judge following idtrack's decided rules (revision sets of judge `rules`) and the correction `copies_apart_early`,
on invented lines: each revises exactly the lines its rule covers."""
from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as MQ


class FakeExceptions:
    """Ensembl's assembly exceptions of one release: Y's PAR equals X's segment; a patch of chromosome 6."""

    def __init__(self):
        self.r = {"Y": [dict(seqname="Y", start="10001", end="2781479", type="PAR", ref_seqname="X", ref_start="10001",
                             ref_end="2781479")],
                  "HSCHR6_MHC_COX": [dict(seqname="HSCHR6_MHC_COX", start="1", end="4795371", type="HAP",
                                          ref_seqname="6", ref_start="28477797", ref_end="33448354")]}

    def rows(self, rel):
        return self.r

    def project(self, place, rel):
        if place["seqname"] != "Y":
            return place
        return dict(place, seqname="X", projected_from="Y", exception="PAR")


class FakeCtx:
    def __init__(self, places, display=None, tables=None):
        self.places, self._display, self.tables, self._c = places, display or {}, tables or {}, {}

    def _cached(self, key, make):
        if key == "exceptions":
            return FakeExceptions()
        if key not in self._c:
            self._c[key] = make()
        return self._c[key]

    def gtf(self, rel):
        return self.places

    def primary(self, rel):
        return set(self.places)

    def display(self, rel):
        n2g = {}
        for g, (name, _src) in self._display.items():
            n2g.setdefault(name, set()).add(g)
        return n2g, self._display

    def bm(self, rel, group, col=1):
        return self.tables.get(group)

    def name_genes(self, rel, v, synonyms=False):
        out = set()
        for grp, (_, v2g) in self.tables.items():
            if grp != "synonym" or synonyms:
                out |= v2g.get(v, set())
        return out


def place(seq, start, end, strand="+"):
    return dict(seqname=seq, start=str(start), end=str(end), strand=strand, region="primary")


def line(qid, ident, expect, rule="set_exact", **kw):
    return dict(qid=qid, kind=kw.pop("kind", "name_ambiguous_to_set"), id=ident, from_release=kw.pop("from_release", 115),
                to_release=kw.pop("to_release", 115), final_database=None, strategy="all", rule=rule,
                cmp=kw.pop("cmp", "versioned"), expect=expect, broad=None, note=None, **kw)


def test_a_y_copy_of_a_shared_name_goes_apart_beside_its_x_gene_only():
    places = {"ENSG00000000001": place("X", 300000, 320000), "ENSG00000000002": place("Y", 300000, 320000), "ENSG00000000003": place("Y", 300000, 320000),
              "ENSG00000000004": place("5", 1, 100)}
    lines = [line("a", "PLCXD1", ["ENSG00000000001.1", "ENSG00000000002.1"]),          # X gene and its Y copy: the copy apart
             line("b", "OTHER", ["ENSG00000000003.1", "ENSG00000000004.1"]),           # a Y gene beside no X gene it repeats: unchanged
             line("c", "ENSG00000000001", ["ENSG00000000001.1", "ENSG00000000002.1"])]           # asked by an Ensembl id, not a name: unchanged
    out = MQ.rev_par_shared_names("homo_sapiens", lines, FakeCtx(places))
    got = {r["qid"]: r for v in out.values() for r in v}
    assert set(got) == {"a"}
    assert got["a"]["expect"] == ["ENSG00000000001.1"] and got["a"]["apart"] == ["ENSG00000000002.1"] and got["a"]["rule"] == "set"
    assert got["a"]["revised_by"] == ["par_shared_names"]


def test_copies_before_82_go_apart_by_the_assembly_exceptions():
    places = {"ENSG00000000001": place("6", 30000000, 30010000), "ENSG00000000002": place("HSCHR6_MHC_COX", 2000, 12000),
              "ENSG00000000003": place("Y", 300000, 320000)}
    lines = [line("a", "ENSG00000000009", ["ENSG00000000001", "ENSG00000000002"], to_release=60, cmp="base"),
             line("b", "ENSG00000000009", ["ENSG00000000001", "ENSG00000000002"], to_release=90),    # from 82 on rule `judge_copies_apart` reads the GTF
             line("c", "ENSG00000000009", ["ENSG00000000001", "ENSG00000000003"], to_release=60)]    # a PAR is no copy here
    out = MQ.rev_copies_apart_early("homo_sapiens", lines, FakeCtx(places))
    got = {r["qid"]: r for v in out.values() for r in v}
    assert set(got) == {"a"} and got["a"]["expect"] == ["ENSG00000000001"] and got["a"]["apart"] == ["ENSG00000000002"]


def test_missing_values_are_pandas_list_and_the_rule_s_two():
    for s in ("N/A", " nan ", "NULL", "<NA>", "#N/A", "NaT", "unknown", ""):
        assert MQ.is_missing_value(s), s
    for s in ("NA1", "TP53", "-", "nil"):
        assert not MQ.is_missing_value(s), s
    lines = [line("a", "N/A", ["ENSG1.1"], kind="synonym_to_gene"), line("b", "TP53", ["ENSG1.1"], kind="synonym_to_gene")]
    out = MQ.rev_missing_values("sus_scrofa", lines, None)
    got = {r["qid"]: r for v in out.values() for r in v}
    assert set(got) == {"a"} and "missing value" in got["a"]["excluded"]


def test_stray_punctuation_is_read_off_the_ends_only():
    assert MQ.strip_stray_punctuation("A1BG;") == "A1BG" and MQ.strip_stray_punctuation("'AGT'") == "AGT"
    assert MQ.strip_stray_punctuation("“TP53”.") == "TP53" and MQ.strip_stray_punctuation("HSP70B'") == "HSP70B'"
    assert MQ.strip_stray_punctuation("14-3-3") == "14-3-3" and MQ.strip_stray_punctuation(";") == ";"
    syn = ({}, {"alpha-1 antiproteinase,": {"G1", "G2"}, "alpha-1 antiproteinase": {"G2", "G3"}})
    ctx = FakeCtx({g: place("12", 1, 2) for g in ("G1", "G2", "G3")}, tables={"synonym": syn})
    lines = [dict(line("a", "alpha-1 antiproteinase,", ["G1", "G2"], kind="synonym_to_gene", rule="set"), cmp="base")]
    out = MQ.rev_stray_punctuation("mus_musculus", lines, ctx)
    (r,) = [r for v in out.values() for r in v]
    assert r["expect"] == ["G1", "G2", "G3"] and r["rule"] == "set_exact"


def test_one_case_strings_read_as_the_mixed_case_authority_s_symbol():
    display = {"G1": ("Narf", "MGI Symbol"), "G2": ("Acad8", "MGI Symbol"), "G3": ("Mup1", "MGI Symbol")}
    mgi = ({"G1": {"Narf"}, "G2": {"Acad8"}, "G3": {"Mup1"}}, {"Narf": {"G1"}, "Acad8": {"G2"}, "Mup1": {"G3"}})
    ctx = FakeCtx({g: place("1", 1, 2) for g in display}, display=display, tables={"mgi": mgi})
    assert MQ.mixed_case_symbols(ctx, 115) is mgi
    lines = [dict(line("a", "narf", ["G2.1"], kind="synonym_to_gene", rule="set"), cmp="base"),
             dict(line("b", "NARF", ["G2.1"], kind="synonym_to_gene", rule="set"), cmp="base")]
    low = {r["qid"]: r for v in MQ.rev_lower_case_symbols("mus_musculus", lines, ctx).values() for r in v}
    up = {r["qid"]: r for v in MQ.rev_capitals_mgi("mus_musculus", lines, ctx).values() for r in v}
    assert set(low) == {"a"} and low["a"]["expect"] == ["G1"] and "G2.1" in low["a"]["broad"]
    assert set(up) == {"b"} and up["b"]["expect"] == ["G1"]
    capitals = FakeCtx(ctx.places, display={g: (n.upper(), "HGNC Symbol") for g, (n, _) in display.items()},
                       tables={"hgnc": ({}, {n.upper(): gs for n, gs in mgi[1].items()})})
    assert MQ.mixed_case_symbols(capitals, 115) is None  # an authority writing capitals: the rule does not apply


def test_the_published_judge_leaves_out_every_revised_judge():
    assert set(H.REVISED_JUDGES) == {"78", "rules"} and MQ.RULES_JUDGE <= set(MQ.REVISIONS)
    assert "copies_apart_early" in MQ.REVISIONS and "copies_apart_early" not in MQ.RULES_JUDGE
    # the rules built after next's first round are published after every earlier revision, as judge `rules`; the
    # support of a name's sources alone (a correction of the corrected judge, `78`) after them
    assert list(MQ.REVISIONS)[-3:] == ["no_mito_nuclear_step", "mirbase_genes", "support_names"]
    assert "support_names" not in MQ.RULES_JUDGE
    assert {"no_mito_nuclear_step", "mirbase_genes"} <= MQ.RULES_JUDGE


def test_the_support_of_a_name_counts_its_name_sources_and_of_a_synonym_its_synonyms():
    """`support_names`: a synonym is no source of a name (idtrack reads the names that carry a string, rule
    `names_before_synonyms`); a string read as a synonym is counted by its synonyms."""
    class C:
        T = 115

        @staticmethod
        def support(rel, gene, v, synonyms=False):
            return 2 + (1 if synonyms and gene == G2 else 0)  # two name sources each; G2 carries the string as a synonym too

    G1, G2 = "ENSG00000000001", "ENSG00000000002"

    def line(kind, ident, support=None):
        rec = dict(qid=f"{kind}|{ident}", kind=kind, id=ident, strategy="best", rule="set_exact", expect=[f"{G1}.1", f"{G2}.1"],
                   from_release=100, to_release=100)
        if support:
            rec["support"] = support
        return rec
    out = MQ.rev_support_names("homo_sapiens", [line("name_ambiguous_to_set", "N1", {f"{G1}.1": 2, f"{G2}.1": 3}),
                                                line("synonym_to_gene", "S1")], C)
    # the name: equal support without the synonym, so the line loses it; the synonym: its synonym source counts
    name = out["name_ambiguous_to_set"][0]
    assert "support" not in name and name["revised_by"] == ["support_names"]
    assert out["synonym_to_gene"][0]["support"] == {f"{G1}.1": 2, f"{G2}.1": 3}
    assert not MQ.rev_support_names("homo_sapiens", [name, out["synonym_to_gene"][0]], C)  # run again: nothing
