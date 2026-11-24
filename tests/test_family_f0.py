"""F0, the common machinery: the query line in its typed form and back, views, variants and units."""
import pytest

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L


def test_a_line_survives_its_typed_row():
    line = dict(qid="a1|X|56|115", family="a", stratum="a1_versioned", kind="a1", id="ENSSSCG00000000001.3",
                from_release=56, to_release=115, final_database=None, strategy="both", rule="set", cmp="versioned",
                expect=["ENSSSCG00000000002.1"], broad=None, era="Sscrofa9", dist=59)
    row = L.to_row(line)
    assert row["extra"] == C.dumps({"era": "Sscrofa9", "dist": 59}) and row["h"] == C.crc(line["qid"])
    back = L.from_row(row)
    assert back["era"] == "Sscrofa9" and back["expect"] == line["expect"] and back["from_release"] == 56


def test_a_line_needs_a_string_id():
    with pytest.raises(ValueError):
        L.to_row(dict(qid="q", kind="k", id=5))


def test_views_variants_units():
    assert L.views(dict(strategy="both")) == ["best", "all"] and L.views(dict(strategy="all")) == ["all"]
    assert L.variant(dict(from_release=None), "best") == "best,default"
    assert L.unit_of(dict(id="ENSG00000000003.15")) == "ENSG00000000003"
    assert L.unit_of(dict(id="ENSG00000000003.15", unit="lineage7")) == "lineage7"
    assert L.judged_by(dict(revised_by=["places", "support"])) == "revised:places+support"


def test_shards_depend_on_the_qid_only():
    assert C.shard_of("q1", 7) == C.crc("q1") % 7
    assert C.base_id("LRG_12.1") == "LRG_12" and C.base_id("TP53") == "TP53"


class _Record:
    content_digest = "abc"
    package_commit = "0f88cf1"
    label_names = {"RFAM": "RFAM label"}


class _GraphFormat3:
    """A graph of format 3: its record, and the database names it records as renamed."""
    record = _Record()
    renamed_databases = {"NCBI gene": "NCBI gene (formerly Entrezgene)"}
    graph = {}


class _GraphBefore:
    """A graph of d4c2dda or 71690ab: the manifest in graph.graph, no record, no renamed names."""
    graph = {"manifest": {"content_digest": "def", "package_commit": "71690ab"}}


def test_the_adapter_reads_the_record_of_either_graph_format():
    from idtrack_reproducibility.families import adapter as A
    assert A.graph_manifest(_GraphFormat3()) == {"content_digest": "abc", "package_commit": "0f88cf1"}
    assert A.graph_manifest(_GraphBefore()) == {"content_digest": "def", "package_commit": "71690ab"}
    assert A.renamed_databases(_GraphFormat3()) == {"NCBI gene": "NCBI gene (formerly Entrezgene)"}
    assert A.renamed_databases(_GraphBefore()) == {}


def test_a_database_the_graph_records_as_renamed_is_asked_under_the_graph_s_name():
    from idtrack_reproducibility.families import adapter as A
    ad = A.Adapter.__new__(A.Adapter)  # no idtrack: only the name mapping
    ad.renamed = A.renamed_databases(_GraphFormat3())
    assert ad.database("NCBI gene") == "NCBI gene (formerly Entrezgene)"
    assert ad.database("EntrezGene") == "EntrezGene" and ad.database(None) is None
    ad.renamed = A.renamed_databases(_GraphBefore())
    assert ad.database("NCBI gene") == "NCBI gene"


def test_a_line_asking_a_database_s_labels_is_asked_under_the_graph_s_name_of_the_labels():
    """Rule `one_kind_per_database`: a graph that records a name of their own for a database's labels is asked under it
    by a line that says it asks the labels; a graph without the record, and a line that asks the accessions, under the
    database's name. Fails if a label line asks the accession rows on a graph that names its labels apart."""
    from idtrack_reproducibility.families import adapter as A
    ad = A.Adapter.__new__(A.Adapter)
    ad.renamed, ad.label_names = {}, A.label_names(_GraphFormat3())
    assert ad.from_database({"from_database": "RFAM", "from_database_rows": "labels"}) == "RFAM label"
    assert ad.from_database({"from_database": "RFAM"}) == "RFAM"
    ad.label_names = A.label_names(_GraphBefore())
    assert ad.label_names == {} and ad.from_database({"from_database": "RFAM", "from_database_rows": "labels"}) == "RFAM"


def test_lines_asked_again_win_over_their_earlier_answers(tmp_path, monkeypatch):
    """runner `--again`: a revision set's lines are asked again whatever the folder holds, into `again_<set>_*` parts that
    every reader takes after the shard parts. Fails if the stale answer wins or the lines are skipped as answered."""
    from types import SimpleNamespace
    from idtrack_reproducibility.families import runner as R
    from idtrack_reproducibility.families import scorer as S
    (tmp_path / "answers").mkdir()
    rows = [dict(qid=q, view="best", outcome="answer", got=["old"]) for q in ("a", "b")]
    R.write_part(str(tmp_path), 0, 1, rows)
    R.write_part(str(tmp_path), 0, 1, [dict(qid="a", view="best", outcome="answer", got=["new"])], "again_s_r1")
    files = [p.rsplit("/", 1)[1] for p in C.answer_files(str(tmp_path))]
    assert files == ["shard_0_of_1.0.parquet", "again_s_r1_0_of_1.0.parquet"]
    got = S.read_answers(str(tmp_path), 0, 1).set_index("qid")["got"]
    assert list(got["a"]) == ["new"] and list(got["b"]) == ["old"]
    monkeypatch.setattr(L, "_qid_column", lambda path: __import__("pyarrow").array(["a", "b"]))
    monkeypatch.setattr(L, "_set_file", lambda d, name: "x")
    monkeypatch.setattr(L, "query_dir", lambda f, o: "d")
    a = SimpleNamespace(out=str(tmp_path), again="s_r1", family="f", organism="o")
    todo, prefix = R.to_ask(a, [{"qid": "a"}, {"qid": "b"}, {"qid": "c"}], 0, 1)
    assert prefix == "again_s_r1" and [q["qid"] for q in todo] == ["b"]  # a: asked again already; c: not in the set
    todo, prefix = R.to_ask(SimpleNamespace(out=str(tmp_path), again=None), [{"qid": "a"}, {"qid": "c"}], 0, 1)
    assert prefix == "shard" and [q["qid"] for q in todo] == ["c"]


def test_causes_written_again_under_the_current_names(tmp_path):
    import gzip
    from idtrack_reproducibility.families import causes as K
    p = tmp_path / "causes.tsv.gz"
    with gzip.open(p, "wt") as fo:
        fo.write("organism\tqid\tview\tkind\tstratum\tunit\tcause\tclass\tdetail\n")
        fo.write("o\tq1\tbest\tk\ts1\tu1\texternal_bridge\tz\td\n")
        fo.write("o\tq2\tall\tk\ts1\tu2\texternal_bridge\tz\td\n")
    K.summary_only(str(tmp_path))
    with gzip.open(p, "rt") as fh:
        assert [ln.split("\t")[7] for ln in fh.read().splitlines()[1:]] == ["b", "b"]
    summ = (tmp_path / "causes_summary.tsv").read_text().splitlines()
    assert summ[1].split("\t") == ["s1", "external_bridge", "b", "rule `external_bridge`", "2", "2"]


def test_lrg_gene_ids_lose_their_version_as_ensembl_ids_do():
    from idtrack_reproducibility.harness import score as S
    assert C.base_id("LRG_10_g1.1") == "LRG_10_g1" and C.base_id("LRG_1t1.2") == "LRG_1t1"
    assert C.base_id("LRG_12.1") == "LRG_12" and C.base_id("ENSG00000141510.18") == "ENSG00000141510"
    assert C.base_id("LRG1.2") == "LRG1.2"  # a gene name, not an LRG id
    assert S.norm(["LRG_10_g1.1"], "base") == {"LRG_10_g1"}


def test_a_wrong_answer_under_a_rule_decided_after_the_tree_has_its_cause():
    from idtrack_reproducibility.families import causes as K
    c = K.Ctx(dict(kind="name_ambiguous_to_set", query='{"id": "CSF2RA", "to_release": 115}', got="ENSG00000198223.1,ENSG00000292357.1",
                   expected="ENSG00000198223.1", rule="revised:places+par_shared_names",
                   _before=dict(expect=["ENSG00000198223.1", "ENSG00000292357.1"])), None, None)
    assert "par_shared_names" in K.t_rule_decided_after_tree(c)
    c.r["rule"] = "revised:places"  # a correction of the judge, not a decided rule
    assert K.t_rule_decided_after_tree(c) is None


def test_an_outcome_limitation_can_name_the_judge_s_reason_or_the_answer_s_status(tmp_path):
    from idtrack_reproducibility.families import limits_outcomes as LO
    p = tmp_path / "map.tsv"
    p.write_text("limitation\tgroup\ttitle\twhat\tuser_action\tfixed_by\trule\tclass\tstrata\twhen\n"
                 "L-x\tg\tt\tw\tu\tf\tr\ta\t^old_gene_to_db$\treason:not expected: .*\\b(XM|NG)_\n"
                 "L-y\tg\tt\tw\tu\tf\tr\ta\t^travel_split$\tstatus:^no_target$\n")
    m = LO.load_map(str(p))
    assert LO.classify(m, "old_gene_to_db", {"reason": "not expected: NM_1,XM_2"})["limitation"] == "L-x"
    assert LO.classify(m, "old_gene_to_db", {"reason": "not expected: NM_1"}) is None
    assert LO.classify(m, "travel_split", {"status": "no_target", "reason": "no answer"})["limitation"] == "L-y"
    assert LO.classify(m, "travel_split", {"status": "not_found"}) is None
    bad = tmp_path / "bad.tsv"
    bad.write_text("limitation\tgroup\ttitle\twhat\tuser_action\tfixed_by\trule\tclass\tstrata\twhen\nL\tg\tt\tw\tu\tf\tr\ta\t.\tmatch:x\n")
    with pytest.raises(SystemExit):
        LO.load_map(str(bad))


def test_a_revision_of_sets_beyond_a_string_arrays_limit_is_published(tmp_path, monkeypatch):
    """`publish(revises=...)` checks that every qid of a revision set is published in a set it revises. An organism's
    sets can hold more qid text than one pyarrow `string` array (32-bit offsets, 2 GiB: human C's 61 million qids
    overflowed it when the check joined every revised set's qids). The limit is lowered here to 40 bytes for `string`
    arrays joined from chunks, as pyarrow raises it; the check must hold the revised sets' qids set by set, in
    `large_string`, and still refuse a qid published nowhere. Fails on the joined `string` check."""
    import pyarrow as pa
    real = pa.chunked_array

    class Joined:
        def __init__(self, arrays):
            self.c = real(arrays)

        def combine_chunks(self):
            if self.c.type == pa.string() and sum(len(x.as_py()) for x in self.c) > 40:
                raise pa.ArrowInvalid("offset overflow while concatenating arrays")
            return self.c.combine_chunks()
    monkeypatch.setattr(pa, "chunked_array", lambda arrays, *a, **k: Joined(arrays))
    monkeypatch.setattr(L, "query_dir", lambda family, organism: str(tmp_path / family / organism))

    def line(q):
        return dict(qid=q, family="x", stratum="s", kind="k", id=q, expect=[])
    for s in ("x1_v1", "x2_v1"):
        L.publish("x", "o", s, [line(f"{s}|qid{i}") for i in range(6)], "base", __file__)
    L.publish("x", "o", "x_r1", [line("x1_v1|qid1"), line("x2_v1|qid5")], "revision", __file__, revises=["x1_v1", "x2_v1"])
    assert L.manifest("x", "o", "x_r1")["n"] == 2
    with pytest.raises(SystemExit, match="new"):
        L.publish("x", "o", "x_r2", [line("x1_v1|qid1"), line("nowhere")], "revision", __file__, revises=["x1_v1", "x2_v1"])
