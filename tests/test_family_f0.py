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
