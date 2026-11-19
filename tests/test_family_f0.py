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
