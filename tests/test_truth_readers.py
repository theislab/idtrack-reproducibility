"""The truth readers and the history walk on an invented history (conftest.truth_store)."""
import math

from idtrack_reproducibility.harness import harness_lib as H
from idtrack_reproducibility.harness import make_queries as M

from conftest import ORG

A, B, C, D, F, B0 = (f"ENSSSCG0000000000{i}" for i in (1, 2, 3, 4, 6, 9))


def test_gtf_readers(truth_store):
    assert H.gtf_releases(ORG) == [1, 2, 3]
    g = H.gtf_genes(ORG, 3)
    assert set(g) == {B, D, F}
    assert g[B]["gene_version"] == "2"
    assert H.gtf_genes(ORG, 4) is None
    assert H.gtf_assemblies(ORG) == {1: "Sscrofa1", 2: "Sscrofa1", 3: "Sscrofa1"}


def test_core_tables_keep_gene_events_only(truth_store):
    core, ev, sessions = H.core_tables(ORG, 3)
    assert core == "sus_scrofa_core_3_1"
    assert len(ev) == 6  # the transcript event is left out
    assert sessions[2]["old_release"] == "2" and sessions[2]["new_release"] == "3"


def test_set_names_in_publication_order(tmp_path):
    for n in ("v1", "v2", "v10", "v9"):
        (tmp_path / f"{n}.jsonl.gz").write_bytes(b"")
    assert H.set_names(str(tmp_path)) == ["v1", "v2", "v9", "v10"]  # a plain sort puts v10 before v2


def test_strict_walk_follows_a_split(truth_store):
    t = M.Truth(ORG, 3)
    assert t.walk_fwd(A, 1, 3, strict=True) == {B, D}
    assert t.walk_back(D, 3, 1, strict=True) == {A}


def test_self_event_keeps_only_the_id_in_the_strict_walk(truth_store):
    t = M.Truth(ORG, 3)
    assert t.walk_fwd(A, 1, 2, strict=True) == {A}


def test_an_id_that_lives_on(truth_store):
    """B is in the GTF on both sides of session 2, which records only B -> F: the strict walk leaves B unless
    `living` (rule `id_lives_on`) keeps it."""
    t = M.Truth(ORG, 3)
    assert t.walk_fwd(B, 2, 3, strict=True) == {F}
    assert t.walk_fwd(B, 2, 3, strict=True, living=True) == {B}
    assert t.walk_fwd(C, 1, 3, strict=True, living=True) == {B}
    assert t.walk_back(B, 3, 1, strict=True) == {A}
    assert t.walk_back(B, 3, 1, strict=True, living=True) == {C}


def test_broad_walk_adds_the_branches(truth_store):
    t = M.Truth(ORG, 3)
    assert t.walk_fwd(B, 2, 3, strict=False, living=True) == {B, F}


def test_walk_series_equals_single_walks(truth_store):
    t = M.Truth(ORG, 3)
    series = t.walk_series(A, 1, [1, 2, 3], strict=True)
    assert series == {1: {A}, 2: {A}, 3: t.walk_fwd(A, 1, 3, strict=True)}
    back = t.walk_series(D, 3, [1, 2], strict=True)
    assert back[1] == t.walk_back(D, 3, 1, strict=True)


def test_versions_of_the_release(truth_store):
    t = M.Truth(ORG, 3)
    assert t.vid(B, 3) == f"{B}.2" and t.vid(B, 2) == f"{B}.1"
    assert t.cmp_for(3) == "versioned"


def test_wilson_interval():
    lo, hi = H.wilson(50, 100)
    assert round(lo, 4) == 0.4038 and round(hi, 4) == 0.5962
    assert all(math.isnan(x) for x in H.wilson(0, 0))  # no interval without lines
    lo, hi = H.wilson(10, 10)
    assert round(hi, 9) == 1.0 and 0.72 < lo < 0.73
