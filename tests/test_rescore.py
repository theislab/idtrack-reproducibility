"""Scoring the stored scorings again under the judge of today (rescore.py), and what it rests on: the sets a scoring
reads today, the lines bounded to the sets it read, and the classifier's order of rows."""
from idtrack_reproducibility import rescore as RS
from idtrack_reproducibility.families import causes as K
from idtrack_reproducibility.families import lines as L


def fake_sets(monkeypatch, sets):
    """sets: [(name, revises, judge)] in publication order."""
    monkeypatch.setattr(L, "set_names", lambda family, organism: [n for n, _, _ in sets])
    man = {n: dict(name=n, revises=r, judge=j) for n, r, j in sets}
    monkeypatch.setattr(L, "manifest", lambda family, organism, name: man[name])


def test_a_scoring_reads_its_sets_and_every_later_revision_of_them(monkeypatch):
    fake_sets(monkeypatch, [("d2_v1", [], None), ("d2_v2", [], None), ("d2_r1", ["d2_v1", "d2_v2"], None),
                            ("d2_r2", ["d2_v1", "d2_v2", "d2_r1"], None), ("d3_v1", [], None)])
    # a run that asked d2_v1 only: the revisions apply, bounded to d2_v1's qids (d2_r1 also revises d2_v2); a scoring
    # with sets added is always bounded to the qids of the sets it read
    assert RS.current_sets("d", "o", ["d2_v1"], False) == (["d2_v1", "d2_r1", "d2_r2"], ["d2_v1"])
    assert RS.current_sets("d", "o", ["d2_v1", "d2_v2"], False) == (["d2_v1", "d2_v2", "d2_r1", "d2_r2"], ["d2_v1", "d2_v2"])
    assert RS.current_sets("d", "o", ["d3_v1"], False) == (["d3_v1"], None)


def test_a_decided_revision_enters_only_a_decided_scoring_after_the_others(monkeypatch):
    fake_sets(monkeypatch, [("e2", [], None), ("e2_dec", ["e2"], "decided"), ("e2_dec2", ["e2", "e2_dec"], "decided"),
                            ("e2_r1", ["e2"], "now"), ("b6_v1", [], None), ("b6_dec", ["b6_v1"], "b_v1")])
    assert RS.current_sets("e", "o", ["e2"], False) == (["e2", "e2_r1"], ["e2"])
    assert RS.current_sets("e", "o", ["e2", "e2_dec"], True) == (["e2", "e2_r1", "e2_dec", "e2_dec2"], ["e2", "e2_dec"])
    assert RS.current_sets("b", "o", ["b6_v1"], False) == (["b6_v1"], None)


def test_lines_bounded_to_the_sets_a_run_asked(data_root):
    line = dict(kind="k", from_release=None, to_release=115, final_database=None, strategy="best", rule="set",
                cmp="versioned", broad=None)
    L.publish("t", "sus_scrofa", "t_v1", [dict(line, qid="a", id="A", expect=["G1.1"])], "v1", __file__)
    L.publish("t", "sus_scrofa", "t_v2", [dict(line, qid="b", id="B", expect=["G2.1"])], "v2", __file__)
    L.publish("t", "sus_scrofa", "t_r1", [dict(line, qid="a", id="A", expect=["G3.1"]),
                                          dict(line, qid="b", id="B", expect=["G4.1"])], "r1", __file__,
              revises=["t_v1", "t_v2"])
    got, _, _ = L.load("t", "sus_scrofa", ["t_v1", "t_r1"], within=["t_v1"])
    assert [(q["qid"], q["expect"]) for q in got] == [("a", ["G3.1"])]
    got, _, _ = L.load("t", "sus_scrofa", ["t_v1", "t_r1"])
    assert sorted(q["qid"] for q in got) == ["a", "b"]


def test_the_classifier_meets_rows_release_by_release_and_answers_in_their_order():
    rows = [dict(query='{"from_release": 90, "to_release": 115}'), dict(query='{"from_release": null, "to_release": 115}'),
            dict(query='{"from_release": 60, "to_release": 115}'), dict(query='{"from_release": 90, "to_release": 100}')]
    assert K.release_order(rows) == [1, 2, 3, 0]
