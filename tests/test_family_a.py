"""Family A (history travel, exhaustive): its report's axes and the exceptions filter, on invented inputs."""
from idtrack_reproducibility.families.a import gen_a as GA
from idtrack_reproducibility.families.a import report_a as RA


def test_distance_buckets_are_symmetric_and_cover_every_distance():
    assert RA.bucket(0) == "0" and RA.bucket(-1) == "1" and RA.bucket(1) == "1"
    assert RA.bucket(5) == "2-5" and RA.bucket(-40) == "21-40" and RA.bucket(67) == "41+"


class FakeExceptions:
    def __init__(self, rows):
        self._rows, self.cache = rows, {}

    def rows(self, rel):
        return self._rows

    def project(self, place, rel):
        return dict(place, projected=tuple(sorted(self.cache[rel])))


def test_rows_without_coordinates_are_left_out_and_remembered():
    good = dict(start="1", end="10", ref_start="5", ref_end="14")
    bad = dict(start="na", end="na", ref_start="5", ref_end="14")
    ex = GA.ExceptionsA(FakeExceptions({"HAP1": [good, bad], "HAP2": [bad]}))
    assert ex.rows(3) == {"HAP1": [good]}
    assert ex.na[3] == {"HAP1", "HAP2"} and ex.has(3)
    assert ex.project(dict(seqname="HAP1"), 3)["projected"] == ("HAP1",)  # projected with the filtered rows only


def test_the_a8_build_runs_the_tree_s_own_api(tmp_path, monkeypatch):
    from idtrack_reproducibility import config as CFG
    from idtrack_reproducibility.families.a import build_a8
    for name, marker in (("old", "_api.py"), ("new", "api.py")):
        (tmp_path / name / "idtrack").mkdir(parents=True)
        (tmp_path / name / "idtrack" / marker).write_text("")
    cfg = tmp_path / "idtrack-repro.toml"
    cfg.write_text('[data]\nroot = "data"\n[trees.old]\nsource = "old"\npython = "py311"\n'
                   '[trees.new]\nsource = "new"\npython = "py313"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(cfg))
    monkeypatch.delenv("IDTRACK_REPRO_DATA_ROOT", raising=False)
    CFG.reset()
    try:
        old = build_a8.command("old", "sus_scrofa", 110, str(tmp_path / "out"))
        new = build_a8.command("new", "sus_scrofa", 110, str(tmp_path / "out"))
    finally:
        CFG.reset()
    assert old[0] == "py311" and "build_graph(organism_name='sus_scrofa', snapshot_release=110)" in old[2]
    assert new[0] == "py313" and "API.build('sus_scrofa', 110," in new[2]
    assert str(tmp_path / "old") in old[2] and str(tmp_path / "new") in new[2]


def test_a0_summaries_get_the_reasons_of_this_package(tmp_path):
    import json
    from idtrack_reproducibility.families.a import a0_fidelity as A0
    assert A0.design_of(A0.BY_DESIGN["history"], "birth:node_missing:id_in_no_gtf").startswith("builder:")
    assert A0.design_of(A0.BY_DESIGN["names"], "display:other_spelling").startswith("rule `spelling_steps`")
    assert A0.design_of(A0.BY_DESIGN["transcripts"], "tl:graph_extra") == ""
    (tmp_path / "summary.tsv").write_text("check\tclass\tby_design\tcount\tids\texample\n"
                                          "names\tdisplay:other_spelling\tan earlier text\t18\t18\t\n"
                                          "transcripts\ttl:graph_extra\t\t39\t1\t\n")
    (tmp_path / "summary.json").write_text(json.dumps({"organism": "homo_sapiens"}))
    assert A0.resummarize(str(tmp_path)) == 2
    rows = (tmp_path / "summary.tsv").read_text().splitlines()
    assert rows[1].split("\t")[2].startswith("rule `spelling_steps`") and rows[1].split("\t")[3] == "18"
    assert rows[2].split("\t")[2] == "" and "by_design_written" in json.loads((tmp_path / "summary.json").read_text())
