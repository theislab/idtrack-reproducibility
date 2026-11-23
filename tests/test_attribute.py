"""The attribution of a decision's effect: moved answers between consecutive harness scorings, per step, on invented
runs (no data, no idtrack)."""
import gzip
import os

from idtrack_reproducibility.harness import attribute as A


def scored(folder, outcomes):
    os.makedirs(folder, exist_ok=True)
    with gzip.open(os.path.join(folder, "outcomes.tsv.gz"), "wt") as fo:
        fo.write("qid\tkind\tvariant\toutcome\tgot\treason\n")
        for qid, outcome in outcomes.items():
            fo.write(f"{qid}\tk\t*\t{outcome}\t-\t-\n")


def test_each_moved_answer_is_counted_at_the_step_that_moved_it(tmp_path, capsys):
    root = str(tmp_path)
    scored(f"{root}/base/sus_scrofa", {"q1": "wrong", "q2": "undetermined", "q3": "right", "q4": "right"})
    scored(f"{root}/part1/sus_scrofa", {"q1": "right", "q2": "undetermined", "q3": "right", "q5": "wrong"})
    scored(f"{root}/part2/sus_scrofa", {"q1": "right", "q2": "right", "q3": "wrong", "q5": "right"})
    scored(f"{root}/base/mus_musculus", {"q1": "right"})
    scored(f"{root}/part1/mus_musculus", {"q1": "undetermined"})  # no part2 of mouse: skipped with a note
    runs = [(n, f"{root}/{n}/{{organism}}") for n in ("base", "part1", "part2")]
    rows = A.attribute(runs, ["sus_scrofa", "mus_musculus"])
    assert rows == [("base->part1", "sus_scrofa", "wrong -> right", 1),  # q4, q5: scored in one run only
                    ("base->part1", "mus_musculus", "right -> undetermined", 1),
                    ("part1->part2", "sus_scrofa", "right -> wrong", 1),
                    ("part1->part2", "sus_scrofa", "undetermined -> right", 1),
                    ("part1->part2", "sus_scrofa", "wrong -> right", 1)]
    assert "mus_musculus skipped in part1->part2" in capsys.readouterr().err


def test_the_table_has_one_row_per_step_organism_and_move(tmp_path):
    scored(f"{tmp_path}/a/sus_scrofa", {"q1": "wrong"})
    scored(f"{tmp_path}/b/sus_scrofa", {"q1": "right"})
    out = tmp_path / "attribution.tsv"
    assert A.main(["--run", f"a={tmp_path}/a/{{organism}}", "--run", f"b={tmp_path}/b/{{organism}}",
                   "--organisms", "sus_scrofa", "--out", str(out)]) == 0
    assert out.read_text() == "step\torganism\tmove\tanswers\na->b\tsus_scrofa\twrong -> right\t1\n"
