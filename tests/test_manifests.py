"""Data manifests: one row per file with its source and rebuild, and the check against the data root."""
import os

from idtrack_reproducibility import manifests as MAN


def setup(tmp_path, data_root):
    repo = tmp_path / "repo"
    (repo / "manifests").mkdir(parents=True)
    (repo / "manifests" / "areas.toml").write_text(
        '[[area]]\nname = "truth"\nfolders = ["truth"]\nexclude = [\'\\.tmp$\']\nrebuild = "fetch"\n'
        "rules = [['/rest/', \"none: a service's answer of a date\"]]\n")
    t = data_root / "truth" / "sus_scrofa"
    (t / "gtf").mkdir(parents=True)
    (t / "rest").mkdir()
    (t / "gtf" / "genes_3.tsv.gz").write_bytes(b"abc")
    (t / "rest" / "archive.json").write_text("{}")
    (t / "gtf" / "x.tmp").write_text("scratch")
    (t / "SOURCES.tsv").write_text("what\turl\tbytes\tmd5\tfetched\ngtf/3\thttps://example.org/g.gtf.gz\t3\tm\t2026-11-11T00:00:00\n")
    return str(repo), t


def test_build_lists_every_file_with_source_and_rebuild(tmp_path, data_root):
    repo, t = setup(tmp_path, data_root)
    area = MAN.areas(repo)[0]
    rows = {r["path"]: r for r in MAN.build(repo, "truth", area)}
    assert set(rows) == {"truth/sus_scrofa/gtf/genes_3.tsv.gz", "truth/sus_scrofa/rest/archive.json",
                         "truth/sus_scrofa/SOURCES.tsv"}
    g = rows["truth/sus_scrofa/gtf/genes_3.tsv.gz"]
    assert g["source"] == "https://example.org/g.gtf.gz" and g["time"] == "2026-11-11T00:00:00" and g["bytes"] == 3
    assert rows["truth/sus_scrofa/rest/archive.json"]["rebuild"].startswith("none:")
    assert os.path.exists(os.path.join(repo, "manifests", "truth.tsv"))


def test_check_finds_changed_and_missing_files(tmp_path, data_root):
    repo, t = setup(tmp_path, data_root)
    MAN.build(repo, "truth", MAN.areas(repo)[0])
    problems, skipped = [], []
    MAN.check(repo, problems, skipped)
    assert problems == []
    (t / "gtf" / "genes_3.tsv.gz").write_bytes(b"abd")  # same size, other content
    os.remove(t / "rest" / "archive.json")
    MAN.check(repo, problems, skipped)
    assert sorted(p[2].split(":")[0] for p in problems) == ["missing in the data root", "sha256 differs"]
