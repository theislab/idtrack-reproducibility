"""Fixtures: a small invented truth store and data root, so that the tests need neither the real data nor idtrack."""
import gzip
import os

import pytest

from idtrack_reproducibility import config as CFG

ORG = "sus_scrofa"


def write_gz(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fo:
        if header:
            fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")


GENE_COLS = ["gene_id", "gene_version", "gene_name", "gene_biotype", "seqname", "start", "end", "strand", "region"]


@pytest.fixture
def data_root(tmp_path, monkeypatch):
    """An empty data root read through the configuration (no configuration file)."""
    monkeypatch.delenv("IDTRACK_REPRO_CONFIG", raising=False)
    for loc in CFG.LOCATIONS:
        monkeypatch.delenv(f"IDTRACK_REPRO_{loc.upper()}", raising=False)
    monkeypatch.setenv("IDTRACK_REPRO_DATA_ROOT", str(tmp_path / "data"))
    monkeypatch.chdir(tmp_path)
    CFG.reset()
    yield tmp_path / "data"
    CFG.reset()


@pytest.fixture
def truth_store(data_root):
    """Three releases of an invented organism history (written under the pig's name):

    release 1: A, C, B0          release 2: A, B, B0          release 3: B, D, F
    session 1 (1 -> 2): A -> A (self), C -> B (C retired into the new B)
    session 2 (2 -> 3): A -> D, A -> B (A split into B and D, no self event); B0 -> F; B lives on, no event of its own;
                        B -> F as well although B is in the GTF of 2 and of 3 (an id that lives on beside an event)
    """
    t = data_root / "truth" / ORG
    genes = {
        1: [("ENSSSCG00000000001", 1, "GA", "protein_coding", "1", 100, 200, "+", "primary"),
            ("ENSSSCG00000000003", 1, "GC", "protein_coding", "1", 1000, 1200, "+", "primary"),
            ("ENSSSCG00000000009", 1, "GB0", "protein_coding", "2", 100, 300, "-", "primary")],
        2: [("ENSSSCG00000000001", 1, "GA", "protein_coding", "1", 100, 200, "+", "primary"),
            ("ENSSSCG00000000002", 1, "GB", "protein_coding", "1", 1000, 1250, "+", "primary"),
            ("ENSSSCG00000000009", 1, "GB0", "protein_coding", "2", 100, 300, "-", "primary")],
        3: [("ENSSSCG00000000002", 2, "GB", "protein_coding", "1", 1000, 1250, "+", "primary"),
            ("ENSSSCG00000000004", 1, "GD", "protein_coding", "1", 100, 180, "+", "primary"),
            ("ENSSSCG00000000006", 1, "GF", "protein_coding", "3", 5, 50, "+", "primary")],
    }
    for r, rows in genes.items():
        write_gz(str(t / "gtf" / f"genes_{r}.tsv.gz"), GENE_COLS, rows)
    core = "sus_scrofa_core_3_1"
    write_gz(str(t / "core" / f"{core}.mapping_session.txt.gz"), None, [
        (1, "sus_scrofa_core_1_1", "sus_scrofa_core_2_1", "1", "2", "Sscrofa1", "Sscrofa1", "2026-01-01 00:00:00"),
        (2, "sus_scrofa_core_2_1", "sus_scrofa_core_3_1", "2", "3", "Sscrofa1", "Sscrofa1", "2026-02-01 00:00:00")])
    write_gz(str(t / "core" / f"{core}.stable_id_event.txt.gz"), None, [
        ("ENSSSCG00000000001", 1, "ENSSSCG00000000001", 1, 1, "gene", 1.0),
        ("ENSSSCG00000000003", 1, "ENSSSCG00000000002", 1, 1, "gene", 0.9),
        ("ENSSSCG00000000001", 1, "ENSSSCG00000000004", 1, 2, "gene", 0.8),
        ("ENSSSCG00000000001", 1, "ENSSSCG00000000002", 2, 2, "gene", 0.3),
        ("ENSSSCG00000000009", 1, "ENSSSCG00000000006", 1, 2, "gene", 0.9),
        ("ENSSSCG00000000002", 1, "ENSSSCG00000000006", 1, 2, "gene", 0.2),
        ("ENSSSCG00000000002", 1, "ENSSSCG00000000002", 2, 2, "transcript", 1.0)])
    with open(t / "SOURCES.tsv", "w") as fo:
        fo.write("what\turl\tbytes\tmd5\tfetched\n")
        for r in (1, 2, 3):
            fo.write(f"gtf/{r}\thttps://example.org/release-{r}/gtf/sus_scrofa/Sus_scrofa.Sscrofa1.{r}.gtf.gz\t1\tx\t"
                     f"2026-11-11T00:00:00\n")
    return t
