"""A tree's graphs: the table farm of a build, graph file names, the difference of two graphs, the history check, the
stored copy and the tree's own verify, on small invented cases (no data, no idtrack)."""
import json
import os
import sys
import zipfile

import numpy as np
import pandas as pd
import pytest

from idtrack_reproducibility.graphs import build as GB
from idtrack_reproducibility.graphs import common as G
from idtrack_reproducibility.graphs import diff as D
from idtrack_reproducibility.graphs import judge as J
from idtrack_reproducibility.graphs import stored as S


def touch(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def test_graph_file_names_of_both_generations_are_read_and_others_refused(tmp_path):
    assert G.graph_of("graph_sus_scrofa_assembly111_min56_max115.npz") == {
        "organism": "sus_scrofa", "assembly": 111, "first": 56, "snapshot": 115}
    assert G.graph_of("/x/graph_homo_sapiens_assembly38_min48_max115_narrow.npz")["snapshot"] == 115
    with pytest.raises(SystemExit):
        G.graph_of("graph_sus_scrofa_115.npz")
    assert G.graph_in(str(tmp_path), "sus_scrofa") is None
    touch(tmp_path / "graph_sus_scrofa_assembly111_min56_max115.npz")
    touch(tmp_path / "graph_mus_musculus_assembly39_min48_max115.npz")
    assert G.graph_in(str(tmp_path), "sus_scrofa").endswith("graph_sus_scrofa_assembly111_min56_max115.npz")
    touch(tmp_path / "graph_sus_scrofa_assembly111_min56_max110.npz")
    with pytest.raises(SystemExit):
        G.graph_in(str(tmp_path), "sus_scrofa")


def test_the_farm_links_ensembls_tables_copies_the_indices_and_leaves_the_derived_tables_out(tmp_path):
    src, repo = tmp_path / "stored", tmp_path / "build"
    tables = src / "idtrack_tables" / "sus_scrofa_assembly-111"
    touch(tables / "ens100_mysql_gene__cols=gene_id+stable_id.h5", "raw")
    touch(tables / "ens100_external-database_gene.h5", "derived")
    touch(tables / "ens100_relationcurrent.h5", "derived")
    touch(src / "assembly_chain" / "sus_scrofa" / "Sscrofa10.2_to_Sscrofa11.1.chain.gz", "chain")
    touch(src / "ensembl_https_index_sus_scrofa_assembly-111.json", "{}")
    touch(src / "graph_sus_scrofa_assembly111_min56_max115.npz", "a graph")
    # the stored tables may themselves be links (a farm of a farm): the links point to the resolved files
    os.symlink(tables / "ens100_mysql_gene__cols=gene_id+stable_id.h5", tables / "ens101_mysql_gene__cols=gene_id+stable_id.h5")
    farmed = G.farm(str(src), str(repo))
    assert farmed["counts"] == {"linked": 3, "copied": 1, "left_out": 3, "kept": 0}
    raw = repo / "idtrack_tables" / "sus_scrofa_assembly-111" / "ens101_mysql_gene__cols=gene_id+stable_id.h5"
    assert os.path.islink(raw) and os.readlink(raw) == os.path.realpath(tables / "ens100_mysql_gene__cols=gene_id+stable_id.h5")
    assert os.path.islink(repo / "assembly_chain" / "sus_scrofa" / "Sscrofa10.2_to_Sscrofa11.1.chain.gz")
    index = repo / "ensembl_https_index_sus_scrofa_assembly-111.json"
    assert index.is_file() and not os.path.islink(index)
    assert not (repo / "idtrack_tables" / "sus_scrofa_assembly-111" / "ens100_external-database_gene.h5").exists()
    assert not list(repo.glob("graph_*"))
    assert G.farm(str(src), str(repo))["counts"]["kept"] == 4  # a second farm keeps what is there
    assert G.changed_sources(farmed["sources"]) == []
    touch(tables / "ens100_mysql_gene__cols=gene_id+stable_id.h5", "written through a link")
    assert G.changed_sources(farmed["sources"]) == [os.path.realpath(tables / "ens100_mysql_gene__cols=gene_id+stable_id.h5")]


def test_a_linked_folder_is_refused(tmp_path):
    touch(tmp_path / "elsewhere" / "ens100_mysql_gene__cols=gene_id.h5")
    os.makedirs(tmp_path / "stored" / "idtrack_tables")
    os.symlink(tmp_path / "elsewhere", tmp_path / "stored" / "idtrack_tables" / "sus_scrofa_assembly-111")
    with pytest.raises(SystemExit):
        G.farm(str(tmp_path / "stored"), str(tmp_path / "build"))


def test_the_content_digest_is_read_from_the_record_else_from_the_manifest():
    class Record:
        content_digest = "abc"

    class Typed:
        record = Record()

    class Old:
        record = None
        graph = {"manifest": {"content_digest": "def"}}

    assert GB.content_digest(Typed()) == "abc" and GB.content_digest(Old()) == "def"


def arrays(nodes, links, history):
    """The arrays of a graph file: nodes [(name, type)], links [(external, gene, [(database, mask)])],
    history [(source, target, old release, new release)]."""
    node_types = sorted({t for _, t in nodes})
    names = [n for n, _ in nodes]
    databases = sorted({db for _, _, conns in links for db, _ in conns})
    e_kind, e_src, e_dst, e_ref, c_ptr, c_db, c_asm, c_mask, h_old, h_new = [], [], [], [], [0], [], [], [], [], []
    for ext, gene, conns in links:
        e_kind.append(1), e_src.append(names.index(ext)), e_dst.append(names.index(gene)), e_ref.append(len(c_ptr) - 1)
        for db, mask in conns:
            c_db.append(databases.index(db)), c_asm.append(111), c_mask.append([mask, 0])
        c_ptr.append(len(c_db))
    for a, b, old, new in history:
        e_kind.append(0), e_src.append(names.index(a)), e_dst.append(names.index(b)), e_ref.append(len(h_old))
        h_old.append(old), h_new.append(new)
    document = {"node_types": node_types, "databases": databases, "record": {}}
    return {"document": np.frombuffer(json.dumps(document).encode(), dtype=np.uint8),
            "names": np.frombuffer("\x00".join(names).encode(), dtype=np.uint8),
            "ntype": np.array([node_types.index(t) for _, t in nodes]), "e_kind": np.array(e_kind),
            "e_src": np.array(e_src), "e_dst": np.array(e_dst), "e_ref": np.array(e_ref), "c_ptr": np.array(c_ptr),
            "c_db": np.array(c_db, dtype=np.int64), "c_asm": np.array(c_asm, dtype=np.int64),
            "c_mask": np.array(c_mask, dtype=np.uint64).reshape(-1, 2), "h_old": np.array(h_old), "h_new": np.array(h_new)}


def test_the_difference_of_two_graphs_counts_nodes_links_by_database_and_history(tmp_path):
    genes = [("G1.1", "gene"), ("G2.1", "gene"), ("G2.2", "gene")]
    old = arrays(genes + [("X1", "external"), ("X2", "external"), ("---", "external")],
                 [("X1", "G1.1", [("RefSeq_dna", 0b110)]), ("X2", "G2.1", [("HGNC", 0b10)]), ("---", "G1.1", [("HGNC", 0b1)])],
                 [("G2.1", "G2.2", 100, 101)])
    new = arrays(genes + [("X1", "external"), ("X2", "external"), ("X3 {ECO:1}", "external")],
                 [("X1", "G1.1", [("RefSeq_mRNA", 0b110)]), ("X2", "G2.1", [("HGNC", 0b110)]),
                  ("X3 {ECO:1}", "G2.2", [("UniProt", 0b1)])],
                 [("G2.1", "G2.2", 100, 101), ("G1.1", "G2.2", 101, 102)])
    np.savez(tmp_path / "old.npz", **old)
    result = D.diff(D.read(str(tmp_path / "old.npz")), D.from_arrays(new))
    assert result["nodes"]["added"] == {"external: UniProt evidence text": 1}
    assert result["nodes"]["removed"] == {"external: no letter or digit": 1}
    assert result["links_moved_between_databases"] == {"RefSeq_dna -> RefSeq_mRNA": 1}
    assert result["links"] == {"HGNC": {"other_releases": 1, "only_old": 1}, "UniProt": {"only_new": 1}}
    assert (result["history"]["n_added"], result["history"]["n_removed"]) == (1, 0)
    assert result["history"]["added"] == [["G1.1", "G2.2", 101, 102]]


def events(rows):
    return pd.DataFrame(rows, columns=["old_stable_id", "old_version", "new_stable_id", "new_version", "old_release",
                                       "new_release"])


def test_the_history_check_tells_held_placed_absent_missing_and_invented_edges():
    truth = events([("A", "1", "B", "1", "100", "101"),    # held
                    ("A", "1", "C", "2", "100", "101"),    # placed with the version of the release: C.3
                    ("A", "1", "Z", "1", "100", "101"),    # Z is in no table of the graph
                    ("B", "1", "C", "3", "101", "102")])   # both in the graph, no edge: missing
    pairs = {("A.1", "B.1"), ("A.1", "C.3"), ("C.3", "D.1")}
    placed = {("A.1", "C.2"): ("A.1", "C.3")}
    result = J.compare_history(truth, pairs, placed, {"A.1", "B.1", "C.3", "D.1"})
    assert (result["held"], result["placed_with_the_version_of_the_release"], result["ids_in_no_table"]) == (1, 1, 1)
    assert result["missing"] == 1 and result["missing_examples"] == [["B.1", "C.3", "101", "102"]]
    assert result["edges_without_event"] == 1 and result["edges_without_event_examples"] == [["C.3", "D.1"]]


def test_gene_events_keep_gene_events_between_different_ids_of_the_graphs_window():
    ev = pd.DataFrame([["A", "1", "B", "1", "1", "gene", "1"], ["A", "1", "A", "2", "1", "gene", "1"],
                       ["A", "1", "B", "0", "1", "gene", "1"], ["T", "1", "U", "1", "1", "transcript", "1"],
                       ["C", "1", "D", "1", "2", "gene", "1"], ["E", "1", None, None, "1", "gene", "1"]])
    sessions = pd.DataFrame([["1", "a", "b", "100", "101", "x", "x", "t"], ["2", "a", "b", "90", "91", "x", "x", "t"]])
    kept = J.gene_events(ev, sessions, 95, 115)
    assert kept[["old_stable_id", "new_stable_id"]].values.tolist() == [["A", "B"]]


def test_the_stored_copy_holds_the_same_members_uncompressed_and_loads_the_same(tmp_path):
    np.savez_compressed(tmp_path / "graph.npz", a=np.arange(1000), b=np.frombuffer(b"names\x00x", dtype=np.uint8))
    copy = S.stored_copy(str(tmp_path / "graph.npz"), str(tmp_path / "stored"))
    with zipfile.ZipFile(copy) as z:
        assert {i.compress_type for i in z.infolist()} == {zipfile.ZIP_STORED}
    assert S.same_members(str(tmp_path / "graph.npz"), copy)
    with np.load(copy) as x, np.load(tmp_path / "graph.npz") as y:
        assert all(np.array_equal(x[k], y[k]) for k in y.files)
    assert S.sha256(copy) == S.sha256(S.stored_copy(str(tmp_path / "graph.npz"), str(tmp_path / "again")))
    with pytest.raises(SystemExit):
        S.stored_copy(copy, str(tmp_path / "stored"))


def test_verify_is_the_trees_own_command_or_not_available(tmp_path):
    old = tmp_path / "old"
    touch(old / "idtrack" / "__init__.py", "")
    assert J.verify_here(sys.executable, str(old), "g.npz", None)["status"].startswith("not available")
    tree = tmp_path / "tree"
    touch(tree / "idtrack" / "__init__.py", "")
    touch(tree / "idtrack" / "__main__.py", 'import sys\nchoices = ["verify"]\nprint(" ".join(sys.argv[1:]))\n'
          'sys.exit(0 if sys.argv[2] == "good.npz" else 1)\n')
    good = J.verify_here(sys.executable, str(tree), "good.npz", str(tmp_path))
    assert good["status"] == "holds" and good["output"] == [f"verify good.npz --repository {tmp_path}"]
    assert J.verify_here(sys.executable, str(tree), "bad.npz", None)["status"] == "fails"


FAKE_IDTRACK = '''import os


class _Record:
    content_digest = "digest-of-the-new-build"


class _Graph:
    record = _Record()


class API:
    graph = _Graph()

    @classmethod
    def build(cls, organism, snapshot, repository):
        name = f"graph_{organism}_assembly111_min56_max{snapshot}.npz"
        with open(os.path.join(repository, name), "wb") as fh:
            fh.write(b"the new build")
        also = os.environ.get("FAKE_IDTRACK_ALSO_WRITE")
        if also:  # another process puts a graph file into the tree's folder meanwhile
            with open(os.path.join(also, name), "wb") as fh:
                fh.write(b"someone else's graph")
        return cls()
'''


@pytest.fixture
def fake_tree(data_root, tmp_path, monkeypatch):
    """A configured tree `t` whose idtrack is a stand-in that builds a small file (no idtrack, no tables)."""
    from idtrack_reproducibility import config as CFG
    src = tmp_path / "tree"
    touch(src / "idtrack" / "__init__.py", FAKE_IDTRACK)
    touch(src / "idtrack" / "api.py", "")
    os.makedirs(tmp_path / "tables")
    cfg = tmp_path / "idtrack-repro.toml"
    cfg.write_text(f'[trees.t]\ncommit = "abc"\nsource = "{src}"\npython = "{sys.executable}"\n')
    monkeypatch.setenv("IDTRACK_REPRO_CONFIG", str(cfg))
    monkeypatch.setattr(sys, "path", list(sys.path))
    monkeypatch.delitem(sys.modules, "idtrack", raising=False)
    CFG.reset()
    yield data_root / "graphs" / "t" / "sus_scrofa", str(tmp_path / "tables")
    sys.modules.pop("idtrack", None)


def test_a_build_records_the_file_it_wrote_and_a_resubmission_records_it_again(fake_tree):
    out, tables = fake_tree
    r = GB.build_here("t", "sus_scrofa", 115, tables)
    graph = out / "graph_sus_scrofa_assembly111_min56_max115.npz"
    assert graph.read_bytes() == b"the new build" and r["content_digest"] == "digest-of-the-new-build"
    assert json.load(open(out / "build.json"))["graph_file_sha256"] == GB.sha256(str(graph))
    sys.modules.pop("idtrack", None)
    assert GB.build_here("t", "sus_scrofa", 115, tables)["graph_file_sha256"] == GB.sha256(str(graph))


def test_an_existing_graph_that_is_not_the_builds_is_never_described_by_its_record(fake_tree):
    out, tables = fake_tree
    graph = out / "graph_sus_scrofa_assembly111_min56_max115.npz"
    touch(graph, "a graph made elsewhere")
    with pytest.raises(SystemExit, match="not this build's graph.*--keep-existing"):
        GB.build_here("t", "sus_scrofa", 115, tables)
    assert graph.read_text() == "a graph made elsewhere" and not (out / "build.json").exists()
    assert GB.build_here("t", "sus_scrofa", 115, tables, keep_existing=True) is None
    assert graph.read_text() == "a graph made elsewhere" and not (out / "build.json").exists()


def test_a_graph_file_that_appears_during_the_build_is_not_recorded_either(fake_tree, monkeypatch):
    out, tables = fake_tree
    os.makedirs(out)
    monkeypatch.setenv("FAKE_IDTRACK_ALSO_WRITE", str(out))
    with pytest.raises(SystemExit, match="not this build's graph"):
        GB.build_here("t", "sus_scrofa", 115, tables)
    assert (out / "graph_sus_scrofa_assembly111_min56_max115.npz").read_bytes() == b"someone else's graph"
    assert not (out / "build.json").exists()
