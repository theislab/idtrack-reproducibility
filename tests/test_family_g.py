"""Family G (other tools): the client's rate, backoff and resume on a fake server, the adapters on small invented
responses, the rule of equal truth on an invented history, the sampling's allocation and which lines a tool can pose.
No network, no data root beyond a temporary one."""
import collections
import json

import pytest

from idtrack_reproducibility import records
from idtrack_reproducibility.families.g import adapters as A
from idtrack_reproducibility.families.g import client as K
from idtrack_reproducibility.families.g import likeforlike as LL
from idtrack_reproducibility.families.g import sets as GS
from idtrack_reproducibility.families.g import tools as T


class FakeResponse:
    def __init__(self, status, content, headers=None):
        self.status_code, self.content, self.headers = status, content, headers or {}


class FakeSession:
    def __init__(self, script):
        self.script, self.calls = list(script), []

    def request(self, method, url, **kw):
        self.calls.append((method, url, kw.get("params"), kw.get("json"), kw.get("data")))
        r = self.script.pop(0)
        if isinstance(r, Exception):
            raise r
        return r


class Clock:
    def __init__(self):
        self.t, self.slept = 0.0, []

    def __call__(self):
        return self.t

    def sleep(self, s):
        self.slept.append(round(s, 3))
        self.t += s


def client_with(script):
    sess, clock = FakeSession(script), Clock()
    return K.Client(session_factory=lambda: sess, sleep=clock.sleep, clock=clock), sess, clock


def test_one_host_is_asked_at_its_interval(tmp_path):
    c, sess, clock = client_with([FakeResponse(200, b"[]"), FakeResponse(200, b"[1]")])
    st = K.RawStore(str(tmp_path), "ensembl_rest", "116", "sus_scrofa")
    c.fetch(st, "GET", "https://rest.ensembl.org/a", qids=["q1"], endpoint="a")
    c.fetch(st, "GET", "https://rest.ensembl.org/b", qids=["q2"], endpoint="b")
    # two requests to one host: the second waits the host's interval (7.5 per second: 0.133 s)
    assert clock.slept == [round(K.HOSTS["ensembl_rest"].interval, 3)]


def test_a_refused_request_is_waited_out_for_its_retry_after_then_asked_again(tmp_path):
    c, sess, clock = client_with([FakeResponse(429, b"slow down", {"Retry-After": "12.5"}), FakeResponse(200, b"{}")])
    st = K.RawStore(str(tmp_path), "mygene", "x", "sus_scrofa")
    row = c.fetch(st, "POST", "https://mygene.info/v3/query", data={"q": "A"}, qids=["q1"], endpoint="query")
    assert 12.5 in clock.slept and row["http_status"] == 200 and not row["error"] and len(sess.calls) == 2


def test_a_failing_server_is_retried_then_recorded_with_its_error(tmp_path):
    c, sess, clock = client_with([FakeResponse(500, b"")] * (len(K.RETRY_WAITS) + 1))
    st = K.RawStore(str(tmp_path), "biomart", "115", "sus_scrofa")
    row = c.fetch(st, "GET", "https://sep2025.archive.ensembl.org/biomart/martservice", params={"type": "registry"})
    assert row["error"] == "HTTP 500" and [s for s in clock.slept if s in K.RETRY_WAITS] == list(K.RETRY_WAITS)


def test_a_stored_answer_is_not_asked_again_and_the_index_is_the_records(tmp_path):
    c, sess, clock = client_with([FakeResponse(200, b'[{"id":"X"}]')])
    st = K.RawStore(str(tmp_path), "ensembl_rest", "116", "sus_scrofa")
    a = c.fetch(st, "POST", "https://rest.ensembl.org/archive/id", json_body={"id": ["X"]}, qids=["q1"], endpoint="archive/id")
    st2 = K.RawStore(str(tmp_path), "ensembl_rest", "116", "sus_scrofa")  # a new pass reads the store
    b = c.fetch(st2, "POST", "https://rest.ensembl.org/archive/id", json_body={"id": ["X"]}, qids=["q1"], endpoint="archive/id")
    assert a["request_id"] == b["request_id"] and len(sess.calls) == 1
    assert st2.answered("archive/id") == {"q1"} and st2.body(b) == b'[{"id":"X"}]'
    st2.compile()
    raw = records.load_raw(str(tmp_path))
    assert list(raw.columns) == records.RAW_COLUMNS and raw.iloc[0]["qids"].tolist() == ["q1"]


def test_a_version_call_is_asked_again_by_every_pass():
    assert K.request_id("GET", "u", salt="p1") != K.request_id("GET", "u", salt="p2") != K.request_id("GET", "u")


def test_hosts_of_one_service_share_its_rate():
    assert K.host_group("https://sep2025.archive.ensembl.org/biomart") == K.host_group("https://dec2021.archive.ensembl.org/x")
    with pytest.raises(ValueError):
        K.host_group("https://example.org/")


def test_retry_after_reads_seconds_and_dates():
    assert K.retry_after("3", 60) == 3.0 and K.retry_after(None, 60) == 60 and K.retry_after("garbage", 7) == 7


def store_with(tmp_path, tool, rows):
    st = K.RawStore(str(tmp_path / tool), tool, "v", "sus_scrofa")
    for i, (endpoint, qids, body, status) in enumerate(rows):
        st.put(dict(request_id=f"r{i}", tool=tool, tool_version="v", organism="sus_scrofa", endpoint=endpoint,
                    method="GET", url="u", params={}, body=None, qids=qids, asked_at=f"2026-11-24T00:00:0{i}",
                    http_status=status, seconds=1.0, bytes=None, sha256=None, file=None, error=None),
               body if isinstance(body, bytes) else json.dumps(body).encode(), "json", {})
    return st


def line(qid, ident, task, **kw):
    return dict(dict(qid=qid, id=ident, task=task, to_release=115, from_release=None, source_family="b",
                     source_stratum="b1_x", expect=[], cmp="base"), **kw)


def test_rest_archive_answers_a_current_id_itself_and_a_retired_one_its_replacements(tmp_path):
    st = store_with(tmp_path, "ensembl_rest", [("archive/id", ["q1", "q2", "q3"], [
        {"id": "ENSSSCG00000000001", "latest": "ENSSSCG00000000001.4", "is_current": "1", "possible_replacement": []},
        {"id": "ENSSSCG00000000002", "latest": "ENSSSCG00000000002.1", "is_current": "",
         "possible_replacement": [{"stable_id": "ENSSSCG00000000005", "score": 1}]}], 200)])
    ls = [line("q1", "ENSSSCG00000000001", "a"), line("q2", "ENSSSCG00000000002.1", "a"), line("q3", "ENSSSCG00000000003", "a")]
    out = A.rest(st, ls, "archive/id")
    assert out["q1"].got == ["ENSSSCG00000000001"] and out["q1"].status == "current"
    assert out["q2"].got == ["ENSSSCG00000000005"] and out["q3"].got == [] and out["q3"].status == "not_found"
    assert out["q1"].seconds == pytest.approx(1 / 3)


def test_gprofiler_none_and_na_are_no_conversion(tmp_path):
    st = store_with(tmp_path, "gprofiler", [("convert", ["q1", "q2", "q3"], {"result": [
        {"incoming": "TP53", "converted": "ENSSSCG00000000009"}, {"incoming": "XX", "converted": "None"},
        {"incoming": "YY", "converted": "N/A"}]}, 200)])
    out = A.gprofiler(st, [line("q1", "TP53", "c"), line("q2", "XX", "c"), line("q3", "YY", "c")])
    assert out["q1"].got == ["ENSSSCG00000000009"] and out["q2"].got == [] and out["q3"].got == []


def test_biomart_reads_rows_by_the_filter_value_and_needs_its_success_line(tmp_path):
    good = b"Gene stable ID\tGene name\nENSSSCG00000000001\tTP53\nENSSSCG00000000002\tTp53\n[success]\n"
    ids = b"Gene stable ID\tGene stable ID version\nENSSSCG00000000007\tENSSSCG00000000007.3\n[success]\n"
    st = store_with(tmp_path, "biomart", [("martservice:external_gene_name", ["q1"], good, 200),
                                          ("martservice:external_gene_name", ["q2"], b"Gene stable ID\tGene name\nE\tX\n", 200),
                                          ("martservice:ensembl_gene_id", ["q3"], ids, 200)])
    out = A.biomart(st, [line("q1", "TP53", "c"), line("q2", "X", "c"), line("q3", "ENSSSCG00000000007", "b")], "sus_scrofa")
    assert out["q1"].got == ["ENSSSCG00000000001", "ENSSSCG00000000002"]  # the mart's filter matches without case
    assert out["q2"].outcome == "error"  # no [success]: the response may be cut
    assert out["q3"].got == ["ENSSSCG00000000007"]


def test_mygene_leaves_out_other_annotations_of_the_species(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "reference_genes", lambda o, r: {"ENSSSCG00000000001"})
    st = store_with(tmp_path, "mygene", [("query:symbol", ["q1", "q2"], [
        {"query": "A", "ensembl.gene": ["ENSSSCG00000000001", "ENSSSCG00015000001"]},
        {"query": "B", "ensembl.gene": "ENSSSCG00015000002"}], 200)])
    out = A.mygene(st, [line("q1", "A", "c"), line("q2", "B", "c")], "sus_scrofa", 116)
    assert out["q1"].got == ["ENSSSCG00000000001"] and out["q2"].got == [] and out["q2"].status == "other_annotation_only"


def test_hgnc_takes_the_first_field_of_the_chain_that_finds_a_record(tmp_path):
    st = store_with(tmp_path, "hgnc", [("fetch/symbol", ["q1"], {"response": {"numFound": 0, "docs": []}}, 200),
                                       ("fetch/prev_symbol", ["q1"], {"response": {"numFound": 1, "docs": [
                                           {"symbol": "NEW", "ensembl_gene_id": "ENSG00000000001"}]}}, 200)])
    out = A.hgnc(st, [line("q1", "OLD", "d")])
    assert out["q1"].got == ["ENSG00000000001"] and out["q1"].status == "prev_symbol"


def test_bridgedb_reads_the_ensembl_mappings(tmp_path, monkeypatch):
    st = store_with(tmp_path, "bridgedb", [("xrefsBatch/L", ["q1", "q2"],
                                            b"100\tEntrez Gene\tQ:NM_1,En:ENSSSCG00000000003,H:X\n200\tEntrez Gene\tN/A\n", 200)])
    ls = [line("q1", "100", "e", source_family="c", source_stratum="c1_default"),
          line("q2", "200", "e", source_family="c", source_stratum="c1_default")]
    out = A.bridgedb(st, ls, "sus_scrofa")
    assert out["q1"].got == ["ENSSSCG00000000003"] and out["q2"].got == []


def test_a_line_is_posed_at_the_tools_release_or_not_supported():
    rest, mart = T.TOOLS["ensembl_rest"], T.TOOLS["biomart"]
    to115 = line("q", "TP53", "c")
    atR = line("q", "TP53", "c", from_release=90, to_release=90)
    assert rest.posable(to115, "sus_scrofa") == (116, "")
    assert rest.posable(atR, "sus_scrofa")[0] is None  # a current-only service cannot answer at 90
    assert mart.posable(to115, "sus_scrofa") == (115, "")  # BioMart's archive of 115
    assert mart.posable(dict(atR, from_release=106, to_release=106), "sus_scrofa") == (106, "")
    assert T.TOOLS["hgnc"].posable(to115, "sus_scrofa")[0] is None  # human only
    assert T.TOOLS["gprofiler"].posable(line("q", "ENS", "a"), "sus_scrofa")[0] is None  # history is not its task
    assert mart.posable(line("q", "A,B", "c"), "sus_scrofa")[0] is None  # a comma cannot be a filter value


def test_gprofiler_archives_are_read_from_its_listing():
    listing = [{"date": "2018-10-02", "url": "https://biit.cs.ut.ee/gprofiler_archive3/e94_eg41_p11"},
               {"date": "2017-10-19", "url": "https://biit.cs.ut.ee/gprofiler_archive3/r1741_e90_eg37/web"},
               {"date": "2025-01-01", "url": "https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19"}]
    assert T.gprofiler_archives(listing) == {94: "https://biit.cs.ut.ee/gprofiler_archive3/e94_eg41_p11",
                                             113: "https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19"}


class FakeTruth:
    """Releases 114, 115, 116; gene G1 lives through, G2 is retired at 116 (into G3), G4 is born at 116."""

    def __init__(self):
        self.present = {"G1": {114, 115, 116}, "G2": {114, 115}, "G3": {114, 115, 116}, "G4": {116}}
        self.genes = {114: {}, 115: {}, 116: {}}
        self.succ = collections.defaultdict(dict, {"s116": {"G2": {"G3"}}})

    def sessions_between(self, lo, hi):
        return ["s116"] if lo < 116 <= hi else []

    def walk_series(self, gid, a, targets, strict):
        out = {}
        for t in targets:
            x = {gid}
            if a < 116 <= t:
                x = set().union(*(self.succ["s116"].get(g, {g}) for g in x))
            out[t] = {g for g in x if t in self.present.get(g, ())}
        return out


def release(names=None, acc=None, synonyms=None, has_synonyms=True):
    r = LL.Release.__new__(LL.Release)
    r.names = collections.defaultdict(set, names or {})
    r.acc = collections.defaultdict(set, acc or {})
    r.synonyms = collections.defaultdict(set, synonyms or {})
    r.has_xref, r.has_synonyms, r.files = True, has_synonyms, []
    return r


def test_the_rule_of_equal_truth():
    tr = FakeTruth()
    a115, a116 = release({"x": {"G1"}}), release({"x": {"G1"}})
    assert LL.judge(line("q", "G1", "a", from_release=114, expect=["G1"]), tr, a115, a116, 116) == (True, "equal")
    eq, why = LL.judge(line("q", "G2", "a", from_release=114, expect=["G2"]), tr, a115, a116, 116)
    assert not eq and "not the same" in why  # the expected gene is retired between 115 and 116
    eq, why = LL.judge(line("q", "X", "c", expect=["G1"]), tr, a115, release({"x": {"G1", "G4"}}), 116)
    assert not eq and "names other genes" in why  # the string names a new gene at 116
    assert LL.judge(line("q", "X", "c", expect=["G1"]), tr, a115, a116, 116)[0]
    eq, why = LL.judge(line("q", "X", "d", expect=["G1"]), tr, a115, release({"x": {"G1"}}, has_synonyms=False), 116)
    assert not eq and "no synonyms" in why  # a synonym line needs the synonyms of both releases
    eq, why = LL.judge(line("q", "NM_1.2", "e", expect=["G1"]), tr, release(acc={"NM_1": {"G1"}}),
                       release(acc={"NM_1": {"G3"}}), 116)
    assert not eq and "accession" in why
    # the same rule toward an earlier release (a tool at 114)
    assert LL.judge(line("q", "G1", "b", from_release=115, expect=["G1"]), tr, a115, a115, 114) == (True, "equal")


def test_allocation_is_equal_over_cells_and_takes_small_cells_whole():
    k = GS.allocate({"a": 5, "b": 1000, "c": 1000}, 100)
    assert k["a"] == 5 and k["b"] + k["c"] == 95 and abs(k["b"] - k["c"]) <= 1
    assert sum(GS.allocate({"a": 3, "b": 4}, 100).values()) == 7
    small, large = GS.allocate({"a": 50, "b": 500}, 40), GS.allocate({"a": 50, "b": 500}, 400)
    assert all(small[c] <= large[c] for c in small)  # the web lines are a part of the local ones


def test_a_g_line_is_judged_at_the_gene():
    src = dict(qid="a1|X|90|115", family="a", stratum="A1", kind="a1", id="ENSSSCG00000000001", from_release=90,
               to_release=115, cmp="versioned", expect=["ENSSSCG00000000002.3", "ENSSSCG00000000002.4"], broad=None,
               set="a1_v1", support={"ENSSSCG00000000002.3": 2})
    g = GS.to_g(src, "a1", "a", "web")
    assert g["cmp"] == "base" and g["expect"] == ["ENSSSCG00000000002"] and g["task"] == "a"
    assert g["source_cmp"] == "versioned" and g["stratum"] == "g_a1_a1" and g["support"] == {"ENSSSCG00000000002": 2}
    assert GS.task_of(dict(src, from_release=115), "a") == "b"


def test_the_scorer_keeps_a_line_a_tool_cannot_pose_as_not_supported():
    from idtrack_reproducibility.families import scorer as S
    ln = dict(qid="q", id="X", kind="k", rule="set", cmp="base", expect=["G1"], strategy="both", to_release=115)
    ans = dict(outcome="not_supported", got=[], error="the tool has no data for the organism")
    assert S.score_line(ln, "all", ans, None) == ("not_supported", "the tool has no data for the organism")


def test_the_best_reading_removes_a_version_a_tool_does_not_take():
    assert T.value(line("q", "ENSSSCG00000000007.3", "b"), "best") == "ENSSSCG00000000007"
    assert T.value(line("q", "ENSSSCG00000000007.3", "b"), "stated") == "ENSSSCG00000000007.3"
    assert T.value(line("q", "NM_001001258.1", "e"), "best") == "NM_001001258"
    assert T.value(line("q", "A0A024BTL2.45", "e"), "best") == "A0A024BTL2"
    assert T.value(line("q", "ABC.1", "c"), "best") == "ABC.1"  # a name is asked as it stands
    assert T.value(line("q", "A5GFW5-2", "e"), "best") == "A5GFW5-2"  # an isoform is not a version
    ls = [line("q1", "ENSSSCG00000000007.3", "b"), line("q2", "ENSSSCG00000000008", "b")]
    assert [q["qid"] for q in T.in_form(ls, "best")] == ["q1"] and len(T.in_form(ls, "stated")) == 2


def test_each_form_is_read_from_its_own_requests(tmp_path, monkeypatch):
    st = store_with(tmp_path, "gprofiler", [
        ("convert@stated", ["q1", "q2"], {"result": [{"incoming": "ENSSSCG00000000007.3", "converted": "None"},
                                                     {"incoming": "ENSSSCG00000000008", "converted": "ENSSSCG00000000008"}]}, 200),
        ("convert@best", ["q1"], {"result": [{"incoming": "ENSSSCG00000000007", "converted": "ENSSSCG00000000007"}]}, 200)])
    monkeypatch.setattr(K.RawStore, "of", classmethod(lambda cls, tool, version, organism: st))
    ls = [line("q1", "ENSSSCG00000000007.3", "b"), line("q2", "ENSSSCG00000000008", "b")]
    best, stated = A.answers("gprofiler", "sus_scrofa", "v", ls, form="best"), A.answers("gprofiler", "sus_scrofa", "v", ls, form="stated")
    assert best["q1"].got == ["ENSSSCG00000000007"] and stated["q1"].got == []
    assert best["q2"].got == stated["q2"].got == ["ENSSSCG00000000008"]  # an unversioned line shares the stated request


def test_a_line_of_another_release_is_kept_only_by_the_rule_of_equal_truth():
    import pandas as pd
    from idtrack_reproducibility.families.g import report as GR
    o = pd.DataFrame(dict(qid=["q1", "q2", "q3", "q4"], release=[116, 116, 115, None]))
    lines = {"q1": dict(to_release=115), "q2": dict(to_release=115), "q3": dict(to_release=115), "q4": dict(to_release=115)}
    lf = {116: {"q1": (True, "equal"), "q2": (False, "a gene of the truth is not the same at 116: G")}}
    kept, why = GR.kept_mask(o, lines, lf)
    assert kept == [True, False, True, True] and why[1] == "a gene of the truth is not the same at 116"


def test_each_subject_is_posed_by_its_own_rule():
    from idtrack_reproducibility.families.g import score as GSc
    sym = line("q", "TP53", "c")
    assert GSc.posable("ensembl_rest_lookup", line("q", "ENS", "b"), "homo_sapiens")[0] is None
    assert GSc.posable("gget", line("q", "ENSG00000000003", "b", from_release=115), "homo_sapiens") == (116, "")
    assert GSc.posable("orgdb", sym, "sus_scrofa")[0] is None  # org.Ss.eg.db has no Ensembl genes
    assert GSc.posable("ensdb", sym, "sus_scrofa") == (113, "")
    assert GSc.tier_of("seurat") == "web" and GSc.tier_of("pyensembl") == "local"


def stored_files(folder):
    import os
    return {f: open(os.path.join(folder, f), "rb").read() for f in os.listdir(folder)
            if f not in ("index.parquet",) and not f.endswith(".tmp")}


def test_failed_refused_and_retried_requests_keep_every_stored_response(tmp_path):
    """Every response stays stored as it came, failed and refused ones too, across retries and later passes; a later
    pass only adds (the store's files before are all there after, byte for byte; index.jsonl only grows)."""
    folder = str(tmp_path / "store")
    script = [FakeResponse(500, b"server error 1")] * (len(K.RETRY_WAITS) + 1) + \
        [FakeResponse(429, b"refused", {"Retry-After": "1"}), FakeResponse(400, b'{"error":"Too Many Hits"}'),
         FakeResponse(200, b'[{"query":"A"}]'), FakeResponse(200, b'[{"query":"B"}]')]
    c, sess, clock = client_with(script)
    st = K.RawStore(folder, "mygene", "v", "sus_scrofa")
    a = c.fetch(st, "POST", "https://mygene.info/v3/query", data={"q": "A,B"}, qids=["q1", "q2"], endpoint="query@best")
    assert a["error"] == "HTTP 500"
    before = stored_files(folder)
    # the next pass asks the same request again (an error is not an answer): refused, then a 400 of the whole batch
    st2 = K.RawStore(folder, "mygene", "v", "sus_scrofa")
    b = c.fetch(st2, "POST", "https://mygene.info/v3/query", data={"q": "A,B"}, qids=["q1", "q2"], endpoint="query@best")
    assert b["http_status"] == 400 and not K.usable(b) and st2.answered("query@best") == set()
    after = stored_files(folder)
    assert all(after.get(f) == v or f in ("index.jsonl", "headers.jsonl") for f, v in before.items())
    assert after["index.jsonl"].startswith(before["index.jsonl"])  # the index only grows
    # a pass of other batches answers the lines; the refused batch's response is still stored
    for q, t in (("q1", "A"), ("q2", "B")):
        c.fetch(st2, "POST", "https://mygene.info/v3/query", data={"q": t}, qids=[q], endpoint="query@best")
    final = stored_files(folder)
    assert set(after) <= set(final) and all(final[f] == v or f in ("index.jsonl", "headers.jsonl") for f, v in after.items())
    assert st2.answered("query@best") == {"q1", "q2"}
    bodies = {final[r["file"]] for r in (json.loads(x) for x in final["index.jsonl"].decode().splitlines()) if r.get("file")}
    assert {b"server error 1", b'{"error":"Too Many Hits"}', b'[{"query":"A"}]'} <= bodies


def test_an_interrupted_download_keeps_its_part(tmp_path):
    class Stream:
        status_code, headers = 200, {"Last-Modified": "x"}

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def iter_content(self, n):
            yield b"partial bytes"
            raise ConnectionError("cut")

    class Sess:
        def get(self, *a, **k):
            return Stream()

    c = K.Client(session_factory=Sess, sleep=lambda s: None)
    st = K.RawStore(str(tmp_path), "pyensembl", "115", "sus_scrofa")
    row = K.download(st, "https://ftp.ensembl.org/pub/x.gtf.gz", client=c)
    assert row["error"].startswith("ConnectionError") and row["file"].endswith(".part")
    assert open(tmp_path / row["file"], "rb").read() == b"partial bytes"


def test_scoring_reads_the_store_with_the_network_off(tmp_path, monkeypatch):
    """The adapters read stored responses only: with every socket refused, the answers come from the store."""
    import socket

    def refuse(*a, **k):
        raise OSError("network off")
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    st = store_with(tmp_path, "gprofiler", [("convert@stated", ["q1"], {"result": [
        {"incoming": "TP53", "converted": "ENSSSCG00000000009"}]}, 200)])
    monkeypatch.setattr(K.RawStore, "of", classmethod(lambda cls, tool, version, organism: st))
    out = A.answers("gprofiler", "sus_scrofa", "v", [line("q1", "TP53", "c")], form="stated")
    assert out["q1"].got == ["ENSSSCG00000000009"]
