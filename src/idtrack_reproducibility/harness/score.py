"""Score a clone of idtrack on a graph against the fixed query sets of the accuracy harness.

Usage:
  score.py --clone <clone> --graph <graph file or folder> --organism <name> --out <folder>
           [--sets v1,v2] [--kinds a,b] [--timeout 30] [--release 115] [--shard i/n]
  score.py ... --submit 16      submit 16 shard jobs and a merge job (the configured runner), then exit
  score.py ... --merge          only score the results already in <out>/results/ (all shards)

The graph is loaded through the clone's public API (`idtrack.API(local_repository=...)`, `build_graph`). A graph
file is linked into `<out>/repo/` under its own name, so nothing is written beside the original; a folder is used
as the local repository as it is. Every query runs `API.convert_identifier` once with a time limit (`--timeout`
seconds; over it the query is `timeout`). Results are appended to `<out>/results/shard_<i>_of_<n>.jsonl` (resumable).
When all shards are complete (or without --shard), the results are scored and written:

  scores.tsv      kind, variant (`*` = all variants of the kind; else `<strategy>,<from|default>`), n, right, wrong, undetermined, timeout, error, rate, ci_low, ci_high
                  rate = right / (right + wrong + error); undetermined and timeout are left out of the rate
  outcomes/part_0_of_1.parquet   the record (docs/records.md): one row per query with the query, the truth, the
                  answer and the verdict, in the columns of the common scorer (families/scorer.py, `OUT_COLS`)
  outcomes.tsv.gz one line per query: qid, kind, variant, outcome, reason, got, seconds (derived from the record; what
                  compare_scores.py and travel_report.py read)
  wrong.tsv       every wrong answer: qid, kind, variant, query, expected, broad, got, reason
  errors.tsv      every exception: qid, kind, variant, query, error
  timeouts.tsv    every query over the limit
  undetermined.tsv  every answer that Ensembl's data neither confirm nor refute, with the reason
  run.json        the record's version, the commit of the checkout that answered (from the answers' load records),
                  the graph, sets and their sha256, time limit, times
"""
import argparse
import datetime
import gc
import glob
import gzip
import json
import os
import re
import signal
import socket
import subprocess
import sys
import time
import traceback
import zlib

from idtrack_reproducibility.harness import harness_lib as H

ENS = re.compile(rf"^({H.STABLE_ID_FORM})\.\d+$")


class QueryTimeout(BaseException):  # not an Exception: the package must not swallow it
    pass


def _alarm(signum, frame):
    raise QueryTimeout()


def load_queries(org, sets, kinds, judge="revised", within=None):
    """The latest line of every qid over the sets named (default: every published set, in order). With judge
    `published`, the revision sets of the corrected judge and of the judge following idtrack's decided rules (manifest
    `judge` `78` or `rules`: harness_lib.REVISED_JUDGES) are left out: the judge as it was before. `within` (set
    names): only the qids these sets publish (the kinds of a run, scored with revision sets that also revise the
    travel matrix's lines)."""
    d = H.query_dir(org)
    if not sets:
        sets = H.set_names(d)
        if judge == "published":
            sets = [x for x in sets if json.load(open(os.path.join(d, f"{x}.manifest.json"))).get("judge")
                    not in H.REVISED_JUDGES]
    by_qid, man, bound = {}, {}, set()
    for s in sets:
        p = os.path.join(d, f"{s}.jsonl.gz")
        m = json.load(open(os.path.join(d, f"{s}.manifest.json")))
        man[s] = m["sha256"]
        with gzip.open(p, "rt") as fh:
            for line in fh:
                rec = json.loads(line)
                if kinds and rec["kind"] not in kinds:
                    continue
                rec["set"] = s
                # A set that revises the truth of published queries (make_queries.py --copies-apart) repeats their
                # qid: the line of the later set is the truth, at the place of the earlier line. A query is asked
                # once, and its stored answer is scored under whichever sets are named with --sets.
                by_qid[rec["qid"]] = rec
                if within and s in within:
                    bound.add(rec["qid"])
    if within:
        by_qid = {q: r for q, r in by_qid.items() if q in bound}
    return list(by_qid.values()), man


def variant(rec):
    return f"{rec['strategy']},{'from' if rec['from_release'] is not None else 'default'}"


# ---------------- running ----------------

def targets_of(res, rec):
    """The ids the API answered, in the requested database. One place to adapt when the result format changes."""
    if res.get("no_corresponding") or res.get("no_conversion"):
        return []
    tids = res.get("target_id") or []
    db = rec.get("final_database")
    out = []
    for t in tids:
        if isinstance(t, dict):  # per-target records
            if db and t.get("database") not in (None, db):
                continue
            if t.get("no_target") or t.get("fallback"):
                continue
            out.append(str(t.get("id")))
        else:
            out.append(str(t))
    if db and (res.get("no_target") is True or res.get("final_database") not in (None, db)):
        return []
    return out


def load_api(clone, graph, org, release, out):
    sys.path.insert(0, clone)
    import idtrack  # noqa: E402

    if os.path.isdir(graph):
        repo = graph
    else:
        repo = os.path.join(out, "repo")
        os.makedirs(repo, exist_ok=True)
        link = os.path.join(repo, os.path.basename(graph))
        try:  # shards start at the same time: the first one makes the link
            os.symlink(os.path.abspath(graph), link)
        except FileExistsError:
            pass
        if os.path.realpath(link) != os.path.realpath(graph):
            raise SystemExit(f"{link} points to {os.path.realpath(link)}, not to {graph}: use another --out")
    before = set(os.listdir(repo))
    api = idtrack.API(local_repository=repo)
    t0 = time.time()
    api.build_graph(organism_name=org, snapshot_release=release)
    gc.collect()
    gc.freeze()  # a full collection over the loaded graph pauses for seconds; it would be counted as query time
    new = sorted(set(os.listdir(repo)) - before)
    return api, dict(idtrack_file=idtrack.__file__, pythonhashseed=os.environ.get("PYTHONHASHSEED"), load_seconds=round(time.time() - t0, 1), repo=repo,
                     files_written_in_repo=new)


def run(a, queries):
    i, n = (int(x) for x in a.shard.split("/"))
    # a query's shard depends only on its qid, so adding sets or kinds never moves a query to another shard
    mine = [q for q in queries if zlib.crc32(q["qid"].encode()) % n == i]
    rdir = os.path.join(a.out, "results")
    os.makedirs(rdir, exist_ok=True)
    path = os.path.join(rdir, f"shard_{i}_of_{n}.jsonl")
    done = set()  # results of every earlier run in this folder, whatever its sharding
    for p in glob.glob(os.path.join(rdir, "shard_*.jsonl")):
        with open(p) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue  # a line cut by an interrupted run is rerun
                done.add(r["qid"])
    todo = [q for q in mine if q["qid"] not in done]
    sys.stderr.write(f"shard {i}/{n}: {len(mine)} queries, {len(todo)} to run\n")
    if not todo:
        return
    api, info = load_api(a.clone, a.graph, a.organism, a.release, a.out)
    json.dump(info, open(os.path.join(rdir, f"load_{i}_of_{n}.json"), "w"), indent=1)
    signal.signal(signal.SIGALRM, _alarm)
    t_all = time.time()
    with open(path, "a") as fo:
        for k, rec in enumerate(todo):
            r = dict(qid=rec["qid"])
            t0 = time.time()
            try:
                signal.setitimer(signal.ITIMER_REAL, a.timeout)
                try:
                    res = api.convert_identifier(rec["id"], from_release=rec["from_release"], to_release=rec["to_release"],
                                                 final_database=rec["final_database"], strategy=rec["strategy"])
                finally:
                    signal.setitimer(signal.ITIMER_REAL, 0)
                r["got"] = sorted(set(targets_of(res, rec)))
                if "alternative_locus_genes" in res:  # the copies reported apart (trees that have rule `copies_apart`)
                    r["apart"] = sorted(str(x) for x in res["alternative_locus_genes"])
                r["flags"] = {f: bool(res.get(f)) for f in ("no_corresponding", "no_conversion", "no_target") if f in res}
                r["final_database"] = res.get("final_database")
            except QueryTimeout:
                r["timeout"] = a.timeout
            except Exception as e:  # noqa: BLE001  a crash on a valid query is a result
                r["error"] = f"{type(e).__name__}: {e}"[:500]
                r["traceback"] = traceback.format_exc()[-1200:]
            r["seconds"] = round(time.time() - t0, 4)
            fo.write(json.dumps(r, sort_keys=True) + "\n")
            fo.flush()
            if k % 2000 == 0:
                sys.stderr.write(f"{k}/{len(todo)} {time.time() - t_all:.0f}s\n")
    sys.stderr.write(f"shard {i}/{n} done in {time.time() - t_all:.0f}s\n")


# ---------------- scoring ----------------

class Places:
    """Gene places at a release, from Ensembl's GTF (for the locus rule); with `projected`, a place on a haplotype,
    patch, alternative locus or Y PAR put on its reference sequence (assembly_map.Exceptions)."""

    def __init__(self, org):
        self.org, self.cache, self.asm, self._exc = org, {}, None, None

    def core(self, rel, gid):
        """The place of a gene of the core database that the GTF of the release lacks (`core_places/<release>.tsv`,
        fetch_tables.py core_places), or None."""
        key = ("core", rel)
        if key not in self.cache:
            p = os.path.join(H.truth_dir(self.org), "core_places", f"{rel}.tsv")
            self.cache[key] = {r["gene_id"]: r for r in H.read_tsv(p)} if os.path.exists(p) else {}
        return self.cache[key].get(gid)

    def projected(self, rel, g):
        if g is None:
            return None
        if self._exc is None:
            from idtrack_reproducibility.harness import assembly_map
            self._exc = assembly_map.Exceptions(self.org)
        return self._exc.project(g, rel)

    def assembly(self, rel):
        if self.asm is None:
            self.asm = H.gtf_assemblies(self.org)
        return self.asm.get(rel)

    def get(self, rel, gid):
        if rel not in self.cache:
            self.cache[rel] = H.gtf_genes(self.org, rel) or {}
        return self.cache[rel].get(gid)


def norm(ids, cmp):
    if cmp == "base":
        return {ENS.sub(lambda m: m.group(1), x) if ENS.match(x) else x for x in ids}
    return set(ids)


def judged_by(rec):
    """The judge of a line: `published`, or `revised:<corrections of the corrected judge that changed it>`."""
    return "revised:" + "+".join(rec["revised_by"]) if rec.get("revised_by") else "published"


def _place_rule(loc, got, B, rec, places, mapped=None):
    """The locus rule on one assembly: every returned gene outside `broad` overlaps the place -> undetermined.

    Lines revised by `exceptions` (the corrected judge, 2026-11-18) compare places on the reference sequence: the query's
    place as the line gives it (`locus_projected`, or `locus_mapped` from it) and each returned gene's place at the
    target release projected through Ensembl's assembly exceptions (`Places.projected`)."""
    from idtrack_reproducibility.harness import assembly_map
    rest = got - (B or set())
    far = []
    project = bool(rec.get("exceptions"))
    moved = bool(loc.get("projected_from"))
    if loc.get("unplaced"):
        return "undetermined", f"Ensembl has no successor; the query's place has no place on the reference sequence " \
                               f"({loc['unplaced']})"
    open_ = []
    for t in sorted(rest):
        gid = ENS.sub(lambda m: m.group(1), t) if ENS.match(t) else t
        g = places.get(rec["to_release"], gid)
        if g is None and rec.get("core_places"):  # a gene of the core that the GTF lacks (the corrected judge, core_places)
            g = places.core(rec["to_release"], gid)
        if project:
            g = places.projected(rec["to_release"], g)
            moved = moved or bool(g and g.get("projected_from"))
            if g and g.get("unplaced"):
                open_.append(t)  # the answer's own place cannot be put on the reference: not judged by place
                continue
        if not assembly_map.overlaps(loc, g):
            far.append(t)
    where = f" (place mapped {mapped['via']}, {mapped['aligned']:.0%} aligned)" if mapped else ""
    if moved:
        where += " (places on haplotypes, patches or PARs projected on the reference sequence)"
    if not far:
        return "undetermined", "Ensembl has no successor; every returned gene overlaps the query's place" + where + (
            f" or has no place on the reference sequence ({','.join(open_)})" if open_ else "")
    return "wrong", f"Ensembl has no successor; not at the query's place{where}: {','.join(far)}"


def score_one(rec, res, places):
    """Score one stored answer; then, for a right answer, the copies the truth expects to be reported apart."""
    outcome, why = score_answer(rec, res, places)
    if outcome != "right" or not rec.get("apart"):
        return outcome, why
    if "apart" not in res:  # an answer stored by a tree, or by a run, that has no list of copies: nothing to check
        return outcome, APART_UNCHECKED
    lacking = sorted(set(rec["apart"]) - norm(res["apart"], rec["cmp"]))
    if lacking:
        return "wrong", f"copies not reported apart: {','.join(lacking[:5])}"
    return outcome, why


APART_UNCHECKED = "copies not checked: the stored answer has no list of them"


def score_answer(rec, res, places):
    if res is None:
        return "missing", "not run"
    if "timeout" in res:
        return "timeout", f"over {res['timeout']} s"
    if "error" in res:
        return "error", res["error"]
    got = norm(res["got"], rec["cmp"])
    E = set(rec["expect"])
    B = set(rec["broad"]) if rec.get("broad") is not None else None
    if rec.get("spellings"):  # the corrected judge: two spellings of one value of one database are one value
        one = {v: grp[0] for grp in rec["spellings"] for v in grp}
        got, E = {one.get(x, x) for x in got}, {one.get(x, x) for x in E}
        B = {one.get(x, x) for x in B} if B is not None else None
    exact = rec["strategy"] == "all" or rec["rule"] == "set_exact"
    if exact:
        if got == E:
            return "right", ""
        sup = rec.get("support")
        if sup and rec["strategy"] == "best" and got and got < E:
            # rule `support_ranks_shared_name` (the corrected judge): the returned genes carry the name in strictly more name sources than
            # every expected gene left out; the truth cannot decide between the best-supported gene and the set
            if min(sup.get(g, 0) for g in got) > max(sup.get(g, 0) for g in E - got):
                return "undetermined", "rule `support_ranks_shared_name`: the returned genes carry the name in more name sources than " \
                                       "every expected gene left out"
        if B is not None and E <= got <= B:
            return "undetermined", "within the wider set that Ensembl's history also allows"
    else:
        if (E and got and got <= E) or (not E and not got):
            return "right", ""
        if B is not None and got and got <= B:
            return "undetermined", "within the wider set that Ensembl's history also allows"
    if not got:
        return "wrong", "no answer" + (" (id not found)" if res.get("flags", {}).get("no_corresponding") else "")
    if not E and rec.get("locus"):
        loc = rec["locus"]
        if places.assembly(rec["to_release"]) != loc["assembly"]:
            m = rec.get("locus_mapped")  # the corrected judge: the place mapped through Ensembl's assembly chain files
            if m and "unmapped" not in m:
                return _place_rule(m, got, B, rec, places, mapped=m)
            why = f"; {m['unmapped']}" if m else ""
            return "undetermined", f"Ensembl has no successor; places not comparable ({loc['assembly']} vs {places.assembly(rec['to_release'])}){why}"
        return _place_rule(rec.get("locus_projected") or loc if rec.get("exceptions") else loc, got, B, rec, places)
    extra = sorted(got - E)
    missing = sorted(E - got)
    return "wrong", (f"not expected: {','.join(extra[:5])}" if extra else "") + (" " if extra and missing else "") + (
        f"missing: {','.join(missing[:5])}" if missing else "")


def write_tsv(path, header, rows):
    with open(path + ".tmp", "w") as fo:
        fo.write("\t".join(header) + "\n")
        for r in rows:
            fo.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
    os.replace(path + ".tmp", path)


def status_of(r):
    """The answer's status in the record's words (docs/records.md) from the d4c2dda answer's flags."""
    if not r or "got" not in r:
        return None
    f = r.get("flags") or {}
    if f.get("no_conversion"):
        return "no_conversion"
    if f.get("no_corresponding"):
        return "not_found"
    if f.get("no_target"):
        return "no_target"
    return "converted"


def harness_key(a, man):
    """What makes a scoring of this scorer (as `scorer.scoring_key`): the organism, judge, sets with their checksums,
    kinds, limit, the record version and the code of the judge; `rescore` leaves a scoring of the same key as it is."""
    import hashlib
    from idtrack_reproducibility.families import lines as L
    from idtrack_reproducibility.families import scorer as SC
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc = dict(scorer="harness", organism=a.organism, judge=a.judge, sets=man, kinds=sorted(a.kinds.split(",")) if a.kinds else None,
               within=a.within, limit=a.limit, record_version=SC.RECORD_VERSION,
               code={m: L.sha256(os.path.join(here, m)) for m in SC.JUDGE_MODULES})
    return hashlib.sha256(json.dumps(doc, sort_keys=True).encode()).hexdigest()


def merge_and_score(a, queries, man):
    res = {}
    for p in sorted(glob.glob(os.path.join(a.out, "results", "shard_*.jsonl"))):
        with open(p) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                res[r["qid"]] = r
    places = Places(a.organism)
    tab = {}
    wrong, errors, timeouts, und = [], [], [], []
    missing = 0
    apart = dict(expected=0, checked=0, unchecked=0)
    stored, a.out = a.out, (a.scores_to or a.out)  # the tables of this scoring; the stored answers stay where they are
    os.makedirs(os.path.join(a.out, "outcomes"), exist_ok=True)
    import pyarrow as pa
    import pyarrow.parquet as pq
    from idtrack_reproducibility.families import lines as FL
    from idtrack_reproducibility.families import scorer as SC
    record, rows_path = [], os.path.join(a.out, "outcomes", "part_0_of_1.parquet")
    with gzip.open(os.path.join(a.out, "outcomes.tsv.gz.tmp"), "wt") as fo:
        fo.write("qid\tkind\tvariant\toutcome\treason\tgot\tseconds\trule\n")
        for rec in queries:
            r = res.get(rec["qid"])
            o, why = score_one(rec, r, places)
            if o == "missing":
                missing += 1
                continue
            if rec.get("apart"):
                apart["expected"] += 1
                if why == APART_UNCHECKED:
                    apart["unchecked"] += 1
                elif "apart" in r:
                    apart["checked"] += 1
            v = variant(rec)
            got = ",".join(r.get("got", [])) if r else ""
            rule = judged_by(rec)
            fo.write("\t".join([rec["qid"], rec["kind"], v, o, why.replace("\t", " "), got, str(r.get("seconds", "")), rule]) + "\n")
            base = dict(qid=rec["qid"], view=rec["strategy"], family="harness", stratum=rec["kind"], kind=rec["kind"], variant=v,
                        unit=FL.unit_of(rec), set=rec.get("set"),
                        rule=rule)
            record.append(SC._row(base, rec, o, why, list(r.get("got") or []) if r else [], r.get("seconds") if r else None,
                                  1 if r else 0, dict(status=status_of(r), match=None,
                                                      error_type=(r.get("error") or "").split(":", 1)[0] or None if r else None)))
            for key in ((rec["kind"], "*"), (rec["kind"], v)):
                c = tab.setdefault(key, dict(right=0, wrong=0, undetermined=0, timeout=0, error=0))
                c[o] += 1
            qtxt = H.dumps({k: rec[k] for k in ("id", "from_release", "to_release", "final_database", "strategy")})
            if o == "wrong":
                wrong.append([rec["qid"], rec["kind"], v, qtxt, ",".join(rec["expect"]), ",".join(rec["broad"] or []), got, why, rule])
            elif o == "error":
                errors.append([rec["qid"], rec["kind"], v, qtxt, why])
            elif o == "timeout":
                timeouts.append([rec["qid"], rec["kind"], v, qtxt, r.get("seconds")])
            elif o == "undetermined":
                und.append([rec["qid"], rec["kind"], v, qtxt, ",".join(rec["expect"]), got, why, rule])
    os.replace(os.path.join(a.out, "outcomes.tsv.gz.tmp"), os.path.join(a.out, "outcomes.tsv.gz"))
    pq.write_table(pa.Table.from_pylist(record, schema=SC.out_schema()), rows_path + ".tmp", compression="zstd")
    os.replace(rows_path + ".tmp", rows_path)
    n_record = len(record)
    del record
    rows = []
    for (k, v), c in sorted(tab.items(), key=lambda kv: (kv[0][0], kv[0][1] != "*", kv[0][1])):
        n = sum(c.values())
        den = c["right"] + c["wrong"] + c["error"]
        lo, hi = H.wilson(c["right"], den)
        rate = c["right"] / den if den else float("nan")
        rows.append([k, v, n, c["right"], c["wrong"], c["undetermined"], c["timeout"], c["error"], f"{rate:.4f}",
                     f"{lo:.4f}", f"{hi:.4f}"])
    write_tsv(os.path.join(a.out, "scores.tsv"),
              ["kind", "variant", "n", "right", "wrong", "undetermined", "timeout", "error", "rate", "ci_low", "ci_high"], rows)
    write_tsv(os.path.join(a.out, "wrong.tsv"), ["qid", "kind", "variant", "query", "expected", "broad", "got", "reason", "rule"], wrong)
    write_tsv(os.path.join(a.out, "errors.tsv"), ["qid", "kind", "variant", "query", "error"], errors)
    write_tsv(os.path.join(a.out, "timeouts.tsv"), ["qid", "kind", "variant", "query", "seconds"], timeouts)
    write_tsv(os.path.join(a.out, "undetermined.tsv"), ["qid", "kind", "variant", "query", "expected", "got", "reason", "rule"], und)
    loads = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(stored, "results", "load_*.json")))]
    from idtrack_reproducibility import config as CFG
    from idtrack_reproducibility.families import common as FC
    prev_path = os.path.join(a.out, "run.json")
    prev = json.load(open(prev_path)) if os.path.exists(prev_path) else {}
    # the checkout that answered (the answers' load records name it), not the one this scoring was started with
    clone = os.path.dirname(os.path.dirname(loads[0]["idtrack_file"])) if loads and loads[0].get("idtrack_file") else a.clone
    commit, dirty = FC.git_commit(clone) if os.path.isdir(clone) else ("", False)
    if not commit and prev.get("commit"):  # that checkout is gone: the commit the earlier scoring of these answers read
        commit, dirty = prev["commit"], prev.get("clone_modified")
    if loads and loads[0].get("repo"):
        a.graph = loads[0]["repo"] if os.path.exists(loads[0]["repo"]) else a.graph
    g = os.path.realpath(a.graph)
    info = dict(record_version=SC.RECORD_VERSION, key=harness_key(a, man), family="harness", n_rows=n_record,
                clone=CFG.in_root(clone, None), commit=commit, clone_modified=bool(dirty), graph=CFG.in_root(a.graph, None),
                graph_real=CFG.in_root(g, None), graph_bytes=(os.path.getsize(g) if os.path.isfile(g) else None),
                graph_mtime=(datetime.datetime.fromtimestamp(os.path.getmtime(g)).isoformat(timespec="seconds")),
                organism=a.organism, release=a.release, sets=man, timeout=a.timeout, n_queries=len(queries),
                n_excluded=getattr(a, "n_excluded", 0), excluded=getattr(a, "excluded", {}), judge=a.judge,
                n_missing=missing, scored_at=datetime.datetime.now().isoformat(timespec="seconds"),
                answers=CFG.in_root(stored), run=CFG.in_root(stored), copies_apart=apart,
                load=[{k: x.get(k) for k in ("load_seconds", "files_written_in_repo", "pythonhashseed", "python")} for x in loads],
                query_seconds_total=round(sum(r.get("seconds", 0) for r in res.values()), 1))
    json.dump(info, open(os.path.join(a.out, "run.json"), "w"), indent=1)
    print(open(os.path.join(a.out, "scores.tsv")).read())
    if missing:
        print(f"WARNING: {missing} queries have no result yet (shards incomplete)")


def submit(a):
    from idtrack_reproducibility import config as CFG
    from idtrack_reproducibility import jobs
    logs = os.path.join(a.out, "logs")
    common = (f"--clone {a.clone} --graph {a.graph} --organism {a.organism} --out {a.out} --timeout {a.timeout} "
              f"--release {a.release} --judge {a.judge}" + (f" --sets {a.sets}" if a.sets else "") + (f" --kinds {a.kinds}" if a.kinds else ""))
    mod = "idtrack_reproducibility.harness.score"
    py = CFG.get().python_for(a.clone)
    mem = a.mem or jobs.memory_for(a.organism)
    tag = a.tag or f"harness_{a.organism}"
    ids = [jobs.submit(f"{tag}_{i}", jobs.module(py, mod, f"{common} --shard {i}/{a.submit}"), logs, cpus=2, mem_gb=mem,
                       time=a.time) for i in range(a.submit)]
    m = jobs.submit(f"{tag}_merge", jobs.module(CFG.get().python, mod, f"{common} --merge"), logs, cpus=1, mem_gb=16,
                    time="02:00:00", after=ids)
    print("shards", " ".join(ids), "merge", m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sets", default=None, help="comma-separated query sets (default: every published set)")
    ap.add_argument("--kinds", default=None)
    ap.add_argument("--within", default=None, help="only the qids these sets publish (default: every line of the sets)")
    ap.add_argument("--judge", default="revised", choices=["revised", "published"],
                    help="without --sets: `published` leaves out the revision sets of the corrected judge")
    ap.add_argument("--limit", type=int, default=None, help="smoke test: only the first N queries of each kind")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--release", type=int, default=None, help="snapshot release of the graph (default: from its name)")
    ap.add_argument("--shard", default=None)
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--scores-to", default=None, help="with --merge: write the tables of this scoring into another "
                    "folder (the stored answers of --out are scored again, e.g. under other --sets, and stay as they are)")
    ap.add_argument("--submit", type=int, default=0)
    ap.add_argument("--tag", default=None, help="with --submit: the jobs' name (default harness_<organism>)")
    ap.add_argument("--mem", type=int, default=None)
    ap.add_argument("--time", default="12:00:00")
    a = ap.parse_args()
    a.out = os.path.abspath(a.out)
    a.scores_to = os.path.abspath(a.scores_to) if a.scores_to else None
    a.clone = os.path.abspath(a.clone)
    if a.release is None:
        m = re.search(r"max(\d+)", os.path.basename(os.path.realpath(a.graph)) if os.path.isfile(a.graph) else
                      " ".join(os.listdir(a.graph)))
        a.release = int(m.group(1)) if m else H.ORGANISMS[a.organism][2]
    os.makedirs(a.out, exist_ok=True)
    if a.submit:
        return submit(a)
    a.within = a.within.split(",") if a.within else None
    queries, man = load_queries(a.organism, a.sets.split(",") if a.sets else None, set(a.kinds.split(",")) if a.kinds else None,
                                a.judge, a.within)
    n_all = len(queries)
    a.excluded = {}
    for q in queries:
        why = "harness_lib.EXCLUDED" if (a.organism, q["kind"], q["final_database"]) in H.EXCLUDED else (
            "excluded by " + judged_by(q) if q.get("excluded") else None)
        if why:
            a.excluded[f"{q['kind']}: {why}"] = a.excluded.get(f"{q['kind']}: {why}", 0) + 1
    queries = [q for q in queries if (a.organism, q["kind"], q["final_database"]) not in H.EXCLUDED and not q.get("excluded")]
    a.n_excluded = n_all - len(queries)
    if a.limit:
        per = {}
        queries = [q for q in queries if per.setdefault(q["kind"], []).append(1) or len(per[q["kind"]]) <= a.limit]
    if not a.merge:
        if a.shard is None:
            a.shard = "0/1"
        run(a, queries)
        if a.shard != "0/1":
            return
    merge_and_score(a, queries, man)


if __name__ == "__main__":
    main()
