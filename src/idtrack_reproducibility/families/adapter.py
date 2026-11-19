"""The adapter (F0): ask one query line of any family through one idtrack tree, whatever its api generation.

    a = Adapter(clone, graph_file, organism, release, workdir)   # loads the graph once (gc.freeze after)
    answers = a.ask(line, explain=False)                          # one answer per view of the line

Generations (`common.generation`):

- `d4c2dda`: `API(local_repository=...).build_graph(organism, snapshot_release)`, then `API.convert_identifier`
  once per view (`best`, `all`): a line with `strategy: both` costs two calls. The ids answered are read by the
  harness's own `score.targets_of` (nothing when `no_corresponding` / `no_conversion`; with a `final_database`,
  nothing when `no_target` or another database).
- `convert` (71690ab and later): `API.load(graph)`, then one `API.convert` call; `result.targets(view)` for each
  view, with the same reading (nothing for `not_found`, `unknown_version`, `no_conversion`; with a
  `final_database`, nothing for `no_target` or another database).

An answer (one per view; plain data, the runner's Parquet row): `qid`, `view`, `outcome` (`answer`, `timeout`,
`error`), `got` (the ids, sorted), `apart` (copies reported apart; None when the tree has no such list), `status`
(`converted`, `tie`, `no_target`, `no_conversion`, `not_found`, `unknown_version`; for d4c2dda derived from its
flags), `match`, `tie`, `from_release` (used), `fallback`, `final_database` (of the result), `flags`
(`no_corresponding`, `no_conversion`, `no_target` as d4c2dda names them), `error_type`, `error_mro` (class names,
for the rule `expect_error`), `error` (message), `traceback`, `seconds` (of the call that gave this view), `calls`
(calls made for the line), and with `explain` `path` (JSON: per candidate of the view `ensembl_gene`, `target_id`,
`entry_gene`, `external_jump`, `edges` [source, target, kind, release], `edge_dbs` [`<database>@<assembly>`, ...]: the
format the cause classifier reads, for both generations).

Budget: `timeout` seconds per call. d4c2dda: SIGALRM at `timeout` (the harness's limit); 71690ab on:
`budget_seconds=timeout` and SIGALRM at `timeout + 30` as a backstop. A `SearchBudgetExceededError` is the outcome
`timeout` (principle 5: over the budget), as the alarm is.
"""
import gc
import json
import math
import os
import signal
import sys
import time
import traceback

from idtrack_reproducibility.families import common as C


class QueryTimeout(BaseException):  # not an Exception: the package must not swallow it
    pass


def _alarm(signum, frame):
    raise QueryTimeout()


def decode_id(line):
    """The identifier argument of a line (`id_codec`: str, json, bytes_hex, float_nan)."""
    codec = line.get("id_codec") or "str"
    v = line["id"]
    if codec == "str":
        return v
    if codec == "json":
        return json.loads(v)
    if codec == "bytes_hex":
        return bytes.fromhex(v)
    if codec == "float_nan":
        return math.nan
    raise ValueError(f"unknown id_codec {codec}")


def link_graph(graph, workdir):
    """Link the graph file into <workdir>/repo/ (nothing is written beside the original); return the link."""
    repo = os.path.join(workdir, "repo")
    os.makedirs(repo, exist_ok=True)
    link = os.path.join(repo, os.path.basename(graph))
    try:  # shards start at the same time: the first one makes the link
        os.symlink(os.path.abspath(graph), link)
    except FileExistsError:
        pass
    if os.path.realpath(link) != os.path.realpath(graph):
        raise SystemExit(f"{link} points to {os.path.realpath(link)}, not to {graph}: use another run folder")
    return link


class Adapter:
    def __init__(self, clone, graph, organism, release, workdir, timeout=30.0):
        self.clone, self.organism, self.release, self.timeout = os.path.abspath(clone), organism, release, timeout
        self.gen = C.generation(self.clone)
        sys.path.insert(0, self.clone)
        import idtrack  # noqa: E402

        self.idtrack = idtrack
        if os.path.isdir(graph):
            raise SystemExit("--graph is the graph file (it is linked into the run folder)")
        link = link_graph(graph, workdir)
        repo = os.path.dirname(link)
        before = set(os.listdir(repo))
        t0 = time.time()
        if self.gen == "d4c2dda":
            self.api = idtrack.API(local_repository=repo)
            self.api.build_graph(organism_name=organism, snapshot_release=release)
            self.graph = self.api.track.graph
            from idtrack._track import SearchBudgetExceededError as Budget
        else:
            self.api = idtrack.API.load(link)
            self.graph = self.api.graph
            from idtrack.errors import SearchBudgetExceededError as Budget
        self.budget_error = Budget
        manifest = (getattr(self.graph, "graph", None) or {}).get("manifest") or {}
        gc.collect()
        gc.freeze()  # a full collection over the loaded graph pauses for seconds; it would be counted as query time
        signal.signal(signal.SIGALRM, _alarm)
        self.info = dict(generation=self.gen, idtrack_file=idtrack.__file__, pythonhashseed=os.environ.get("PYTHONHASHSEED"),
                         python=sys.version.split()[0], load_seconds=round(time.time() - t0, 1), graph=os.path.abspath(graph),
                         graph_real=os.path.realpath(graph), files_written_in_repo=sorted(set(os.listdir(repo)) - before),
                         content_digest=manifest.get("content_digest"), package_commit=manifest.get("package_commit"))

    # ---------------- one line ----------------

    def ask(self, line, explain=False):
        vs = ["best", "all"] if (line.get("strategy") or "best") == "both" else [line.get("strategy") or "best"]
        if self.gen == "d4c2dda":
            out = [self._ask_d4c2dda(line, v, explain) for v in vs]
        else:
            out = self._ask_convert(line, vs, explain)
        for a in out:
            a["calls"] = len(vs) if self.gen == "d4c2dda" else 1
        return out

    def _call(self, fn):
        """Run fn() under the budget; returns (result, None) or (None, failure dict)."""
        limit = self.timeout if self.gen == "d4c2dda" else self.timeout + 30
        t0 = time.perf_counter()
        try:
            signal.setitimer(signal.ITIMER_REAL, limit)
            try:
                res = fn()
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            return res, None, time.perf_counter() - t0
        except QueryTimeout:
            return None, dict(outcome="timeout", error_type="alarm", error=f"over {limit} s"), time.perf_counter() - t0
        except self.budget_error as e:
            return None, dict(outcome="timeout", error_type=type(e).__name__, error=str(e)[:500]), time.perf_counter() - t0
        except Exception as e:  # noqa: BLE001  a crash on a valid query is a result
            return None, dict(outcome="error", error_type=type(e).__name__, error_mro=[k.__name__ for k in type(e).__mro__],
                              error=f"{type(e).__name__}: {e}"[:500], traceback=traceback.format_exc()[-1200:]), \
                time.perf_counter() - t0

    @staticmethod
    def _blank(line, view):
        return dict(qid=line["qid"], view=view, outcome="answer", got=[], apart=None, status=None, match=None, graph_id=None, tie=None,
                    from_release=None, fallback=None, final_database=None, flags=None, error_type=None, error_mro=None,
                    error=None, traceback=None, seconds=None, calls=None, path=None)

    def _ask_d4c2dda(self, line, view, explain):
        from idtrack_reproducibility.harness import score  # the harness's reading of d4c2dda's answer (score.targets_of)

        kw = dict(from_release=line.get("from_release"), to_release=line.get("to_release"),
                  final_database=line.get("final_database"), strategy=view)
        if explain:
            kw["explain"] = True
        if line.get("from_database") is not None:
            kw["from_database"] = line["from_database"]
        ident = decode_id(line)
        res, fail, sec = self._call(lambda: self.api.convert_identifier(ident, **kw))
        a = self._blank(line, view)
        a["seconds"] = round(sec, 6)
        if fail:
            a.update(fail)
            return a
        a["got"] = sorted(set(score.targets_of(res, line)))
        if "alternative_locus_genes" in res:
            a["apart"] = sorted(str(x) for x in res["alternative_locus_genes"])
        a["flags"] = {f: bool(res.get(f)) for f in ("no_corresponding", "no_conversion", "no_target") if f in res}
        if res.get("no_corresponding"):
            a["status"] = "unknown_version" if res.get("known_versions") else "not_found"
        elif res.get("no_conversion"):
            a["status"] = "no_conversion"
        elif res.get("no_target"):
            a["status"] = "no_target"
        elif res.get("tie"):
            a["status"] = "tie"
        else:
            a["status"] = "converted"
        a.update(match=res.get("match"), graph_id=res.get("graph_id"), tie=bool(res.get("tie")), from_release=res.get("from_release"),
                 fallback=res.get("from_release_fallback"), final_database=res.get("final_database"))
        if explain:
            a["path"] = json.dumps(self._paths_d4c2dda(res.get("the_path") or []), default=str)
        return a

    def _ask_convert(self, line, views, explain):
        kw = dict(from_release=line.get("from_release"), to_release=line.get("to_release"),
                  final_database=line.get("final_database"), budget_seconds=self.timeout)
        if explain:
            kw["explain"] = True
        if line.get("from_database") is not None:
            kw["from_database"] = line["from_database"]
        ident = decode_id(line)
        res, fail, sec = self._call(lambda: self.api.convert(ident, **kw))
        out = []
        db = line.get("final_database")
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(fail)
                out.append(a)
                continue
            got = list(res.targets(v))
            if res.status in ("not_found", "unknown_version", "no_conversion"):
                got = []
            if db and (res.no_target or res.final_database not in (None, db)):
                got = []
            a["got"] = sorted(set(str(x) for x in got))
            a["apart"] = sorted(str(x) for x in res.alternative_locus_genes)
            a["flags"] = dict(no_corresponding=res.status in ("not_found", "unknown_version"),
                              no_conversion=res.status == "no_conversion", no_target=bool(res.no_target))
            a.update(status=res.status, match=res.match, graph_id=res.graph_id, tie=bool(res.tie), from_release=res.from_release,
                     fallback=res.from_release_fallback, final_database=res.final_database)
            if explain:
                a["path"] = json.dumps(self._paths_convert(res, v), default=str)
            out.append(a)
        return out

    # ---------------- explained paths, in one format for both generations ----------------

    def _edge_dbs_d4c2dda(self, edges):
        G = self.graph
        dbs = []
        for e in edges:
            a, b = e[0], e[1]
            rel = e[3] if len(e) > 3 else None
            if a is None or b is None:
                dbs.append(None)
                continue
            ens_a, ens_b = bool(C.ENS_ID.match(str(a))), bool(C.ENS_ID.match(str(b)))
            ext, ens = (b, a) if ens_a and not ens_b else ((a, b) if ens_b and not ens_a else (None, None))
            if ext is None:
                dbs.append(None)
                continue
            found = set()
            try:
                for t, c in G.xref_edges(ext):
                    if t == ens:
                        for d, asm, rs in c.items:
                            if rel is None or rel in rs:
                                found.add(f"{d}@{asm}")
            except Exception as ex:  # noqa: BLE001
                found.add(f"?{type(ex).__name__}")
            dbs.append(sorted(found))
        return dbs

    def _paths_d4c2dda(self, the_path, keep=40):
        out, seen = [], {}
        for p in the_path:
            k = str(p.get("ensembl_gene"))
            if seen.get(k, 0) >= 2:
                continue
            seen[k] = seen.get(k, 0) + 1
            q = dict(ensembl_gene=p.get("ensembl_gene"), target_id=p.get("target_id"), entry_gene=p.get("entry_gene"),
                     external_jump=p.get("external_jump"), external_gap=p.get("external_gap"), score=p.get("score"),
                     edges=p.get("edges") or [])
            q["edge_dbs"] = self._edge_dbs_d4c2dda(q["edges"])
            out.append(q)
            if len(out) >= keep:
                break
        return out

    def _paths_convert(self, res, view, keep=40):
        steps = res.path.steps if res.path is not None else ()
        out, seen = [], {}
        for c in res.view(view):
            if seen.get(c.gene, 0) >= 2:
                continue
            seen[c.gene] = seen.get(c.gene, 0) + 1
            idx = list(c.path) + list(c.final_path)
            edges, dbs = [], []
            for i in idx:
                s = steps[i]
                edges.append([s.source, s.target, s.kind, s.release if s.release is not None else s.new_release])
                dbs.append([f"{d}@{s.assembly}" for d in s.databases] if s.databases else None)
            out.append(dict(ensembl_gene=c.gene, target_id=c.target, entry_gene=c.entry_gene, external_jump=c.bridges,
                            external_gap=c.bridge_gap, score=c.score, edges=edges, edge_dbs=dbs))
            if len(out) >= keep:
                break
        return out
