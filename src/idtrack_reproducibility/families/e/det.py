"""E7 determinism and concurrency: one set of lines answered under many configurations, compared byte for byte.

Usage (the tree's Python; one configuration per process, its hash seed set before the interpreter starts):
  det.py run --organism <o> --clone <clone> --graph <graph> --out <run>/<organism> --config <name> [--sets ...]
  det.py submit --organism <o> --clone <clone> --graph <graph> --out <run>/<organism> [--configs a,b]
  det.py compare --out <run>/<organism>          -> e7_report.tsv, e7_differences.jsonl.gz

The lines: every line of the harness's published sets (`lines.load("harness", ...)`, judge revised) and family E's
`convert` lines of sets e1, e2, e8 (strings only), each asked with its own arguments.

The canonical answer of a line (docs/methods.md, E, E7): from 71690ab on `ConversionResult.to_dict()` of one
`API.convert` call (both views are in it), on d4c2dda the result dicts of `convert_identifier` for the line's
strategy, written as `json.dumps(sort_keys=True)`; an error as its class and message, except the budget error,
whose message holds seconds (`SearchBudgetExceededError` alone). Searches are bounded by states, not seconds
(`budget_seconds` 10^6 on 71690ab; d4c2dda's own budget of states), so that a budget cannot make an answer depend on
timing. A configuration stores, per qid, the sha256 of the canonical answer; the reference (`ref`: hash seed 0,
one process, file order, single calls) also stores the canonical answers, so that every difference is shown.

Configurations (each compared with `ref` line by line):
  seed{0,1,2,R}_p{1,2,8}   PYTHONHASHSEED 0, 1, 2 or a random seed (recorded); 1, 2 or 8 processes forked after one
                           load of the graph (rule `process_workers`'s processes over one graph; the graph as loaded, not
                           memory-mapped at 71690ab), lines dealt by crc32(qid) % n
  threads8                 8 threads sharing one API
  shuffled                 one process, the lines in a seeded random order
  batch1, batch100, batch10000   `convert_many` (d4c2dda: `convert_identifier_multiple`) over lines of equal arguments
                           (a blank identifier is asked alone: a batch that holds one is refused whole, as documented)
  duplicates               batches of 100 in which every identifier comes twice, the second copy later
  junk                     batches of 100 with 25 identifiers that are nothing interleaved (E8: an added item changes
                           no other)
  after_raise              each call preceded by a call that raises (an empty string; an unknown release; on 71690ab a
                           budget of one state)
  explain                  `explain=True`: the answer without its `path` field
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
import time
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

BIG = 10 ** 6
CONFIGS_CONVERT = ["seed0_p1", "seed1_p1", "seed2_p1", "seedR_p1", "seed0_p2", "seed1_p2", "seed2_p2", "seedR_p2",
                   "seed0_p8", "seed1_p8", "seed2_p8", "seedR_p8", "threads8", "shuffled", "batch1", "batch100",
                   "batch10000", "duplicates", "junk", "after_raise", "explain"]
JUNK = ["zz_nothing_%d_zz" % i for i in range(25)]


def the_lines(organism, sets=None):
    out = []
    hl, _, _ = L.load("harness", organism)
    for q in hl:
        out.append(dict(qid="harness|" + q["qid"], id=q["id"], from_release=q.get("from_release"), to_release=q.get("to_release"),
                        final_database=q.get("final_database"), strategy=q.get("strategy") or "best"))
    names = set(L.set_names("e", organism))
    for s in (sets or ["e1", "e2", "e8"]):
        if s not in names:
            continue
        el, _, _ = L.load("e", organism, [s])
        for q in el:
            if (q.get("call") or "convert") != "convert" or q.get("kwargs") or (q.get("id_codec") or "str") != "str" \
                    or not isinstance(q.get("id"), str) or len(q["id"]) > 2000:
                continue
            out.append(dict(qid=q["qid"], id=q["id"], from_release=q.get("from_release"), to_release=q.get("to_release"),
                            final_database=q.get("final_database"), from_database=q.get("from_database"),
                            strategy=q.get("strategy") or "best"))
    return out


class Asker:
    def __init__(self, clone, graph, organism, workdir):
        sys.path.insert(0, clone)
        from idtrack_reproducibility.families import adapter as F0
        self.gen = C.generation(clone)
        import idtrack
        link = F0.link_graph(graph, workdir)
        if self.gen == "d4c2dda":
            self.api = idtrack.API(local_repository=os.path.dirname(link))
            self.api.build_graph(organism_name=organism, snapshot_release=C.release_of_graph(graph))
        else:
            self.api = idtrack.API.load(link)
        import gc
        gc.collect()
        gc.freeze()

    def kw(self, q, explain=False):
        k = {x: q.get(x) for x in ("from_release", "to_release", "final_database", "from_database") if q.get(x) is not None}
        if self.gen == "d4c2dda":
            k["strategy"] = q["strategy"] if q["strategy"] in ("best", "all") else "best"
        else:
            k["budget_seconds"] = BIG
        if explain:
            k["explain"] = True
        return k

    @staticmethod
    def err(e):
        n = type(e).__name__
        return {"raised": n} if n == "SearchBudgetExceededError" else {"raised": n, "message": str(e)[:2000]}

    def plain(self, r, drop_path=False):
        """The canonical answer: the result without its explained path (None unless asked for; `explain` compares the
        rest)."""
        d = dict(r) if self.gen == "d4c2dda" else r.to_dict()
        d.pop("the_path" if self.gen == "d4c2dda" else "path", None)
        for c in d.get("candidates") or ():  # the indices into the explained path (filled only with explain)
            c.pop("path", None)
            c.pop("final_path", None)
        return d

    def single(self, q, explain=False):
        try:
            if self.gen == "d4c2dda":
                return self.plain(self.api.convert_identifier(q["id"], **self.kw(q, explain)), explain)
            return self.plain(self.api.convert(q["id"], **self.kw(q, explain)), explain)
        except Exception as e:  # noqa: BLE001
            return self.err(e)

    def batch(self, qs, ids):
        """One batch call over the identifiers `ids` with the arguments of qs[0]; results in order (or one error)."""
        k = self.kw(qs[0])
        try:
            if self.gen == "d4c2dda":
                res = self.api.convert_identifier_multiple(list(ids), verbose=False, **k)
            else:
                res = self.api.convert_many(list(ids), **k)
        except Exception as e:  # noqa: BLE001
            return [self.err(e)] * len(ids)
        out = []
        for r in res:
            d = self.plain(r)
            # an item whose conversion raised: written as the single call's error would be
            e = d.get("error") if isinstance(d, dict) else None
            if e and (self.gen == "d4c2dda" or d.get("status") == "error"):
                out.append({"raised": e["type"]} if e["type"] == "SearchBudgetExceededError" else
                           {"raised": e["type"], "message": e["message"][:2000]})
            else:
                out.append(d)
        return out

    def raiser(self, i):
        try:
            if i % 3 == 0:
                self.api.convert("") if self.gen != "d4c2dda" else self.api.convert_identifier("")
            elif i % 3 == 1:
                self.api.convert("TP53", to_release=10 ** 9) if self.gen != "d4c2dda" else \
                    self.api.convert_identifier("TP53", to_release=10 ** 9)
            elif self.gen != "d4c2dda":
                self.api.convert("ENSG00000141510", budget_states=1, budget_seconds=BIG)
            else:
                self.api.convert_identifier("")
        except Exception:  # noqa: BLE001
            pass


def canon(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"), default=str)


def dig(x):
    return hashlib.sha256(canon(x).encode()).hexdigest()[:32]


def group_key(q):
    return canon([q.get("from_release"), q.get("to_release"), q.get("final_database"), q.get("from_database"), q["strategy"]])


def run_config(a):
    cfg = a.config
    qs = the_lines(a.organism, a.sets)
    qs.sort(key=lambda q: q["qid"])
    if a.limit:
        qs = [q for q in qs if zlib.crc32(q["qid"].encode()) % a.limit == 0]  # smoke runs: one line in `limit`
    outdir = os.path.join(a.out, "e7", cfg)
    os.makedirs(outdir, exist_ok=True)
    asker = Asker(a.clone, a.graph, a.organism, a.out)
    t0 = time.time()
    procs = 1
    if cfg.startswith("seed"):
        procs = int(cfg.split("_p")[1])
    results = {}
    if cfg == "ref" or (cfg.startswith("seed") and procs == 1):
        for q in qs:
            results[q["qid"]] = asker.single(q)
    elif cfg.startswith("seed"):
        r_, pids = [], []
        for w in range(procs):
            rd, wr = os.pipe()
            pid = os.fork()
            if pid == 0:  # the child: its share of the lines, its digests to the pipe
                os.close(rd)
                part = {q["qid"]: dig(asker.single(q)) for q in qs if zlib.crc32(q["qid"].encode()) % procs == w}
                data = json.dumps(part).encode()
                with os.fdopen(wr, "wb") as fo:
                    fo.write(data)
                os._exit(0)
            os.close(wr)
            r_.append(rd)
            pids.append(pid)
        for rd in r_:
            with os.fdopen(rd, "rb") as fi:
                results.update({k: ("digest", v) for k, v in json.loads(fi.read().decode()).items()})
        for pid in pids:
            os.waitpid(pid, 0)
    elif cfg == "threads8":
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(8) as pool:
            for q, r in zip(qs, pool.map(asker.single, qs)):
                results[q["qid"]] = r
    elif cfg == "shuffled":
        order = qs[:]
        random.Random(20261118).shuffle(order)
        for q in order:
            results[q["qid"]] = asker.single(q)
    elif cfg == "explain":
        for q in qs:
            results[q["qid"]] = asker.single(q, explain=True)
    elif cfg == "after_raise":
        for i, q in enumerate(qs):
            asker.raiser(i)
            results[q["qid"]] = asker.single(q)
    else:
        groups = {}
        for q in qs:
            if not q["id"].strip():
                # a batch holding an empty identifier is refused as a whole (documented): these lines are asked alone
                results[q["qid"]] = asker.single(q)
                continue
            groups.setdefault(group_key(q), []).append(q)
        size = {"batch1": 1, "batch100": 100, "batch10000": 10000, "duplicates": 100, "junk": 100}[cfg]
        rng = random.Random(20261118)
        for g in groups.values():
            for i in range(0, len(g), size):
                chunk = g[i:i + size]
                ids = [q["id"] for q in chunk]
                if cfg == "duplicates":
                    ids2 = ids + ids[::-1]
                    res = asker.batch(chunk, ids2)
                    n = len(ids)
                    for j, q in enumerate(chunk):
                        a1, a2 = res[j], res[2 * n - 1 - j]
                        results[q["qid"]] = a1 if canon(a1) == canon(a2) else {"duplicate_differs": [a1, a2]}
                elif cfg == "junk":
                    mixed, where = [], []
                    for j, x in enumerate(ids):
                        if j % 4 == 0:
                            mixed.append(JUNK[rng.randrange(len(JUNK))])
                        where.append(len(mixed))
                        mixed.append(x)
                    res = asker.batch(chunk, mixed)
                    for j, q in enumerate(chunk):
                        results[q["qid"]] = res[where[j]]
                else:
                    res = asker.batch(chunk, ids)
                    for q, r in zip(chunk, res):
                        results[q["qid"]] = r
    secs = time.time() - t0
    import pyarrow as pa
    import pyarrow.parquet as pq
    rows = []
    for qid, r in results.items():
        if isinstance(r, tuple) and r[0] == "digest":
            rows.append(dict(qid=qid, digest=r[1], answer=None))
        else:
            rows.append(dict(qid=qid, digest=dig(r), answer=canon(r) if cfg == "ref" else None))
    pq.write_table(pa.Table.from_pylist(rows), os.path.join(outdir, "digests.parquet.tmp"), compression="zstd")
    os.replace(os.path.join(outdir, "digests.parquet.tmp"), os.path.join(outdir, "digests.parquet"))
    if cfg != "ref":
        # the answers of the lines that differ from the reference (asked again here, as this configuration gave them)
        ref = os.path.join(a.out, "e7", "ref", "digests.parquet")
        if os.path.exists(ref):
            rd = dict(zip(*[pq.read_table(ref, columns=["qid", "digest"]).column(c).to_pylist() for c in ("qid", "digest")]))
            diff = [x for x in rows if rd.get(x["qid"]) not in (None, x["digest"])]
            with gzip.open(os.path.join(outdir, "differences.jsonl.gz"), "wt") as fo:
                for x in diff[:5000]:
                    r = results[x["qid"]]
                    fo.write(canon(dict(qid=x["qid"], answer=None if isinstance(r, tuple) else r)) + "\n")
    json.dump(dict(config=cfg, organism=a.organism, n=len(rows), seconds=round(secs, 1), generation=asker.gen,
                   pythonhashseed=os.environ.get("PYTHONHASHSEED"), host=os.uname().nodename, slurm_job=os.environ.get("SLURM_JOB_ID"),
                   clone=a.clone, commit=C.git_commit(a.clone)[0], graph=os.path.realpath(a.graph), finished=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(outdir, "info.json"), "w"), indent=1)
    sys.stderr.write(f"{cfg}: {len(rows)} lines in {secs:.0f} s\n")


def seed_of(cfg):
    if cfg in ("ref",) or not cfg.startswith("seed"):
        return "0"
    s = cfg[4:].split("_")[0]
    return str(random.SystemRandom().randrange(3, 2 ** 32 - 1)) if s == "R" else s


CONFIGS_D4C2DDA = ["seed0_p1", "seed1_p1", "seed2_p1", "seedR_p1", "seed0_p2", "seed1_p2", "seed2_p2", "seedR_p2", "threads8",
                   "shuffled", "batch1", "batch100", "batch10000", "duplicates", "junk", "after_raise", "explain"]


def submit(a):
    """One job per organism and tree that runs the reference, then every configuration, one process after the other
    (each with its own hash seed); configurations already done are skipped, so a resubmission continues.
    d4c2dda's graphs take 11 to 50 GB per process: there 1 and 2 processes only (8 forked processes of the human
    graph would need about 400 GB)."""
    gen = C.generation(a.clone)
    cfgs = a.configs.split(",") if a.configs else (CONFIGS_CONVERT if gen == "convert" else CONFIGS_D4C2DDA)
    logs = os.path.join(a.out, "e7", "logs")
    os.makedirs(logs, exist_ok=True)
    from idtrack_reproducibility import config as CFG
    from idtrack_reproducibility import jobs
    py = C.python_for(a.clone)
    me = "-m idtrack_reproducibility.families.e.det"
    common = f"--organism {a.organism} --clone {a.clone} --graph {a.graph} --out {a.out}" + (f" --sets {a.sets_arg}" if a.sets_arg else "")
    if gen == "d4c2dda":
        mem = 2 * jobs.memory_for(a.organism) + 8
        cpus = 8
    else:
        mem = 48
        cpus = 9
    script = os.path.join(logs, f"bundle_{time.strftime('%Y%m%d_%H%M%S')}.sh")
    with open(script, "w") as fo:
        fo.write("#!/bin/bash\n")
        for cfg in ["ref"] + [c for c in cfgs if c != "ref"]:
            done = os.path.join(a.out, "e7", cfg, "digests.parquet")
            seed = seed_of(cfg)
            fo.write(f"[ -f {done} ] || PYTHONHASHSEED={seed} {py} {me} run {common} --config {cfg}" + (" || exit 1" if cfg == "ref" else "") + "\n")
        fo.write(f"{CFG.get().python} {me} compare --out {a.out}\n")
    os.chmod(script, 0o700)
    jid = jobs.submit(f"e7_{gen}_{a.organism}", f"bash {script}", logs, cpus=cpus, mem_gb=mem, time=a.time)
    print(jid)
    json.dump(dict(job=jid, script=script, configs=cfgs), open(script + ".json", "w"))


def compare(a):
    import pyarrow.parquet as pq
    root = os.path.join(a.out, "e7")
    ref = pq.read_table(os.path.join(root, "ref", "digests.parquet")).to_pandas().set_index("qid")
    rows = []
    from idtrack_reproducibility.families import report as RP
    with gzip.open(os.path.join(a.out, "e7_differences.jsonl.gz"), "wt") as fo:
        for cfg in sorted(os.listdir(root)):
            p = os.path.join(root, cfg, "digests.parquet")
            if cfg in ("ref", "logs") or not os.path.exists(p):
                continue
            info = json.load(open(os.path.join(root, cfg, "info.json")))
            d = pq.read_table(p, columns=["qid", "digest"]).to_pandas().set_index("qid")
            j = d.join(ref[["digest"]].rename(columns={"digest": "digest_ref"}), how="outer")
            n = int(j["digest_ref"].notna().sum())
            missing = int(j["digest"].isna().sum())
            same = int((j["digest"] == j["digest_ref"]).sum())
            differ = n - same - missing
            lo, hi = RP.wilson(same, same + differ) if same + differ else (float("nan"), float("nan"))
            kinds = {}
            dfile = os.path.join(root, cfg, "differences.jsonl.gz")
            if differ and os.path.exists(dfile):
                for raw in gzip.open(dfile, "rt"):
                    x = json.loads(raw)
                    k = x["qid"].split("|")[1] if x["qid"].startswith("harness|") else x["qid"].split("|")[1]
                    kinds[k] = kinds.get(k, 0) + 1
                    fo.write(canon(dict(config=cfg, qid=x["qid"], got=x["answer"], ref=json.loads(ref.at[x["qid"], "answer"])
                                        if isinstance(ref.at[x["qid"], "answer"], str) else None)) + "\n")
            rows.append([cfg, info.get("generation"), info.get("pythonhashseed"), n, same, differ, missing,
                         f"{same / (same + differ):.6f}" if same + differ else "", f"{lo:.6f}", f"{hi:.6f}", info.get("seconds"),
                         ";".join(f"{k}:{v}" for k, v in sorted(kinds.items(), key=lambda kv: -kv[1])[:10])])
    hdr = ["config", "generation", "pythonhashseed", "n", "identical", "different", "missing", "share_identical", "wilson_low",
           "wilson_high", "seconds", "different_by_kind"]
    with open(os.path.join(a.out, "e7_report.tsv"), "w") as fo:
        fo.write("\t".join(hdr) + "\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    print("\n".join("\t".join(str(x) for x in r) for r in [hdr] + rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["run", "submit", "compare"])
    ap.add_argument("--organism")
    ap.add_argument("--clone")
    ap.add_argument("--graph")
    ap.add_argument("--out", required=True)
    ap.add_argument("--config")
    ap.add_argument("--configs")
    ap.add_argument("--sets", default=None)
    ap.add_argument("--time", default="24:00:00")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    a.out = os.path.abspath(a.out)
    a.sets_arg = a.sets
    a.sets = a.sets.split(",") if a.sets else None
    if a.clone:
        a.clone = os.path.abspath(a.clone)
    if a.graph:
        a.graph = os.path.abspath(a.graph)
    if a.what == "run":
        run_config(a)
    elif a.what == "submit":
        submit(a)
    else:
        compare(a)


if __name__ == "__main__":
    main()
