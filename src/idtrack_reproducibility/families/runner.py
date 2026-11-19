"""The runner (F0): ask a family's query lines through one tree and one graph, sharded and resumable.

Usage:
  runner.py --family <f> --organism <o> --clone <clone> --graph <graph file> --out <run folder>/<organism>
            [--sets a,b] [--strata x,y] [--kinds k] [--timeout 30] [--release 115] [--shard i/n] [--limit N]
  runner.py ... --submit 16 [--score [--score-parts P]]   16 shard jobs (and, with --score, a scoring job after them;
                                        with --score-parts, P partition jobs and a merge job), then exit; then
  runner.py ... --submit 0 --score --score-parts P [--scores-to <folder>] [--judge published]
                                        a rescoring of the answers already in <out> (scoring jobs only)

- Jobs go through the configured runner (`idtrack_reproducibility.jobs`): names `<prefix><family>_<run>_<organism>_s<i>of<n>`,
  `..._score`, `..._score<p>of<P>`, `..._merge<P>` (`--tag` replaces `<family>_<run>_<organism>`), one name per work
  unit, kept by a resubmission; every step resumes (a shard asks only what has no answer, a scoring is whole).

- A line belongs to shard `crc32(qid) % n` (as the harness's score.py); a shard skips every qid that has an answer in
  any answer file of the folder, whatever its sharding, so a resubmission (any n, or after a new set is published)
  asks only what is missing.
- The graph is loaded once per shard (adapter.Adapter: the file linked into `<out>/repo/`), `gc.freeze()` after;
  `PYTHONHASHSEED=0` (the runner restarts itself with it when it is not set), a per-call budget (`--timeout`).
- Answers (adapter.py: one row per line and view) go to `<out>/answers/shard_<i>_of_<n>.<k>.parquet`, a new part
  every `--flush` lines or 5 minutes (an interrupted shard loses at most one part; those lines are asked again).
- The tree's Python is the configuration's (`[trees.<name>] python`, common.python_for); the scoring jobs run in this
  package's Python. Jobs go through the configured runner (jobs.submit).
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time

from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L


def answer_schema():
    import pyarrow as pa
    return pa.schema([
        ("qid", pa.string()), ("h", pa.uint32()), ("view", pa.string()), ("outcome", pa.string()),
        ("got", pa.list_(pa.string())), ("apart", pa.list_(pa.string())), ("status", pa.string()),
        ("match", pa.string()), ("tie", pa.bool_()), ("from_release", pa.int64()), ("fallback", pa.bool_()),
        ("final_database", pa.string()), ("flags", pa.string()), ("error_type", pa.string()),
        ("error_mro", pa.list_(pa.string())), ("error", pa.string()), ("traceback", pa.string()),
        ("seconds", pa.float64()), ("calls", pa.int32()), ("path", pa.string())])


def answered(out, shard=None):
    """qids with an answer in any answer file of the folder (with shard=(i, n): only those of that shard, read
    through the answers' `h` column, so a shard of a set of hundreds of millions of lines holds only its own)."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    done = set()
    for p in glob.glob(os.path.join(out, "answers", "*.parquet")):
        t = pq.read_table(p, columns=["qid", "h"])
        if shard is not None:
            t = t.filter(pa.array(t.column("h").to_numpy() % shard[1] == shard[0]))
        done.update(t.column("qid").to_pylist())
    return done


def write_part(out, i, n, rows):
    import pyarrow as pa
    import pyarrow.parquet as pq
    d = os.path.join(out, "answers")
    k = len(glob.glob(os.path.join(d, f"shard_{i}_of_{n}.*.parquet")))
    while os.path.exists(os.path.join(d, f"shard_{i}_of_{n}.{k}.parquet")):
        k += 1
    for r in rows:
        r["h"] = C.crc(r["qid"])
        if isinstance(r.get("flags"), dict):
            r["flags"] = C.dumps(r["flags"])
        if r.get("from_release") is not None:
            r["from_release"] = int(r["from_release"])
    path = os.path.join(d, f"shard_{i}_of_{n}.{k}.parquet")
    pq.write_table(pa.Table.from_pylist(rows, schema=answer_schema()), path + ".tmp", compression="zstd")
    os.replace(path + ".tmp", path)


def run(a):
    i, n = (int(x) for x in a.shard.split("/"))
    mine, shas, excluded = L.load(a.family, a.organism, a.sets, a.strata, a.kinds, shard=(i, n), limit=a.limit)
    os.makedirs(os.path.join(a.out, "answers"), exist_ok=True)
    done = answered(a.out, (i, n))
    todo = [q for q in mine if q["qid"] not in done]
    sys.stderr.write(f"shard {i}/{n}: {len(mine)} lines, {len(todo)} to ask\n")
    if not todo:
        return
    from idtrack_reproducibility.families import adapter
    ad = adapter.Adapter(a.clone, a.graph, a.organism, a.release, a.out, timeout=a.timeout)
    info = dict(ad.info, shard=a.shard, sets=shas, host=os.uname().nodename, slurm_job=os.environ.get("SLURM_JOB_ID"),
                started=time.strftime("%Y-%m-%d %H:%M:%S"), n_lines=len(mine), n_todo=len(todo))
    json.dump(info, open(os.path.join(a.out, "answers", f"load_{i}_of_{n}.json"), "w"), indent=1)
    buf, t_all, t_flush = [], time.time(), time.time()
    for k, q in enumerate(todo):
        buf.extend(ad.ask(q))
        if len(buf) >= a.flush or time.time() - t_flush > 300:
            write_part(a.out, i, n, buf)
            buf, t_flush = [], time.time()
        if k % 5000 == 0:
            sys.stderr.write(f"{k}/{len(todo)} {time.time() - t_all:.0f}s\n")
    if buf:
        write_part(a.out, i, n, buf)
    sys.stderr.write(f"shard {i}/{n} done in {time.time() - t_all:.0f}s\n")


def submit(a):
    """The shards, and with --score the scoring after them (one job, or --score-parts partitions and a merge), through
    the configured runner (jobs.submit; job names `<prefix><tag>_<step>`, a resubmission of a step keeps its name; a
    scoring job's dependency is `afterany`: the scorer reports lines without answers as missing). `--submit 0 --score`
    is a rescoring of the answers already there, into `--scores-to`. Writes `<logs>/submit_<time>.json`."""
    from idtrack_reproducibility import config as CFG
    from idtrack_reproducibility import jobs
    logs = os.path.join(a.out if a.submit else (a.scores_to or a.out), "logs")  # a rescoring logs beside its scores
    os.makedirs(logs, exist_ok=True)
    sel = ((f" --sets {','.join(a.sets)}" if a.sets else "") + (f" --strata {','.join(a.strata)}" if a.strata else "")
           + (f" --kinds {','.join(a.kinds)}" if a.kinds else "") + (f" --limit {a.limit}" if a.limit else ""))
    common = (f"--family {a.family} --organism {a.organism} --clone {a.clone} --graph {a.graph} --out {a.out} "
              f"--timeout {a.timeout} --release {a.release} --flush {a.flush}" + sel)
    py = C.python_for(a.clone)
    mem = a.mem or jobs.memory_for(a.organism)
    tag = a.tag or f"{a.family}_{os.path.basename(os.path.dirname(a.out))}_{a.organism}"
    ids = [jobs.submit(f"{tag}_s{i}of{a.submit}", jobs.module(py, "idtrack_reproducibility.families.runner",
                                                              f"{common} --shard {i}/{a.submit}"),
                       logs, cpus=a.cpus, mem_gb=mem, time=a.time) for i in range(a.submit)]
    print("shards", " ".join(ids))
    if a.score:
        scargs = (f"--family {a.family} --organism {a.organism} --run {a.out} --judge {a.judge}" + sel
                  + (f" --scores-to {a.scores_to}" if a.scores_to else ""))

        def job(name, args, dep):
            return jobs.submit(f"{tag}_{name}", jobs.module(CFG.get().python, "idtrack_reproducibility.families.scorer",
                                                           f"{scargs} {args}"),
                               logs, cpus=2, mem_gb=a.score_mem, time=a.score_time, after=dep)
        if a.score_parts > 1:
            parts = [job(f"score{p}of{a.score_parts}", f"--parts {a.score_parts} --part {p}", ids)
                     for p in range(a.score_parts)]
            m = job(f"merge{a.score_parts}", f"--parts {a.score_parts} --merge", parts)
            print("score parts", " ".join(parts), "merge", m)
            ids += parts + [m]
        else:
            m = job("score", "", ids)
            print("score", m)
            ids.append(m)
    json.dump(dict(jobs=ids, submitted=time.strftime("%Y-%m-%d %H:%M:%S"), args=vars(a)),
              open(os.path.join(logs, f"submit_{time.strftime('%Y%m%d_%H%M%S')}.json"), "w"), indent=1)
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--clone", required=True)
    ap.add_argument("--graph", required=True, help="the graph file (linked into <out>/repo/)")
    ap.add_argument("--out", required=True, help="the run folder of this organism: <families>/<family>/runs/<run>/<organism>")
    ap.add_argument("--sets", default=None)
    ap.add_argument("--strata", default=None)
    ap.add_argument("--kinds", default=None)
    ap.add_argument("--limit", type=int, default=None, help="smoke run: the first N lines of each stratum")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--release", type=int, default=None)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--flush", type=int, default=20000)
    ap.add_argument("--submit", type=int, default=0)
    ap.add_argument("--score", action="store_true", help="with --submit: a scoring job after the shards")
    ap.add_argument("--score-mem", type=int, default=32)
    ap.add_argument("--score-parts", type=int, default=1, help="with --score: partitions scored by parallel jobs, then a merge")
    ap.add_argument("--scores-to", default=None, help="with --score: the scored folder (default: the run folder)")
    ap.add_argument("--judge", default="revised", choices=["revised", "published"])
    ap.add_argument("--cpus", type=int, default=2)
    ap.add_argument("--mem", type=int, default=None)
    ap.add_argument("--time", default="12:00:00")
    ap.add_argument("--score-time", default="12:00:00")
    ap.add_argument("--tag", default=None, help="the run's part of the job names (default <family>_<run>_<organism>)")
    a = ap.parse_args()
    if os.environ.get("PYTHONHASHSEED") != "0" and not (a.submit or a.score):
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, "-m", __spec__.name] + sys.argv[1:] if __spec__ else [sys.executable] + sys.argv)
    a.out, a.clone, a.graph = os.path.abspath(a.out), os.path.abspath(a.clone), os.path.abspath(a.graph)
    a.scores_to = os.path.abspath(a.scores_to) if a.scores_to else None
    for k in ("sets", "strata", "kinds"):
        setattr(a, k, getattr(a, k).split(",") if getattr(a, k) else None)
    a.release = a.release or C.release_of_graph(a.graph)
    os.makedirs(a.out, exist_ok=True)
    if a.submit or a.score:  # --submit 0 --score: a rescoring of answers already there (the scoring jobs only)
        return submit(a)
    run(a)


if __name__ == "__main__":
    main()
