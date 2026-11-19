"""The regression corpus's steps as subcommands (docs/corpus.md, "Making and comparing").

    python -m idtrack_reproducibility.corpus <subcommand> --tree <tree> --organism <organism> [options] [--submit]

- `make`: draw a corpus with a tree of the d4c2dda generation (`make_corpus.py`) into a new file.
- `run`: ask the corpus through a tree (`run_corpus.py`) into `<corpus>/runs/<tree>/<organism>/<tag>.jsonl.gz`.
- `harmonize`: the HarmonizeFeatures part (`run_harmonize.py`, one process per configuration, the parts joined)
  into `<tag>.harmonize.jsonl.gz` of the run folder.
- `compare`: a run (or a reference) against another tree's reference: `compare_<tag>_vs_<other>.{txt,json}`.
- `reference`: make a tree's reference answers: the corpus twice in two processes under the hash seeds 0 and 1
  (`run1`, `run2`), their comparison (`determinism.compare.{txt,json}`), the calls over the cap of `run1` asked again
  with a long cap and merged (`timeouts.txt`, `slow.jsonl.gz`, `merged.jsonl.gz`), the harmonize part
  (`harmonize.jsonl.gz`), then both copied to `<corpus>/reference/<organism>.<tree>.jsonl.gz` and
  `.harmonize.jsonl.gz`, and a summary `reference.json`; with `--compare-with <tree>`, both compared with that tree's
  reference. It refuses a checkout with local changes or at another commit than the configured one (a reference is
  named by its commit), and never overwrites a reference. `--no-network` leaves out the calls that need Ensembl over
  the network (`resolve_organism`), for a machine or a day without it; `reference.json` says so.
- `summarize`: answers and exceptions of a run by entry and stratum (`summarize_results.py`).

`--submit` runs the same subcommand as one job through the configured runner (`jobs.submit`, named
`corpus_<subcommand>_<tree>_<organism>`, log in the run folder's `logs/`), memory from `[jobs.memory_gb]`. Without it
the subcommand runs here; the steps that ask idtrack run in the tree's Python (`[trees.<tree>] python`) in their own
processes, with `PYTHONHASHSEED` set per step.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time

from idtrack_reproducibility import config as CFG
from idtrack_reproducibility import jobs
from idtrack_reproducibility.corpus import _common as C
from idtrack_reproducibility.corpus import compare as CMP


# ---------------------------------------------------------------------------------------------------------------------
# Helpers that read answer files (no idtrack).


def timed_out(path: str) -> list[str]:
    """The qids of the calls that hit the cap in an answer file."""
    out = []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            r = json.loads(ln)
            if "timeout" in r:
                out.append(r["qid"])
    return out


def merge(base: str, rerun: str, out: str) -> int:
    """Write `base` with the records of `rerun` in place of the same qids, unless the rerun also hit a cap not larger
    than the base's. The meta line of the base is kept with the rerun's meta as `meta_rerun`. Returns the number of
    records replaced."""
    re_, re_meta = {}, {}
    with gzip.open(rerun, "rt") as fh:
        for ln in fh:
            r = json.loads(ln)
            if "meta" in r or "meta_end" in r:
                re_meta.update(r)
                continue
            re_[r["qid"]] = r
    n_rep = 0
    with gzip.open(base, "rt") as fh, gzip.open(out + ".partial", "wt") as fo:
        for ln in fh:
            r = json.loads(ln)
            if "meta" in r:
                fo.write(json.dumps({**r, "meta_rerun": re_meta}) + "\n")
                continue
            if "qid" in r and r["qid"] in re_:
                n = re_[r["qid"]]
                if not ("timeout" in n and n["timeout"] <= r.get("timeout", 0)):
                    r = n
                    n_rep += 1
            fo.write(json.dumps(r, sort_keys=True) + "\n")
    os.replace(out + ".partial", out)
    return n_rep


def join_gzip(parts: list[str], out: str) -> None:
    """Concatenate gzip files (gzip members concatenate into one valid stream)."""
    with open(out + ".partial", "wb") as fo:
        for p in parts:
            with open(p, "rb") as fi:
                shutil.copyfileobj(fi, fo)
    os.replace(out + ".partial", out)


# ---------------------------------------------------------------------------------------------------------------------
# Steps that ask a tree (each in its own process of the tree's Python).


def tree_process(tree: str, module: str, args: list[str], hash_seed: int, log: str) -> None:
    t = C.tree(tree)
    env = dict(os.environ, **jobs.job_env({"PYTHONHASHSEED": str(hash_seed)}))
    env["PYTHONPATH"] = CFG.SRC + (os.pathsep + os.environ["PYTHONPATH"] if os.environ.get("PYTHONPATH") else "")
    cmd = [t["python"], "-m", module] + args
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "a") as fo:
        fo.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} PYTHONHASHSEED={hash_seed} {' '.join(cmd)}\n")
        fo.flush()
        r = subprocess.run(cmd, stdout=fo, stderr=subprocess.STDOUT, env=env)
    if r.returncode:
        raise SystemExit(f"{module} failed (exit {r.returncode}); log: {log}")


def run_file(tree: str, organism: str, tag: str) -> str:
    return os.path.join(C.run_dir(tree, organism), f"{tag}.jsonl.gz")


def do_run(tree, organism, tag="run1", hash_seed=0, timeout=60, extra=(), corpus=None):
    out = run_file(tree, organism, tag)
    args = ["--tree", tree, "--organism", organism, "--out", out, "--timeout", str(timeout)] + list(extra)
    if corpus:
        args += ["--corpus", corpus]
    tree_process(tree, "idtrack_reproducibility.corpus.run_corpus", args, hash_seed,
                 os.path.join(C.run_dir(tree, organism), "logs", f"{tag}.log"))
    return out


def do_harmonize(tree, organism, tag="harmonize", hash_seed=0, corpus=None, exclude_slow=None):
    from idtrack_reproducibility.corpus.run_harmonize import CONFIGURATIONS

    d = C.run_dir(tree, organism)
    out = os.path.join(d, f"{tag}.jsonl.gz")
    parts = []
    for cfg in CONFIGURATIONS:
        part = os.path.join(d, f"{tag}.{cfg}.jsonl.gz")
        args = ["--tree", tree, "--organism", organism, "--out", part, "--only", cfg,
                "--workdir", os.path.join(d, "harmonize")]
        if corpus:
            args += ["--corpus", corpus]
        if exclude_slow:
            args += ["--exclude-slow", exclude_slow]
        tree_process(tree, "idtrack_reproducibility.corpus.run_harmonize", args, hash_seed,
                     os.path.join(d, "logs", f"{tag}.log"))
        parts.append(part)
    join_gzip(parts, out)
    for p in parts:
        os.remove(p)
    return out


def write_compare(a, b, corpus, stem):
    rep = CMP.compare(a, b, corpus)
    with open(stem + ".json", "w") as fh:
        json.dump(rep, fh, indent=1)
    txt = CMP.text(rep)
    with open(stem + ".txt", "w") as fh:
        fh.write(txt)
    return txt.splitlines()[0]


def clean_checkout(tree: str) -> str:
    """The checkout's HEAD commit; refuses local changes or a commit other than the configured one."""
    t = C.tree(tree)
    st = subprocess.run(["git", "-C", t["source"], "status", "--porcelain"], capture_output=True, text=True, check=True)
    if st.stdout.strip():
        raise SystemExit(f"tree {tree}: the checkout has local changes; a reference is named by a commit")
    head = subprocess.run(["git", "-C", t["source"], "rev-parse", "HEAD"], capture_output=True, text=True,
                          check=True).stdout.strip()
    if t["commit"] and not head.startswith(t["commit"]) and not t["commit"].startswith(head):
        raise SystemExit(f"tree {tree}: the checkout is at {head}, the configuration names {t['commit']}")
    return head


def do_reference(tree, organism, timeout=60, slow_cap=900, compare_with=None, no_network=False):
    head = clean_checkout(tree)
    ref = C.reference_file(organism, tree)
    href = C.reference_file(organism, tree, harmonize=True)
    for p in (ref, href):
        if os.path.exists(p):
            raise SystemExit(f"{p} exists: a reference is never overwritten")
    d = C.run_dir(tree, organism)
    corpus = C.corpus_file(organism)
    summary = {"tree": tree, "commit": head, "organism": organism, "corpus": os.path.basename(corpus),
               "corpus_sha256": C.sha256(corpus), "timeout": timeout, "slow_cap": slow_cap, "network": not no_network,
               "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    extra = ["--no-network"] if no_network else []
    run1 = do_run(tree, organism, "run1", 0, timeout, extra)
    run2 = do_run(tree, organism, "run2", 1, timeout, extra)
    summary["determinism"] = write_compare(run1, run2, corpus, os.path.join(d, "determinism.compare"))
    slow = timed_out(run1)
    with open(os.path.join(d, "timeouts.txt"), "w") as fh:
        fh.write("".join(q + "\n" for q in slow))
    merged = os.path.join(d, "merged.jsonl.gz")
    summary["timeouts"] = len(slow)
    if slow:
        slow_run = do_run(tree, organism, "slow", 0, slow_cap, extra + ["--qids", os.path.join(d, "timeouts.txt")])
        summary["replaced_by_slow_pass"] = merge(run1, slow_run, merged)
        summary["still_over_cap"] = len(timed_out(merged))
    else:
        shutil.copyfile(run1, merged)
    harm = do_harmonize(tree, organism, "harmonize", 0)
    os.makedirs(os.path.dirname(ref), exist_ok=True)
    shutil.copyfile(merged, ref)
    shutil.copyfile(harm, href)
    summary["reference"] = os.path.relpath(ref, C.root())
    summary["reference_harmonize"] = os.path.relpath(href, C.root())
    if compare_with:
        summary[f"compare_{compare_with}"] = write_compare(C.reference_file(organism, compare_with), ref, corpus,
                                                          os.path.join(d, f"compare_vs_{compare_with}"))
        summary[f"compare_{compare_with}_harmonize"] = write_compare(
            C.reference_file(organism, compare_with, harmonize=True), href, None,
            os.path.join(d, f"compare_vs_{compare_with}.harmonize"))
    summary["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    C.write_json(os.path.join(d, "reference.json"), summary)
    print(json.dumps(summary, indent=1))
    return summary


# ---------------------------------------------------------------------------------------------------------------------


def submit(a, argv):
    """The same subcommand as one job (without --submit)."""
    rest = [x for x in argv if x != "--submit"]
    d = C.run_dir(a.tree, a.organism)
    logs = os.path.join(d, "logs")
    os.makedirs(logs, exist_ok=True)
    cmd = jobs.module(CFG.get().python, "idtrack_reproducibility.corpus", " ".join(_quote(x) for x in rest))
    jid = jobs.submit(f"corpus_{a.cmd}_{a.tree}_{a.organism}", cmd, logs, cpus=a.cpus,
                      mem_gb=a.mem or jobs.memory_for(a.organism), time=a.time)
    print(jid)
    return jid


def _quote(x):
    import shlex
    return shlex.quote(x)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="python -m idtrack_reproducibility.corpus", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--tree", required=True)
        p.add_argument("--organism", required=True)
        p.add_argument("--submit", action="store_true", help="run as one job through the configured runner")
        p.add_argument("--cpus", type=int, default=2)
        p.add_argument("--mem", type=int, default=None, help="GB (default: [jobs.memory_gb] of the organism)")
        p.add_argument("--time", default="08:00:00")

    p = sub.add_parser("make", help="draw a corpus into a new file")
    common(p)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=20261111)
    p.add_argument("--tables", action="store_true")
    p = sub.add_parser("run", help="ask the corpus through a tree")
    common(p)
    p.add_argument("--tag", default="run1")
    p.add_argument("--hash-seed", type=int, default=0)
    p.add_argument("--timeout", type=int, default=60)
    p.add_argument("--corpus", default=None)
    p.add_argument("--no-network", action="store_true")
    p.add_argument("--shuffle", type=int, default=None)
    p.add_argument("--entries", default="")
    p.add_argument("--limit", type=int, default=None)
    p = sub.add_parser("harmonize", help="the HarmonizeFeatures part")
    common(p)
    p.add_argument("--tag", default="harmonize")
    p.add_argument("--hash-seed", type=int, default=0)
    p.add_argument("--corpus", default=None)
    p.add_argument("--exclude-slow", default=None)
    p = sub.add_parser("compare", help="a run against another tree's reference")
    common(p)
    p.add_argument("--tag", default="run1", help="the run file <tag>.jsonl.gz of the run folder, or `reference`")
    p.add_argument("--against", required=True, help="the tree whose reference is the base")
    p.add_argument("--harmonize", action="store_true", help="compare the harmonize part (<tag> default `harmonize`)")
    p = sub.add_parser("reference", help="make a tree's reference answers")
    common(p)
    p.add_argument("--timeout", type=int, default=60)
    p.add_argument("--slow-cap", type=int, default=900)
    p.add_argument("--compare-with", default=None, help="a tree whose reference the new one is compared with")
    p.add_argument("--no-network", action="store_true", help="leave out resolve_organism (needs Ensembl over the network)")
    p = sub.add_parser("summarize", help="answers and exceptions of a run by entry and stratum")
    common(p)
    p.add_argument("--tag", default="run1")
    a = ap.parse_args(argv)
    if a.submit:
        submit(a, argv)
        return 0
    if a.cmd == "make":
        tree_process(a.tree, "idtrack_reproducibility.corpus.make_corpus",
                     ["--tree", a.tree, "--organism", a.organism, "--out", os.path.abspath(a.out), "--seed", str(a.seed)]
                     + (["--tables"] if a.tables else []), 0,
                     os.path.join(C.run_dir(a.tree, a.organism), "logs", "make.log"))
    elif a.cmd == "run":
        extra = (["--no-network"] if a.no_network else []) + (["--shuffle", str(a.shuffle)] if a.shuffle is not None
                                                               else []) + (["--entries", a.entries] if a.entries else []) \
            + (["--limit", str(a.limit)] if a.limit else [])
        print(do_run(a.tree, a.organism, a.tag, a.hash_seed, a.timeout, extra, a.corpus))
    elif a.cmd == "harmonize":
        print(do_harmonize(a.tree, a.organism, a.tag, a.hash_seed, a.corpus, a.exclude_slow))
    elif a.cmd == "compare":
        d = C.run_dir(a.tree, a.organism)
        if a.harmonize:
            tag = a.tag if a.tag != "run1" else "harmonize"
            b = C.reference_file(a.organism, a.tree, True) if tag == "reference" else os.path.join(d, f"{tag}.jsonl.gz")
            base, corpus = C.reference_file(a.organism, a.against, True), None
        else:
            tag = a.tag
            b = C.reference_file(a.organism, a.tree) if tag == "reference" else run_file(a.tree, a.organism, tag)
            base, corpus = C.reference_file(a.organism, a.against), C.corpus_file(a.organism)
        os.makedirs(d, exist_ok=True)
        print(write_compare(base, b, corpus, os.path.join(d, f"compare_{tag}_vs_{a.against}")))
    elif a.cmd == "reference":
        do_reference(a.tree, a.organism, a.timeout, a.slow_cap, a.compare_with, a.no_network)
    elif a.cmd == "summarize":
        from idtrack_reproducibility.corpus import summarize_results
        f = run_file(a.tree, a.organism, a.tag)
        summarize_results.main([f, "--corpus", C.corpus_file(a.organism),
                                "--json", os.path.join(C.run_dir(a.tree, a.organism), f"{a.tag}.summary.json")])
    return 0


if __name__ == "__main__":
    sys.exit(main())
