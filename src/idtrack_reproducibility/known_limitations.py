"""`idtrack-repro known-limitations`: the known-limitations table of every family, from the current scored runs.

Each family's builder measures its rows from the causes and reports of its runs (docs/methods.md, "Known
limitations"); this command runs them on the runs `results.toml` names and writes `results/limitations/<family>_<tree>.tsv`
and their union `results/limitations/limitations.tsv` (the common columns, `limitations.HEADER`), each with a
`.provenance` file (its sources and their sha256) that `idtrack-repro results` reads. No number is typed.

The builders: the harness's sets (`f0` runs; F0's map of the residual causes), `astar` (`limitations_astar`), `a`
(`limits_a`: the causes with A's map plus the undetermined answers by reason), `e` (`limits_e`, judge `now`: what users
get from the tree).
"""
import argparse
import csv
import datetime
import os
import subprocess
import sys
import tempfile

from . import config as CFG
from . import results as RES

BUILDERS = {
    # family in results.toml -> (family in the table, command template; {runs} the scored folders, {out} the table)
    "f0": ("f0", "{py} -m idtrack_reproducibility.families.limitations --family f0 --scored {runs} --out {out}"),
    "astar": ("astar", "{py} -m idtrack_reproducibility.families.astar.limitations_astar {run} --out {out}"),
    "a": ("a", "{py} -m idtrack_reproducibility.families.a.limits_a --out {out} {runs}"),
    "e": ("e", "{py} -m idtrack_reproducibility.families.e.limits_e --runs {runs} --judge now --out {out}"),
    "c": ("c", "{py} -m idtrack_reproducibility.families.limitations --family c --scored {runs} --out {out}"),
    "d": ("d", "{py} -m idtrack_reproducibility.families.d.limits_d --out {out} {runs}"),
}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idtrack-repro known-limitations", description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", default=None, help="comma-separated families")
    a = ap.parse_args(argv)
    from .check import ROOT
    reg = RES.registry(ROOT)
    data = CFG.get().root
    out_dir = os.path.join(ROOT, "results", "limitations")
    os.makedirs(out_dir, exist_ok=True)
    union, prov_all = [], []
    now = datetime.datetime.now().isoformat(timespec="seconds")
    for e in reg.get("run", []):
        fam = e["family"]
        if fam not in BUILDERS or (a.only and fam not in a.only.split(",")) or "/runs/" not in e["folder"]:
            continue
        if fam == "astar" and e["tree"] != "d4c2dda":
            continue  # A*'s builder reads every organism's run of one name (its transcript runs have no causes)
        name, tmpl = BUILDERS[fam]
        runs = [os.path.join(data, e["folder"].format(organism=o)) for o in e["organisms"]]
        runs = [r for r in runs if os.path.exists(os.path.join(r, "causes.tsv.gz")) or
                os.path.exists(os.path.join(r, "scores_now", "causes.tsv.gz"))]  # a run whose answers were classified
        if not runs:
            continue
        out = os.path.join(out_dir, f"{name}_{e['tree']}.tsv")
        run_name = e["folder"].split("/runs/")[1].split("/")[0]
        cmd = tmpl.format(py=sys.executable, runs=" ".join(runs), out=out, run=run_name)
        subprocess.run(cmd, shell=True, check=True)
        rows = RES.read_tsv(out)
        union += rows
        srcs = []
        for r in runs:
            for f in ("causes.tsv.gz", "scores_now/causes.tsv.gz", "report.tsv", "scores_now/report.tsv",
                      "undetermined.tsv.gz", "e_summary.tsv"):
                p = os.path.join(r, f)
                if os.path.exists(p):
                    srcs.append(p)
        prov = [dict(file=os.path.relpath(out, ROOT), sha256="", generated_by=cmd.replace(sys.executable, "python")
                     .replace(data + os.sep, "<data>/").replace(ROOT + os.sep, ""), source=os.path.relpath(p, data), source_sha256=RES.sha256(p),
                     scored_at="", commit=e["commit"], graph_digest="", generated_at=now) for p in srcs]
        RES.write_tsv(out + ".provenance", RES.PROV_COLS, prov)
        prov_all += [dict(r, file="results/limitations/limitations.tsv") for r in prov]
        print(f"{name} {e['tree']}: {len(rows)} rows")
    from .families import limitations as LIM
    path = os.path.join(out_dir, "limitations.tsv")
    with open(path + ".tmp", "w") as fo:
        fo.write("\t".join(LIM.HEADER) + "\n")
        for r in union:
            fo.write("\t".join(str(r.get(h, "")) for h in LIM.HEADER) + "\n")
    os.replace(path + ".tmp", path)
    RES.write_tsv(path + ".provenance", RES.PROV_COLS, prov_all)
    print(f"{len(union)} rows -> results/limitations/limitations.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
