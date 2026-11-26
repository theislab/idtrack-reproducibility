"""`idtrack-repro known-limitations`: the known-limitations table of every family, from the current scored runs.

Each family's builder measures its rows from the causes and reports of its runs (docs/methods.md, "Known
limitations"); this command runs them on the runs `results.toml` names and writes `results/limitations/<family>_<tree>.tsv`
and their union `results/limitations/limitations.tsv` (the common columns, `limitations.HEADER`), each with a
`.provenance` file (its sources and their sha256) that `idtrack-repro results` reads. No number is typed.

The builders: the harness's sets (`f0` runs; F0's map of the residual causes), `astar` (`limitations_astar`), `a`
(`limits_a`: the causes with A's map plus the undetermined answers by reason; A0's summaries with `--a0`), `e`
(`limits_e`, judge `now`: what users get from the tree), `c` (F0's builder with C's map), `d` (`limits_d`); and
`f` (`limits_f`: F's causes with its map, the undetermined answers by reason and the checks that do not
hold, from each run's scorings `scores_gene` and `scores_symbol`); `g` (`limits_g`: idtrack's wrong lines of G's read
sample, by known limitation and G's map of new causes; another tool's runs have no rows); `b` (`limits_b`: B's map of the wrong answers'
classes and B2's names whose old gene ended, from `scores_now` and `report_b.tsv`); and `outcomes` for a registry entry that says `limits = "outcomes"` (`limits_outcomes`: the class of a wrong answer read
from its stratum and the answer, with the family's `limitations_map_outcomes.tsv`; the scored folders of the entry's
`scored`), for sets whose wrong answers need no cause classification.
"""
import argparse
import csv
import datetime
import glob
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
    "a0": ("a0", "{py} -m idtrack_reproducibility.families.a.limits_a --a0 --out {out} {runs}"),
    "e": ("e", "{py} -m idtrack_reproducibility.families.e.limits_e --runs {runs} --judge now --out {out}"),
    "c": ("c", "{py} -m idtrack_reproducibility.families.limitations --family c --scored {runs} --map "
               "{src}/idtrack_reproducibility/families/c/limitations_map_c.tsv --out {out}"),
    "d": ("d", "{py} -m idtrack_reproducibility.families.d.limits_d --out {out} {runs}"),
    "f": ("f", "{py} -m idtrack_reproducibility.families.f.limits_f --out {out} {runs}"),
    "b": ("b", "{py} -m idtrack_reproducibility.families.b.limits_b --out {out} {runs}"),
    "g": ("g", "{py} -m idtrack_reproducibility.families.g.limits_g --out {out} {runs}"),
    "outcomes": ("{family}", "{py} -m idtrack_reproducibility.families.limits_outcomes --family {family} --map "
                             "{src}/idtrack_reproducibility/families/{family}/limitations_map_outcomes.tsv --out {out} {runs}"),
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
    union, prov_all, seen = [], [], set()
    now = datetime.datetime.now().isoformat(timespec="seconds")
    for e in reg.get("run", []):
        fam = e["family"]
        if e.get("subject"):
            continue  # another tool's runs (family G): the table states idtrack's limitations, read from idtrack's runs
        if fam == "a" and "/a0/" in e["folder"]:
            fam = "a0"  # A0, graph fidelity: its own builder over its summaries
        elif "/runs/" not in e["folder"] and not (e["folder"].startswith("modest/") and e.get("limits")):
            continue  # (a part of a modest run is read where its entry names a builder of rows from the outcomes)
        fam = e.get("limits", fam)
        if fam not in BUILDERS or (fam == "astar" and e["tree"] != "d4c2dda"):
            continue  # (A*'s builder reads every organism's run of one name; its transcript runs have no causes)
        name, tmpl = BUILDERS[fam]
        name = name.format(family=e["family"])
        out = os.path.join(out_dir, f"{name}_{e['tree']}.tsv")
        runs = [os.path.join(data, e["folder"].format(organism=o)) for o in e["organisms"]]
        if fam == "outcomes":  # the scored folders of the entry (each with run.json, report.tsv and outcomes)
            runs = [os.path.dirname(os.path.join(r, (e.get("scored") or ["run.json"])[0])) for r in runs]
            runs = [r for r in runs if os.path.exists(os.path.join(r, "run.json"))]
        runs = [r for r in runs if fam == "outcomes" or os.path.exists(os.path.join(r, "causes.tsv.gz")) or
                os.path.exists(os.path.join(r, "scores_now", "causes.tsv.gz")) or
                (fam == "f" and glob.glob(os.path.join(r, "scores_*", "causes.tsv.gz"))) or
                (fam == "a0" and os.path.exists(os.path.join(r, "summary.tsv"))) or
                (fam == "g" and os.path.exists(os.path.join(r, "run.json")))]  # a run whose answers were classified
        # (G: its answers are classified by the reading of its sample, report/reading.tsv, a source below)
        if not runs or out in seen:
            continue
        seen.add(out)
        if a.only and e["family"] not in a.only.split(","):  # another family's table, kept as it is in the union
            if os.path.exists(out):
                union += RES.read_tsv(out)
                if os.path.exists(out + ".provenance"):
                    prov_all += [dict(r, file="results/limitations/limitations.tsv") for r in RES.read_tsv(out + ".provenance")]
            continue
        run_name = e["folder"].split("/runs/")[1].split("/")[0] if "/runs/" in e["folder"] else e["tree"]
        cmd = tmpl.format(py=sys.executable, runs=" ".join(runs), out=out, run=run_name, src=CFG.SRC, family=e["family"])
        subprocess.run(cmd, shell=True, check=True)
        rows = RES.read_tsv(out)
        union += rows
        srcs = []
        for r in runs:
            for f in ("causes.tsv.gz", "scores_now/causes.tsv.gz", "report.tsv", "scores_now/report.tsv",
                      "undetermined.tsv.gz", "e_summary.tsv", "summary.tsv") + (("run.json",) if fam in ("outcomes", "g") else ()) + \
                    (tuple(f"scores_{c}/{x}" for c in ("gene", "symbol") for x in ("causes.tsv.gz", "report.tsv",
                                                                                  "undetermined.tsv.gz")) + ("f_checks.tsv",)
                     if fam == "f" else ()):
                p = os.path.join(r, f)
                if os.path.exists(p):
                    srcs.append(p)
        if fam == "g":
            srcs.append(os.path.join(data, "families", "g", "report", "reading.tsv"))
        prov = [dict(file=os.path.relpath(out, ROOT), sha256="", generated_by=cmd.replace(sys.executable, "python")
                     .replace(data + os.sep, "<data>/").replace(CFG.SRC + os.sep, "").replace(ROOT + os.sep, ""), source=os.path.relpath(p, data), source_sha256=RES.sha256(p),
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
