"""Running work as jobs: one interface, two runners (`local` and `slurm`), chosen by the configuration.

`submit(what, command, ...)` returns a job id. With the `slurm` runner the command is an `sbatch --wrap` job in the
configured partition and QOS, named `<name_prefix><what>` (a rerun of the same work keeps its name), its log in
`<log_dir>/<name>_<id>.out`; dependencies are `afterany` (a preempted or failed predecessor does not hold the queue:
each step checks its inputs). With the `local` runner the command runs at once in this process's environment, its
output in the same log file, and dependencies are met by order; the id is `local:<n>`.

Every job's environment: `PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 IDTRACK_ENSEMBL_ACCESS=https`,
`PYTHONPATH` with this package's source (so a tree's own Python, which does not have the package installed, imports
it), and `IDTRACK_REPRO_CONFIG` (the configuration file of the submitting process).
"""
import os
import shlex
import subprocess

from . import config as CFG

_local_n = [0]


def job_env(extra=None):
    cfg = CFG.get()
    env = {"PYTHONDONTWRITEBYTECODE": "1", "PYTHONUNBUFFERED": "1", "PYTHONHASHSEED": "0",
           "IDTRACK_ENSEMBL_ACCESS": "https", "PYTHONPATH": CFG.SRC}
    if cfg.path:
        env["IDTRACK_REPRO_CONFIG"] = cfg.path
    for k in ("IDTRACK_REPRO_DATA_ROOT",) + tuple(f"IDTRACK_REPRO_{x.upper()}" for x in CFG.LOCATIONS):
        if os.environ.get(k):
            env[k] = os.path.abspath(os.environ[k])
    env.update(extra or {})
    return env


def module(python, name, args=""):
    """The command that runs a module of this package with a given Python (`python -m <module> <args>`)."""
    return f"{python} -m {name} {args}".strip()


def submit(what, command, log_dir, cpus=2, mem_gb=8, time="04:00:00", after=None, env=None):
    """Run `command` (a shell command line) as a job; returns its id."""
    cfg = CFG.get()
    name = f"{cfg.name_prefix}{what}"
    os.makedirs(log_dir, exist_ok=True)
    e = job_env(env)
    if cfg.runner == "local":
        _local_n[0] += 1
        jid = f"local:{_local_n[0]}"
        log = os.path.join(log_dir, f"{name}_{jid.replace(':', '')}.out")
        with open(log, "w") as fo:
            r = subprocess.run(command, shell=True, stdout=fo, stderr=subprocess.STDOUT, env=dict(os.environ, **e))
        if r.returncode:
            raise SystemExit(f"{name} failed (exit {r.returncode}); log: {log}")
        return jid
    exports = " ".join(f"{k}={shlex.quote(v)}" for k, v in e.items())
    cmd = ["sbatch", "--parsable", "-p", cfg.slurm["partition"]]
    if cfg.slurm["qos"]:
        cmd += ["--qos", cfg.slurm["qos"]]
    cmd += ["-c", str(cpus), "--mem", f"{int(mem_gb)}G", "-t", time, f"--job-name={name}",
            f"--output={log_dir}/%x_%j.out"] + cfg.slurm["extra"]
    deps = [j for j in (after or []) if j and not str(j).startswith("local:")]
    if deps:
        cmd.append(f"--dependency=afterany:{':'.join(str(j) for j in deps)}")
    cmd.append(f"--wrap=export {exports}; {command}")
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    return out.split(";")[0]


def memory_for(organism, default=24):
    return CFG.get().memory_gb.get(organism, default)
