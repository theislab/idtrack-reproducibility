"""E9 the web layer: idtrack-web's api, served from local graphs, against the package's own answers, and its limits.

Usage: web_e.py --out <run folder> [--n 10000] [--organisms a,b]    (envs/webapp: idtrack-web's environment, whose
       idtrack is the commit it pins: 71690ab)

No container runs on this cluster: the api is served as its SETUP.md section 8 runs it locally (`idtrack-web serve`:
gunicorn, two forked workers, a directory store of links to the graphs made by `idtrack-web publish --link`), on
127.0.0.1, and asked over HTTP from this process, which holds the same graph files through `idtrack.API.load`.

1. Answers (decision W2: nothing rounded or reworded): per organism up to `--n` lines of every family with
   published sets (the harness's sets, family E's `convert` lines, A*'s if present), seeded; each asked as
   `GET /organisms/{o}/convert` with its arguments, and in batches of 100 lines of equal arguments as
   `POST .../convert/batch`; the JSON must equal `API.convert(...).to_dict()` (resp. `convert_many`), every float
   with the same bits. An identifier the api refuses by its documented validation (empty, longer than
   `identifier_max_length`, control characters) must be refused with `400 invalid_identifier`; an argument the
   package refuses must be refused with the documented code.
2. Limits (SETUP.md section 6 defaults; a second server with a small `response_max_bytes`, a third with the
   default rate): `batch_max` (1000 accepted, 1001 -> 413 `too_many_identifiers`), `identifier_max_length` (256
   accepted, 257 -> 400 `invalid_identifier`), `body_max_bytes` (1 MiB + 1 -> 413 `payload_too_large`),
   `response_max_bytes` (-> 422 `response_too_large`, never a cut answer), the rate (burst 30, then 429
   `rate_limited` with `Retry-After`), an unknown parameter (400 `unknown_parameter`), a release as a JSON string (400
   `invalid_request`), an unknown organism (404 `unknown_organism`), a body that is not JSON (415).

Writes `<out>/e9/report.tsv` and `<out>/e9/differences.jsonl`.
"""
import argparse
import collections
import glob
import json
import math
import os
import random
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
from idtrack_reproducibility.families import common as C
from idtrack_reproducibility.families import lines as L

def web_command():
    """The idtrack-web command (configuration `[tools] idtrack_web`; default `idtrack-web` on the PATH)."""
    from idtrack_reproducibility import config as CFG
    return CFG.get().raw.get("tools", {}).get("idtrack_web", "idtrack-web")
P = "/api/v1"


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def floats_of(d):
    if isinstance(d, float):
        return [d]
    if isinstance(d, dict):
        return [f for v in d.values() for f in floats_of(v)]
    if isinstance(d, list):
        return [f for v in d for f in floats_of(v)]
    return []


def same(web, ref):
    return web == ref and [f.hex() for f in floats_of(web)] == [f.hex() for f in floats_of(ref)] and all(
        math.isfinite(f) for f in floats_of(web))


class Server:
    def __init__(self, store, logs, name, **env):
        self.port = free_port()
        e = dict(os.environ, IDTRACK_WEB_SOURCE="directory", IDTRACK_WEB_DIRECTORY=store, IDTRACK_WEB_HOST="127.0.0.1",
                 IDTRACK_WEB_PORT=str(self.port), IDTRACK_WEB_WORKERS="2", IDTRACK_WEB_LOG_LEVEL="WARNING")
        e.update({f"IDTRACK_WEB_{k.upper()}": str(v) for k, v in env.items()})
        self.log = open(os.path.join(logs, f"server_{name}.log"), "w")
        self.p = subprocess.Popen([web_command(), "serve"], env=e, stdout=self.log, stderr=subprocess.STDOUT)
        self.url = f"http://127.0.0.1:{self.port}{P}"
        import requests
        t0 = time.time()
        while time.time() - t0 < 900:
            try:
                if requests.get(self.url + "/health", timeout=5).status_code == 200:
                    return
            except Exception:  # noqa: BLE001
                pass
            if self.p.poll() is not None:
                raise SystemExit(f"server {name} exited: {self.p.returncode}; see {self.log.name}")
            time.sleep(2)
        raise SystemExit(f"server {name} did not start")

    def stop(self):
        self.p.terminate()
        try:
            self.p.wait(60)
        except subprocess.TimeoutExpired:
            self.p.kill()


def params(q):
    p = {"identifier": q["id"]}
    for k in ("from_release", "to_release", "final_database", "from_database"):
        if q.get(k) is not None:
            p[k] = str(q[k])
    return p


def sample_lines(org, n, rng):
    fams = []
    hl, _, _ = L.load("harness", org)
    fams.append(("harness", hl))
    for fam in ("e", "astar", "a"):
        d = L.query_dir(fam, org)
        jsonl = [s for s in L.set_names(fam, org) if not s.endswith("_dec") and os.path.exists(os.path.join(d, s + ".jsonl.gz"))]
        if jsonl:  # (this environment has no pyarrow: the sets written as jsonl.gz)
            ln, _, _ = L.load(fam, org, jsonl)
            ln = [q for q in ln if (q.get("call") or "convert") == "convert" and not q.get("kwargs")
                  and (q.get("id_codec") or "str") == "str" and isinstance(q.get("id"), str) and len(q["id"]) < 5000]
            fams.append((fam, ln))
    out = []
    for fam, ln in fams:
        ln.sort(key=lambda q: q["qid"])
        for q in rng.sample(ln, min(n, len(ln))):
            out.append((fam, q))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--organisms", default=",".join(C.ORGANISMS))
    a = ap.parse_args()
    import requests
    out = os.path.join(os.path.abspath(a.out), "e9")
    os.makedirs(out, exist_ok=True)
    store = os.path.join(out, "store")
    from idtrack_reproducibility import config as CFG
    graphs = [glob.glob(os.path.join(CFG.get().graphs, o, "graph_*.npz"))[0] for o in a.organisms.split(",")]
    if not os.path.exists(os.path.join(store, "idtrack")):
        cmd = [web_command(), "publish", "--directory", store, "--link"]
        for g in graphs:
            cmd += ["--graph", g]
        subprocess.run(cmd, check=True, stdout=open(os.path.join(out, "publish.log"), "w"), stderr=subprocess.STDOUT)
    import idtrack
    apis = {}
    repo = os.path.join(out, "repo")
    os.makedirs(repo, exist_ok=True)
    for o, g in zip(a.organisms.split(","), graphs):
        link = os.path.join(repo, os.path.basename(g))
        if not os.path.exists(link):
            os.symlink(g, link)
        apis[o] = idtrack.API.load(link)
    srv = Server(store, out, "main", rate_per_minute=0, search_seconds=30, request_seconds=120)
    lim = json.loads(requests.get(srv.url + "/limits").text)
    stats = collections.Counter()
    diff = open(os.path.join(out, "differences.jsonl"), "w")
    rng = random.Random(20261118)
    max_len = lim.get("identifier_max_length", 256)
    S = requests.Session()
    try:
        for o, api in apis.items():
            batches = collections.defaultdict(list)
            for fam, q in sample_lines(o, a.n, rng):
                k = (o, fam)
                ident = q["id"]
                kw = {x: q.get(x) for x in ("from_release", "to_release", "final_database", "from_database") if q.get(x) is not None}
                bad_id = (not ident.strip()) or len(ident) > max_len or any(ord(c) < 32 or 0x7F <= ord(c) < 0xA0 for c in ident)
                try:
                    ref = api.convert(ident, budget_seconds=30, **kw).to_dict()
                    ref_err = None
                except Exception as e:  # noqa: BLE001
                    ref, ref_err = None, type(e).__name__
                r = S.get(f"{srv.url}/organisms/{o}/convert", params=params(q), timeout=300)
                stats[k + ("asked",)] += 1
                try:
                    body = r.json()
                except Exception:  # noqa: BLE001
                    body = {"_raw": r.text[:300]}
                if bad_id:
                    ok = r.status_code == 400 and (body.get("error") or {}).get("code") in ("invalid_identifier", "missing_value")
                    stats[k + ("refused_by_validation", "ok" if ok else "wrong")] += 1
                elif ref_err is not None:
                    ok = r.status_code in (400, 404, 422) and "error" in body
                    stats[k + ("package_error", "ok" if ok else "wrong")] += 1
                else:
                    ok = r.status_code == 200 and same(body, ref)
                    stats[k + ("equal", "ok" if ok else "wrong")] += 1
                    if not bad_id and not q.get("from_database"):
                        batches[json.dumps(kw, sort_keys=True)].append((fam, q, ref))
                if not ok:
                    diff.write(json.dumps(dict(organism=o, family=fam, qid=q["qid"], status=r.status_code, web=body, ref=ref,
                                               ref_error=ref_err), default=str)[:20000] + "\n")
            # batches of 100 of equal arguments
            for kwj, items in batches.items():
                kw = json.loads(kwj)
                for i in range(0, len(items), 100):
                    chunk = items[i:i + 100]
                    ids = [q["id"] for _, q, _ in chunk]
                    ref = [x.to_dict() for x in api.convert_many(ids, budget_seconds=30, **kw)]
                    r = S.post(f"{srv.url}/organisms/{o}/convert/batch", json=dict(identifiers=ids, **kw), timeout=600)
                    body = r.json() if r.headers.get("content-type", "").startswith("application/json") else None
                    res = body.get("results") if isinstance(body, dict) else body
                    ok = r.status_code == 200 and isinstance(res, list) and len(res) == len(ref) and all(same(w, x) for w, x in zip(res, ref))
                    stats[(o, "batch", "equal", "ok" if ok else "wrong")] += 1
                    if not ok:
                        diff.write(json.dumps(dict(organism=o, family="batch", args=kw, status=r.status_code,
                                                   web=(str(body)[:3000]), n=len(ids)), default=str) + "\n")
        # ---- limits ----
        o = list(apis)[0]
        api = apis[o]
        some = sorted(n for n in api.graph.nodes if C.ENS_ID.match(n) and "." in n)[:2000]

        def check(name, resp, status, code=None):
            try:
                b = resp.json()
            except Exception:  # noqa: BLE001
                b = {}
            got_code = (b.get("error") or {}).get("code") if isinstance(b, dict) else None
            ok = resp.status_code == status and (code is None or got_code == code)
            stats[("limits", name, "ok" if ok else "wrong")] += 1
            if not ok:
                diff.write(json.dumps(dict(limit=name, status=resp.status_code, code=got_code, expected=[status, code],
                                           body=str(b)[:500])) + "\n")
        bm = lim.get("batch_max", 1000)
        check("batch_max_at", S.post(f"{srv.url}/organisms/{o}/convert/batch", json={"identifiers": some[:bm]}, timeout=600), 200)
        check("batch_max_beyond", S.post(f"{srv.url}/organisms/{o}/convert/batch", json={"identifiers": some[:bm] + ["x"]}), 413, "too_many_identifiers")
        check("identifier_length_at", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "A" * max_len}), 200)
        check("identifier_length_beyond", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "A" * (max_len + 1)}), 400, "invalid_identifier")
        check("identifier_control", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "TP53\x00"}), 400, "invalid_identifier")
        check("identifier_empty", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": ""}), 400)
        bb = lim.get("body_max_bytes", 1048576)
        big = '{"identifiers": ["A"], "pad": "' + "x" * (bb + 1 - len('{"identifiers": ["A"], "pad": ""}')) + '"}'
        check("body_beyond", S.post(f"{srv.url}/organisms/{o}/convert/batch", data=big.encode(), headers={"Content-Type": "application/json"}), 413, "payload_too_large")
        check("unknown_parameter", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "x", "nope": "1"}), 400, "unknown_parameter")
        check("release_as_string_in_body", S.post(f"{srv.url}/organisms/{o}/convert/batch", json={"identifiers": ["x"], "to_release": "115"}), 400, "invalid_request")
        check("unknown_organism", S.get(f"{srv.url}/organisms/zebra/convert", params={"identifier": "x"}), 404, "unknown_organism")
        check("not_json", S.post(f"{srv.url}/organisms/{o}/convert/batch", data=b"identifiers=x", headers={"Content-Type": "text/plain"}), 415, "unsupported_media_type")
        check("unknown_release", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "x", "to_release": "999"}), 400, "unknown_release")
        check("unknown_database", S.get(f"{srv.url}/organisms/{o}/convert", params={"identifier": "x", "final_database": "Nope"}), 400, "unknown_database")
        srv.stop()
        small = Server(store, out, "small_response", rate_per_minute=0, response_max_bytes=65536)
        check("response_small_ok", requests.post(f"{small.url}/organisms/{o}/convert/batch", json={"identifiers": some[:3]}), 200)
        check("response_beyond", requests.post(f"{small.url}/organisms/{o}/convert/batch", json={"identifiers": some[:1000], "explain": True}),
              422, "response_too_large")
        small.stop()
        rate = Server(store, out, "rate", rate_per_minute=120, rate_burst=30)
        codes = [requests.get(f"{rate.url}/organisms/{o}/final-databases").status_code for _ in range(45)]
        n_ok = sum(c == 200 for c in codes[:30])
        stats[("limits", "rate_burst_accepted", "ok" if n_ok == 30 else "wrong")] += 1
        r = requests.get(f"{rate.url}/organisms/{o}/final-databases")
        check("rate_beyond", r, 429, "rate_limited")
        stats[("limits", "rate_retry_after_header", "ok" if r.headers.get("Retry-After") else "wrong")] += 1
        h = requests.get(f"{rate.url}/health")
        stats[("limits", "health_not_rate_limited", "ok" if h.status_code == 200 else "wrong")] += 1
        rate.stop()
    finally:
        try:
            srv.stop()
        except Exception:  # noqa: BLE001
            pass
        diff.close()
    with open(os.path.join(out, "report.tsv"), "w") as fo:
        fo.write("key\tcount\n")
        for k, v in sorted(stats.items(), key=lambda kv: str(kv[0])):
            fo.write("|".join(str(x) for x in k) + f"\t{v}\n")
    json.dump(dict(limits=lim, idtrack=idtrack.__file__, finished=time.strftime("%Y-%m-%d %H:%M:%S")), open(os.path.join(out, "info.json"), "w"), indent=1)
    print(open(os.path.join(out, "report.tsv")).read())


if __name__ == "__main__":
    main()
