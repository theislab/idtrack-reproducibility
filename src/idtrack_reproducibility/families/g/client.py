"""The client of family G: one per host, at its polite rate, with backoff, every response stored as fetched before it is
parsed (docs/families/g.md, "The client"; docs/records.md, "Another tool").

- **Hosts** (`HOSTS`): one session per host group, with one connection (a pool of one, blocking), and the interval
  between two requests: half the host's documented rate, one request per second where none is documented. Hosts of one
  service share a group (every archive of Ensembl's BioMart is one group), so a pass asks a service through one
  connection at its rate whatever host name it uses. A pass of a host is one process (one job per host): one
  connection per host over all jobs together.
- **Backoff**: a 429 or 503 is waited out for its `Retry-After` (seconds or an HTTP date; without it 60 s, doubling to
  15 minutes) and asked again, for at most `MAX_WAIT` seconds in all; other 5xx answers, time-outs and connection errors
  wait 30 s, 2, 8 and 30 minutes, then the request is recorded with its error, and the next pass asks it again. Every
  other status (200, and the 400 or 404 a service answers for an unknown identifier) is an answer and is stored.
- **Store** (`RawStore`): `<families>/g/raw/<tool>/<version>/<organism>/`: the body of each response as it came (after
  the transfer's content encoding is undone) in a file named by its sha256 (`<sha256>.<ext>`, so a body that comes twice
  is kept once), then the request's row (`records.RAW_COLUMNS`) appended to `index.jsonl` and its response headers to
  `headers.jsonl`; `compile()` writes `index.parquet` from `index.jsonl` (the last row of a request id wins), the
  index `records.load_raw` reads.
- **Resume**: a request id is the sha256 of the method, URL, parameters and body; a request whose id has a stored row
  without error is not asked again (`fetch` returns that row). `answered(endpoint)` gives the qids that a stored
  response of an endpoint holds, so that a pass made of other batches (a pilot's) asks only the lines without one.
- Nothing of the person running it is sent: the User-Agent is a generic client name (`USER_AGENT`), with no contact.
"""
import datetime
import email.utils
import hashlib
import json
import os
import time
import urllib.parse

from idtrack_reproducibility.families import common as C

USER_AGENT = "gene-id-benchmark-client/1.0"
MAX_WAIT = 6 * 3600  # the longest a refused request is waited out in all before it is recorded with its error
RETRY_WAITS = (30, 120, 480, 1800)  # waits after a 5xx other than 503, a time-out or a connection error
TIMEOUT = (30, 600)  # connect, read (BioMart answers a query within its five minutes)


class Host:
    """A host group: its polite interval between requests and why (the documented rate, quoted in g.md)."""

    def __init__(self, name, per_second, why):
        self.name, self.interval, self.why = name, 1.0 / per_second, why


# The documented rates (docs/families/g.md, "The tools", with the quotes and their addresses), halved; one request per
# second where none is documented.
HOSTS = {
    "ensembl_rest": Host("ensembl_rest", 7.5, "Ensembl REST: 15 requests per second (55,000 per hour) documented; half"),
    "ensembl_biomart": Host("ensembl_biomart", 1.0, "Ensembl BioMart (every archive host): no rate documented"),
    "gprofiler": Host("gprofiler", 1.0, "g:Profiler: no rate documented"),
    "mygene": Host("mygene", 1.0, "MyGene.info: 'one per second' for massive queries (its blog); no number in its documents"),
    "hgnc": Host("hgnc", 5.0, "HGNC REST: 10 requests per second documented; half"),
    "uniprot": Host("uniprot", 1.0, "UniProt REST: no rate documented"),
    "bridgedb": Host("bridgedb", 1.0, "BridgeDb web service: no rate documented"),
    "ensembl_ftp": Host("ensembl_ftp", 1.0, "Ensembl's FTP site over HTTPS (pyensembl's downloads): no rate documented"),
}


def host_group(url):
    h = urllib.parse.urlparse(url).hostname or ""
    if h == "rest.ensembl.org":
        return "ensembl_rest"
    if h.endswith("archive.ensembl.org") or h in ("www.ensembl.org", "useast.ensembl.org", "asia.ensembl.org"):
        return "ensembl_biomart"
    if h == "ftp.ensembl.org":
        return "ensembl_ftp"
    if h == "biit.cs.ut.ee":
        return "gprofiler"
    if h.endswith("mygene.info"):
        return "mygene"
    if h == "rest.genenames.org":
        return "hgnc"
    if h == "rest.uniprot.org":
        return "uniprot"
    if h.endswith("bridgedb.org"):
        return "bridgedb"
    raise ValueError(f"{url}: a host family G does not ask (add it to HOSTS with its documented rate)")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def request_id(method, url, params=None, body=None, salt=None):
    """The request's id: the sha256 of its method, URL, parameters and body (and `salt`: the pass of a version call,
    which every pass asks again)."""
    doc = json.dumps([method.upper(), url, params or {}, body] + ([salt] if salt else []), sort_keys=True,
                     separators=(",", ":"), default=str)
    return hashlib.sha256(doc.encode()).hexdigest()[:32]


def retry_after(value, default):
    """Seconds to wait from a `Retry-After` header (seconds, possibly fractional, or an HTTP date)."""
    if not value:
        return default
    try:
        return max(0.0, float(value))
    except ValueError:
        pass
    try:
        when = email.utils.parsedate_to_datetime(value)
        return max(0.0, (when - datetime.datetime.now(when.tzinfo)).total_seconds())
    except (TypeError, ValueError):
        return default


def usable(row):
    """A stored response that answers its lines: no error recorded, and not a batch the service refused as a whole
    (a status of 400 or more to a request of several lines: MyGene.info's "Too Many Hits", BioMart's 405). A single
    line's 400 or 404 is the service's answer for that input and is kept as one. An unusable response stays stored;
    the next pass asks its lines again under a new request."""
    if row.get("error"):
        return False
    st = row.get("http_status")
    return not (st is not None and st >= 400 and len(row.get("qids") or []) > 1)


def raw_dir(tool, version, organism):
    return os.path.join(C.family_dir("g"), "raw", tool, str(version), organism)


class RawStore:
    """The raw responses of one tool, version and organism (module docstring, "Store")."""

    def __init__(self, folder, tool, version, organism):
        self.folder, self.tool, self.version, self.organism = folder, tool, str(version), organism
        os.makedirs(folder, exist_ok=True)
        self.index_path = os.path.join(folder, "index.jsonl")
        self.rows = {}
        if os.path.exists(self.index_path):
            with open(self.index_path) as fh:
                for raw in fh:
                    if raw.strip():
                        r = json.loads(raw)
                        self.rows[r["request_id"]] = r

    @classmethod
    def of(cls, tool, version, organism):
        return cls(raw_dir(tool, version, organism), tool, version, organism)

    def done(self, rid):
        r = self.rows.get(rid)
        return r if r and usable(r) else None

    def answered(self, endpoint):
        """qids that a usable stored response of this endpoint holds (`usable`)."""
        out = set()
        for r in self.rows.values():
            if r["endpoint"] == endpoint and usable(r):
                out.update(r.get("qids") or [])
        return out

    def put(self, row, content, ext, headers):
        if content is not None:
            sha = hashlib.sha256(content).hexdigest()
            name = f"{sha}.{ext}"
            path = os.path.join(self.folder, name)
            if not os.path.exists(path):
                tmp = f"{path}.{os.getpid()}.{time.time_ns()}.tmp"  # one name per write: no temp is written twice
                with open(tmp, "wb") as fh:
                    fh.write(content)
                try:
                    os.link(tmp, path)  # a new name for the bytes; fails where the response exists, never replaces it
                except FileExistsError:
                    pass  # the same bytes (the name is their sha256) were stored meanwhile by another write
                # the temp name goes, never a response: its bytes are the response stored under `path`
                os.unlink(tmp)
            row.update(bytes=len(content), sha256=sha, file=name)
        with open(self.index_path, "a") as fh:
            fh.write(json.dumps(row, sort_keys=True, default=str) + "\n")
        with open(os.path.join(self.folder, "headers.jsonl"), "a") as fh:
            fh.write(json.dumps({"request_id": row["request_id"], "headers": headers or {}}, sort_keys=True) + "\n")
        self.rows[row["request_id"]] = row
        return row

    def body(self, row):
        with open(os.path.join(self.folder, row["file"]), "rb") as fh:
            return fh.read()

    def compile(self):
        """index.parquet from index.jsonl (records.RAW_COLUMNS; the last row of each request id)."""
        import pyarrow as pa
        import pyarrow.parquet as pq
        from idtrack_reproducibility import records
        cols = records.RAW_COLUMNS
        rows = []
        for r in self.rows.values():
            d = {c: r.get(c) for c in cols}
            for c in ("params", "body"):
                if d[c] is not None and not isinstance(d[c], str):
                    d[c] = json.dumps(d[c], sort_keys=True, separators=(",", ":"))
            d["qids"] = list(d["qids"] or [])
            rows.append(d)
        rows.sort(key=lambda d: (d["asked_at"] or "", d["request_id"]))
        path = os.path.join(self.folder, "index.parquet")
        pq.write_table(pa.Table.from_pylist(rows, schema=records.raw_schema()), path + ".tmp")
        os.replace(path + ".tmp", path)
        return len(rows)


class Client:
    """Asks the hosts of family G (module docstring). `sleep` and `session_factory` are for tests."""

    def __init__(self, session_factory=None, sleep=time.sleep, clock=time.monotonic, log=None):
        self._sessions, self._last = {}, {}
        self._session_factory = session_factory or _session
        self.sleep, self.clock = sleep, clock
        self.log = log or (lambda *a: None)
        self.asked = 0

    def session(self, group):
        if group not in self._sessions:
            self._sessions[group] = self._session_factory()
        return self._sessions[group]

    def _wait_turn(self, group):
        last = self._last.get(group)
        if last is not None:
            gap = HOSTS[group].interval - (self.clock() - last)
            if gap > 0:
                self.sleep(gap)
        self._last[group] = self.clock()

    def fetch(self, store, method, url, *, params=None, data=None, json_body=None, qids=(), endpoint=None, ext="json",
              headers=None, version=None, salt=None, allow_redirects=True):
        """Ask one request through the store (a stored answer is returned, not asked again); the stored row."""
        body = json_body if json_body is not None else data
        rid = request_id(method, url, params, body, salt)
        old = store.done(rid)
        if old is not None:
            return old
        group = host_group(url)
        hdrs = {"User-Agent": USER_AGENT}
        hdrs.update(headers or {})
        waited, refused_wait, attempt = 0.0, 60.0, 0
        row = dict(request_id=rid, tool=store.tool, tool_version=version or store.version, organism=store.organism,
                   endpoint=endpoint or url, method=method.upper(), url=url, params=params or {}, body=body,
                   qids=list(qids), asked_at=None, http_status=None, seconds=None, bytes=None, sha256=None, file=None,
                   error=None)
        while True:
            self._wait_turn(group)
            row["asked_at"] = now()
            t0 = time.monotonic()
            try:
                r = self.session(group).request(method.upper(), url, params=params, data=data, json=json_body,
                                                headers=hdrs, timeout=TIMEOUT, allow_redirects=allow_redirects)
            except Exception as e:  # noqa: BLE001  (time-outs and connection errors)
                status, resp, err = None, None, f"{type(e).__name__}: {e}"
            else:
                status, resp, err = r.status_code, r, None
            row["seconds"] = round(time.monotonic() - t0, 3)
            self.asked += 1
            if status in (429, 503):
                wait = retry_after(resp.headers.get("Retry-After"), refused_wait)
                refused_wait = min(refused_wait * 2, 900.0)
                if waited + wait > MAX_WAIT:
                    row.update(http_status=status, error=f"refused ({status}) for {int(waited)} s; given up for now")
                    return store.put(row, resp.content, ext, dict(resp.headers))
                self.log("refused", group, status, f"waiting {wait:.0f} s")
                self.sleep(wait)
                waited += wait
                continue
            if status is None or status >= 500:
                if attempt < len(RETRY_WAITS):
                    self.log("failed", group, status, err, f"waiting {RETRY_WAITS[attempt]} s")
                    self.sleep(RETRY_WAITS[attempt])
                    attempt += 1
                    continue
                row.update(http_status=status, error=err or f"HTTP {status}")
                return store.put(row, resp.content if resp is not None else None, ext,
                                 dict(resp.headers) if resp is not None else {})
            row.update(http_status=status)
            return store.put(row, resp.content, ext, dict(resp.headers))


def _session():
    import requests
    from requests.adapters import HTTPAdapter
    s = requests.Session()
    ad = HTTPAdapter(pool_connections=1, pool_maxsize=1, pool_block=True, max_retries=0)
    s.mount("https://", ad)
    s.mount("http://", ad)
    return s


def download(store, url, endpoint=None, qids=(), version=None, client=None, chunk=1 << 20):
    """A file of a local tool (a GTF of pyensembl's): streamed into the store under its sha256, its row with URL, size,
    sha256 and the time; the md5 and the server's Last-Modified in the headers' line. Returns the stored row.

    Each attempt streams into a file of its own (`download_<request id>_<time>.part`); a complete file gets its
    sha256 name as a second name (never over an existing file) and the part name goes; an interrupted or failed
    download keeps its part under its own name and is recorded with its error (the next pass downloads again)."""
    client = client or Client()
    rid = request_id("GET", url)
    old = store.done(rid)
    if old is not None:
        return old
    group = host_group(url)
    client._wait_turn(group)
    row = dict(request_id=rid, tool=store.tool, tool_version=version or store.version, organism=store.organism,
               endpoint=endpoint or url, method="GET", url=url, params={}, body=None, qids=list(qids), asked_at=now(),
               http_status=None, seconds=None, bytes=None, sha256=None, file=None, error=None)
    t0 = time.monotonic()
    part = f"download_{rid}_{time.time_ns()}.part"
    tmp = os.path.join(store.folder, part)
    sha, md5, n, headers = hashlib.sha256(), hashlib.md5(), 0, {}
    try:
        with client.session(group).get(url, headers={"User-Agent": USER_AGENT}, stream=True, timeout=TIMEOUT) as r:
            row["http_status"] = r.status_code
            headers = dict(r.headers)
            if r.status_code != 200:
                row.update(error=f"HTTP {r.status_code}", seconds=round(time.monotonic() - t0, 3))
                return store.put(row, r.content, "txt", headers)
            with open(tmp, "xb") as fh:  # "x": a part is never written twice
                for b in r.iter_content(chunk):
                    fh.write(b)
                    sha.update(b)
                    md5.update(b)
                    n += len(b)
    except Exception as e:  # noqa: BLE001  (the part, as far as it came, stays under its own name)
        row.update(error=f"{type(e).__name__}: {e}", seconds=round(time.monotonic() - t0, 3), bytes=n,
                   file=part if os.path.exists(tmp) else None)
        return store.put(row, None, "part", headers)
    ext = url.rsplit("/", 1)[-1].split(".", 1)[-1] if "." in url.rsplit("/", 1)[-1] else "bin"
    name = f"{sha.hexdigest()}.{ext}"
    try:
        os.link(tmp, os.path.join(store.folder, name))  # never over an existing file
    except FileExistsError:
        pass  # the same bytes were stored before under their sha256
    os.unlink(tmp)  # the part's name only: its bytes are the file stored under `name`
    headers["x-md5"] = md5.hexdigest()
    row.update(seconds=round(time.monotonic() - t0, 3), bytes=n, sha256=sha.hexdigest(), file=name)
    return store.put(row, None, ext, headers)
