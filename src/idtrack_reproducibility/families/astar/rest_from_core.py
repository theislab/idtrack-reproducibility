"""A*'s second reading where Ensembl REST does not answer: REST `archive/id` read from the tables it serves.

REST `archive/id` (release 116) answers from the release-116 core database: an id in the `gene` table is current
(`release` 116, no `possible_replacement`); an id that is not is the old side of the events that removed it, and
`possible_replacement` lists the new ids of the session that removed it. rest.ensembl.org has not answered from
this cluster since 2026-11-18 ~15:30 while ftp.ensembl.org does; for the ids REST had not answered before
(`archive_astar.json` lacks them), this script builds the same records from the FTP MySQL dumps of the 116 core
(`stable_id_event`, `mapping_session`, `gene`; fetched into the truth store's `core/` with a SOURCES line each) and
writes them to `rest/archive_astar_core116.json` (`_meta.source` says so). `truth.Rest` reads them after the REST
records, and every line they judge names its source (`rest_source`).

Usage: rest_from_core.py <organism> --ids <file>
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
import os
import re
import sys

from idtrack_reproducibility.harness import harness_lib as H

FTP = "https://ftp.ensembl.org/pub/release-116/mysql"


def core116(org):
    html = H.http_get(FTP + "/").decode()
    dbs = sorted({d for d in re.findall(rf'href="({org}_core_116_[^"/]+?)/?"', html) if not d.endswith(".gz")})
    db = max(dbs, key=lambda d: int(re.search(r"_(\d+)[a-z]*$", d).group(1)))  # human also serves the _37 archive core
    paths = {}
    for t in ("stable_id_event", "mapping_session", "gene"):
        rp = f"core/{db}.{t}.txt.gz"
        p = os.path.join(H.truth_dir(org), rp)
        if not os.path.exists(p):
            url = f"{FTP}/{db}/{t}.txt.gz"
            data = H.http_get(url)
            with open(p + ".tmp", "wb") as fo:
                fo.write(data)
            os.replace(p + ".tmp", p)
            H.add_source(org, rp, url, len(data), hashlib.md5(data).hexdigest(),
                         datetime.datetime.now().isoformat(timespec="seconds") + " (bench A*: REST's tables at 116)")
        paths[t] = p
    return db, paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("organism")
    ap.add_argument("--ids", required=True)
    a = ap.parse_args()
    org = a.organism
    db, p = core116(org)
    cols = None
    current = {}
    with gzip.open(p["gene"], "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if cols is None:
                cols = (12, 13) if len(f) >= 14 else None
            current[f[12]] = f[13]
    sess = {}
    with gzip.open(p["mapping_session"], "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            o, n = H.rel_int(f[3]), H.rel_int(f[4])
            sess[int(f[0])] = ((n if n > o else o + 1) - 1, f[7])  # the last release before the session takes effect
    removed = collections.defaultdict(lambda: collections.defaultdict(set))  # id -> session -> new ids
    with gzip.open(p["stable_id_event"], "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[5] != "gene" or f[0] == "\\N" or f[0] == f[2]:
                continue
            removed[f[0]][int(f[4])].add(None if f[2] == "\\N" else f[2])
    ids = {x.strip() for x in open(a.ids) if x.strip()}
    out_p = os.path.join(H.truth_dir(org), "rest", "archive_astar_core116.json")
    out = json.load(open(out_p)) if os.path.exists(out_p) else {}
    out.pop("_meta", None)
    todo = set(ids)
    while todo:
        nxt = set()
        for x in sorted(todo):
            if x in out:
                continue
            if x in current:
                out[x] = dict(id=x, release="116", is_current="1", possible_replacement=[], source="core116")
                continue
            ss = removed.get(x)
            if not ss:
                out[x] = {"_error": "no record in the 116 core (neither current nor removed by an event)"}
                continue
            # the session that removed it last (REST's latest archive entry); the release before it takes effect is the
            # id's last release (make_queries.Truth.effective_release)
            s = max(ss, key=lambda k: sess.get(k, (0, ""))[0:2])
            reps = sorted(y for y in ss[s] if y)
            out[x] = dict(id=x, release=str(sess.get(s, (0,))[0]), is_current="", source="core116",
                          possible_replacement=[dict(stable_id=y) for y in reps])
            nxt |= {y for y in reps if y not in out}
        todo = nxt
    out["_meta"] = dict(source=f"{FTP}/{db} (stable_id_event, mapping_session, gene), read as REST archive/id reads them",
                        built=datetime.datetime.now().isoformat(timespec="seconds"))
    json.dump(out, open(out_p + ".tmp", "w"), sort_keys=True)
    os.replace(out_p + ".tmp", out_p)
    H.add_source(org, "rest/archive_astar_core116.json", out["_meta"]["source"], os.path.getsize(out_p),
                 hashlib.md5(open(out_p, "rb").read()).hexdigest(), out["_meta"]["built"])
    print(org, len(out) - 1, "records", collections.Counter("error" if "_error" in v else v.get("is_current") or "retired"
                                                           for k, v in out.items() if k != "_meta"))


if __name__ == "__main__":
    main()
