"""Family E's adapter: F0's adapter (`../adapter.py`) widened to every call and argument family E asks.

A line of family E is F0's line (docs/families/f0.md, "The query line") plus:

- `call`: `convert` (default; F0's reading of the answer), `convert_many` (a batch; the item `pos` is scored, or the
  call as a whole when `pos` is None), `infer_source`, `classify` (a batch converted, then classified), `api` (any
  other method of the API: `method`, `args`, `kwargs`).
- `kwargs`: arguments of the call, each an *encoded value* (below); they replace the line's typed fields
  (`from_release`, `to_release`, `final_database`, `from_database`), so a release can be `"115"`, `115.0`, `True`,
  NaN or a NumPy integer.
- `id_codec`: F0's codecs (`str`, `json`, `bytes_hex`, `float_nan`) plus `enc` (`id` is the JSON of an encoded value).
- `batch` (a list of encoded values) and `batch_container` (`list` default, `tuple`, `np_object`, `np_str`,
  `pd_series`, `pd_series_string`, `pd_index`, `pd_categorical`, `generator`, `set`, `str`: the batch argument as
  that container; `str` passes the first item alone as the batch), `pos`.
- `generations`: the api generations the line is asked on (default both; a line whose argument only one generation
  has names it). A line not asked is left out by `run_e.py` and counted in its log.

An *encoded value* is a JSON value (a string, number, boolean, null, list) or `{"$": codec, "v": value}` with codec
`float` (`"nan"`, `"inf"`, `"-inf"` or a number), `np` (`[dtype, value]`: `np.int64(115)`, `np.float64(nan)`,
`np.str_("TP53")`), `bytes_hex`, `pd_na`, `pd_nat`, `decimal`, `tuple`, `object` (a plain `object()`).

The answer row is F0's (`../adapter.py`). Added: `flags` (JSON) carries `digest` (sha256 of the canonical whole
result: `to_dict()` from 71690ab on, the result dict on d4c2dda, both as `json.dumps(sort_keys=True)`; for a
batch, of item `pos`; for `infer_source`, of the returned list), `status` of the item, `n_items`; the error classes
of a d4c2dda error are tagged `d4c2dda:<class>` in `error_mro` beside the plain names (d4c2dda raises the built-in
`TypeError` and `ValueError` where 71690ab raises its own classes that derive from them; an `expect_error` of
`MissingValueError|d4c2dda:ValueError` is then strict on both trees); an item of a batch whose conversion raised
(status `error`) is the outcome `error` with the item's error class and its classes (looked up in `idtrack.errors`).
"""
import decimal
import hashlib
import json
import math
import os
import sys
import time

from idtrack_reproducibility.families import adapter as F0
from idtrack_reproducibility.families import common as C

GENERATIONS = ("d4c2dda", "convert")
# d4c2dda's names of the infer_source modes (71690ab renamed two)
MODE_D4C2DDA = {"release": "ensembl_release", "assembly_release": "assembly_ensembl_release"}


def decode(v):
    """An encoded value (module docstring) as the Python value passed to idtrack."""
    if not isinstance(v, dict) or "$" not in v:
        return v
    c, x = v["$"], v.get("v")
    if c == "float":
        return float(x)
    if c == "np":
        import numpy as np
        dtype, val = x
        if dtype == "str_":
            return np.str_(val)
        if dtype == "bool_":
            return np.bool_(val)
        return getattr(np, dtype)(float(val) if isinstance(val, str) else val)
    if c == "bytes_hex":
        return bytes.fromhex(x)
    if c == "pd_na":
        import pandas as pd
        return pd.NA
    if c == "pd_nat":
        import pandas as pd
        return pd.NaT
    if c == "decimal":
        return decimal.Decimal(x)
    if c == "tuple":
        return tuple(decode(y) for y in x)
    if c == "object":
        return object()
    if c == "list":
        return [decode(y) for y in x]
    if c == "dict":
        return {k: decode(y) for k, y in x.items()}
    raise ValueError(f"unknown codec {c}")


def decode_id(line):
    if (line.get("id_codec") or "str") == "enc":
        return decode(json.loads(line["id"]))
    return F0.decode_id(line)


def container(items, kind):
    kind = kind or "list"
    if kind == "list":
        return list(items)
    if kind == "tuple":
        return tuple(items)
    if kind == "set":
        return set(items)
    if kind == "generator":
        return (x for x in items)
    if kind == "str":
        return items[0]
    import numpy as np
    if kind == "np_object":
        return np.array(items, dtype=object)
    if kind == "np_str":
        return np.array(items, dtype=str)
    import pandas as pd
    if kind == "pd_series":
        return pd.Series(items, dtype=object)
    if kind == "pd_series_string":
        return pd.Series(items, dtype="string")
    if kind == "pd_index":
        return pd.Index(items)
    if kind == "pd_categorical":
        return pd.Series(pd.Categorical(items))
    raise ValueError(f"unknown batch_container {kind}")


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str, allow_nan=True)


def digest(obj):
    return hashlib.sha256(canonical(obj).encode()).hexdigest()[:32]


class AdapterE(F0.Adapter):
    """F0's adapter with every call of family E (module docstring)."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        if self.gen == "convert":
            import idtrack.errors as errors
            self.errors = errors
        else:
            self.errors = None

    # ---------------- helpers ----------------

    def _fail_tag(self, fail):
        """Tag d4c2dda's error classes (`d4c2dda:<class>`) beside the plain names."""
        if fail and fail.get("outcome") == "error" and self.gen == "d4c2dda" and fail.get("error_mro"):
            fail["error_mro"] = list(fail["error_mro"]) + [f"d4c2dda:{x}" for x in fail["error_mro"]]
        return fail

    def _mro_of(self, name):
        """The classes of an error named by a batch item (idtrack.errors, else a built-in)."""
        import builtins
        cls = getattr(self.errors, name, None) if self.errors is not None else None
        cls = cls or getattr(builtins, name, None)
        if isinstance(cls, type):
            mro = [k.__name__ for k in cls.__mro__]
        else:
            mro = [name]
        if self.gen == "d4c2dda":
            mro += [f"d4c2dda:{x}" for x in mro]
        return mro

    def _kwargs(self, line, keys=("from_release", "to_release", "final_database", "from_database")):
        kw = {k: line.get(k) for k in keys if line.get(k) is not None}
        for k, v in (line.get("kwargs") or {}).items():
            kw[k] = decode(v)
        return kw

    # ---------------- one line ----------------

    def ask(self, line, explain=False):
        call = line.get("call") or "convert"
        if call == "convert" and not line.get("kwargs") and (line.get("id_codec") or "str") != "enc":
            out = super().ask(line, explain)
            for a in out:
                self._fail_tag(a)
            return out
        vs = ["best", "all"] if (line.get("strategy") or "best") == "both" else [line.get("strategy") or "best"]
        fn = getattr(self, f"_ask_{call}")
        out = fn(line, vs, explain)
        for a in out:
            a["calls"] = a.get("calls") or 1
        return out

    def _fill_convert(self, a, res, line, view):
        """Fill an answer row from one ConversionResult (71690ab on) for one view (F0's reading)."""
        db = line.get("final_database") or (line.get("kwargs") or {}).get("final_database")
        got = list(res.targets(view))
        if res.status in ("not_found", "unknown_version", "no_conversion", "error"):
            got = []
        # F0's reading (nothing for no_target or another database); a spelling of the database (E5) is the database
        if isinstance(db, str) and (res.no_target or (res.final_database is not None
                                                      and res.final_database.casefold() != db.strip().casefold())):
            got = []
        a["got"] = sorted(set(str(x) for x in got))
        a["apart"] = sorted(str(x) for x in res.alternative_locus_genes)
        a.update(status=res.status, match=res.match, graph_id=res.graph_id, tie=bool(res.tie),
                 from_release=res.from_release, fallback=res.from_release_fallback, final_database=res.final_database)
        a["flags"] = dict(no_corresponding=res.status in ("not_found", "unknown_version"),
                          no_conversion=res.status == "no_conversion", no_target=bool(res.no_target),
                          digest=digest(res.to_dict()), reason=res.reason)
        if res.status == "error" and res.error is not None:
            a.update(outcome="error", error_type=res.error.type, error_mro=self._mro_of(res.error.type),
                     error=f"{res.error.type}: {res.error.message}"[:500])

    def _fill_d4c2dda(self, a, res, line, view):
        from idtrack_reproducibility.harness import score
        rec = dict(line)
        if line.get("kwargs") and "final_database" in line["kwargs"]:
            db = decode(line["kwargs"]["final_database"])
            # a spelling of the database (E5) is the database the result names
            if isinstance(db, str) and isinstance(res.get("final_database"), str) and \
                    res["final_database"].casefold() == db.strip().casefold():
                db = res["final_database"]
            rec["final_database"] = db
        a["got"] = sorted(set(score.targets_of(res, rec)))
        if "alternative_locus_genes" in res:
            a["apart"] = sorted(str(x) for x in res["alternative_locus_genes"])
        if res.get("no_corresponding"):
            st = "unknown_version" if res.get("known_versions") else "not_found"
        elif res.get("no_conversion"):
            st = "no_conversion"
        elif res.get("no_target"):
            st = "no_target"
        elif res.get("error"):
            st = "error"
        elif res.get("tie"):
            st = "tie"
        else:
            st = "converted"
        a.update(status=st, match=res.get("match"), graph_id=res.get("graph_id"), tie=bool(res.get("tie")),
                 from_release=res.get("from_release"), fallback=res.get("from_release_fallback"),
                 final_database=res.get("final_database"))
        a["flags"] = {f: bool(res.get(f)) for f in ("no_corresponding", "no_conversion", "no_target")}
        a["flags"]["digest"] = digest(res)
        if res.get("error"):
            e = res["error"]
            a["got"] = []
            a.update(outcome="error", error_type=e.get("type"), error_mro=self._mro_of(e.get("type")),
                     error=f"{e.get('type')}: {e.get('message')}"[:500])

    def _ask_convert(self, line, views, explain):  # noqa: D401 - overrides F0's (lines with kwargs or encoded ids)
        if line.get("call") in (None, "convert") and not line.get("kwargs") and (line.get("id_codec") or "str") != "enc":
            return super()._ask_convert(line, views, explain)
        kw = self._kwargs(line)
        ident = decode_id(line)
        out = []
        if self.gen == "d4c2dda":
            for v in views:
                k2 = dict(kw)
                k2.setdefault("strategy", v)
                for x in ("budget_states", "budget_seconds", "allow_bridges"):
                    k2.pop(x, None)
                res, fail, sec = self._call(lambda: self.api.convert_identifier(ident, **k2))
                a = self._blank(line, v)
                a["seconds"] = round(sec, 6)
                if fail:
                    a.update(self._fail_tag(fail))
                else:
                    self._fill_d4c2dda(a, res, line, v)
                a["calls"] = len(views)
                out.append(a)
            return out
        kw.setdefault("budget_seconds", self.timeout)
        kw.pop("strategy", None)
        res, fail, sec = self._call(lambda: self.api.convert(ident, **kw))
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(fail)
            else:
                self._fill_convert(a, res, line, v)
            out.append(a)
        return out

    def _batch(self, line):
        return container([decode(x) for x in line["batch"]], line.get("batch_container"))

    def _ask_convert_many(self, line, views, explain):
        kw = self._kwargs(line)
        pos = line.get("pos")
        n_items = len(line["batch"])
        out = []
        if self.gen == "d4c2dda":
            for x in ("budget_states", "budget_seconds", "allow_bridges", "progress"):
                kw.pop(x, None)
            kw.setdefault("verbose", False)
            for v in views:
                k2 = dict(kw)
                k2.setdefault("strategy", v)
                batch = self._batch(line)
                res, fail, sec = self._call(lambda: self.api.convert_identifier_multiple(batch, **k2))
                a = self._blank(line, v)
                a["seconds"] = round(sec, 6)
                if fail:
                    a.update(self._fail_tag(fail))
                elif pos is None:
                    a["got"] = []
                    a["flags"] = dict(digest=digest(res), n_items=len(res))
                else:
                    self._fill_d4c2dda(a, res[pos], line, v)
                    a["flags"]["n_items"] = len(res)
                a["calls"] = len(views)
                out.append(a)
            return out
        kw.setdefault("budget_seconds", self.timeout)
        batch = self._batch(line)
        res, fail, sec = self._call(lambda: self.api.convert_many(batch, **kw))
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(fail)
            elif pos is None:
                a["got"] = []
                a["status"] = ",".join(sorted({r.status for r in res}))
                a["flags"] = dict(digest=digest([r.to_dict() for r in res]), n_items=len(res))
            else:
                self._fill_convert(a, res[pos], line, v)
                a["flags"]["n_items"] = len(res)
            out.append(a)
        _ = n_items
        return out

    def _ask_infer_source(self, line, views, explain):
        kw = {k: decode(v) for k, v in (line.get("kwargs") or {}).items()}
        batch = self._batch(line)
        if self.gen == "d4c2dda":
            if "mode" in kw and isinstance(kw["mode"], str):
                kw["mode"] = MODE_D4C2DDA.get(kw["mode"], kw["mode"])
            fn = lambda: self.api.infer_identifier_source(batch, **kw)  # noqa: E731
        else:
            fn = lambda: self.api.infer_source(batch, **kw)  # noqa: E731
        res, fail, sec = self._call(fn)
        out = []
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(self._fail_tag(fail))
            else:
                a["got"] = sorted(canonical(x) for x in res)
                a["status"] = "inferred"
                a["flags"] = dict(digest=digest(res), n=len(res))
            out.append(a)
        return out

    def _ask_classify(self, line, views, explain):
        kw = {k: decode(v) for k, v in (line.get("kwargs") or {}).items()}
        batch = self._batch(line)
        strategy = kw.pop("strategy", "best")
        out = []

        def fn():
            if self.gen == "d4c2dda":
                return self.api.classify_multiple_conversion(self.api.convert_identifier_multiple(batch, verbose=False,
                                                                                                  strategy="best"))
            return self.api.classify(self.api.convert_many(batch), strategy=strategy)

        res, fail, sec = self._call(fn)
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(self._fail_tag(fail))
            else:
                a["got"] = sorted(f"{k}:{len(x)}" for k, x in res.items())
                a["status"] = "classified"
                a["flags"] = dict(digest=digest({k: len(x) for k, x in res.items()}))
            out.append(a)
        return out

    def _ask_api(self, line, views, explain):
        method = line["method"] if self.gen == "convert" else (line.get("method_d4c2dda") or line["method"])
        args = [decode(x) for x in (line.get("args") or [])]
        kw = {k: decode(v) for k, v in (line.get("kwargs") or {}).items()}
        fn = getattr(self.api, method)
        res, fail, sec = self._call(lambda: fn(*args, **kw))
        out = []
        for v in views:
            a = self._blank(line, v)
            a["seconds"] = round(sec, 6)
            if fail:
                a.update(self._fail_tag(fail))
            else:
                plain = res.to_dict() if hasattr(res, "to_dict") else (
                    [x.to_dict() if hasattr(x, "to_dict") else x for x in res] if isinstance(res, (list, tuple)) else res)
                a["got"] = [] if res is None else sorted({canonical(x) for x in (plain if isinstance(plain, list) else [plain])})[:50]
                a["status"] = "answered"
                a["flags"] = dict(digest=digest(plain))
            out.append(a)
        return out


def asked_on(line, gen):
    g = line.get("generations")
    return not g or gen in g


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


_ = math
