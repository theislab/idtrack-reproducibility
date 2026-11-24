"""Family F's `.h5ad` files written and read with h5py alone, in AnnData's on-disk format (anndata's "Elements" of
encoding version 0.1.0 and 0.2.0: `anndata`, `dataframe`, `csr_matrix`, `array`, `string-array`, `categorical`,
`dict`, `string`, `numeric-scalar`), so that the generator and the scorer need no anndata (the package's environment
has none; the harmonizer's environment reads the files with anndata itself, and its run checks that what it read is
what the generator wrote).

    write(path, X, obs, var, layers={}, obsm={}, obsp={}, varm={}, uns={})
      X, layers, obsp: scipy CSR matrices or NumPy arrays; obs, var: {"index": [...], "columns": {name: column}} with a
      column a NumPy array (numbers, booleans) or {"categories": [...], "codes": [...]} (a categorical) or a list of
      strings; obsm, varm: NumPy arrays; uns: a dict of strings, numbers, NumPy arrays and dicts.
    read(path)      the same structure back (what the scorer compares with the harmonizer's outputs)
"""
import numpy as np


def _str_dtype():
    import h5py
    return h5py.string_dtype(encoding="utf-8")


def _attrs(node, kind, version):
    node.attrs["encoding-type"] = kind
    node.attrs["encoding-version"] = version


def _strings(group, name, values):
    d = group.create_dataset(name, data=np.array([str(v) for v in values], dtype=object), dtype=_str_dtype())
    _attrs(d, "string-array", "0.2.0")
    return d


def _packed(group, name, a):
    """A dataset, gzip-compressed when it is large (anndata reads either)."""
    a = np.asarray(a)
    kw = dict(compression="gzip", compression_opts=4) if a.ndim and a.size > 1000 else {}
    return group.create_dataset(name, data=a, **kw)


def _array(group, name, values):
    d = _packed(group, name, values)
    _attrs(d, "array", "0.2.0")
    return d


def _matrix(group, name, m):
    import scipy.sparse as sp
    if sp.issparse(m):
        m = sp.csr_matrix(m)
        g = group.create_group(name)
        _attrs(g, "csr_matrix", "0.1.0")
        g.attrs["shape"] = np.array(m.shape, dtype=np.int64)
        index = np.int32 if m.nnz < 2**31 else np.int64  # one dtype for both index arrays, as scipy requires
        _packed(g, "data", m.data)
        _packed(g, "indices", m.indices.astype(index))
        _packed(g, "indptr", m.indptr.astype(index))
        return g
    return _array(group, name, m)


def _frame(group, name, frame):
    g = group.create_group(name)
    _attrs(g, "dataframe", "0.2.0")
    g.attrs["_index"] = "_index"
    cols = list(frame.get("columns", {}))
    g.attrs.create("column-order", data=np.array(cols, dtype=object), dtype=_str_dtype())
    _strings(g, "_index", frame["index"])
    for c, v in frame.get("columns", {}).items():
        if isinstance(v, dict):
            cg = g.create_group(c)
            _attrs(cg, "categorical", "0.2.0")
            cg.attrs["ordered"] = False
            _array(cg, "codes", np.asarray(v["codes"], dtype=np.int8))
            _strings(cg, "categories", v["categories"])
        elif isinstance(v, list) or (isinstance(v, np.ndarray) and v.dtype == object):
            _strings(g, c, v)
        else:
            _array(g, c, v)
    return g


def _dict(group, name, d):
    g = group.create_group(name)
    _attrs(g, "dict", "0.1.0")
    for k, v in d.items():
        if isinstance(v, dict):
            _dict(g, k, v)
        elif isinstance(v, str):
            ds = g.create_dataset(k, data=v, dtype=_str_dtype())
            _attrs(ds, "string", "0.2.0")
        elif isinstance(v, (bool, int, float, np.integer, np.floating, np.bool_)):
            ds = g.create_dataset(k, data=v)
            _attrs(ds, "numeric-scalar", "0.2.0")
        else:
            _matrix(g, k, v)
    return g


def write(path, X, obs, var, layers=None, obsm=None, obsp=None, varm=None, uns=None):
    """Write one AnnData file (see the module docstring); written beside the path and renamed."""
    import os

    import h5py
    tmp = path + ".tmp"
    with h5py.File(tmp, "w") as f:
        _attrs(f, "anndata", "0.1.0")
        _matrix(f, "X", X)
        _frame(f, "obs", obs)
        _frame(f, "var", var)
        for name, d in (("layers", layers), ("obsm", obsm), ("obsp", obsp), ("varm", varm), ("varp", {})):
            g = f.create_group(name)
            _attrs(g, "dict", "0.1.0")
            for k, v in (d or {}).items():
                _matrix(g, k, v)
        _dict(f, "uns", uns or {})
    os.replace(tmp, path)


def _read(node):
    import h5py
    import scipy.sparse as sp
    kind = node.attrs.get("encoding-type", "")
    kind = kind.decode() if isinstance(kind, bytes) else kind
    if kind == "csr_matrix":
        return sp.csr_matrix((node["data"][()], node["indices"][()], node["indptr"][()]), shape=tuple(node.attrs["shape"]))
    if kind == "string-array":
        return [x.decode() if isinstance(x, bytes) else str(x) for x in node[()]]
    if kind == "string":
        v = node[()]
        return v.decode() if isinstance(v, bytes) else str(v)
    if kind == "numeric-scalar":
        return node[()].item()
    if kind == "categorical":
        return {"categories": _read(node["categories"]), "codes": node["codes"][()].tolist()}
    if kind == "dataframe":
        order = [x.decode() if isinstance(x, bytes) else str(x) for x in node.attrs["column-order"]]
        return {"index": _read(node[node.attrs["_index"]]), "columns": {c: _read(node[c]) for c in order}}
    if kind == "dict" or isinstance(node, h5py.Group):
        return {k: _read(node[k]) for k in node}
    return node[()]


def read(path):
    """{X, obs, var, layers, obsm, obsp, varm, uns} of a file written by `write` (or by anndata in the same encodings)."""
    import h5py
    with h5py.File(path, "r", locking=False) as f:
        return {k: _read(f[k]) for k in ("X", "obs", "var", "layers", "obsm", "obsp", "varm", "uns") if k in f}


# ---------------- the forms both sides compare (the run in the harmonizer's environment, the scorer here) ----------------

def digest(m):
    """A short content digest of an array or a sparse matrix: its shape, its dtype and its values (sparse: as CSR with
    sorted indices, the index arrays as int64, so that the index dtype a library chose does not count)."""
    import hashlib

    import scipy.sparse as sp
    h = hashlib.sha256()
    if sp.issparse(m):
        m = sp.csr_matrix(m)
        m.sort_indices()
        h.update(f"csr{m.shape}|{m.data.dtype.str}|".encode())
        for a in (m.data, m.indices.astype(np.int64), m.indptr.astype(np.int64)):
            h.update(np.ascontiguousarray(a).tobytes())
    else:
        a = np.asarray(m)
        h.update(f"dense{a.shape}|{a.dtype.str}|".encode())
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()[:24]


def column_digests(m, rows=None):
    """The digest of each column of a matrix (its values as float64, in row order; `rows` a slice of the rows): what a
    feature holds in X or a layer, compared whatever the matrix's kind and dtype (counts are exact in float32)."""
    import hashlib

    import scipy.sparse as sp
    if rows is not None:
        m = m[rows]
    if sp.issparse(m):
        m = sp.csc_matrix(m)
        n = m.shape[0]
        out = []
        for j in range(m.shape[1]):
            v = np.zeros(n, dtype=np.float64)
            s, e = m.indptr[j], m.indptr[j + 1]
            v[m.indices[s:e]] = m.data[s:e]
            out.append(hashlib.sha256(v.tobytes()).hexdigest()[:16])
        return out
    a = np.asarray(m, dtype=np.float64)
    return [hashlib.sha256(np.ascontiguousarray(a[:, j]).tobytes()).hexdigest()[:16] for j in range(a.shape[1])]


def zero_digest(n):
    import hashlib
    return hashlib.sha256(np.zeros(n, dtype=np.float64).tobytes()).hexdigest()[:16]


def as_text(v):
    """A table value as text, the same on both sides: missing values (None, NaN, pandas' NA and NaT) as "", booleans as
    True/False, integers without a point, floats by repr."""
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return "" if v != v else repr(float(v))
    try:
        import pandas as pd
        if v is pd.NA or v is pd.NaT:
            return ""
    except ImportError:  # pragma: no cover
        pass
    return str(v)


def frame_values(frame):
    """{column: [values as text]} of a frame written by `write` and read by `read` (categoricals as their labels)."""
    out = {}
    for c, v in frame.get("columns", {}).items():
        if isinstance(v, dict):
            cats = v["categories"]
            out[c] = ["" if k < 0 else as_text(cats[k]) for k in v["codes"]]
        else:
            out[c] = [as_text(x) for x in (v.tolist() if hasattr(v, "tolist") else v)]
    return out


def plain(x):
    """A value of `uns` as plain JSON data (arrays as lists, numbers as Python numbers)."""
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if hasattr(x, "tolist"):
        return x.tolist()
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x
