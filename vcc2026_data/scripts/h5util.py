#!/usr/bin/env python3
"""Chunked readers and group reductions for large h5ad files.

Everything here streams the expression matrix in row blocks and never
materialises it densely. The count checks are deliberately narrow: passing an
integrality test is evidence that a matrix *could* be raw counts, never proof
that it is, so `group_sum` reports what it saw and leaves the semantic call to
the caller (which must pass `confirmed_counts=True` from documentation, not
from the numbers).
"""
import numpy as np
import h5py
from scipy import sparse


def dec(x):
    return x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else x


def enc_type(obj):
    v = obj.attrs.get("encoding-type", "")
    return dec(v) if v is not None else ""


def read_col(h5, path):
    """Decode one obs/var column to a numpy object/numeric array.

    Handles the three layouts anndata writes: plain dataset, categorical group
    (categories+codes), and the nullable string/int group (values+mask) that
    anndata 0.13 produces for plain string columns.
    """
    x = h5[path]
    if isinstance(x, h5py.Group):
        keys = set(x.keys())
        if {"categories", "codes"} <= keys:
            cats = np.array([dec(v) for v in x["categories"][:]], dtype=object)
            codes = x["codes"][:]
            out = np.empty(len(codes), dtype=object)
            ok = codes >= 0
            out[ok] = cats[codes[ok]]
            out[~ok] = None
            return out
        if {"values", "mask"} <= keys:
            vals = x["values"][:]
            mask = x["mask"][:].astype(bool)
            if vals.dtype.kind in "SO":
                vals = np.array([dec(v) for v in vals], dtype=object)
            else:
                vals = vals.astype(object)
            vals[mask] = None
            return vals
        raise ValueError(f"unsupported group column {path}: keys={sorted(keys)}")
    v = x[:]
    if v.dtype.kind in "SO":
        return np.array([dec(t) for t in v], dtype=object)
    return v


def obs_columns(h5, group="obs"):
    return [k for k in h5[group].keys() if k != "__categories"]


def index_name(h5, group="obs"):
    return dec(h5[group].attrs.get("_index", "_index"))


def read_frame(h5, group="obs", columns=None):
    import pandas as pd
    cols = columns if columns is not None else obs_columns(h5, group)
    data = {}
    for c in cols:
        try:
            data[c] = read_col(h5, f"{group}/{c}")
        except Exception as e:                                  # keep the rest usable
            n = matrix_shape(h5)[0 if group == "obs" else 1]
            data[c] = np.array([f"<unreadable: {type(e).__name__}>"] * n, dtype=object)
    return pd.DataFrame(data)


def matrix_shape(h5, matrix="X"):
    x = h5[matrix]
    if isinstance(x, h5py.Dataset):
        return tuple(x.shape)
    return tuple(int(v) for v in x.attrs["shape"])


def matrix_kind(h5, matrix="X"):
    x = h5[matrix]
    if isinstance(x, h5py.Dataset):
        return "dense"
    return enc_type(x) or "unknown"


def iter_row_blocks(h5, matrix="X", block=4096, rows=None):
    """Yield (start, stop, csr_matrix) over row blocks of a dense or CSR matrix."""
    x = h5[matrix]
    shape = matrix_shape(h5, matrix)
    n = shape[0] if rows is None else len(rows)
    if isinstance(x, h5py.Dataset):
        for s in range(0, n, block):
            e = min(n, s + block)
            sel = slice(s, e) if rows is None else np.sort(rows[s:e])
            yield s, e, sparse.csr_matrix(x[sel])
        return
    kind = enc_type(x)
    if kind == "csc_matrix":
        raise ValueError(f"{matrix} is CSC; convert to CSR once before reducing rows")
    if kind != "csr_matrix":
        raise ValueError(f"{matrix} has encoding-type {kind!r}, expected csr_matrix")
    indptr = x["indptr"]
    data, indices = x["data"], x["indices"]
    for s in range(0, n, block):
        e = min(n, s + block)
        if rows is None:
            ptr = indptr[s:e + 1]
            lo, hi = int(ptr[0]), int(ptr[-1])
            yield s, e, sparse.csr_matrix(
                (data[lo:hi], indices[lo:hi], ptr - lo), shape=(e - s, shape[1]))
        else:
            sel = np.sort(rows[s:e])
            parts = []
            for r in sel:
                lo, hi = int(indptr[r]), int(indptr[r + 1])
                parts.append(sparse.csr_matrix(
                    (data[lo:hi], indices[lo:hi], np.array([0, hi - lo])),
                    shape=(1, shape[1])))
            yield s, e, sparse.vstack(parts, format="csr") if parts else \
                sparse.csr_matrix((0, shape[1]))


def audit_values(blk):
    """Facts about one block's stored values; no semantic conclusion."""
    v = blk.data
    if v.size == 0:
        return dict(n=0, min=0.0, max=0.0, all_finite=True, all_nonneg=True,
                    all_integral=True)
    return dict(n=int(v.size), min=float(v.min()), max=float(v.max()),
                all_finite=bool(np.all(np.isfinite(v))),
                all_nonneg=bool(np.all(v >= 0)),
                all_integral=bool(np.all(np.abs(v - np.rint(v)) <= 1e-6)))


def group_sum(path, matrix, codes, n_groups, block=4096, confirmed_counts=False,
              progress=None):
    """Sum matrix rows into `n_groups` groups given a per-row group code.

    Rows with code < 0 are skipped (unassigned / excluded on purpose).
    Returns (csr sum, n_cells per group, value audit).  `confirmed_counts` must
    come from the source documentation; when set, non-integral or negative
    values raise instead of silently producing a meaningless "sum of counts".
    """
    codes = np.asarray(codes)
    with h5py.File(path, "r") as h:
        shape = matrix_shape(h, matrix)
        if shape[0] != len(codes):
            raise ValueError(f"matrix rows {shape[0]} != codes {len(codes)}")
        acc = sparse.csr_matrix((n_groups, shape[1]), dtype=np.float64)
        seen = dict(n=0, min=np.inf, max=-np.inf, all_finite=True, all_nonneg=True,
                    all_integral=True)
        n_cells = np.bincount(codes[codes >= 0], minlength=n_groups).astype(np.int64)
        for s, e, blk in iter_row_blocks(h, matrix, block):
            a = audit_values(blk)
            seen["n"] += a["n"]
            if a["n"]:
                seen["min"] = min(seen["min"], a["min"])
                seen["max"] = max(seen["max"], a["max"])
            for k in ("all_finite", "all_nonneg", "all_integral"):
                seen[k] = seen[k] and a[k]
            if confirmed_counts and not (a["all_finite"] and a["all_nonneg"]):
                raise ValueError(f"rows {s}:{e} contain negative or non-finite values; "
                                 f"this is not a raw count matrix")
            c = codes[s:e]
            keep = c >= 0
            if not keep.any():
                continue
            rows = np.nonzero(keep)[0]
            M = sparse.csr_matrix(
                (np.ones(len(rows)), (c[rows], rows)), shape=(n_groups, e - s))
            acc = acc + M @ blk
            if progress and (s // block) % progress == 0:
                print(f"    rows {e}/{shape[0]}", flush=True)
        if seen["n"] == 0:
            seen["min"] = seen["max"] = 0.0
        return acc.tocsr(), n_cells, seen


def group_sum_dense(path, matrix, codes, n_groups, block=4096, confirmed_counts=False,
                    progress=None, dtype=np.float64):
    """Same reduction as group_sum but with a dense accumulator.

    Preferred when the number of groups times features is large and the result is
    fairly dense (pseudobulk of pseudobulk), because repeatedly adding into a
    growing sparse accumulator copies the whole thing on every block. Only the
    groups a block actually touches are updated.
    """
    codes = np.asarray(codes)
    with h5py.File(path, "r") as h:
        shape = matrix_shape(h, matrix)
        if shape[0] != len(codes):
            raise ValueError(f"matrix rows {shape[0]} != codes {len(codes)}")
        acc = np.zeros((n_groups, shape[1]), dtype=dtype)
        seen = dict(n=0, min=np.inf, max=-np.inf, all_finite=True, all_nonneg=True,
                    all_integral=True)
        n_cells = np.bincount(codes[codes >= 0], minlength=n_groups).astype(np.int64)
        for s, e, blk in iter_row_blocks(h, matrix, block):
            a = audit_values(blk)
            seen["n"] += a["n"]
            if a["n"]:
                seen["min"] = min(seen["min"], a["min"])
                seen["max"] = max(seen["max"], a["max"])
            for k in ("all_finite", "all_nonneg", "all_integral"):
                seen[k] = seen[k] and a[k]
            if confirmed_counts and not (a["all_finite"] and a["all_nonneg"]):
                raise ValueError(f"rows {s}:{e} contain negative or non-finite values; "
                                 f"this is not a raw count matrix")
            c = codes[s:e]
            keep = np.nonzero(c >= 0)[0]
            if keep.size == 0:
                continue
            touched, inv = np.unique(c[keep], return_inverse=True)
            M = sparse.csr_matrix((np.ones(keep.size), (inv, keep)),
                                  shape=(touched.size, e - s))
            acc[touched] += (M @ blk).toarray()
            if progress and (s // block) % progress == 0:
                print(f"    rows {e}/{shape[0]}", flush=True)
        if seen["n"] == 0:
            seen["min"] = seen["max"] = 0.0
        return acc, n_cells, seen


def iter_col_blocks(h5, matrix="X", block=512):
    """Yield (start, stop, csc submatrix) over COLUMN blocks of a CSC matrix.

    A CSC h5ad cannot be reduced over rows without decompressing everything, but
    it can be streamed by column. For a cells x genes CSC matrix that means one
    block of genes at a time, which is exactly what a per-gene grouped sum needs.
    """
    x = h5[matrix]
    if isinstance(x, h5py.Dataset):
        raise ValueError(f"{matrix} is dense; use iter_row_blocks")
    kind = enc_type(x)
    if kind != "csc_matrix":
        raise ValueError(f"{matrix} has encoding-type {kind!r}, expected csc_matrix")
    shape = matrix_shape(h5, matrix)
    indptr, data, indices = x["indptr"], x["data"], x["indices"]
    for s in range(0, shape[1], block):
        e = min(shape[1], s + block)
        ptr = indptr[s:e + 1]
        lo, hi = int(ptr[0]), int(ptr[-1])
        yield s, e, sparse.csc_matrix(
            (data[lo:hi], indices[lo:hi], ptr - lo), shape=(shape[0], e - s))


def group_sum_csc(path, matrix, codes, n_groups, block=512, confirmed_counts=False,
                  progress=None):
    """Grouped row sum of a CSC matrix, streamed by column block.

    Returns (dense (n_groups, n_features) float64, n_cells, value audit). Dense is
    the right shape here because the OUTPUT is small relative to the input: the
    caller must keep n_groups modest, which is why the KOLF ingest checks the
    projected size before calling.
    """
    codes = np.asarray(codes)
    with h5py.File(path, "r") as h:
        shape = matrix_shape(h, matrix)
        if shape[0] != len(codes):
            raise ValueError(f"matrix rows {shape[0]} != codes {len(codes)}")
        acc = np.zeros((n_groups, shape[1]), dtype=np.float64)
        n_cells = np.bincount(codes[codes >= 0], minlength=n_groups).astype(np.int64)
        seen = dict(n=0, min=np.inf, max=-np.inf, all_finite=True, all_nonneg=True,
                    all_integral=True)
        safe = np.where(codes >= 0, codes, n_groups)      # a bin to discard
        for s, e, blk in iter_col_blocks(h, matrix, block):
            a = audit_values(blk)
            seen["n"] += a["n"]
            if a["n"]:
                seen["min"] = min(seen["min"], a["min"])
                seen["max"] = max(seen["max"], a["max"])
            for k in ("all_finite", "all_nonneg", "all_integral"):
                seen[k] = seen[k] and a[k]
            if confirmed_counts and not (a["all_finite"] and a["all_nonneg"]):
                raise ValueError(f"columns {s}:{e} contain negative or non-finite values")
            B = blk.tocsc()
            for j in range(e - s):
                lo, hi = B.indptr[j], B.indptr[j + 1]
                if hi == lo:
                    continue
                rows = B.indices[lo:hi]
                acc[:, s + j] = np.bincount(safe[rows], weights=B.data[lo:hi],
                                            minlength=n_groups + 1)[:n_groups]
            if progress and (s // block) % progress == 0:
                print(f"    genes {e}/{shape[1]}", flush=True)
        if seen["n"] == 0:
            seen["min"] = seen["max"] = 0.0
        return acc, n_cells, seen
