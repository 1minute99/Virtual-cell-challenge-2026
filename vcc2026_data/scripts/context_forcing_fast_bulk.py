#!/usr/bin/env python3
"""Exact, fused count reconstruction for context-forcing bulk diagnostics.

This helper preserves ``harness_pred.reconstructed_blocks`` sampling, NumPy
random-number order, float64 arithmetic order, clipping, and stochastic
rounding.  It only removes repeated work: the normalized control base is built
once, reconstruction diagnostics are accumulated during reconstruction, and a
strict (``fastmath=False``) Numba loop forms counts and bulk sums together.

It is intentionally standalone.  Existing scorers do not import it.
"""
from __future__ import annotations

from dataclasses import dataclass
import pathlib
import sys
from typing import Iterable, Iterator

import numpy as np
from scipy import sparse

import simple_residual_eval as E

# ``scripts/coverage.py`` is a project helper, not the third-party package that
# Numba probes optionally.  Hide this directory only while importing Numba so
# that the local module cannot be mistaken for that optional dependency.
_HERE = pathlib.Path(__file__).resolve().parent
_PATH = sys.path[:]
try:
    sys.path[:] = [p for p in sys.path if pathlib.Path(p or ".").resolve() != _HERE]
    from numba import njit
finally:
    sys.path[:] = _PATH


@njit(cache=True, fastmath=False, error_model="numpy")
def _round_and_sum_target(base, exp_row, libraries, random_values):
    """Return every rounded count plus fused sums and raw-count diagnostics."""
    n_cells, n_genes = base.shape
    counts = np.empty((n_cells, n_genes), dtype=np.float64)
    gene_sums = np.zeros(n_genes, dtype=np.float64)
    cell_sums = np.zeros(n_cells, dtype=np.float64)
    invalid = 0
    clipped = 0
    clipped_amount = 0.0
    maximum = 0.0
    for cell in range(n_cells):
        library = libraries[cell]
        cell_sum = 0.0
        for gene in range(n_genes):
            # Keep the operation order of harness_pred.reconstructed_blocks.
            raw = (base[cell, gene] * exp_row[gene] - 1.0) * library / 1.0e4
            if not np.isfinite(raw):
                invalid += 1
                counts[cell, gene] = np.nan
                continue
            if raw < 0.0:
                clipped += 1
                clipped_amount += -raw
                wanted = 0.0
            else:
                wanted = raw
            floor = np.floor(wanted)
            value = floor + (1.0 if random_values[cell, gene] < wanted - floor else 0.0)
            counts[cell, gene] = value
            gene_sums[gene] += value
            cell_sum += value
            if value > maximum:
                maximum = value
        cell_sums[cell] = cell_sum
    return counts, gene_sums, cell_sums, invalid, clipped, clipped_amount, maximum


@dataclass
class _Prepared:
    changed: np.ndarray
    unchanged: np.ndarray
    exponentiated_delta: np.ndarray
    base: np.ndarray
    libraries: np.ndarray
    unchanged_coo: sparse.coo_matrix
    rng: np.random.Generator


def _prepare(controls: sparse.csr_matrix, delta: np.ndarray, cells: int, seed: int) -> _Prepared:
    libraries = np.asarray(controls.sum(1)).ravel().astype(np.float64)
    if libraries.min() <= 0:
        raise ValueError("a control cell has zero counts")
    changed = np.nonzero(np.abs(delta).max(0) > 0)[0]
    unchanged = np.setdiff1d(np.arange(controls.shape[1]), changed, assume_unique=True)
    # Advanced column indexing produces Fortran layout; each kernel consumes a
    # target row, so make rows contiguous without changing any values.
    exponentiated_delta = np.ascontiguousarray(
        np.exp(delta[:, changed].astype(np.float64))
    )
    rng = np.random.default_rng(seed)
    pick = rng.choice(controls.shape[0], size=cells, replace=False)
    block = controls[pick]
    selected_libraries = libraries[pick]
    dense = np.asarray(block[:, changed].todense(), dtype=np.float64)
    # This is the target-independent left side of the original expression.
    base = dense / selected_libraries[:, None] * 1.0e4 + 1.0
    unchanged_coo = block[:, unchanged].tocoo()
    return _Prepared(
        changed=changed,
        unchanged=unchanged,
        exponentiated_delta=exponentiated_delta,
        base=base,
        libraries=selected_libraries,
        unchanged_coo=unchanged_coo,
        rng=rng,
    )


def _blocks_with_diagnostics(prepared: _Prepared) -> Iterator[tuple]:
    for target_index, exp_row in enumerate(prepared.exponentiated_delta):
        random_values = prepared.rng.random(prepared.base.shape)
        result = _round_and_sum_target(
            prepared.base, exp_row, prepared.libraries, random_values
        )
        yield target_index, *result


def reconstructed_blocks(
    controls: sparse.csr_matrix,
    delta: np.ndarray,
    cells: int = 400,
    seed: int = 0,
) -> Iterator[tuple[int, np.ndarray, sparse.coo_matrix, np.ndarray, np.ndarray]]:
    """Yield the original block interface with count-for-count parity."""
    prepared = _prepare(controls, delta, cells, seed)
    for target_index, counts, *_ in _blocks_with_diagnostics(prepared):
        yield (
            target_index,
            counts,
            prepared.unchanged_coo,
            prepared.changed,
            prepared.unchanged,
        )


def reconstructed_bulk(
    controls: sparse.csr_matrix,
    delta: np.ndarray,
    targets: Iterable[str],
    cells: int,
    seed: int,
):
    """Return the same labels, normalized bulk, and diagnostics as the frozen helper."""
    targets = list(targets)
    if len(targets) != len(delta):
        raise ValueError("target count differs from delta rows")
    prepared = _prepare(controls, delta, cells, seed)
    sums = np.zeros((len(targets) + 1, controls.shape[1]), dtype=np.float64)

    unchanged_data = prepared.unchanged_coo.data.astype(np.float64)
    unchanged_gene_sums = np.bincount(
        prepared.unchanged_coo.col,
        weights=unchanged_data,
        minlength=len(prepared.unchanged),
    )
    unchanged_cell_sums = np.bincount(
        prepared.unchanged_coo.row,
        weights=unchanged_data,
        minlength=cells,
    )
    unchanged_invalid = int(np.count_nonzero(~np.isfinite(unchanged_data)))
    unchanged_negative = int(np.count_nonzero(unchanged_data < 0))
    unchanged_noninteger = int(
        np.count_nonzero(np.abs(unchanged_data - np.rint(unchanged_data)) > 0)
    )
    unchanged_maximum = float(np.nanmax(unchanged_data)) if unchanged_data.size else 0.0

    invalid = negative = noninteger = 0
    clipped = 0
    clipped_amount = 0.0
    target_totals = []
    minimum_cell_total = np.inf
    maximum_cell_total = 0.0
    maximum_count = unchanged_maximum

    for (
        target_index,
        counts,
        changed_gene_sums,
        changed_cell_sums,
        target_invalid,
        target_clipped,
        target_clipped_amount,
        target_maximum,
    ) in _blocks_with_diagnostics(prepared):
        invalid += int(target_invalid) + unchanged_invalid
        negative += unchanged_negative
        noninteger += unchanged_noninteger
        clipped += int(target_clipped)
        clipped_amount += float(target_clipped_amount)
        maximum_count = max(maximum_count, float(target_maximum))
        sums[target_index, prepared.changed] = changed_gene_sums
        sums[target_index, prepared.unchanged] = unchanged_gene_sums
        cell_totals = changed_cell_sums + unchanged_cell_sums
        minimum_cell_total = min(minimum_cell_total, float(cell_totals.min(initial=np.inf)))
        maximum_cell_total = max(maximum_cell_total, float(cell_totals.max(initial=0)))
        target_totals.append(float(cell_totals.sum()))
        del counts

    if invalid:
        raise ValueError(f"reconstruction produced {invalid} non-finite count values")
    if negative or noninteger:
        raise ValueError(
            "invalid reconstructed counts: "
            f"nonfinite={invalid}, negative={negative}, noninteger={noninteger}"
        )
    if maximum_count > 2**24:
        raise ValueError("reconstruction has integer counts outside exact float32 range")

    sums[-1] = np.asarray(controls.astype(np.float64).sum(axis=0)).ravel()
    labels = np.asarray([*targets, E.CONTROL])
    order = np.argsort(labels)
    diagnostics = {
        "invalid_value_count": 0,
        "clipped_value_count": clipped,
        "clipped_count_amount": clipped_amount,
        "changed_gene_count": int(len(prepared.changed)),
        "noninteger_value_count": 0,
        "negative_value_count": 0,
        "predicted_count_total": float(sum(target_totals)),
        "per_target_count_total": dict(zip(targets, target_totals)),
        "control_count_total": float(sums[-1].sum()),
        "minimum_predicted_cell_total": float(minimum_cell_total),
        "maximum_predicted_cell_total": float(maximum_cell_total),
    }
    return labels[order], E.bulk_lognorm_means(sums[order], E.BULK_TARGET_SUM), diagnostics
