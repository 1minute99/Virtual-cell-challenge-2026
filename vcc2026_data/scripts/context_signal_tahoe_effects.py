#!/usr/bin/env python3
"""Raw-count aggregation primitives for the staged Tahoe pilot.

Groups must be frozen from metadata (sample x cell line) before expression reads.
Full-library totals are preserved even when only a fixed gene subset is retained.
This module performs no download and does not infer biological replicates.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def remove_cls(offsets, gene_tokens, counts):
    """Remove only the documented leading (token1, value-2) marker in each cell.

    All other negative/nonintegral values remain errors in CountAccumulator.
    This is not a transformation or clipping of expression measurements.
    """
    offsets = np.asarray(offsets, dtype=np.int64)
    genes, values = np.asarray(gene_tokens), np.asarray(counts)
    lo, hi = int(offsets[0]), int(offsets[-1])
    genes, values = genes[lo:hi], values[lo:hi]
    offsets = offsets - lo
    if np.any(np.diff(offsets) < 1):
        raise ValueError("Publisher expression row lacks the leading CLS marker")
    starts = offsets[:-1]
    if not (np.all(genes[starts] == 1) and np.all(values[starts] == -2)):
        raise ValueError("Unexpected CLS representation; do not silently remove expression")
    keep = np.ones(len(genes), dtype=bool)
    keep[starts] = False
    return offsets - np.arange(len(offsets)), genes[keep], values[keep]


@dataclass
class CountAccumulator:
    n_groups: int
    gene_lookup: np.ndarray
    n_genes: int
    confirmed_raw_counts: bool

    def __post_init__(self):
        if not self.confirmed_raw_counts:
            raise ValueError("Count semantics need publisher/schema evidence before raw-count aggregation")
        self.gene_lookup = np.asarray(self.gene_lookup, dtype=np.int64)
        self.sums = np.zeros((self.n_groups, self.n_genes), dtype=np.float64)
        self.full_library_sums = np.zeros(self.n_groups, dtype=np.float64)
        self.cells = np.zeros(self.n_groups, dtype=np.int64)
        self.zero_library_cells = np.zeros(self.n_groups, dtype=np.int64)
        self.total_cells_seen = 0
        self.total_cells_kept = 0

    def add_flat(self, group_ids, offsets, gene_tokens, counts):
        groups = np.asarray(group_ids, dtype=np.int64)
        offsets = np.asarray(offsets, dtype=np.int64)
        gene_tokens = np.asarray(gene_tokens, dtype=np.int64)
        counts = np.asarray(counts, dtype=np.float64)
        if offsets.shape != (len(groups) + 1,) or np.any(np.diff(offsets) < 0):
            raise ValueError("Invalid ragged offsets")
        # Sliced Arrow lists may retain offsets into the parent values array.
        lo, hi = int(offsets[0]), int(offsets[-1])
        gene_tokens, counts = gene_tokens[lo:hi], counts[lo:hi]
        offsets = offsets - lo
        if len(gene_tokens) != len(counts) or offsets[-1] != len(counts):
            raise ValueError("Gene/count ragged arrays disagree")
        if not np.isfinite(counts).all() or np.any(counts < 0):
            raise ValueError("Non-finite or negative expression is not a raw count input")
        if np.any(np.abs(counts - np.rint(counts)) > 1e-5):
            raise ValueError("Non-integral expression: revisit publisher count semantics")
        if np.any(gene_tokens < 0) or np.any(gene_tokens >= len(self.gene_lookup)):
            raise ValueError("Gene vocabulary does not cover observed tokens")
        if np.any(groups >= self.n_groups):
            raise ValueError("Group mapping exceeds frozen metadata group table")
        self.total_cells_seen += len(groups)
        kept_cells = groups >= 0
        if not np.any(kept_cells):
            return
        lengths = np.diff(offsets)
        cell_tokens = np.repeat(np.arange(len(groups)), lengths)
        libraries = np.bincount(cell_tokens, weights=counts, minlength=len(groups))
        self.full_library_sums += np.bincount(groups[kept_cells], weights=libraries[kept_cells], minlength=self.n_groups)
        self.cells += np.bincount(groups[kept_cells], minlength=self.n_groups)
        zeros = kept_cells & (libraries == 0)
        self.zero_library_cells += np.bincount(groups[zeros], minlength=self.n_groups)
        self.total_cells_kept += int(kept_cells.sum())
        unique, local_codes = np.unique(groups[kept_cells], return_inverse=True)
        local = np.full(len(groups), -1, dtype=np.int64)
        local[kept_cells] = local_codes
        token_groups = np.repeat(local, lengths)
        measured = self.gene_lookup[gene_tokens]
        keep = (token_groups >= 0) & (measured >= 0)
        index = token_groups[keep] * self.n_genes + measured[keep]
        local_sums = np.bincount(index, weights=counts[keep], minlength=len(unique) * self.n_genes)
        self.sums[unique] += local_sums.reshape(len(unique), self.n_genes)

    def logcp(self):
        cp = np.divide(self.sums * 1e4, self.full_library_sums[:, None],
                       out=np.zeros_like(self.sums), where=self.full_library_sums[:, None] > 0)
        return np.log1p(cp)

    def save(self, path, groups, genes, metadata):
        import json
        np.savez_compressed(path, sum_counts=self.sums, full_library_sums=self.full_library_sums,
            cell_counts=self.cells, zero_library_cells=self.zero_library_cells,
            groups=np.asarray(groups), genes=np.asarray(genes),
            metadata=np.asarray(json.dumps(metadata)), total_cells_seen=self.total_cells_seen,
            total_cells_kept=self.total_cells_kept)


def matched_effect(accumulator, treated_group, control_groups):
    """One sample effect against pooled controls matched externally on plate/line/time."""
    controls = np.asarray(control_groups, dtype=np.int64)
    if not len(controls):
        raise ValueError("No matched controls")
    control_sum = accumulator.sums[controls].sum(0)
    control_library = float(accumulator.full_library_sums[controls].sum())
    treatment_library = float(accumulator.full_library_sums[treated_group])
    if control_library <= 0 or treatment_library <= 0:
        raise ValueError("Empty treatment/control count library")
    return (np.log1p(accumulator.sums[treated_group] * 1e4 / treatment_library) -
            np.log1p(control_sum * 1e4 / control_library))


def check_count_primitives():
    """Small scientific checks: excluded-gene library mass, cell filtering, chunk invariance."""
    # Gene1 is outside the evaluation subset but must remain in full-library totals.
    lookup = np.array([0, -1, 1])
    groups = np.array([0, 1, -1, 2])
    offsets = np.array([0, 3, 6, 9, 9])
    genes = np.tile(np.arange(3), 3)
    counts = np.array([2, 3, 5, 5, 5, 10, 999, 999, 999], dtype=float)
    a = CountAccumulator(3, lookup, 2, True)
    a.add_flat(groups, offsets, genes, counts)
    np.testing.assert_array_equal(a.sums, [[2, 5], [5, 10], [0, 0]])
    np.testing.assert_array_equal(a.full_library_sums, [10, 20, 0])
    np.testing.assert_array_equal(a.cells, [1, 1, 1])
    np.testing.assert_array_equal(a.zero_library_cells, [0, 0, 1])
    expected = np.log1p([2500, 5000]) - np.log1p([2000, 5000])
    np.testing.assert_allclose(matched_effect(a, 1, [0]), expected, atol=1e-12)
    b = CountAccumulator(3, lookup, 2, True)
    b.add_flat(groups[:1], offsets[:2], genes, counts)
    b.add_flat(groups[1:], offsets[1:], genes, counts)
    np.testing.assert_array_equal(a.sums, b.sums)
    np.testing.assert_array_equal(a.full_library_sums, b.full_library_sums)
    of, ge, co = remove_cls(np.array([0, 3, 5]), np.array([1, 3, 5, 1, 8]),
                            np.array([-2, 4, 6, -2, 9]))
    np.testing.assert_array_equal(of, [0, 2, 3])
    np.testing.assert_array_equal(ge, [3, 5, 8])
    np.testing.assert_array_equal(co, [4, 6, 9])
    try:
        remove_cls([0, 1], [3], [-2])
    except ValueError:
        pass
    else:
        raise AssertionError("An undocumented negative value was silently removed")
    return dict(excluded_gene_library_preserved=True, unwanted_cells_excluded=True,
                zero_library_cells_recorded=True, split_batch_aggregation_identical=True,
                documented_cls_removed=True, unexpected_negative_marker_rejected=True)


if __name__ == "__main__":
    import json
    print(json.dumps(check_count_primitives(), indent=2))
