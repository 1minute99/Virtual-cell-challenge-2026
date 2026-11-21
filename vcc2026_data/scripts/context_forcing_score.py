#!/usr/bin/env python3
"""Persistent installed-PDS worker for the context-forcing follow-up.

Requests contain frozen, full official-axis deltas.  Reconstruction delegates to
``simple_residual_eval._reconstructed_bulk`` and scoring delegates to
``simple_residual_eval._score_bulk``; this module defines no replacement metric.
"""
from __future__ import annotations

import argparse
import contextlib
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import inspect
import json
import os
import pathlib
import sys
import time
from typing import Iterable, Mapping

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_key, "2")

import anndata as ad
import numpy as np
from scipy import sparse

import simple_residual_eval as E
import vcclib as V
from cell_eval2.metrics import discrimination as PDS


ROOT = pathlib.Path(__file__).resolve().parent.parent
REPORT = ROOT / "reports" / "context_forcing_2026-09-16"
RN_REPORT = ROOT / "reports" / "simple_residual_2026-09-16"
RN_PARTITIONS = RN_REPORT / "independent_controls_1944"
RN_MOCK = RN_REPORT / "mock"
KOLF_REPORT = RN_REPORT / "kolf_full_axis_2026-09-16"
KOLF_REAL = ROOT / "harness" / "kolf_real.h5ad"
N_FEATURES = V.N_FEATURES
OFFICIAL_CELL_LIMIT = 1_000_000
RN_LINES = ("K562", "RPE1", "HepG2", "Jurkat")
LINE_ALIASES = {
    "k562": "K562",
    "rpe1": "RPE1",
    "hepg2": "HepG2",
    "jurkat": "Jurkat",
    "kolf": "KOLF",
    "kolf2.1j": "KOLF",
    "kolf2_1j": "KOLF",
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def sparse_sha256(value: sparse.spmatrix) -> str:
    matrix = value.tocsr()
    digest = hashlib.sha256()
    digest.update(np.asarray(matrix.shape, dtype=np.int64).tobytes())
    for array in (matrix.indptr, matrix.indices, matrix.data):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(np.asarray(contiguous.shape, dtype=np.int64).tobytes())
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _implementation_hashes_cached() -> tuple[tuple[str, str], ...]:
    return tuple(sorted({
        "context_forcing_score": sha256(pathlib.Path(__file__)),
        "simple_residual_eval": sha256(pathlib.Path(E.__file__)),
        "installed_discrimination": sha256(pathlib.Path(inspect.getfile(PDS.discrimination_score))),
    }.items()))


def implementation_hashes() -> dict[str, str]:
    return dict(_implementation_hashes_cached())


@lru_cache(maxsize=1)
def _gene_map() -> V.GeneMap:
    return V.GeneMap()


def _write_new_json(path: pathlib.Path, value: Mapping) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temporary, path)


def canonical_line(line: str) -> str:
    key = str(line).strip().lower()
    if key not in LINE_ALIASES:
        raise ValueError(f"unknown line {line!r}; expected RN4 or KOLF")
    return LINE_ALIASES[key]


def random_derangements(n: int, count: int = 32) -> list[list[int]]:
    """Fixed seed-0..count-1 derangements, independent of reconstruction seed."""
    if count < 0:
        raise ValueError("null_repeats must be nonnegative")
    if count and n < 2:
        raise ValueError("random nulls require at least two targets")
    identity = np.arange(n)
    result: list[list[int]] = []
    for seed in range(count):
        rng = np.random.default_rng(seed)
        while True:
            order = rng.permutation(n)
            if np.all(order != identity):
                result.append([int(x) for x in order])
                break
    return result


@dataclass
class ControlData:
    line: str
    targets: tuple[str, ...]
    genes: tuple[str, ...]
    official_indices: np.ndarray
    query_controls: sparse.csr_matrix
    reference_controls: sparse.csr_matrix
    provenance: dict
    input_hashes: dict


@dataclass
class LineData(ControlData):
    real_bulk: tuple[np.ndarray, np.ndarray]

    def controls_only(self) -> ControlData:
        fields = {
            name: getattr(self, name)
            for name in ControlData.__dataclass_fields__
        }
        return ControlData(**fields)


def _validate_control_data(data: ControlData) -> None:
    if data.targets != tuple(sorted(data.targets)) or len(data.targets) != len(set(data.targets)):
        raise ValueError(f"{data.line}: targets must be sorted and unique")
    if len(data.genes) != len(set(data.genes)):
        raise ValueError(f"{data.line}: genes must be unique")
    indices = np.asarray(data.official_indices, dtype=np.int64)
    if indices.shape != (len(data.genes),):
        raise ValueError(f"{data.line}: official index shape mismatch")
    if np.any(indices < 0) or np.any(indices >= N_FEATURES) or len(np.unique(indices)) != len(indices):
        raise ValueError(f"{data.line}: invalid official feature indices")
    for name, matrix in (("query", data.query_controls), ("reference", data.reference_controls)):
        if not sparse.isspmatrix_csr(matrix):
            raise ValueError(f"{data.line}: {name} controls must be CSR")
        if matrix.shape[1] != len(data.genes):
            raise ValueError(f"{data.line}: {name} control gene dimension mismatch")
        if matrix.shape[0] == 0:
            raise ValueError(f"{data.line}: {name} control split is empty")
        values = matrix.data
        if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values != np.rint(values)):
            raise ValueError(f"{data.line}: {name} controls are not finite nonnegative counts")


def _control_profile(controls: sparse.csr_matrix) -> np.ndarray:
    sums = np.asarray(controls.astype(np.float64).sum(axis=0)).ravel()
    if not np.all(np.isfinite(sums)) or np.any(sums < 0) or sums.sum() <= 0:
        raise ValueError("control split has an invalid native count denominator")
    return E.bulk_lognorm_means(sums[None, :], E.BULK_TARGET_SUM)[0]


def _subset_truth(
    real_bulk: tuple[np.ndarray, np.ndarray],
    targets: list[str],
    reference_controls: sparse.csr_matrix,
) -> tuple[np.ndarray, np.ndarray]:
    labels = [str(x) for x in real_bulk[0]]
    if len(labels) != len(set(labels)):
        raise ValueError("real pseudobulk labels are not unique")
    positions = {label: i for i, label in enumerate(labels)}
    wanted = [*targets, E.CONTROL]
    missing = sorted(set(wanted) - set(positions))
    if missing:
        raise ValueError(f"real pseudobulk is missing labels: {missing[:5]}")
    values = np.asarray(real_bulk[1], dtype=np.float64)[[positions[x] for x in wanted]].copy()
    values[-1] = _control_profile(reference_controls)
    order = np.argsort(np.asarray(wanted))
    return np.asarray(wanted)[order], values[order]


def _score(pred_bulk, real_bulk, genes: tuple[str, ...]) -> tuple[float, dict[str, float]]:
    per_target = E._score_bulk(pred_bulk, real_bulk, list(genes))
    values = {str(key): float(value) for key, value in per_target.items()}
    return float(np.mean(list(values.values()))), values


def _permuted_groups(
    pred_bulk: tuple[np.ndarray, np.ndarray], targets: list[str], order: Iterable[int]
) -> tuple[np.ndarray, np.ndarray]:
    labels, values = pred_bulk
    labels = np.asarray(labels)
    values = np.asarray(values)
    positions = {str(label): i for i, label in enumerate(labels)}
    target_rows = np.asarray([positions[target] for target in targets], dtype=np.int64)
    permutation = np.asarray(list(order), dtype=np.int64)
    result = values.copy()
    result[target_rows] = values[target_rows[permutation]]
    return labels.copy(), result


def _relabel_score_matrix(
    pred_bulk: tuple[np.ndarray, np.ndarray],
    real_bulk: tuple[np.ndarray, np.ndarray],
    genes: tuple[str, ...],
    expected_identity_per_target: Mapping[str, float],
) -> dict:
    """Compute all relabel scores from the installed scorer's own primitives.

    The expensive cosine distance matrix is independent of the label derangement.
    This follows the installed ``discrimination_score`` panel/cosine branch and
    verifies its diagonal against a direct ``E._score_bulk`` call before use.
    """
    real_perts = np.asarray(real_bulk[0]).astype(str)
    pred_perts = np.asarray(pred_bulk[0]).astype(str)
    real_means = np.asarray(real_bulk[1], dtype=np.float64)
    pred_means = np.asarray(pred_bulk[1], dtype=np.float64)
    real_keys, real_effects = PDS.delta(real_means, real_perts, E.CONTROL)
    control_hits = np.flatnonzero(real_perts == E.CONTROL)
    if control_hits.size != 1:
        raise ValueError("real bulk must contain exactly one control row")
    predicted_mask = pred_perts != E.CONTROL
    pred_keys = pred_perts[predicted_mask]
    pred_effects = pred_means[predicted_mask] - real_means[control_hits[0]]
    real_keys = np.asarray(real_keys).astype(str)
    pred_keys = np.asarray(pred_keys).astype(str)
    if set(pred_keys) != set(real_keys):
        raise ValueError("fast relabel matrix requires matching prediction/reference panels")
    panel_columns = PDS.resolve_panel_columns(real_keys, np.asarray(genes, dtype=str))
    pred_reduced = PDS.panel_reduced(pred_effects, panel_columns)
    real_reduced = PDS.panel_reduced(real_effects, panel_columns)
    distances = PDS.pairwise_full(pred_reduced, real_reduced, "cosine")
    denominator = len(real_keys) - 1
    scores = np.empty(distances.shape, dtype=np.float64)
    for source_row in range(distances.shape[0]):
        for destination_column in range(distances.shape[1]):
            rank = PDS._match_rank(distances[source_row], destination_column, "midrank")
            scores[source_row, destination_column] = (
                1.0 if denominator <= 0 else 1.0 - rank / denominator
            )
    pred_position = {target: i for i, target in enumerate(pred_keys)}
    real_position = {target: i for i, target in enumerate(real_keys)}
    identity = {
        target: float(scores[pred_position[target], real_position[target]])
        for target in expected_identity_per_target
    }
    maximum_error = max(
        abs(identity[target] - float(expected_identity_per_target[target]))
        for target in identity
    )
    if maximum_error > 1e-12:
        raise ValueError(f"installed relabel-matrix identity parity failed: {maximum_error}")
    return {
        "pred_keys": pred_keys.tolist(),
        "real_keys": real_keys.tolist(),
        "scores": scores,
        "identity_per_target": identity,
        "identity_pds": float(np.mean(list(identity.values()))),
        "maximum_identity_per_target_error": float(maximum_error),
        "implementation": "installed discrimination.delta/resolve_panel_columns/panel_reduced/pairwise_full/_match_rank",
    }


def _score_relabel_order(
    matrix: Mapping, targets: list[str], order: Iterable[int]
) -> tuple[float, dict[str, float]]:
    pred_position = {target: i for i, target in enumerate(matrix["pred_keys"])}
    real_position = {target: i for i, target in enumerate(matrix["real_keys"])}
    scores = np.asarray(matrix["scores"], dtype=np.float64)
    permutation = list(order)
    per_target = {
        destination: float(
            scores[pred_position[targets[source]], real_position[destination]]
        )
        for destination, source in zip(targets, permutation)
    }
    return float(np.mean(list(per_target.values()))), per_target


def prepare_delta(
    delta: np.ndarray,
    input_targets: Iterable[str],
    data: ControlData,
    scale: float,
) -> tuple[np.ndarray, dict]:
    if not np.isfinite(scale):
        raise ValueError("scale must be finite")
    raw_targets = [str(x) for x in input_targets]
    gene_map = _gene_map()
    canonical = [gene_map.canon(value) for value in raw_targets]
    supported_by_canonical: dict[str, str] = {}
    for evaluation_label in data.targets:
        key = gene_map.canon(evaluation_label)
        previous = supported_by_canonical.get(key)
        if previous is not None and previous != evaluation_label:
            raise ValueError(
                f"{data.line}: frozen target labels collide after canonicalization: "
                f"{previous!r}, {evaluation_label!r}"
            )
        supported_by_canonical[key] = evaluation_label
    resolved = [
        raw if raw in data.targets else supported_by_canonical.get(canon)
        for raw, canon in zip(raw_targets, canonical)
    ]
    unsupported = [raw for raw, value in zip(raw_targets, resolved) if value is None]
    if unsupported:
        raise ValueError(f"targets are unavailable for {data.line}: {unsupported[:5]}")
    resolved = [str(value) for value in resolved]
    if len(resolved) != len(set(resolved)):
        raise ValueError("target aliases collide onto the same frozen evaluation label")
    if len(resolved) < 2:
        raise ValueError("at least two targets are required")
    value = np.asarray(delta)
    expected = (len(resolved), N_FEATURES)
    if value.shape != expected:
        raise ValueError(f"delta shape {value.shape} is invalid; expected {expected}")
    if not np.all(np.isfinite(value)):
        raise ValueError("delta contains non-finite values")
    input_position = {target: i for i, target in enumerate(resolved)}
    sorted_targets = sorted(resolved)
    rows = np.asarray([input_position[target] for target in sorted_targets], dtype=np.int64)
    local = np.asarray(value[rows][:, data.official_indices], dtype=np.float32)
    local *= float(scale)
    aliases = [
        {"input": raw, "canonical": canon, "evaluation_label": evaluation}
        for raw, canon, evaluation in zip(raw_targets, canonical, resolved)
        if raw != canon or raw != evaluation
    ]
    return local, {
        "input_targets": raw_targets,
        "canonical_input_targets": canonical,
        "resolved_source_targets": resolved,
        "sorted_targets": sorted_targets,
        "sorted_from_input_rows": rows.tolist(),
        "aliases": aliases,
    }


class ContextForcingScoreWorker:
    """Cache source material once per line and serve score/export requests."""

    def __init__(self, report_root: pathlib.Path = REPORT, cache_root: pathlib.Path | None = None):
        self.report_root = pathlib.Path(report_root).resolve()
        self.cache_root = pathlib.Path(cache_root or self.report_root / "score_cache").resolve()
        self._ensure_under_report(self.cache_root)
        self._control_cache: dict[str, ControlData] = {}
        self._line_cache: dict[str, LineData] = {}

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        self.close()

    def close(self) -> None:
        self._control_cache.clear()
        self._line_cache.clear()

    def _ensure_under_report(self, path: pathlib.Path) -> pathlib.Path:
        resolved = pathlib.Path(path).resolve()
        try:
            resolved.relative_to(self.report_root)
        except ValueError as error:
            raise ValueError(f"new outputs must stay under {self.report_root}: {resolved}") from error
        return resolved

    def _load_controls(self, line: str) -> ControlData:
        line = canonical_line(line)
        if line in self._control_cache:
            return self._control_cache[line]
        if line in self._line_cache:
            data = self._line_cache[line].controls_only()
            self._control_cache[line] = data
            return data
        if line in RN_LINES:
            folder = RN_MOCK / line / "panel300"
            manifest_path = folder / "manifest.json"
            controls_path = folder / "controls.npz"
            partition_path = RN_PARTITIONS / f"{line}.partition.json"
            manifest = json.loads(manifest_path.read_text())
            partition = json.loads(partition_path.read_text())
            controls = sparse.load_npz(controls_path).tocsr()
            query_indices = np.asarray(partition["query_indices"], dtype=np.int64)
            reference_indices = np.asarray(partition["reference_indices"], dtype=np.int64)
            if np.intersect1d(query_indices, reference_indices).size:
                raise ValueError(f"{line}: query/reference controls overlap")
            data = ControlData(
                line=line,
                targets=tuple(str(x) for x in manifest["targets"]),
                genes=tuple(str(x) for x in manifest["genes"]),
                official_indices=np.asarray(manifest["official_feature_indices"], dtype=np.int64),
                query_controls=controls[query_indices].tocsr().copy(),
                reference_controls=controls[reference_indices].tocsr().copy(),
                provenance={
                    "kind": "RN frozen disjoint controls",
                    "panel": str(folder.resolve()),
                    "partition": str(partition_path.resolve()),
                    "query_rows": int(len(query_indices)),
                    "reference_rows": int(len(reference_indices)),
                    "row_disjoint": True,
                },
                input_hashes={
                    "manifest": sha256(manifest_path),
                    "controls": sha256(controls_path),
                    "partition": sha256(partition_path),
                },
            )
        else:
            axis_path = KOLF_REPORT / "axis.json"
            partition_path = KOLF_REPORT / "partition.json"
            axis = json.loads(axis_path.read_text())
            partition = json.loads(partition_path.read_text())
            query_rows = np.asarray(partition["query_rows"], dtype=np.int64)
            reference_rows = np.asarray(partition["reference_rows"], dtype=np.int64)
            if np.intersect1d(query_rows, reference_rows).size:
                raise ValueError("KOLF: query/reference controls overlap")
            backed = ad.read_h5ad(KOLF_REAL, backed="r")
            try:
                genes = tuple(str(x) for x in backed.var_names)
                if genes != tuple(str(x) for x in axis["genes"]):
                    raise ValueError("KOLF gene axis differs from frozen axis")
                query = sparse.csr_matrix(backed.X[query_rows]).copy()
                reference = sparse.csr_matrix(backed.X[reference_rows]).copy()
            finally:
                backed.file.close()
            data = ControlData(
                line="KOLF",
                targets=tuple(str(x) for x in axis["targets"]),
                genes=genes,
                official_indices=np.asarray(axis["official_indices"], dtype=np.int64),
                query_controls=query,
                reference_controls=reference,
                provenance={
                    "kind": "KOLF frozen disjoint controls",
                    "source": str(KOLF_REAL.resolve()),
                    "axis": str(axis_path.resolve()),
                    "partition": str(partition_path.resolve()),
                    "query_rows": int(len(query_rows)),
                    "reference_rows": int(len(reference_rows)),
                    "row_disjoint": True,
                    "export_only_reads_control_rows": True,
                },
                input_hashes={
                    "axis": sha256(axis_path),
                    "partition": sha256(partition_path),
                    "query_counts_csr": sparse_sha256(query),
                    "reference_counts_csr": sparse_sha256(reference),
                },
            )
        _validate_control_data(data)
        self._control_cache[line] = data
        return data

    def _load_line(self, line: str) -> LineData:
        line = canonical_line(line)
        if line in self._line_cache:
            return self._line_cache[line]
        controls = self._load_controls(line)
        real_path = KOLF_REAL if line == "KOLF" else RN_MOCK / line / "panel300" / "real.h5ad"
        real = ad.read_h5ad(real_path)
        try:
            if tuple(str(x) for x in real.var_names) != controls.genes:
                raise ValueError(f"{line}: real gene axis differs from frozen control axis")
            observed = sorted(set(str(x) for x in real.obs["target"].astype(str)) - {E.CONTROL})
            if observed != list(controls.targets):
                raise ValueError(f"{line}: real target labels differ from frozen target axis")
            real_bulk = E.pseudobulk_bulk_lognorm(
                real, "target", bulk_target_sum=E.BULK_TARGET_SUM
            )
        finally:
            del real
        data = LineData(
            **{name: getattr(controls, name) for name in ControlData.__dataclass_fields__},
            real_bulk=(np.asarray(real_bulk[0]), np.asarray(real_bulk[1], dtype=np.float64)),
        )
        data.input_hashes = {**data.input_hashes, "real": sha256(real_path)}
        self._line_cache[line] = data
        return data

    def export_query_features(
        self, line: str, output_path: pathlib.Path | None = None
    ) -> dict:
        line = canonical_line(line)
        output = self._ensure_under_report(
            pathlib.Path(output_path or self.cache_root / "query_features" / f"{line}.query.npz")
        )
        sidecar = output.with_suffix(output.suffix + ".json")
        data = self._load_controls(line)
        if output.exists() or sidecar.exists():
            if not output.is_file() or not sidecar.is_file():
                raise FileExistsError(f"incomplete existing query export: {output}")
            receipt = json.loads(sidecar.read_text())
            if (
                receipt.get("npz_sha256") != sha256(output)
                or receipt.get("line") != line
                or receipt.get("input_hashes") != data.input_hashes
                or receipt.get("implementation_hashes") != implementation_hashes()
            ):
                raise FileExistsError(f"existing query export identity mismatch: {output}")
            return receipt
        local_sums = np.asarray(data.query_controls.astype(np.float64).sum(axis=0)).ravel()
        denominator = float(local_sums.sum())
        if denominator <= 0 or not np.isfinite(denominator):
            raise ValueError(f"{line}: invalid query-control native denominator")
        full_sums = np.zeros(N_FEATURES, dtype=np.float64)
        full_sums[data.official_indices] = local_sums
        available = np.zeros(N_FEATURES, dtype=bool)
        available[data.official_indices] = True
        profile = np.zeros(N_FEATURES, dtype=np.float64)
        profile[data.official_indices] = np.log1p(local_sums / denominator * 1e4)
        official_genes = np.asarray(V.official_features(), dtype=str)
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_name(output.name + ".tmp")
        with temporary.open("wb") as handle:
            np.savez_compressed(
                handle,
                control_log1p_cp10k=profile,
                available_mask=available,
                native_sum_counts=full_sums,
                native_denominator=np.asarray(denominator, dtype=np.float64),
                official_genes=official_genes,
            )
        os.replace(temporary, output)
        receipt = {
            "format_version": 1,
            "status": "complete",
            "line": line,
            "output_path": str(output),
            "npz_sha256": sha256(output),
            "query_cells": int(data.query_controls.shape[0]),
            "native_measured_genes": int(len(data.genes)),
            "official_features": N_FEATURES,
            "native_denominator": denominator,
            "available_coordinates": int(available.sum()),
            "control_log1p_cp10k_array_sha256": array_sha256(profile),
            "available_mask_array_sha256": array_sha256(available),
            "native_sum_counts_array_sha256": array_sha256(full_sums),
            "normalization": "log1p(CP10000) using all available native measured-axis query-control counts",
            "response_rows_consulted": 0,
            "source_only": True,
            "provenance": data.provenance,
            "input_hashes": data.input_hashes,
            "implementation_hashes": implementation_hashes(),
        }
        _write_new_json(sidecar, receipt)
        return receipt

    def score(self, request: Mapping) -> dict:
        started = time.monotonic()
        line = canonical_line(str(request["line"]))
        delta_path = pathlib.Path(str(request["delta_path"])).resolve()
        output_path = self._ensure_under_report(pathlib.Path(str(request["output_path"])))
        cells = int(request.get("cells", 400))
        seed = int(request.get("seed", 0))
        scale = float(request.get("scale", 3.0))
        null_repeats = int(request.get("null_repeats", 32))
        if cells <= 0:
            raise ValueError("cells must be positive")
        if null_repeats < 0:
            raise ValueError("null_repeats must be nonnegative")
        if not delta_path.is_file():
            raise FileNotFoundError(delta_path)
        data = self._load_line(line)
        loaded_at = time.monotonic()
        if cells > data.query_controls.shape[0]:
            raise ValueError(
                f"cells must be in 1..{data.query_controls.shape[0]} for {line}; got {cells}"
            )
        delta = np.load(delta_path, allow_pickle=False)
        targets = request.get("targets")
        input_targets = list(data.targets) if targets is None else [str(x) for x in targets]
        scaled_local, target_audit = prepare_delta(delta, input_targets, data, scale)
        sorted_targets = target_audit["sorted_targets"]
        normalized_request = {
            "command": "score",
            "line": line,
            "delta_path": str(delta_path),
            "output_path": str(output_path),
            "targets": target_audit["input_targets"],
            "cells": cells,
            "seed": seed,
            "scale": scale,
            "null_repeats": null_repeats,
        }
        delta_file_hash = sha256(delta_path)
        delta_array_hash = array_sha256(delta)
        implementations = implementation_hashes()
        prepared_at = time.monotonic()
        if output_path.exists():
            saved = json.loads(output_path.read_text())
            same = (
                saved.get("request") == normalized_request
                and saved.get("delta", {}).get("file_sha256") == delta_file_hash
                and saved.get("delta", {}).get("array_sha256") == delta_array_hash
                and saved.get("source", {}).get("input_hashes") == data.input_hashes
                and saved.get("implementation_hashes") == implementations
            )
            if not same:
                raise FileExistsError(f"existing score output identity mismatch: {output_path}")
            return saved
        pred_labels, pred_values, diagnostics = E._reconstructed_bulk(
            data.query_controls, scaled_local, sorted_targets, cells, seed
        )
        reconstructed_at = time.monotonic()
        pred_bulk = (np.asarray(pred_labels), np.asarray(pred_values, dtype=np.float64))
        truth = _subset_truth(data.real_bulk, sorted_targets, data.reference_controls)
        pds, per_target = _score(pred_bulk, truth, data.genes)
        identity_scored_at = time.monotonic()

        orders = random_derangements(len(sorted_targets), null_repeats)
        null_values: list[float] = []
        null_per_target: list[dict] = []
        relabel_matrix = None
        if null_repeats:
            relabel_matrix = _relabel_score_matrix(pred_bulk, truth, data.genes, per_target)
            parity_pds = float(relabel_matrix["identity_pds"])
            parity_error = max(
                abs(pds - parity_pds),
                float(relabel_matrix["maximum_identity_per_target_error"]),
            )
            parity_method = (
                "direct installed _score_bulk versus diagonal of one installed-primitives "
                "relabel score matrix"
            )
        else:
            parity_pds = float(np.mean(list(per_target.values())))
            parity_error = abs(pds - parity_pds)
            parity_method = "reported PDS versus mean of installed _score_bulk per-target values"
        if parity_error > 1e-12:
            raise ValueError(f"identity PDS parity failed: error={parity_error}")
        matrix_built_at = time.monotonic()
        for null_seed, order in enumerate(orders):
            null_pds, null_per = _score_relabel_order(
                relabel_matrix, sorted_targets, order
            )
            null_values.append(null_pds)
            null_per_target.append(
                {
                    "seed": null_seed,
                    "source_target_by_destination": {
                        target: sorted_targets[source]
                        for target, source in zip(sorted_targets, order)
                    },
                    "per_target": null_per,
                }
            )
        nulls_derived_at = time.monotonic()
        null_array = np.asarray(null_values, dtype=np.float64)
        count_valid = bool(
            diagnostics["invalid_value_count"] == 0
            and diagnostics["negative_value_count"] == 0
            and diagnostics["noninteger_value_count"] == 0
        )
        gene_set = set(data.genes)
        resolved = [target for target in sorted_targets if target in gene_set]
        unresolved = [target for target in sorted_targets if target not in gene_set]
        result = {
            "format_version": 1,
            "status": "complete",
            "request": normalized_request,
            "line": line,
            "targets": sorted_targets,
            "target_order": "sorted canonical request targets",
            "target_audit": target_audit,
            "coordinate_coverage": {
                "official_features": N_FEATURES,
                "native_measured_genes": int(len(data.genes)),
                "requested_targets": int(len(sorted_targets)),
                "resolved_target_gene_coordinates": int(len(resolved)),
                "unresolved_target_gene_coordinates": int(len(unresolved)),
                "resolved_targets": resolved,
                "unresolved_targets": unresolved,
            },
            "delta": {
                "path": str(delta_path),
                "shape": list(delta.shape),
                "dtype": str(delta.dtype),
                "file_sha256": delta_file_hash,
                "array_sha256": delta_array_hash,
                "scaled_local_array_sha256": array_sha256(scaled_local),
            },
            "pds": pds,
            "per_target": per_target,
            "identity_parity": {
                "method": parity_method,
                "recomputed_pds": parity_pds,
                "absolute_error": parity_error,
                "atol": 1e-12,
            },
            "random_null": {
                "kind": "uniform random permutation conditioned on no fixed points; reconstructed groups relabeled without reconstruction",
                "repeats": null_repeats,
                "seeds": list(range(null_repeats)),
                "seed_independent_of_reconstruction_and_training": True,
                "orders": orders,
                "pds": null_values,
                "mean": None if not null_values else float(null_array.mean()),
                "std": None if not null_values else float(null_array.std()),
                "minimum": None if not null_values else float(null_array.min()),
                "maximum": None if not null_values else float(null_array.max()),
                "pds_minus_mean": None if not null_values else float(pds - null_array.mean()),
                "per_target": null_per_target,
                "matrix_reuse": bool(null_repeats),
                "matrix_implementation": None if relabel_matrix is None else relabel_matrix["implementation"],
            },
            "count_diagnostics": {
                **diagnostics,
                "finite_nonnegative_integral": count_valid,
                "official_per_cell_limit": OFFICIAL_CELL_LIMIT,
                "official_per_cell_limit_valid": bool(
                    count_valid
                    and diagnostics["maximum_predicted_cell_total"] <= OFFICIAL_CELL_LIMIT
                ),
            },
            "prediction_bulk_sha256": array_sha256(pred_bulk[1]),
            "truth_bulk_sha256": array_sha256(truth[1]),
            "normalization": f"installed pseudobulk log1p CP{int(E.BULK_TARGET_SUM)}",
            "source": {
                "query_controls_used_for_reconstruction": True,
                "reference_controls_used_for_truth_control_origin": True,
                "query_reference_disjoint": True,
                "provenance": data.provenance,
                "input_hashes": data.input_hashes,
            },
            "implementation_hashes": implementations,
            "phase_seconds": {
                "line_load": loaded_at - started,
                "delta_validation_and_hashing": prepared_at - loaded_at,
                "reconstruction": reconstructed_at - prepared_at,
                "installed_identity_pds": identity_scored_at - reconstructed_at,
                "installed_relabel_matrix": matrix_built_at - identity_scored_at,
                "derive_random_nulls": nulls_derived_at - matrix_built_at,
                "total_before_write": nulls_derived_at - started,
            },
            "interpretation": "Installed raw PDS diagnostic only; no full six-member score or official submission.",
        }
        _write_new_json(output_path, result)
        return result


def dispatch(worker: ContextForcingScoreWorker, request: Mapping) -> tuple[dict, bool]:
    try:
        command = str(request.get("command", "score"))
        if command == "quit":
            return {"ok": True, "result": {"status": "bye"}}, True
        if command == "score":
            return {"ok": True, "result": worker.score(request)}, False
        if command == "export_query":
            result = worker.export_query_features(
                str(request["line"]),
                None if request.get("output_path") is None else pathlib.Path(str(request["output_path"])),
            )
            return {"ok": True, "result": result}, False
        raise ValueError(f"unknown command {command!r}")
    except Exception as error:  # Protocol errors are data, not worker termination.
        return {
            "ok": False,
            "error": {"type": type(error).__name__, "message": str(error)},
        }, False


def _requests_from_file(path: pathlib.Path) -> list[dict]:
    text = path.read_text()
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        return [json.loads(line) for line in text.splitlines() if line.strip()]
    if isinstance(value, dict):
        return [value]
    if isinstance(value, list) and all(isinstance(item, dict) for item in value):
        return value
    raise ValueError("request file must contain one JSON object, an object list, or JSONL")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=pathlib.Path)
    parser.add_argument("--report-root", type=pathlib.Path, default=REPORT)
    parser.add_argument("--cache-root", type=pathlib.Path)
    args = parser.parse_args(argv)
    worker = ContextForcingScoreWorker(args.report_root, args.cache_root)
    try:
        if args.request is not None:
            requests = _requests_from_file(args.request)
            for request in requests:
                with contextlib.redirect_stdout(sys.stderr):
                    reply, should_quit = dispatch(worker, request)
                print(json.dumps(reply, sort_keys=True, allow_nan=False), flush=True)
                if should_quit:
                    break
            return 0
        for line in sys.stdin:
            if not line.strip():
                continue
            try:
                request = json.loads(line)
                if not isinstance(request, dict):
                    raise ValueError("each JSONL request must be an object")
                with contextlib.redirect_stdout(sys.stderr):
                    reply, should_quit = dispatch(worker, request)
            except Exception as error:
                reply = {
                    "ok": False,
                    "error": {"type": type(error).__name__, "message": str(error)},
                }
                should_quit = False
            print(json.dumps(reply, sort_keys=True, allow_nan=False), flush=True)
            if should_quit:
                break
        return 0
    finally:
        worker.close()


if __name__ == "__main__":
    raise SystemExit(main())
