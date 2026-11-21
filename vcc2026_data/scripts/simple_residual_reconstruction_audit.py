#!/usr/bin/env python3
"""Measure mock-score distortion introduced by count reconstruction and normalization."""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib

import anndata as ad
import numpy as np

import simple_residual_eval as E
from cell_eval2.competition import BULK_TARGET_SUM
from cell_eval2.prep import pseudobulk_bulk_lognorm


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _index(labels) -> dict[str, int]:
    values = [str(x) for x in labels]
    if len(values) != len(set(values)):
        raise ValueError("pseudobulk labels are not unique")
    return {label: i for i, label in enumerate(values)}


def _ordered_bulk(bulk, labels: list[str]) -> np.ndarray:
    names, values = bulk
    positions = _index(names)
    missing = sorted(set(labels) - set(positions))
    if missing:
        raise ValueError(f"pseudobulk is missing labels: {missing[:5]}")
    return np.asarray(values, dtype=np.float64)[[positions[x] for x in labels]]


def _pds(predicted: np.ndarray, real_bulk, targets: list[str], genes: list[str]) -> float:
    labels = np.asarray([*targets, E.CONTROL])
    order = np.argsort(labels)
    per_target = E._score_bulk((labels[order], predicted[order]), real_bulk, genes)
    return float(np.mean(list(per_target.values())))


def _center_rows(value: np.ndarray) -> np.ndarray:
    return value - value.mean(axis=1, keepdims=True)


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    x = np.asarray(left, dtype=np.float64).ravel()
    y = np.asarray(right, dtype=np.float64).ravel()
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    if nx == 0 or ny == 0:
        return float("nan")
    return float(np.dot(x, y) / (nx * ny))


def audit_reconstruction(real, prediction, delta, scale: float = 3.0) -> dict:
    """Return direct, reconstructed-CP10K, and reconstructed-CP50K diagnostics."""
    real = pathlib.Path(real)
    folder = real if real.is_dir() else real.parent
    real_path = folder / "real.h5ad" if real.is_dir() else real
    manifest_path = folder / "manifest.json"
    prediction = pathlib.Path(prediction)
    delta = pathlib.Path(delta)
    manifest = json.loads(manifest_path.read_text())
    targets = [str(x) for x in manifest["targets"]]
    genes = [str(x) for x in manifest["genes"]]
    indices = np.asarray(manifest["official_feature_indices"], dtype=np.int64)
    source_delta = np.load(delta, allow_pickle=False)
    if source_delta.shape == (len(targets), len(genes)):
        local_delta = np.asarray(source_delta, dtype=np.float64)
    elif source_delta.shape == (len(targets), E.V.N_FEATURES):
        local_delta = np.asarray(source_delta[:, indices], dtype=np.float64)
    else:
        raise ValueError(
            f"delta shape {source_delta.shape} does not match local or official mock axes"
        )
    if not np.all(np.isfinite(local_delta)) or not np.isfinite(scale):
        raise ValueError("delta and scale must be finite")
    intended = local_delta * float(scale)

    real_data = ad.read_h5ad(real_path)
    predicted_data = ad.read_h5ad(prediction)
    if list(real_data.var_names.astype(str)) != genes:
        raise ValueError("real gene order differs from the prepared manifest")
    if list(predicted_data.var_names.astype(str)) != genes:
        raise ValueError("prediction gene order differs from the prepared manifest")
    observed = sorted(set(predicted_data.obs["target"].astype(str)) - {E.CONTROL})
    if observed != targets:
        raise ValueError("prediction target labels differ from the prepared manifest")

    real10 = pseudobulk_bulk_lognorm(real_data, "target", bulk_target_sum=10_000.0)
    real50 = pseudobulk_bulk_lognorm(real_data, "target", bulk_target_sum=BULK_TARGET_SUM)
    pred10 = pseudobulk_bulk_lognorm(predicted_data, "target", bulk_target_sum=10_000.0)
    pred50 = pseudobulk_bulk_lognorm(predicted_data, "target", bulk_target_sum=BULK_TARGET_SUM)
    labels = [*targets, E.CONTROL]
    real10_values = _ordered_bulk(real10, labels)
    pred10_values = _ordered_bulk(pred10, labels)
    pred50_values = _ordered_bulk(pred50, labels)
    control10 = real10_values[-1]
    direct10_values = np.vstack([control10[None, :] + intended, control10])

    direct_pds = _pds(direct10_values, real10, targets, genes)
    reconstructed10_pds = _pds(pred10_values, real10, targets, genes)
    reconstructed50_pds = _pds(pred50_values, real50, targets, genes)
    realized = pred10_values[:-1] - pred10_values[-1]
    centered_intended = _center_rows(intended)
    centered_realized = _center_rows(realized)
    intended_norm = float(np.linalg.norm(centered_intended))
    realized_norm = float(np.linalg.norm(centered_realized))
    ratios = np.divide(
        np.linalg.norm(centered_realized, axis=1),
        np.linalg.norm(centered_intended, axis=1),
        out=np.full(len(targets), np.nan, dtype=np.float64),
        where=np.linalg.norm(centered_intended, axis=1) > 0,
    )
    correlations = np.asarray(
        [_correlation(x, y) for x, y in zip(centered_intended, centered_realized)],
        dtype=np.float64,
    )

    return {
        "format_version": 1,
        "targets": len(targets),
        "genes": len(genes),
        "scale": float(scale),
        "pds": {
            "intended_direct_cp10k": direct_pds,
            "reconstructed_cp10k": reconstructed10_pds,
            "reconstructed_cp50k": reconstructed50_pds,
            "reconstruction_minus_direct_cp10k": reconstructed10_pds - direct_pds,
            "cp50k_minus_cp10k_reconstructed": reconstructed50_pds - reconstructed10_pds,
        },
        "centered_delta": {
            "intended_frobenius_norm": intended_norm,
            "realized_cp10k_frobenius_norm": realized_norm,
            "realized_to_intended_norm_ratio": (
                realized_norm / intended_norm if intended_norm > 0 else float("nan")
            ),
            "pearson": _correlation(centered_intended, centered_realized),
            "per_target_norm_ratio_mean": float(np.nanmean(ratios)),
            "per_target_norm_ratio_median": float(np.nanmedian(ratios)),
            "per_target_pearson_mean": float(np.nanmean(correlations)),
            "per_target_pearson_median": float(np.nanmedian(correlations)),
        },
        "direct_log_expression": {
            "negative_coordinate_count": int(np.count_nonzero(direct10_values[:-1] < 0)),
            "minimum": float(direct10_values[:-1].min()),
            "note": "The intended diagnostic is real-control log1p(CP10K) plus scale*delta without feasibility clipping.",
        },
        "axes": {
            "target_order": targets,
            "gene_order": genes,
            "target_gene_exclusion": "installed panel exclusion on this same measured gene axis",
        },
        "inputs": {
            "real": str(real_path),
            "real_sha256": _sha256(real_path),
            "prediction": str(prediction),
            "prediction_sha256": _sha256(prediction),
            "delta": str(delta),
            "delta_sha256": _sha256(delta),
        },
        "limitations": [
            "The mock retains only the unique canonical gene intersection of the source screen and official feature axis.",
            "Real treated cells are capped at 200 per target and controls at 800; predictions use 400 sampled controls per target.",
            "Direct CP10K is a log-space pseudobulk diagnostic and is not an official score or a feasible count reconstruction.",
            "Reconstructed CP50K uses the installed scorer normalization but remains a partial-panel mock diagnostic.",
        ],
    }


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real", type=pathlib.Path, required=True)
    parser.add_argument("--prediction", type=pathlib.Path, required=True)
    parser.add_argument("--delta", type=pathlib.Path, required=True)
    parser.add_argument("--scale", type=float, default=3.0)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args(argv)
    if args.out.exists():
        raise FileExistsError(f"refusing to overwrite {args.out}")
    result = audit_reconstruction(args.real, args.prediction, args.delta, args.scale)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["pds"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
