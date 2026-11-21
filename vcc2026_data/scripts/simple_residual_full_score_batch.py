#!/usr/bin/env python3
"""Run immutable simple-residual mock candidates through installed cell-eval2 serially."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import simple_residual_reconstruction_audit as R
E = R.E


RAW_SIX = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: pathlib.Path, value) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def thread_environment(base: dict[str, str] | None = None) -> dict[str, str]:
    result = dict(os.environ if base is None else base)
    for key in (
        "NUMBA_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        result[key] = "2"
    return result


def score_command(executable, prediction, real, outdir) -> list[str]:
    return [
        str(executable), "run", "-ap", str(prediction), "-ar", str(real),
        "--preset", "vcc2026", "--pert-col", "target", "--control", E.CONTROL,
        "--input-type", "counts", "-o", str(outdir),
        "--set", "de.backend=pdex", "--set", "num_threads=2",
        "--set", "gather_threads=2",
    ]


def require_pds_parity(full_score: float, companion: float, atol: float = 1e-12) -> float:
    error = abs(float(full_score) - float(companion))
    if not np.isfinite(error) or error > atol:
        raise RuntimeError(
            f"full-score PDS parity failed: {full_score} != {companion}; "
            f"absolute error {error} exceeds {atol}"
        )
    return error


def normalize_runs(runs, mock_root: pathlib.Path, out_root: pathlib.Path) -> list[dict]:
    if not isinstance(runs, list) or not runs:
        raise ValueError("plan must contain a non-empty run list")
    normalized = []
    destinations = set()
    for item in runs:
        if not isinstance(item, dict):
            raise ValueError("every run must be an object")
        line = str(item["line"])
        name = str(item["name"])
        if line not in E.LINE_FILES:
            raise ValueError(f"unknown mock line: {line}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise ValueError(f"unsafe run name: {name!r}")
        delta = pathlib.Path(item["delta"]).resolve()
        if not delta.is_file():
            raise FileNotFoundError(delta)
        digest = sha256(delta)
        declared = item.get("delta_file_sha256")
        if declared is not None and str(declared) != digest:
            raise ValueError(
                f"declared delta SHA differs for {line}/{name}: {declared} != {digest}"
            )
        destination = pathlib.Path(out_root).resolve() / line / name
        if destination in destinations:
            raise ValueError(f"duplicate output destination: {destination}")
        destinations.add(destination)
        normalized.append(
            {
                "line": line,
                "name": name,
                "delta": delta,
                "delta_sha256": digest,
                "real": pathlib.Path(mock_root).resolve() / line / "panel300",
                "out": destination,
            }
        )
    return normalized


def _completed(out: pathlib.Path) -> bool:
    marker = out / "COMPLETE"
    receipt = out / "receipt.json"
    if not marker.is_file() or not receipt.is_file():
        return False
    state = json.loads(marker.read_text())
    return state.get("receipt_sha256") == sha256(receipt)


def _freeze(run: dict, companion: dict, command: list[str], resource_path: pathlib.Path) -> None:
    out = run["out"]
    score = out / "cell_eval2"
    agg = pd.read_csv(score / "agg_results.csv", index_col=0)
    metric_agg = pd.read_csv(score / "metric_aggregation.csv").set_index("metric")
    raw = {
        metric: {
            "value": float(agg.loc["mean", metric]),
            "aggregation": str(metric_agg.loc[metric, "agg"]),
            "n_used": int(metric_agg.loc[metric, "n_used"]),
            "n_rows": int(metric_agg.loc[metric, "n_rows"]),
            "n_nan": int(metric_agg.loc[metric, "n_nan"]),
        }
        for metric in RAW_SIX
    }
    pds_parity_error = require_pds_parity(
        raw["pds_cosine"]["value"], float(companion["pds_raw"])
    )
    if sha256(run["delta"]) != run["delta_sha256"]:
        raise RuntimeError("immutable delta changed during scoring")
    audit_path = out / "reconstruction_audit.json"
    audit = R.audit_reconstruction(run["real"], out / "prediction.h5ad", run["delta"], 3.0)
    _write_json(audit_path, audit)
    manifest = json.loads((run["real"] / "manifest.json").read_text())
    artifact_names = [
        "prediction.h5ad", "prediction.meta.json", "reconstruction_audit.json",
        "cell_eval2/agg_results.csv", "cell_eval2/metric_aggregation.csv",
        "cell_eval2/results.csv", "cell_eval2/run_meta.json", "cell_eval2/run_params.yaml",
    ]
    artifacts = {
        rel: {"bytes": (out / rel).stat().st_size, "sha256": sha256(out / rel)}
        for rel in artifact_names
    }
    diagnostics = companion["diagnostics"]
    receipt = {
        "format_version": 1,
        "created_utc": utcnow(),
        "status": "complete",
        "purpose": "Installed cell-eval2 raw six-member partial-panel mock score; no calibrated overall and no official submission.",
        "line": run["line"],
        "name": run["name"],
        "seed": 0,
        "cells_per_target": 400,
        "scale": 3.0,
        "delta": {"path": str(run["delta"]), "sha256": run["delta_sha256"]},
        "raw_six_order": RAW_SIX,
        "raw_six_metrics": raw,
        "raw_pds_companion": {
            "candidate_pds": float(companion["pds_raw"]),
            "same_candidate_deranged_pds": float(companion["shuffle_pds_raw"]),
            "pds_minus_deranged": float(companion["pds_raw"] - companion["shuffle_pds_raw"]),
            "installed_full_score_pds_matches_cached": True,
            "absolute_parity_error": pds_parity_error,
            "parity_atol": 1e-12,
        },
        "reconstruction": {
            "finite_nonnegative_integral_counts": True,
            "maximum_predicted_cell_total": diagnostics["maximum_predicted_cell_total"],
            "minimum_predicted_cell_total": diagnostics["minimum_predicted_cell_total"],
            "official_per_cell_count_limit": 1_000_000,
            "official_per_cell_count_limit_valid": diagnostics["maximum_predicted_cell_total"] < 1_000_000,
            "aggregate_predicted_count_total": diagnostics["predicted_count_total"],
            "aggregate_total_exceeds_1e6": diagnostics["predicted_count_total"] > 1_000_000,
            "aggregate_total_note": "The 1,000,000 official limit is per cell, not the aggregate across predicted cells.",
            "clipped_value_count": diagnostics["clipped_value_count"],
            "clipped_count_amount": diagnostics["clipped_count_amount"],
        },
        "coverage_and_metric_limits": {
            "measured_genes": len(manifest["genes"]),
            "resolved_target_gene_coordinates": len(manifest["target_gene_coordinates_resolved"]),
            "unresolved_target_gene_coordinates": len(manifest["target_gene_coordinates_unresolved"]),
            "unresolved_targets": manifest["target_gene_coordinates_unresolved"],
            "lfc_nmae_n_used": raw["de_wilcoxon_lfc_nmae"]["n_used"],
            "direction_reach_n_used": raw["de_wilcoxon_direction_reach_raw"]["n_used"],
        },
        "reconstruction_distortion": audit,
        "scorer": {
            "command": command,
            "thread_environment": {
                key: "2" for key in (
                    "NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                )
            },
            "resource_log": resource_path.name,
            "exit_status": 0,
        },
        "artifact_hashes": artifacts,
        "interpretation": "All metrics and distortion values are offline partial-panel mock diagnostics.",
    }
    _write_json(out / "receipt.json", receipt)
    _write_json(
        out / "COMPLETE",
        {
            "created_utc": utcnow(),
            "receipt": "receipt.json",
            "receipt_bytes": (out / "receipt.json").stat().st_size,
            "receipt_sha256": sha256(out / "receipt.json"),
            "status": "complete",
        },
    )


def execute_run(run: dict, python: pathlib.Path) -> None:
    out = run["out"]
    if _completed(out):
        print(f"skip complete {run['line']}/{run['name']}", flush=True)
        return
    if out.exists():
        raise FileExistsError(f"refusing incomplete/non-empty destination: {out}")
    if not (run["real"] / "real.h5ad").is_file():
        raise FileNotFoundError(run["real"] / "real.h5ad")
    out.mkdir(parents=True)
    evaluator = E.MockEvaluator(run["real"], cells=400, scale=3.0, seed=0)
    delta = np.load(run["delta"], allow_pickle=False)
    companion = evaluator.score(delta)
    evaluator.write_prediction(delta, out / "prediction.h5ad")
    del evaluator, delta
    gc.collect()
    scorer = pathlib.Path(python).with_name("cell-eval2")
    command = score_command(
        scorer, out / "prediction.h5ad", run["real"] / "real.h5ad", out / "cell_eval2"
    )
    resource_path = out / "resource.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(resource_path), *command]
    with (out / "score.stdout.log").open("w") as stdout, (out / "score.stderr.log").open("w") as stderr:
        completed = subprocess.run(
            timed, stdout=stdout, stderr=stderr, env=thread_environment(), check=False
        )
    _write_json(
        out / "full_score_invocation.json",
        {"command": command, "returncode": completed.returncode, "raw_pds_only_companion": companion},
    )
    if completed.returncode:
        raise RuntimeError(f"cell-eval2 failed for {run['line']}/{run['name']}")
    _freeze(run, companion, command, resource_path)
    print(f"complete {run['line']}/{run['name']}", flush=True)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=pathlib.Path, required=True)
    parser.add_argument("--mock-root", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    args = parser.parse_args(argv)
    value = json.loads(args.plan.read_text())
    runs = normalize_runs(value.get("runs", value) if isinstance(value, dict) else value,
                          args.mock_root, args.out_root)
    for run in runs:
        execute_run(run, pathlib.Path(sys.executable))


if __name__ == "__main__":
    main()
