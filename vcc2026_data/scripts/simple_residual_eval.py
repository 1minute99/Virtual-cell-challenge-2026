#!/usr/bin/env python3
"""Bounded partial-panel count reconstruction and raw-PDS mock evaluation.

This utility deliberately reports the installed scorer's raw PDS only.  A prepared
panel contains capped, metadata-sampled treated cells and controls from one public
Replogle/Nadig screen; it is not a complete challenge reconstruction and carries no
calibrated six-member score.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import pathlib
import subprocess
import sys
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import h5util
import harness_pred as H
import vcclib as V

# Avoid scripts/coverage.py shadowing the third-party coverage package during scorer imports.
_HERE = pathlib.Path(__file__).resolve().parent
sys.path[:] = [p for p in sys.path if pathlib.Path(p or ".").resolve() != _HERE]
from cell_eval2.competition import BULK_TARGET_SUM
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.prep import bulk_lognorm_means, pseudobulk_bulk_lognorm

ROOT = _HERE.parent
CONTROL = H.CONTROL
LINE_FILES = {
    "K562": ROOT / "raw/scperturb/ReplogleWeissman2022_K562_gwps.h5ad",
    "RPE1": ROOT / "raw/scperturb/ReplogleWeissman2022_rpe1.h5ad",
    "HepG2": ROOT / "raw/scperturb/NadigOConner2024_hepg2.h5ad",
    "Jurkat": ROOT / "raw/scperturb/NadigOConner2024_jurkat.h5ad",
}


def _write_json(path: pathlib.Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _source_metadata(path: pathlib.Path, gm: V.GeneMap) -> tuple[np.ndarray, np.ndarray, str]:
    with h5py.File(path, "r") as h:
        raw_targets = h5util.read_col(h, "obs/gene")
        targets = np.asarray(
            [CONTROL if str(x) == CONTROL else gm.canon(str(x)) for x in raw_targets],
            dtype=object,
        )
        gene_key = "var/gene_name" if "gene_name" in h["var"] else (
            "var/" + h5util.index_name(h, "var")
        )
        genes = np.asarray([gm.canon(str(x)) for x in h5util.read_col(h, gene_key)], dtype=object)
    return targets, genes, gene_key


def _unique_gene_intersection(source_genes: np.ndarray, gm: V.GeneMap):
    official_raw = list(V.official_features())
    official_canon = [gm.canon(g) for g in official_raw]
    source_count = Counter(source_genes.tolist())
    official_count = Counter(official_canon)
    source_pos = {g: i for i, g in enumerate(source_genes) if source_count[g] == 1}
    official_pos = {g: i for i, g in enumerate(official_canon) if official_count[g] == 1}
    genes = sorted(set(source_pos) & set(official_pos))
    return (
        genes,
        np.asarray([source_pos[g] for g in genes], dtype=np.int64),
        np.asarray([official_pos[g] for g in genes], dtype=np.int64),
        {
            "source_canonical_duplicates_excluded": int(sum(n for n in source_count.values() if n > 1)),
            "official_canonical_duplicates_excluded": int(sum(n for n in official_count.values() if n > 1)),
        },
    )


def _read_dense_selected(
    path: pathlib.Path, rows: np.ndarray, columns: np.ndarray, block: int = 256
) -> sparse.csr_matrix:
    """Read only requested dense rows, in bounded blocks, and then retain selected columns."""
    rows = np.asarray(rows, dtype=np.int64)
    if len(np.unique(rows)) != len(rows):
        raise ValueError("sampled source rows are not unique")
    order = np.argsort(rows)
    sorted_rows = rows[order]
    chunks = []
    with h5py.File(path, "r") as h:
        x = h["X"]
        if not isinstance(x, h5py.Dataset):
            raise ValueError(f"{path}: X must be dense for bounded selected-row streaming")
        for lo in range(0, len(sorted_rows), block):
            selected = sorted_rows[lo : lo + block]
            part = np.asarray(x[selected])[:, columns]
            if not np.all(np.isfinite(part)):
                raise ValueError("selected raw counts contain non-finite values")
            if np.any(part < 0):
                raise ValueError("selected raw counts contain negative values")
            if np.any(np.abs(part - np.rint(part)) > 1e-6):
                raise ValueError("selected raw matrix is not integer counts")
            if part.size and part.max() > 2**24:
                raise ValueError("selected counts exceed exact float32 integer range")
            chunks.append(sparse.csr_matrix(np.rint(part).astype(np.float32)))
    sorted_matrix = sparse.vstack(chunks, format="csr")
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    return sorted_matrix[inverse]


def prepare_real(
    line: str,
    targets: list[str],
    out: pathlib.Path,
    real_cells: int = 200,
    control_cells: int = 800,
    seed: int = 20260916,
) -> dict:
    """Prepare a capped mock reference without selecting cells from expression values."""
    if line not in LINE_FILES:
        raise ValueError(f"unknown line {line!r}; choose one of {sorted(LINE_FILES)}")
    if real_cells <= 0 or real_cells > 200:
        raise ValueError("real_cells must be in 1..200")
    if control_cells <= 0 or control_cells > 800:
        raise ValueError("control_cells must be in 1..800")
    out = pathlib.Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    source = pathlib.Path(LINE_FILES[line])
    if not source.is_file():
        raise FileNotFoundError(source)

    gm = V.GeneMap()
    canonical = [gm.canon(str(t)) for t in targets]
    if len(canonical) != len(set(canonical)):
        raise ValueError("requested targets collide after canonical gene mapping")
    canonical = sorted(canonical)
    if len(canonical) < 2:
        raise ValueError("at least two unique targets are required for a deranged-shuffle PDS")
    source_targets, source_genes, gene_key = _source_metadata(source, gm)
    counts = Counter(source_targets.tolist())
    missing = {t: counts.get(t, 0) for t in canonical if counts.get(t, 0) < 25}
    if missing:
        raise ValueError(f"every requested target requires >=25 cells; insufficient: {missing}")
    control_rows = np.flatnonzero(source_targets == CONTROL)
    if len(control_rows) < 400:
        raise ValueError(f"control pool requires >=400 cells; found {len(control_rows)}")

    genes, source_columns, official_indices, collision_audit = _unique_gene_intersection(
        source_genes, gm
    )
    missing_coordinates = sorted(set(canonical) - set(genes))

    rng = np.random.default_rng(seed)
    sampled_target_rows = []
    sampled_counts = {}
    for target in canonical:
        available = np.flatnonzero(source_targets == target)
        take = min(real_cells, len(available))
        chosen = np.sort(rng.choice(available, size=take, replace=False))
        sampled_target_rows.append(chosen)
        sampled_counts[target] = int(take)
    n_controls = min(control_cells, len(control_rows))
    if n_controls < 400:
        raise ValueError(
            f"requested control cap yields {n_controls}; control_cells must permit at least 400"
        )
    sampled_controls = np.sort(rng.choice(control_rows, size=n_controls, replace=False))
    selected_rows = np.concatenate([*sampled_target_rows, sampled_controls])
    selected = _read_dense_selected(source, selected_rows, source_columns)
    n_treated = sum(sampled_counts.values())
    controls = selected[n_treated:].copy().tocsr()
    sparse.save_npz(out / "controls.npz", controls)

    labels = np.concatenate(
        [
            *[np.full(sampled_counts[t], t, dtype=object) for t in canonical],
            np.full(n_controls, CONTROL, dtype=object),
        ]
    )
    obs = pd.DataFrame(
        {
            "target": pd.Categorical(labels, categories=[*canonical, CONTROL]),
            "source_row": selected_rows,
        },
        index=[f"real_{i}" for i in range(len(labels))],
    )
    real = ad.AnnData(X=selected, obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    real.write_h5ad(out / "real.h5ad", compression="gzip")

    metadata = {
        "format_version": 1,
        "line": line,
        "source": str(source.relative_to(ROOT) if source.is_relative_to(ROOT) else source),
        "source_shape": [int(len(source_targets)), int(len(source_genes))],
        "source_gene_column": gene_key,
        "targets": canonical,
        "genes": genes,
        "official_feature_indices": official_indices.tolist(),
        "target_source_cell_counts": {t: int(counts[t]) for t in canonical},
        "target_sample_cell_counts": sampled_counts,
        "target_gene_coordinates_resolved": sorted(set(canonical) & set(genes)),
        "target_gene_coordinates_unresolved": missing_coordinates,
        "source_control_cells": int(len(control_rows)),
        "sampled_control_cells": int(n_controls),
        "seed": int(seed),
        "sample_caps": {"real_cells": int(real_cells), "control_cells": int(control_cells)},
        "mock_partial_coverage": True,
        "limitations": [
            "Only genes in the unique canonical intersection of this screen and the 18,533 official features are retained.",
            "Targets listed in target_gene_coordinates_unresolved have no unique retained transcript coordinate; installed panel exclusion leaves those coordinates unresolved and warns when scoring.",
            "Treated cells are capped at 200 per target and controls at 800; this is a mock partial-coverage panel.",
            "Cells are sampled from metadata indices only; treated expression is never used to choose cells or controls.",
            "The resulting PDS is an offline raw diagnostic, not an official or calibrated six-member score.",
        ],
        **collision_audit,
    }
    _write_json(out / "manifest.json", metadata)
    _write_json(out / "genes.json", genes)
    _write_json(out / "targets.json", canonical)
    return metadata


def _reconstruction_diagnostics(
    controls: sparse.csr_matrix, delta: np.ndarray, cells: int, seed: int
) -> dict:
    totals = np.asarray(controls.sum(1)).ravel().astype(np.float64)
    if np.any(totals <= 0):
        raise ValueError("a control cell has zero counts")
    changed = np.nonzero(np.abs(delta).max(axis=0) > 0)[0]
    rng = np.random.default_rng(seed)
    pick = rng.choice(controls.shape[0], size=cells, replace=False)
    dense = np.asarray(controls[pick][:, changed].todense(), dtype=np.float64)
    lib = totals[pick][:, None]
    invalid = 0
    clipped = 0
    clipped_amount = 0.0
    with np.errstate(over="ignore", invalid="ignore"):
        expd = np.exp(delta[:, changed].astype(np.float64))
        for row in expd:
            raw = ((dense / lib * 1e4 + 1.0) * row - 1.0) * lib / 1e4
            invalid += int(np.count_nonzero(~np.isfinite(raw)))
            negative = raw < 0
            clipped += int(np.count_nonzero(negative))
            clipped_amount += float((-raw[negative]).sum())
    if invalid:
        raise ValueError(f"reconstruction produced {invalid} non-finite count values")
    return {
        "invalid_value_count": 0,
        "clipped_value_count": clipped,
        "clipped_count_amount": clipped_amount,
        "changed_gene_count": int(len(changed)),
    }


def _reconstructed_bulk(
    controls: sparse.csr_matrix,
    delta: np.ndarray,
    targets: Iterable[str],
    cells: int,
    seed: int,
):
    targets = list(targets)
    diagnostics = _reconstruction_diagnostics(controls, delta, cells, seed)
    sums = np.zeros((len(targets) + 1, controls.shape[1]), dtype=np.float64)
    target_totals = []
    min_cell_total = np.inf
    max_cell_total = 0.0
    invalid = noninteger = negative = 0
    for ti, new, unchanged_coo, changed, unchanged in H.reconstructed_blocks(
        controls, delta, cells, seed, "stochastic"
    ):
        values = np.concatenate([new.ravel(), unchanged_coo.data.astype(np.float64)])
        invalid += int(np.count_nonzero(~np.isfinite(values)))
        negative += int(np.count_nonzero(values < 0))
        noninteger += int(np.count_nonzero(np.abs(values - np.rint(values)) > 0))
        if values.size and np.nanmax(values) > 2**24:
            raise ValueError("reconstruction has integer counts outside exact float32 range")
        sums[ti, changed] = new.sum(axis=0, dtype=np.float64)
        sums[ti, unchanged] = np.bincount(
            unchanged_coo.col,
            weights=unchanged_coo.data.astype(np.float64),
            minlength=len(unchanged),
        )
        cell_totals = new.sum(axis=1, dtype=np.float64) + np.bincount(
            unchanged_coo.row,
            weights=unchanged_coo.data.astype(np.float64),
            minlength=cells,
        )
        min_cell_total = min(min_cell_total, float(cell_totals.min(initial=np.inf)))
        max_cell_total = max(max_cell_total, float(cell_totals.max(initial=0)))
        target_totals.append(float(cell_totals.sum()))
    if invalid or negative or noninteger:
        raise ValueError(
            "invalid reconstructed counts: "
            f"nonfinite={invalid}, negative={negative}, noninteger={noninteger}"
        )
    sums[-1] = np.asarray(controls.astype(np.float64).sum(axis=0)).ravel()
    labels = np.asarray([*targets, CONTROL])
    order = np.argsort(labels)
    diagnostics.update(
        {
            "noninteger_value_count": 0,
            "negative_value_count": 0,
            "predicted_count_total": float(sum(target_totals)),
            "per_target_count_total": {t: v for t, v in zip(targets, target_totals)},
            "control_count_total": float(sums[-1].sum()),
            "minimum_predicted_cell_total": float(min_cell_total),
            "maximum_predicted_cell_total": float(max_cell_total),
        }
    )
    return labels[order], bulk_lognorm_means(sums[order], BULK_TARGET_SUM), diagnostics


def _score_bulk(pred_bulk, real_bulk, genes: list[str]) -> dict[str, float]:
    return discrimination_score(
        pred_bulk=pred_bulk,
        real_bulk=real_bulk,
        distance="cosine",
        rank_denominator="n-1",
        tie_policy="midrank",
        exclude_target_gene=True,
        exclusion_scope="panel",
        control_source="real",
        genes=np.asarray(genes),
    )


def _derangement(n: int) -> np.ndarray:
    if n < 2:
        raise ValueError("at least two targets are required for a derangement")
    # A fixed cyclic shift is a deterministic derangement for every n >= 2.
    return np.roll(np.arange(n), 1)


class MockEvaluator:
    """Evaluate target-by-gene deltas on one prepared partial panel."""

    def __init__(self, real_path, cells: int = 400, scale: float = 3.0, seed: int = 0):
        path = pathlib.Path(real_path)
        self.folder = path if path.is_dir() else path.parent
        self.real_path = self.folder / "real.h5ad" if path.is_dir() else path
        manifest_path = self.folder / "manifest.json"
        controls_path = self.folder / "controls.npz"
        if not self.real_path.is_file() or not manifest_path.is_file() or not controls_path.is_file():
            raise FileNotFoundError(
                f"prepared panel requires real.h5ad, manifest.json, controls.npz in {self.folder}"
            )
        self.manifest = json.loads(manifest_path.read_text())
        self.targets = list(self.manifest["targets"])
        self.genes = list(self.manifest["genes"])
        if self.targets != sorted(self.targets) or len(self.targets) != len(set(self.targets)):
            raise ValueError("manifest targets must be sorted and unique")
        if len(self.targets) < 2:
            raise ValueError("at least two targets are required")
        if len(self.genes) != len(set(self.genes)):
            raise ValueError("manifest genes must be unique")
        if not np.isfinite(scale):
            raise ValueError("scale must be finite")
        self.cells = int(cells)
        self.scale = float(scale)
        self.seed = int(seed)
        self.controls = sparse.load_npz(controls_path).tocsr()
        if self.cells <= 0 or self.cells > self.controls.shape[0]:
            raise ValueError(
                f"cells must be in 1..{self.controls.shape[0]}; got {self.cells}"
            )
        if self.controls.shape[1] != len(self.genes):
            raise ValueError("controls gene dimension does not match manifest")
        real = ad.read_h5ad(self.real_path)
        if list(real.var_names.astype(str)) != self.genes:
            raise ValueError("real var order does not match manifest genes")
        observed = sorted(set(real.obs["target"].astype(str)) - {CONTROL})
        if observed != self.targets:
            raise ValueError("real target labels do not match manifest target order/set")
        self._real_bulk = pseudobulk_bulk_lognorm(
            real, "target", bulk_target_sum=BULK_TARGET_SUM
        )
        indices = np.asarray(self.manifest["official_feature_indices"], dtype=np.int64)
        if indices.shape != (len(self.genes),) or np.any(indices < 0) or np.any(indices >= V.N_FEATURES):
            raise ValueError("manifest official_feature_indices are invalid")
        if len(np.unique(indices)) != len(indices):
            raise ValueError("manifest official_feature_indices must be unique")
        self._official_indices = indices

    def _delta(self, delta) -> np.ndarray:
        value = np.asarray(delta)
        local_shape = (len(self.targets), len(self.genes))
        official_shape = (len(self.targets), V.N_FEATURES)
        if not np.all(np.isfinite(value)):
            raise ValueError("delta contains non-finite values")
        if value.shape == local_shape:
            local = value
        elif value.shape == official_shape:
            local = value[:, self._official_indices]
        else:
            raise ValueError(
                f"delta shape {value.shape} is invalid; expected {local_shape} or {official_shape}"
            )
        return np.asarray(local, dtype=np.float32) * self.scale

    def score(self, delta: np.ndarray, seed: int | None = None) -> dict:
        use_seed = self.seed if seed is None else int(seed)
        scaled = self._delta(delta)
        pred_labels, pred_values, diagnostics = _reconstructed_bulk(
            self.controls, scaled, self.targets, self.cells, use_seed
        )
        per_target = _score_bulk((pred_labels, pred_values), self._real_bulk, self.genes)
        order = _derangement(len(self.targets))
        sh_labels, sh_values, sh_diagnostics = _reconstructed_bulk(
            self.controls, scaled[order], self.targets, self.cells, use_seed
        )
        shuffle_per_target = _score_bulk((sh_labels, sh_values), self._real_bulk, self.genes)
        return {
            "pds_raw": float(np.mean(list(per_target.values()))),
            "per_target": per_target,
            "shuffle_pds_raw": float(np.mean(list(shuffle_per_target.values()))),
            "shuffle_per_target": shuffle_per_target,
            "shuffle_source_target": {
                target: self.targets[int(source)] for target, source in zip(self.targets, order)
            },
            "seed": use_seed,
            "cells": self.cells,
            "scale": self.scale,
            "bulk_target_sum": BULK_TARGET_SUM,
            "diagnostics": diagnostics,
            "shuffle_diagnostics": sh_diagnostics,
            "interpretation": "Raw installed-scorer PDS on a capped partial-coverage mock panel; not a calibrated composite or official score.",
        }

    def _prediction_matrix(self, scaled: np.ndarray, seed: int):
        rows, cols, values = [], [], []
        for ti, new, unchanged_coo, changed, unchanged in H.reconstructed_blocks(
            self.controls, scaled, self.cells, seed, "stochastic"
        ):
            rr, cc = np.nonzero(new)
            offset = ti * self.cells
            rows.append(np.concatenate([unchanged_coo.row, rr]) + offset)
            cols.append(np.concatenate([unchanged[unchanged_coo.col], changed[cc]]))
            values.append(
                np.concatenate(
                    [unchanged_coo.data.astype(np.float32), new[rr, cc].astype(np.float32)]
                )
            )
        matrix = sparse.csr_matrix(
            (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
            shape=(len(self.targets) * self.cells, len(self.genes)),
        )
        matrix.data = np.rint(matrix.data).astype(np.float32)
        return matrix

    def write_prediction(self, delta: np.ndarray, out: pathlib.Path, seed: int | None = None):
        use_seed = self.seed if seed is None else int(seed)
        out = pathlib.Path(out)
        if out.exists():
            raise FileExistsError(f"refusing to overwrite {out}")
        out.parent.mkdir(parents=True, exist_ok=True)
        scaled = self._delta(delta)
        # Validate and record the exact same reconstruction before serializing it.
        _, _, diagnostics = _reconstructed_bulk(
            self.controls, scaled, self.targets, self.cells, use_seed
        )
        predicted = self._prediction_matrix(scaled, use_seed)
        matrix = sparse.vstack([predicted, self.controls], format="csr")
        labels = np.concatenate(
            [
                np.repeat(self.targets, self.cells),
                np.full(self.controls.shape[0], CONTROL, dtype=object),
            ]
        )
        result = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(
                {"target": pd.Categorical(labels, categories=[*self.targets, CONTROL])},
                index=[f"pred_{i}" for i in range(matrix.shape[0])],
            ),
            var=pd.DataFrame(index=pd.Index(self.genes)),
        )
        result.write_h5ad(out, compression="gzip")
        metadata = {
            "seed": use_seed,
            "cells": self.cells,
            "scale": self.scale,
            "shape": list(result.shape),
            "diagnostics": diagnostics,
            "mock_partial_coverage": True,
        }
        _write_json(out.with_suffix(".meta.json"), metadata)
        return metadata


def _read_targets(path: pathlib.Path) -> list[str]:
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text())
        if not isinstance(value, list):
            raise ValueError("JSON targets file must contain a list")
        return [str(x) for x in value]
    frame = pd.read_csv(path)
    if "target" not in frame:
        raise ValueError(f"CSV targets file must have a 'target' column; got {list(frame)}")
    return [str(x) for x in frame["target"]]


def _add_eval_arguments(parser):
    parser.add_argument("--real", type=pathlib.Path, required=True)
    parser.add_argument("--delta", type=pathlib.Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--scale", type=float, default=3.0)
    parser.add_argument("--cells", type=int, default=400)


def _full_score_command(executable, prediction, real, outdir) -> list[str]:
    """Build the explicit installed-scorer boundary for generated count h5ads."""
    return [
        str(executable),
        "run",
        "-ap",
        str(prediction),
        "-ar",
        str(real),
        "--preset",
        "vcc2026",
        "--pert-col",
        "target",
        "--control",
        CONTROL,
        "--input-type",
        "counts",
        "-o",
        str(outdir),
        "--set",
        "de.backend=pdex",
    ]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="prepare one capped real-data mock panel")
    prepare.add_argument("--line", required=True, choices=sorted(LINE_FILES))
    prepare.add_argument("--targets-file", type=pathlib.Path, required=True)
    prepare.add_argument("--out", type=pathlib.Path, required=True)
    prepare.add_argument("--real-cells", type=int, default=200)
    prepare.add_argument("--control-cells", type=int, default=800)
    prepare.add_argument("--seed", type=int, default=20260916)
    score = sub.add_parser("score", help="compute raw installed-scorer PDS and shuffle PDS")
    _add_eval_arguments(score)
    score.add_argument("--out", type=pathlib.Path, required=True)
    full = sub.add_parser("full-score", help="write prediction and invoke cell-eval2 separately")
    _add_eval_arguments(full)
    full.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args(argv)

    if args.command == "prepare":
        metadata = prepare_real(
            args.line,
            _read_targets(args.targets_file),
            args.out,
            real_cells=args.real_cells,
            control_cells=args.control_cells,
            seed=args.seed,
        )
        print(json.dumps(metadata, indent=2))
        return

    evaluator = MockEvaluator(args.real, cells=args.cells, scale=args.scale, seed=args.seed)
    delta = np.load(args.delta, allow_pickle=False)
    if args.command == "score":
        result = evaluator.score(delta)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        if args.out.exists():
            raise FileExistsError(f"refusing to overwrite {args.out}")
        _write_json(args.out, result)
        print(json.dumps(result, indent=2))
        return

    if args.out.exists():
        raise FileExistsError(f"full-score output directory must be new: {args.out}")
    args.out.mkdir(parents=True)
    prediction = args.out / "prediction.h5ad"
    evaluator.write_prediction(delta, prediction)
    score_dir = args.out / "cell_eval2"
    executable = pathlib.Path(sys.executable).with_name("cell-eval2")
    command = _full_score_command(executable, prediction, evaluator.real_path, score_dir)
    completed = subprocess.run(command, check=False)
    _write_json(
        args.out / "full_score_invocation.json",
        {"command": command, "returncode": completed.returncode, "raw_pds_only_companion": evaluator.score(delta)},
    )
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
