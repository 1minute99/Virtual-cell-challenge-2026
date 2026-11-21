#!/usr/bin/env python3
"""Prepare immutable full-local pseudobulk inputs for simple residual models.

This script performs no download and never reads response-derived model outputs.
It exhausts the four local Replogle/Nadig primary screens and the 157 retained
Tahoe raw-count shards.  Earlier pilot artifacts remain read-only.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import pathlib
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

import context_signal_tahoe_effects as TE
import denoising_track_replogle_panel as RN
import denoising_track_models as DM
import h5util as H
import vcclib as V


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "reports/simple_residual_2026-09-16/data"
TAHOE = ROOT / "reports/context_signal_2026-09-15/tahoe"
MIN_CELLS = 25


def write_json(path, value):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_repo_path(path, base=None):
    path = pathlib.Path(path)
    if not path.is_absolute():
        path = (pathlib.Path.cwd() if base is None else pathlib.Path(base)) / path
    return str(path.resolve().relative_to(ROOT.resolve()))


def unicode_array(values):
    return np.asarray([str(value) for value in values], dtype=str)


def logcp_from_library(sum_counts, full_library_sums):
    counts = np.asarray(sum_counts, dtype=np.float64)
    library = np.asarray(full_library_sums, dtype=np.float64)
    if counts.ndim == 1:
        counts = counts[None]
        scalar = True
    else:
        scalar = False
    if library.ndim == 0:
        library = library[None]
    if library.shape != (len(counts),) or np.any(library <= 0):
        raise ValueError("Every pseudobulk requires one positive full-native-library denominator")
    out = np.log1p(counts * 1e4 / library[:, None])
    return out[0] if scalar else out


def batch_matched_control_sums(control_sums, control_cells, target_batch_cells):
    control_sums = np.asarray(control_sums, dtype=np.float64)
    control_cells = np.asarray(control_cells, dtype=np.float64)
    weights = np.asarray(target_batch_cells, dtype=np.float64)
    if control_sums.ndim != 2 or weights.ndim != 2:
        raise ValueError("Control sums and target-by-batch counts must be matrices")
    if len(control_sums) != len(control_cells) or weights.shape[1] != len(control_cells):
        raise ValueError("Batch axes disagree")
    means = np.divide(control_sums, control_cells[:, None], out=np.zeros_like(control_sums),
                      where=control_cells[:, None] > 0)
    return weights @ means


def per_line_target_indices(targets, counts, minimum=MIN_CELLS):
    targets = np.asarray(targets)
    counts = np.asarray(counts)
    if counts.ndim != 2 or counts.shape[1] != len(targets):
        raise ValueError("counts must be line by target")
    return [np.flatnonzero(row >= minimum) for row in counts]


_FAMILIES = {
    "A549": "A549", "CVCL0023": "A549",
    "HT29": "HT29", "CVCL0320": "HT29",
    "HEPG2": "HEPG2_FAMILY", "HEPG2C3A": "HEPG2_FAMILY",
    "C3A": "HEPG2_FAMILY", "CVCL1098": "HEPG2_FAMILY",
    "K562": "K562", "RPE1": "RPE1", "JURKAT": "JURKAT",
    "MCF7": "MCF7", "HAP1": "HAP1", "BXPC3": "BXPC3",
}


def biological_family(value):
    key = "".join(c for c in str(value).upper() if c.isalnum())
    return _FAMILIES.get(key, key)


def assert_disjoint_biologies(training, heldout):
    left = {biological_family(v) for v in training}
    right = {biological_family(v) for v in heldout}
    overlap = sorted(left & right)
    if overlap:
        raise ValueError(f"biological-group overlap: {overlap}")


def smiles_fingerprint(smiles, n_bits=512):
    """Stable response-free hashed SMILES character-ngram fingerprint."""
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")
    text = "" if smiles is None or pd.isna(smiles) else str(smiles).strip()
    out = np.zeros(int(n_bits), dtype=np.uint8)
    if not text:
        return out
    padded = "^" + text + "$"
    tokens = {padded[i:i+n] for n in range(1, 6)
              for i in range(max(0, len(padded) - n + 1))}
    for token in tokens:
        digest = hashlib.sha256(token.encode("utf8")).digest()
        out[int.from_bytes(digest[:8], "little") % n_bits] = 1
        out[int.from_bytes(digest[8:16], "little") % n_bits] = 1
    return out


def select_tahoe_rows(samples, lines, plates, conditions, cell_counts, minimum=MIN_CELLS):
    samples, lines, plates, conditions, cell_counts = map(np.asarray,
        (samples, lines, plates, conditions, cell_counts))
    if not (len(samples) == len(lines) == len(plates) == len(conditions) == len(cell_counts)):
        raise ValueError("Tahoe group columns have unequal lengths")
    control = np.array(["DMSO_TF" in str(v).upper() for v in conditions])
    keep = (~control) & (cell_counts >= minimum)
    indices = np.flatnonzero(keep)
    return {"indices": indices, "is_repeat": plates[indices] == "plate14"}


def canonical_condition(value):
    parsed = ast.literal_eval(str(value)) if isinstance(value, str) else value
    return json.dumps([[str(a), float(b), str(c)] for a, b, c in parsed],
                      ensure_ascii=False, separators=(",", ":"))


def _collapse_genes(genes, *matrices):
    genes = np.asarray(genes).astype(str)
    order = list(dict.fromkeys(genes.tolist()))
    if len(order) == len(genes):
        return np.asarray(order), matrices
    lookup = {g: i for i, g in enumerate(order)}
    codes = np.asarray([lookup[g] for g in genes])
    collapsed = []
    for matrix in matrices:
        matrix = np.asarray(matrix)
        out = np.zeros(matrix.shape[:-1] + (len(order),), dtype=matrix.dtype)
        for old, new in enumerate(codes):
            out[..., new] += matrix[..., old]
        collapsed.append(out)
    return np.asarray(order), tuple(collapsed)


def _rn_aggregate(path, meta, targets, destination):
    start = time.monotonic()
    n, b, g = len(targets), len(meta["batch_names"]), len(meta["genes"])
    target_map = {target: i for i, target in enumerate(targets)}
    codes = np.asarray([target_map.get(target, -1) for target in meta["target"]], dtype=np.int64)
    codes[meta["control"]] = n + meta["batch"][meta["control"]]
    codes[~meta["allowed"]] = -1
    treated = (codes >= 0) & (codes < n)
    target_batch_counts = np.asarray(sparse.coo_matrix((np.ones(treated.sum(), dtype=np.int64),
        (codes[treated], meta["batch"][treated])), shape=(n, b)).toarray())
    sums = np.zeros((n + b, g), dtype=np.float64)
    with h5py.File(path, "r") as handle:
        matrix = handle["X"]
        chunk_rows = matrix.chunks[0]
        block_rows = chunk_rows * max(1, 8192 // chunk_rows)
        for lo in range(0, len(codes), block_rows):
            hi = min(lo + block_rows, len(codes))
            local = codes[lo:hi]
            valid = local >= 0
            if valid.any():
                block = np.asarray(matrix[lo:hi])
                if not np.isfinite(block).all() or block.min() < 0 or not np.allclose(block, np.rint(block)):
                    raise ValueError(f"{path.name} is not finite nonnegative integer raw counts")
                unique, inverse = np.unique(local[valid], return_inverse=True)
                reducer = sparse.csr_matrix((np.ones(len(inverse)),
                    (inverse, np.arange(len(inverse)))), shape=(len(unique), len(inverse)))
                sums[unique] += reducer @ block[valid]
            if lo == 0 or hi == len(codes) or (lo // block_rows) % 25 == 0:
                print(json.dumps({"dataset": "replogle_nadig", "file": path.name,
                    "cells_scanned": hi, "total_cells": len(codes),
                    "elapsed_seconds": time.monotonic() - start}), flush=True)
    control_sums = sums[n:]
    matched = batch_matched_control_sums(control_sums, meta["control_n"], target_batch_counts)
    target_cells = target_batch_counts.sum(1)
    source_library = sums[:n].sum(1)
    matched_library = matched.sum(1)
    genes, collapsed = _collapse_genes(meta["genes"], sums[:n], matched, control_sums)
    source_counts, matched_counts, batch_control_counts = collapsed
    effect = (logcp_from_library(source_counts, source_library) -
              logcp_from_library(matched_counts, matched_library)).astype(np.float32)
    pooled_counts = batch_control_counts.sum(0)
    pooled_library = float(control_sums.sum())
    pooled_control = logcp_from_library(pooled_counts, pooled_library).astype(np.float32)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, targets=np.asarray(targets), genes=genes,
        effect=effect, pooled_control=pooled_control, counts=target_cells,
        target_cell_counts=target_cells, sum_counts=source_counts,
        source_full_library_sums=source_library,
        matched_control_sums=matched_counts,
        matched_control_full_library_sums=matched_library,
        batch_control_sum_counts=batch_control_counts,
        batch_control_cell_counts=meta["control_n"], batches=meta["batch_names"],
        target_batch_cell_counts=target_batch_counts,
        pooled_control_sum_counts=pooled_counts,
        pooled_control_full_library_sum=np.asarray(pooled_library),
        measured_mask=np.ones(len(genes), dtype=bool),
        metadata=np.asarray(json.dumps({"minimum_target_cells": MIN_CELLS,
            "normalization": "log1p(CP10K) using complete native measured library before gene collapse/restriction",
            "controls": "within-screen batch means weighted by each target's cells per batch",
            "biological_replicate": False, "source": str(path.relative_to(ROOT))}, sort_keys=True)))
    return {"targets": len(targets), "genes": len(genes), "cells": int(target_cells.sum()),
            "seconds": time.monotonic() - start, "bytes": destination.stat().st_size,
            "sha256": sha256(destination)}


def prepare_replogle_nadig(out):
    destination = pathlib.Path(out) / "replogle_nadig"
    gm = V.GeneMap()
    rows = []
    for line, filename in RN.FILES[:4]:
        path = ROOT / "raw/scperturb" / filename
        meta = RN.metadata(path, gm)
        targets = np.asarray(sorted(meta["counts"][meta["counts"] >= MIN_CELLS].index))
        cache = destination / line / "pseudobulk.npz"
        if cache.exists():
            raise RuntimeError(f"Refusing to overwrite immutable output {cache}")
        summary = _rn_aggregate(path, meta, targets, cache)
        rows.append({"line": line, "biological_group": biological_family(line),
            "path": relative_repo_path(cache), **summary})
    write_json(destination / "manifest.json", {"status": "complete", "lines": rows,
        "selection": "all canonical targets with >=25 cells in that line and a batch containing controls; ragged, no intersection",
        "normalization": "raw sums -> full-native-library log1p(CP10K) -> target-batch-matched control subtraction",
        "outer_split_rule": "hold out the complete biological family across every source/study/modality",
        "source_files": [str((ROOT / "raw/scperturb" / f).relative_to(ROOT)) for _, f in RN.FILES[:4]]})
    return rows


def _tahoe_groups():
    coverage = pd.read_csv(TAHOE / "candidate_shard_coverage.csv")
    coverage = coverage.groupby(["sample", "cell_line", "plate"], as_index=False).cells.sum()
    samples = pd.read_parquet(TAHOE / "raw/metadata/sample_metadata.parquet")
    groups = coverage.merge(samples, on=["sample", "plate"], how="left", validate="many_to_one")
    if groups.drug.isna().any():
        raise ValueError("Tahoe local groups lack sample metadata")
    groups["condition"] = groups.drugname_drugconc.map(canonical_condition)
    groups = groups.sort_values(["plate", "sample", "cell_line"], kind="mergesort").reset_index(drop=True)
    groups["group_id"] = np.arange(len(groups))
    return groups


def _tahoe_gene_axis():
    spec = json.loads((TAHOE / "gene_spec.json").read_text())
    genes = np.asarray(spec["genes"])
    metadata = pd.read_parquet(TAHOE / "raw/metadata/gene_metadata.parquet")
    index = {g: i for i, g in enumerate(genes)}
    lookup = np.full(int(metadata.token_id.max()) + 1, -1, dtype=np.int64)
    for row in metadata.itertuples():
        if row.gene_symbol in index:
            lookup[int(row.token_id)] = index[row.gene_symbol]
    return genes, lookup


def _tahoe_line_names():
    frame = pd.read_parquet(TAHOE / "raw/metadata/cell_line_metadata.parquet")
    frame = frame[["Cell_ID_Cellosaur", "cell_name"]].drop_duplicates("Cell_ID_Cellosaur")
    return dict(zip(frame.Cell_ID_Cellosaur.astype(str), frame.cell_name.astype(str)))


def _tahoe_drug_features(row_drugs):
    meta = pd.read_parquet(TAHOE / "raw/metadata/drug_metadata.parquet")
    exact = {str(row.drug): row for row in meta.itertuples()}
    stripped = {str(row.drug).strip(): row for row in meta.itertuples()}
    names = sorted(set(map(str, row_drugs)))
    smiles, features = [], []
    for name in names:
        row = exact.get(name, stripped.get(name.strip()))
        value = "" if row is None or pd.isna(row.canonical_smiles) else str(row.canonical_smiles)
        smiles.append(value)
        features.append(smiles_fingerprint(value))
    lookup = {name: i for i, name in enumerate(names)}
    return (np.asarray(names), np.asarray(smiles), np.stack(features),
            np.asarray([lookup[str(v)] for v in row_drugs], dtype=np.int32))


def prepare_tahoe(out):
    start = time.monotonic()
    destination = pathlib.Path(out) / "tahoe" / "table.npz"
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite immutable output {destination}")
    groups = _tahoe_groups()
    genes, lookup = _tahoe_gene_axis()
    key_to_group = {(r.sample, r.cell_line): int(r.group_id) for r in groups.itertuples()}
    acc = TE.CountAccumulator(len(groups), lookup, len(genes), True)
    plan = json.loads((TAHOE / "candidate_shard_plan.json").read_text())
    shards = plan["shards"] if isinstance(plan, dict) else plan
    for number, shard in enumerate(shards, 1):
        relative = shard["path"]
        path = TAHOE / "raw" / relative
        if not path.exists() or path.stat().st_size != int(shard["bytes"]):
            raise FileNotFoundError(f"Missing or size-mismatched frozen Tahoe shard {path}")
        parquet = __import__("pyarrow.parquet", fromlist=["ParquetFile"]).ParquetFile(path)
        seen = 0
        for batch in parquet.iter_batches(batch_size=256,
                columns=["genes", "expressions", "sample", "cell_line_id", "plate"], use_threads=False):
            samples = batch.column("sample").to_pylist()
            lines = batch.column("cell_line_id").to_pylist()
            plates = batch.column("plate").to_pylist()
            if not all(p == shard["expected_plate"] for p in plates):
                raise ValueError("Tahoe shard plate differs from frozen manifest")
            codes = np.asarray([key_to_group.get((s, line), -1) for s, line in zip(samples, lines)])
            ga, xa = batch.column("genes"), batch.column("expressions")
            offsets = ga.offsets.to_numpy()
            if not np.array_equal(offsets, xa.offsets.to_numpy()):
                raise ValueError("Tahoe gene/count offsets disagree")
            offsets, tokens, values = TE.remove_cls(offsets, ga.values.to_numpy(), xa.values.to_numpy())
            acc.add_flat(codes, offsets, tokens, values)
            seen += len(batch)
        if seen != int(shard["expected_rows"]):
            raise ValueError("Tahoe shard row count differs from frozen manifest")
        if number == 1 or number % 10 == 0 or number == len(shards):
            print(json.dumps({"dataset": "tahoe", "shards": number, "total_shards": len(shards),
                "cells_scanned": acc.total_cells_seen, "elapsed_seconds": time.monotonic() - start}), flush=True)
    np.testing.assert_array_equal(acc.cells, groups.cells.to_numpy())
    if acc.total_cells_seen != int(groups.cells.sum()):
        raise ValueError("Tahoe scanned-cell total differs from metadata coverage")
    selected = select_tahoe_rows(groups["sample"], groups["cell_line"], groups["plate"],
        groups["drug"], acc.cells, MIN_CELLS)
    take = selected["indices"]
    control = groups.drug.eq("DMSO_TF")
    control_keys, control_sums, control_libraries, control_cells = [], [], [], []
    for (plate, line), part in groups[control].groupby(["plate", "cell_line"], sort=True):
        ids = part.group_id.to_numpy(dtype=np.int64)
        control_keys.append((plate, line))
        control_sums.append(acc.sums[ids].sum(0))
        control_libraries.append(acc.full_library_sums[ids].sum())
        control_cells.append(acc.cells[ids].sum())
    control_sums = np.asarray(control_sums)
    control_libraries = np.asarray(control_libraries)
    control_cells = np.asarray(control_cells)
    control_index = {key: i for i, key in enumerate(control_keys)}
    row_control = np.asarray([control_index[(groups.plate[i], groups.cell_line[i])] for i in take])
    n, g = len(take), len(genes)
    effects = np.empty((n, g), dtype=np.float32)
    controls = np.empty((n, g), dtype=np.float32)
    treated_sum_counts = np.empty((n, g), dtype=np.uint32)
    for lo in range(0, n, 256):
        hi = min(n, lo + 256)
        source = acc.sums[take[lo:hi]]
        if source.max(initial=0) > np.iinfo(np.uint32).max or not np.array_equal(source, np.rint(source)):
            raise ValueError("Tahoe raw pseudobulk sums cannot be represented exactly as uint32")
        treated_sum_counts[lo:hi] = source.astype(np.uint32)
        cids = row_control[lo:hi]
        c = logcp_from_library(control_sums[cids], control_libraries[cids])
        x = logcp_from_library(source, acc.full_library_sums[take[lo:hi]])
        controls[lo:hi], effects[lo:hi] = c.astype(np.float32), (x - c).astype(np.float32)
    names = _tahoe_line_names()
    line_ids = unicode_array(groups.cell_line.iloc[take])
    line_names = np.asarray([names.get(v, v) for v in line_ids])
    drug_names, drug_smiles, drug_features, drug_feature_index = _tahoe_drug_features(
        groups.drug.iloc[take].astype(str).to_numpy())
    metadata = {"minimum_treated_cells": MIN_CELLS, "local_shards": len(shards),
        "cells_scanned": int(acc.total_cells_seen), "groups_aggregated": len(groups),
        "normalization": "raw count sum -> log1p(CP10K) with complete native RNA library denominator; fixed 18150-gene challenge intersection retained",
        "controls": "both DMSO_TF wells pooled in raw-count space within exact plate and line",
        "primary_plates": ["plate3", "plate6"], "repeat_plate": "plate14",
        "outer_split_rule": "exclude complete biological family across every plate, repeat, study, and modality",
        "fingerprint": "response-free SHA256-hashed canonical-SMILES character ngrams length1-5, 512 bits; RDKit unavailable",
        "source_revision": "2dc57900b7981cfcf5e211527169a0b006546a95"}
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, effects=effects, controls=controls,
        treated_sum_counts=treated_sum_counts,
        treated_full_library_sums=acc.full_library_sums[take],
        counts=acc.cells[take], treatment_cell_counts=acc.cells[take],
        targets=unicode_array(groups.condition.iloc[take]),
        drug_names_by_row=unicode_array(groups.drug.iloc[take]),
        lines=line_names, line_ids=line_ids,
        biological_groups=np.asarray([biological_family(v) for v in line_names]),
        plates=unicode_array(groups.plate.iloc[take]),
        samples=unicode_array(groups["sample"].iloc[take]),
        is_repeat=selected["is_repeat"], genes=genes,
        measured_mask=np.ones(len(genes), dtype=bool),
        control_profile_index=row_control,
        control_profile_keys=np.asarray([f"{p}|{line}" for p, line in control_keys]),
        control_sum_counts=control_sums.astype(np.uint32),
        control_full_library_sums=control_libraries,
        control_cell_counts=control_cells,
        drug_feature_names=drug_names, drug_smiles=drug_smiles,
        drug_features=drug_features, drug_feature_index=drug_feature_index,
        drug_feature_method=np.asarray(metadata["fingerprint"]),
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)))
    summary = {"status": "complete", "path": relative_repo_path(destination),
        "rows": n, "primary_rows": int((~selected["is_repeat"]).sum()),
        "repeat_rows": int(selected["is_repeat"].sum()),
        "lines": int(len(set(line_ids))), "conditions": int(groups.condition.iloc[take].nunique()),
        "genes": len(genes), "cells_scanned": int(acc.total_cells_seen),
        "bytes": destination.stat().st_size, "sha256": sha256(destination),
        "seconds": time.monotonic() - start, **metadata}
    write_json(destination.parent / "manifest.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=("all", "replogle_nadig", "tahoe"), default="all")
    parser.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    result = {"status": "complete", "dataset": args.dataset, "threads": args.threads}
    with threadpool_limits(args.threads):
        if args.dataset in ("all", "replogle_nadig"):
            result["replogle_nadig"] = prepare_replogle_nadig(args.out)
        if args.dataset in ("all", "tahoe"):
            result["tahoe"] = prepare_tahoe(args.out)
    result["seconds"] = time.monotonic() - started
    write_json(args.out / "complete.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
