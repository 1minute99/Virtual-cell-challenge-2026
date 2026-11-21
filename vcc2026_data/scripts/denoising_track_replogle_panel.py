#!/usr/bin/env python3
"""Consistent count-sum, batch-matched Replogle/Nadig genetic pilot from local files."""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

import denoising_track_models as D
import h5util as H
import vcclib as V

FILES = [
    ("K562", "ReplogleWeissman2022_K562_gwps.h5ad"),
    ("RPE1", "ReplogleWeissman2022_rpe1.h5ad"),
    ("HepG2", "NadigOConner2024_hepg2.h5ad"),
    ("Jurkat", "NadigOConner2024_jurkat.h5ad"),
    ("K562_essential_screen", "ReplogleWeissman2022_K562_essential.h5ad"),
]


def metadata(path, gm):
    with h5py.File(path, "r") as h:
        target = H.read_col(h, "obs/gene")
        mapping = {v: gm.canon(str(v)) for v in set(target) if v != "non-targeting"}
        mapping["non-targeting"] = "non-targeting"
        target = np.asarray([mapping[v] for v in target])
        batch_names, batch = np.unique(H.read_col(h, "obs/batch").astype(str), return_inverse=True)
        gene_key = "var/gene_name" if "gene_name" in h["var"] else "var/" + H.index_name(h, "var")
        genes = np.asarray([gm.canon(str(g)) for g in H.read_col(h, gene_key)])
        control = target == "non-targeting"
        nc = np.bincount(batch[control], minlength=len(batch_names))
        allowed = nc[batch] > 0
        counts = pd.Series(target[allowed & ~control]).value_counts()
        info = dict(path=str(path.relative_to(D.ROOT)), bytes=path.stat().st_size,
                    shape=list(h["X"].shape), chunks=list(h["X"].chunks),
                    n_batches=len(batch_names), n_batches_with_controls=int((nc > 0).sum()),
                    n_control_cells=int(control.sum()), no_control_batch_cells=int((~allowed).sum()),
                    gene_column=gene_key, n_canonical_duplicate_genes=len(genes) - len(set(genes)))
    return dict(target=target, batch=batch, batch_names=batch_names, genes=genes,
                control=control, control_n=nc, allowed=allowed, counts=counts, audit=info)


def logcp(x):
    total = x.sum(-1, keepdims=True)
    return np.log1p(np.divide(x * 1e4, total, out=np.zeros_like(x, dtype=np.float64), where=total > 0))


def aggregate(path, meta, targets, folder):
    start = time.monotonic()
    n, b, g = len(targets), len(meta["batch_names"]), len(meta["genes"])
    target_map = {t: i for i, t in enumerate(targets)}
    codes = np.asarray([target_map.get(t, -1) for t in meta["target"]], dtype=np.int64)
    codes[meta["control"]] = n + meta["batch"][meta["control"]]
    codes[~meta["allowed"]] = -1
    keep = (codes >= 0) & (codes < n)
    target_batch_counts = np.asarray(sparse.coo_matrix((np.ones(keep.sum(), dtype=np.int64),
        (codes[keep], meta["batch"][keep])), shape=(n, b)).toarray())
    sums = np.zeros((n + b, g), dtype=np.float64)
    sample_integrality_error, minimum = 0.0, np.inf
    with h5py.File(path, "r") as h:
        x = h["X"]
        chunk_rows = x.chunks[0]
        block_rows = chunk_rows * max(1, 8192 // chunk_rows)
        for lo in range(0, len(codes), block_rows):
            hi = min(lo + block_rows, len(codes))
            local = codes[lo:hi]
            valid = local >= 0
            if not valid.any():
                continue
            block = np.asarray(x[lo:hi])
            assert np.isfinite(block).all()
            minimum = min(minimum, float(block.min()))
            sample = block[::max(1, len(block) // 64), ::max(1, g // 128)]
            sample_integrality_error = max(sample_integrality_error, float(np.abs(sample - np.rint(sample)).max()))
            assert minimum >= 0 and sample_integrality_error < 1e-6
            unique, inverse = np.unique(local[valid], return_inverse=True)
            reducer = sparse.csr_matrix((np.ones(len(inverse)), (inverse, np.arange(len(inverse)))),
                                        shape=(len(unique), len(inverse)))
            sums[unique] += reducer @ block[valid]
            if lo == 0 or hi == len(codes) or lo // block_rows % 25 == 0:
                print(json.dumps(dict(file=path.name, cells_scanned=hi, total_cells=len(codes),
                                      elapsed_seconds=time.monotonic() - start)), flush=True)
    control_sums = sums[n:]
    control_means = np.divide(control_sums, meta["control_n"][:, None],
                              out=np.zeros_like(control_sums), where=meta["control_n"][:, None] > 0)
    matched_control_sums = target_batch_counts @ control_means
    target_counts = target_batch_counts.sum(1)
    effect = logcp(sums[:n]) - logcp(matched_control_sums)
    effect[target_counts == 0] = np.nan
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / "pseudobulk.npz", targets=targets, genes=meta["genes"],
        sum_counts=sums[:n], matched_control_sums=matched_control_sums,
        control_sum_counts=control_sums, control_cell_counts=meta["control_n"],
        batches=meta["batch_names"], target_batch_cell_counts=target_batch_counts,
        target_cell_counts=target_counts, effect=effect.astype(np.float32),
        pooled_control=logcp(control_sums.sum(0)[None])[0].astype(np.float32),
        batch_control=logcp(control_sums).astype(np.float32))
    audit = dict(**meta["audit"], target_cells_kept=int(target_counts.sum()), n_targets_with_cells=int((target_counts > 0).sum()),
                 n_targets_min25=int((target_counts >= 25).sum()), min_count=float(minimum),
                 sampled_integrality_max_error=sample_integrality_error,
                 elapsed_seconds=time.monotonic() - start, preprocessing="Full measured gene library log1pCP10K sums; batch-matched control mean weighted by target cells per batch",
                 target_counts_and_control_counts_saved=True, sha256=D.sha(folder / "pseudobulk.npz"))
    D.write_json(folder / "audit.json", audit)
    return effect, logcp(control_sums.sum(0)[None])[0], target_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=D.OUT / "replogle_nadig_panel")
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if (args.out / "panel.npz").exists():
        raise RuntimeError("Existing completed panel is immutable")
    gm = V.GeneMap()
    metas = [metadata(D.ROOT / "raw/scperturb" / file, gm) for _, file in FILES]
    excluded = set(json.loads((D.ROOT / "processed/cross_source_residual/crsres128/manifest.json").read_text())["excluded_evaluation_targets"])
    target_sets = [set(m["counts"][m["counts"] >= 25].index) for m in metas[:4]]
    targets = np.asarray(sorted(set.intersection(*target_sets) - excluded))
    common_genes = set.intersection(*[set(m["genes"]) for m in metas])
    common_genes &= {gm.canon(g) for g in V.official_features()}
    genes = np.asarray(sorted(common_genes - set(targets) - excluded))
    assert len(targets) > 500 and len(genes) > 2000
    fingerprints = {str((D.ROOT / "raw/scperturb" / file).relative_to(D.ROOT)):
                    D.sha(D.ROOT / "raw/scperturb" / file) for _, file in FILES}
    spec = D.freeze_spec(args.out / "spec.json", dict(
        primary_lines=[name for name, _ in FILES[:4]], separate_screen=FILES[4][0],
        targets=targets.tolist(), genes=genes.tolist(), min_target_cells=25,
        target_selection="Metadata cell counts and coverage only; all four primary lines >=25 cells in batches with controls; exclude official and KOLF evaluation targets",
        gene_selection="Measured canonical gene intersection across five screens and official features; exclude all benchmark targets plus official/KOLF panel coordinates",
        normalization="Count SUM across cells then log1pCP10K over each full measured feature library before canonical merging and intersection",
        controls="Within-screen batch control mean raw counts mixed using target cell numbers in each batch; batches without controls excluded",
        counts_semantics="scPerturb publication explicitly supplies scRNAseq count matrices; counts checked for finite/nonnegative/integral sampled entries",
        biological_replicates="Batch labels not treated as biological replicates. K562 essential is a separate screen, never an independent context",
        no_context_ids_as_model_inputs=True, source_sha256=fingerprints,
        metadata=[m["audit"] for m in metas],
        provenance=["https://www.nature.com/articles/s41592-023-02144-y",
                    "https://zenodo.org/records/13350497"],
        license="CC BY 4.0 per local registry and public Zenodo release",
        threads=args.threads), ["scripts/denoising_track_replogle_panel.py", "scripts/h5util.py", "scripts/vcclib.py"])
    start = time.monotonic()
    effects, controls, cells = [], [], []
    with threadpool_limits(args.threads):
        for (name, file), meta in zip(FILES, metas):
            folder = args.out / name
            if (folder / "pseudobulk.npz").exists():
                saved = np.load(folder / "pseudobulk.npz")
                np.testing.assert_array_equal(saved["targets"], targets)
                effect, control, count = saved["effect"], saved["pooled_control"], saved["target_cell_counts"]
            else:
                effect, control, count = aggregate(D.ROOT / "raw/scperturb" / file, meta, targets, folder)
            # If canonical aliases collide, average log effects; normalization was already full-library.
            mapping = {g: np.flatnonzero(meta["genes"] == g) for g in genes}
            effects.append(np.stack([effect[:, mapping[g]].mean(1) for g in genes], axis=1))
            controls.append(np.array([control[mapping[g]].mean() for g in genes]))
            cells.append(count)
    np.savez_compressed(args.out / "panel.npz", targets=targets, genes=genes,
                        lines=np.asarray([name for name, _ in FILES[:4]]),
                        effects=np.stack(effects[:4]).astype(np.float32),
                        controls=np.stack(controls[:4]).astype(np.float32),
                        target_cell_counts=np.stack(cells[:4]),
                        k562_essential_screen=effects[4].astype(np.float32),
                        k562_essential_screen_cells=cells[4])
    pd.DataFrame(np.stack(cells, axis=1), index=targets, columns=[name for name, _ in FILES]).to_csv(args.out / "coverage.csv", index_label="target")
    D.write_json(args.out / "summary.json", dict(stage="complete", n_targets=len(targets), n_genes=len(genes),
                 primary_lines=4, separate_k562_screen=True, elapsed_seconds=time.monotonic() - start,
                 panel_sha256=D.sha(args.out / "panel.npz"), no_new_downloads=True,
                 biological_replicates_not_claimed=True))
    print((args.out / "summary.json").read_text(), flush=True)


if __name__ == "__main__":
    main()
