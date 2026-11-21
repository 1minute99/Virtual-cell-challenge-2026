#!/usr/bin/env python3
"""Apply the frozen response operator to a target-aligned baseline delta."""
import argparse
import json
import pathlib
import h5py
import numpy as np
from threadpoolctl import threadpool_limits
import vcclib as V


def transform(delta, model):
    original = np.asarray(delta, dtype=np.float64)
    pred = original - original.mean(0, keepdims=True)
    cols = model["columns"]
    pred[:, cols] += (pred[:, cols] @ model["basis"].astype(np.float64)) @ model["coefficient"].astype(np.float64)
    factor = np.abs(original).sum() / np.abs(pred).sum()
    return (pred * factor).astype(np.float32), float(factor)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default="processed/cross_source_residual/crsres128/model.npz")
    ap.add_argument("--baseline", default="processed/baseline_deltas/b4_sevencontext_consensus_delta.npy")
    ap.add_argument("--out", default="processed/baseline_deltas/b16_learned_residual_delta.npy")
    a = ap.parse_args()
    model = np.load(a.model)
    features, targets = V.official_features(), V.official_targets()
    np.testing.assert_array_equal(model["features"], features)
    assert not (set(model["training_targets"]) & set(targets))
    baseline = np.load(a.baseline)
    assert baseline.shape == (len(targets), len(features))
    assert np.isfinite(baseline).all()
    # Serialization/inference parity on the already evaluated KOLF input.
    kolf = np.load("harness/xl_delta_mt_kolf.npy")
    expected = np.load("harness/xl_delta_crsres128_kolf.npy")
    with threadpool_limits(4):
        check, _ = transform(kolf, model)
        error = float(np.max(np.abs(check-expected)))
        np.testing.assert_allclose(check, expected, atol=5e-7, rtol=5e-5)
        result, factor = transform(baseline, model)
    assert np.isfinite(result).all()
    assert (np.abs(result).max(1) > 0).all()
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out, result)
    manifest = dict(model=a.model, model_sha256=V.sha256(a.model),
        baseline=a.baseline, baseline_sha256=V.sha256(a.baseline),
        output=str(out), output_sha256=V.sha256(out),
        targets=targets, features=features, shape=list(result.shape),
        inference_parity_max_abs_error=error, magnitude_factor=factor,
        l1_ratio=float(np.abs(result).sum()/np.abs(baseline).sum()),
        algorithm="source-only target centering plus frozen signed residual operator; prediction-only global L1 match",
        context_specific=False, recommended_reconstruction_scale=3.0,
        held_out_challenge_response_used=False)
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:manifest[k] for k in ["output","shape","inference_parity_max_abs_error","magnitude_factor","l1_ratio"]}))


if __name__=="__main__":
    main()
