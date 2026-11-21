#!/usr/bin/env python3
"""End-to-end checks on what is actually on disk. Prints PASS/FAIL per check.

Sample checks and exhaustive checks are labelled separately, because "20 groups
recomputed by hand agree" and "every group agrees" are different claims. The
independent recomputation reads the source file directly rather than trusting the
reduction that produced the store.
"""
import json
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import vccdata as VD
import vcclib as V

RESULTS = []


def check(name, scope, fn):
    try:
        ok, detail = fn()
    except Exception as e:                                     # a failure is a result
        ok, detail = False, f"{type(e).__name__}: {e}"
    RESULTS.append(dict(check=name, scope=scope, passed=bool(ok), detail=str(detail)[:400]))
    print(f"[{'PASS' if ok else 'FAIL'}] ({scope}) {name}: {str(detail)[:200]}", flush=True)
    return ok


def main():
    sp = VD.OfficialSpace()

    check("official target and feature files load with the right shape and header",
          "exhaustive",
          lambda: (len(sp.targets) == 300 and len(sp.features) == 18533,
                   f"{len(sp.targets)} targets, {len(sp.features)} features, "
                   f"pert_col={sp.manifest['pert_col']}, "
                   f"control_label={sp.manifest['control_label']}"))

    check("target list and feature list are not the same file",
          "exhaustive",
          lambda: (set(sp.targets) != set(sp.features) and len(set(sp.targets) &
                                                              set(sp.features)) == 300,
                   f"all 300 targets are also readout features, but the lists differ "
                   f"in length ({len(sp.targets)} vs {len(sp.features)})"))

    def ctrl_counts():
        bad = []
        for c in ["A", "B", "C"]:
            d = sp.control_cells(c)
            X = d["X"]
            v = X.data
            if X.shape != (18400, 18533):
                bad.append(f"{c} shape {X.shape}")
            if not np.all(np.isfinite(v)) or np.any(v < 0):
                bad.append(f"{c} has negative or non-finite values")
            if not np.all(np.abs(v - np.rint(v)) < 1e-6):
                bad.append(f"{c} has non-integral values")
            L = np.asarray(X.sum(1)).ravel()
            if L.min() <= 0:
                bad.append(f"{c} has an empty cell")
            if len(set(d["ntc_id"])) != 46:
                bad.append(f"{c} has {len(set(d['ntc_id']))} ntc ids, expected 46")
        return (not bad, "3 contexts x 18,400 cells: counts non-negative, integral, "
                         "finite, no empty cell, 46 NTC ids each" if not bad else "; ".join(bad))
    check("official control cells are usable raw counts", "exhaustive", ctrl_counts)

    keys = VD.list_pb()
    check("pseudobulk stores present", "exhaustive",
          lambda: (len(keys) > 0, f"{len(keys)}: {keys}"))

    for k in keys:
        src = VD.load_pb(k)

        def shapes(src=src, k=k):
            ok = (src.sums.shape[0] == len(src.obs)
                  and src.sums.shape[1] == len(src.var)
                  and len(src.official_to_source) == 18533)
            return ok, (f"{src.sums.shape} vs obs {len(src.obs)} / var {len(src.var)}; "
                        f"{int(src.measured.sum())}/18533 official features measured")
        check(f"{k}: obs/var/matrix alignment", "exhaustive", shapes)

        def sem(src=src, k=k):
            s = src.semantics
            return ("sum" in s.lower() and "umi" in s.lower(),
                    f"count_semantics = {s[:150]}")
        check(f"{k}: count semantics recorded as a sum of UMIs", "exhaustive", sem)

        def values(src=src):
            v = src.sums.data
            return (bool(np.all(v >= 0) and np.all(np.isfinite(v))
                         and np.all(np.abs(v - np.rint(v)) < 1e-6)),
                    f"{src.sums.nnz:,} stored values, min {v.min()}, max {v.max()}, "
                    f"all integral")
        check(f"{k}: stored sums are non-negative finite integers", "exhaustive", values)

        def ncells(src=src):
            n = src.obs.n_cells.to_numpy()
            tot = np.asarray(src.sums.sum(1)).ravel()
            bad = int(((n <= 0) & (tot > 0)).sum())
            return (bad == 0, f"{len(n)} groups, n_cells median {np.median(n):.0f}, "
                              f"{int((n == 0).sum())} empty groups, "
                              f"{bad} groups with counts but no cells")
        check(f"{k}: n_cells agrees with the presence of counts", "exhaustive", ncells)

        def controls(src=src):
            isc = src.obs.is_control.astype(str).str.lower().isin(["true", "1"])
            per = src.obs[isc].groupby("stimulation").size().to_dict() if isc.any() else {}
            return (isc.sum() > 0, f"{int(isc.sum())} control groups; per stimulation {per}")
        check(f"{k}: control groups exist per context", "exhaustive", controls)

        def matched(src=src):
            isc = src.obs.is_control.astype(str).str.lower().isin(["true", "1"]).to_numpy()
            rows = np.nonzero(~isc)[0][:200]
            miss = [int(r) for r in rows if len(src.matched_control_rows(r)) == 0]
            return (len(miss) == 0, f"200 sampled perturbation groups all have a "
                                    f"matched control on the same stimulation"
                                    if not miss else f"{len(miss)} without a control")
        check(f"{k}: matched controls resolve", "sample (200 groups)", matched)

        def dupes(src=src):
            # Use the store's OWN declared group keys, not a fixed list. VCC 2025 groups
            # on `split` as well, and its 48 control groups legitimately recur in each
            # split because the release ships the same control cells with all three -- a
            # fixed key would flag those 96 rows as duplicates when the real duplication
            # is of cells across files, which is recorded in the provenance instead.
            kk = [c for c in src.provenance.get("group_keys", []) if c in src.obs.columns]
            if not kk:
                kk = [c for c in ["target_symbol", "stimulation", "donor",
                                  "technical_batch"] if c in src.obs.columns]
            d = int(src.obs.duplicated(kk).sum())
            return (d == 0, f"no duplicate identity on the declared keys {kk}" if d == 0
                    else f"{d} duplicate rows on {kk}")
        check(f"{k}: group identity is unique", "exhaustive", dupes)

        def official(src=src):
            A, m, s = src.official_view(rows=np.arange(min(8, src.sums.shape[0])))
            unmeasured_all_zero = bool(np.all(A[:, ~m] == 0))
            return (A.shape[1] == 18533 and unmeasured_all_zero,
                    f"official view {A.shape}, unmeasured columns are 0 and masked "
                    f"({int((~m).sum())} of 18,533); semantics: {s[:80]}")
        check(f"{k}: official 18,533-feature view + measured mask", "sample (8 groups)",
              official)

    # independent recomputation straight from the CD4 source file
    def recompute():
        key = "cd4_target_condition_donor"
        d = V.ROOT / "pseudobulk" / key
        if not (d / "sum_counts.npz").exists():
            return False, "store missing"
        src = VD.load_pb(key)
        prov = src.provenance
        srcfile = pathlib.Path(prov["source_file"])
        with h5py.File(srcfile, "r") as h:
            o = pd.DataFrame({c: VD.h5_col(h, f"obs/{c}") for c in
                              ["perturbed_gene_name", "culture_condition", "donor_id",
                               "10xrun_id", "guide_type", "n_cells"]})
            o["n_cells"] = pd.to_numeric(o.n_cells, errors="coerce").fillna(0)
            tgt = np.where(o.guide_type.astype(str).str.lower()
                           .str.contains("non-targeting"), "non-targeting",
                           o.perturbed_gene_name.astype(str))
            rng = np.random.default_rng(0)
            picks = rng.choice(len(src.obs), size=6, replace=False)
            ok = []
            X = h["X"]
            for pi in picks:
                r = src.obs.iloc[pi]
                sel = np.nonzero((tgt == str(r.target_symbol))
                                 & (o.culture_condition.astype(str) == str(r.stimulation))
                                 & (o.donor_id.astype(str) == str(r.donor))
                                 & (o["10xrun_id"].astype(str) == str(r.technical_batch)))[0]
                tot = 0.0
                cells = 0.0
                for rr in sel:
                    lo, hi = int(X["indptr"][rr]), int(X["indptr"][rr + 1])
                    tot += float(X["data"][lo:hi].sum())
                    cells += float(o.n_cells.iloc[rr])
                stored = float(src.sums[pi].sum())
                ok.append((abs(tot - stored) < 1.0, int(r.n_cells) == int(cells),
                           r.target_symbol, tot, stored, int(r.n_cells), int(cells)))
            bad = [x for x in ok if not (x[0] and x[1])]
            return (not bad, f"6 random groups recomputed from {srcfile.name}: totals and "
                             f"n_cells match exactly" if not bad else f"mismatches: {bad}")
    check("independent recomputation of pseudobulk sums from the source file",
          "sample (6 groups)", recompute)

    # effect stores
    ekeys = VD.list_effects()
    check("effect stores present", "exhaustive",
          lambda: (len(ekeys) > 0, f"{len(ekeys)}: {ekeys}"))
    for k in ekeys:
        e = VD.load_effect(k)

        def eshape(e=e):
            return (e.delta.shape[1] == 18533 and e.delta.shape[0] == len(e.pairs),
                    f"delta {e.delta.shape}, {len(e.pairs)} pairs, "
                    f"{int(e.measured.sum())}/18533 measured")
        check(f"{k}: delta shape and pair alignment", "exhaustive", eshape)

        def enan(e=e):
            un = ~e.measured
            all_nan = bool(np.all(np.isnan(e.delta[:, un]))) if un.any() else True
            fin = bool(np.all(np.isfinite(e.delta[:, e.measured])))
            return (all_nan and fin,
                    "unmeasured features are NaN (not 0) and measured features are finite")
        check(f"{k}: unmeasured encoded as NaN, measured finite", "exhaustive", enan)

        def esem(e=e):
            s = e.semantics.lower()
            ref = ("control" in s) or ("non-targeting" in s)
            return ("log1p" in s and "cp10k" in s and ref,
                    f"quantity = {e.semantics[:120]}")
        check(f"{k}: effect quantity is named precisely", "exhaustive", esem)

    # baseline deltas and built submissions
    bd = V.ROOT / "processed/baseline_deltas"
    for f in sorted(bd.glob("*_delta.npy")) if bd.exists() else []:
        def dcheck(f=f):
            D = np.load(f)
            return (D.shape == (300, 18533) and np.all(np.isfinite(D)),
                    f"{D.shape}, {(np.abs(D).max(1) > 0).sum()}/300 targets carry a "
                    f"signal, {(np.abs(D).max(0) > 0).sum()}/18533 genes can move")
        check(f"baseline delta {f.name}", "exhaustive", dcheck)

    # loader batch
    def batch():
        if not keys:
            return False, "no stores"
        b = next(VD.iter_pb_batches([keys[0]], batch=16, targets_only=sp.targets))
        return (b["X"].shape == (16, 18533) and b["measured"].dtype == bool,
                f"batch X {b['X'].shape}, mask {b['measured'].sum()}/18533 measured, "
                f"targets {list(b['obs'].target_symbol[:3])}, "
                f"group_key={b['group_key']}")
    check("training loader returns a batch in official gene order with a mask",
          "sample (1 batch)", batch)

    D = pd.DataFrame(RESULTS)
    out = V.ROOT / "reports/smoke_test.csv"
    D.to_csv(out, index=False)
    n_ok = int(D.passed.sum())
    print(f"\n{n_ok}/{len(D)} checks passed -> {out}")
    (V.ROOT / "reports/smoke_test.json").write_text(json.dumps(dict(
        passed=n_ok, total=len(D),
        failed=[r for r in RESULTS if not r["passed"]]), indent=2))
    return 0 if n_ok == len(D) else 1


if __name__ == "__main__":
    sys.exit(main())
