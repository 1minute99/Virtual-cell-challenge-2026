"""Tests for output/data-analysis/rq3_energy.py (Task 15; plan input/ideas/zero-shot-limits-context-complexity-plan.md,
operational definition 7, H3a). Synthetic cell-level Poisson data with a known low-rank control structure; no real
truth-vs-prediction computation anywhere in this file.

 (i)   zero residual + noise                 -> E[F numerator] ~ 0 (within MC error), also E[denominator] ~ 0
 (ii)  residual planted on the dest control span -> F(U_dest) >> F(U_perm), D_other > 0 with CI excluding 0
 (iii) isotropic residual (true D = 0): shared-baseline variant is biased upward on control-PC subspaces,
       the registered three-way split is not
 (iv)  noise-only negative control (pseudo-targets from controls, dpred = 0) -> numerator-unit D ~ 0
 (v)   P_U is an orthogonal projector; F invariant to rotation of the basis within U; randomized PCs = exact PCs;
       control split deterministic/disjoint/reusable; prediction loader aligns to G-cap
 smoke real H1 controls via common.py -> U_dest, U_perm (k = 8) on G-cap
"""
from pathlib import Path
import sys, time

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import rq3_energy as E  # noqa: E402
import common  # noqa: E402

G, Q, S_LAT = 300, 4, 3.0          # genes, latent control factors, column norm of the loadings (log scale)
KS = (4, 8)


def _loadings(rng):
    W = rng.standard_normal((G, Q))
    return W / np.linalg.norm(W, axis=0) * S_LAT


def _cells(rng, m, W, n, e=None):
    eta = m + rng.standard_normal((n, Q)) @ W.T + (0 if e is None else e)
    X = rng.poisson(np.exp(eta)).astype(float)
    return X, X.sum(1)


def _pb_exp(m, W, e=None):
    lam = np.exp(m + (0 if e is None else e) + 0.5 * (W ** 2).sum(1))
    return np.log1p(1e4 * lam / lam.sum())


class World:
    """Destination panel (controls + T targets) and n_other other panels' controls on the same genes."""

    def __init__(self, seed=0, n_ctrl=3000, T=150, nt=(60, 160), n_other=3):
        rng = np.random.default_rng(seed)
        self.m = rng.normal(1.0, 0.7, G)
        self.W = _loadings(rng)
        self.Qd = np.linalg.qr(self.W)[0]                       # true dest control span (G x Q)
        self.cX, self.cL = _cells(rng, self.m, self.W, n_ctrl)
        base = _pb_exp(self.m, self.W)
        self.tX, self.tL, self.dexp, self.names = [], [], [], []
        for t in range(T):
            e = np.zeros(G); g = rng.choice(G, 10, replace=False); e[g] = rng.choice([-0.7, 0.7], 10)
            X, L = _cells(rng, self.m, self.W, int(rng.integers(*nt)), e)
            self.tX.append(X); self.tL.append(L); self.dexp.append(_pb_exp(self.m, self.W, e) - base); self.names.append(f"t{t}")
        self.dexp = np.stack(self.dexp)
        self.others = [f"o{j}" for j in range(n_other)]
        self.U_other = {}
        for j in self.others:
            X, L = _cells(rng, self.m, _loadings(rng), 2000)
            self.U_other[j] = E.top_pcs(E.cell_log(X, L), max(KS), seed=0)[0]

    def parts(self, rho, seeds=(0, 1, 2), n_perm=5, scheme="three_way", negctrl=False):
        out = []
        for s in seeds:
            split = E.common.split_halves(len(self.cL), s, k=3)
            subs = E.subspaces_for_seed(E.cell_log(self.cX[split[0]], self.cL[split[0]]), self.U_other, s, max(KS), n_perm)
            if negctrl:
                rA, rB, _ = E.negctrl_residuals(self.cX, self.cL, split, [len(x) for x in self.tL], s, n_pseudo=300)
            else:
                hs = [E.target_half_sums(X, L, s, n) for X, L, n in zip(self.tX, self.tL, self.names)]
                rA, rB = E.residuals(hs, self.cX, self.cL, split, self.dexp - rho, scheme=scheme)
            out.append(E.seed_parts_from_residuals(rA, rB, subs, self.others, n_perm, KS))
        return out


def _row(rows, q, k=8):
    return next(r for r in rows if r["quantity"] == q and r["k"] == k)


@pytest.fixture(scope="module")
def world():
    return World(seed=0)


def test_i_zero_residual_numerator_zero(world):
    parts = world.parts(np.zeros_like(world.dexp))
    for key in ("dest", "haar"):
        v = np.concatenate([p[key][:, 1] for p in parts])        # k = 8
        z = v.mean() / (v.std(ddof=1) / np.sqrt(len(v)))
        assert abs(z) < 4, (key, v.mean(), z)
    v = np.concatenate([p["other"][:, :, 1].mean(1) for p in parts])
    assert abs(v.mean() / (v.std(ddof=1) / np.sqrt(len(v)))) < 4
    d = np.concatenate([p["den"] for p in parts])
    assert abs(d.mean() / (d.std(ddof=1) / np.sqrt(len(d)))) < 4
    # and small against the signal scale of test (ii): per-target ||rho||^2 = 4 * 0.3^2 = 0.36
    assert abs(np.concatenate([p["dest"][:, 1] for p in parts]).mean()) < 0.036


def test_ii_planted_on_dest(world):
    rng = np.random.default_rng(1)
    rho = rng.normal(0, 0.3, (len(world.names), Q)) @ world.Qd.T
    rows = E.summarize(world.parts(rho), KS, B=500, others=world.others)
    Fd, Fp, Fo = _row(rows, "F_dest")["point"], _row(rows, "F_perm_mean")["point"], _row(rows, "F_other_mean")["point"]
    assert Fd > 0.6 and Fd > 5 * max(Fp, 1e-3), (Fd, Fp)
    assert Fd > 5 * max(Fo, 1e-3)
    assert _row(rows, "D_other")["lo"] > 0 and _row(rows, "D_perm")["lo"] > 0
    assert abs(_row(rows, "F_haar")["point"] - 8 / G) < 0.05


def test_iii_shared_baseline_bias_removed_by_three_way():
    """Isotropic residual -> true D_other = D_perm = 0. The bias of the shared-baseline variant is +||P_U eps_C||^2 per
    target with ONE baseline-noise draw eps_C, so a single realisation is chi-square-noisy; the claim is about the
    expectation, tested over R independent synthetic replicates (new world, residual and splits each time)."""
    R = 16
    D = {sc: {q: [] for q in ("D_perm", "D_other")} for sc in ("shared", "three_way")}
    for rep in range(R):
        w = World(seed=100 + rep, n_ctrl=300, T=100)   # thirds of 100 cells (RN panels: 800 controls -> 266)
        rho = np.random.default_rng(200 + rep).normal(0, 1 / np.sqrt(G), (len(w.names), G))
        for sc in D:
            rows = E.summarize(w.parts(rho, seeds=(0, 1), n_perm=3, scheme=sc), KS, B=50, others=w.others)
            for q in D[sc]:
                D[sc][q].append(_row(rows, q)["point"])
    z = lambda v: np.mean(v) / (np.std(v, ddof=1) / np.sqrt(len(v)))
    for q in ("D_perm", "D_other"):
        zs, zt = z(D["shared"][q]), z(D["three_way"][q])
        print(f"\n(iii) {q}: shared mean {np.mean(D['shared'][q]):+.4f} (z {zs:+.1f}); "
              f"three-way mean {np.mean(D['three_way'][q]):+.4f} (z {zt:+.1f})")
        assert zs > 4, (q, D["shared"][q])                        # biased upward
        assert abs(zt) < 3, (q, D["three_way"][q])                # bias removed
        assert z(np.subtract(D["shared"][q], D["three_way"][q])) > 3      # paired: shared - three-way > 0


def test_iv_negative_control_noise_only(world):
    parts = world.parts(None, negctrl=True)
    rows = E.summarize(parts, KS, B=500, others=world.others)
    for q in ("N_other", "N_perm"):
        r = _row(rows, q)
        assert abs(r["point"]) < 4 * r["se"], (q, r)
    assert np.isfinite(_row(rows, "D_other")["point"])
    # scaled against a synthetic "observed" run with a planted dest residual, the negative control is < 50 %
    rng = np.random.default_rng(1)
    obs = E.summarize(world.parts(rng.normal(0, 0.3, (len(world.names), Q)) @ world.Qd.T), KS, B=200, others=world.others)
    cmp = E.negctrl_vs_observed(rows, obs, k=8)
    assert abs(cmp["other"]["frac"]) < 0.5 and abs(cmp["perm"]["frac"]) < 0.5, cmp


def test_v_projector_rotation_pcs_split_loader(tmp_path):
    rng = np.random.default_rng(5)
    k = 8
    U = E.haar(G, k, seed=1)
    P = E.projector(U)
    assert np.allclose(P, P.T) and np.allclose(P @ P, P, atol=1e-10) and np.isclose(np.trace(P), k)
    assert np.allclose(U.T @ U, np.eye(k), atol=1e-10)
    rA, rB = rng.standard_normal((50, G)), rng.standard_normal((50, G))
    Qk = np.linalg.qr(rng.standard_normal((k, k)))[0]
    n1, d1 = E.projection_parts(rA, rB, {"u": U}, ks=(k,))
    n2, d2 = E.projection_parts(rA, rB, {"u": U @ Qk}, ks=(k,))
    assert np.allclose(E.F_value(n1["u"], d1), E.F_value(n2["u"], d2))
    assert np.allclose(n1["u"][:, 0], np.einsum("tg,tg->t", rA @ P, rB @ P))
    # randomized SVD subspace == exact SVD subspace on a low-rank + noise matrix
    Y = rng.standard_normal((800, 6)) @ rng.standard_normal((6, G)) * 2 + rng.standard_normal((800, G))
    Ur, _ = E.top_pcs(Y, 6, seed=0)
    Ue = np.linalg.svd(Y - Y.mean(0), full_matrices=False)[2][:6].T
    assert np.linalg.svd(Ur.T @ Ue, compute_uv=False).min() > 0.999
    Up = E.perm_pcs(Y, 6, 2, seed=0)
    assert all(np.allclose(u.T @ u, np.eye(6), atol=1e-8) for u in Up)
    # control split: generated, written, reused identically; obs-row indices converted
    s1, src1 = E.ctrl_split("SYN", 0, 999, write=True, preds_dir=tmp_path)
    s2, src2 = E.ctrl_split("SYN", 0, 999, write=True, preds_dir=tmp_path)
    assert src1.startswith("generated") and src2.startswith("reused")
    assert all(np.array_equal(a, b) for a, b in zip(s1, s2))
    assert sorted(np.concatenate(s1).tolist()) == list(range(999)) and [len(x) for x in s1] == [333] * 3
    rows = np.arange(999) * 3 + 7
    np.savez(tmp_path / "SYN2_ctrl_split_seed1.npz", pc=rows[s1[0]], a=rows[s1[1]], b=rows[s1[2]])
    s3, _ = E.ctrl_split("SYN2", 1, 999, ctrl_rows=rows, preds_dir=tmp_path)
    assert all(np.array_equal(a, b) for a, b in zip(s1, s3))
    # prediction loader: shuffled gene order + extra genes -> aligned to G-cap
    gc = common.g_cap()
    extra = np.array(["NOT_A_GENE_1", "NOT_A_GENE_2"])
    genes = np.concatenate([gc[::-1], extra])
    M = rng.standard_normal((3, len(genes)))
    np.savez(tmp_path / "SYN_prem3.npz", targets=np.array(["a", "b", "c"]), genes=genes, delta=M)
    tg, D, sha = E.load_pred("SYN", "prem3", gc, preds_dir=tmp_path)
    assert list(tg) == ["a", "b", "c"] and D.shape == (3, len(gc)) and np.allclose(D, M[:, : len(gc)][:, ::-1])


H1_TRUTH = common.truth_path("H1")


@pytest.mark.skipif(not H1_TRUTH.exists(), reason="H1 truth not on disk")
def test_smoke_h1_controls_subspaces():
    t0 = time.time()
    C, genes = common.controls("H1")
    L = np.asarray(C.sum(1), float).ravel()
    og = {g: i for i, g in enumerate(common.official_genes())}
    om = common.to_official(genes); first = {}
    for c, o in enumerate(om):
        first.setdefault(o, c)
    gc = common.g_cap()
    cols = np.array([first[og[g]] for g in gc])
    assert len(cols) == 6014
    X = C[:, cols].tocsr()
    split, src = E.ctrl_split("H1", 0, X.shape[0], write=False)
    pc = split[0][L[split[0]] > 0]
    Y = E.cell_log(X[pc], L[pc])
    Ud, s = E.top_pcs(Y, 8, seed=0)
    Up = E.perm_pcs(Y, 8, n_perm=2, seed=0)
    for U in [Ud] + Up:
        assert U.shape == (6014, 8) and np.allclose(U.T @ U, np.eye(8), atol=1e-8)
    # U_dest captures more held-out (C_a) control variance than the permutation null
    Ya = E.cell_log(X[split[1]], L[split[1]]); Ya -= Ya.mean(0)
    ve = lambda U: float(((Ya @ U) ** 2).sum() / (Ya ** 2).sum())
    assert ve(Ud) > ve(Up[0]) and ve(Ud) > 8 / 6014
    print(f"\nH1 smoke: {X.shape[0]} controls, C_pc {len(pc)}, var-explained dest {ve(Ud):.4f} perm {ve(Up[0]):.4f}, "
          f"split {src}, {time.time() - t0:.1f}s")


# ---------------- Task 16 additions: summarize2 / two-level bootstrap ----------------
def _stores(w, rho, seeds=(0, 1, 2), n_perm=5, B=200, two_level=True):
    out = []
    for s in seeds:
        split = E.common.split_halves(len(w.cL), s, k=3)
        subs = E.subspaces_for_seed(E.cell_log(w.cX[split[0]], w.cL[split[0]]), w.U_other, s, max(KS), n_perm)
        hs = [E.target_half_sums(X, L, s, n) for X, L, n in zip(w.tX, w.tL, w.names)]
        rA, rB = E.residuals(hs, w.cX, w.cL, split, w.dexp - rho)
        dA = dB = None
        if two_level:
            rng = np.random.default_rng([E.BOOT_SEED, s, 17])
            dA = E.ctrl_boot_deltas(w.cX, w.cL, split[1], B, rng); dB = E.ctrl_boot_deltas(w.cX, w.cL, split[2], B, rng)
        out.append(E.boot_store(rA, rB, subs, dA, dB, max(KS)))
    return out


def test_vi_summarize2_reproduces_summarize_and_point(world):
    rng = np.random.default_rng(1)
    rho = rng.normal(0, 0.3, (len(world.names), Q)) @ world.Qd.T
    ref = E.summarize(world.parts(rho), KS, B=200, others=world.others)
    st = _stores(world, rho, B=200)
    t1 = E.summarize2(st, KS, B=200)
    t2 = E.summarize2(st, KS, B=200, two_level=True)
    idx = np.arange(0, len(world.names), 3)
    r_ref = E.summarize(world.parts(rho), KS, B=200, idx=idx, others=world.others)
    r_new = E.summarize2(st, KS, B=200, idx=idx)
    for rows_a, rows_b in ((ref, t1), (r_ref, r_new)):
        for r in rows_a:
            m = next(x for x in rows_b if x["quantity"] == r["quantity"] and x["k"] == r["k"])
            for f in ("point", "lo", "hi", "se"):
                assert np.isclose(r[f], m[f], rtol=1e-8, atol=1e-12), (r["quantity"], r["k"], f, r[f], m[f])
    for r in t1:   # two-level row 0 (no resampling) == point estimate
        m = next(x for x in t2 if x["quantity"] == r["quantity"] and x["k"] == r["k"])
        assert np.isclose(r["point"], m["point"], rtol=1e-10, atol=1e-14)


def test_vii_two_level_wider_with_small_control_thirds():
    """Isotropic residual (true D = 0), 100-cell control thirds: the two-level SE of D_perm / D_other exceeds the
    target-only SE, and over R replicates its z-scores are closer to unit spread."""
    R, zt, z2, ratio = 8, {"D_perm": [], "D_other": []}, {"D_perm": [], "D_other": []}, []
    for rep in range(R):
        w = World(seed=300 + rep, n_ctrl=300, T=100)
        rho = np.random.default_rng(400 + rep).normal(0, 1 / np.sqrt(G), (len(w.names), G))
        st = _stores(w, rho, seeds=(0, 1), n_perm=3, B=200)
        a, b = E.summarize2(st, KS, B=200), E.summarize2(st, KS, B=200, two_level=True)
        for q in zt:
            ra, rb = _row(a, q), _row(b, q)
            zt[q].append(ra["point"] / ra["se"]); z2[q].append(rb["point"] / rb["se"]); ratio.append(rb["se"] / ra["se"])
    print(f"\n(vii) se ratio two-level/target median {np.median(ratio):.2f}; sd(z) target "
          f"{np.std(zt['D_perm'] + zt['D_other']):.2f}, two-level {np.std(z2['D_perm'] + z2['D_other']):.2f}")
    assert np.median(ratio) > 1.3
    assert np.std(z2["D_perm"] + z2["D_other"]) < np.std(zt["D_perm"] + zt["D_other"])
