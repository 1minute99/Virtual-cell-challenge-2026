"""RQ3 Task 15 -- split-half residual-energy estimator (operational definition 7; H3a). Library + CLI.

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md -- Task 15 (this file + tests), Task 16 (the real run,
which is NOT done here), H3a and its falsification criterion, operational definition 7. No amendment (A1-A3) touches
RQ3; A3 fixes the truth paths (results/panels.json, used through common.py).

ESTIMATOR (definition 7, unchanged)
  Control cells of the destination panel are split three ways (equal thirds): C_pc (PC estimation), C_a (half-A
  baseline), C_b (half-B baseline). Each target's cells are split at random into halves A / B. In log1p CP10K
  pseudobulk on G-cap (6,014 genes):
      dtruth_A,t = pb(A_t) - pb(C_a),  dtruth_B,t = pb(B_t) - pb(C_b),
      r_A,t = dtruth_A,t - dpred_t,    r_B,t = dtruth_B,t - dpred_t,
      F(U)  = sum_t <P_U r_A,t, P_U r_B,t> / sum_t <r_A,t, r_B,t>,   P_U = U U^T (U orthonormal, G x k).
  A and B use disjoint cells and disjoint baselines, and U_dest is estimated on C_pc, disjoint from both, so the
  noise of r_A and r_B is independent of each other and of U: E[numerator] = sum_t <P_U rho_t, P_U rho_t>, where rho_t
  is the noise-free residual (up to the O(1/n) log-pseudobulk bias, assumption A6).
  Subspaces (all on G-cap, k in {8, 24, 50}; computed once at k_max = 50 and sliced, PCs are nested):
      U_dest   top-k PCs of C_pc (per-cell log1p CP10K, centred across cells, no scaling)
      U_j      top-k PCs of ALL control cells of each other harness panel j (6 for every destination)
      U_perm   top-k PCs of C_pc after an independent per-gene permutation across cells; 20 permutations, F averaged
      U_haar   uniform random k-frame (Gaussian QR) -- reported only
  D_other = F(U_dest) - mean_j F(U_j);  D_perm = F(U_dest) - mean_p F(U_perm,p).
  Each D is averaged over the 5 split seeds (0-4); uncertainty = target bootstrap (B = 2,000; the same resampled
  target set is used for all seeds in a draw, then the seed mean is taken), percentile 95 % CI.
  Strata: quartiles of the predicted amplitude ||dpred_t|| on G-cap (test-time quantity; Q1 = smallest).
  Negative-control mode: 300 pseudo-targets per seed drawn without replacement from C_a u C_b (destination controls
  excluding C_pc), sizes drawn from the panel's real per-target cell counts (capped at half the pool), dpred = 0;
  each pseudo-target's baselines are C_a minus its cells and C_b minus its cells (keeps A, B, C_a', C_b' disjoint).

NOTE ON THE NEGATIVE-CONTROL RATIO (spec ambiguity, flagged to the parent; not resolved here). With pure-control
  pseudo-targets the noise-free residual is 0, so BOTH the numerator and the denominator of F have expectation ~0 and
  the registered ratio-form D is a ratio of two noise terms (unstable, sign-arbitrary). This module therefore reports,
  for both modes, (a) the registered ratio-form D, and (b) the numerator-unit difference per target
  N_other = mean_s [sum_t (num_dest - mean_j num_j) / T] (and N_perm), plus the per-target denominator den/T.
  `negctrl_vs_observed` gives D_neg,scaled = N_neg / (den_obs / T_obs), i.e. the negative-control numerator
  excess expressed on the observed D's scale; Task 16 decides (with the parent) which form enters the H3a 50 % rule.

ASSUMPTIONS (stated; the testable ones are checked in tests/test_rq3_energy.py and the note
knowledge/synthesis/rq3_energy-assumptions.md)
  A1 Truth = the harness truth file (common.py, which="truth"; panels.json). Controls = obs[label] == "non-targeting".
     Control-split indices refer to positions in the file-order list of control rows of that file.
  A2 Pseudobulk pb(S) = log1p(1e4 * sum_{c in S} x_cg / sum_{c in S} L_c), L_c = cell's total UMI over ALL genes of
     its panel (not only G-cap) -- the same rule as common.pseudobulk (sum, then normalise, then log1p).
     Per-cell PCs use log1p(1e4 * x_cg / L_c) on G-cap. Cells with L_c = 0 are dropped.
  A3 G-cap = results/g_cap.json (6,014 official-axis symbols); a panel column is matched through
     common.to_official; if two panel columns map to one official gene, the first is used (1 case on H1).
  A4 dpred (Task 14) is a target x gene pseudobulk delta in log1p CP10K, i.e. pb(pred target) - pb(ctrl), and is
     aligned to G-cap by official symbol (Task 14 writes the truth file's var_names as `genes`, CP10K over all panel
     genes -- consistent with A2). It depends on the destination control MEAN (Task 14 A4); `--ctrl-ab` swaps in
     Task 14's per-seed C_a u C_b rebuild so that dpred is independent of C_pc. Expected file: results/preds/<panel>_{prem3,m3}.npz with keys
     `targets` (T,), `genes` (G,) [optional if G = 18,533 official axis] and one of
     `delta` / `pred_delta` / `dpred` / `D` (T x G). H3a uses prem3.
  A5 Three-way control split rule (used only if results/preds/<panel>_ctrl_split_seed<k>.npz is absent):
     (C_pc, C_a, C_b) = common.split_halves(n_ctrl, seed=k, k=3). Files written by this module carry keys
     C_pc, C_a, C_b, n_ctrl, seed, which, rule. An existing file is reused (keys C_pc/C_a/C_b or pc/a/b or
     idx_pc/idx_a/idx_b; obs-row indices are converted to control positions) and validated (disjoint, in range).
  A6 Log-pseudobulk bias is O(1/n) and common to both halves; at the harness cell counts it is negligible relative to
     the noise of the cross-product (checked on synthetic Poisson data, test (i)).
  A7 Target-level exchangeability for the bootstrap (targets are the resampling unit); split seeds are not
     independent replicates and are averaged, not pooled into the CI.
  A8 Target halves: rng = default_rng([split_seed, 7, crc32(target)]); targets with < 2 cells are dropped.
  A9 PCs by randomized SVD (sklearn, n_oversamples 20, n_iter 7, fixed random_state); checked against exact SVD in
     the tests. The permutation null permutes each gene's column of C_pc independently (rng.permuted, axis 0).
Seeds: split seeds 0-4; PCs random_state 0; permutations default_rng([split_seed, 11]); Haar default_rng([split_seed,
  13]); bootstrap default_rng(20261105); negative-control default_rng([split_seed, 99]).
Environment: ~/miniconda3/envs/ct_transfer; run under nice -n 10, OMP/MKL threads 4.
TASK 16 ADDITIONS (2026-11-05; estimator unchanged; see the log clarification "Pre-data operationalisation for Task 16")
  T16-1 Two-level bootstrap (targets x control cells): each draw b resamples targets (the SAME W as the registered
     target bootstrap) AND, independently per split seed, resamples the cells of C_a and C_b with replacement
     (rng default_rng([20261105, split_seed, 17])). The resampled baseline enters as r*_A,t = r_A,t - dA_b with
     dA_b = pb(C_a*) - pb(C_a); numerator and denominator are expanded exactly (summarize2). C_pc (hence U_dest,
     U_perm) is NOT resampled -- PC-estimation noise is outside both CIs (it varies across the averaged split seeds).
     Row 0 (no resampling) reproduces the point estimate exactly; two_level=False reproduces `summarize` exactly.
  T16-2 space='full': secondary destination-only analysis on the destination's FULL gene axis (truth var_names, all
     columns; Task 14 writes dpred on that axis): U_dest, U_perm, U_haar only (no U_j, no D_other).
  T16-3 load_pred uses allow_pickle=True: Task 14 stores `genes` as an object array (own trusted files).
  T16-4 Negative control: target (pseudo-target) bootstrap only; two-level not defined there (each pseudo-target has
     its own leave-out baselines C_a minus S, C_b minus S).
CLI (Task 16): python rq3_energy.py --panel H1 --mode real --pred prem3   (writes results/rq3/<mode>_<panel>_<pred>.*)
               python rq3_energy.py --collect                              (-> results/rq3_energy.parquet, rq3_negctrl.parquet)
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import argparse, hashlib, json, pathlib, sys, time, zlib
import numpy as np
from scipy import sparse

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402

PREDS = common.RESULTS / "preds"
OUT = common.RESULTS / "rq3"
CACHE = common.RESULTS / "rq3_cache"
KS = (8, 24, 50)
KMAX = 50
SEEDS = (0, 1, 2, 3, 4)
N_PERM = 20
B_BOOT = 2000
BOOT_SEED = 20261105
N_PSEUDO = 300


# ============================== core (pure, tested on synthetic data) ==============================
def cell_log(X, L):
    """Per-cell log1p CP10K. X: (n x G) counts (csr or dense) on the analysis genes; L: (n,) library sizes."""
    s = (1e4 / np.asarray(L, float))
    if sparse.issparse(X):
        Y = sparse.diags(s) @ X.astype(np.float64)
        return np.log1p(Y.toarray())
    return np.log1p(np.asarray(X, float) * s[:, None])


def pb(sums, lib):
    """Pseudobulk log1p CP10K from summed counts (.. x G) and summed library size (..)."""
    sums = np.asarray(sums, float); lib = np.asarray(lib, float)
    return np.log1p(sums * (1e4 / np.maximum(lib, 1e-12))[..., None]) if sums.ndim > 1 else np.log1p(sums * 1e4 / max(lib, 1e-12))


def top_pcs(Y, kmax=KMAX, seed=0):
    """Top-kmax right singular vectors of the column-centred Y (n x G) -> (G x kmax) orthonormal, singular values."""
    from sklearn.utils.extmath import randomized_svd
    Yc = np.asarray(Y, float) - np.asarray(Y, float).mean(0, keepdims=True)
    k = min(kmax, min(Yc.shape) - 1)
    _, s, Vt = randomized_svd(Yc, n_components=k, n_oversamples=20, n_iter=7, random_state=seed)
    return np.ascontiguousarray(Vt.T), s


def perm_pcs(Y, kmax=KMAX, n_perm=N_PERM, seed=0):
    """Top PCs of Y after an independent per-gene permutation across cells, n_perm times."""
    rng = np.random.default_rng([seed, 11])
    return [top_pcs(rng.permuted(Y, axis=0), kmax, seed=0)[0] for _ in range(n_perm)]


def haar(G, kmax=KMAX, seed=0):
    """Uniform random orthonormal G x kmax frame (QR of a Gaussian matrix, sign-fixed)."""
    Q, R = np.linalg.qr(np.random.default_rng([seed, 13]).standard_normal((G, kmax)))
    return Q * np.sign(np.diag(R))


def projector(U):
    return U @ U.T


def projection_parts(rA, rB, subspaces: dict, ks=KS):
    """Per-target F-numerator contributions <P_U r_A,t, P_U r_B,t> for nested k, and the denominator <r_A,t, r_B,t>.
    subspaces: name -> (G x kmax) orthonormal. Returns ({name: (T x len(ks))}, den (T,))."""
    rA = np.asarray(rA, float); rB = np.asarray(rB, float)
    den = np.einsum("tg,tg->t", rA, rB)
    num = {}
    for name, U in subspaces.items():
        cs = np.cumsum((rA @ U) * (rB @ U), axis=1)      # (T x kmax), nested-k cumulative energy
        num[name] = np.stack([cs[:, min(k, U.shape[1]) - 1] for k in ks], axis=1)
    return num, den


def F_value(num_t, den_t, idx=None):
    idx = slice(None) if idx is None else idx
    return np.asarray(num_t)[idx].sum(0) / np.asarray(den_t)[idx].sum()


def group_parts(num: dict, others, perms):
    """Collapse named parts into arrays: dest (T,K), other (T,J,K), perm (T,P,K), haar (T,K)."""
    return dict(dest=num["dest"], other=np.stack([num[f"other:{j}"] for j in others], 1),
                perm=np.stack([num[f"perm:{p}"] for p in range(perms)], 1), haar=num["haar"])


def summarize(seed_parts: list, ks=KS, B=B_BOOT, boot_seed=BOOT_SEED, idx=None, others=None):
    """seed_parts: list over split seeds of dict(dest (T,K), other (T,J,K), perm (T,P,K), haar (T,K), den (T,)).
    Same T (target order) in every seed. Returns list of row dicts (quantity, k, point, boot stats)."""
    S = len(seed_parts); T = seed_parts[0]["den"].shape[0]
    tid = np.arange(T) if idx is None else np.asarray(idx)
    n = len(tid)
    rng = np.random.default_rng(boot_seed)
    W = np.zeros((B + 1, T))
    W[0, tid] = 1.0                                           # row 0 = point estimate
    draws = rng.integers(0, n, size=(B, n))
    for b in range(B):
        W[b + 1] = np.bincount(tid[draws[b]], minlength=T)
    stack = lambda key: np.stack([sp[key] for sp in seed_parts])          # (S, T, ...)
    den = W @ stack("den").T                                              # (B+1, S)
    sd = np.einsum("bt,stk->bsk", W, stack("dest"))
    so = np.einsum("bt,stjk->bsjk", W, stack("other"))
    sp_ = np.einsum("bt,stpk->bspk", W, stack("perm"))
    sh = np.einsum("bt,stk->bsk", W, stack("haar"))
    Fd = sd / den[..., None]; Fo = so / den[..., None, None]; Fp = sp_ / den[..., None, None]; Fh = sh / den[..., None]
    q = {"F_dest": Fd.mean(1), "F_other_mean": Fo.mean(2).mean(1), "F_perm_mean": Fp.mean(2).mean(1),
         "F_haar": Fh.mean(1),
         "D_other": (Fd - Fo.mean(2)).mean(1), "D_perm": (Fd - Fp.mean(2)).mean(1),
         # numerator units per target (used by the negative control; see the module docstring)
         "N_other": ((sd - so.mean(2)) / n).mean(1), "N_perm": ((sd - sp_.mean(2)) / n).mean(1),
         "den_per_target": np.repeat((den / n).mean(1)[:, None], len(ks), 1)}
    if others is not None:
        for jj, j in enumerate(others):
            q[f"F_other:{j}"] = Fo[:, :, jj, :].mean(1)
    rows = []
    for name, v in q.items():                                  # v: (B+1, K)
        for kk, k in enumerate(ks):
            bs = v[1:, kk]
            rows.append(dict(quantity=name, k=k, point=float(v[0, kk]), boot_mean=float(bs.mean()),
                             se=float(bs.std(ddof=1)), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)),
                             frac_le0=float((bs <= 0).mean()), n_targets=int(n), n_seeds=S, B=B,
                             den_sign_ok=bool((den[0] > 0).all())))
    return rows


def negctrl_vs_observed(neg_rows, obs_rows, k=8):
    """D_neg,scaled = N_neg / den_obs_per_target, reported as a fraction of the observed ratio-form D."""
    g = lambda rows, qn: next(r for r in rows if r["quantity"] == qn and r["k"] == k)
    dpt = g(obs_rows, "den_per_target")["point"]
    out = {}
    for which in ("other", "perm"):
        dneg = g(neg_rows, f"N_{which}")["point"] / dpt
        dobs = g(obs_rows, f"D_{which}")["point"]
        out[which] = dict(D_neg_scaled=dneg, D_obs=dobs, frac=dneg / dobs if dobs != 0 else np.nan,
                          D_neg_ratio_form=g(neg_rows, f"D_{which}")["point"])
    return out


# ---------------- deltas ----------------
def target_half_sums(X, L, split_seed, name):
    """Random cell halves of one target -> (sum_A (G,), lib_A, sum_B (G,), lib_B)."""
    n = X.shape[0]
    p = np.random.default_rng([split_seed, 7, zlib.crc32(str(name).encode())]).permutation(n)
    a, b = np.sort(p[: n // 2]), np.sort(p[n // 2:])
    s = lambda ii: np.asarray(X[ii].sum(0), float).ravel()
    return s(a), float(np.sum(L[a])), s(b), float(np.sum(L[b]))


def group_sum(X, L, ii):
    return np.asarray(X[ii].sum(0), float).ravel(), float(np.sum(np.asarray(L)[ii]))


def residuals(half_sums, ctrl_X, ctrl_L, split, dpred, scheme="three_way"):
    """half_sums: list of (sA, lA, sB, lB) per target. scheme 'three_way' (registered): baselines C_a / C_b;
    'shared' (biased reference variant for test iii): one baseline C_a u C_b for both halves."""
    C_pc, C_a, C_b = split
    if scheme == "three_way":
        bA = pb(*group_sum(ctrl_X, ctrl_L, C_a)); bB = pb(*group_sum(ctrl_X, ctrl_L, C_b))
    elif scheme == "shared":
        bA = bB = pb(*group_sum(ctrl_X, ctrl_L, np.concatenate([C_a, C_b])))
    else:
        raise ValueError(scheme)
    sA = np.stack([h[0] for h in half_sums]); lA = np.array([h[1] for h in half_sums])
    sB = np.stack([h[2] for h in half_sums]); lB = np.array([h[3] for h in half_sums])
    return pb(sA, lA) - bA - dpred, pb(sB, lB) - bB - dpred


def negctrl_residuals(ctrl_X, ctrl_L, split, sizes, split_seed, n_pseudo=N_PSEUDO):
    """Pseudo-targets from C_a u C_b, dpred = 0, baselines C_a \\ S and C_b \\ S. Returns rA, rB, used sizes."""
    _, C_a, C_b = split
    pool = np.concatenate([C_a, C_b]); in_a = np.isin(pool, C_a)
    rng = np.random.default_rng([split_seed, 99])
    sz = np.minimum(rng.choice(np.asarray(sizes), n_pseudo, replace=True), len(pool) // 2)
    sa_full, la_full = group_sum(ctrl_X, ctrl_L, C_a); sb_full, lb_full = group_sum(ctrl_X, ctrl_L, C_b)
    rA, rB = [], []
    for i in range(n_pseudo):
        pick = rng.choice(len(pool), sz[i], replace=False)
        S = pool[pick]; Sa, Sb = S[in_a[pick]], S[~in_a[pick]]
        h = target_half_sums(ctrl_X[np.sort(S)], np.asarray(ctrl_L)[np.sort(S)], split_seed, f"pseudo{i}")
        xa, la = group_sum(ctrl_X, ctrl_L, Sa) if len(Sa) else (0.0, 0.0)
        xb, lb = group_sum(ctrl_X, ctrl_L, Sb) if len(Sb) else (0.0, 0.0)
        bA = pb(sa_full - xa, la_full - la); bB = pb(sb_full - xb, lb_full - lb)
        rA.append(pb(h[0], h[1]) - bA); rB.append(pb(h[2], h[3]) - bB)
    return np.stack(rA), np.stack(rB), sz


def subspaces_for_seed(Y_pc, U_other: dict, split_seed, kmax=KMAX, n_perm=N_PERM):
    subs = {"dest": top_pcs(Y_pc, kmax, seed=0)[0]}
    for j, U in U_other.items():
        subs[f"other:{j}"] = U[:, :kmax]
    for p, U in enumerate(perm_pcs(Y_pc, kmax, n_perm, seed=split_seed)):
        subs[f"perm:{p}"] = U
    subs["haar"] = haar(Y_pc.shape[1], kmax, seed=split_seed)
    return subs


def seed_parts_from_residuals(rA, rB, subs, others, n_perm, ks=KS):
    num, den = projection_parts(rA, rB, subs, ks)
    d = group_parts(num, others, n_perm); d["den"] = den
    return d



# ============================== Task 16: two-level bootstrap (targets x control cells) ==============================
def ctrl_boot_deltas(ctrl_X, ctrl_L, idx, B, rng, chunk=250):
    """dA_b = pb(bootstrap resample of control cells idx) - pb(idx), b = 1..B -> (B, G)."""
    idx = np.asarray(idx); n = len(idx)
    Xs = ctrl_X[idx]; Ls = np.asarray(ctrl_L, float)[idx]
    base = pb(np.asarray(Xs.sum(0), float).ravel(), Ls.sum())
    XsT = (Xs.T.tocsr() if sparse.issparse(Xs) else np.asarray(Xs).T)
    out = np.empty((B, Xs.shape[1]))
    for b0 in range(0, B, chunk):
        nb = min(chunk, B - b0)
        M = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(nb)]).astype(np.float64)
        S = np.asarray(XsT @ M.T).T
        out[b0:b0 + nb] = pb(S, M @ Ls) - base
    return out


def boot_store(rA, rB, subs, dA=None, dB=None, kmax=KMAX):
    """Per-target projections a = U^T r_A, b = U^T r_B for every subspace (T, U, kmax), the per-target denominator,
    and (optional) the control-resampling terms padded with a zero row 0 (= no resampling) for summarize2."""
    names = list(subs)
    Ucat = np.concatenate([subs[n][:, :kmax] for n in names], 1)          # (G, U*kmax)
    T = rA.shape[0]; nU = len(names)
    st = dict(names=names, a=(rA @ Ucat).reshape(T, nU, kmax), b=(rB @ Ucat).reshape(T, nU, kmax),
              den_t=np.einsum("tg,tg->t", rA, rB))
    if dA is not None:
        z = lambda M: np.concatenate([np.zeros((1,) + M.shape[1:]), M])
        st.update(x=z((dA @ Ucat).reshape(len(dA), nU, kmax)), y=z((dB @ Ucat).reshape(len(dB), nU, kmax)),
                  PA=z(dA @ rB.T), PB=z(dB @ rA.T), dd=z(np.einsum("bg,bg->b", dA, dB)))
    return st


def summarize2(stores: list, ks=KS, B=B_BOOT, boot_seed=BOOT_SEED, idx=None, two_level=False):
    """Same quantities as `summarize`, from boot_store dicts (one per split seed). two_level=False reproduces
    `summarize` exactly (same W); two_level=True adds the control-cell resampling of C_a / C_b (row b of x/y/PA/PB/dd
    belongs to draw b). Subspace groups are read from the names (dest, other:<j>, perm:<p>, haar)."""
    S = len(stores); T = stores[0]["den_t"].shape[0]
    tid = np.arange(T) if idx is None else np.asarray(idx); n = len(tid)
    rng = np.random.default_rng(boot_seed)
    W = np.zeros((B + 1, T)); W[0, tid] = 1.0
    draws = rng.integers(0, n, size=(B, n))
    for b in range(B):
        W[b + 1] = np.bincount(tid[draws[b]], minlength=T)
    names = stores[0]["names"]
    i_d = names.index("dest"); i_h = names.index("haar")
    others = [x.split(":", 1)[1] for x in names if x.startswith("other:")]
    i_o = [names.index(f"other:{j}") for j in others]
    i_p = [i for i, x in enumerate(names) if x.startswith("perm:")]
    kk = [k - 1 for k in ks]
    NUM, DEN = [], []
    for st in stores:
        a, b = st["a"], st["b"]; U, K = a.shape[1], a.shape[2]
        num = (W @ (a * b).reshape(T, -1)).reshape(B + 1, U, K)
        den = W @ st["den_t"]
        if two_level:
            if "x" not in st:
                raise ValueError("two_level needs control-resampling terms (boot_store with dA, dB)")
            Sa = (W @ a.reshape(T, -1)).reshape(B + 1, U, K); Sb = (W @ b.reshape(T, -1)).reshape(B + 1, U, K)
            x, y = st["x"], st["y"]
            num = num - x * Sb - y * Sa + n * x * y
            den = den - (W * st["PA"]).sum(1) - (W * st["PB"]).sum(1) + n * st["dd"]
        NUM.append(np.cumsum(num, axis=2)[:, :, kk]); DEN.append(den)
    num = np.stack(NUM, 1); den = np.stack(DEN, 1)                         # (B+1, S, U, nk), (B+1, S)
    F = num / den[:, :, None, None]
    Fd, Fh, Fp = F[:, :, i_d], F[:, :, i_h], F[:, :, i_p].mean(2)
    nd, npm = num[:, :, i_d], num[:, :, i_p].mean(2)
    q = {"F_dest": Fd.mean(1), "F_perm_mean": Fp.mean(1), "F_haar": Fh.mean(1), "D_perm": (Fd - Fp).mean(1),
         "N_perm": ((nd - npm) / n).mean(1)}
    if i_o:
        Fo = F[:, :, i_o].mean(2); no = num[:, :, i_o].mean(2)
        q.update({"F_other_mean": Fo.mean(1), "D_other": (Fd - Fo).mean(1), "N_other": ((nd - no) / n).mean(1)})
        for jj, j in zip(i_o, others):
            q[f"F_other:{j}"] = F[:, :, jj].mean(1)
    q["den_per_target"] = np.repeat((den / n).mean(1)[:, None], len(ks), 1)
    rows = []
    for name, v in q.items():
        for c, k in enumerate(ks):
            bs = v[1:, c]
            rows.append(dict(quantity=name, k=k, point=float(v[0, c]), boot_mean=float(bs.mean()),
                             se=float(bs.std(ddof=1)), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)),
                             frac_le0=float((bs <= 0).mean()), n_targets=int(n), n_seeds=S, B=B,
                             den_sign_ok=bool((den[0] > 0).all()), den_boot_frac_pos=float((den[1:].mean(1) > 0).mean())))
    return rows


# ============================== real-data I/O ==============================
class PanelData:
    """Backed access to one panel's truth on G-cap. read(idx) -> (csr n x G-cap counts, (n,) library size)."""

    def __init__(self, panel, which="truth", space="gcap"):
        self.panel, self.which, self.space = panel, which, space
        self.a = common.open_panel(panel, which)
        col, ctrl = common._label_col(panel, which)
        self.lab = self.a.obs[col].astype(str).to_numpy()
        self.ctrl_rows = np.flatnonzero(self.lab == ctrl)
        self.ctrl_label = ctrl
        self.gcap = common.g_cap()
        self.var_names = self.a.var_names.astype(str).to_numpy()
        if space == "full":                                      # T16-2
            self.cols = np.arange(len(self.var_names)); self.genes = self.var_names
            return
        self.genes = self.gcap
        og = {g: i for i, g in enumerate(common.official_genes())}
        want = np.array([og[g] for g in self.gcap])
        om = common.to_official(self.a.var_names.astype(str).to_numpy())
        first = {}
        for c, o in enumerate(om):
            if o >= 0 and o not in first:
                first[o] = c
        missing = [g for g, o in zip(self.gcap, want) if o not in first]
        if missing:
            raise ValueError(f"{panel}: {len(missing)} G-cap genes absent, e.g. {missing[:5]}")
        self.cols = np.array([first[o] for o in want])

    def read(self, idx):
        X = common._rows(self.a, idx)                           # sorts idx
        L = np.asarray(X.sum(1), float).ravel()
        return X[:, self.cols].tocsr(), L

    def targets(self):
        return sorted(set(self.lab) - {self.ctrl_label})

    def close(self):
        self.a.file.close()


def ctrl_split(panel, seed, n_ctrl, which="truth", write=True, ctrl_rows=None, preds_dir=None):
    """Three-way control split (positions in the file-order control list). Reuse <panel>_ctrl_split_seed<k>.npz."""
    d = pathlib.Path(preds_dir) if preds_dir else PREDS
    f = d / f"{panel}_ctrl_split_seed{seed}.npz"
    if f.exists():
        z = np.load(f, allow_pickle=False)
        names = next(t for t in (("C_pc", "C_a", "C_b"), ("pc", "a", "b"), ("idx_pc", "idx_a", "idx_b")) if all(x in z for x in t))
        parts = [np.asarray(z[x]).astype(np.int64) for x in names]
        if max(p.max() for p in parts) >= n_ctrl:                # obs-row indices -> control positions
            if ctrl_rows is None:
                raise ValueError(f"{f}: indices exceed n_ctrl={n_ctrl} and no ctrl_rows given")
            pos = {r: i for i, r in enumerate(ctrl_rows)}
            parts = [np.array([pos[r] for r in p]) for p in parts]
        source = f"reused {f.name}"
    else:
        parts = common.split_halves(n_ctrl, seed, k=3)
        source = "generated common.split_halves(n_ctrl, seed, k=3)"
        if write:
            d.mkdir(parents=True, exist_ok=True)
            tmp = f.with_name(f.name + f".tmp{os.getpid()}.npz")
            np.savez(tmp, C_pc=parts[0], C_a=parts[1], C_b=parts[2], n_ctrl=n_ctrl, seed=seed, which=which,
                     rule="(C_pc, C_a, C_b) = common.split_halves(n_ctrl, seed, k=3); positions in file-order control rows")
            if not f.exists():
                os.replace(tmp, f)
            else:
                os.remove(tmp)
                return ctrl_split(panel, seed, n_ctrl, which, False, ctrl_rows, preds_dir)
    allc = np.concatenate(parts)
    assert len(np.unique(allc)) == len(allc), f"{panel} seed {seed}: control split not disjoint"
    assert allc.min() >= 0 and allc.max() < n_ctrl
    return [np.sort(p) for p in parts], source


def load_pred(panel, kind, gcap, preds_dir=None, full_genes=None):
    """Task 14 prediction -> (targets (T,), dpred (T x |G-cap|), sha256). full_genes: the destination's var_names
    (space='full'); then the prediction's gene axis must equal it and is returned unchanged."""
    f = (pathlib.Path(preds_dir) if preds_dir else PREDS) / f"{panel}_{kind}.npz"
    z = np.load(f, allow_pickle=True)                         # T16-3: `genes` is an object array in Task 14 files
    if full_genes is not None:
        key = next(k for k in ("delta", "pred_delta", "dpred", "D") if k in z)
        g = z["genes"].astype(str)
        assert list(g) == list(full_genes), f"{f}: gene axis differs from the truth var_names"
        return z["targets"].astype(str), np.asarray(z[key], float), hashlib.sha256(f.read_bytes()).hexdigest()
    key = next(k for k in ("delta", "pred_delta", "dpred", "D") if k in z)
    M = np.asarray(z[key], float)
    genes = z["genes"].astype(str) if "genes" in z else (common.official_genes() if M.shape[1] == len(common.official_genes()) else None)
    if genes is None:
        raise ValueError(f"{f}: no `genes` and width {M.shape[1]} is not the official axis")
    og = {g: i for i, g in enumerate(common.official_genes())}
    om = common.to_official(genes)
    col = {}
    for c, o in enumerate(om):
        col.setdefault(o, c)
    want = [og[g] for g in gcap]
    miss = [g for g, o in zip(gcap, want) if o not in col]
    if miss:
        raise ValueError(f"{f}: {len(miss)} G-cap genes missing from the prediction, e.g. {miss[:5]}")
    sha = hashlib.sha256(f.read_bytes()).hexdigest()
    return z["targets"].astype(str), M[:, [col[o] for o in want]], sha


def other_panel_pcs(panel, which="truth", kmax=KMAX):
    """Top-kmax control PCs of ALL control cells of `panel` on G-cap (cached)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"ctrlpcs_{panel}_{which}_k{kmax}.npz"
    gc = common.g_cap()
    if f.exists():
        z = np.load(f, allow_pickle=False)
        if list(z["genes"].astype(str)) == list(gc):
            return z["U"]
    P = PanelData(panel, which)
    try:
        X, L = P.read(P.ctrl_rows)
    finally:
        P.close()
    keep = L > 0
    U, s = top_pcs(cell_log(X[keep], L[keep]), kmax, seed=0)
    np.savez(f, U=U, s=s, genes=gc, n_cells=int(keep.sum()))
    return U


def run(panel, mode="real", pred="prem3", seeds=SEEDS, ks=KS, which="truth", n_perm=N_PERM, B=B_BOOT,
        write_split=True, out_dir=None, log=print, ctrl_ab=False, space="gcap", two_level=True):
    """Task 16 entry point. mode 'real' (needs results/preds/<panel>_<pred>.npz) or 'negctrl'.
    ctrl_ab=True: per split seed use Task 14's <panel>_<pred>_ctrlAB_seed<s>.npz (pre-M3 rebuilt with every
    control-derived quantity from C_a u C_b only, i.e. independent of C_pc); strata use the main file's amplitude.
    space 'gcap' (registered: G-cap, U_dest/U_j/U_perm/U_haar) or 'full' (T16-2: destination's full gene axis, no U_j).
    Rows carry boot = 'target' (registered target bootstrap) and, in real mode, 'two_level' (T16-1)."""
    t0 = time.time()
    out_dir = pathlib.Path(out_dir) if out_dir else OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    others = [p for p in common.PANELS7 if p != panel] if space == "gcap" else []
    U_other = {j: other_panel_pcs(j, which) for j in others}
    log(f"[{panel}] other-panel PCs ready ({len(U_other)}) {time.time()-t0:.0f}s")
    P = PanelData(panel, which, space)
    tag = (pred + ("_ctrlAB" if ctrl_ab else "")) if mode == "real" else "none"
    meta = dict(panel=panel, mode=mode, pred=tag, which=which, space=space, n_genes=int(len(P.cols)), seeds=list(seeds),
                ks=list(ks), n_perm=n_perm, B=B, ctrl_ab=ctrl_ab, two_level=bool(two_level and mode == "real"))
    fg = P.genes if space == "full" else None
    try:
        cX, cL = P.read(P.ctrl_rows)
        if (cL <= 0).any():
            log(f"[{panel}] {int((cL <= 0).sum())} control cells with L=0 kept in the split, dropped from PCs")
        tg_all = P.targets()
        counts = {t: int((P.lab == t).sum()) for t in tg_all}
        if mode == "real":
            ptg, dpred_all, sha = load_pred(panel, pred, P.gcap, full_genes=fg)
            meta["pred_sha256"] = sha
            tg = [t for t in ptg if counts.get(t, 0) >= 2]
            meta["targets_dropped"] = sorted(set(ptg) - set(tg)) + sorted(set(tg_all) - set(ptg))
            sel = [list(ptg).index(t) for t in tg]
            dpred = dpred_all[sel]
            dpred_seed = {s: dpred for s in seeds}
            if ctrl_ab:
                for s in seeds:
                    t2, d2, sha2 = load_pred(panel, f"{pred}_ctrlAB_seed{s}", P.gcap, full_genes=fg)
                    assert list(t2) == list(ptg), f"{panel} seed {s}: ctrlAB target order differs"
                    dpred_seed[s] = d2[sel]; meta[f"pred_ctrlAB_seed{s}_sha256"] = sha2
            halves = {s: [] for s in seeds}
            for i, t in enumerate(tg):
                X, L = P.read(np.flatnonzero(P.lab == t))
                for s in seeds:
                    halves[s].append(target_half_sums(X, L, s, t))
            log(f"[{panel}] {len(tg)} targets read {time.time()-t0:.0f}s")
        stores, splits = [], {}
        for s in seeds:
            split, src = ctrl_split(panel, s, len(P.ctrl_rows), which, write_split, P.ctrl_rows)
            splits[s] = src
            pc = split[0][cL[split[0]] > 0]
            subs = subspaces_for_seed(cell_log(cX[pc], cL[pc]), U_other, s, KMAX, n_perm)
            dA = dB = None
            if mode == "real":
                rA, rB = residuals(halves[s], cX, cL, split, dpred_seed[s])
                if two_level:
                    rng = np.random.default_rng([BOOT_SEED, s, 17])
                    dA = ctrl_boot_deltas(cX, cL, split[1], B, rng); dB = ctrl_boot_deltas(cX, cL, split[2], B, rng)
            else:
                rA, rB, _ = negctrl_residuals(cX, cL, split, list(counts.values()), s)
            stores.append(boot_store(rA, rB, subs, dA, dB, KMAX))
            del dA, dB
            log(f"[{panel}] seed {s} done {time.time()-t0:.0f}s ({src})")
    finally:
        P.close()
    meta["splits"] = splits
    boots = ["target"] + (["two_level"] if (two_level and mode == "real") else [])
    rows = []
    for bt in boots:
        rows += [dict(r, stratum="all", boot=bt) for r in summarize2(stores, ks, B, two_level=(bt == "two_level"))]
    if mode == "real":
        amp = np.linalg.norm(dpred, axis=1)                    # main (non-ctrlAB) file, analysis gene space
        qs = np.quantile(amp, [0.25, 0.5, 0.75])
        lab = np.searchsorted(qs, amp, side="right")             # 0..3
        meta["amp_quartile_edges"] = qs.tolist()
        meta["amp_quartile_n"] = np.bincount(lab, minlength=4).tolist()
        for qi in range(4):
            for bt in boots:
                rows += [dict(r, stratum=f"Q{qi+1}", boot=bt)
                         for r in summarize2(stores, ks, B, idx=np.flatnonzero(lab == qi), two_level=(bt == "two_level"))]
    import pandas as pd
    df = pd.DataFrame(rows); df.insert(0, "panel", panel); df.insert(1, "mode", mode)
    df.insert(2, "pred", tag); df.insert(3, "space", space)
    stem = out_dir / f"{mode}_{panel}_{tag}{'_full' if space == 'full' else ''}"
    df.to_parquet(f"{stem}.parquet", index=False)
    np.savez_compressed(f"{stem}_parts.npz", names=np.array(stores[0]["names"]),
                        **{f"s{s}_{k}": v for s, st in zip(seeds, stores) for k, v in st.items() if k != "names"},
                        targets=np.array(tg if mode == "real" else [f"pseudo{i}" for i in range(N_PSEUDO)]))
    meta["runtime_s"] = time.time() - t0
    pathlib.Path(f"{stem}.json").write_text(json.dumps(meta, indent=1, default=str))
    log(f"[{panel}] wrote {stem}.parquet in {meta['runtime_s']:.0f}s")
    return df


def collect():
    import pandas as pd
    for mode, name in (("real", "rq3_energy"), ("negctrl", "rq3_negctrl")):
        fs = sorted(OUT.glob(f"{mode}_*.parquet"))
        if fs:
            pd.concat([pd.read_parquet(f) for f in fs]).to_parquet(common.RESULTS / f"{name}.parquet", index=False)
            print(f"{name}.parquet <- {len(fs)} files")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel"); ap.add_argument("--mode", default="real", choices=["real", "negctrl"])
    ap.add_argument("--pred", default="prem3"); ap.add_argument("--which", default="truth")
    ap.add_argument("--collect", action="store_true"); ap.add_argument("--ctrl-ab", action="store_true")
    ap.add_argument("--space", default="gcap", choices=["gcap", "full"])
    ap.add_argument("--panels", help="comma list: run several panels sequentially (Task 16 batch)")
    ap.add_argument("--modes", default=None, help="comma list of modes for --panels (default: --mode)")
    a = ap.parse_args()
    if a.collect:
        collect()
    elif a.panels:
        for pn in a.panels.split(","):
            for md in (a.modes or a.mode).split(","):
                run(pn, md, a.pred, which=a.which, ctrl_ab=a.ctrl_ab, space=a.space)
    else:
        run(a.panel, a.mode, a.pred, which=a.which, ctrl_ab=a.ctrl_ab, space=a.space)
