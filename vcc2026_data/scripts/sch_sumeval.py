"""Fast OFFLINE evaluator of the sum-side metrics (vcc2026 definitions) for a pinned-sum prediction:
  mu_pred,p = bulk_lognorm(control_mean * exp(L_p) * n)  (override entries multiplicative), TS = 5e4
  MSE_norm  = sum_p [ ||mu_pred,p - mu_real,p||^2 - r * min(C_pred,p, C_real,p) - C_real,p ]
              / sum_p [ ||mu_real,p - mu_ctrl||^2 - C_real,p - C_ctrl ]        (own gene excluded per row)
              r = min(1, B_pred / sum_p min(C_pred,p, C_real,p)), B_pred = across-perturbation centred SS of mu_pred
  PDS       = mean_p [1 - midrank_p / n], cosine distance of (mu_pred,p - mu_ctrl) to every (mu_real,j - mu_ctrl), panel genes excluded
C_real / C_ctrl from the scorer's real moments cache; C_pred from a reference scored run's pred moments (jk), or C_real.
Deterministic (no emission sampling), seconds per panel. Used to screen sum-side transformations on all six panels."""
import glob, sys
import numpy as np
from cell_eval2.prep import bulk_lognorm_means
R0 = "/home/wonmin-kim/AIVC/vcc2026_data/reports"
TS = 5e4
_REAL = {}
def real_moments(panel):
    if panel not in _REAL:
        c = f"{R0}/transfer_v2_2026-09-23/cache_real_H1" if panel == "H1" else f"{R0}/final_harness_2026-09-24/cache_real_{panel}"
        _REAL[panel] = np.load(glob.glob(f"{c}/pseudobulk_moments_bulk_lognorm-*.moments.npz")[0])
    return _REAL[panel]
def evaluate(panel, targets, genes, ctrl_counts, L, OV, n, cpred=None):
    """targets/genes: panel order; ctrl_counts: CSR controls; L, OV: (n_t x G) sum-side lfc / override (NaN = none)."""
    Rm = real_moments(panel); rp = list(Rm["perts"]); c = "non-targeting"
    ir = np.array([rp.index(t) for t in targets]); mr = Rm["means"][ir]; cr = Rm["jk"][ir]
    m0 = Rm["means"][rp.index(c)]; c0 = float(Rm["jk"][rp.index(c)])
    fm = np.asarray(ctrl_counts.sum(0)).ravel() / ctrl_counts.shape[0]
    e = np.exp(L); j = np.isfinite(OV); e[j] = np.exp(OV[j])
    mp = bulk_lognorm_means(np.round(fm[None, :] * e * n), TS)
    gi = {g: k for k, g in enumerate(genes)}
    own = np.zeros(mp.shape, bool)
    for k, t in enumerate(targets):
        if t in gi:
            own[k, gi[t]] = True
    sq = np.where(own, 0.0, (mp - mr) ** 2).sum(1)
    claim = np.minimum(cpred if cpred is not None else cr, cr)
    mpc = np.where(own, np.nan, mp); B = float(np.nansum((mpc - np.nanmean(mpc, 0)) ** 2))
    r = min(1.0, B / claim.sum())
    num = sq - r * claim - cr
    dist = np.where(own, 0.0, (mr - m0[None, :]) ** 2).sum(1) - cr - c0
    mse = float(num.sum() / dist.sum())
    panel_cols = np.array(sorted({gi[t] for t in targets if t in gi}))
    keep = np.ones(len(genes), bool); keep[panel_cols] = False
    P_ = (mp - m0[None, :])[:, keep]; R_ = (mr - m0[None, :])[:, keep]
    Pn = P_ / (np.linalg.norm(P_, axis=1, keepdims=True) + 1e-12); Rn = R_ / (np.linalg.norm(R_, axis=1, keepdims=True) + 1e-12)
    Dm = 1.0 - Pn @ Rn.T                                   # cosine distance, rows = queries (pred), cols = real
    nT = len(targets); ranks = np.empty(nT)
    for k in range(nT):
        row = Dm[k]; x = row[k]; less = (row < x).sum(); eq = (row == x).sum()
        ranks[k] = less + (eq - 1) / 2.0
    pds = float(np.mean(1.0 - ranks / nT))
    return dict(mse=round(mse, 4), pds=round(pds, 4), r=round(r, 3), B_over_claim=round(B / claim.sum(), 3))
if __name__ == "__main__":
    from cell_eval2.moments import jackknife_correction  # noqa
    sys.argv, arm = ["x"], sys.argv[1] if len(sys.argv) > 1 else "em_mult_share_k3000_du0.5_ag_cpc"
    sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
    import fh_harness as H
    for panel in ["H1", "KOLF"]:
        d, ov, _ = H.arm_delta("blend30_e", panel); r = H.real(panel); ok = r["gmap"] >= 0
        Ld = H.on_panel(d, panel); Od = np.full(Ld.shape, np.nan); Od[:, ok] = ov[:, r["gmap"][ok]]
        print(panel, "offline", evaluate(panel, r["targets"], r["genes"], r["ctrl"], Ld, Od, H.PANELS[panel]["n_cells"]))
