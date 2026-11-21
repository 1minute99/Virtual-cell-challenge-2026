"""Task 20 of plan zero-shot-limits-context-complexity (input/ideas/zero-shot-limits-context-complexity-plan.md).
Appendix: correlation upper bound for TARGET-IDENTITY-MEDIATED drug -> CRISPRi/shRNA transfer, from the frozen
exploratory concordance analysis vcc2026_data/reports/drug_gene_concordance_2026-09-16/ (RESULTS.md, pair_metrics.csv,
repeat_metrics.csv, summary.csv). Note: knowledge/synthesis/appx-drug-gene-bound.md (derivation and caveats).

NOT USED, NOT CITED: the gene side of vcc2026_data/reports/drug_context_pretrain_2026-09-20 (gene context vectors were
scrambled; AGENTS.md "BUG FOUND 2026-09-21", FINDINGS caveat). Nothing from that directory is read here.

DERIVATION (full version in the note). Per (drug d hitting target t, line l), genes as samples, rank-transformed:
  D = s_D + e_D (Tahoe drug response), G = s_G + e_G (genetic knock-down response), e_D, e_G independent of each
  other and of the signals. Write s_G = beta*s_D + u with <u, s_D> = 0. A target-identity-mediated transfer
  predictor is g_hat = a*s_D + w, where w carries no target-specific information about s_G (<w, s_G> = <w, s_D> = 0
  after removing the generic response shared with wrong targets). Cauchy-Schwarz:
      corr(g_hat, s_G) = a*<s_D, s_G> / (||g_hat|| ||s_G||) <= |<s_D, s_G>| / (||s_D|| ||s_G||) = |rho*|,
  equality iff w = 0. rho* (noise-free, target-specific concordance) is estimated by Spearman's disattenuation:
      rho* = Delta / sqrt(R_D * R_G),  Delta = rho(D_t, G_t) - mean_{t' not a target of d} rho(D_t, G_t')
  (matched minus wrong-target concordance = the target-identity part). The bound against the observed (noisy)
  genetic truth is rho* * sqrt(R_G); the share of reproducible genetic-response variance transferable this way is rho*^2.

ASSUMPTIONS (stated; checks in results/appx_drug_gene.json["checks"]):
  D1  Classical attenuation holds for the Spearman correlations used (rank-space approximation of the Pearson result).
  D2  Reliabilities: Jiang genetic response = mean of Rep1 and Rep2 (ANALYSIS_PLAN.md), so R_G = Spearman-Brown of the
      Rep1-vs-Rep2 rank correlation r_G (genetic_repeat_rho): R_G = 2 r_G / (1 + r_G). Tahoe = primary treatment only, so
      R_D = r_D (primary-vs-repeat rank correlation, drug_repeat_rho). r_D is measured on the drugs with an archived
      repeat (7 drugs on the Jiang ~16.6k-gene axis, 26 on the LINCS 922-gene axis) and assumed representative.
  D3  LINCS shRNA has no biological replicate in these files -> R_G unidentified; reported as a curve over R_G and with
      the drug-side-only correction (R_G = 1), which is a FLOOR of the bound, not the bound.
  D4  Shared non-targeting controls between replicates would inflate r_G and r_D, which makes rho* too SMALL (the
      bound anti-conservative). Not identifiable from these files; stated, not corrected.
  D5  Same aggregation as the source analysis: pair -> mean over targets within drug x line -> mean over lines within
      drug -> equal weight per drug. Raw effects are primary (reliabilities were measured on raw effects only); the
      generic-removed transform reuses the raw reliabilities (sensitivity only).
  D6  Bootstrap: B = 2,000 percentile; drugs resampled within each comparison (Delta and r_G jointly, same draw),
      repeat drugs resampled independently for r_D (shared across the four Jiang stimuli within a draw). Draws with
      R_D*R_G <= 0 are counted and dropped. Drug draws ignore shared-target dependence (as in the source analysis).
Seed 20261105."""
import json, pathlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED, B = 20261105, 2000
SRC = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/reports/drug_gene_concordance_2026-09-16")
OUT = pathlib.Path("/home/wonmin-kim/AIVC/output/data-analysis/results")
pm = pd.read_csv(SRC / "pair_metrics.csv"); rm = pd.read_csv(SRC / "repeat_metrics.csv"); sm = pd.read_csv(SRC / "summary.csv")
JI = ["jiang_IFNB", "jiang_INS", "jiang_TGFB", "jiang_TNFA"]; LI = "lincs_full_local_CGS"
checks = {"inputs_not_read": "reports/drug_context_pretrain_2026-09-20 (invalid gene side)"}

def by_drug(df, cols):
    """pair -> mean over targets within drug x line -> mean over lines within drug (D5)."""
    return df.groupby(["drug", "cell_line"])[cols].mean().groupby("drug")[cols].mean()

# ---- D5 reproduce the frozen summary --------------------------------------------------------------------------
rep = {}
for (ds, mode), g in pm.groupby(["dataset", "mode"]):
    v = by_drug(g, ["spearman", "wrong_spearman", "spearman_advantage"]).mean()
    ref = sm[(sm.dataset == ds) & (sm["mode"] == mode) & (sm.subset == "all_annotated_inhibitors")].iloc[0]
    rep[f"{ds}|{mode}"] = float(abs(v.spearman_advantage - ref.spearman_advantage) + abs(v.spearman - ref.spearman))
checks["D5_max_abs_diff_vs_summary_csv"] = max(rep.values())

# ---- reliabilities --------------------------------------------------------------------------------------------
sb = lambda r: 2 * r / (1 + r)
rD_pool = {}
jd = rm[rm.dataset.isin(JI) & rm.drug_repeat_rho.notna()].groupby(["drug", "cell_line"]).drug_repeat_rho.mean()
rD_pool["jiang_axis"] = jd.groupby("drug").mean()
ld = rm[(rm.dataset == LI) & rm.drug_repeat_rho.notna()].groupby(["drug", "cell_line"]).drug_repeat_rho.mean()
rD_pool["lincs_axis"] = ld.groupby("drug").mean()
checks["D2_n_repeat_drugs"] = {k: int(len(v)) for k, v in rD_pool.items()}
checks["D2_rD_spread_across_jiang_stimuli_same_drug_line"] = float(rm[rm.dataset.isin(JI) & rm.drug_repeat_rho.notna()]
                                                                  .groupby(["drug", "cell_line"]).drug_repeat_rho.std().max())

# genetic r_G per drug (same aggregation), joined to Delta per drug, raw mode (reliabilities measured raw)
rG = {ds: by_drug(rm[rm.dataset == ds], ["genetic_repeat_rho"]).genetic_repeat_rho for ds in JI}
checks["D2_rG_all_pairs_present"] = {ds: bool(rm[rm.dataset == ds].genetic_repeat_rho.notna().all()) for ds in JI}
checks["D3_lincs_genetic_repeat_available"] = bool(rm[rm.dataset == LI].genetic_repeat_rho.notna().any())

# D1/D2 model check: averaging two genetic replicates should raise the cross-modal Delta by sqrt(2/(1+r_G)) on average
mc = rm[rm.dataset.isin(JI)].copy()
mc["adv_single"] = (mc.cross_gene_rep1_advantage + mc.cross_gene_rep2_advantage) / 2
mcd = pd.concat([by_drug(mc[mc.dataset == ds], ["cross_primary_advantage", "adv_single", "cross_primary", "cross_gene_rep1", "cross_gene_rep2"]).mean()
                 for ds in JI], axis=1).T.assign(dataset=JI)
r_bar = float(np.mean([rG[ds].mean() for ds in JI]))
checks["D1_attenuation_model_check"] = {
    "expected_ratio_avg_over_single": float(np.sqrt(2 / (1 + r_bar))),
    "observed_ratio_Delta_avg_over_single_pooled": float(mcd.cross_primary_advantage.mean() / mcd.adv_single.mean()),
    "observed_ratio_matched_rho_pooled": float(mcd.cross_primary.mean() / ((mcd.cross_gene_rep1 + mcd.cross_gene_rep2) / 2).mean()),
    "note": "pooled over the 4 Jiang stimuli; small signals, so the ratio is noisy; consistent if near the expectation"}

# ---- point estimates + bootstrap ------------------------------------------------------------------------------
rng = np.random.default_rng(SEED)
def drug_table(ds, mode):
    t = by_drug(pm[(pm.dataset == ds) & (pm["mode"] == mode)], ["spearman", "wrong_spearman", "spearman_advantage"])
    if ds in rG: t = t.join(rG[ds].rename("rG"), how="left")
    return t

res, boots = [], {}
for mode in ["raw", "generic_removed"]:
    rDj, rDl = rD_pool["jiang_axis"].to_numpy(), rD_pool["lincs_axis"].to_numpy()
    idx_j = rng.integers(len(rDj), size=(B, len(rDj))); idx_l = rng.integers(len(rDl), size=(B, len(rDl)))
    RDj_b, RDl_b = rDj[idx_j].mean(1), rDl[idx_l].mean(1)
    jiang_star_b = []
    for ds in JI + [LI]:
        t = drug_table(ds, mode); n = len(t); ii = rng.integers(n, size=(B, n))
        dl, rho = t.spearman_advantage.to_numpy(), t.spearman.to_numpy()
        Dl_b, rho_b = dl[ii].mean(1), rho[ii].mean(1)
        if ds in JI:
            RD, RD_b = rDj.mean(), RDj_b
            rg = t.rG.to_numpy(); RG, RG_b = sb(rg.mean()), sb(rg[ii].mean(1))
        else:
            RD, RD_b = rDl.mean(), RDl_b; RG, RG_b = np.nan, np.full(B, np.nan)
        def ci(x):
            x = x[np.isfinite(x)]; return [float(np.quantile(x, .025)), float(np.quantile(x, .975))] if len(x) else [np.nan, np.nan]
        prod_b = RD_b * (RG_b if ds in JI else 1.0); ok = prod_b > 0
        star_b = np.where(ok, Dl_b / np.sqrt(np.where(ok, prod_b, 1)), np.nan)
        row = {"comparison": ds, "mode": mode, "n_drugs": n, "rho_matched": float(rho.mean()), "Delta_obs": float(dl.mean()),
               "Delta_obs_ci": ci(Dl_b), "r_D": float(RD), "r_D_ci": ci(RD_b), "R_D": float(RD), "n_repeat_drugs_rD": int(len(rDj) if ds in JI else len(rDl))}
        if ds in JI:
            star = dl.mean() / np.sqrt(RD * RG)
            row.update(r_G_single=float(rg.mean()), R_G=float(RG), R_G_ci=ci(RG_b), rho_star=float(star), rho_star_ci=ci(star_b),
                       rho_star_sq=float(star ** 2), bound_vs_observed_truth=float(star * np.sqrt(RG)),
                       bound_vs_observed_truth_ci=ci(star_b * np.sqrt(RG_b)),
                       matched_total_corrected=float(rho.mean() / np.sqrt(RD * RG)), dropped_draws=int((~ok).sum()))
            jiang_star_b.append(star_b)
        else:
            star = dl.mean() / np.sqrt(RD)
            row.update(r_G_single=None, R_G=None, rho_star_drugside_only=float(star), rho_star_drugside_only_ci=ci(star_b),
                       rho_star_vs_R_G={f"{g:.2f}": float(star / np.sqrt(g)) for g in [1.0, 0.5, 0.3, 0.2, 0.1]},
                       matched_total_corrected_drugside_only=float(rho.mean() / np.sqrt(RD)), dropped_draws=int((~ok).sum()))
        res.append(row); boots[(ds, mode)] = star_b
    sb_mean = np.nanmean(np.vstack(jiang_star_b), axis=0)
    pt = float(np.mean([r["rho_star"] for r in res if r["mode"] == mode and r["comparison"] in JI]))
    res.append({"comparison": "jiang_4stim_equal", "mode": mode, "rho_star": pt, "rho_star_ci": [float(np.nanquantile(sb_mean, .025)), float(np.nanquantile(sb_mean, .975))],
                "rho_star_sq": pt ** 2, "note": "equal-weight mean of the four stimulus-specific rho*; one shared r_D draw"})

# max over Jiang comparisons of the upper CI = the most permissive bound these data support
raw_j = [r for r in res if r["mode"] == "raw" and r["comparison"] in JI]
summary = {
  "plan": "input/ideas/zero-shot-limits-context-complexity-plan.md", "task": 20, "seed": SEED, "B": B,
  "jiang_modality": "CRISPRi (dCas9-KRAB-MeCP2; sgRNAs from the Dolcetto CRISPRi library) -- GEO GSE281048 description; Jiang et al. 2025 Nat Cell Biol 27:505-517, doi:10.1038/s41556-025-01622-z (bibkey jiang-2025-pathway-signatures)",
  "lincs_modality": "shRNA knock-down (LINCS L1000 CGS, 96 h)",
  "results": res,
  "headline_raw": {
     "jiang_4stim_equal_rho_star": next(r for r in res if r["mode"] == "raw" and r["comparison"] == "jiang_4stim_equal"),
     "max_point_rho_star_jiang": max(r["rho_star"] for r in raw_j),
     "max_upper_ci_rho_star_jiang": max(r["rho_star_ci"][1] for r in raw_j),
     "max_upper_ci_share_of_reproducible_variance": max(r["rho_star_ci"][1] for r in raw_j) ** 2},
  "checks": checks,
}
(OUT / "appx_drug_gene.json").write_text(json.dumps(summary, indent=1, default=float))
flat = []
for r in res:
    f = {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
    for k, v in r.items():
        if isinstance(v, list): f[k + "_lo"], f[k + "_hi"] = v
    flat.append(f)
pd.DataFrame(flat).to_csv(OUT / "appx_drug_gene.csv", index=False)

# ---- figure: observed Delta vs noise-corrected rho*, raw mode, one axis ----------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 3.4))
lab = {"jiang_IFNB": "Jiang IFNB", "jiang_INS": "Jiang INS", "jiang_TGFB": "Jiang TGFB", "jiang_TNFA": "Jiang TNFA", "jiang_4stim_equal": "Jiang, 4 stimuli", LI: "LINCS shRNA*"}
rows = [r for r in res if r["mode"] == "raw"]
for k, r in enumerate(rows):
    y = len(rows) - k
    if "Delta_obs" in r:
        ax.errorbar(r["Delta_obs"], y + .15, xerr=[[r["Delta_obs"] - r["Delta_obs_ci"][0]], [r["Delta_obs_ci"][1] - r["Delta_obs"]]], fmt="o", ms=5, color="#2a78d6", lw=1.5, capsize=0, label="observed Δ" if k == 0 else None)
    st = r.get("rho_star", r.get("rho_star_drugside_only")); c = r.get("rho_star_ci", r.get("rho_star_drugside_only_ci"))
    ax.errorbar(st, y - .15, xerr=[[st - c[0]], [c[1] - st]], fmt="s", ms=5, color="#eb6834", lw=1.5, capsize=0, label="noise-corrected ρ* (bound)" if k == 0 else None)
ax.axvline(0, color="#52514e", lw=.8); ax.set_yticks([len(rows) - k for k in range(len(rows))]); ax.set_yticklabels([lab[r["comparison"]] for r in rows], fontsize=8)
ax.set_xlabel("target-identity concordance (Spearman, matched − wrong target)"); ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="lower right"); ax.text(0.01, -0.28, "*LINCS: drug-side correction only (no shRNA replicate) = floor of the bound", transform=ax.transAxes, fontsize=7, color="#52514e")
fig.tight_layout(); fig.savefig(OUT / "appx_drug_gene.png", dpi=150)
print(json.dumps({"headline_raw": summary["headline_raw"], "checks": checks}, indent=1, default=float))
for r in res: print(r["comparison"], r["mode"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k in ("Delta_obs", "R_D", "R_G", "rho_star", "rho_star_ci", "rho_star_drugside_only", "rho_star_drugside_only_ci", "bound_vs_observed_truth", "rho_star_vs_R_G")})
