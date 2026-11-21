---
title: "Huu Loc & Van 2026 — Closed-form, control-aware cell quota for perturbation directions"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: huuloc-2026-cell-quota
depth: standard
tags: [sample-size, angular-error, control-pool, tahoe, crispri, noise-floor, sequencing-depth, h1, h4]
relations:
  - target: directional-cell-quota
    type: supports
    confidence: extracted
    because: "«We derive a closed-form, anisotropic cell quota for the treated cells required to hold the angular error of an estimated perturbation direction at a chosen tolerance» (p. 1)."
  - target: tahoe-100m
    type: supports
    confidence: extracted
    because: "Characterises Tahoe-100M's noise regime: «55% of 56,827 conditions are control-pool-limited and only 11% are over-sampled at a 5.7-degree tolerance» (p. 1; Table 1, p. 8)."
  - target: rq2_params-assumptions
    type: supports
    confidence: extracted
    because: "Independent support for the project's T8 finding that plate × line DMSO noise does not average out over drugs: «One pooled DMSO_TF centroid of about 3,113 cells per (plate × line) serves a median of 94 conditions, and its own sampling noise sets a median floor …» (p. 7). The source says nothing about b_c or C*."
  - target: zero-shot-transfer-ceiling
    type: supports
    confidence: inferred
    because: "Shows per-condition measurement noise (treated and control sampling) bounding what any response can be resolved to (pp. 7, 10–11); the source states no prediction ceiling."
  - target: control-subspace-projection
    type: supports
    confidence: ambiguous
    because: "«drug displacements do not concentrate on the high-variance axes» of the within-condition covariance (p. 14), in line with Molina & Zhang's residual-outside-control-PCs reading; but Σ is pooled within-condition and diagonal on chemical data, not a control-only PCA, and no transfer residual is analysed."
---

# How many cells resolve a perturbation direction? A closed-form, control-aware cell quota for single-cell perturbation screens

## Bibliographic Details

Huu Loc, Tran Thai, and Nguyen Thi Van. 2026. "How many cells resolve a perturbation direction? A closed-form,
control-aware cell quota for single-cell perturbation screens." bioRxiv preprint, posted 5 October 2026.
doi:10.64898/2026.10.02.754788. CC BY 4.0. Code: github.com/huulocmedvnu/sample_sufficiency_calculator (p. 15).
Read from `<library>/pdf/huuloc-2026-cell-quota.pdf` (18 pp.; printed page numbers = PDF page numbers).
Citekey: @huuloc-2026-cell-quota.

**Name order.** The bib entry follows Crossref (family "Huu Loc"). The PDF's author biography calls the first author
"Dr. Tran" (p. 15), so the family name is most likely **Tran**. The bibkey is left unchanged here (it is a join key);
see the report to the caller.

**Not in the PDF.** Supplemental Notes S1–S3, Supplemental Methods and Tables S1–S4 are "available online" (p. 15) and
were not read. The acknowledgments state that Claude (via Claude Code) drafted and revised the text and code under the
authors' direction (p. 15).

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract | 1 |
| ✓ F1 | Introduction | 1–2 |
| ✓ F1 | Methods: The angular-error model and the cell quota | 2–3 |
| ✓ F1 | Methods: Datasets and preprocessing | 3–4 |
| ✓ F1 | Methods: Units, controls and regime classification | 4 |
| ✓ F1 | Methods: Downsample-and-measure protocol | 4 |
| — | Methods: Verification and software | 4–5 |
| — | Results: Overview of the framework (Fig. 1, p. 6) | 5–6 |
| ✓ F1 | Results: Only perpendicular noise rotates a direction (Fig. 2) | 6–7 |
| ✓ F1 | Results: The sufficiency spectrum of Tahoe-100M is set by its control pool (Fig. 3, Table 1) | 7–9 |
| ✓ F1 | Results: Tolerance and control-pool size both act on the floor (Table 2) | 9 |
| ✓ F1 | Results: Consequences for screen design — compression, replication, budget and reliability (Table 3) | 9–10 |
| ✓ F1 | Results: Genetic perturbation: detection-limited genome-wide, depth-limited for essential genes (Table 4, Fig. 4) | 10–11 |
| — | Results: The angular-error law holds on every atlas where the signal is strong enough (Fig. 5, Table 5) | 12–13 |
| ✓ F1 | Discussion (incl. seven limitations and three extensions) | 13–14 |
| — | Conclusion | 14–15 |
| — | Supplemental Material; Data and Code Availability; Acknowledgments; Author Biographies | 15–16 |
| — | References (39) | 16–18 |

## Kernthesen

### 1. Direction-based screens need a sample-size rule for angular precision, and none exists (pp. 1–2)

Atlases are used to place an unannotated perturbation next to known ones by the **direction** of its centroid
displacement from matched controls in a shared embedding (connectivity-map reasoning, pp. 1–2). «Two mechanistically
distinct perturbations are separable only when each direction is resolved more tightly than the angle between them,
so angular precision is the binding constraint and raw cell count is only a proxy for it» (p. 2). Conventional power
calculations assume isotropic noise. scPower and bulk replicate planning address detection, «and neither accounts for
the control arm» (p. 2).

### 2. ⭐ Only noise perpendicular to the effect rotates its direction (pp. 2–3; Fig. 2, p. 7)

With v̂ = μ̂_t − μ̂_c ∼ N(v, S), S = Σ_t/n_t + Σ_c/n_c, the tangent projector P = I − uuᵀ gives the exact Jacobian
P/m of the normalisation map. A first-order delta expansion in the high-SNR regime ρ² = m²/(uᵀSu) ≫ 1 yields
E[θ²] = tr(PSP)/m² + O(ρ⁻⁴) = (tr S − uᵀSu)/m² + O(ρ⁻⁴) (pp. 2–3). Noise along the effect changes only its length.
For isotropic noise tr(PΣP) = σ²(d − 1). The first-order rule has relative error O(θ*²), below 1 % at θ* = 0.1 rad
(p. 3, Supp. Note S1).

### 3. ⭐⭐ A closed-form two-arm cell quota (p. 3)

Setting E[θ²] = θ*² with a common Σ gives n_t* = 1 / (m²θ*²/tr(PΣP) − 1/n_c) (p. 3). It reduces to
2·tr(PΣP)/(m²θ*²) for equal arms and to tr(PΣP)/(m²θ*²) as n_c → ∞. A Laurent–Massart bound gives a confidence quota
(Supp. Note S2). For subsampling from N₀ acquired cells the variance carries the finite-population factor
(1/n − 1/N₀) (p. 3). Since n* ∝ 1/m², «halving the magnitude quadruples the requirement» (p. 7).

### 4. ⭐⭐ The control pool sets a magnitude floor below which no treated depth suffices (p. 3; p. 9)

The quota's denominator vanishes at m_min = √(tr(PΣP)/(n_c θ*²)). «below it the control's own sampling noise rotates
the direction by more than θ⋆, n⋆t = ∞, and the only remedies are a larger control pool or a looser tolerance»
(p. 3). Because m_min ∝ 1/(θ*√n_c), a condition below the floor can be freed only by enlarging the shared control or
coarsening the tolerance (p. 9). «A shared vehicle or non-targeting control must scale its pool with the number of
conditions it references, or the screen buys treated precision it cannot realize» (p. 9).

### 5. The deterministic core is machine-checked; the probabilistic part is simulation-checked (p. 4)

SymPy confirms the Jacobian P/m (residual 2.8 × 10⁻¹⁷). Monte-Carlo with K = 50,000 reproduces E[θ²] to 0.13 % and the
second-order correction to < 0.5 %. The tail bound's empirical exceedance is 0.34 % at δ = 10 %. The Jacobian and the
trace identity tr(PΣP) = tr Σ − uᵀΣu are proven in Lean 4/Mathlib; «The probabilistic results are not formalized»
(p. 4).

### 6. The quota depends on an effective noise dimension well below the ambient dimension (p. 7)

d_eff = (tr PΣP)² / tr((PΣP)²) is «well below d − 1 (12.2 in the Monte-Carlo check)» (p. 7). With Tahoe's calibration
(d = 50, σ² = 0.9567, θ* = 0.1 rad) the large-control limit is n* = (d − 1)σ²/(m²θ*²) = 4,688/m². A strong cytotoxic
with m = 15 needs about 20 cells, a median condition about 3,900 (p. 7).

### 7. ⭐⭐ Most Tahoe-100M conditions are control-pool-limited against their plate vehicle (pp. 1, 7–8; Fig. 3, Table 1)

Over all 56,827 drug × dose × line conditions at θ* = 0.1 rad: 55.4 % pool-limited, 5.6 % not detectable, 25.2 %
treated-depth-limited, 2.4 % "ghost" (n* > 50,000), and only 11.4 % over-sampled (p. 7; Table 1, p. 8). Among finite
quotas the median is 3,622 cells against a median acquired depth of 1,296 (p. 7). The plate-shared DMSO_TF vehicle of
about 3,113 cells per plate × line serves a median of 94 conditions, and its median floor m_min = 1.22 exceeds the
median bias-corrected magnitude 1.096 (p. 7). The conclusion holds under a coarser marginal variance (Supp. Note S3).

### 8. ⭐ Effect magnitude is set jointly by compound and cell line, so quotas cannot come from mechanism labels (pp. 8–9)

240 of 379 drugs have median magnitude below the floor; only 9 are over-sampled in most conditions (Supp. Table S2,
p. 8). The strongest displacers are harringtonine, homoharringtonine, bortezomib and panobinostat (median m ≈ 3.4–5.6,
n* ≈ 125–340); 4EGI-1 and (S)-crizotinib (m ≤ 1.3) are under-sampled or pool-limited almost everywhere (p. 8).
Homoharringtonine at 5 µM reaches m = 14.9 in NCI-H460 but only m = 3.0 in NCI-H661, whose plate vehicle holds about
210 cells (m_min = 4.7), so it is pool-limited there (pp. 8–9). The over-sampled fraction rises from 6.9 % at the
lowest dose to about 13.5 % at the higher doses (pp. 8–9).

### 9. ⭐ The control reference is part of the estimand; a per-line-mean reference removes the floor but measures something else (p. 9)

Referencing each condition to its line's mean over about 94 treatments removes the floor (m_min ≈ 0.05) and gives
21.7 % over-sampled with no pool limitation. «That reference measures “difference from the average drug” and not
“difference from untreated”, and it claims a control no real experiment has, so we report it only as a sensitivity
analysis» (p. 9). Downsampling the 6,467 conditions predicted over-sampled against the vehicle confirms 99.6 % within
1.05 θ*; the 25 exceptions are strong responders, 23 in NCI-H460, «which points to a within-condition covariance for
strong responders that departs from the atlas-wide Σ» (p. 9).

### 10. Loosening the tolerance frees pool-limited conditions without extra cells (p. 9; Table 2)

From θ* = 0.10 to 0.30 rad the median floor falls 1.22 → 0.41, pool-limited 55.4 % → 7.2 %, over-sampled 11.4 % →
68.1 %; not-detectable stays 5.6 % because detectability does not depend on the angular target (Table 2, p. 9). «On
this atlas, enlarging the shared DMSO pool would resolve more conditions than enlarging any treated arm» (p. 9).

### 11. Quota-guided allocation saves cells only where the control is large (pp. 9–10; Table 3)

Compression to the quota leaves an over-sampled EmeraldBay centroid stable (distance ≤ 0.6 %, neighbour cosines
≤ 0.005) but distorts an under-sampled one (distance +19 %, cosines −0.10 to −0.12) (pp. 9–10, Supp. Table S4). Under
the per-line-mean reference, adaptive allocation cuts a 600-condition screen from 6.6 M to 2.4 M treated cells (2.8×;
Table 3, p. 10). Against the real vehicle, adaptive and flat allocation give the same drug-similarity graph (k-NN
Jaccard within 0.006), and mechanism-of-action recovery does not beat the majority class even at full depth
(accuracy 0.062 vs 0.142) (p. 10).

### 12. ⭐ At a stringent tolerance most similarity edges rest on an under-resolved direction (p. 10)

A direction is resolved only if it clears the vehicle floor and N₀ ≥ n*: 11.4 % of conditions at 0.1 rad, 44.8 % at
0.2 rad, 68.1 % at 0.3 rad. «only 7.5% of k = 10 neighbour edges at 0.1 rad connect two resolved directions (38.6% and
64.0% at 0.2 and 0.3 rad)» (p. 10).

### 13. Tahoe has almost no replication (p. 10)

«Replicates are largely absent in Tahoe-100M (86.3% of drug–dose combinations occupy one plate)», and among finite
quotas the shortfall between quota and acquired depth is about 2.8-fold (p. 10).

### 14. ⭐⭐ Per-cell noise variance is portable across platforms and modalities (pp. 10, 13; Table 1)

In one fixed embedding the within-condition variance is σ² = 0.957 (Tahoe), 0.917 (EmeraldBay), 0.913 (Orion HCT116)
and 1.033 (Orion HEK293T) — «in the same band … despite a different perturbation type, platform and laboratory. What
transfers is the measurement noise» (p. 10). Large-pool screens share an operating constant of 4,500–9,100 per m²
(p. 10). σ² rises to 1.50 / 1.87 on the essential-gene TRADE screens «because strong and diverse responses spread
variance across many components» (p. 11). σ² depends on the embedding's normalisation scale «and must be
re-estimated for each platform and pipeline» (p. 13).

### 15. ⭐⭐ Genome-wide CRISPRi is detection-limited; essential-gene CRISPRi is depth-limited (pp. 10–11; Table 4, Fig. 4)

«Only 14.7% (HCT116) and 35.3% (HEK293T) of knockdowns clear the detection floor, because a genome-wide library is
dominated by weak non-essential knockdowns (median bias-corrected m 0.14 and 0.35) acquired at only 155–204 cells»
(p. 10). Among detectable knockdowns 0.0 % and 0.1 % are over-sampled, at a median deficit of 35× and 17× (p. 10;
Table 4, p. 11). The strongest are core essential machinery (ribosomal proteins, ribosome biogenesis, Pol II
subunits, splicing factors, NUP93, EIF2S3; pp. 10–11). TRADE essential-gene knockdowns are strong (median m 1.84 /
2.25) and mostly detectable (70.6 % / 60.9 %) but sequenced at only 48–85 cells per gene, so 1.2 % / 4.4 % of
detectable knockdowns are over-sampled (p. 11). «a screen's regime is set the moment its library and its depth are
chosen» (p. 13).

### 16. The angular-error law holds wherever the signal is strong enough (p. 12; Fig. 5, Table 5)

Downsampling every group with ≥ 400 cells, the parameter-free slope tr(PΣP)/m² is recovered at a median
realized/predicted ratio of 0.94 (Tahoe, 1,790 of 56,195 groups with ρ² ≥ 3) and 0.98 (EmeraldBay, 25 of 958), median
R² ≈ 0.999 (p. 12; Table 5, p. 13). Below ρ² = 3 the deficit grows monotonically with 1/ρ² (Spearman 0.57–0.95). Orion
has no group with ρ² ≥ 3 (max 1.3 and 2.0), so its median ratios 0.15 / 0.30 reflect the magnitude distribution «and
not the law» (p. 12). TRADE knockdowns with m > 4 recover the slope (ratio 1.02 / 0.96, R² = 0.996; p. 12).

### 17. Effect sizes are bias-corrected for the sampling floor before any quota is computed (p. 4; p. 13)

m² → max(0, m²_raw − tr(S)) with tr(S) = tr(Σ)(1/n_t + 1/n_c); a condition is detectable when m_raw/√tr(S) > 1.5
(p. 4). «Because raw displacement is inflated by sampling noise at shallow depth … the quota subtracts the
sampling-floor term from the squared magnitude on every dataset before it is applied» (p. 13).

### 18. Covariance is anisotropic, but drug displacements avoid its high-variance axes (pp. 13–14)

The within-condition covariance spans 3.6–5-fold across components on the chemical atlases, «yet tr(PΣP) varies by
only 1.5–3% across conditions because drug displacements do not concentrate on the high-variance axes» (pp. 13–14).

### 19. Directional precision is a different functional from detection power (p. 14)

Detection is governed by the Hotelling non-centrality vᵀΣ⁻¹v, which rewards low-noise directions; direction by
tr(PΣP), which grows with perpendicular noise, «so a screen can be well powered to detect a shift yet under-powered to
orient it» (p. 14). The authors propose the quota as «the sample-size standard for direction-based single-cell
screens» (p. 14).

### 20. ⭐ The noise model omits between-well, plate and replicate variance, so it understates requirements (p. 14)

The seven limitations (p. 14): (1) «Between-well biological variation, plate effects and cell-cycle composition do
not shrink with cell number and are not included, so n⋆ understates the requirement where they dominate»; (2) the
downsampling test resamples the cells that define its reference and «does not test the independence of replicate
wells», and the over-sampled test shares the control centroid, so «a replicate-well test of m_min remains to be
done»; (3) first-order only, counts for weak detectable conditions are orders of magnitude; (4) θ* = 0.1 rad is a
convention; (5) chemical spectra use a diagonal Σ and EmeraldBay vs Tahoe use different references (per-line mean vs
vehicle), so their rows are not directly comparable; (6) only the pseudobulk centroid direction is certified;
(7) only the deterministic algebra is machine-checked.

## Figures and tables

- Fig. 1 (p. 6): study overview with per-dataset σ² and headline result.
- Fig. 2 (p. 7): geometry of the angular error (along-signal vs perpendicular noise).
- Fig. 3 (p. 8): **Tahoe sufficiency spectrum** by regime, and finite quotas on a log axis.
- **Table 1 (p. 8): sufficiency spectrum across six screens** (n_c, σ², large-control bound, regime %, median m and depth).
- Table 2 (p. 9): pool limit vs tolerance on Tahoe.
- Table 3 (p. 10): budget comparison for a 600-condition screen.
- **Table 4 (p. 11): genetic sufficiency spectrum** (detectable %, deficit, pool-limited, held-out slope).
- Fig. 4 (p. 11): magnitude distributions and regime splits across screens.
- Fig. 5 (p. 12) / Table 5 (p. 13): held-out downsample-and-measure test.

## What this source does not address

- **No prediction model and no ceiling for prediction.** It is a sampling-design paper: it scores no predictor, gives
  no correlation or reproducibility ceiling, and uses no split-half estimate. Noise-corrected prediction ceilings are
  the project's [[zero-shot-transfer-ceiling]]; split-half reliability ceilings are in @wang-2026-reliable-perturbations
  (not yet ingested).
- **No context axis.** Cells per condition, not cell lines, are the sample size. Lines appear only as part of the
  condition and as a modulator of magnitude (p. 9). Context counts are the subject of @dibaeinia-2026-context-scale and
  of the project's H2 (not addressed here).
- **No control-PC analysis or correction.** "Control-aware" means the control-arm sampling term 1/n_c and the floor
  m_min (p. 3), not projection onto control principal components (see [[control-subspace-projection]]).
- **No replicate-well or between-plate variance** (p. 14). Tahoe has almost no replicates to estimate it (p. 10).
- **Gene space and metrics.** Everything is in a 50-PC embedding of 2,000 HVGs with angular error, k-NN Jaccard and
  MoA accuracy (pp. 3, 10). No gene-level DE, PDS or VCC-style composite.
- **No STATE/linear-baseline style models, no hESC or KO screens.** Genetic data are CRISPRi only (Orion, TRADE).

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 4, 7, 8, 9, 13, 14, 15, 17, 18, 20.** The hypotheses are those of the plan at
`input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3 at the end).

### What this focus needs beyond the Kernthesen

1. **Datasets, exactly** (pp. 3–4). Tahoe-100M: 100,648,790 cells, 379 small molecules × 3 doses, **50** cancer lines,
   14 96-well plates; 95,624,334 cells pass QC; 56,827 conditions, median N₀ = 1,296; 2–3 DMSO_TF wells per plate,
   pooled to one vehicle centroid per plate × line (median n_c ≈ 3,113). EmeraldBay: same Vevo Mosaic platform, read
   at 5 days, 1,831,648 cells, 52 lines, 91 treatments + 2 vehicles, 4,912 conditions, median 240 cells, referenced to
   the per-line mean (n_c ≈ 2.3 × 10⁴). X-Atlas/Orion (@huang-2025-x-atlas-orion): genome-wide dual-sgRNA CRISPRi,
   HCT116 and HEK293T, 18,903 genes, about 8 M cells, pooled NTC 1.7–2.2 × 10⁵. TRADE essential-gene screens
   (@nadig-2025-transcriptome-wide, GSE264667): Jurkat and HepG2, 2,393 genes, about 2.5 × 10⁵ cells per line, median
   48–85 cells per gene, NTC 11,514 / 4,380. Genetic knockdowns need ≥ 25 cells (p. 4).
2. **Metrics, exactly.** Angular error θ = arccos(ûᵀu) at θ* = 0.1 rad (5.73°); five regimes (over / under / ghost /
   pool-limited / not detectable; p. 4); held-out slope ratio with R² and Spearman deficit vs 1/ρ² (p. 4); k = 10
   neighbour-edge resolution and k-NN Jaccard; MoA label accuracy vs majority class (p. 10). Embedding: normalise to
   10⁴, log1p, 2,000 HVGs, PCA to 50 (p. 3). Σ is the pooled within-condition covariance from single cells; diagonal
   for chemical, full 50 × 50 for genetic (p. 3).
3. **What "control-aware" means here — and what it is not.** The control enters as sampling variance Σ_c/n_c and the
   floor m_min (p. 3). The projector P = I − uuᵀ removes the **effect** direction from the noise; it does not deflate
   control PCs. There is no analogue of M3 (`--sum-defl 8,0.75,0.5`, plan line 161).
4. **The project's Tahoe target sits in the pool-limited regime (a reader's link).** The H2a/H2b response is "matched
   to the plate × line DMSO pool" ([[rq2_params-assumptions]]), the reference whose noise Kernthese 7 shows to dominate.
   The project's T8 finding that b's empirical split-half noise is ~37× the independent-noise model (log, 2026-11-05)
   has the same mechanism the source quantifies: one vehicle centroid serves about 94 conditions (p. 7). Tiny pools
   occur in both (NCI-H661 vehicle about 210 cells, p. 9; the project's plate14 × CVCL_1571 pool with 1-cell wells).
5. **Noise not modelled (relevant to H2b's σ²_ε).** Between-well, plate and cell-cycle variance «do not shrink with
   cell number and are not included» (p. 14), and 86.3 % of Tahoe drug–dose combinations sit on one plate (p. 10). A
   cell-split-half σ²_ε has the same blind spot. This is one more reason why an independent-noise linear-Gaussian
   formula need not calibrate on Tahoe (reader's inference; the source says nothing about context models).

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **partial (noise-floor ingredient only)** | A per-condition sampling-noise floor (bias correction m² − tr(S), detection at SNR > 1.5) and the share of conditions whose response is not resolvable at all. Genome-wide Orion HCT116: 85.3 % of knockdowns not detectable at median 155 cells (Table 1). The harness panel OR_HCT116X uses Orion, so a low Var_reproducible and a low, noisy r_max are to be expected there (reader inference). No prediction, no reproducible-variance ceiling, no split-half, no correlation bound. | p. 4, p. 10, Table 1 p. 8, Table 4 p. 11 |
| **H2a** (Tahoe: context correction negative at C = 12; C*_emp ≥ 24) | **none** | No learned model and no context-count axis. Indirect: most Tahoe conditions are noise-dominated against the plate vehicle (55.4 % pool-limited), consistent with the very small residual-MSE gains (~10⁻⁴) the project measures. The source does not bear on the 37.4 crossing itself. | p. 7, Table 1 p. 8 |
| **H2b** (C*_theory within ×2 of C*_emp) | **none (adjacent evidence on the noise model)** | No break-even theory. Two adjacent points: (i) control noise shared by ~94 conditions per plate × line (p. 7) and unmodelled between-well / plate variance (p. 14) contradict the independent-noise assumption behind C* = d + 1 + d·σ²_res/(R²·Var(i)); this supports dropping the quantitative Proposition 2 and is a citable mechanism for the T8 S1emp finding. (ii) Its d_eff = (tr PΣP)²/tr((PΣP)²) has the participation-ratio form of the project's d, but on noise, not on descriptors (p. 7). | p. 7, p. 10, p. 14 |
| **H2c** (genetic C*_min > C_public) | **none** | No count of public lines. Adjacent: genome-wide CRISPRi knockdowns are mostly below the detection floor (14.7–35.3 % detectable, median m 0.14–0.35), so σ²_ε/Var(i) is large on genome-wide data, which pushes C*_min up. Given H2b's calibration failure this can only be a qualitative remark. | p. 10, Table 4 p. 11 |
| **H3a** (transfer residual on destination top-8 control PCs) | **none** | "Control-aware" ≠ control-PC projection (focus item 3). The only related statement: drug displacements do not concentrate on the high-variance within-condition axes (pp. 13–14), which concerns observed responses, a pooled Σ, chemical data, and no transfer residual. | p. 3, pp. 13–14 |
| **H3b** (M3 beats five corrections on PDS) | **none (framing only)** | No correction is compared. The per-line-mean reference (p. 9) is a common-response centring like the Systema-style alternative; the authors call it a different estimand («difference from the average drug») and report it only as a sensitivity. | p. 9 |
| **H4** (adding H1GW beats the largest modelling gain) | **partial (depth / control-pool framing only)** | Genome-wide CRISPRi directions are under-resolved per gene (0.0–0.1 % over-sampled, deficit 17–35× at 155–204 cells); essential-gene screens at 48–85 cells per gene are depth-limited; a screen's regime is fixed by its library and depth; enlarging a shared control can beat enlarging treated arms. So a same-line genome-wide screen (H1GW, @pan-2026-global-cell-state) is likely to add value through strong essential knockdowns and shared structure, not through resolving weak ones (reader inference). No prediction, no adding a screen to a predictor, no composite score. | pp. 9–11, p. 13, Table 4 p. 11 |

**Positioning.** The source is the **cells axis** complement to the project's **contexts axis**: it shows that the
per-condition response itself is often unresolvable (Tahoe: control pool; genome-wide CRISPRi: effect size and depth).
It does not change C*_emp ≈ 37 or the H2a verdict. It strengthens the case for dropping the quantitative
Proposition 2: the Tahoe noise is block-shared through the plate vehicle and has unmodelled plate / well components,
so an independent-noise σ²_ε is misspecified there.

### Direct quotes (for the draft)

> "One pooled DMSO_TF centroid of about 3,113 cells per (plate × line) serves a median of 94 conditions, and its own
> sampling noise sets a median floor m_min = 1.22 that lies above the median bias-corrected effect magnitude of 1.096
> (raw 1.14). For most conditions the control alone rotates the direction by more than θ⋆, so no treated depth resolves
> them." (p. 7; subscripts normalised from the PDF text layer)

> "Because magnitude is set jointly by compound and context, quotas must be calibrated per condition from the measured
> magnitude and cannot be inferred from mechanism labels." (p. 9)

> "Only 14.7% (HCT116) and 35.3% (HEK293T) of knockdowns clear the detection floor, because a genome-wide library is
> dominated by weak non-essential knockdowns (median bias-corrected m 0.14 and 0.35) acquired at only 155–204 cells."
> (p. 10)

> "Between-well biological variation, plate effects and cell-cycle composition do not shrink with cell number and are
> not included, so n⋆ understates the requirement where they dominate" (p. 14)

> "That reference measures “difference from the average drug” and not “difference from untreated”, and it claims a
> control no real experiment has, so we report it only as a sensitivity analysis." (p. 9)

### Examples & illustrations (for later drafting)

- Homoharringtonine 5 µM: m = 14.9, n* = 21 in NCI-H460 (N₀ = 6,060, about 290-fold surplus) but m = 3.0 and
  pool-limited in NCI-H661, whose plate vehicle has about 210 cells (m_min = 4.7) (pp. 8–10). One drug, two lines, two
  regimes: a concrete line × perturbation interaction in magnitude, and a concrete control-pool failure.
- Orion HCT116 vs TRADE HepG2: weak knockdowns at moderate depth (m 0.14, 155 cells) vs strong knockdowns at very
  shallow depth (m 2.25, 48 cells); both fail the quota for opposite reasons (p. 13; Table 1, p. 8).
- MoA recovery on Tahoe's similarity graph: accuracy 0.062 against a majority-class baseline of 0.142 even at full depth
  (p. 10).

### Boundary of this focus

This focus does not work up the verification suite, the held-out law test (Kernthesen 5, 16) or the budget
arithmetic beyond its regime logic. A reader looking for a context-count curve, a prediction ceiling or a
control-PC correction will not find one in this source.

## Mentioned entities

- Concepts: [[directional-cell-quota]], [[zero-shot-transfer-ceiling]], [[control-subspace-projection]]
- Datasets: [[tahoe-100m]], [[replogle-nadig-crispri-collection]] (the TRADE HepG2 / Jurkat screens)
- Project notes: [[rq2_params-assumptions]]

## Connections

- Introduces / supports: [[directional-cell-quota]] (p. 1, p. 3).
- Supports (characterises the noise regime of): [[tahoe-100m]] (pp. 7–10).
- Supports (independent mechanism for the T8 shared-DMSO-noise finding): [[rq2_params-assumptions]] (p. 7).
- Supports (inferred; measurement noise bounds resolvable signal): [[zero-shot-transfer-ceiling]].
- Supports (ambiguous; displacements avoid high-variance axes, p. 14): [[control-subspace-projection]].
- Shares data with [[molina-2026-response-decomposition]] (Tahoe, Nadig HepG2/Jurkat) but neither cites the other; no
  stance.
- Cites (not yet ingested): @zhang-2026-tahoe, @huang-2025-x-atlas-orion, @nadig-2025-transcriptome-wide,
  @replogle-2022-genome-scale.
