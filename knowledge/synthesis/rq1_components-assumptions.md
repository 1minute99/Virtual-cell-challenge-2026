---
title: "RQ1 Task 12 — noise-corrected variance components and ceilings on the 7 harness panels: assumptions and results"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ1 Task 12: variance components, measurement ceiling and r_max

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 12, operational definitions 4 and 6, H1 (frozen)
and Amendment A3. **H1 is not computed here** (Task 13). No prediction file was read and no scorer was run.

Files:
- Script (library + CLI): `output/data-analysis/rq1_components.py`.
- Extraction of the panel targets that the Task 9 cache lacks: `output/data-analysis/rq1_gene_extract.py`, output in
  `results/rq1_gene_cache/`.
- Tests: `output/data-analysis/tests/test_rq1_components.py` (7 pass).
- Results:
  - `results/rq1_components.parquet` is in long format: one row per panel × other-line set × quantity, with the
    estimate, the target-bootstrap CI, the gene-block jackknife CI and the seed range.
  - `results/rq1_components.json`.
  - `results/rq1_components.{png,svg}`.
  - `results/rq1_components/<panel>_pack.npz` holds the arrays Task 13 needs: the destination truth halves `dA`/`dB`
    and the other-line target profiles `qA`/`qB` for each of the 5 split seeds, plus `mu`, `b`, `a_full`, `delta_full`
    and the target and gene lists.
- Log: `results/logs/rq1_components.log`.

Upstream: [[rq2_params-assumptions]] (the split-half estimator idea and the A1 depth rule),
[[rq3_build_preds-assumptions]] (the panel conventions), and the Task 9 cache `results/rq2_gene_cache/` (read-only).

## Model and estimator

δ_{c,t} = μ + a_t + b_c + i_{c,t} + ε (per gene). The space is the response as log1p CP10K(target) − log1p
CP10K(control) of the same half, on G = 6,009 genes: G_cap ∩ the gene axis of the Task 9 cache.

For destination c0 and its panel targets T (each measured with ≥ 20 cells in ≥ 2 OTHER lines):

| Symbol | Estimate |
|---|---|
| q_t | Unweighted two-way additive fit δ_{k,t} = c + α_t + β_k over the other lines only. The fitted table is completed and averaged over the lines k, which makes it gauge-invariant. |
| μ | mean_t q_t |
| a_t | q_t − μ (taken from the other lines' truth) |
| b_c0 | mean_t δ_{c0,t} − μ |
| i | δ − μ − b − a |

All variances are per element (target × gene) and noise-free. Each one is ⟨X^A, X^B⟩ over cell-level split halves.
Every line has its own **separate control halves**.

- For the Task 9 lines, the halves are Task 9's.
- For the new extraction, the control split is aligned to Task 9. `verify.json` confirms the control sums are
  bit-identical.
- For the destination, the panel truth file is split with seed `[20261112, s, panel, group]`.

The quantities:
- Var_rep = Var(a) + Var(b) + Var(i) + 2C_ai. This identity is exact, with a residual of about 1e-18.
- σ²_ε = ¼·mean(δ^A − δ^B)² is the destination noise at full panel depth (A1).
- The measurement ceiling is √(Var_rep / (Var_rep + σ²_ε)). The split-half r with Spearman–Brown agrees to 3
  decimals.
- The other ceilings:
  - r_max = √[(Var(a) + Var(b))/Var_rep] (the registered formula).
  - r_oracle is the noise-free r of μ + a + b itself.
  - r_proj is the r of the best-rescaled oracle αa + b.
  - r_max_orth = √[(Va + Vb)/(Va + Vb + Vi)] (the Molina & Zhang convention).

These are computed in the primary space, which is μ-centred. Pearson-space versions are in the parquet.

Uncertainty is estimated two ways:
- **Target bootstrap.** The Gram matrices are averaged over 5 seeds, then targets are resampled with B = 2,000. Every
  centring is recomputed for each draw. The other-line fit q is held fixed.
- **Delete-1 jackknife over 20 random gene blocks.** It captures the control noise shared by all targets, which the
  target bootstrap cannot see. The test `test_gene_jackknife_tracks_control_noise_in_b` checks this: the jackknife SE
  is within ×0.6–1.6 of the replicate SD, while the bootstrap SE is less than half of it.

## Assumptions

- **B1.** The destination truth is the harness truth file, which is what is scored. Its library is the total over all
  panel genes, as in Tasks 14 and 15.
- **B2.** The other lines are the cell-level pseudobulks of the Task 9 cache (universe U) plus `rq1_gene_cache` (474
  panel targets outside U, for VCC_H1, CD4, KOLF, K562, HCT116 and HEK293T). The library is the native cell total.
- **B3.** Genes are G_cap ∩ the cache axis. The first panel column is used for each official gene.
- **B4.**
  - FENG is treated as one line, with per-donor responses weighted by 1/(1/n_t + 1/n_c).
  - For the KOLF destination, the kolf_2 and kolf_3 donors are dropped (A3).
  - CD4 halves are unit-level, as in Task 9's E4.
- **B5.** Noise is independent between halves. The O(1/n) bias of the log pseudobulk is not removed.
- **B6.** The other-line sets are:
  - **pack** (primary): the panel's ERA pack families, i.e. the pipeline's own sources.
  - **harness6** (sensitivity): the other 6 harness lines.
  - **rn3** (RN panels only): the other 3 of K562, RPE1, HepG2 and Jurkat.
  - **rn4_insample** (cite-check only): the destination is inside the fit, as in Molina & Zhang.
- **B7.** The additive fit is unweighted.
- **B8 (finding, see below).** The registered r_max formula assumes a ⊥ i. That fails when a comes from the other lines.

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Synthetic known components (balanced, control noise shared within a line) | within 5 % | Var_a, Var_b, Var_i, Var_rep, σ²_ε, r_max, ceiling all within 5 % | pass |
| Unbalanced design vs noise-free plug-in | within 5 % | pass; 6 extra seeds: Var_a, Var_i, Var_rep, r_max ≤ 1.8 %; Var_b up to 7 % (control-noise sampling, see jackknife) | pass |
| Identity Var_rep = Va + Vb + Vi + 2C_ai | exact | residual ≤ 1e-17 on real data | pass |
| Quadratic-form bootstrap = brute-force recomputation | equal | rtol 1e-8 | pass |
| Count-level (Poisson cells → log pseudobulk) | Var_rep, r_max within 5 %; σ²_ε = full-depth replicate variance within 5 % | pass | pass |
| Control-split alignment of the new extraction with Task 9 | bit-identical CTRL sums, all 5 seeds | identical (VCC_H1, CD4, KOLF, K562, HCT116, HEK293T) | pass |
| No prediction or scorer input read | none | grep test | pass |

## Results

Run: unit `aivc-paper-t12-components`, 243 s, peak RSS 4.9 GB (log `results/logs/paper-t12-components.log`).
Setup: 5 split seeds, B = 2,000 (seed 20261113) and 6,009 genes. Variances are per element, in units of
10⁻³ (log1p CP10K)². CIs are 95 % target-bootstrap intervals. Gene-block jackknife CIs, the seed range and the
Pearson-space versions are in the parquet.

| panel | set | T | lines | Var_a | Var_b | Var_i | 2C_ai | Var_rep | σ²_ε | ceiling_meas | r_max (reg.) | r_oracle | r_proj | r_max_orth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H1 | pack | 300 | 9 | 0.19 | 0.07 | 1.13 | 0.08 | 1.48 | 0.43 | 0.880 [0.860, 0.896] | 0.423 [0.395, 0.462] | 0.490 [0.465, 0.521] | 0.492 [0.465, 0.526] | 0.436 [0.412, 0.467] |
| H1 | harness6 | 298 | 6 | 0.36 | 0.11 | 1.24 | -0.19 | 1.52 | 0.43 | 0.883 [0.864, 0.898] | 0.555 [0.518, 0.605] | 0.442 [0.420, 0.468] | 0.447 [0.429, 0.470] | 0.523 [0.500, 0.553] |
| OR_HCT116X | pack | 232 | 3 | -0.21 | 0.01 | -0.18 | 0.58 | 0.20 | 1.53 | 0.338 [0.293, 0.374] | nan [nan, nan] | nan [nan, nan] | nan [0.024, 0.194] | 0.725 [0.683, 0.778] |
| OR_HCT116X | harness6 | 226 | 3 | 0.31 | 0.02 | 0.37 | -0.50 | 0.21 | 1.53 | 0.345 [0.303, 0.384] | 1.276 [1.136, 1.458] | 0.326 [0.278, 0.386] | 0.409 [0.364, 0.482] | 0.690 [0.667, 0.713] |
| KOLF | pack | 278 | 5 | 0.02 | 0.05 | 1.20 | 0.02 | 1.28 | 2.71 | 0.567 [0.554, 0.579] | 0.232 [0.164, 0.299] | 0.265 [0.194, 0.362] | 0.269 [0.187, 0.598] | 0.234 [0.167, 0.297] |
| KOLF | harness6 | 247 | 3 | 0.13 | 0.05 | 1.28 | -0.19 | 1.27 | 2.69 | 0.566 [0.552, 0.580] | 0.375 [0.349, 0.413] | 0.180 [0.143, 0.238] | 0.217 [0.182, 0.274] | 0.350 [0.328, 0.384] |
| RN_K562 | pack | 300 | 8 | 0.95 | 0.24 | 1.35 | -0.68 | 1.85 | 2.73 | 0.636 [0.608, 0.661] | 0.800 [0.752, 0.856] | 0.571 [0.546, 0.597] | 0.582 [0.561, 0.606] | 0.684 [0.667, 0.702] |
| RN_K562 | harness6 | 300 | 6 | 1.26 | 0.35 | 1.56 | -1.20 | 1.96 | 2.73 | 0.647 [0.620, 0.671] | 0.906 [0.856, 0.962] | 0.567 [0.544, 0.593] | 0.594 [0.574, 0.619] | 0.713 [0.698, 0.729] |
| RN_K562 | rn3 | 300 | 3 | 1.97 | 0.75 | 1.94 | -2.30 | 2.36 | 2.73 | 0.681 [0.657, 0.704] | 1.073 [1.023, 1.126] | 0.619 [0.597, 0.644] | 0.679 [0.658, 0.706] | 0.764 [0.753, 0.775] |
| RN_K562 | rn4_insample | 300 | 4 | 1.51 | 0.45 | 1.09 | -0.98 | 2.06 | 2.73 | 0.656 [0.631, 0.679] | 0.974 [0.933, 1.020] | 0.730 [0.713, 0.747] | 0.742 [0.728, 0.757] | 0.802 [0.792, 0.812] |
| RN_HepG2 | pack | 300 | 8 | 0.73 | 0.54 | 2.60 | 0.60 | 4.47 | 3.69 | 0.740 [0.706, 0.770] | 0.534 [0.512, 0.560] | 0.660 [0.642, 0.678] | 0.668 [0.648, 0.690] | 0.574 [0.553, 0.597] |
| RN_HepG2 | harness6 | 300 | 6 | 0.96 | 0.45 | 2.59 | 0.38 | 4.38 | 3.69 | 0.736 [0.703, 0.766] | 0.566 [0.545, 0.593] | 0.643 [0.623, 0.664] | 0.645 [0.624, 0.668] | 0.593 [0.575, 0.613] |
| RN_HepG2 | rn3 | 300 | 3 | 1.45 | 0.34 | 2.74 | -0.26 | 4.27 | 3.69 | 0.732 [0.698, 0.762] | 0.648 [0.618, 0.687] | 0.600 [0.579, 0.623] | 0.600 [0.580, 0.623] | 0.628 [0.612, 0.649] |
| RN_HepG2 | rn4_insample | 300 | 4 | 1.56 | 0.19 | 1.53 | 0.84 | 4.12 | 3.69 | 0.726 [0.693, 0.755] | 0.652 [0.629, 0.681] | 0.808 [0.793, 0.821] | 0.810 [0.794, 0.824] | 0.730 [0.717, 0.746] |
| RN_Jurkat | pack | 300 | 8 | 0.89 | 0.26 | 1.70 | -0.37 | 2.47 | 3.87 | 0.624 [0.598, 0.650] | 0.680 [0.643, 0.723] | 0.570 [0.544, 0.599] | 0.574 [0.550, 0.600] | 0.634 [0.616, 0.655] |
| RN_Jurkat | harness6 | 300 | 6 | 1.19 | 0.34 | 1.93 | -0.91 | 2.56 | 3.87 | 0.631 [0.605, 0.656] | 0.775 [0.732, 0.824] | 0.545 [0.521, 0.573] | 0.559 [0.538, 0.583] | 0.665 [0.648, 0.684] |
| RN_Jurkat | rn3 | 300 | 3 | 1.82 | 0.67 | 2.20 | -1.80 | 2.88 | 3.87 | 0.653 [0.628, 0.677] | 0.928 [0.883, 0.977] | 0.592 [0.570, 0.615] | 0.626 [0.608, 0.648] | 0.728 [0.716, 0.742] |
| RN_Jurkat | rn4_insample | 300 | 4 | 1.51 | 0.35 | 1.24 | -0.53 | 2.57 | 3.87 | 0.632 [0.606, 0.656] | 0.851 [0.816, 0.890] | 0.730 [0.711, 0.748] | 0.732 [0.715, 0.750] | 0.775 [0.763, 0.786] |
| RN_RPE1 | pack | 300 | 8 | 0.63 | 3.72 | 5.30 | 0.94 | 10.59 | 3.77 | 0.859 [0.847, 0.870] | 0.641 [0.614, 0.672] | 0.710 [0.691, 0.732] | 0.730 [0.710, 0.751] | 0.671 [0.646, 0.700] |
| RN_RPE1 | harness6 | 300 | 6 | 0.84 | 3.50 | 5.28 | 0.75 | 10.37 | 3.77 | 0.856 [0.844, 0.868] | 0.647 [0.620, 0.678] | 0.703 [0.684, 0.726] | 0.712 [0.693, 0.734] | 0.672 [0.647, 0.700] |
| RN_RPE1 | rn3 | 300 | 3 | 1.19 | 3.18 | 5.33 | 0.35 | 10.05 | 3.77 | 0.853 [0.840, 0.864] | 0.660 [0.631, 0.691] | 0.686 [0.665, 0.709] | 0.688 [0.667, 0.711] | 0.672 [0.648, 0.699] |
| RN_RPE1 | rn4_insample | 300 | 4 | 1.62 | 1.80 | 2.99 | 2.26 | 8.67 | 3.77 | 0.835 [0.821, 0.847] | 0.628 [0.606, 0.653] | 0.835 [0.827, 0.844] | 0.864 [0.853, 0.876] | 0.730 [0.713, 0.750] |


### Per-panel summary (primary set = pack)

| Panel | Targets | Cells per target | Measurement ceiling | r_max (registered) | r_oracle | r_proj |
|---|---|---|---|---|---|---|
| H1 | 300 | 300 | 0.880 [0.860, 0.896] | 0.423 [0.395, 0.462] | 0.490 | 0.492 |
| OR_HCT116X | 232 / 263 | 160 | 0.338 [0.293, 0.374] | not estimable (Var_a < 0) | — | — |
| KOLF | 278 / 279 | 200 | 0.567 [0.554, 0.579] | 0.232 [0.164, 0.299] | 0.265 | 0.269 |
| RN_K562 | 300 | 189 | 0.636 [0.608, 0.661] | 0.800 [0.752, 0.856] | 0.571 | 0.582 |
| RN_HepG2 | 300 | 54 | 0.740 [0.706, 0.770] | 0.534 [0.512, 0.560] | 0.660 | 0.668 |
| RN_Jurkat | 300 | 96 | 0.624 [0.598, 0.650] | 0.680 [0.643, 0.723] | 0.570 | 0.574 |
| RN_RPE1 | 300 | 85 | 0.859 [0.847, 0.870] | 0.641 [0.614, 0.672] | 0.710 | 0.730 |

### Findings that matter for Task 13 (flagged; no hypothesis was touched)

1. **The registered r_max formula is not a correlation when a comes from the other lines.**
   - **Why.** Var_rep = Va + Vb + Vi + 2C_ai, and C_ai is far from 0. It is negative when the destination responds
     more weakly than the other lines' average (amplitude slope (Va + C_ai)/Va: RN_K562 0.64, RN_Jurkat 0.79). It is
     positive when the destination responds more strongly (HepG2 1.41, RPE1 1.74, H1 1.22, KOLF 1.47).
   - **Consequence.** √[(Va + Vb)/Var_rep] then over- or understates what the oracle μ + a + b attains:
     - It exceeds 1 for OR_HCT116X with harness6 (1.28), RN_K562 with rn3 (1.07) and RN_K562 with rn4_insample
       (0.97).
     - It is undefined when Va < 0 (OR_HCT116X, pack).
   - **Alternatives.** r_oracle is the noise-free r of the oracle predictor itself. r_proj is the best-rescaled oracle,
     i.e. "the correlation achievable by an oracle that knows a_t and b_c" (definition 6 in words). Both are bounded and
     consistent with the definition's wording.
   - **For H1 specifically.**
     - With the registered formula, an oracle that rescaled a could beat r_max: the r_proj CI (0.465–0.526) lies
       above the point estimate of r_max (0.423), although the r_max CI (0.395–0.462) overlaps it. Task 13's
       "any method's r > r_max with CI excluding r_max → estimator invalid" rule would be at risk.
     - **Which ceiling H1 uses is a decision for the main session, logged before Task 13 looks at method r.**
2. **The OR_HCT116X panel has almost no reproducible signal in the μ-centred space.**
   - Var_rep is 0.20e-3 against σ²_ε of 1.53e-3, so reliability is 0.11 and the measurement ceiling 0.34.
   - Its 263 targets are non-essential, and K562 and CD4 barely respond to them. The conserved component is not
     distinguishable from 0 (Va = −0.21e-3).
   - **No transfer ceiling is estimable on this panel.**
3. **KOLF.** Va is about 0 (0.02e-3) with the pack lines (K562, HCT116, HEK293T, CD4 and 13 VCC_H1 targets). The
   ceiling is 0.57; r_max is 0.23, with a wide r_proj CI. Nearly all of KOLF's reproducible variance is
   interaction (frac_i ≈ 0.93).
4. **H1.**
   - The measurement ceiling is high (0.88; 300 cells per target).
   - The transferable share is small:
     - μ-centred Va 0.19e-3 and Vb 0.07e-3, against Vi 1.13e-3.
     - r_max 0.42, r_oracle 0.49, r_proj 0.49 (pack).
     - With harness6: r_max 0.56, r_oracle 0.44.
   - The choice of formula therefore moves the H1 ceiling by about ±0.07.
5. **The other-line set matters.**
   - *Pack* = 8–9 heterogeneous lines, including CD4, HEK293T and FENG.
   - *RN-only (rn3)* gives a larger Var_a and a larger |C_ai|.
   - The primary choice (pack = the pipeline's own sources) is defensible, but the ceiling is not unique.
6. **The target bootstrap misses the control noise in b.**
   - The control noise is shared by all targets of a line, so the target bootstrap cannot see it. The gene-block
     jackknife can, and on real data its CIs are narrower than the bootstrap's for every quantity reported.
   - So for the reported quantities the target bootstrap is the conservative interval. Both are in the parquet.

### Cite-check: Molina & Zhang (2026), Fig. 1c

- **Their numbers.** Fig. 1c gives b 27.8 %, a 29.4 %, i 23.5 % and noise 19.3 %. This implies
  r_max = √((27.8 + 29.4)/80.7) = √0.7088 = **0.842**. Supp. Note 1's 27.9 / 29.4 / 23.5 give 0.842. **The arithmetic
  in `input/notes/overlap-molina-zhang.md` is correct.**
- **Their setup.** The decomposition is in-sample over 4 lines (the destination is inside the fit; pooled over
  lines, C_ai cancels), on 2,000 HVGs and 638 balanced perturbations.
- **Our closest analogue.** rn4_insample on the 4 RN panels, in the MZ convention (C_ai dropped):

| Panel | r_max_orth |
|---|---|
| K562 | 0.80 |
| HepG2 | 0.73 |
| Jurkat | 0.78 |
| RPE1 | 0.73 |

  These are 0.04–0.11 below 0.84, which is consistent given 6,009 genes rather than 2,000 HVGs and the 300 panel
  targets.
- **What this means for the paper.**
  - **0.84 is an in-sample, line-pooled quantity, not a zero-shot ceiling.** Removing the destination from the fit
    (rn3) lowers the comparable r_max_orth to 0.63–0.76. The pipeline's own sources (pack) give 0.57–0.68 on the RN
    panels and 0.44 on H1.
  - The paper should cite 0.84 as "implied by MZ's in-sample decomposition". It should state that the
    leave-destination-out ceiling is lower, and by how much.

## Deviations (logged)

- **D1. FENG was not extracted for the 474 panel targets outside U.** Only 9 of the 474 have a ≥ 20-cell donor group,
  and the reader takes 1 h. FENG therefore enters only through the Task 9 cache (U targets): 14 H1 targets and 37 RN
  targets. FENG is absent for KOLF (A3 donors dropped; no KOLF panel target in U).
- **D2. The first extraction unit was stopped on the main session's memory request** (about 13 GB while reading
  HCT116).
  - HCT116 and HEK293T were re-extracted after Task 9 finished, with a row-batched reader (`extract_lean`).
  - Its output is bit-identical to Task 9's control sums, which also validates the lean reader.
  - VCC_H1, CD4, KOLF and K562 come from the first unit; they finished before the stop and were verified.
- **D3. Quantities added beyond the registered r_max.** r_oracle, r_proj and r_max_orth; the Pearson-space variant;
  the gene-block jackknife; the rn3, harness6 and rn4_insample sets.
- **D4. The registered r_max is ill-posed** (finding 1); not resolved here.
- **D5. The μ-centred measurement ceiling depends on the other-line set** (through μ). It moves by ≤ 0.011 between
  pack and harness6, but by up to 0.045 with rn3 (RN_K562: 0.636 vs 0.681). The Pearson-space ceiling is
  destination-only.

## What stays open

- Which ceiling H1 and Task 13 use: the registered r_max, r_oracle or r_proj. This must be decided and logged before
  any method r is computed.
- The noise-corrected r for predictions in Task 13 must use truth halves that are independent of the prediction's own
  dependence on the controls (Task 14 notes that predictions use all destination controls). The saved `dA`/`dB` use
  this task's destination split, not Task 14's C_a/C_b.
- The log-pseudobulk O(1/n) bias at half depth is not corrected (RN_HepG2: 27 cells per half).
- The OR_HCT116X transfer ceiling is not estimable with the available lines.

## Evidence base

- `results/rq1_components.{parquet,json,png}`
- `results/rq1_gene_cache/verify.json`
- `results/rq2_gene_cache/` (Task 9)
- `results/panels.json`, `results/g_cap.json`
- `input/notes/overlap-molina-zhang.md` (Fig. 1c numbers)
- [[rq2_params-assumptions]], [[rq3_build_preds-assumptions]]

