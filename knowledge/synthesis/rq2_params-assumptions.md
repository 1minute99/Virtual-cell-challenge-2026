---
title: "RQ2 Task 8 — Proposition-2 parameters on Tahoe (d, σ²_ε, Var(i), R²): assumptions and results"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ2 Task 8: Proposition-2 parameters on Tahoe

Script: `output/data-analysis/rq2_params.py` (library + CLI; Task 9 reuses `participation_ratio`, `descriptor_pcs`,
`half_assignment`, `filter_min_lines`, `halves_stats`, `aggregate`, `c_star_*`). Tests:
`output/data-analysis/tests/test_rq2_params.py` (18 pass). Results: `output/data-analysis/results/rq2_params_tahoe.json`,
`results/rq2_h2b.json`. Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md`, Task 8, definitions 3–4, H2b. Upstream:
[[rq2_tahoe_curve-assumptions]] (descriptor).

## Status: complete (2026-11-05, after the shards were restored)

The previous run stopped because all 157 frozen Tahoe parquet shards (15.2 GiB, revision 2dc57900…) had been
deleted from `vcc2026_data/reports/context_signal_2026-09-15/tahoe/raw/data/`. With the main session's approval they
were re-downloaded by `output/code/restore_tahoe_shards.py`. The script took the exact list (path, bytes, sha256)
from `tahoe/candidate_shard_plan.json` and cross-checked every entry against the HF LFS listing
(`source_documents/tahoe_data_directory.json`). It used no credentials, and all 157 of 157 sha256 checks passed.
The receipt is `tahoe/raw/restore_receipt_2026-11-05.json`. The full split-half run then completed in 818 s, with
a peak RSS of 9.5 GB.

**Bug fixed during the run (no change to any estimator).** The first full run aborted in `tahoe_half_responses`.
Control pool plate14 × CVCL_1571 has 2 DMSO wells of 1 cell each, so half A gets floor(1/2) = 0 cells and its
library is 0. That pool is used by no row; the 144 pools that are used have at least 30 cells, and at least 15 in
half A. The script had computed control pseudobulks for all 150 pools. It now computes them only for the pools
that rows reference. Regression test: `test_half_responses_ignore_unused_empty_control_pool`.

## What each estimate means

| Quantity | Meaning | Estimator |
|---|---|---|
| d | effective dimension of the 48 Tahoe context descriptors | (Σλ)²/Σλ² of the centred 48 × 5,955 `profile_adj` |
| Var(i) | noise-free variance of the line × drug interaction | ⟨i^A, i^B⟩ / residual d.o.f., with i^h = the residual of an additive (μ + a_t + b_c) fit made within half h |
| σ²_ε | noise variance of the full-depth target | ½·½‖i^A − i^B‖² / residual d.o.f. (the second ½ rescales half depth to full depth) |
| R² | share of Var(i) that is linear in the top-k descriptor PCs | fit on A and score on B per drug, ⟨(H_t − J_t) i^A, i^B⟩, pooled with the classical adjusted-R² d.o.f. |
| C*_theory, C*_min | Proposition 2 | d + 1 + d·σ²_res/(R²·Var(i)); d + 1 + d·σ²_ε/Var(i) |

Units: per (line, drug, gene) element, in (log1p CP10K)², on the model's 18,150-gene target axis. Var(i) and σ²_ε
share the same denominator, so every ratio that enters C* is free of the denominator.

## Assumptions

- **A1 depth rescaling.** Sampling variance scales as 1/cells for both the treated and the control pool.
- **A2 independent halves.** The halves are disjoint cells, and the DMSO wells are split as well.
- **A3 linear-Gaussian model.** This is the model of Proposition 2: a linear map from PCs to i.
- **A4 unweighted additive fit.** The two-way fit per gene is unweighted. Its leverage is identical across genes
  because they share the same observation mask.
- **A5 bootstrap scope.** The bootstrap (B = 2,000) runs over drugs × split seeds. d and the additive fit are held
  fixed; d's context jackknife is reported separately.
- **A6 gene axis.** All 18,150 target genes are used.
- **Response definition.** The response is the model target (dcp_prepare): log1p CP10K with the full native
  library as denominator, matched to the plate × line DMSO pool. Rows are combined with weight
  1/(1/n_t + 1/n_c).

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Descriptor = Task 6 `context_descriptors()` | identical | identical (48 × 5,955) | pass |
| Metadata row table = model `tahoe.npz` | 13,011 rows; same ctx/pert/repeat | identical; weight max rel. diff 5.9e-8; 4,431,233 cells | pass |
| Fold/subset helpers = `dcp_drug_lodo` | identical | identical | pass |
| Synthetic recovery (10 % tolerance) | Var(i) 1.0, σ²_ε 0.5, R² 0.6, d exact | 1.007, 0.500, 0.603 (adj.; raw 0.638); d to 1e-8; C*_theory 12.92 vs 13.00 | pass |
| Noise invariance of Var(i), R² | same at σ² 0.1 and 2.0 | within 10 % / 0.06 | pass |
| A1 (cell-level Poisson simulation) | split-half σ²_ε = replicate noise | +1.3 % | pass |
| Reader on a synthetic Tahoe-schema shard | sums exact, halves floor/ceil | exact | pass |
| Blinding grep | no curve or scaling file named | none | pass |
| Cell-level source on disk | 157 shards, frozen sha256 | 157/157 sha256-verified (restored) | pass |
| Cells scanned = frozen group table | every group's cell count | `finish_check` passed (4,431,233 cells) | pass |
| Full-depth rebuild = model target `tahoe.npz` | (ctx, drug) means equal | 204 of 8,568 checked, max abs diff 3.8e-7 | pass |
| Split-seed stability | small spread | C*_theory per seed 228.2–232.2; C*_min 80.45–80.72 | pass |

## Numbers obtained

### Split-half estimates

Setup: 5 split seeds; 185 drugs measured in ≥ 3 lines; 48 lines; 8,568 (line, drug) cells; 18,150 genes;
additive rank 232. CIs are percentile bootstrap intervals over drugs × split seeds, with B = 2,000 and seed
20261109.

| Quantity | Estimate | 95 % CI |
|---|---|---|
| Var(i) | 1.497e-3 | 1.143e-3 – 1.894e-3 |
| σ²_ε (full depth) | 4.302e-3 | 3.923e-3 – 4.689e-3 |
| σ²_ε / Var(i) | 2.87 | 2.28 – 3.70 |
| R²_adj, k = ⌈d⌉ = 21 (raw 0.650) | 0.347 | 0.327 – 0.369 |
| R²_adj, k = ⌊d⌋ = 20 (raw 0.630) | 0.338 | 0.318 – 0.358 |
| R²_adj, k = 2⌈d⌉ = 42 (176 drugs; raw 0.966) | 0.575 | 0.529 – 0.618 |
| **C*_theory, k = 21 (primary)** | **230.3** | **192.2 – 287.0** |
| C*_theory, k = 20 | 236.9 | 197.6 – 295.8 |
| C*_theory, k = 42 | 139.4 | 115.1 – 176.0 |
| **C*_min** | **80.6** | **68.5 – 97.6** |

Reading the table:

- **Noise.** Measurement noise is about 2.9 times the reproducible interaction variance, per element at full
  depth.
- **R².** Only about a third of the interaction is linear in the top 21 descriptor PCs after the
  degrees-of-freedom correction. The raw-to-adjusted gap is large because k = 21 is close to n_t ≈ 46.
- **C* against the pre-registered grid.** Both C* values lie far above the largest pre-registered C, which is 40.
  These numbers are facts. They are not an H2b verdict: the verdict needs C*_emp from Task 7.


- **d, primary (full descriptor): 20.56.** Leave-one-context-out range 19.94–21.38; jackknife SE 1.81. The top
  eigenvalue shares are 0.143, 0.101 and 0.053, so the spectrum is flat.
- **d, secondary (network input z):**
  - All 48 contexts: 20.56, equal to the primary as expected for a rank-preserving map.
  - C = 40 fold pools: 18.88 on the 40 training contexts (range 18.0–19.4) and 17.52 with the 8 held-out contexts
    projected.
  - Descriptive d of the C training descriptors: 4.6, 5.7, 8.6, 10.7, 14.4, 17.1 and 18.9 for C = 6, 8, 12, 16, 24,
    32 and 40. This is bounded by C − 1 and is not the population d.
- **PC counts for R²:** k = ⌈d⌉ = 21 is primary; ⌊d⌋ = 20 and 2⌈d⌉ = 42 are sensitivity analyses. With k = 42,
  only drugs measured in ≥ 44 lines are usable.
- **Structural floor (found by the previous agent; a fact about the formula).** For any σ²_ε/Var(i) ≥ 0 and any
  R² ∈ (0, 1], C*_min ≥ d + 1 ≈ 21.6 and C*_theory ≥ d/R² + 1 ≥ d + 1 ≈ 21.6. Two consequences:
  - If C*_emp is censored at > 40, the censored-case H2b rejection rule (C*_theory < 20) cannot fire, whatever
    the noise estimates are.
  - If C*_emp ≤ 43, a ratio C*_theory/C*_emp below 0.5 is impossible.

  With the estimates above (C*_theory 230, CI lower bound 192), there are two cases for the uncensored
  comparison. If C*_emp ≤ 40, the ratio is at least 4.8 at the CI lower bound, outside [0.5, 2]. A ratio inside
  [0.5, 2] would need C*_emp ≥ 115. The verdict is still Task 7's to give.

## Secondary sensitivities S1–S3 (pre-specified 2026-11-05; reported next to the primary, never replacing it)

Run: `rq2_params.py --dataset tahoe --sensitivity`. Results sit under `sensitivity` in `rq2_params_tahoe.json` and
`sensitivity_secondary` in `rq2_h2b.json`. The primary outputs are byte-identical to the primary run apart from the
`created` field. The sufficient statistics are cached in `results/logs/rq2_params_tahoe_extra_stats.npz`.

### Definitions

- **b_c.** The line coefficients of the same unweighted additive fit, computed per half and centred across lines.
- **V_b.** The noise-free line variance, ⟨(I − J) b^A, b^B⟩ / ((L − 1) G).
- **R²_b.** The adjusted R² of b on [1, top-k PCs] across the 48 lines, fitted on half A and scored on half B as
  for i.
- **Noise terms.** n_i = σ²_ε, and n_b = σ²_ε / n̄ with n̄ = 178.5 drugs per line, as pre-specified.

### Formulas

- **S1 (model target b + i).** The context gains of the two components add, so the residual and explained
  variances of Proposition 2 are summed:
  - C*_S1 = d + 1 + d·[(1 − R²_b)V_b + n_b + (1 − R²_i)V_i + n_i] / [R²_b V_b + R²_i V_i]
  - C*_S1,min = d + 1 + d (n_b + n_i)/(V_b + V_i)

  This is equivalent to a mixture of the component noise-to-signal ratios, weighted by each component's
  explained (or total) variance.
- **S2 (line component alone).** C*_b = d + 1 + d[(1 − R²_b)V_b + n_b]/(R²_b V_b) and
  C*_b,min = d + 1 + d n_b/V_b.
- **S3 (pre-specified line bootstrap).** Resample the lines with replacement, B = 1,000, seed 20261110. Each draw
  recomputes the centring, PCs, d, k = ⌈d⌉ and every component from per-line sufficient statistics. The additive
  fit and the half-specific b and i are held fixed (A5).

### Findings that change how these sensitivities can be read (flagged; no hypothesis is affected)

1. **The pre-specified b-noise model is ~37× too small.** The split-half noise of b is 9.0e-4 (jackknife CI
   6.4–11.7e-4), while σ²_ε/n̄ gives 2.4e-5. All drugs of a plate × line share one DMSO control pseudobulk, so
   control noise does not average out over drugs; it lands in b. The variants S1emp and S2emp use the empirical
   noise and are reported as flagged extras. In the same terms, the review's "Var(b) ≈ 1.37e-3" is about the
   noise-free V_b (5.6e-4) plus this noise (9.0e-4); that is my reading, not the review's.
2. **The S3 line bootstrap is invalid here.** Every point estimate lies outside its S3 percentile CI. With
   replacement only ~63 % of lines are distinct, so the rank-capped d drops to 11.6–17.0 against a point of 20.6,
   and duplicated rows inflate the adjusted R² to 0.55–0.72 against 0.35. The S3 CIs are kept in the JSON as
   specified, together with an `S3_validity` flag. The flagged replacement **S3b** is a delete-1 jackknife over
   lines: d, PCs, k and all components are recomputed without the deleted line, and the CI is point ± 1.96 SE.

### Results (Tahoe, 185 drugs × 48 lines, 5 split seeds)

All CIs are S3b line-jackknife CIs, which propagate d.

| Quantity | Point | S3b 95 % CI |
|---|---|---|
| d | 20.56 | 17.0 – 24.1 |
| V_i | 1.50e-3 | 1.40 – 1.59e-3 |
| σ²_ε | 4.30e-3 | 3.32 – 5.29e-3 |
| R²_adj, i | 0.347 | 0.296 – 0.398 |
| V_b | 5.62e-4 | 5.21 – 6.02e-4 |
| R²_adj, b | 0.299 | 0.232 – 0.365 |
| C*_theory (primary quantity) | 230.3 | 161.7 – 298.9 |
| C*_min (primary quantity) | 80.6 | 61.0 – 100.3 |
| **S1 C*_theory (b + i)** | **191.9** | **134.9 – 248.9** |
| S1 C*_min | 64.8 | 49.4 – 80.1 |
| **S2 C*_b** | **72.8** | **52.0 – 93.7** |
| S2 C*_b,min | 22.4 | 18.7 – 26.2 |
| S1emp C*_theory (flagged) | 218.1 | 150.9 – 285.4 |
| S1emp C*_min (flagged) | 73.5 | 55.2 – 91.9 |
| S2emp C*_b (flagged) | 180.4 | 113.0 – 247.9 |
| S2emp C*_b,min (flagged) | 54.6 | 39.9 – 69.2 |

The S3 bootstrap CIs, which are not valid here (see finding 2), are C*_theory 65.6–117.5 and C*_min 43.4–71.5.

**Reading.** Including the line component lowers C*_theory by about 17 %, to 192. If b's real noise is used, the
drop is only 5 %. Every C* except the noise-free bound S2 C*_b,min stays well above the largest C on the grid,
which is 40.

## Evidence base

`vcc2026_data/data/processed/drug_context_pretrain/{contexts,tahoe}.npz`; metadata in
`vcc2026_data/reports/context_signal_2026-09-15/tahoe/` (`candidate_shard_plan.json`, `candidate_shard_coverage.csv`,
`raw/metadata/*.parquet`, `gene_spec.json`); helpers `vcc2026_data/scripts/{dcp_data,simple_residual_prepare,context_signal_tahoe_effects}.py` (read-only).

## What stays open

- **H2b verdict.** It waits for C*_emp from Task 7, which is not read here (blinding).
- **Limitations:**
  - The estimate depends on the linear-Gaussian theory.
  - It covers only the Tahoe 5 µM subset (48 lines, 185 drugs).
  - The descriptor choice is fixed: the DMSO pseudobulk on 5,955 genes. d is large (≈ 21 of a possible 47), so
    the d + 1 floor dominates.
  - Noise is heteroscedastic across genes and cell counts, whereas A4 is exact only for homoscedastic noise.
  - With k = 21 on n_t ≈ 46 contexts, the adjusted R² carries a large degrees-of-freedom correction.
