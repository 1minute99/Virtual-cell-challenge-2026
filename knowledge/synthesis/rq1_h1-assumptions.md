---
title: "RQ1 Task 13 — methods vs transfer ceilings and the H1 verdict: assumptions, checks and results"
type: synthesis
created: 2026-11-05
updated: 2026-11-05
status: draft
author: llm
---

# RQ1 Task 13: methods vs ceilings, H1

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1 and its falsification criterion, operational
definition 6, Task 13). Decision taken before any method r was computed (log, 2026-11-05):
- H1 is evaluated as registered on the H1 panel, with r_max = √[(Var a + Var b)/Var_rep].
- The plan's validity check is applied.
- r_proj and the per-panel ill-posedness of r_max are exploratory.

Upstream: [[rq1_components-assumptions]] (Task 12 components, halves and arrays) and
[[rq3_build_preds-assumptions]] (Task 14 predictions).

Files:
- Script: `output/data-analysis/rq1_ceiling.py`, with the subcommands `analyse`, `build`, `real` and `collect`.
- Results:
  - `results/rq1_ceiling.csv` (long format: panel × predictor × quantity, with 95 % CIs).
  - `results/rq1_h1.json`.
  - `results/rq1_ceiling.{png,svg}`.
  - `results/rq1_ceiling/<panel>_point_ci.json`.
  - `results/rq1_ceiling/build.json` and `H1_real_arms.npz` (the real-scorer packages).
- Logs: `results/logs/rq1_ceiling_analyse.log` and `results/logs/paper-t13-queue.log`.
- Queue: `output/data-analysis/rq1_ceiling_queue.sh`, run as unit `aivc-paper-t13-queue`.

## Estimator

The estimator is the same as Task 12's:
- **Space.** μ-centred response on the 6,009 genes and the Task 12 pack target set.
- **Halves.** The 5 split seeds of the destination truth (dA/dB) and of the other-line additive profile q (qA/qB).
- **Inner products.** Every product is a cross-half product.
  - Gram matrices of [δ; q; predictors] are symmetrised and averaged over the seeds.
  - The truth's measurement noise is removed.
  - A deterministic prediction keeps its own noise, which counts as method error.
- **Correlation.** r_P = S(P,δ)/√(S(P,P)·S(δ,δ)), where S(x,y) = Σ_t w_t⟨x_t^A − μ^A, y_t^B − μ^B⟩.
- **Definition 6.**
  - Attained fraction: r_m3 / r_max.
  - Share: E_tr/(E_tr + Var i), where E_tr = ‖(μ + â + b̂) − P‖²/N with the oracle's noise removed.
- **Bootstrap.** Targets are resampled with B = 2,000, using the same multinomial draws as Task 12 (seed 20261113).
  All centring is recomputed in each draw. The ceilings and the method r share the draws, so comparisons are paired.

## Assumptions (stated in the script header)

- **C1.** Truth, targets, genes, halves, μ, q, b and Var(i) are exactly Task 12's.
- **C2.** Predictions are aligned by label and by `panel_cols`. KOLF uses the kolf-fixed main build.
- **C3.** Predictions use all destination controls.
  - The fm-split diagnostic (`m3_fmsplit`) re-emits the pipeline delta with the control mean of the same truth half.
  - Budget and ew-cpm dependence is ignored (≤ 2 % of ‖L‖).
- **C4. Validity check.**
  - *Primary:* a method's own 95 % CI lies entirely above the r_max point estimate.
  - *Sensitivity:* the paired CI of r_method − r_max lies entirely above 0.
  - *Methods:* m3, pre-M3 and xmean. The architectures are unavailable.
  - *Diagnostics, not methods:* xmean_nf, xmean_naive, m3_fmsplit, and the oracle μ + â + b̂ (the ceiling's own
    predictor).
- **C5.** The verdict uses point estimates, as the rule is written. CIs are reported.
  - Sensitivity: the share after the optimal rescaling α = C/V.
- **C6. Real-scorer runs.**
  - They use the Task 17 emission (`rq3_corrections.emit_share`) and the same command (PDS only).
  - The oracle package sets the 6,009 Task 12 genes to expm1(o + log1p CP10K(fm)); 8,024 of 1.8 M entries are
    clipped at 0.
  - All other genes keep the pipeline's profile. The result is renormalised and carries the pipeline's library.
  - The emitted oracle delta has cosine 0.988 with o.
- **C7. Source-noise overlap.** Added after the first run showed a negative xmean share, and before any verdict was
  written.
  - q and every source-built prediction use the same source cells.
  - A deterministic prediction's cross terms with q therefore carry shared noise, which biases E_tr (and so the
    share) downward.
  - **xmean:** corrected exactly with its own halves. The uncorrected value is reported as `xmean_naive`.
  - **m3:** cannot be corrected, because no half-built prediction exists. The bound is
    share_adj = E_tr + κ·bias(xmean), with κ = the LS loading of m3 on the noisy xmean (κ_ls 0.92 on H1) or the norm
    ratio (κ_amp 1.30).
  - This is a plausible bound, not a strict one.

## Checks

| Check | Expected | Got | Pass |
|---|---|---|---|
| Var_rep vs Task 12 (7 panels) | rel < 1e-6 | ≤ 2.6e-10 | yes |
| r_oracle vs Task 12 | abs < 1e-6 | ≤ 6e-11 (OR: n/a) | yes |
| Bootstrap r_max CI vs Task 12 | identical draws | ≤ 7e-11 | yes |
| q rebuilt from the cached responses vs Task 12 q (full depth and seed-0 half) | equal | ≤ 6e-8 | yes |
| Stored m3 delta = delta recomputed from lfc, Od and fm | equal | ≤ 4.2e-7 | yes |
| Control mean = stored fm | equal | rel ≤ 4e-8 | yes |
| Shared-control effect (m3_fmsplit vs m3), H1 r | negligible | 0.3804 vs 0.3804 | yes |
| Pipeline package lfc = Task 17 M3_dest lfc | equal | 2e-7 | yes |
| Assumption: prediction shares no noise with q | — | **fails** for source-built predictions | see C7 |

## Results: H1 panel (registered)

| quantity | estimate [95 % CI] |
|---|---|
| r_max (registered) | 0.423 [0.395, 0.462] |
| r, pipeline M1 + M2 + M3 | 0.380 [0.348, 0.411] |
| r, pre-M3 | 0.382 [0.351, 0.412] |
| r, cross-line mean transfer (noise-corrected cross terms) | 0.289 [0.260, 0.316] |
| r, cross-line mean, own source noise removed (diagnostic) | 0.449 [0.414, 0.486] |
| r, oracle μ + â + b̂ (diagnostic) | 0.490 [0.465, 0.521] |
| **attained fraction r_m3 / r_max** | **0.899 [0.788, 1.004]** |
| **remaining-error share (registered formula)** | **0.137 [0.115, 0.166]** |
| share, source-noise bound κ_ls / κ_amp | 0.392 [0.363, 0.428] / 0.459 [0.426, 0.499] |
| share after optimal rescaling (α = 0.50) | 0.004 [−0.013, 0.029] |
| measurement ceiling | 0.880 [0.860, 0.896] |

**Validity check** (primary, own CI above r_max = 0.423). No method's CI lies above r_max:

| predictor | role | r CI | paired CI of r − r_max |
|---|---|---|---|
| m3 | method | [0.348, 0.411] | [−0.098, +0.001] |
| pre-M3 | method | [0.351, 0.412] | [−0.096, +0.003] |
| xmean | method | [0.260, 0.316] | [−0.188, −0.092] |
| xmean_nf | diagnostic | [0.414, 0.486] | [−0.036, +0.082] (not above) |
| oracle | diagnostic | [0.465, 0.521] | [+0.016, +0.114] (**above**) |

The ceiling's own predictor beats the registered formula. That is the Task 12 C_ai finding: a and i correlate when a
comes from the other lines.

**H1 verdict (registered rule): CONFIRMED.** The fraction is 0.899 ≥ 0.7, the share is 0.137 < 0.5, and the validity
check passes for methods. Caveats that must accompany it:
1. **If the oracle predictor is counted as a "method", H1 is untestable** (its CI excludes r_max).
2. The share is biased low (C7). With the source-noise bound it is 0.39–0.46, so the point estimate is still < 0.5,
   but the κ_amp upper CI reaches 0.499.
3. The fraction's CI spans 0.79–1.00.
4. The pipeline's amplitude is about 2× the least-squares optimum (α = 0.50). After rescaling, almost no error is
   left in the transferable components.

**Exploratory (r_proj = 0.492 [0.465, 0.526]).**
- Fraction 0.773 [0.722, 0.810], share unchanged.
- No method and no diagnostic exceeds r_proj (the oracle's paired CI is [−0.005, 0.000]).
- The registered rule applied with r_proj gives **confirmed**.

## Results: other panels (secondary, descriptive)

| panel | r_max | r_m3 | fraction | share (reg.) | share adj κ_ls | registered rule | r_proj rule (exploratory) |
|---|---|---|---|---|---|---|---|
| OR_HCT116X | not estimable | 0.231 | — | not estimable (Var i < 0) | — | untestable | untestable |
| KOLF | 0.232 | 0.058 | 0.250 | 0.725 | 0.810 | rejected (model-limited) | rejected |
| RN_K562 | 0.800 | 0.459 | 0.574 | 0.340 | 0.562 | inconclusive (→ rejected with adj) | confirmed |
| RN_HepG2 | 0.534 | 0.472 | 0.885 | 0.587 | 0.668 | inconclusive | inconclusive |
| RN_Jurkat | 0.680 | 0.437 | 0.643 | 0.565 | 0.674 | rejected (model-limited) | inconclusive |
| RN_RPE1 | 0.641 | 0.327 | 0.510 | 0.521 | 0.577 | rejected (model-limited) | rejected |

- **Validity on the other panels.**
  - No method exceeds r_max on any panel.
  - The oracle exceeds r_max on HepG2 and RPE1.
  - xmean_nf exceeds it on HepG2 (0.598 [0.577, 0.619] vs 0.534).
  - The verdicts are therefore fragile wherever r_max < r_oracle: HepG2, RPE1 and H1.
- **Cross-line mean transfer.** It beats the pipeline in r on HepG2 (0.484) and RPE1 (0.354).
- **KOLF.** The pipeline is close to no signal in the μ-centred space (r 0.058).
- **Architectures.** Not available or not comparable.
  - Only score JSONs remain in `arch_screen_2026-09-22{,_HepG2}`; the `<arm>.npy` arrays are deleted.
  - The targets do match the RN_K562/HepG2 panels (300/300).
  - But the arrays were on the official 18,533-gene axis in CP50000, against disjoint replacement controls.

## Real-scorer runs (H1, PDS; queued)

- **Unit:** `aivc-paper-t13-queue`.
- **Start condition:** it waits for `aivc-paper-t18-queue`, `aivc-paper-t10-queue` and `aivc-aivc-paper-t17-rest`
  to exit.
- **Order:** it builds the packages under `flock results/.celleval.lock`, then runs oracle and pipeline one after the
  other under the lock.
- **Results:** they go to `rq1_h1.json` → `real_scorer`, including oracle − pipeline per-target PDS with a paired
  bootstrap CI.
- **Reproducibility check:** the pipeline run repeats Task 17 M3_dest (0.91013).
- These runs do not change the H1 verdict.

## What stays open

- **C7 for m3 cannot be removed without rebuilding the pipeline on source halves.** The share is therefore a lower
  value with a plausible upper bound, not an unbiased estimate.
- **Whether the oracle counts as a "method" in the validity check** is for the main session to decide (and log).
  - Counted: H1 is untestable.
  - Not counted: H1 is confirmed.
- **The verdict depends on r_max's ill-posedness** (Task 12 finding 1). The exploratory r_proj reading agrees
  (confirmed).

## Evidence base

- `results/rq1_ceiling.csv`, `results/rq1_h1.json` and `results/rq1_ceiling/*_point_ci.json`.
- Task 12: `results/rq1_components/<panel>_pack.npz` and `results/rq1_components.json`.
- Task 14: `results/preds/<panel>_{m3,prem3}.npz`.
- Task 17: `results/rq3_corr/real/M3_dest/result.json`.
- [[rq1_components-assumptions]], [[rq3_build_preds-assumptions]]
