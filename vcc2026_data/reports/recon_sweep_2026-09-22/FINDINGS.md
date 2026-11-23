# Is the delta scale calibration or metric tuning? — 2026-09-22

Motivated by the public leaderboard API (`reports/leaderboard_2026-09-22.json`): our raw expression-MSE
ratio is 14.7 against 0.58–0.79 for the top 30 (anything < 1 earns credit), and our composite gap
to the leader decomposes as MSE 0.072, PDS 0.048, NMAE 0.040, reach 0.033, fidelity 0.018,
Jaccard 0.005. Code `scripts/recon_sweep.py` (λ = 1 reproduces the frozen reconstruction count for
count, verified), `recon_sweep_lines.sh`, `recon_sweep_refit.sh`; results `results.json`,
`amplitude_slopes.json`. Official `cell-eval2 vcc2026` preset, seed 0, five held-out lines, the
existing leave-one-line-out deltas (refit_v2 and the S4 model film_v2 at its nested alpha).

## 1. Cell dispersion (λ) — negative
Shrinking predicted cells toward the control mean before applying the delta lowers the composite
at every scale on K562 (λ 0.5: −0.03 to −0.07; λ 0.25: −0.05 to −0.11): every gene becomes
"significant", fidelity and Jaccard fall. Dropped.

## 2. Scale, five lines, λ = 1

| scale | refit_v2 mean | S4 mean | S4 − refit | MSE ratio | PDS | fidelity | reach | NMAE | Jaccard |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.135 | 0.149 | +0.014 | 1.34 | 0.711 | 0.42 | 0.46 | 0.853 | 0.055 |
| 1.5 | 0.227 | 0.243 | +0.016 | 1.56 | 0.748 | 0.57 | 0.47 | 0.825 | 0.049 |
| 2.0 | 0.271 | 0.289 | +0.018 | 1.94 | 0.770 | 0.62 | 0.47 | 0.812 | 0.047 |
| **2.5** | **0.285** | **0.304** | +0.019 | 2.58 | 0.786 | 0.62 | 0.47 | 0.814 | 0.045 |
| 3.0 (submitted) | 0.281 | 0.300 | +0.019 | 3.38 | 0.797 | 0.62 | 0.47 | 0.827 | 0.044 |

Per-line optimum is identical for refit and S4: K562 2.0, HepG2 2.5, Jurkat 2.5, RPE1 3.0,
KOLF 3.0. The five-line optimum is 2.5, **+0.004 over 3.0** for both models; K562's +0.041 at
scale 2 is a one-line effect. S4's margin over the refit is +0.014 to +0.019 at every scale, so
the residual's value does not depend on the scale.

## 3. Empirical amplitude on the same held-out lines (`amplitude_slopes.json`)

| line | slope true~pred (S4) | r | RMS pred/true | variance-matching multiplier | composite-optimal scale |
|---|---|---|---|---|---|
| K562 | 0.38 | 0.28 | 0.74 | 1.35 | 2.0 |
| RPE1 | 0.92 | 0.29 | 0.32 | 3.2 | 3.0 |
| HepG2 | 0.82 | 0.40 | 0.49 | 2.0 | 2.5 |
| Jurkat | 0.43 | 0.24 | 0.55 | 1.8 | 2.5 |
| KOLF | 0.07 | 0.03 | 0.43 | 2.3 | 3.0 |

Two notions of calibration disagree. The least-squares amplitude (slope of truth on prediction)
is **< 1 on every line, mean 0.5**, because correlation with the truth is only 0.24–0.40: in
expression-error terms our deltas are already too large, and the MSE ratio is > 1 at every scale
(1.34 at scale 1). The variance-matching multiplier is **1.3–3.2, mean ~2**: cross-biology
averaging shrinks magnitude by that much. The composite's preferred scale tracks the second, line
by line (Spearman 0.9 over five lines; K562 lowest, RPE1 highest), not the first.

## Verdict
The scale is **amplitude restoration for the Wilcoxon-based direction metrics, not error
calibration**: it moves fidelity/reach/Jaccard, costs expression MSE at every value, and its
optimum differs by a factor 1.5 across lines in step with each line's amplitude deficit. **No
single global scale is supported by the data**; 2.5 vs 3.0 is +0.004 on five lines, inside seed
noise (earlier paired-seed sd ≈ 0.01), and the leaderboard's x1/x3/x4 ordering already bracketed
this. Expected effect of moving S4 from scale 3 to 2.5 on the official composite: at most a few
thousandths, sign uncertain. **Do not spend a slot on the scale.** What the sweep does establish:
S4's residual is worth +0.014–0.019 over the refit at any scale, and the expression-MSE gap to the
leaders (0.072) cannot be closed by scaling — it needs predictions that are correlated with the
truth (r 0.24–0.40 now), i.e. better deltas, not a different reconstruction.
