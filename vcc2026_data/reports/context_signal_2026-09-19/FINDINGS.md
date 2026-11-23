# Does context carry transferable signal for perturbation response? — 2026-09-19

A measurement, not a model. Residual = observed response − leave-one-biology-out refit baseline
(`training_anchors` fits each biology's baseline from the OTHER biologies, so the residual is what
the baseline could not explain without seeing that biology). Contexts come from the corrected
matched-control pipeline: 16,263 distinct contexts vs 205 before.

## Verdict

1. **Measurable context → response signal?** Statistically detectable, practically nil. A
   cross-validated linear map from context deviation to residual deviation explains
   **R² ≈ 0.001** (0.1% of variance). The larger distance-based correlations are a sampling-noise
   artifact (see the confound below).
2. **Does it generalise to unseen biology?** Weakly and inconsistently. Nearest-context retrieval
   beats a random source by +0.044 cosine on average, but by −0.088 on K562, the largest and
   best-powered fold. Not a reliable channel.
3. **Which datasets have usable context variation?** Only Replogle-Nadig, which alone has
   per-target matched controls (K562 9,729 / Jurkat 2,201 / RPE1 2,122 / HepG2 2,010 distinct
   contexts). Every other source has exactly ONE control profile. But even RN's variation is
   dominated by sampling noise.
4. **Strong enough to justify another context-conditioned model?** **No.**
5. **Recommendation: abandon the context branch.** Details below.

## Test A — cross-source pairs, source identity held constant

For a fixed source pair, per-target matched controls make context distance vary across targets, so
this asks whether context predicts response *within* a pair.

| RN pair | n | raw Spearman | partial (controls residual magnitude) | z vs null |
|---|---|---|---|---|
| HepG2/Jurkat | 1,901 | 0.057 | 0.023 | 2.52 |
| HepG2/K562 | 1,980 | 0.052 | 0.008 | 2.24 |
| HepG2/RPE1 | 1,849 | 0.078 | 0.070 | 3.44 |
| Jurkat/K562 | 2,171 | 0.059 | 0.012 | 2.77 |
| Jurkat/RPE1 | 2,003 | 0.090 | 0.003 | 4.44 |
| K562/RPE1 | 2,089 | 0.128 | **−0.009** | 5.95 |
| **mean** | | **0.077** | **0.018** | 3.56 |

Significant raw correlation in all six pairs, but **77% vanishes** under the partial correlation,
and the strongest raw pair turns negative. Power is also poor: context distance varies by only
~0.3% within a pair (CV 0.0014–0.0048).

## Test A2 — within a single source (removes source identity by construction)

| source | n | CV R² | null R² | z | pairwise ρ | z |
|---|---|---|---|---|---|---|
| K562 | 9,729 | 0.0010 | −0.0003 | 9.62 | 0.240 | 211.8 |
| RPE1 | 2,122 | −0.0004 | −0.0011 | 0.60 | 0.138 | 97.7 |
| HepG2 | 2,010 | 0.0000 | −0.0025 | 3.13 | −0.002 | −1.8 |
| Jurkat | 2,201 | 0.0007 | −0.0011 | 2.56 | 0.099 | 92.8 |

The predictive test and the distance test disagree by two orders of magnitude, which is the tell.

## The confound that resolves the disagreement

Per-target context deviation is largely a proxy for how few cells that target had:

| source | ρ(‖context deviation‖, 1/effective n) |
|---|---|
| K562 | **0.685** |
| RPE1 | 0.612 |
| HepG2 | 0.632 |
| Jurkat | 0.572 |

A target measured in few cells gets a noisy matched control (large context deviation) AND a noisy
residual (large residual deviation). Neither causes the other. That shared factor inflates every
distance-based correlation — Test A's raw 0.077, and Test A2's ρ up to 0.240 — while the
cross-validated predictive map, which cannot exploit it, reports R² ≈ 0.001. HepG2's ρ of −0.002
against K562's 0.240 is the inconsistency you would expect from an artifact, not a mechanism.

## Test B — nearest-context retrieval on held-out biology

Borrow the same target's residual from the training source whose context is nearest.

| held out | targets | nearest ctx | random src | source mean | nearest − random |
|---|---|---|---|---|---|
| K562 | 10,104 | −0.224 | −0.136 | −0.461 | **−0.088** |
| RPE1 | 2,119 | 0.022 | −0.148 | −0.503 | +0.170 |
| HepG2 | 2,010 | −0.142 | −0.146 | −0.498 | +0.004 |
| Jurkat | 2,199 | −0.034 | −0.123 | −0.483 | +0.089 |

Mean +0.044, but the sign flips on the largest fold, and every absolute cosine is at or below
zero: borrowing another source's residual is uninformative about the target's true residual. (The
strongly negative "source mean" column is expected by construction — residuals are centred against
a cross-source anchor — and is not evidence about context.)

## Why the context branch cannot work on this data

The usable context variation is (a) confined to one source family, (b) ~0.3% of the between-source
scale, and (c) 57–69% explained by sampling noise. What remains predicts ~0.1% of residual
variance. For reference, the seed noise floor on the scoring harness is ±0.010 PDS; a channel this
small cannot move it.

This also explains the companion modelling result (`reports/context_model_2026-09-19`): the context
residual failed to beat "predict nothing" in 9 of 10 held-out-biology fits, and its best scored
configuration tied the refit baseline at +0.0003.

## Recommendation

**Abandon the context branch.** Redirect to the baseline and the perturbation representation:

1. **Submit the refit baseline.** It is the best local predictor (0.8005 vs b19's 0.7958) and has
   never been scored officially. One upload settles whether the neural residual helps at all.
2. **Attack expression accuracy, not context.** We score 0 on that metric (raw 15.5 against a
   context-mean baseline of ~1.10); it is 38% of the gap to the leaderboard leaders and nothing
   context-related touches it.
3. **If any context work continues, it needs different data, not a different model** — control
   profiles that vary within a source for reasons other than cell count, i.e. genuine biological
   replicates, which this dataset does not have (only 2.4% of rows have any replicate at all).
