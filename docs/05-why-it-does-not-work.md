# 05 — Why it does not work

This file collects the explanations behind the negative results in `04-what-did-not-work.md`. Part A covers the data
and the task; Part B covers the scorer properties that decided which changes could pay at all. Measurement labels
follow `04-what-did-not-work.md` (O official leaderboard, R real local scorer, X offline or proxy). The post-challenge
analyses (November 2026) followed a pre-registered plan; exploratory results are marked. M1, M2 and M3 are the fusion
transfer, the agreement amplitude and the control-PC deflation (`01-challenge-and-final-system.md`).

## Part A — Empirical and theoretical reasons

### A1. The noise ceiling is low, and lowest where the score is decided

**Few replicates, low reproducibility.** Only 2.4 % of the rows in the training table had any true replicate (the
Jiang screen); KOLF's three batches screened disjoint gene sets. Where reproducibility could be measured, two
replicates of the same line and stimulus agreed at median r = 0.209 (range 0.04–0.41).

**Most of a line's deviation is noise.** For each training biology, the deviation of its responses from the
cross-biology mean was split into signal and sampling noise by fitting ‖δ_p‖² = a + b / n_eff,p across perturbations
and taking the intercept as signal. The reproducible share averaged 0.277: HepG2 0.681, CD4 0.565, Jurkat 0.344, FENG
0.199, HEK293T 0.192, K562 0.180, RPE1 0.165, HCT116 0.163, KOLF 0.007. A model asked to learn line-specific
deviations from public data was therefore fitting about 72 % noise on average and more than 99 % on KOLF; 94 % of the
deviation was perturbation-specific rather than a generic line shift.

Sources: `vcc2026_data/reports/data_diagnosis_2026-09-19/FINDINGS.md`; `vcc2026_data/reports/basal_context_2026-09-20/FINDINGS.md`.

**Destination ceilings.** Split-half estimates on the evaluation panels (6,009 genes, noise-corrected) give the
correlation a perfect predictor could reach against the measured truth:

| Panel | Cells per target | Reliability | Measurement ceiling |
|---|---|---|---|
| H1 | 300 | 0.775 | 0.880 [0.860, 0.896] |
| RN_RPE1 | 85 | 0.737 | 0.859 |
| RN_HepG2 | 54 | 0.548 | 0.740 |
| RN_K562 | 189 | 0.404 | 0.636 |
| RN_Jurkat | 96 | 0.390 | 0.624 |
| KOLF | 200 | 0.322 | 0.567 |
| OR_HCT116X | 160 | 0.114 | 0.338 |

Sources: `output/data-analysis/results/rq1_components.json`; `knowledge/synthesis/rq1_components-assumptions.md`.

**Cross-line agreement for the same target is far below the ceiling.** Target-centred correlations of the same
target's response: CD4 resting vs 8 h stimulated (the same cells and donors) 0.270; RN K562 vs RPE1 0.163, HepG2 vs
Jurkat 0.174 (same lab); Orion HCT116 vs HEK293T 0.063; across sources 0.013–0.071. Two independent K562 screens agreed
at cosine 0.28, i.e. one measurement was about 30 % signal; noise attenuated cross-line agreement 8–12×. Holding out
one context instead of a whole biology inflated apparent transfer about six-fold (0.1332 vs 0.0218), because the CD4
states are near-replicates.

Sources: `vcc2026_data/reports/data_diagnosis_2026-09-19/FINDINGS.md`; `AGENTS.md`.

**Weak responders decide PDS and do not transfer.** On three RN lines (300 shared targets), the within-line
split-half cosine was 0.108 for the weakest quartile of targets against a cross-line cosine of 0.013–0.025, a
transferable fraction of 0.12–0.23; for the strongest quartile the fraction was 0.49–0.74. On H1 the gap between our
PDS and the split-half ceiling (0.895 vs 0.963) sat 66 % in the weakest quartile (0.686 vs 0.866). These targets have
responses that are measurable inside a line but barely shared between lines.

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§11, §12).

### A2. What can and cannot be transferred

The response of target t in context c is decomposed as δ_{c,t} = μ + a_t + b_c + i_{c,t} + ε: a shared response, a
conserved target component, a context component, a context × target interaction and noise (Molina & Zhang 2026,
`molina-2026-response-decomposition`). For an unseen line, μ and a_t can be estimated from other lines; b_c and
i_{c,t} would have to be inferred from the line's control cells. On the H1 panel, with the pipeline's nine source
lines as "other lines", the noise-free variances were (×10⁻³ per element):

| Var(a) | Var(b) | Var(i) | 2 Cov(a, i) | Reproducible total | Noise σ²_ε |
|---|---|---|---|---|---|
| 0.19 | 0.07 | 1.13 | 0.08 | 1.48 | 0.43 |

The interaction carried 76 % of the reproducible variance on H1 and about 93 % on KOLF. The noise-free correlation
of an oracle that knows μ + a + b was 0.490 [0.465, 0.521] (best-rescaled 0.492 [0.465, 0.526]), against a measurement
ceiling of 0.880. The shipped recipe (M1 + M2 + M3) reached r = 0.380 [0.348, 0.411], i.e. 0.77 [0.72, 0.81] of the
bounded ceiling. This last ratio is exploratory: the registered ceiling formula assumes a ⊥ i, which fails when a is
estimated from other lines, so the registered hypothesis was declared untestable. Molina & Zhang's in-sample
decomposition implies a ceiling of about 0.84; with the destination left out the comparable values were 0.63–0.76
(three RN lines) and 0.44 on H1 with the pipeline's own sources.

Sources: `output/data-analysis/results/rq1_components.json`; `output/data-analysis/results/rq1_h1.json`;
`knowledge/synthesis/rq1_components-assumptions.md`; `knowledge/synthesis/results-summary.md`.

**Context effects are real but not predictable from controls.** If a line's response was known from an independent
batch, prediction improved by +0.1953 cosine over the context-free prediction (569/649 targets; a random-split control
gave 0.95–1.00× the same gain, so it is not technical leakage), worth +0.27 to +0.48 composite. On a strict
unseen-context benchmark the best control-based context model recovered 22.3 % of the oracle headroom in cosine, and
on the harness about none of it in PDS (`04-what-did-not-work.md`, (b)). Molina & Zhang likewise report that the
interaction is not recovered zero-shot.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `knowledge/synthesis/related-work-positioning.md`.

**Within the transferable part, the pipeline was near what the sources allow.** The sum side reached 80–91 % of the
best per-target non-negative mix of our own sources (a restricted-model bound, not an information ceiling), and that
mix reached only about 45–60 % of the truth's noise ceiling on the strong responders. On RN_Jurkat the recipe beat
oracle transfer that used the other RN lines' measured responses for the same targets (PDS 0.9078 vs 0.8974 for the
mean of three lines).

Sources: `AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§12, §20.7–20.8).

### A3. The number of contexts

**In-challenge evidence.** The context branch of the learned residual acted as a source lookup: on the KOLF fold there
were 58 distinct control profiles across 92,029 rows, so the "context input" was a categorical label. With nine
training biologies at most eight context vectors were available to learn a map from. Adding 31 real contexts (CD4
donor × state, HipSci lines) left discrimination and PDS flat from 6 to 35 contexts, but those additions came from two
cell types.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/basal_context_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md`.

**Post-challenge learning curve (pre-registered, confirmed).** On Tahoe-100M data (Zhang et al. 2026,
`zhang-2026-tahoe`; 48 drug-treated cancer lines), a FiLM model's context pathway was compared with a matched
constant-descriptor ablation while the number of training contexts C varied (3 seeds × 6 folds × 7 values of C,
cluster bootstrap). The gain is in held-out residual-MSE units; negative means the context pathway hurts:

| C | 6 | 8 | 12 | 16 | 24 | 32 | 40 |
|---|---|---|---|---|---|---|---|
| Gain (×10⁻⁴) | −3.93 | −4.86 | −4.20 [−5.15, −3.32] | −2.80 | −1.08 | −0.32 | +0.14 [−0.06, +0.35] |

Every interval up to C = 32 lies below zero. The break-even is C* = 37.4 [33.3, > 40]; 8.75 % of bootstrap draws were
censored above 40. The registered expectation (harm at small C, break-even at 24 or more contexts) was confirmed.

Sources: `output/data-analysis/results/rq2_h2a.json`; `knowledge/synthesis/theory-context-complexity.md`.

**Public genetic data sit below that.** The source registry contains 11 distinct cell types with a usable genome- or
essential-scale CRISPRi Perturb-seq response (K562, HCT116, HEK293T, primary CD4+ T, KOLF2.1J, the HipSci iPSC pool,
HepG2, Jurkat, hTERT-RPE1, DLD-1, H1), or 25 if HipSci donor lines are counted separately. The genetic analogue of the
minimum break-even was 16.45 [16.09, 16.82]: above 11, below 25. Because the same estimator overestimated the Tahoe
break-even 2.2-fold, this is weak support only.

Sources: `output/data-analysis/results/rq2_h2c.json`; `output/data-analysis/results/families.json`;
`knowledge/synthesis/theory-context-complexity.md`.

**Identifiability.** If the learned correction is linear in a descriptor x_c of the control profile, î_{c,t} = W_t x_c,
and the C training descriptors are centred, their span has rank at most C − 1. Two maps that agree on that span fit
the training data equally well and differ for any new line with a component outside it, so at most C − 1 descriptor
directions can be learned, whatever the model capacity. For a nonlinear map the function is constrained only at C
points, and under smoothness the minimax rate is the nonparametric rate with sample size C (Stone 1982,
`stone-1982-optimal-rates`). This is the context-count analogue of the environment-count result for invariant risk
minimisation (Rosenfeld et al. 2021, `rosenfeld-2021-risks-irm`) and is cited as such, not as a new theorem.

**A closed-form break-even was ×6 too pessimistic.** Under a linear-Gaussian model with OLS, a d-dimensional
descriptor and explained fraction R² of the interaction variance, the expected gain over predicting no interaction is
R²·Var(i) − σ²_res·d / (C − d − 1), giving C* = d + 1 + d·σ²_res / (R²·Var(i)) (Breiman & Freedman 1983,
`breiman-1983-how-many-variables`; Rosset & Tibshirani 2020, `rosset-2020-random-x`; Dicker 2013,
`dicker-2013-equivariant-prediction`). With Tahoe parameters estimated blind to the curve (d 20.6, noise/interaction
ratio 2.87, R²_adj 0.35) it predicts C* = 230 [192, 287], 6.2 times the measured 37.4; even its R² = 1 lower bound
(80.6) exceeds the measurement, so the pre-registered test rejected the formula. A plausible but untested reading: the
learner is regularised and stops early, predicts through a shared rank-32 basis, and much of its target is the
per-line component b_c, which is cheap to learn. What survives is the identifiability argument and the measured curve.

Sources: `output/data-analysis/results/rq2_h2b.json`; `knowledge/synthesis/theory-context-complexity.md`.

**More source families help the transfer only slightly.** During the challenge, random subsets of source families with
every part rebuilt gave PDS ≈ 0.796 + 0.048 ln k on H1, about +0.003 composite per added family. With the learned S4
part held fixed, the post-challenge curves were shallower (+0.004 to +0.017 PDS per ln k); no leave-one-family-out
interval excluded 0, and the per-context value on H1 was +0.00086 [−0.00045, +0.00255] composite. The S4-free curve
was 1.6–15× steeper (4–8× on four of the seven panels), so most of what an extra family adds was already carried by
S4; most of the remaining value came from K562. Small public screens added nothing because they overlap a 300-target
panel at chance.

Sources: `AGENTS.md`; `output/data-analysis/results/rq2_source_curve.json`; `knowledge/synthesis/rq2_source_curve-assumptions.md`.

The literature agrees: context-conditioned models lose to a cross-line mean when trained on three lines (Piao et al.
2026, `piao-2026-adaptive-contexts`); scalar-readout learning curves cross a context-free baseline at
dataset-dependent points, about 6 lines on LINCS and 170 on GDSC (Ji et al. 2024, `ji-2024-universal-prediction`);
context coverage rather than scale has been argued to bind (Dibaeinia et al. 2026, `dibaeinia-2026-context-scale`).

Sources: `knowledge/synthesis/related-work-positioning.md`.

### A4. Same-line data is worth more than modelling

On the H1 panel (real scorer, 300 targets × 300 cells), replacing transferred responses by those of H1GW, a
genome-scale CRISPRi screen of the same line from another lab (Pan et al. 2026, `pan-2026-global-cell-state`),
raised the composite from 0.3151 to 0.3841 at weight 1 (+0.0690; 0.3887 at weight 0.5, 0.3956 at 0.75). The largest
single modelling gains were M2 at +0.0065 on this panel and the best seed-mean architecture at +0.0075 (group (a) of
`04-what-did-not-work.md`): ratios of about 10.6 and 9.2. During the challenge H1GW reached a strong-target cosine of
0.49 against our 0.34, and about 70 % of its value was on the DE side. Related lines helped less and could hurt
(KOLF ← H1GW +0.021 at weight 0.5; OR_HCT116 ← DLD-1 +0.005 at 0.5 and −0.035 at 1). No usable public same-line screen
existed for the validation lines (Jurkat, HeLa, CAL-33) or for the three final lines.

Sources: `output/data-analysis/results/rq4_h4.json`; `output/data-analysis/results/rq2_arch/rq2_arch_seeds_summary.json`;
`knowledge/synthesis/results-summary.md`; `AGENTS.md`; `vcc2026_data/reports/final_2026-10-27/README.md`.

### A5. Where the remaining error lives

**On the destination's own control axes (pre-registered, confirmed).** On H1, the top 8 principal components of the
destination's control cells (8 of 6,014 gene dimensions) held 0.0842 [0.0711, 0.0983] of the energy of the transfer
residual. The same number of PCs from other lines' controls held 0.0659, from gene-permuted controls 0.0060, and a
random subspace 0.0014 (≈ k/G). The excess over other lines' PCs was 0.0182 [0.0092, 0.0267] and over permuted
controls 0.0782 [0.0653, 0.0922]; a pseudo-target negative control reached ≤ 2 % of these. The excess over permuted
controls was positive on 7/7 panels, and the share was largest for weakly predicted targets.

**But the margin is small once PC-estimation noise is matched (exploratory).** The registered design built the
destination PCs and the other lines' PCs from different numbers of cells. At matched cell numbers the H1 excess over
other lines halved to 0.0093 [0.0003, 0.0176], and at a common 267 or 800 cells its interval included 0. At matched
numbers it was positive on 7/7 panels; the negative values on RN_K562 and RN_Jurkat in the registered run were an
artefact of the mismatch.

Sources: `output/data-analysis/results/rq3_h3a.json`; `output/data-analysis/results/rq3_h3a_matchedn.json`;
`knowledge/synthesis/rq3_h3a-assumptions.md`; `knowledge/synthesis/rq3_matchedn-exploratory.md`.

On strong targets the truth carried about 2–6 times more energy on the destination's top-24 control PCs than our prediction
(H1 0.188 vs 0.070; KOLF 0.119 vs 0.021), yet amplifying the prediction along those axes lowered PDS on every panel.
The location of the error was known; its content was not.

Sources: `vcc2026_data/reports/final_2026-10-27/README.md`.

**The M3 correction against alternatives was inconclusive.** Offline, over 7 panels, deflating the weaker half of
targets along the destination's control PCs (M3) beat every alternative on mean ΔPDS:

| M3 minus | ΔPDS (7-panel mean) | Panels > 0 |
|---|---|---|
| no correction | +0.0131 [+0.0101, +0.0160] | 7/7 |
| other lines' control PCs | +0.0058 [+0.0039, +0.0076] | 7/7 |
| permuted controls | +0.0131 [+0.0101, +0.0159] | 7/7 |
| PCs of the prediction | +0.0065 [+0.0036, +0.0095] | 7/7 |
| Systema-style mean-response removal (cf. Viñas Torné et al. 2026, `vinastorne-2026-systema`) | +0.0145 [+0.0108, +0.0179] | 6/7 |

On H1 with the real scorer, however, the Systema-style arm beat M3 (PDS 0.9178 vs 0.9101), so the strict criterion —
the same sign on H1 against all five — failed and the verdict is inconclusive. At matched cell numbers the margin
over other lines' PCs was unchanged (+0.0058 [+0.0042, +0.0075]), with H1 the smallest (+0.0005).

Sources: `output/data-analysis/results/rq3_h3b.json`; `output/data-analysis/results/rq3_h3b_matchedn.json`;
`knowledge/synthesis/rq3_matchedn-exploratory.md`.

**The leaderboard gap was the direction of strong targets.** At the end of validation, the 0.30–0.31 tier led our
best entry by +0.049 to +0.058 composite, of which the MSE member supplied +0.023 to +0.028. From the MSE formula,
their implied energy-weighted cosine with the truth was about 0.4 against our 0.2, with the amplitude at the budget
kink for both. That is the population where same-line data helps (A4) and transfer saturates (A1–A2).

Sources: `vcc2026_data/reports/final_2026-10-27/README.md`.

## Part B — Scorer mechanics that shaped what could work

**B1. PDS is a discrimination rank.** For each target, PDS ranks the cosine distance from the prediction to its own
truth among all 300 targets' truths, after every panel target gene is removed. The on-target knockdown therefore
counts for nothing, and a generic all-target response raises ordinary cosine (about 6.6×) while lowering PDS. On the
KOLF harness a target-shuffled delta kept 79 % of the composite but only 3 % of PDS, so PDS was the one member that
discriminated targets. Cosine and PDS repeatedly moved in opposite directions (HepG2: +0.0294 cosine, −0.0727 PDS),
and PDS-only development rankings reversed on the full preset three times.

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `AGENTS.md`; `vcc2026_data/reports/drug_context_pretrain_2026-09-20/FINDINGS.md`.

**B2. The bulk comparator depends only on group sums.** The scorer's pseudobulk is log1p(TS · S_g / Σ_g S_g), computed
from each target's summed counts; PDS and MSE are deterministic functions of the 300 group totals. Nothing done to the
distribution of counts across a target's 400 cells can move them. The four DE members (NMAE, fidelity, reach,
Jaccard) come from Wilcoxon tests on the emitted cells, so they can be moved with the sums pinned — which is why the
emission work paid (`03-what-raised-the-score.md`) and why no sum-side "composition mismatch" existed to exploit.
Fidelity is k / max(n_pred, n_conf), so over-calling only dilutes it slowly while under-calling collapses it, and the
right number of calls depends on the test data's power (median 145 real DE genes per target on H1, 27 on KOLF).

Sources: `AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §16.2, §25.4).

**B3. PDS and MSE lose in disjoint target populations.** On H1 (S16 structure, exact offline), the per-target cosine
with the truth had median 0.085. The strongest three deciles (cosine 0.32 / 0.26 / 0.18) held 97 % of the MSE
numerator and were already at PDS 0.996 / 0.986 / 0.985; the weakest three (cosine 0.025 / 0.013 / 0.008) had PDS
0.861 / 0.738 / 0.672; 126 of 300 targets were already ranked first. No single lever raises both. Each reduces to the
per-target correlation, which is bounded by the sources in both populations (A1).

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§25.4).

**B4. The MSE member has a kink.** Its numerator is ‖μ_pred − μ_real‖² − r · min(C_pred, C_real) − C_real, with a
rebate r = min(1, B_pred / claim), where B_pred is the predicted response energy and the per-target claim is capped by
the truth's own noise. Below the kink extra energy is free and MSE falls with ⟨prediction, truth⟩; above it MSE rises
like a cliff. The kink is set by the test data's noise; in units of the base sum amplitude it was about 1.0 on H1,
1.0–1.1 on A/B/C, 1.1 on K562, 1.6 on Jurkat, 1.8 on RPE1, 2.0 on HepG2 and 2.8 on KOLF. On H1 our plug-in error was
worse than predicting the control (25,410 vs 24,503), so the whole MSE credit was rebate; approximately
MSE ≈ 1 − 2κρ√(B / claim) (H1: 0.748 predicted, 0.766 measured). Consequences:

- per-target amplitude is coupled to the global budget: a greedy per-target oracle that zeroed 125 targets dropped the
  rebate from 0.89 to 0.287 and lost;
- budget changes are panel-specific: the same +7 % was −0.181 on OR_X and +0.0075 on RN_Jurkat; both official probes
  away from the shipped margin lost (S17u −0.0462, S17d −0.0073);
- extra per-cell dispersion cannot raise the claim, because the cap already binds on every target (263/263 on OR_X).

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §8, §24.6, §25.1, §25.4); `AGENTS.md`;
`output/data-analysis/results/official_record.csv`.

**B5. The offline evaluator is exact only for budget-invariant changes.** Two budget-invariant amplitude arms on H1
reproduced the real scorer's ΔPDS and ΔMSE to four decimals (+0.0023 / −0.0065 and +0.0038 / −0.0118). A
budget-changing arm sign-flipped MSE (offline −0.0211, real +0.0055), because the evaluator uses the real panel's
jackknife noise in place of the emitted cells' own while the cap binds. Earlier, matching candidates to the base
budget instead of their own had biased offline MSE by a constant +0.0103. Orion's offline MSE was unreliable in
general: the kink-factor probe claimed +0.0768 where the scorer gave −0.181.

Sources: `AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§17, §24.6, §25.5–25.6, §25.10).

**B6. The J proxy is blind to NMAE.** J = 2.23 ΔPDS − 1.05 ΔMSE (≈ 6 × Δcomposite) was derived for arms that leave the
DE side alone. A per-target amplitude arm moves NMAE through the composition pin that couples the sums to the DE log
fold changes. On OR_X the arm `tzb0.25` had J +0.0398 and a real composite change of −0.0108, the whole gap being NMAE
(+0.0373). On H1, where NMAE moved by −0.0001, J/6 = +0.0020 matched the real +0.0022. Earlier amplitude screens had
also used J; the shipped M2 survived because the official NMAE behaved like H1's, not because J was right.

Sources: `AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§25.8).

**B7. Harness panels are not the official task.** The RN panels are essential-gene screens on about 8,000 genes with a
shared response of 35–66 % of the signal, against about 1 % on A/B/C. The KOLF harness mispredicted the sign of the
uploads it was used to justify (B11, B15). Only H1's NMAE behaved like the official one; on shallow panels one target
moved NMAE by −0.018 and the emission seed alone by −0.012. For source-family composition the stem-line H1 panel gave
the wrong sign (dropping KOLF: H1 J −0.0114, Orion J +0.0368, official +0.0061). Until 2026-10-06 the harness also
used three source families (HEPG2, RNJURKAT, RPE1) that the shipped specification lacked; scoring the shipped set cost
−0.0023 (H1) and −0.0013 (OR_X), one candidate contributor to harness over-prediction.

Sources: `AGENTS.md`; `vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`
(§18.3, §25.19).

**B8. Local gains shrink on the leaderboard.** The one valid case for a learned context model, S4 against the refit
S1, went from +0.0247 locally (15/15 seeds) to +0.00131 officially, a ratio of 0.053 with the sign kept. Other
transfers ranged widely: S10c matched its harness gain (+0.0247), S8 kept 63 %, S12 about a third, S11 kept its
sum-side gain but not its DE-side fixes (+0.00296 against +0.0044 to +0.0064 forecast), S10f went from +0.006 to
−0.0015, and S15 exceeded its forecast (+0.0106 against +0.005). Hence no uploads below an expected +0.005 from
2026-09-27.

Sources: `output/data-analysis/results/e4_official.json`; `knowledge/synthesis/e4_official-assumptions.md`; `AGENTS.md`.

## Lessons

- **Measure the noise ceiling and the transferable share before modelling.** On these panels the interaction
  dominates the reproducible variance and is not transferable zero-shot; capacity cannot buy what the data do not
  contain (A1–A2).
- **Count contexts, not rows.** The effective sample size for a context map is the number of distinct contexts; with
  about 11 public cell types a learned context pathway sits on the harmful side of a curve whose break-even was
  about 37 on drug data (A3).
- **Same-line or closely related data first.** One same-line screen was worth about nine times the largest modelling
  gain; related lines help less and can hurt at full weight (A4).
- **Model only what the scorer reads.** PDS and MSE see group sums; DE members see cells. PDS and MSE are decided in
  different target populations, so a single reshaping lever trades one for the other (B2–B3).
- **Use an evaluator only where its exactness condition holds.** Price PDS offline; price budget or amplitude changes
  on the real scorer; never decide an amplitude arm on an NMAE-blind proxy (B4–B6).
- **Controls must test identity, not magnitude.** Wrong-context, shuffled-context, zero-context and label-shuffled
  controls each caught a false positive; look at the worst line, not only the mean (`04-what-did-not-work.md`, (b)).
- **Treat suspiciously clean or exactly zero results as possible bugs.** Scrambled control profiles invalidated two
  uploads, the foundation-model context comparison and the gene side of the drug-pretraining experiment.
- **Expect an order-of-magnitude shrinkage from local to official** and adopt a change only with a mechanism and a
  margin well above that (B8).
