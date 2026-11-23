# DE-set classifier — how many genes move, and which (2026-09-24)

Question (user): leaving everything else aside, can a classifier predict the number of significantly moving genes
per knockdown and which genes they are — or, with the count fixed, only which genes?

## Design
- Labels: real cell-eval2 Wilcoxon DE tables (pdex, per-pert BH FDR < 0.05, 5 CPM gate; on-target gene excluded,
  as in `de_sig_jaccard`). Six contexts: H1 (VCC2025, Flex, 300 cells/target), KOLF (200 cells/target, 279 official
  targets), RN K562/HepG2/Jurkat/RPE1 (the same 300 essential-gene targets in four lines).
- Leave-one-context-out. The held-out context's own source is removed from every row's features; each training
  context's own source is removed from its own rows; the label prior excludes the row's own context.
- Features (all available for an unseen official context): per-source Poisson z / log fold change of g under t,
  aggregated over 12-13 public screens (`dec_sources.py`: RN x4, Orion HCT116/HEK293T, CD4 x3, KOLF, Feng x3, H1);
  CRISPRi TSS distance; the context's control expression of g and of t; source-side gene responsiveness; label prior;
  source knockdown strength, on-target KD, DepMap mean effect.
- Models: sklearn HistGradientBoosting (full; generic = no pair features) and single-feature rankings; the full model
  on target-shuffled test features isolates the target-specific part. Count: HGB regressor on target-level features,
  log1p(count) centred within context (the context level is set separately).
- Code: `scripts/dec_sources.py`, `scripts/dec_eval.py` (`python -P scripts/dec_eval.py [CTX]` re-runs one fold and
  saves scores). Results: `results.json`, `scores_H1.npz`, logs in `logs/`. Runtime 6.5 min for all folds.

## 1. Which genes, count given (R-precision = hits in top-n / n, n = true count; mean over targets with n >= 1)

| held-out | base rate | full | generic (no pair) | full, targets shuffled | expression only |
|---|---|---|---|---|---|
| H1 | 0.083 | **0.207** | 0.187 | 0.186 | 0.173 |
| KOLF | 0.005 | 0.128 | 0.194 | 0.121 | **0.311** |
| RN K562 | 0.016 | **0.222** | 0.156 | 0.131 | 0.169 |
| RN HepG2 | 0.050 | **0.266** | 0.195 | 0.178 | 0.189 |
| RN Jurkat | 0.023 | **0.242** | 0.152 | 0.124 | 0.153 |
| RN RPE1 | 0.076 | **0.335** | 0.285 | 0.257 | 0.277 |

Per-target AUROC of the full model 0.75 (H1) to 0.85; precision of the top 10 (normalised by min(10, n)):
H1 0.39, KOLF 0.19, RN 0.41-0.59. By target strength (H1): weak targets (<= 10 DE, 76/300) 0.05, 11-200 DE 0.10,
> 200 DE 0.35.

- **Most of what the model gets right is target-generic** on H1: shuffling the targets keeps 0.186 of 0.207. The
  target-specific part is clear only when the same targets were measured in same-lab sibling lines (RN: shuffled
  0.12-0.26 -> full 0.22-0.34; top-10 precision roughly doubles).
- **Weak targets are unpredictable** (H1 R-precision 0.05; they are also where the scorer's empty-set convention
  bites).
- **KOLF reverses**: expression alone (0.311) beats every model. With 200 cells and low depth, significance is mostly
  power, so the called genes are simply the most highly expressed; models trained on higher-power contexts (H1 median
  145 DE, RPE1 355) rank mid-expression genes too high. Which genes pass the test depends on cells x depth as much as
  on biology, so a classifier has to be calibrated on a context with the target platform's power (for the official
  panel: H1 Flex).

## 2. How many genes (per-target count)

| held-out | true median | Spearman count model | Spearman sum of probabilities | log-MAE model (oracle level) | log-MAE constant |
|---|---|---|---|---|---|
| H1 | 145 | 0.19 | 0.39 | 2.23 | 2.27 |
| KOLF | 28 | 0.11 | 0.16 | 1.62 | 1.04 |
| RN K562 | 21 | 0.51 | 0.44 | 1.56 | 1.88 |
| RN HepG2 | 11.5 | 0.62 | 0.64 | 1.97 | 2.52 |
| RN Jurkat | 22.5 | 0.56 | 0.55 | 1.61 | 2.02 |
| RN RPE1 | 355 | 0.66 | 0.64 | 1.90 | 2.00 |

The relative count is partly predictable only with same-lab sibling lines (Spearman 0.5-0.66); on H1 and KOLF it is
weak (0.16-0.39) and the log error is no better than a constant. The absolute level (median 11-355) is set by power,
not by the model.

## 3. Set agreement (the scorer's `sig_jaccard`, mean over all 300 H1 targets)

| count strategy (ranking = full model) | Jaccard |
|---|---|
| oracle k (true count per target) | **0.151** |
| oracle k, expression ranking | 0.127 |
| count model + oracle context level | 0.046 |
| count model + training-context level (fully automatic) | 0.019 |
| sum of probabilities (automatic) | 0.035 |
| constant k = 145 / 1,000 / 3,000 / 5,000 | 0.040 / 0.079 / 0.098 / 0.100 |
| **S6 recipe on the H1 harness (cell-eval2)** | **0.066** |

With the true counts the classifier would more than double the current Jaccard; with any count we can actually
predict it is worse than S6. Calling ~3,000 genes for every target scores 0.098 because the mean is carried by the 84
targets with > 1,000 DE genes (Jaccard 0.26 there, ~0 on weak targets). That is a Jaccard-only observation: fidelity
(k / max(n_pred, n_conf)) and reach need correct direction and punish over-calling, so it was not tested as a
recipe.

## Conclusion
- "Which genes" is moderately learnable: with the count given, ~21 % (H1) to ~34 % (RPE1) of the top-n are right
  (base rates 0.5-8 %), but on the official-like H1 most of it is generic (expression and responsive genes), not
  target-specific. Target-specific gene choice needs the same target measured in closely related lines.
- "How many" is the bottleneck: counts span 0-6,929 per target, the relative count is weakly predictable (H1 0.39) and
  the absolute level is a property of the assay's power. Knowing the true counts is worth +0.085 Jaccard on H1;
  predicting them loses to the current emission.
- Direction (needed by fidelity and reach) is outside this classifier.
- Not uploaded, nothing changed in the submission pipeline.
