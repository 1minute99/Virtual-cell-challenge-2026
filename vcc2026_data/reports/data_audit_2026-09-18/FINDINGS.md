# Data audit — 2026-09-18

All numbers are local reconstructed-count PDS or raw `cell-eval2` metrics on previously studied
held-out panels. Not official scores. No submission was made.

## 1. Two claims from the first pass, corrected

**RETRACTED: "we train in CP10K and are scored in CP50K, and that is why expression accuracy is
zero."** `harness_pred.reconstructed_blocks` applies the delta as
`((CP10K + 1) * exp(delta) - 1) * L / 1e4`, the exact inverse of log1p(CP10K). A perfect CP10K
delta therefore reproduces the cells and scores perfectly. Confirmed by the oracle run below.
What survives is only that the metric's space weights low-abundance genes ~3x more than our
training loss does. That is an objective-weighting mismatch, not a bug.

**RETRACTED: "dropping the x3 scale would recover expression accuracy."** Measured on KOLF with
our own delta: at scale 1 the expression error is still 2.63 where 1.00 already means "predicted
no change", so it scores zero at either scale, while scale 1 is worse on the other five metrics.
x3 is the right choice for the delta we have.

## 2. The pipeline is sound; the delta is what is weak

Oracle test: real measured KOLF effects fed through our own reconstruction and scorer.

| metric | no-skill | replicate | oracle x1 | oracle x3 | our delta x1 | our b19 x3 |
|---|---|---|---|---|---|---|
| pds | 0.500 | 0.9515 | **1.0000** | 1.0000 | 0.5582 | 0.5659 |
| expression (lower better) | 1.000 | ~-0.91 | **-2.853** | 3.169 | 2.633 | 5.154 |
| log-FC (lower better) | 1.000 | 0.392 | **0.222** | 1.352 | 0.971 | 0.963 |
| fidelity | 0.513 | 0.807 | **0.985** | 0.962 | 0.252 | 0.499 |
| reach | 0.079 | 0.965 | **1.000** | 1.000 | 0.085 | 0.101 |
| jaccard | 0.030 | 0.401 | 0.065 | 0.012 | 0.010 | 0.012 |

- No structural defect in reconstruction, gene mapping, target alignment or normalization.
- Scaling a GOOD delta by 3 is destructive: expression -2.85 -> +3.17, log-FC 0.222 -> 1.352.
- One real ceiling: even a perfect mean delta scores 0.065 jaccard against 0.401 for a replicate,
  because we shift real control cells rather than generating perturbed heterogeneity. That is the
  only place single-cell modeling would pay, and it is worth at most 1/6 of the score.

## 3. THE data problem: evaluation panels were cutting the training gene axis

`simple_residual_assemble.common_axis` intersects EVERY source axis:

| axis | official genes |
|---|---|
| xl TRAINING sources jointly | **17,847** |
| rn_K562 / rn_RPE1 / rn_Jurkat / rn_HepG2 evaluation panels | 8,248 / 8,748 / 8,881 / 9,624 |
| plus Jiang panels and the Tahoe drug cache | — |
| **resulting training axis** | **5,592** |

The genetic model never uses Tahoe's 12,922 drug rows, yet Tahoe constrained it. The official
contexts measure all 18,533 features, so ~70% of the scored readout got no learned component.

**Mean transfer was never affected.** `simple_residual_baselines.empirical` reads the per-source
effect stores on the full official axis with proper masking. Measured identical in both pipelines
(0.7816 on four lines, difference 0.0000). Only the refit ridge and the neural residual were
restricted.

## 4. Biologies matter far more than genes

Narrow arm = identical rows to the new table, restricted back to 5,592 genes, i.e. the control
that isolates the cost of dropping the RN/Jiang rows (9% of rows, three extra biologies).

| line | refit frozen | refit narrow | b19 frozen | b19 narrow |
|---|---|---|---|---|
| K562 | 0.8821 | 0.7977 | 0.8835 | 0.7767 |
| RPE1 | 0.8107 | 0.7465 | 0.8103 | 0.7404 |
| HepG2 | 0.8846 | 0.7888 | 0.8564 | 0.7777 |
| Jurkat | 0.8596 | 0.7890 | 0.8597 | 0.7756 |
| KOLF | 0.5658 | **0.5683** | 0.5689 | **0.5712** |
| five-line mean | **0.8005** | 0.7381 | 0.7958 | 0.7283 |

Cost of dropping those biologies: -0.062 on the baseline, -0.068 on the full model.

**KOLF is the exception and the warning.** It is the only line where losing them costs nothing.
KOLF is iPSC, furthest from every source, and the closest local analogue to the unseen official
contexts. Extra biologies help on lines that resemble our sources and do not help on the line that
does not. The five-line mean is therefore flattering, and matches the previously measured
cross-biology transfer ceiling.

**The loss falls on the BASELINE, not the neural residual** (d(refit) ≈ d(b19) per line), so the
baseline is what must see the masked rows. Fitting it on dense sources only would reproduce the
entire loss while looking like a fix.

## 5. Artifacts built

- `table_fullaxis` — 84,752 xl rows x **17,282** genes. Verified against the frozen table: same
  rows, targets, bios, sources, weights, and **bit-identical values on all 5,592 shared genes**,
  old axis a strict subset, fully dense. See `table_fullaxis_verification.json`.
- `table_fullaxis_masked` — **103,295** rows (exactly the frozen genetic row count) x 17,282
  genes, every biology restored, 9.1% of entries NaN. Jiang covers ~16,200 genes; only the four RN
  panels are genuinely narrow.
- `fullaxis_masked_data.py` — mask-aware layer. Baseline path imputes unmeasured entries at the
  target mean (neutral for the pooling that follows); training path never imputes and fits each
  row's rank-32 latent by least squares over its measured genes only, sharing one pseudo-inverse
  per distinct mask. Basis fitted on fully measured rows.
- `score_delta_full.py` — scores any delta at any scale through the frozen six-metric path.

## 6. Open

The masked baseline has NOT produced a number yet; it was OOM-killed once and stopped once, and
now runs from `fullaxis_sequence.sh` after sweep 2, serialized and memory-gated. The target to
beat is the frozen five-line refit mean of **0.8005**, which remains the best local predictor,
ahead of b19 at 0.7958.
