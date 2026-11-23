# Stronger residual sweep and official b19 — complete

**b19 is published at 0.13646991123621474, rank 193 at retrieval.**
This is +0.007276868 (5.63%) above preserved b16's overall score.
The new package contains a nonzero context-conditioned residual and is
not a copy of b16/b18. b19 is the best recorded official overall result in this
project. b16 remains preserved as the reference checkpoint and has better official
PDS in this comparison.

The higher overall result is mainly an LFC-NMAE improvement. **Official PDS falls**,
and the controlled local alpha-zero baseline remains better on average. This is
a whole-system official improvement, not evidence that stronger residuals or
context information consistently improve perturbation discrimination. b19 also
changes the source-refitted baseline rank and training population relative to
historical b16, so its official score cannot isolate the causal residual effect.

## Actual official result

- Entry: **`pamfw9b1glTBvlXvEyQA`**, name **b19**, status **published**.
- Uploaded 3,515,402,240 bytes with server checksum verification.
- Published score retrieved at `2026-09-17T03:17:29.791724+00:00`.
- Validation panel `vcc2026-val-1`, anchor `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`;
  both match b16's recorded result.
- Rank 193 is the observed rank at retrieval, not a guaranteed current rank.

All values below are returned official calibrated score components, not local
estimates. Overall is exactly their arithmetic mean.

| Metric | Preserved b16 | New b19 | Difference |
| --- | ---: | ---: | ---: |
| Overall | 0.129193043 | 0.136469911 | +0.007276868 |
| PDS | 0.554670505 | 0.543156464 | -0.011514041 |
| MSE | 0.000000000 | 0.000000000 | +0.000000000 |
| LFC NMAE | 0.030130366 | 0.075023684 | +0.044893318 |
| Direction fidelity | 0.010708046 | 0.018796652 | +0.008088607 |
| Direction reach | 0.187439871 | 0.188644931 | +0.001205060 |
| Significant-gene Jaccard | -0.007790531 | -0.006802264 | +0.000988267 |

For completeness, the server's underlying raw metrics are:

| Raw metric | b16 | b19 |
| --- | ---: | ---: |
| pds_cosine | 0.751222594 | 0.746086957 |
| expr_mse_unbiased_capped_norm | 17.250293613 | 15.516872987 |
| de_wilcoxon_lfc_nmae | 0.981584189 | 0.955234555 |
| de_wilcoxon_direction_fidelity_yield_raw | 0.515907900 | 0.518357642 |
| de_wilcoxon_direction_reach_raw | 0.245284204 | 0.246322983 |
| de_wilcoxon_sig_jaccard | 0.027528719 | 0.027894231 |

Receipts: [actual official response](submission/receipt.json),
[new score ledger](submission/official_scores_new.csv),
[exact comparison data](OFFICIAL_COMPARISON.json). The historical submission
ledger and b16/b18 records were retained; the new result is recorded separately.
Only **one** new entry was uploaded. No hyperparameters, final fit or package were
changed after seeing the official score.

## Stronger residual sweep

The fixed-model sweep used alpha **1, 1.5, 2, 3 and 4**, with zero as a benchmark.
It completed 120 source-inner count checks, 90 three-seed outer outcomes and
20 matched context controls. Only the added residual was multiplied; baseline,
subspace, learned checkpoint and scale-3 reconstruction stayed fixed. This reuses
the preceding trained models; it is not described as 90 new training fits.

All five source-only whole-line selections chose **alpha 1**. The standard MLP
was selected for four lines; RPE1 selected the strict interaction plus contrast
weight 4. The frozen frequency rule chose standard MLP / alpha 1 for full training
and the explicitly requested nonzero-residual official evaluation. Alpha zero
remained an evaluated benchmark and was ineligible for that submission.

Local raw PDS, three paired seeds averaged within each line:

| Alpha | Five-line mean PDS | KOLF PDS | Valid outcomes |
| --- | ---: | ---: | ---: |
| 0 | 0.800538 | 0.565788 | 15/15 |
| 1 | 0.794715 | 0.568947 | 15/15 |
| 1.5 | 0.790091 | 0.571465 | 15/15 |
| 2 | 0.785289 | 0.574495 | 14/15 |
| 3 | — | 0.582167 | 12/15 |
| 4 | — | 0.584307 | 11/15 |

The alpha-2 mean includes an ineligible count-cap diagnostic. Alpha 3/4 have
incomplete five-line means because of unscorable HepG2 outputs; they are not
imputed. Five stronger outputs are scored but exceed the count cap, and three
are explicitly unscorable. The selected alpha-1 outputs and all context controls
are valid. KOLF benefits from larger corrections, but the other lines generally
worsen, and KOLF alpha 4 still trails b16's local 0.596134.

For practical full-axis comparisons, local five-line PDS is 0.736193 for mean
transfer, 0.754619 for historical b16, 0.800538 for the source-refitted alpha-zero
baseline, and 0.794715 for the selected alpha-one residual. These local raw PDS
numbers are not comparable in scale to the official calibrated overall score.

Correct context beats the separately trained no-context residual in only 2/5
lines; the zero-query control is better on average. This does not establish a
robust incremental context benefit. See [all per-line results and controls](LOCAL_RESULTS.md),
[paired null comparisons](paired_null_comparisons.csv), and the prior
[context-forcing/direct-model report](../context_forcing_2026-09-16/FINAL_REPORT.md).
The latter's direct model without empirical mean transfer remains worse locally,
including a comparison on identical response genes; it was not submitted here.

## Selected full-data model

The new checkpoint is a small **24,192-parameter MLP**, context dimension 16,
hidden width 64, 64 gene features and a **rank-32 source-residual basis**. It adds
an alpha-1 correction to a newly source-refitted rank-32 empirical denoiser.
Baseline coordinates are an input; energy penalty is 0.5, weight decay 0.03,
learning rate 0.003 and batch size 256. Seed 0 received one exhaustive coverage
pass (404 batches) and 2,500 biology-balanced sampled batches: 103,295 unique
permitted genetic rows, 18,711 targets, 743,295 examples seen in total.

Data reuse Replogle/Nadig, Jiang, Orion, CD4, Feng, H1 and KOLF assets.
**No drug rows** enter this fit. The preceding genetic-only/joint-Tahoe/drug-
pretraining comparison found no stable drug benefit and remains in its
[completed report](../simple_residual_2026-09-16/independent_verification_1937/FINAL_REPORT.md).
This scalar-strength continuation adds no new claim about drug transfer.

All 300 official targets have source measurements. Query information is only
18,400 supplied controls per official A/B/C context. Training anchors exclude the
example's destination biology; the residual basis, centering and preprocessing
use sources only. Output outside the residual subspace retains the baseline.
A/B/C corrections are nonzero and distinct, with norms 5.6713/6.3964/5.9549 against
baseline norm 40.9288. The fixed training-probe MSE decreases from 0.576515 to
0.459220; this is training fit, not unseen-context validation. The model did not
simply produce a zero residual.

Historical b16 is the preserved rank-128 checkpoint, not this newly trained
expanded-data model. Its model and package remain unchanged.

## Artifacts and verification

- [New checkpoint](final/model.pt), [basis/preprocessing](final/coordinates.npz),
  [baseline operator](final/source_baseline.npz), [training receipt](final/training_receipt.json).
- [Frozen fit manifest and artifact hashes](final/COMPLETE.json),
  [source-inner decision](final_choice.json), [pre-upload freeze](SUBMISSION_FREEZE.json).
- [Submitted b19 package](submission/b19.vcc), [cell matrix](submission/b19.h5ad).
  Package SHA256: `ca8d60b1cd8f5f39b215cfde2cd4973f2fc2047967463384760ff0c125f92f47`.
- [Official preparation/validation and published receipt](submission/receipt.json).
  The matrix has 360,000 cells, 18,533 genes and 400 cells per target/context.
- [Eight targeted scientific tests](final_scientific_tests.log),
  [independent full-fit review](STRENGTH_FINAL_ADAPTER_REVIEW.json),
  [independent complete-grid review](STRENGTH_OUTER_FINAL_REVIEW.json),
  [physical package review](PACKAGE_REVIEW.json),
  [local-report review](LOCAL_REPORT_REVIEW.json).
- [Deliverable index](DELIVERABLES.md) and [restart status](research_status.json).

The official validator passed, reconstruction matched independent cell checks,
and live authorization/allowance checks passed immediately before upload.
The allowance endpoint reported `limit_reached=false`; it did not disclose a
numeric remaining quota. Current requirements were checked against the
[official VCC guide](https://vcc-cli-wiki.virtualcellchallenge.org/) and
[CLI reference](https://vcc-cli-wiki.virtualcellchallenge.org/cli-reference).

The local test still cannot exactly reproduce hidden official contexts, target
panels, gene coverage or control availability. Previously studied outer panels
are not untouched external replication. All targets here have source effects;
this is not an unseen-target result. Local count-reconstructed PDS and the actual
official six-metric score are reported separately. Invalid outputs, startup/CSV
corrections, original logs and hash-preservation evidence remain archived as
linked in the deliverable index. No training, upload or score poll remains pending.
