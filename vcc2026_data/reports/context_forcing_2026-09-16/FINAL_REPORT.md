# Context forcing and prediction without mean transfer — completed 2026-09-16

The models learned on the source data. The previous b18 score equalled b16 because
the submitted package deliberately reused b16 with the new context residual
disabled. It was not evidence that optimization produced zero weights. Historical
b16 itself contains a learned signed denoiser on top of empirical mean transfer.

Both requested follow-ups have now been implemented, trained and evaluated:
stronger residual learning through several mechanisms, and direct prediction with
no empirical mean-transfer input or output. Neither produces a stable improvement
over the appropriate baseline on the held-out lines in this experiment. Direct
prediction is substantially worse even when every method is restricted to the
same genes. This is evidence about the tested small models, features, training
objectives and data; it does not establish that transferable context cannot be
learned.

## Completed experiment

There are **158 distinct fits: 110 source-inner pilot fits and 48 full-data fits**.
The primary full-data comparisons use five whole-line holdouts and three paired
training/reconstruction seeds. Separately trained context controls use seed 0.
All **210 full count-score records**, 80 inner count-score records and 220 inner
native-score records are complete. The 50 full model records include two explicitly
labelled aliases of the context-free K562 direct model, not two extra fits.

The frozen [protocol](protocol.json) and [source-inner selections](selection.json)
precede all full outer evaluation. Outer lines are K562, RPE1, HepG2, Jurkat and
KOLF. K562/RPE1/HepG2/Jurkat are called RN4 below. The input table contains 103,295
genetic pseudobulk rows over 18,711 targets, reusing the existing Replogle/Nadig,
Jiang, Orion, CD4, Feng, H1 and KOLF preparation assets. Every study and stimulus
from the held-out biological line is excluded. Full fits use all remaining rows:

| Held-out line | Permitted training rows | Residual setting | Direct setting |
| --- | ---: | --- | --- |
| K562 | 93,145 | Standard MLP | Gene-only MLP |
| RPE1 | 101,173 | Strict context interaction + contrast weight 4 | Gene + context MLP |
| HepG2 | 101,285 | Standard MLP | Gene + context MLP |
| Jurkat | 101,094 | Standard MLP | Gene + context MLP |
| KOLF | 92,029 | Standard MLP | Gene + context MLP + contrast weight 4 |

This follow-up uses genetic data to isolate the context and mean-transfer
mechanisms. The completed genetic-versus-Tahoe comparison remains in the
[previous report](../simple_residual_2026-09-16/independent_verification_1937/FINAL_REPORT.md).
No new drug-transfer result is claimed here.

## What was changed and tested

All new networks are small: 64 gene-feature inputs, a 16-dimensional context
bottleneck, hidden width 64 and a rank-32 output basis. Residual bases are fitted
to permitted-source residuals; direct bases are fitted to permitted-source
perturbation effects. No held-out treated response enters these fits.

Seven residual candidates test the following mechanisms:

- Standard MLP with baseline coordinates as an input, output-energy penalty 0.5
  and AdamW weight decay 0.03.
- Relaxed penalties: no output-energy penalty and weight decay 0.001.
- Removal of baseline coordinates from the network input, keeping the additive
  empirical baseline outside the network.
- Baseline-input dropout of 0.5.
- A strict perturbation-by-context interaction with no target-only output head;
  its correction is exactly zero at zero centered context.
- That strict interaction with same-target, cross-context contrast supervision
  at weights 1 and 4.

Four direct candidates predict the perturbation effect from gene features alone,
gene features plus context, or gene plus context with contrast weights 1 and 4.
**Direct models never receive an empirical baseline as an input, target offset or
output addition.** Their final expression prediction is control plus the learned
effect. Effects outside the learned gene coordinates are zero, not filled with
mean transfer.

Residual strength is fixed at **alpha = 1** for this experiment. Alpha zero is
scored separately and cannot erase the residual arm. All candidates receive an
exhaustive coverage pass, then 800 pilot or 2,500 full biology-balanced batches
of 256 at learning rate 0.003. All full fit receipts confirm every permitted row
was seen. These are complete permitted populations, not the pilot target subset.

The pilot population contains 1,368 targets, selected independently of response
quality. Eleven candidates on ten excluded-line pairs produce 110 fits, each
evaluated in both query directions. Native decoded gene-space PDS shortlists two
candidates per route and outer line; installed count-reconstructed PDS across
all four inner lines selects one. Inner count checks use 64 targets and 64 cells
per target. Counts and the native gene-space metric are reported separately in
the [generated tables](generated_summary.md). Candidate averages over different
count-stage shortlists are not a uniform candidate ranking.

The standard MLP wins the residual inner selection for four of five outer lines;
strict interaction with the stronger contrast loss wins for RPE1. None of the
more forceful mechanisms establishes a broadly better transferable correction.

## Held-out count results

These are **local raw PDS values**, with larger values better and a no-change
value of 0.5. They are not the calibrated six-component official overall score.
Means average three seeds within each line and then give each line equal weight.

The direct comparison below gives every candidate the same response support:
5,592 common genes for RN4 and the 5,583 measured members of that set for KOLF.
All candidates use the same native count normalization and reconstruction.
This restricts where predicted effects can be nonzero; the scorer retains the
full native gene axis and count denominator rather than renormalizing a sliced
assay.

| Method, same gene support | RN4 mean | KOLF | All five lines |
| --- | ---: | ---: | ---: |
| No change | 0.500000 | 0.500000 | 0.500000 |
| Empirical mean transfer | 0.771064 | 0.558732 | 0.728597 |
| Historical b16 denoiser | 0.770710 | **0.597208** | 0.736009 |
| Source-refitted baseline, alpha zero | **0.840262** | 0.570120 | **0.786234** |
| Baseline + selected residual, alpha one | 0.833098 | 0.583731 | 0.783225 |
| Direct learned effect, no empirical mean | 0.588944 | 0.516477 | 0.574450 |

The direct model loses to the baseline with the same gene support in **15/15
paired line-seed comparisons**. Removing the empirical route therefore does not
improve prediction in this setup. It does not make the model incapable of
training: source-objective losses improve, as documented below.

For the practical full native-gene comparison, empirical methods retain their
source effects outside the learned subspace:

| Method, full native scoring axis | RN4 mean | KOLF | All five lines |
| --- | ---: | ---: | ---: |
| Empirical mean transfer | 0.781570 | 0.554683 | 0.736193 |
| Historical b16 denoiser | 0.794240 | **0.596134** | 0.754619 |
| Source-refitted baseline, alpha zero | **0.859225** | 0.565788 | **0.800538** |
| Baseline + selected residual, alpha one | 0.851158 | 0.568947 | 0.794715 |

The full-axis residual changes the five-line mean by **−0.005822**, and wins only
6/15 paired line-seed comparisons. Per-line details show why the pooled result
must not be described as a universal loss or a context success:

| Held-out line | Refit, alpha zero | Residual, alpha one | Difference | Residual seed wins |
| --- | ---: | ---: | ---: | ---: |
| K562 | 0.882070 | 0.883530 | +0.001460 | 1/3 |
| RPE1 | 0.810684 | 0.805006 | −0.005678 | 0/3 |
| HepG2 | 0.884560 | 0.856429 | −0.028131 | 0/3 |
| Jurkat | 0.859588 | 0.859666 | +0.000078 | 2/3 |
| KOLF | 0.565788 | 0.568947 | +0.003159 | 3/3 |

KOLF improves slightly relative to the refit in all three paired seeds, but the
result remains **0.027187 below b16**. No new model is promoted over b16.

Each prediction also receives 32 random target derangements after cell
reconstruction. The [generated tables](generated_summary.md) report
`(PDS − candidate null)/(refit PDS − refit null)`, paired within line and seed.
These per-line ratios weight low-baseline-excess lines differently from the raw
PDS mean. For example, the common-support residual's mean ratio is 1.017314 while
its raw mean PDS is lower than the refit. That ratio is not evidence of a stable
pooled improvement; the frozen primary criterion is raw PDS.

## Did the networks actually learn, and did context help?

Replaying seed-matched initial and fitted networks on every source training row
shows improved biology-balanced main MSE in **30/30 primary full fits**. Relative
to initialization, residual MSE reductions range from **11.73% to 42.20%**;
direct-model reductions range from **11.87% to 23.38%**. The fixed-pair data-loss
estimate also improves in all 30. This audit concerns training fit, not accuracy
on an unseen biology. It excludes decoupled weight decay from the reported
data-loss estimate. See [the exact source-objective audit](FULL_SOURCE_OBJECTIVE_AUDIT.json).

Residual outputs are substantial: the mean correction norm is **19.5–66.8%** of
the baseline norm on the common coordinates, depending on held-out line. They
did not simply collapse to zero.

Context attribution is less convincing. The following table uses **matched
seed-0 primary predictions and seed-0 controls**, never three-seed primaries
against one-seed controls. Residuals use full support; direct models use common
support, so compare controls within a column.

| Context treatment | Residual PDS, five lines | Direct PDS, five lines |
| --- | ---: | ---: |
| Correct query context | 0.795414 | 0.573623 |
| Separately trained without context | 0.793970 | 0.572024 |
| Separately trained with shuffled context assignments | 0.742342 | 0.560575 |
| Same fitted model, zero query context | **0.799576** | 0.559305 |
| Same fitted model, fixed swapped query context | 0.785168 | 0.569095 |

Correct context beats the separately trained no-context residual in only **2/5
lines**. Among the four direct folds that actually selected a context model, it
wins **2/4**. K562 selected gene-only direct prediction, so its context controls
are explicit identical aliases. RPE1's strict residual has zero correction by
construction in the no-context control; that control is not a learned gene-only
residual. Zero and swapped queries can also be outside the training distribution.

For KOLF specifically, all four seed-0 full-axis residual context controls score
above its correct-context residual. Its small gain over the refit consequently
does not establish that the destination's context information caused the gain.
For RPE1 direct prediction, a wrong query context scores above the correct one.
There is sensitivity to context and some useful source fitting, but no consistent
incremental context benefit across these controls.

The [prior-checkpoint diagnosis](PRIOR_DIAGNOSIS.md) independently found all 60
earlier full residual checkpoints nonzero. Only 7/20 earlier inner settings chose
alpha zero. It also identified an exact limitation of the earlier additive
linear architecture: target-centering cancels its query-context term. The new
interaction models avoid that cancellation. None of those facts changes why
b18 equalled b16: the submitted package itself was identical.

## Evaluation integrity and practical limits

Training anchors exclude the training example's entire destination biology.
Response bases, centering panels, context preprocessing and model selection use
permitted sources only. Compatible same-target contrasts stay within matched
assay/stimulus families; Jiang stimuli are never subtracted from one another.
Within each compatible family, source label centering uses the same training
target panel, preventing different target coverage from creating artificial
context contrasts. Response-poisoning, baseline-poisoning, source-coverage and
outside-subspace checks are recorded in the independent reviews.

Query expression uses only the frozen query-control partition, disjoint from
evaluation reference controls: 400/400 cells for RN lines and 400/11,600 for
KOLF. Nine common query genes absent from KOLF are imputed using the current
fold's source-control mean before context projection; they remain absent from
the KOLF scoring axis. All actual missing-coordinate and target-alias checks
are recorded. No held-out responses supply those missing input values.

Every full candidate uses 400 reconstructed cells per target, scale 3, the same
installed CP50K PDS scorer and the same paired reconstruction seeds. All 210
output-validity records pass, including finite, nonnegative reconstructed counts.
Negative signed effects are allowed; negative output counts are not. The exact
count accelerator was checked against the frozen reconstruction on actual data;
bulk counts and scores match. Its floating-point diagnostic for the total amount
clipped can differ at roundoff scale. All 82 receipts produced before the
acceleration change were preserved; 128 additional receipts bind the accelerated
worker and helper hashes. See [the rollout audit](FAST_ROLLOUT_AUDIT.json).

Important limits remain:

- These are previously studied biological panels, not untouched external
  replication. Source-only nested selection prevents this run's held-out
  responses from selecting its configurations; it cannot undo earlier project
  exposure to these panels.
- All scored targets have permitted source measurements. This is an unseen-cell
  context test, not a wholly unseen-target benchmark.
- Direct prediction is limited to a rank-32 basis on the common genes and the
  reused gene features. Its KOLF basis retains only **2.919%** of the centered
  observed effect energy, versus 22.8–45.4% on the RN lines. This is an
  evaluation-only projection diagnostic, not repeatable biological variance
  explained, and was never used for fitting or selection.
- This follow-up evaluates reconstructed PDS, not all six official score
  members. No new official score has been estimated or claimed. The prior
  completed study retains its full-six and repeat-based diagnostics.
- The small search can reject these particular attempts to remove or force the
  residual path; it cannot reject every direct or context model. Stronger
  context dependence by construction does not guarantee useful transfer.

## Artifacts, verification and submission status

The [deliverable index](DELIVERABLES.md) links code, frozen settings, checkpoints,
predictions, all scores and verification records. Full fold artifacts are under
[`full/folds/`](full/folds/): each trained model has `model.pt` and `receipt.json`,
each fold has source preparation and anchors, and `predictions/` contains the
saved full-axis effect arrays. [Combined scores](full/count_scores.csv) and
[diagnostics](full/diagnostics.csv) are the measured source for the tables.

Forty scientific tests passed: 29 data/model/runner/dispatch tests, eight scorer
tests and three accelerated-reconstruction tests. Independent audits checked
all 110 pilot fits, all 48 unique full checkpoints, all 210 score receipts, all
140 prediction arrays, numerical inference reproduction and report arithmetic.
The five frozen historical references, including the b16 checkpoint/package
and prior submission receipt/report/status, retain their original hashes.

No new official upload was made in this follow-up. The preserved published
result remains **b18, ID `Ax90uDSJhyMEFBacq5WB`, overall
0.12919304282260788**, from the package byte-identical to b16. That number is an
actual historical official result, not a score for any model in this report.
There is no pending fit, score retrieval or submission, and no replacement
deployment checkpoint was selected. The recommendation remains to keep the
empirical baseline and preserved b16 while treating the new context/direct
results as completed negative or mixed experiments.
