# Pooled-control PCA context on top of the refit baseline: harness verdict

Date: 2026-09-20. Scripts: `ctx_harness_eval.py`, `ctx_harness_report.py`.
Baseline to beat: frozen five-line refit mean **0.8005**. Seeds 0/1/2, scale 3, 400 cells.

## Answer: no. It does not beat the refit baseline.

| arm | K562 | RPE1 | HepG2 | Jurkat | KOLF | five-line | vs refit |
|---|---|---|---|---|---|---|---|
| refit (alpha=0, validation) | 0.8821 | 0.8107 | 0.8846 | 0.8596 | 0.5658 | **0.8005** | +0.0000 |
| context both, alpha=0.5 | 0.8854 | 0.8343 | 0.8528 | 0.8628 | 0.5681 | 0.8007 | +0.0001 |
| context both, alpha=1 | 0.8878 | 0.8522 | 0.8119 | 0.8661 | 0.5702 | 0.7976 | -0.0029 |
| context both, alpha=2 | 0.8906 | 0.8759 | 0.7302 | 0.8707 | 0.5683 | 0.7871 | -0.0134 |
| context both, alpha=3 | 0.8908 | 0.8828 | 0.6745 | 0.8729 | 0.5602 | 0.7762 | -0.0243 |
| context both, alpha=4 | 0.8880 | 0.8736 | 0.6395 | 0.8742 | 0.5497 | 0.7650 | -0.0356 |
| interaction (`s`) only, alpha=1 | 0.8832 | 0.8330 | 0.8677 | 0.8640 | 0.5578 | 0.8011 | +0.0006 |
| generic shift (`g`) only, alpha=1 | 0.8873 | 0.8350 | 0.8216 | 0.8615 | 0.5750 | 0.7961 | -0.0045 |
| **wrong context, alpha=1** | 0.8751 | 0.8059 | 0.8927 | 0.8599 | 0.5708 | **0.8009** | **+0.0003** |
| zero context, alpha=1 | 0.8822 | 0.8104 | 0.8847 | 0.8597 | 0.5658 | 0.8006 | +0.0000 |

The best arm is +0.0006, inside the ±0.010 seed floor. The untuned natural choice, alpha=1 with
both terms, is -0.0029. **Wrong context scores 0.8009, statistically tied with correct context and
better than the full correct-context model.** That single row is the result.

Validation: the alpha=0 arm reproduces the stored refit run to 1e-5 per seed (0.881962 vs
0.881951 on K562 seed 0), so the baseline path is faithful and the comparison is like-for-like.

## Paired per-seed differences

Seeds are shared across arms, so pairing removes most of the harness seed noise.

| arm | mean | sd | min | max | wins |
|---|---|---|---|---|---|
| both, alpha=0.5 | +0.0001 | 0.0185 | -0.0321 | +0.0238 | 12/15 |
| both, alpha=1 | -0.0029 | 0.0389 | -0.0730 | +0.0416 | 12/15 |
| both, alpha=3 | -0.0243 | 0.1000 | -0.2107 | +0.0724 | 9/15 |
| interaction only | +0.0006 | 0.0136 | -0.0171 | +0.0224 | 9/15 |
| generic shift only | -0.0045 | 0.0313 | -0.0633 | +0.0246 | 12/15 |
| wrong context | +0.0003 | 0.0059 | -0.0072 | +0.0083 | 8/15 |
| zero context | +0.0000 | 0.0002 | -0.0004 | +0.0002 | 10/15 |

Note the shape: 12 of 15 wins at alpha=1 with a *negative* mean. Four lines gain a little and one
loses a lot. Counting folds won would have reported this as a success.

## HepG2 is the failure, and it is not a small one

| line | refit | ctx alpha=1 | ctx alpha=3 | deviation/refit norm |
|---|---|---|---|---|
| K562 | 0.8821 | +0.0058 | +0.0087 | 0.140 |
| RPE1 | 0.8107 | +0.0415 | +0.0721 | 0.347 |
| Jurkat | 0.8596 | +0.0065 | +0.0134 | 0.119 |
| KOLF | 0.5658 | +0.0044 | -0.0056 | 0.290 |
| **HepG2** | 0.8846 | **-0.0727** | **-0.2101** | **0.487** |

HepG2 gets the largest predicted correction of any line and it is confidently wrong. On HepG2
wrong context *beats* correct context (+0.0081), the reverse of every other line. Raising alpha
amplifies the error faster than it amplifies the other four lines' gains, which is why the
five-line mean falls monotonically while four of five lines individually improve.

## Honest alpha selection makes it worse, not better

Reading the best alpha off the five-line mean is selection on the evaluation set. Choosing alpha
for each line from the *other four* lines only:

| line | alpha* | PDS | refit | gain |
|---|---|---|---|---|
| K562 | 0.5 | 0.8854 | 0.8821 | +0.0033 |
| RPE1 | 0.5 | 0.8343 | 0.8107 | +0.0237 |
| HepG2 | 3 | 0.6745 | 0.8846 | **-0.2101** |
| Jurkat | 0.5 | 0.8628 | 0.8596 | +0.0032 |
| KOLF | 0.5 | 0.5681 | 0.5658 | +0.0023 |
| **nested five-line** | | **0.7650** | 0.8005 | **-0.0355** |

Nested selection is worse than any fixed alpha because HepG2's optimum is the opposite of every
other line's, so the other four lines actively mis-select for it. A method whose hyperparameter
cannot be chosen out-of-sample is not deployable.

## Why the cosine benchmark said the opposite

`reports/context_benchmark_2026-09-20` measured cosine against a batch-held-out replicate with the
cross-biology mean as the baseline, and reported correct context recovering 22.3% of oracle
headroom, winning 3/3 folds, beating wrong context by +0.0399. Every one of those conclusions
fails to transfer:

| | cosine benchmark | harness PDS |
|---|---|---|
| HepG2, correct context | +0.0294 (42% of headroom) | -0.0727 |
| correct vs wrong | +0.0399 for correct | +0.0003 for **wrong** |
| aggregate | +0.0268 | -0.0029 |

Same line, same context representation, same fold discipline, opposite sign. This is exactly the
failure mode recorded in the standing note that row cosine is not PDS: a generic, well-aimed
all-target response inflates row cosine while doing nothing for, or actively harming,
discrimination. I did not expect it to bite here because the ablations looked clean, and the
ablations looked clean *in the wrong metric*.

## What this overturns

The recommendation written six hours earlier in `reports/context_benchmark_2026-09-20/FINDINGS.md`
said "continue the context-conditioned model" on the strength of the cosine result, while flagging
score conversion as unverified. It is now verified and it does not convert. That recommendation is
withdrawn.

Standing record, unchanged: **the refit baseline at 0.8005 remains the best local predictor**, and
it has still never been scored officially.

## What is salvageable

1. **The interaction term, not the generic shift.** `s_only` is +0.0006 and `g_only` is -0.0045.
   The generic shift is the dangerous half: one vector applied to every target, which wrecks
   discrimination across the whole panel when the predicted direction is wrong. Any future attempt
   should drop `g` entirely. But `s_only` at +0.0006 is still indistinguishable from zero.
2. **Four of five lines do gain**, modestly and consistently across seeds. The blocker is one
   catastrophic fold, not an absence of signal everywhere.
3. **A per-line trust gate is the only obvious route.** Something that detects "the predicted
   correction for this line is large and unsupported" and falls back to refit would keep the four
   gains and drop HepG2's loss. That is a real research direction, but it needs a criterion
   validated out-of-sample, and with nine biologies there is very little to validate on.

## Recommendation

Do not pursue the context residual further on the current data without first solving fold-level
robustness. The evidence across three experiments now points the same way: the honest
cross-biology signal is small, the sample is eight training context vectors, and every gain
measured in a surrogate metric has failed to survive contact with PDS. Prioritise adding
independent biologies over further architecture work, and score the refit baseline officially.
