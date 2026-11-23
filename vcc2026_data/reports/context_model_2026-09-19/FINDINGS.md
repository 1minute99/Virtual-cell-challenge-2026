# Context-aware pseudobulk residual on the refit baseline — 2026-09-19

`prediction = refit_baseline + alpha * residual(context, perturbation)`

Local reconstructed-count PDS, leave-one-biology-out over five lines, seeds 0/1/2, frozen scoring
harness. No submission. Configs in `fits.json`, all scores in `scores.json`.

## Verdict against the two success criteria

1. **Beat the refit baseline: NO.** Best configuration 0.8008 vs 0.8005 (+0.0003), and performance
   falls monotonically as alpha rises. The optimum is alpha -> 0, i.e. the baseline itself.
2. **Genuine context dependence: PARTIAL, and confounded.** Substituting another line's context
   costs 0.042, which is large. But a gene-shuffled context SCORES BETTER than the correct one, so
   much of the ablation spread tracks residual magnitude rather than context accuracy.

Both criteria must hold. They do not.

## Setup

Architecture (24k-parameter class, deliberately small):
  perturbation branch  Linear(64 ESM + 32 baseline latent -> 64) + GELU + dropout 0.1
  context encoder      Linear(32 control coords -> 16) + GELU + dropout 0.1
  interaction          element-wise gate: p * W_g(c)   (explicit, not concatenation)
  head                 Linear(64 -> 32), ZERO-INITIALISED so training starts at the baseline
  output               rank-32 response space, decoded by the existing basis

Training: AdamW, lr 3e-3, weight decay 0.03, batch 256, max 20,000 steps, early stopping on a
HELD-OUT TRAINING BIOLOGY (patience 12 checks at 250-step intervals). Every fold-specific quantity
-- refit baseline, anchors, response basis, control basis -- refitted on training biologies only.

The enabling data fix: the Replogle-Nadig effects were computed against per-target batch-weighted
controls, but the assembler paired them with the line's POOLED control (median L2 error 0.646 over
9,624 genes). Restoring the matched control raised distinct contexts from 205 to 16,263, so context
is no longer a synonym for source. Without this the experiment was structurally impossible.

## Results

| configuration | K562 | RPE1 | HepG2 | Jurkat | KOLF | mean |
|---|---|---|---|---|---|---|
| refit baseline | 0.8821 | 0.8107 | 0.8846 | 0.8596 | 0.5658 | **0.8005** |
| mean transfer | 0.7604 | 0.8073 | 0.7472 | 0.8114 | 0.5547 | **0.7362** |
| ctx mean a=0.25 | 0.8861 | 0.8019 | 0.8835 | 0.8665 | 0.5661 | **0.8008** |
| ctx mean a=0.5 | 0.8878 | 0.7926 | 0.8802 | 0.8686 | 0.5675 | **0.7993** |
| ctx mean a=1.0 | 0.8782 | 0.7729 | 0.8708 | 0.8554 | 0.5727 | **0.7900** |
| ctx mean a=2.0 | 0.7612 | 0.6891 | 0.8516 | 0.7642 | 0.5845 | **0.7301** |
| ctx pop a=0.25 | 0.8841 | 0.8065 | 0.8846 | 0.8583 | 0.5652 | **0.7997** |
| ctx pop a=0.5 | 0.8859 | 0.8021 | 0.8833 | 0.8551 | 0.5662 | **0.7985** |
| ctx pop a=1.0 | 0.8868 | 0.7863 | 0.8798 | 0.8462 | 0.5743 | **0.7947** |
| ctx pop a=2.0 | 0.8804 | 0.7337 | 0.8695 | 0.8241 | 0.5808 | **0.7777** |

| ablation (alpha 1.0) | mean encoder | population encoder |
|---|---|---|
| correct | 0.7900 | 0.7947 |
| wrong | 0.7484 | 0.7642 |
| shuffled | 0.7954 | 0.7971 |
| zero | 0.7731 | 0.7915 |

## The training diagnostic that decided it

Because the head is zero-initialised, the validation loss at step 0 IS the refit baseline's loss.
Logging it turned out to be the single most informative number: in **9 of 10 fits the model never
beat "predict no residual at all"** on the held-out validation biology, with deterioration from 1%
(K562) to 46% (HepG2). Training loss declined in every fit, so without this reference the curves
would have looked like healthy convergence.

The one exception was KOLF/mean (1.5360 vs 1.5470, +0.7%). Its sibling KOLF/pop missed (1.5509),
which argues noise around a threshold rather than a real effect.

## Mean vs population context

Inconsistent and within noise: population won on RPE1, HepG2 and Jurkat folds, lost on K562 and
KOLF. Neither cleared the zero-residual reference anywhere except KOLF/mean. **Important caveat:**
only 16,062 of 92,029 rows (17%) carry a genuine population; the rest repeat a single profile
because those sources really do have one control each. So this is not a clean head-to-head.

## The one result that cuts the other way

On KOLF the residual HELPS and helps MORE with larger alpha: 0.5658 (baseline) -> 0.5675 -> 0.5727
-> 0.5845 at alpha 0.25/0.5/1.0/2.0, i.e. **+0.019**. The other four lines all degrade with alpha.
KOLF is the line least like the training sources and the closest analogue to the scored contexts.
This is one line and should not be over-read, but it is monotone and it is the opposite of the
aggregate.

## Failure analysis

The residual is harmful on source-like lines and helpful on the most distant one. The ablation
pattern (shuffled > correct > zero > wrong) is largely ordered by how LARGE a residual each variant
produces, not by how accurate it is: at alpha 1 the residual subtracts value, so anything that
shrinks it scores better. That is the signature of a correction with the wrong direction rather
than the wrong scale.

## Implementation defects found (all fixed; they bear on confidence)

1. Inner validation split matched zero rows (biology names canonicalise to upper case). Caught only
   because the split size was asserted -- would otherwise have silently disabled early stopping.
2. The population encoder was never TRAINED on sets, only evaluated on them; the two arms returned
   bit-identical validation losses, which is how it surfaced.
3. float32 sampling weights rejected by numpy ("probabilities do not sum to 1").
4. Set size read as the control feature width for the population encoder.
Four defects, all in the population encoder path, three of which failed loudly and one silently.
