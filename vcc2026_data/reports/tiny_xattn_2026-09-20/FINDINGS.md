# Perturbation-conditioned cross-attention for unseen-context prediction: large and tiny

Date: 2026-09-20. Scripts: `xattn_context_model.py`, `xattn_harness_eval.py` (large),
`tiny_xattn_eval.py` (tiny). Results: `xattn_context_2026-09-20/`, `tiny_xattn_2026-09-20/`.
Baseline: frozen five-line refit mean **0.8005**, seed-noise floor ±0.010.

## Verdict

**Stop this direction.** Correct context beats a real wrong biology's context on **1 of 5** lines,
and that one line is KOLF, whose deviation is 99.3% noise. Reducing capacity 28-fold fixed the
overfitting symptom completely and did not change the sign of the result.

## The comparison

| arm | params | five-line | vs refit | paired | sd | wins |
|---|---|---|---|---|---|---|
| **PCA interaction (ridge)** | — | **0.8030** | **+0.0024** | +0.0024 | 0.0107 | 9/15 |
| STATE mean-pool (ridge) | — | 0.8012 | +0.0007 | +0.0007 | 0.0019 | 9/15 |
| STATE mean-pool (neural) | 277,760 | 0.8012 | +0.0007 | +0.0007 | 0.0006 | 13/15 |
| **tiny wrong context (real biology)** | 9,728 | 0.8010 | +0.0005 | +0.0005 | 0.0009 | 12/15 |
| refit baseline | 0 | 0.8005 | +0.0000 | — | — | — |
| tiny zero context | 9,728 | 0.8005 | +0.0000 | +0.0000 | 0.0000 | 0/15 |
| tiny correct, alpha 0.25 | 9,728 | 0.8002 | -0.0003 | -0.0003 | 0.0003 | 2/15 |
| large wrong context | 277,760 | 0.8002 | -0.0003 | -0.0003 | 0.0006 | 5/15 |
| **TINY cross-attention, correct** | **9,728** | **0.7999** | **-0.0007** | -0.0007 | 0.0012 | 4/15 |
| **LARGE cross-attention, correct** | **277,760** | **0.7992** | **-0.0013** | -0.0013 | 0.0014 | 5/15 |
| tiny32 (32 tokens, rank 32) | 10,240 | 0.7991 | -0.0014 | -0.0014 | 0.0008 | 0/15 |
| tiny corrupted context | 9,728 | 0.7988 | -0.0018 | -0.0018 | 0.0017 | 3/15 |

Validation: alpha=0 reproduces the stored refit **bit-identically** on all five lines and all three
seeds (difference 0.00e+00), in both runs. The zero-context arm is exactly +0.0000 on 0/15,
confirming no leakage.

Nothing clears the ±0.010 seed floor in either direction.

## The required context-specificity criterion: FAILED

| line | correct | wrong (real biology) | corrupted | correct > wrong | correct > corrupted |
|---|---|---|---|---|---|
| K562 | -0.0010 | +0.0002 | +0.0003 | **no** | **no** |
| RPE1 | -0.0022 | +0.0015 | -0.0005 | **no** | **no** |
| HepG2 | -0.0012 | +0.0011 | -0.0043 | **no** | yes |
| Jurkat | -0.0001 | +0.0005 | -0.0018 | **no** | yes |
| KOLF | +0.0012 | -0.0009 | -0.0025 | yes | yes |

Correct context beats a real different biology's context on **1 of 5** folds. Across all 15 paired
seeds, wrong context wins 12/15 against refit while correct context wins 4/15.

The single pass is KOLF, the line whose perturbation deviation was measured at **99.3% noise**
(`basal_context_2026-09-20`). It is also the line where the tiny model draws by far the largest
correction, 31.8% of the refit norm against 3.5-9.0% elsewhere. A large confident correction on
the noisiest line is a warning, not a result.

## Capacity was a real problem, and fixing it was not enough

| | large | tiny |
|---|---|---|
| trainable parameters | 277,760 | **9,728** |
| of which 2058-d supervised token projections | 263,424 (95%) | 0 |
| validation beats zero residual | **7 of 20** | **10 of 10** |
| five-line PDS vs refit | -0.0013 | -0.0007 |
| correct beats real wrong context | 0/5 folds | 1/5 folds |

Replacing the supervised 2058->64 key/value projections with a **frozen PCA fitted on training
biologies' control cells only** (68.3-69.9% variance kept at 64 dimensions) removed the
memorisation entirely: every tiny training now beats a zero residual on a held-out *training*
biology, where the large model managed it in a third of cases.

The large model's failure signature was unambiguous. On K562 its train loss fell 0.0755 to 0.0695
while held-out-training-biology validation **rose** 0.1808 to 0.2047 from the first checkpoint, a
3x gap. The tiny model still overfits, more slowly: train 0.1343 to 0.1256 while validation rises
0.2970 to 0.3158 after step 200, a 2.2x gap.

So capacity explains the *validation* failure and roughly half the PDS deficit. It does not
explain the sign.

## What is actually happening

The learned residual has genuine **average** validity but its **context-specific** component is
anti-correlated with what PDS rewards. Three facts pin this down:

1. Feeding the model *any* plausible biology's context helps slightly (wrong context +0.0005,
   12/15 wins).
2. Feeding it the *correct* biology's context hurts (-0.0007, 4/15).
3. Feeding it a *corrupted* context hurts most (-0.0018).

A model that had learned nothing would land on the zero-context arm, at exactly +0.0000. Instead
the ordering is wrong ≈ mean-pool > zero > correct > corrupted. The model is using context; it is
using it in a direction that does not survive transfer to an unseen biology.

## Attention is not collapsed; the mechanism works

| line | tiny entropy / uniform | max weight | cross-perturbation spread |
|---|---|---|---|
| K562 | 2.631 / 2.773 | 0.145 | 0.057 |
| RPE1 | 2.586 / 2.773 | 0.200 | 0.049 |
| HepG2 | 2.763 / 2.773 | 0.081 | 0.008 |
| Jurkat | 2.546 / 2.773 | 0.219 | 0.041 |
| KOLF | 2.765 / 2.773 | 0.079 | 0.014 |

Different perturbations do attend to different prototypes; attention never fully collapses to a
single pattern. The large model was more concentrated still (K562 entropy 2.389 of 3.466, max
weight 0.383, spread 0.202) and scored worse. **More attention structure correlated with worse
transfer**, which is what memorisation looks like.

HepG2 is the near-collapse case (entropy 2.763 of 2.773, spread 0.008) and is also the fold where
early stopping fires at step 100 with the smallest possible margin over zero.

## Ablations within the architecture

| ablation | five-line | vs refit |
|---|---|---|
| cross-attention | 0.7992 | -0.0013 |
| **mean pooling instead** | **0.8012** | **+0.0007** |
| multiplicative interaction | 0.7992 | -0.0013 |
| additive instead | 0.7999 | -0.0006 |

**Mean pooling beats cross-attention by +0.0020**, and is the only STATE arm that is positive.
The perturbation-conditioned selection of which control cells matter is exactly the component
that costs performance.

## The five questions

1. **Was the previous failure mostly caused by excessive capacity?** Partly, not mostly. Capacity
   fully explains the validation failure (7/20 to 10/10 beat-zero) and about half the PDS deficit
   (-0.0013 to -0.0007). It does not explain why correct context is worse than wrong context,
   which persists unchanged at 28x smaller.
2. **Does frozen unsupervised STATE compression help?** Yes, measurably: it removes memorisation
   and halves the deficit. It does not make the model useful.
3. **Does tiny cross-attention learn real context-specific signal?** No. It learns *something*
   (it beats zero-residual validation on 10/10 and beats the zero-context arm's behaviour), but
   the context-specific part does not transfer: 1 of 5 folds pass, on the noisiest line.
4. **Does correct context beat wrong context on real PDS?** **No.** -0.0007 against +0.0005, and
   wrong context wins 12/15 paired seeds against correct context's 4/15.
5. **Continue or stop?** **Stop.** This was pre-declared the final architecture test for the
   direction, and it fails every clause of the decision rule.

## Which failure mode, specifically

Of the candidates named in advance:

- **Cross-attention overfitting** — confirmed for the large model, eliminated in the tiny one,
  and the result did not flip. Contributing cause, not root cause.
- **Prototype compression** — ruled out. Attention is structured, not collapsed, and 16 vs 32
  tokens changes nothing (-0.0007 vs -0.0014).
- **Low-rank response limitation** — ruled out. Rank 16 and rank 32 both fail; rank 32 is worse.
- **Uninformative STATE cell embeddings** — partly. A PCA of pooled controls (+0.0024) beats every
  STATE arm, so STATE is not adding usable information over a linear summary of the same cells.
- **Insufficient context diversity** — the binding constraint. Six training biologies per fold,
  and both the ESM arms result (the context branch was a 58-value source lookup) and the cell-FM
  result (a 600M model ties PCA) point the same way.
- **The interaction residual is simply not predictable from untreated RNA at this sample size** —
  the conclusion the evidence supports.

## Recommendation

Do not build a larger Transformer. The failure is not architectural: a 9,728-parameter model with
frozen unsupervised compression, clean early stopping and non-collapsed attention lands in the
same place as a 277,760-parameter one. The refit baseline at **0.8005** remains the best local
predictor and has still never been scored officially. The highest-value investment is more
independent biologies, not more interaction capacity.
