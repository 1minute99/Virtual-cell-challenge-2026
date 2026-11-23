# Strict unseen-context benchmark — can a context representation recover oracle headroom?

Date: 2026-09-20. Scripts: `context_benchmark.py`, `context_benchmark_decomp.py`.
Data: `reports/basal_context_2026-09-20/`, truth from `reports/oracle_2026-09-20/`.

## Design

Every fold removes **all** perturbation data for the held-out biology. At test time the model
sees only that biology's pooled control profile. The perturbation mean `mu`, the response basis,
the context basis, the generic map and the bilinear interaction map are fitted on training
biologies only.

Predictions are scored against the oracle experiment's **batch-held-out** replicate B, so the
held-out biology's perturbation data never touches training, and the three quantities are on
one scale:

```
universal  cos(mu,                   truth_B)   no context
context    cos(mu + delta_predicted, truth_B)
oracle     cos(effect_A,             truth_B)   empirical ceiling, replicate A of same line
recovered  = (context - universal) / (oracle - universal)
```

All context arms are projected to the same dimensionality (8) before the interaction head and
share the same perturbation representation, response space, architecture, ridge penalty and
training budget. Ablations per arm: correct context, wrong context (a training biology's
vector), shuffled coordinates, and zeros.

## Headline result — context beats no-context

n = 649 held-out (biology, target) pairs across HepG2, Jurkat, RPE1.

| arm | correct | 95% CI on gain | wrong | shuffled | zero | gain vs no-context | oracle gain recovered |
|---|---|---|---|---|---|---|---|
| no-context (`mu`) | 0.2377 | — | — | — | — | — | — |
| `none` (sanity) | 0.2396 | [+0.0006, +0.0031] | 0.2396 | 0.2396 | 0.2396 | +0.0019 | 1.6% |
| `raw` | 0.2631 | [+0.0215, +0.0294] | 0.1999 | 0.2260 | 0.2386 | +0.0254 | 21.2% |
| **`pca`** | **0.2645** | **[+0.0239, +0.0299]** | 0.2246 | 0.2304 | 0.2376 | **+0.0268** | **22.3%** |
| oracle ceiling | 0.3578 | — | — | — | — | +0.1201 | 100% |

Bootstrap: 4,000 resamples over pairs.

### Correct context beats wrong and shuffled

| arm | correct − wrong | correct − shuffled | pairs won vs wrong |
|---|---|---|---|
| `raw` | +0.0633 | +0.0371 | 465 / 649 |
| `pca` | +0.0399 | +0.0342 | 499 / 649 |

The `none` arm gains 1.6% and its four ablations are bit-identical, which is the intended
no-leak sanity check: with a zero context vector the pipeline cannot express biology-specific
structure and does not.

### Per-biology

Correct-context gain over no-context, and the same as a fraction of that fold's own headroom:

| arm | HepG2 (headroom +0.0694) | Jurkat (+0.0841) | RPE1 (+0.1940) | folds won |
|---|---|---|---|---|
| `none` | +0.0000 (0%) | +0.0020 (2%) | +0.0030 (2%) | 3/3 |
| `raw` | +0.0266 (38%) | −0.0019 (−2%) | +0.0541 (28%) | 2/3 |
| `pca` | +0.0294 (42%) | +0.0084 (10%) | +0.0449 (23%) | 3/3 |

PCA context is the only arm positive on all three held-out biologies. Jurkat is the weak fold
for both real-context arms.

## Ranked comparison

1. **`pca`** — +0.0268, 22.3% of oracle headroom, 3/3 folds, beats wrong context by +0.0399.
2. **`raw`** — +0.0254, 21.2%, 2/3 folds, beats wrong context by +0.0633 but is negative on Jurkat.
3. **`none`** — +0.0019, 1.6%, no ablation contrast. Behaves exactly as a null arm should.

Raw and PCA are within each other's confidence intervals on the aggregate gain. PCA is
preferred on fold consistency, not on aggregate size.

## Cell foundation models were not benchmarked

The requested comparison against GB.Cell, UCE, scFoundation and scGPT could not be run, for
reasons of artifact availability rather than compute:

| model | status |
|---|---|
| GB.Cell | no such repository on Hugging Face |
| scGPT | canonical repo returns 401 unauthenticated |
| UCE (`minwoosun/UCE-100m`, `minwoosun/uce-650m`) | weights + config downloadable, **no vocabulary/tokenizer or inference stack** |
| scFoundation (`genbio-ai/scFoundation`) | same: weights + config only |
| Geneformer (`ctheodoris/Geneformer`) | ungated and complete, not in the requested list |
| scLong | 12 GB of weights already local at `checkpoints/sclong/` |

Embedding cells with UCE or scFoundation requires reconstructing their gene vocabulary and
tokenisation exactly; guessing it would produce a representation that is not the published
model and would make the comparison meaningless. This is a real gap in the benchmark, and the
decision below rests only on the arms that did run.

## Decision

The user's rule was: continue if context representations beat no-context and correct context
beats wrong/shuffled on held-out biologies; stop and prioritise more biologies if everything
fails despite large oracle headroom.

**Both conditions are met by the arms that ran.** Raw and PCA context beat no-context with
confidence intervals excluding zero, and correct context beats both wrong and shuffled context.
So the evidence favours continuing to develop the context-conditioned model.

### This refutes the prediction made before the run

Before running, the stated expectation was: *"the context arms will fail to beat no-context,
because the constraint is eight training context vectors rather than their coordinates."*
They did not fail. The prediction was wrong and the earlier recommendation to abandon context
does not survive this benchmark or the oracle result.

The reconciliation with the earlier null is that `reports/basal_context_2026-09-20/` measured
R² of the deviation itself and found it negative, while this benchmark measures cosine of the
total prediction. Both hold at once: the predicted deviation is a poor magnitude estimate, yet
adding it still improves directional alignment of the total. The earlier null was
over-generalised from a single narrower statistic.

## Follow-up: is the gain real, or an artefact of a weak baseline?

The headline no-context baseline is the plain cross-biology mean `mu`, which is mean transfer.
Our deployed predictor is the refit baseline, which beats mean transfer. A gain measured only
over `mu` could therefore be absorbed by a better context-independent baseline. Two checks,
both fitted on training biologies only (`context_benchmark_decomp.py`, same folds, same truth).

### Stronger no-context baselines

| baseline | cosine | vs plain `mu` |
|---|---|---|
| `mu` | 0.2377 | — |
| `mu` × best global scale | 0.2377 | +0.0000 |
| `mu` projected onto a rank-32 target basis | 0.2533 | +0.0156 (13% of headroom) |

The fitted global scale came out at exactly 1.0 on all three folds, so the context gain is not
a magnitude correction. Low-rank denoising is the mechanism by which the refit baseline beats
mean transfer, and it recovers 13% of headroom on its own.

### The gain survives the stronger baseline

| baseline | + context | gain | 95% CI | pairs won |
|---|---|---|---|---|
| `mu` | 0.2377 → 0.2646 | +0.0269 | [+0.0238, +0.0300] | 499/649 |
| `mu` × scale | 0.2377 → 0.2646 | +0.0269 | [+0.0238, +0.0301] | 499/649 |
| `mu` rank-32 denoised | 0.2533 → 0.2761 | **+0.0227** | [+0.0187, +0.0267] | 454/649 |

About 84% of the context gain is additive on top of the denoising that makes the refit baseline
better than mean transfer, and it stays positive on **all three** held-out biologies
(HepG2 +0.0123, Jurkat +0.0090, RPE1 +0.0448). Denoising and context are largely orthogonal:
together they reach 0.2761 against a 0.3578 ceiling, or 32% of the headroom.

### What kind of signal is it?

| component | cosine | gain | 95% CI | recovered |
|---|---|---|---|---|
| generic per-biology shift `g` only | 0.2538 | +0.0161 | [+0.0134, +0.0186] | 13.4% |
| target × context interaction `s` only | 0.2514 | +0.0137 | [+0.0119, +0.0156] | 11.4% |
| both | 0.2646 | +0.0269 | [+0.0238, +0.0300] | 22.4% |

The two components are nearly additive, so roughly half the gain is a per-biology offset that
applies to every target alike, and half is a genuine target-specific interaction. The
interaction term is the part that a context-conditioned model is supposed to buy, and it is
positive with a confidence interval excluding zero. It is also the smaller half, which is worth
remembering before attributing the full 22% to context-specific biology.

## Recommendation

1. **Continue the context-conditioned model.** Both of the user's conditions are met, and the
   gain survives a baseline strengthened the same way the refit baseline is strengthened.
2. **Use PCA context over raw projection.** Same aggregate gain, positive on 3/3 folds rather
   than 2/3.
3. **Fit the generic shift and the interaction as separate terms.** They are additive here and
   the cheap generic shift delivers over half the gain on its own.
4. **Treat 22% as the realistic near-term target, not 100%.** 78% of oracle headroom is
   unrecovered. Eight training context vectors is a hard sample-size constraint on how much
   context structure can be learned, so adding independent biologies remains the highest-value
   data investment even though the "stop and add biologies" branch of the decision rule was not
   triggered.
5. **Not yet established: that this converts into leaderboard score.** Everything here is cosine
   against a batch-held-out replicate. The mapping from cosine to the six-metric composite is
   not linear, and PDS is a discrimination rank that row cosine is known to mislead on. Scoring
   the context predictions through cell-eval2 is the necessary next step before any submission.

## Supporting control: the oracle is not a leakage artefact

The random-split control finished alongside this benchmark. Splitting cells at random instead of
by batch gives essentially the same oracle gain, so shared technical noise is not inflating it:

| line | batch-split gain | random-split gain | inflation |
|---|---|---|---|
| HepG2 | +0.1736 | +0.1707 | 0.98× |
| Jurkat | +0.1423 | +0.1348 | 0.95× |
| RPE1 | +0.2676 | +0.2665 | 1.00× |

Aggregate batch-split oracle 0.3812 vs universal 0.1859, +0.1953, winning 569/649 pairs.

## How much is 22% of headroom worth in score units?

The oracle phase-2 run scores the same arms through cell-eval2 and normalises each metric against
the published anchors (0 = context-mean baseline, 1 = split-half replicate). All three lines
complete.

| line | arm | composite | pds | expr | lfc nmae | fidelity | reach | jaccard |
|---|---|---|---|---|---|---|---|---|
| HepG2 | oracle | 0.5875 | 0.886 | 0.283 | 0.747 | 0.488 | 0.978 | 0.143 |
| HepG2 | universal | 0.2175 | 0.709 | 0.043 | 0.125 | 0.000 | 0.354 | 0.074 |
| HepG2 | wrong | 0.1512 | 0.467 | 0.000 | 0.117 | 0.000 | 0.288 | 0.035 |
| HepG2 | none | 0.0710 | 0.000 | 0.000 | 0.000 | 0.000 | 0.058 | 0.368 |
| Jurkat | oracle | 0.4956 | 0.824 | 0.000 | 0.760 | 0.426 | 0.928 | 0.036 |
| Jurkat | universal | 0.2273 | 0.555 | 0.000 | 0.174 | 0.000 | 0.371 | 0.263 |
| Jurkat | wrong | 0.1641 | 0.298 | 0.000 | 0.178 | 0.000 | 0.292 | 0.217 |
| Jurkat | none | 0.0796 | 0.000 | 0.000 | 0.022 | 0.000 | 0.079 | 0.377 |
| RPE1 | oracle | 0.6626 | 0.922 | 0.451 | 0.784 | 0.607 | 0.999 | 0.213 |
| RPE1 | universal | 0.1870 | 0.626 | 0.000 | 0.091 | 0.000 | 0.196 | 0.209 |
| RPE1 | wrong | 0.1695 | 0.451 | 0.000 | 0.139 | 0.000 | 0.211 | 0.216 |
| RPE1 | none | 0.0344 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.206 |

Composite headroom from no-context to oracle, and what a linear 22% recovery would be worth:

| line | no-context | oracle | headroom | wrong vs no-context | 22% of headroom |
|---|---|---|---|---|---|
| HepG2 | 0.2175 | 0.5875 | +0.3700 | -0.0663 | +0.0825 |
| Jurkat | 0.2273 | 0.4956 | +0.2683 | -0.0632 | +0.0598 |
| RPE1 | 0.1870 | 0.6626 | +0.4755 | -0.0175 | +0.1060 |

Wrong context scores *below* no-context on all three lines, matching the cosine ablation. The
margin is small on RPE1 (-0.0175), which is the line where the cosine context gain was largest;
the two do not have to agree, since composite is dominated by pds and reach rather than row cosine.

For scale, our published entry is 0.1365 and the top entry is 0.3304, so a linear 22% would be a
large move. **That extrapolation is not evidence.** The cosine-to-composite mapping is not linear,
and PDS is a discrimination rank that row cosine is known to mislead on. Note that the `none` arm
scores 0.03-0.08 composite while sitting at 0.24 row cosine, which is exactly that failure mode.
It is offered only to show the question is worth the cost of answering properly, which means
scoring actual context-model predictions through cell-eval2.
