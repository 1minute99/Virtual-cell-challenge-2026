# Does the stable basal state of a cell line predict its perturbation-response deviation? — 2026-09-20

Re-run of the context question with the definition corrected. The previous analysis used
per-target matched-control fluctuations WITHIN a line; those are ~60% sampling noise and are not
what the challenge asks. Here context is the robust POOLED control profile of a whole biology,
built from all of its untreated cells.

## Setup

9 biologies x 11,058 shared perturbations x 5,592 genes (>=800 targets per biology, each target in
>=4 biologies). Decomposition, with `mu` from TRAINING biologies only:

    effect[i,p] = mu[p] + delta[i,p]        mu = mean transfer / refit behaviour
    delta[i,p]  = g[i]  + s[i,p]            g = generic line shift, s = perturbation-specific

`s` is the only part the scorer can reward -- the metric deletes the generic all-target response.
Models are linear / low-rank / bilinear only: a ridge map from basal latent to `g`, and a bilinear
(basal latent x perturbation latent) map to `s`. Leave-one-biology-out, plus a stricter
leave-one-LAB-out that breaks the rn={K562,RPE1,HepG2,Jurkat} and orion={HCT116,HEK293T} groupings.

## The two controls that frame everything

**1. 94% of the deviation is perturbation-specific.** Mean generic share 0.058 (range 0.003-0.253).
So the half that context might plausibly explain is the small half.

**2. Only ~28% of the deviation is reproducible at all.** Fitting
`||delta_p||^2 = a + b/effective_n_p` across perturbations and taking the intercept as the
noise-free part:

| biology | signal fraction | | biology | signal fraction |
|---|---|---|---|---|
| HepG2 | 0.681 | | K562 | 0.180 |
| CD4 | 0.565 | | RPE1 | 0.165 |
| Jurkat | 0.344 | | HCT116 | 0.163 |
| FENG | 0.199 | | **KOLF** | **0.007** |
| HEK293T | 0.192 | | **mean** | **0.277** |

Deviation magnitude correlates 0.40-0.91 with 1/effective_n. **KOLF's deviation is 99.3% noise** --
which retrospectively explains the unexplained "+0.019 on KOLF" anomaly from the companion
modelling run: there is essentially nothing reproducible there to predict.

## Results

Leave-one-biology-out, context rank sweep (mean over 9 folds):

| context rank | generic: ctx | wrong | zero | ctx>wrong | specific: ctx | zero |
|---|---|---|---|---|---|---|
| 1 | -0.062 | -0.473 | -0.010 | 7/9 | -0.0006 | 0.0006 |
| 2 | -0.068 | -0.316 | -0.009 | 7/9 | -0.0016 | 0.0007 |
| 3 | -0.148 | -0.376 | -0.010 | 7/9 | -0.0024 | 0.0007 |

Leave-one-LAB-out (rank 2): generic ctx **+0.078**, wrong -0.244, zero -0.007;
ctx>wrong 6/9, ctx>zero 5/9. Specific ctx +0.0010 vs zero +0.0006, ctx>zero 3/9.
Per-fold generic R2 ranges from -0.046 (K562) to +0.373 (HEK293T) -- extremely variable.

## Verdict

**1. How much variance is explained?** Essentially none. On the perturbation-specific 94%, R2 is
-0.002 (LOBO) to +0.001 (lab-out) against a noise ceiling of 0.28 -- i.e. context captures none of
the reproducible deviation. On the generic 6%, R2 is -0.062 to +0.078 depending on protocol,
never reliably above zero.

**2. Does it generalise to a completely held-out biology?** No. Correct context beats predicting
nothing in 5 of 9 folds -- a coin flip -- under both protocols and every rank.

**3. Correct vs shuffled / wrong / no context?** Mixed, and this is the one positive.
**Correct context beats WRONG context 7/9 (LOBO, all three ranks) and 6/9 (lab-out)**, a
reproducible ordering: the basal profile does carry directional information about how a line's
responses shift. But correct context does NOT beat ZERO context. Directional signal exists; it is
worth less than doing nothing.

**4. Large enough to improve the refit baseline?** No. The seed noise floor on the scoring harness
is +-0.010 PDS; a channel that cannot beat a zero prediction on its own target variable cannot move it.

**5. Continue or abandon?** **Abandon the context branch for now** -- but for a specific, stated
reason rather than "context is useless". Reducing the context to a SINGLE number per cell line
still fails to beat zero, so this is not an under-regularisation problem that a better model fixes.
The binding constraints are in the data:
  * 9 biologies means at most 8 context vectors to learn any context->response map from;
  * 72% of the deviation is sampling noise;
  * 94% of what remains is perturbation-specific, where context shows nothing at all.

**What would change the answer:** more independent biologies (the map is estimated from 8 points),
and deeper coverage per condition (the target variable is mostly noise). Neither is a modelling
change. Until then, redirect to (a) submitting the refit baseline -- best local predictor 0.8005 vs
b19's 0.7958, still never scored officially; and (b) expression accuracy, where we score 0 and
which is 38% of the leaderboard gap.

**Correction to the earlier report:** the previous conclusion ("context signal is sampling noise,
R2~0.001") was measured on within-line matched-control fluctuations and was the wrong test. This
analysis uses pooled basal states as intended. The conclusion is similar but the reasoning is
different, and one finding is new: correct context reliably beats wrong context on the generic
shift, so the signal is real but too weak to use.
