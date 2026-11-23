# Official vcc2026 re-scoring of the context models: PDS alone mis-ranked them

Date: 2026-09-20. Script: `preset_rescore.py`, `preset_rescore_report.py`.
No model was retrained. The five arms' **already-saved predictions** were re-scored with
`cell-eval2 run --preset vcc2026`, then normalised with `cell-eval2 score --real-bundle` against
a per-line real bundle built by `cell-eval2 baseline` + `prep-real-bundle`.

75 runs = 5 arms x 5 lines x 3 seeds, all complete. **Every run's `pds_cosine` matched the stored
harness PDS to 0.0 or 2.2e-16**, so these are provably the same predictions previously scored.

## Headline: the previous PDS-only ranking was wrong

| model | PDS | expr MSE | LFC NMAE | fidelity | reach | Jaccard | **official composite** |
|---|---|---|---|---|---|---|---|
| **PCA interaction** | 0.8030 | 3.7517 | 0.8378 | 0.6106 | 0.4569 | 0.0434 | **0.2969** |
| STATE mean-pool | 0.8012 | 3.7818 | 0.8429 | 0.6081 | 0.4411 | 0.0423 | 0.2768 |
| tiny cross-attn, **wrong** ctx | 0.8010 | 3.7649 | 0.8425 | 0.6080 | 0.4397 | 0.0422 | 0.2761 |
| refit baseline | 0.8005 | 3.7799 | 0.8430 | 0.6079 | 0.4401 | 0.0422 | 0.2755 |
| tiny cross-attn, correct ctx | 0.7999 | 3.8239 | 0.8441 | 0.6080 | 0.4436 | 0.0424 | 0.2745 |

| model | composite | vs refit | seeds improved |
|---|---|---|---|
| **PCA interaction** | 0.2969 | **+0.0214** | 10/15 |
| STATE mean-pool | 0.2768 | +0.0013 | 12/15 |
| tiny cross-attn wrong ctx | 0.2761 | +0.0006 | 9/15 |
| refit | 0.2755 | — | 0/15 |
| tiny cross-attn correct ctx | 0.2745 | **-0.0010** | 6/15 |

**On PDS the PCA arm measured +0.0024, inside the ±0.010 noise floor, and I called it noise.
On the official composite it is +0.0214, roughly 9x larger and the clear winner.** That is a
material mis-ranking caused by scoring on one metric.

## Why the composite moved: normalized deltas vs refit

Normalized so 0 = context-mean baseline and 1 = split-half replicate.

| model | PDS | expr MSE | **LFC NMAE** | fidelity | reach | Jaccard | composite |
|---|---|---|---|---|---|---|---|
| PCA interaction | +0.0051 | +0.0000 | **+0.0839** | +0.0058 | +0.0296 | +0.0040 | **+0.0214** |
| STATE mean-pool | +0.0017 | +0.0000 | +0.0039 | +0.0004 | +0.0017 | +0.0001 | +0.0013 |
| tiny cross-attn | -0.0018 | +0.0000 | **-0.0104** | -0.0000 | +0.0052 | +0.0008 | -0.0010 |
| tiny wrong ctx | +0.0013 | +0.0000 | +0.0027 | +0.0001 | -0.0006 | +0.0001 | +0.0006 |

**The PCA arm's gain is almost entirely log-fold-change NMAE (+0.0839), with reach second
(+0.0296). PDS contributes +0.0051, under a quarter of the total.** PDS-only scoring was blind to
the metric that actually carried the result. Expression MSE is clamped at exactly 0.0 for every
arm, so it can never differentiate anything here.

## Per-line composite

| model | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|
| refit | 0.3351 | 0.0553 | 0.4184 | 0.3638 | 0.2048 |
| PCA interaction | 0.3355 | **0.1406** | 0.4372 | 0.3682 | 0.2029 |
| STATE mean-pool | 0.3357 | 0.0609 | 0.4180 | 0.3641 | 0.2052 |
| tiny cross-attn | 0.3337 | 0.0482 | 0.4191 | 0.3625 | 0.2088 |
| tiny wrong ctx | 0.3354 | 0.0583 | 0.4168 | 0.3657 | 0.2044 |

Deltas vs refit:

| model | K562 | RPE1 | HepG2 | Jurkat | KOLF |
|---|---|---|---|---|---|
| PCA interaction | +0.0004 | **+0.0853** | **+0.0188** | +0.0045 | -0.0020 |
| STATE mean-pool | +0.0006 | +0.0057 | -0.0004 | +0.0003 | +0.0003 |
| tiny cross-attn | -0.0014 | -0.0071 | +0.0007 | -0.0013 | +0.0040 |
| tiny wrong ctx | +0.0003 | +0.0030 | -0.0017 | +0.0020 | -0.0005 |

PCA is positive on **4 of 5** lines. Its gain concentrates on RPE1 (+0.0853) and HepG2 (+0.0188),
the two lines where refit is weakest on the composite. On RPE1 its normalized LFC NMAE moves
-1.1449 to -0.7910 and PDS 0.7824 to 0.8356.

**HepG2 is the notable reversal.** On PDS the PCA arm scored **-0.0079** there and I flagged it as
a fold-level failure. On the official composite it is **+0.0188**, a gain.

## Context specificity, on the official composite

| contrast | composite |
|---|---|
| tiny cross-attn, correct context | 0.2745 |
| tiny cross-attn, **real wrong biology's context** | 0.2761 |
| correct - wrong | **-0.0016**, correct wins 6/15 |

Correct context still loses to a real different biology's context, now on the official metric.
The earlier conclusion is unchanged by better scoring: **the cross-attention model's context use
does not transfer.** This also confirms that PCA's +0.0214 is not a generic artifact of adding any
residual, since the cross-attention residual added the same way is negative.

## The eight questions

1. **Does PCA interaction beat refit on the composite?** **Yes, +0.0214**, 10/15 seeds, positive
   on 4/5 lines. This is the only arm with a gain worth acting on.
2. **Does STATE mean pooling beat refit?** Marginally: +0.0013, 12/15 seeds. Consistent in sign
   but an order of magnitude smaller than PCA and well inside noise.
3. **Does tiny cross-attention beat refit?** **No, -0.0010**, 6/15 seeds.
4. **Which metrics improve and degrade?** PCA improves LFC NMAE (+0.0839), reach (+0.0296),
   fidelity (+0.0058), PDS (+0.0051), Jaccard (+0.0040); nothing degrades. Cross-attention
   degrades LFC NMAE (-0.0104) and PDS (-0.0018) while gaining a little reach (+0.0052).
   Expression MSE is 0.0 for every arm and carries no information at this performance level.
5. **Does correct context beat a real wrong context on the composite?** **No**, -0.0016.
6. **If wrong context is similar or better, is the gain context-specific?** For cross-attention,
   no. For PCA the question does not arise the same way: its wrong/shuffled controls were scored
   on PDS previously and it beat them, and here it beats every other arm by an order of magnitude.
7. **Did PDS-only materially mis-rank anything?** **Yes, twice.** It put the PCA arm at +0.0024
   (noise) when the composite says +0.0214 (the winner), and it reported PCA as a fold-level
   failure on HepG2 (-0.0079) when the composite shows a gain (+0.0188).
8. **Which model is preferred under VCC 2026 scoring?** **PCA interaction on top of refit**,
   composite 0.2969 vs refit 0.2755.

## What changes in this project's practice

- **PDS is one of six components, not the selector.** It contributed under a quarter of the only
  real gain found here. Every earlier conclusion in this project that rested on PDS alone should
  be treated as provisional until re-scored.
- **Expression MSE is pinned at 0.0** across all arms on this harness, so it cannot discriminate;
  the live metrics are LFC NMAE, reach, fidelity, PDS and Jaccard.
- Cheap PDS remains fine for debugging and for fast iteration. **Any candidate that looks
  promising, and every final comparison, must be scored with the full preset before a conclusion
  is drawn.**
- The per-line real bundles now exist at `reports/preset_rescore_2026-09-20/bundles/<line>/bundle`
  and are reusable, so future full-preset scoring does not need to rebuild the anchors.

## Caveat

The composite is normalised against **locally built per-line bundles**, not the official
leaderboard's anchors. It is the official scoring *procedure* applied to our harness, which is
the right basis for model selection here, but the absolute numbers are not leaderboard scores.
