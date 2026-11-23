# Two authorised submissions: refit baseline, and refit + PCA interaction

Date: 2026-09-20. Authorised by the user ("submit the best two model", then "Refit + PCA
interaction"). Two uploads = the full daily quota. Build: `submit_build.py`.

## Results

| | model | entry | overall | rank |
|---|---|---|---|---|
| previous best | b19 | `pamfw9b1glTBvlXvEyQA` | 0.136470 | 193 at submission |
| **S1** | **b20-wk7f3a — refit baseline** | `3CFK0NSTG9Fw6SfjD3IC` | **0.137665** | **231** |
| S2 | b21-wk2e9c — refit + PCA interaction | `hlhJKNFjK8UHE1dq2rsV` | 0.135702 | 241 |

| metric | b19 | S1 refit | S2 +PCA | S1 - b19 | S2 - S1 |
|---|---|---|---|---|---|
| overall | 0.136470 | **0.137665** | 0.135702 | **+0.001195** | **-0.001963** |
| pds | 0.543156 | 0.545108 | 0.533791 | +0.001952 | **-0.011317** |
| lfc nmae | 0.075024 | 0.075883 | 0.074323 | +0.000859 | -0.001560 |
| fidelity | 0.018797 | 0.018175 | 0.018405 | -0.000621 | +0.000230 |
| reach | 0.188645 | 0.193821 | 0.194632 | +0.005176 | +0.000811 |
| jaccard | -0.006802 | -0.006997 | -0.006937 | -0.000195 | +0.000060 |
| expr mse | 0 | 0 | 0 | 0 | 0 |

## Two findings, one good and one that matters more

### 1. The refit baseline is the project's new best official score

**0.137665 vs b19's 0.136470, +0.001195.** It had been the best local predictor for days and had
never been submitted. Locally it beat b19 by +0.0047 five-line PDS; officially it beats it by
+0.0020 PDS and +0.0012 overall. **The direction transferred and the magnitude shrank by roughly
4x.** That is the opposite of b15, where the local signal pointed the wrong way entirely, so the
harness does carry real but attenuated predictive weight.

### 2. The local composite gain did NOT transfer, and reversed

The official vcc2026 re-scoring predicted the PCA interaction arm at **+0.0214 composite** over
refit on our harness. Officially it is **-0.001963 overall**: a loss.

| | local composite | official overall |
|---|---|---|
| refit | 0.2755 | 0.137665 |
| refit + PCA interaction | 0.2969 (**+0.0214**) | 0.135702 (**-0.0020**) |

**The sign flipped.** And the metric that carried the local gain is the one that broke: PDS fell
-0.011317, by far the largest single movement in either submission. Locally the PCA arm's gain was
mostly log-FC NMAE (+0.0839 normalised) with PDS contributing +0.0051; officially NMAE improved
(+0.001560 vs refit, in the right direction) while PDS collapsed and took the composite with it.

So the correction I made earlier -- moving from PDS-only to the full preset -- was right in method
and still produced a candidate that lost. Full-preset scoring on our harness is necessary but was
not sufficient.

## Why this most likely happened

The honest answer is that we do not know, but the evidence narrows it:

1. **Contexts A/B/C are genuinely unseen biologies with no measured perturbations, and our seven
   training biologies are few.** Every context result in this project has been capped by the same
   constraint. The interaction was fitted on 7 biologies and asked to extrapolate to 3 new ones.
2. **The local harness lines are not exchangeable with A/B/C.** The local gain was concentrated on
   RPE1 (+0.0853) and HepG2 (+0.0188), two of five lines. Nothing guarantees A/B/C resemble those.
3. **Our local composite uses locally-built anchors**, not the leaderboard's. The procedure is
   official; the normalisation constants are not. A +0.0214 local composite is not a +0.0214
   leaderboard prediction, and this pair is now direct evidence of that gap.

## What is established now

- **Best official model: the refit baseline, 0.137665, entry `3CFK0NSTG9Fw6SfjD3IC`.**
- The harness predicts leaderboard *direction* for a large, structural change (refit vs b19) and
  failed to predict it for a small context-conditioned correction (+PCA).
- Context-conditioned residuals have now failed on the real leaderboard, not just locally. Across
  this project: cosine benchmarks said +22% of oracle headroom, harness PDS said noise, the full
  official preset said +0.0214, and the leaderboard says -0.0020. **Each time the measurement got
  closer to the real objective, the gain shrank.** That ordering is itself the result.

## Build provenance

Both packages: 360,000 cells x 18,533 genes, 120,000 per context, 400 cells per target, scale 3,
stochastic rounding, seed 0, targets verified against the official list, counts preserved.
Distinct SHA256 (`9e9e8b3d...` and `b9d3e2dd...`) -- checked explicitly because b18 was silently
byte-identical to b16 and wasted a slot.

S1 is one context-independent delta reused for A/B/C. S2 carries three distinct per-context
deltas (pairwise L2 0.15-0.38), with the interaction applied to 285 of 300 official targets; the
other 15 are absent from the basal panel and fall back to pure refit. Context axis = the 6,413-gene
common axis intersected with the official features = 6,092 genes. Interaction only; no generic
per-biology shift.
