# Reliable-gene subsets do not improve transfer — 2026-09-19

Question: does restricting the model to the most reliably perturbation-responsive genes improve
held-out cell-line performance? Answer: **no, at any size tested.** Local reconstructed-count PDS
on the frozen evaluation path, five held-out lines, seeds 0/1/2. No submission.

## Setup

Universe is the current frozen 5,592-gene table, so this is a clean subset test against the
standing baseline. Per held-out line, genes are ranked on that fold's permitted rows ONLY
(destination biology already excluded) by `between-target variance x cross-biology ICC`:
responsive AND reproducible. Genes outside the subset keep the empirical cross-source mean,
exactly as unmodelled genes are treated today. Every size has a SIZE-MATCHED RANDOM control.
Validation: the all-genes arm reproduces the frozen refit at 0.8005 / K562 0.8821 exactly.

## Result 1 — no subset improves (the question asked)

Ranked subset minus all 5,592 genes:

| K | K562 | RPE1 | HepG2 | Jurkat | KOLF | mean | lines | paired |
|---|---|---|---|---|---|---|---|---|
| 500 | -0.0589 | +0.0091 | -0.0703 | -0.0168 | -0.0102 | -0.0294 | 1/5 | 3/15 |
| 1000 | -0.0412 | +0.0091 | -0.0465 | -0.0094 | -0.0081 | -0.0192 | 1/5 | 3/15 |
| 2000 | -0.0224 | +0.0035 | -0.0216 | -0.0035 | -0.0064 | -0.0101 | 1/5 | 3/15 |
| 3000 | -0.0121 | +0.0001 | -0.0107 | -0.0014 | -0.0011 | -0.0050 | 1/5 | 2/15 |
| 5000 | -0.0028 | -0.0010 | -0.0021 | +0.0009 | +0.0022 | -0.0006 | 2/5 | 6/15 |

Every size fails the project gate (gain>=.001, >=3/5 lines, >=9/15 paired). The deficit shrinks
monotonically toward zero as K approaches the full set: **more genes is simply better.** RPE1 is
the lone exception and is also the line with the weakest full baseline (0.8107), so it has the
most room for restriction to help — a confound, not support.

**KOLF, the closest analogue to the unseen contexts:** -0.0102 / -0.0081 / -0.0064 / -0.0011 at
500-3000, and +0.0022 only at K=5000, which is 89% of the full set. No useful subset.

## Result 2 — but the ranking itself works (worth keeping)

Ranked minus size-matched random:

| K | K562 | RPE1 | HepG2 | Jurkat | KOLF | mean |
|---|---|---|---|---|---|---|
| 500 | +0.0425 | +0.0080 | +0.0498 | +0.0215 | **-0.0011** | +0.0241 |
| 1000 | +0.0327 | +0.0069 | +0.0627 | +0.0178 | **+0.0004** | +0.0241 |
| 2000 | +0.0168 | -0.0026 | +0.0537 | +0.0110 | +0.0026 | +0.0163 |
| 3000 | +0.0045 | -0.0049 | +0.0357 | +0.0047 | +0.0074 | +0.0095 |

Given that you must pick 500 genes, the reliability statistic picks much better ones than chance
(+0.024 mean). So the metric measures something real. It just cannot pay for the information
destroyed by discarding the other 5,092 genes.

**The important caveat: the ranking does NOT help on KOLF** (-0.0011 at 500, +0.0004 at 1000)
while helping strongly on K562/HepG2/Jurkat. Cross-biology reliability estimated from our 14
source biologies transfers to lines that resemble those sources and NOT to the one that does not.
This is the same pattern as every other result in this project.

Top-ranked genes are tubulins, mitochondrial genes and proliferation markers (MDM2, STMN1, PCLAF,
CENPF) — a coherent fitness/growth axis, which is real biology but evidently source-specific in
its detail.

## Conclusion

Do not pursue gene selection. The model wants MORE output coordinates, not fewer, and the
reliability ranking's usefulness is confined to lines that resemble the training sources. Best
local predictor remains the frozen five-line refit mean of **0.8005** (b19 0.7958).
