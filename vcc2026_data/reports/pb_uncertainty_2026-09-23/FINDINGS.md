# Pseudobulk construction and within-condition uncertainty — 2026-09-23 (NEGATIVE, stopped by user)

Upstream-data experiment on the extended table (table_v2) + refit_v2 (LeanSourceBaseline rank 32),
five leave-one-line-out holdouts, seed 0, scale 3 / lambda 1. Everything downstream held fixed.
Code: `scripts/pb_cellstats.py` (cell-level cache), `pb_arms2.py` (arms), `pb_confidence.py`,
`pb_uncertainty_cd4.py`, `pb_score.sh`. Caches: `cells/{feng,kolf,orion}/`, `cd4/`, `arms2/`.

## Coverage
- Cell-level rebuild: Feng targeted lines (all 7,435 rows), KOLF (all 11,266), Orion (1,079
  evaluation-target rows only; other 34,360 Orion rows unchanged). = 2,301 / 8,446 (27 %) of the
  genetic rows of evaluation targets. The rebuilt current estimator reproduces the table exactly
  (corr 1.0, max |diff| 1.2e-7) for all three; group cell counts match the stores.
- No cell-level data: CD4 (author guide x run pseudobulk), RN, pooled Feng, xl/feng, Jiang (per-rep
  pseudobulk); H1 has cells but was not rebuilt (22 eval rows). Summaries unchanged.
- SE: cells = analytic within-experiment cell sampling (ratio estimator, half-count Poisson floor,
  shared control term); CD4 = guide x run spread (111,080 rows; 12,463 single-row targets get none);
  Jiang = replicate pair. Never cross-biology spread.

## Cheap delta evaluation (five-line mean vs base; every arm 0/5 lines on corr)
| arm | d corr | d mse_ratio | d sign prec | d strong recall | d prec@top5% |
|---|---|---|---|---|---|
| weight (gene-wise, clip 1/4..4 of cell-count weight) | -0.0006 | +0.002 | -0.0002 | -0.0003 | -0.0007 |
| umean | -0.0017 | +0.005 | -0.0004 | -0.0008 | -0.0005 |
| median | -0.0293 | +0.119 | -0.0039 | -0.0089 | -0.0370 |
| winsor99 / winsor995 | -0.0017 | +0.004 | -0.0004 | -0.0007 | -0.0008 |
| trim1 / trim5 | -0.0024 / -0.0020 | +0.006 | -0.0007 | -0.0014 | -0.0019 / -0.0013 |
| nested robust (chosen on other 4 lines: umean/winsor99) | -0.0017 | +0.005 | -0.0004 | -0.0008 | -0.0006 |
| nested + weight | -0.0022 | +0.007 | -0.0005 | -0.0008 | -0.0011 |

Median collapses effects to exactly zero: Feng 42 %, KOLF 84 %, Orion 34 % of (row, gene).
Winsorised ~= unweighted mean (corr to current 0.90-0.93 for both): the difference from the current
effect is removing library weighting, not tail capping. Winsor99 is a no-op for n < 100 (75 % of Feng groups).
CD4 multipliers hit a cap 43 % of the time (variance from ~2 guide rows, ~1 df); cell-based
multipliers sit 0.65-1.5 (5-95 %), i.e. the cell-count weight was already ~inverse-variance.

## Official vcc2026 composite, seed 0 (delta vs base)
| arm | K562 | RPE1 | HepG2 | Jurkat | KOLF | mean |
|---|---|---|---|---|---|---|
| base (= refit_v2) | 0.3481 | 0.0447 | 0.4270 | 0.3772 | 0.2101 | 0.2814 |
| weight | -0.0006 | -0.0000 | -0.0014 | -0.0003 | -0.0011 | -0.0007 |
| nested robust | -0.0003 | +0.0011 | -0.0014 | -0.0018 | +0.0004 | -0.0004 |
| nested + weight | -0.0000 | +0.0007 | (stopped) | | | |

Gains attributable to robust aggregation: none; to uncertainty weighting: none; to both: none.

## Directional confidence (nested logistic, prediction = refit_v2)
Within-experiment z features (u_z, u_sig, u_t, u_nb; ~2.9 of ~6.8 biologies per target have SE):
over |pred| alone AUC +0.0056 (5/5), top-20 % |truth| AUC +0.0137 (5/5); over the cross-biology
features of confidence_gate AUC +0.0014 (5/5), top-20 % +0.0036 (5/5), sign precision in top 10 %
+0.0025 (4/5). Real but small; confidence shrinkage already lost officially (confidence_gate_2026-09-23).

## Conclusion
Pseudobulk construction and within-condition reliability weighting are not a lever for this
pipeline: the current library-weighted pseudobulk with cell-count weights is already close to the
efficient estimator where cells exist, and the sources without cells (CD4, RN) dominate the
evaluation targets. Close the direction; no leaderboard slot.
Ops note: a one-pass Feng store hit 40 GB and oomd killed it (and Firefox) at 09:02; stream instead.
recon_sweep.py needs `--threads 2` here (numba max 2), which also explains the 08:25 denoise failure.
