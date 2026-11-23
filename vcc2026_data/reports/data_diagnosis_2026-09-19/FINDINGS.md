# Second data diagnosis: is the data okay to use? — 2026-09-19

The first audit asked whether the data was PROCESSED correctly. This one asks whether the
measurements themselves are good enough to learn transfer from. Different question, different
answer. All figures are target-specific correlations: each source is centred over its own
targets first, so the shared generic response cannot inflate them.

## 1. Processing: verified correct, repeatedly

- Oracle test: real measured effects fed through our own reconstruction score pds 1.000 and beat
  a replicate on 5 of 6 metrics. No bug in reconstruction, gene mapping, target alignment or
  normalization.
- Wide-axis rebuild reproduced the frozen table bit-identically on all 5,592 shared genes.
- The gene-subset run reproduced the frozen refit exactly (0.8005 five-line, K562 0.8821).
- Official/mock target coordinates (600) are excluded from training everywhere.

Nothing found here is a preprocessing defect.

## 2. We cannot measure our own noise — only 2.4% of rows have a replicate

| source | true same-condition replicate? |
|---|---|
| Jiang (2,481 rows) | yes, Rep1/Rep2 |
| KOLF ALPHA/BETA/GAMMA | **no** — the three batches screen DISJOINT gene sets (0, 0 and 1 shared targets) |
| Orion, CD4, Feng, VCC2025, Replogle-Nadig | no |

So for 97.6% of the table there is no way to separate measurement noise from real biological
difference. That is the single biggest blind spot in the dataset.

## 3. Where reproducibility IS measurable, it is low

Jiang Rep1 vs Rep2, same line, same stimulus: **median r = 0.209** (range 0.04-0.41).
Roughly 80% of that signal is not reproducible. Jiang also has the lowest cell counts, so this
understates reliability elsewhere — but it is the only direct measurement we have.

## 4. Cross-biology agreement, and what modulates it

| comparison | n targets | r | mean abs effect |
|---|---|---|---|
| CD4 Rest vs Stim8hr — the SAME cells, different state | 11,873 | **0.270** | — |
| CD4 Rest vs Stim48hr | 11,817 | 0.187 | — |
| RN K562 vs RPE1 (same source) | 2,089 | 0.163 | 0.099 |
| RN HepG2 vs Jurkat (same source) | 1,901 | 0.174 | 0.099 |
| Orion HCT116 vs HEK293T (same source) | 17,498 | 0.063 | 0.058 |
| Orion HCT116 vs RN K562 (cross source) | 9,123 | 0.071 | 0.063 |
| CD4 Rest vs RN K562 (cross source) | 8,456 | 0.030 | 0.063 |
| Orion vs Feng / CD4 / KOLF (cross source) | 1,795-11,434 | 0.013-0.017 | — |

**The decisive number is the first row.** The same CD4 cells, same donors, same protocol, merely
8 hours of stimulation apart, agree at only 0.270 — barely above Jiang's technical-replicate
ceiling of 0.209. If a state change in identical cells nearly exhausts reproducibility, then
cross-CELL-TYPE agreement of 0.03-0.07 is close to the biological ceiling, not a modelling
failure.

RN pairs agree better (0.163-0.174) than Orion pairs (0.063). Restricting Orion to RN-screened
targets moves it only 0.063 -> 0.085, so target composition explains part of it. The rest is
amplitude: on the SAME targets RN reports effects 1.6x larger than Orion (0.099 vs 0.061), i.e.
the sources are not on a common scale. (Per-source amplitude calibration was tried previously and
scored -0.003 pds, so this is a known, already-rejected lever.)

## Verdict

**The data is safe to use and is not mispreprocessed.** It is also not short of rows: 103,295
genetic rows over 18,711 targets.

**But it is close to exhausted as a source of cross-biology signal.** The transferable component
is r ~ 0.02-0.07, the reproducibility ceiling where we can measure it is ~0.21, and the same cells
in two states only reach 0.27. Every method tried against this — longer training, learning rates,
capacity, branch structure, seed ensembling, the full gene axis, reliability-ranked gene subsets —
has returned between -0.07 and +0.001, with the gains concentrated on lines that resemble the
training sources and absent on KOLF.

**The honest bottleneck is not the pipeline and not the gene axis. It is that we have ~8 biological
contexts, no replicates to calibrate them, and a signal that barely survives a change of cell
state — let alone a change of cell type.**

### What would actually move this
1. More independent biological contexts, especially ones resembling the scored contexts
   (T-cell and keratinocyte). This is a data-acquisition problem, not an engineering one.
2. Replicate structure for the large sources, so noise can be separated from biology and the
   models can be weighted by real reliability rather than cell counts.
3. Noise-aware shrinkage rather than more capacity: with r ~ 0.05, the optimal predictor is
   heavily shrunk, which is consistent with the measured MSE-optimal scale of 0.10-0.64 while we
   submit at 3.
