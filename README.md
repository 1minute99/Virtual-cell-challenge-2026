# Virtual Cell Challenge 2026 — participation record

This repository is our record of taking part in the Arc Institute **Virtual Cell Challenge 2026**: predicting the
transcriptional response to CRISPRi knockdowns in cell lines a model has never seen (zero-shot), scored by
`cell-eval2 --preset vcc2026`. It documents **what we tried, why most of it did not work, and which methods actually
raised the score** — with the original lab notebooks, the code of the final system, and the follow-up analyses that
tested our explanations.

> Private record. Large data, official challenge data and model checkpoints are not included. Other participants'
> names in the archived notebooks are replaced by "[another team]".

## At a glance

| | |
|---|---|
| Official submissions | 38 (36 in the validation phase, 2 in the final phase) |
| Best validation score | **0.2545** (S18) |
| Standing at the end of validation | **129 of 1,635 teams (top 7.9 %)**; gap to the top-10 threshold 0.0649 |
| Best standing during validation | 34 of 1,244 (top 2.7 %, 2026-09-30) |
| Final-phase entry | F2 (uploaded 2026-10-28); final scores not yet published |

Details: [`docs/02-submission-history.md`](docs/02-submission-history.md).

**The final system** (F2) predicted each knockdown's response in a new cell line by **transferring measured responses
from public Perturb-seq screens of other lines**, not by a learned cell model:

- **M1 — multi-source fold-change fusion.** Each public screen gave a fold change for the target. These were fused
  across sources and applied to the destination's own control cells. A small learned part ("S4") and a co-essentiality
  neighbour blend were mixed in.
- **M2 — agreement amplitude.** Targets on which the sources agreed got larger predicted responses (rank^0.5).
- **M3 — control-PC deflation.** For the weaker half of targets, the part of the prediction lying along the
  destination's own top control-cell axes was shrunk (k = 8, η = 0.75).
- **MSE budget at the "kink".** The total predicted response energy was matched to the point where the scorer's MSE
  member stops rewarding more amplitude.
- **Emission and DE side.** Count-level emission pinned the group sums and created the differential-expression
  signal inside cells (dual-moment emission, cpc composition pin).

Details: [`docs/01-challenge-and-final-system.md`](docs/01-challenge-and-final-system.md).

## What raised the score

The largest official gains, by size:

| Change | Submissions | Official Δ | Attribution |
|---|---|---|---|
| Six-source mean + stochastic count rounding | B2 → B6 | +0.0750 | bundle of 2 |
| Response scaled ×3 | B6 → B8 | +0.0447 | single (notes) |
| Dual-moment emission (sums pinned, DE created within cells) | S8 → S10b | +0.0348 | bundle of 2 |
| cpc composition pin (DE members) | S10b → S10c | +0.0247 | single |
| **M2** agreement amplitude | S14fk2 → S15 | +0.0106 | single |

Twenty changes raised the official score; the full ranked list, including M3 (which could only be measured inside a
bundle, S16 → S18 +0.0070) and the small ride-alongs, is in
[`docs/03-what-raised-the-score.md`](docs/03-what-raised-the-score.md). The single largest *local* lever — a same-line
screen from another lab (+0.069 composite on an H1 panel, about nine times the largest modelling gain) — never reached an
official entry, because no such screen existed for the validation or final cell lines.

## What did not work

- **Learned architectures.** All 17 architectures in a joint screen lost 0.04–0.22 composite to the transfer baseline,
  mostly through PDS. Rerunning six of them with two more seeds confirmed it: none beats the baseline on the two-line
  mean. The best one helped K562 and hurt HepG2.
- **Context-conditioned models.** Training on more contexts did not help (flat from 6 to 35 contexts). On drug data,
  models told contexts apart long before that made their predictions any better.
- **Foundation models and embeddings.** STATE, scLong, UCE, ESM-2 and PIE gave no usable gain; random or zero gene
  features matched or beat the embeddings.
- **Drug data as a signal.** Tahoe-100M context models gained nothing even on drug data, the one valid drug→gene test
  was null, and LINCS sources were rejected. After correcting for reliability, the drug–CRISPRi concordance explains
  about 0.3 % of the reproducible response.
- **Other sources and learned weighting.** VIPerturb-seq, DepMap beyond the neighbour blend, and response-space kNN did
  not help, and learned source weighting failed six times.
- **Amplitude and metric tricks.** Amplitude reshaping beyond M2, and moving the budget off the kink, traded PDS for MSE
  or lost officially.
- **The final "big-lever" search.** Every remaining public lever screened before the final phase was worth ≤ +0.005.

Catalogue with numbers and evidence: [`docs/04-what-did-not-work.md`](docs/04-what-did-not-work.md).

## Why

- **The noise ceiling is low.** On average only about 28 % of a line's deviation from the cross-line mean was
  reproducible (0.7 % on KOLF). The same target's response agreed poorly across cell lines.
- **The interaction cannot be transferred.** The part of the response specific to a line × target combination cannot
  be learned zero-shot. Weak responders, which decide PDS, are exactly the targets that transfer worst.
- **Too few contexts.** A learned context correction needs many distinct cell contexts. In a follow-up experiment on
  48 drug-treated lines it was harmful below about 24 contexts and broke even only near 37, while public genetic
  Perturb-seq covers about 11 distinct cell types.
- **The scorer reads sums.** PDS and MSE depend only on each target's summed counts and are decided in different
  target populations (MSE in the strongest responders, PDS in the weakest). One reshaping lever cannot raise both, and
  the MSE member has a kink set by the test data's noise.
- **Local gains shrank officially.** The one valid learned-context case shrank about 19-fold, from +0.0247 locally to
  +0.0013 officially.

Explanations, scorer mechanics and lessons: [`docs/05-why-it-does-not-work.md`](docs/05-why-it-does-not-work.md).

## Follow-up analyses

After the challenge we tested these explanations with a pre-registered set of analyses (plan:
`input/ideas/zero-shot-limits-context-complexity-plan.md`; summary: `knowledge/synthesis/results-summary.md`). Their
outcomes:

| Question | Outcome |
|---|---|
| Does a learned context correction need many contexts? (48 Tahoe drug lines) | yes — harmful at 12, break-even ≈ 37 [33.3, > 40] |
| Does a linear-Gaussian formula predict that break-even? | no — ×6 too pessimistic (230 vs 37) |
| Are public genetic data below the break-even? | for distinct cell types (11) yes; counting donor lines (25) no |
| Does transfer error sit on the destination's own control axes? | yes, consistently across panels, but small on H1 |
| Does M3 beat five alternative corrections? | inconclusive (wins on the 7-panel mean, loses on one panel) |
| Is one same-line screen worth more than modelling? | yes — about 9× the largest modelling gain |
| How close is the pipeline to the zero-shot ceiling? | untestable as registered; exploratory ≈ 0.77 of a bounded ceiling |

The figures are `output/article/figures/fig1.pdf` … `fig6.pdf`, with captions in `CAPTIONS.md`. The analysis code is
`output/data-analysis/` and its tests are `pytest output/data-analysis/tests`.

## Repository layout

```
docs/                    the write-up (start here: 01 → 05)
AGENTS.md                the day-by-day lab notebook kept during the challenge (newest first)
RESULTS.md, VCC2026_gene_perturbation_report.md   early results and the public-data survey
vcc2026_data/reports/    archived notebooks: PROJECT_HISTORY.md (Sept 5–20), FINDINGS/README files per
                         experiment, official_scores.csv, F1/F2 upload receipts (list: docs/archive_files.txt)
vcc2026_data/scripts/    the challenge pipeline (modules used by the final system and the analyses)
external/era-aivc/       offline evaluator modules from google-research/era (Apache-2.0)
output/data-analysis/    follow-up analysis code, tests and result summaries
output/article/figures/  figures 1–6
output/code/             input manifest with sha256, conda environments, data-restore script
knowledge/, input/ideas/ notes from the follow-up analyses (theory, results summary, analysis plan)
```

## Notes

- **Data.** Inputs are public screens (Replogle 2022, Nadig 2025, X-Atlas/Orion, Tahoe-100M, GSE295214, KOLF2.1J,
  HipSci, CD4 T-cell CRISPRi, DLD-1, DepMap), listed with sha256 in `output/code/INPUTS.md`. Official challenge data
  are not redistributed.
- **Paths.** Scripts use absolute paths of the original workstation.
- **Licence.** Not yet chosen. `external/era-aivc` keeps its Apache-2.0 licence.
