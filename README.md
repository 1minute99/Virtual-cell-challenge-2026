# Data-limited, not model-limited — zero-shot perturbation prediction in unseen cell lines

Research code, results and the pre-registered plan behind a preprint in preparation (target: arXiv cs.LG,
late November 2026). The work grew out of our participation in the Arc Institute **Virtual Cell Challenge 2026**
(zero-shot CRISPRi response prediction in unseen cell lines).

**Question.** When responses to CRISPRi knockdowns in an unseen cell line are predicted from public screens in other
lines, is accuracy limited by model capacity, or by the data — the measurement-noise ceiling, the line × perturbation
interaction that cannot be transferred, and the number of cell contexts available to learn from?

> Status: analysis complete; manuscript being drafted. This repository is private until the preprint is posted.

## Pre-registration

The hypotheses, falsification criteria, stop criterion and task list were frozen before the confirmatory analyses ran:

- Design: `input/ideas/zero-shot-limits-context-complexity-design.md` (approved 2026-11-05)
- Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` — signed off 2026-11-05
  (sha256 at sign-off `59503cbae05f2c5a210a7ba8bbf4f64a38ce7578075ab8012486ebc2dc9ccf18`); three pre-data amendments
  (A1–A3) appended the same day (current sha256 `7adfc134514eea6ba9add95379f6da9e0b7ce62064e9eb5079c118b33e87114d`).
  The original hypothesis text above the amendments is unchanged.
- Every deviation, clarification and independent review is recorded in the project log (kept outside this repository)
  and summarised in `knowledge/synthesis/results-summary.md`.

## Results at a glance

| Hypothesis | Verdict | Key numbers |
|---|---|---|
| H2a — a learned context correction (Tahoe, 48 drug-treated lines) hurts at small C and breaks even at ≥ 24 contexts | **confirmed** | context-pathway gain at C = 12: −4.2e-4 [−5.2, −3.3]e-4; break-even **C\* ≈ 37.4** [33.3, > 40] |
| H2b — a linear-Gaussian break-even formula predicts C\* within ×2 | **rejected** | C\*_theory 230 vs 37.4 (×6 too pessimistic) |
| H2c — public genetic data lie below the break-even | **holds for distinct cell types, not donor lines** | C\*_min 16.45 vs 11 cell types / 25 donor lines |
| H1 — the transfer pipeline reaches ≥ 0.7 of the zero-shot ceiling | **untestable** (registered ceiling internally inconsistent) | exploratory: 0.77 of a bounded ceiling |
| H3a — transfer error concentrates on the destination's own control axes | **confirmed** | D_other 0.018 [0.009, 0.027]; halves at matched n, positive on 7/7 panels |
| H3b — destination control-PC deflation beats five alternative corrections | **inconclusive** | beats all five on the 7-panel mean; loses to a Systema-style correction on the H1 real scorer |
| H4 — one same-line screen is worth more than the largest modelling gain | **confirmed** | +0.069 composite vs +0.0075 (≈ ×9) |

Figures: `output/article/figures/fig1.pdf` … `fig6.pdf` (captions in `output/article/figures/CAPTIONS.md`).

## Repository layout

```
input/ideas/            design and pre-registered plan (earlier, superseded design/plan kept for the record)
input/notes/            external gates (challenge rules, arXiv endorsement), Molina & Zhang overlap note
input/bibliography/     literature guide and acquisition worklist (no PDFs)
knowledge/              research wiki: theory, positioning, results summary, per-task assumption notes, source pages
output/data-analysis/   analysis code (one script per plan task), tests, figure code, top-level result summaries
output/article/         figures and captions (manuscript to follow)
output/code/            input manifest (INPUTS.md, inputs.sha256), import-closure tool and list, data-restore script, conda envs
output/bibtex/          references.bib (entries verified against Crossref/arXiv)
vcc2026_data/scripts/   the challenge pipeline modules imported by the analysis (import closure only)
external/era-aivc/      evaluator modules from google-research/era (Apache-2.0, LICENSE included)
scripts/                job runner (run_detached.sh), wiki lint and graph tools
```

## Reproducing

1. **Environments:** `output/code/envs/vcc2026.yml` (pipeline, scorer) and `output/code/envs/ct_transfer.yml`
   (Tahoe/genetic extraction, architecture screen). The official scorer is `cell-eval2 --preset vcc2026`
   (`--set de.backend=pdex`).
2. **Inputs:** all inputs are listed with sha256 in `output/code/INPUTS.md` / `inputs.sha256` (440 files, ~55 GB).
   They are public datasets (Replogle 2022, Nadig 2025, X-Atlas/Orion, Tahoe-100M, GSE295214 H1 CRISPRi, KOLF2.1J,
   HipSci/Feng, CD4 T-cell CRISPRi, DLD-1) plus derived tables. The Tahoe cell-level shards are restored with
   `output/code/restore_tahoe_shards.py` (immutable revision `2dc57900…`, sha256-checked).
   **Official Virtual Cell Challenge data are not redistributed**; only derived statistics are reported.
3. **Run order** follows the plan's task numbers (`output/data-analysis/`): inputs `task04_freeze_inputs.py`;
   RQ2 `rq2_tahoe_curve.py` → `rq2_h2a.py`, `rq2_params.py`, `rq2_gene_*.py`, `rq2_arch_replicate.sh`,
   `rq2_source_curve*.py`; RQ1 `rq1_components.py` → `rq1_ceiling.py`; RQ3 `rq3_build_preds.py` → `rq3_energy.py` /
   `rq3_h3a.py` → `rq3_corrections.py` (+ matched-n sensitivities); RQ4 `rq4_sameline.py`; record `e4_official.py`;
   appendix `appx_drug_gene.py`; figures `figures.py`. Tests: `pytest output/data-analysis/tests`.
4. **Compute:** CPU only; real-scorer runs take ~25 min and ~15 GB each and must run one at a time
   (`scripts/run_detached.sh`, a shared `flock`). Paths in the scripts are absolute to the original workstation and
   will be made relative for the public release.

## Licence and citation

Licence to be decided before the public release; `external/era-aivc` keeps its Apache-2.0 licence. A citation entry
will be added when the preprint is posted.
