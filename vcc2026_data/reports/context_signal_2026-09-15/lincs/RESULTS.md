# LINCS genetic context-transfer pilot — completed 2026-09-15

**This pilot does not establish transferable, perturbation-specific destination-context prediction.** A source-only denoiser remains the strongest PDS model. The small additive ridge gain comes entirely from a destination-specific vector shared by every perturbation. Interaction models fail the matched context controls or consistency tests. Expanding the GP response space from 16 to all 924 evaluation genes does not change that conclusion.

The results concern shRNA L1000 signatures and external CCLE baseline expression. They are a separate benchmark from Jiang and the official challenge. Historical b16/Jiang artifacts and scores were preserved; no official reconstruction or upload was performed.

## Data, provenance and endpoint

- **GSE106127 ordinary CGS**, eight lines at 96 hours: A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7 and PC3. VCAP has a different 120-hour endpoint and is excluded.
- Metadata provides 3,931 target symbols at 96 hours: 3,801 occur in at least three and at least five lines; **2,632 occur in all eight**. A fixed target-name SHA256 order selected **200** from that fully shared set before the response matrix was opened. Original source target IDs remain fixed; all 200 resolve to approved HGNC targets. See [coverage summary](coverage_summary.csv), [coverage matrix](../coverage_tables/lincs_shrna_cgs_96.csv) and [selection](panel_selection.json).
- The downloaded CGS and Level4 files contain **978 directly measured landmarks**. All 54 landmark coordinates corresponding to the 200 screened target genes are excluded from response features and evaluation, leaving **924**. The complete selected 978-gene signatures and Level4 instances are retained separately. Inferred genes and globally PC-corrected PRIME signatures are not used. Symbol mappings and exclusions are explicit in [landmark_gene_mapping.csv](landmark_gene_mapping.csv) and [target_gene_mapping.json](target_gene_mapping.json).
- Context is **DepMap/CCLE 23Q4 log2(TPM+1) baseline RNA-seq**, uniquely matched through the local publisher `Model.csv`. The feature axis has 19,191 columns after collision-safe alias joins and duplicate-symbol means. Missing values are retained until training-only imputation. This is an **external untreated assay**, not a culture- and batch-matched LINCS control. Context PCA fits source lines only: six independent components in an outer fold, five in an inner fold. A single training-derived scale preserves expression distances; full-rank per-PC whitening would artificially equalize training distances.
- CGS is a consensus over shRNA signatures in **plate-population z-score space**, not log1p(CP10K) perturbed-minus-control expression. This endpoint mismatch matters when interpreting agreement and future cross-modality transfer. The genetic subset reuses experiments from larger LINCS releases and is not an independent cohort. [Publisher data levels](https://clue.io/connectopedia/data_levels), [CGS construction](https://clue.io/connectopedia/what_is_a_cgs), [GSE106127 files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106127/suppl/).

Acquisition totaled **1,706,623,288 bytes**: 141,593,307-byte CGS; 1,106,966,866-byte ordinary Level4; 449,782,440-byte CCLE expression; and small metadata/checksum files. Publisher SHA512 checksums were verified for LINCS and the publisher MD5 for CCLE; SHA256s, URLs and retrieval dates are in [acquisition_manifest.json](acquisition_manifest.json). Public NIH-funded LINCS data have no redistribution restriction; CCLE is the cited CC BY 4.0 DepMap release. [LINCS policy](https://clue.io/connectopedia/data_redistribution), [DepMap 23Q4 release](https://doi.org/10.25452/figshare.plus.24667905.v2).

## Frozen experiment

Every outer test holds out an entire cell line, while perturbation identities are observed in source lines. Features, centering, response PCA, context PCA, imputation, denoiser fitting and hyperparameter selection use permitted source lines only. No line IDs are model inputs.

The four primary models use 16-dimensional response features and predict a correction to source mean transfer. The simple denoiser is a separately reported source-only anchor, using the unchanged rank128/alpha1 operator from the Jiang engine, refitted on allowed LINCS sources. It is not the official b16 checkpoint applied to a different gene axis.

Each family receives three regularization settings and residual strengths 0/.25/.5/1, selected by inner whole-line PDS. Comparisons include matched no-context, independently selected no-context, context-only, every available wrong source context, zero context, five training-context derangements, five training-label derangements and 32 prediction-target derangements. The learned metric also has two extra fixed-setting seeds. Closed-form deterministic models do not gain independent robustness evidence by repeating an identical seed.

[pilot_spec.json](pilot_spec.json) freezes code/data hashes and settings. Fresh fitting with a poisoned outer response matrix leaves **all 84 inner scores per family, selected settings, unshrunk predictions and the denoiser bit-identical**. The [synthetic preflight](preflight.json) recovers an implanted held-out bilinear context rule, checks context-free query invariance and verifies exact agreement with the existing scorer. The primary run completed **32 outer model-folds in 34.0 seconds**, including the poisoning audit. Cached inner fits make individual family runtimes unsuitable as a hardware benchmark.

## Primary results

Numbers are macro-averages over eight held-out lines. Wrong context averages all permitted wrong query contexts. GP parameter counts include data-dependent posterior coefficients; only two kernel scalars are shared hyperparameters.

| Model | Params | PDS | Matched no-context | Wrong context | Zero context | Correct − no-context | Line wins |
|---|---:|---:|---:|---:|---:|---:|---:|
| Mean transfer | 0 | 0.920107 | — | — | — | — | — |
| Source-only denoiser | rank128 operator | **0.953222** | — | — | — | — | — |
| Additive ridge | 352 | 0.937073 | 0.936106 | 0.932423 | 0.936106 | +0.000967 | 5/8 |
| Bilinear | 1,888 | 0.935176 | 0.936106 | 0.925230 | 0.934384 | −0.000930 | 4/8 |
| Learned metric | 92 | 0.919786 | 0.920107 | 0.913108 | 0.921746 | −0.000320 | 5/8 |
| Shared GP, rank16 response | 22,402 | 0.920091 | 0.920107 | 0.917089 | 0.920383 | −0.000016 | 5/8 |

The line-bootstrap 95% intervals for every correct-minus-no-context gain span zero. For ridge the interval is approximately [−0.000983, +0.002842]. The learned metric's seed standard deviation, 0.000765, exceeds the magnitude of its negative primary gain. Prediction-target shuffle PDS is approximately 0.5 in this raw pseudobulk rank statistic; this is distinct from reconstructed-cell official scores.

| Held-out line | Mean | Denoiser | Ridge | Bilinear | Metric | GP |
|---|---:|---:|---:|---:|---:|---:|
| A375 | 0.943894 | 0.968769 | 0.953040 | 0.955829 | 0.944372 | 0.944271 |
| A549 | 0.935503 | 0.975477 | 0.951935 | 0.951910 | 0.930729 | 0.930905 |
| HA1E | 0.913090 | 0.949146 | 0.929899 | 0.919573 | 0.905854 | 0.910729 |
| HCC515 | 0.819045 | 0.868894 | 0.842965 | 0.835477 | 0.821307 | 0.820553 |
| HEPG2 | 0.936709 | 0.971884 | 0.963342 | 0.957864 | 0.936332 | 0.938693 |
| HT29 | 0.910302 | 0.943492 | 0.927538 | 0.929698 | 0.913668 | 0.912010 |
| MCF7 | 0.956985 | 0.975754 | 0.971357 | 0.974523 | 0.959975 | 0.956985 |
| PC3 | 0.945327 | 0.972362 | 0.956508 | 0.956533 | 0.946055 | 0.946583 |

All models lose to the source-only denoiser on every line. The denoiser improves PDS but reduces mean response cosine from 0.368928 to 0.300605 and increases RMSE from 0.702789 to 0.732351. It is a discrimination improvement, not a general expression-error improvement. Cosine, Pearson, RMSE, MAE, magnitude error, deviation R² and every null are in [the complete CSV](lincs_genetic_transfer_pilot.csv), also copied to the requested [top-level deliverable](../lincs_genetic_transfer_pilot.csv).

**Why ridge is not the desired context result:** its correct-minus-no-context prediction is exactly a target-independent vector in every line. Removing the mean across perturbation targets gives identical centered PDS, **0.951743**, and identical perturbation-specific deviation R², **0.014769**, for correct and no-context models. Thus the extra PDS gain does not encode a perturbation × destination interaction. The numerical energy check is saved in [context_information_checks.csv](context_information_checks.csv). Bilinear adds the interaction but lowers centered PDS by 0.002597 relative to its matched no-context model.

## Replication and context-signal headroom

The primary paper's methods classify the shRNA triplicates as **technical replicates**. Therefore neither collapsed CGS nor X-coded Level4 groups establish biological replication. Reproducible shRNA seed effects are a documented confound. [Primary study, Data curation and Results](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2003213).

Metadata-only filtering gives **82 targets** with identical hairpin reagent sets in X2 and X3 across all eight lines, and **16** with all X1/X2/X3 groups. Technical reruns are averaged within reagent/group, then the same reagents are averaged uniformly in each group. No CGS response-derived weight is reused. All instance, reagent, plate, batch and grouping provenance is retained. These are evaluations of frozen primary predictions: **zero refits or model selection** using the replicate answers.

| Model | X2/X3 mean context gain, 82 targets | Lines winning both groups | Disjoint-reagent mean gain, 65 targets | Lines winning both halves |
|---|---:|---:|---:|---:|
| Ridge | +0.001270 | 5/8 | +0.000225 | 4/8 |
| Bilinear | −0.000847 | 2/8 | −0.002990 | 1/8 |
| Metric | +0.000499 | 4/8 | −0.000240 | 2/8 |
| GP | −0.000047 | 4/8 | +0.000180 | 3/8 |

The stronger reagent check uses **65 targets with distinct seed-family reagents shared across all eight lines**. Target-specific hash order assigns seed families to disjoint halves before their values are read. Every reagent is first averaged over X2/X3. Ridge, metric and GP gains have opposite signs in the two reagent halves; bilinear loses in both. Different reagents are not renamed biological replicates.

For the same 65 targets, with **identical pooled mean and identical denominator**, an approximate cross-half energy decomposition after removing the shared all-target shift is:

| Agreement comparison | Shared across lines | Cell-specific, cross-half repeatable | Half disagreement |
|---|---:|---:|---:|
| Same-reagent technical X2/X3 | 14.98% | 16.38% | 68.64% |
| Distinct seed-family reagent halves | 2.52% | 7.98% | 89.49% |

These are signed cross-product estimators under independent residual-error assumptions; negative estimates would be retained. The identities sum exactly to measured pooled-mean energy. **The reagent disagreement is not pure measurement noise**: it contains differential reagent activity, off-target effects, sampling error and other heterogeneity. Shared assay/plate effects may remain in repeatable components. The remaining cell-specific cross-reagent component is real headroom for further investigation, but its predictability from the external untreated assay was not established here. See [matched decomposition](matched_reagent_technical_comparison/variance_components.csv) and [reagent check](disjoint_reagent_check/replicate_consistency.csv).

![LINCS context controls and repeatability](figures/lincs_context_diagnostics.png)

## Full-response GP sensitivity and next decision

A separately prefrozen sensitivity removed the primary 16D response bottleneck, applying the shared context kernel to **all 924 response genes**. It repeats nested selection, outer-response poisoning and all wrong/zero/training shuffle controls. PDS is **0.920345**, versus no-context **0.920107** and zero-context **0.920320**: gains of only +0.000239 and +0.000025. Centered PDS falls by 0.000166, and perturbation-specific deviation R² is −0.001707. Five of eight lines win raw PDS. This does not support a context success claim or identify rank16 as the principal explanation. The sensitivity took 6.0 seconds and preserves the primary run. [Specification and controls](full_response_gp/spec.json), [scores](full_response_gp/fold_metrics.csv).

LINCS adds many shared targets but only **eight** matched contexts, compared with Jiang's six. It is not yet the many-dozen-context test. Important limits are the external baseline culture mismatch, z-signature endpoint, cancer-heavy line set, reagent heterogeneity, incomplete biological replication and the primary compact feature sweep. The results justify the planned Tahoe pilot with matched controls and many more independent lines; they do not justify a larger Transformer or a drug→gene transfer claim. Biological-replicate and perturbation-specific-context success gates remain unmet.

## Reproduction

New scripts are `scripts/context_signal_lincs_data.py`, `context_signal_lincs_models.py`, `context_signal_lincs.py`, `context_signal_lincs_checks.py`, `context_signal_lincs_replicates.py`, `context_signal_lincs_fullgp.py`, `context_signal_lincs_reagent_check.py` and `context_signal_lincs_report.py`. Use the `ct_transfer` Python environment. Source snapshots for the primary run are under `lincs/sources/`; subsequent sensitivities have separate specs and outputs. Cached complete outer folds are reused only when the primary panel and source hashes match.

The canonical response panel is `processed/context_signal_2026-09-15/lincs_shrna_cgs_panel.npz`. Full 978-coordinate signatures, selected Level4 instances, X-group panels and the disjoint-reagent panel are saved beside it. Downloading the large parent LINCS Level3/4/5 collections was unnecessary.
