# B19 optimization: completed local study

Completed 2026-09-17T08:01:09.194507+00:00. No new official submission.

The frozen outer promotion gate **failed**. The source-selected pipeline changes mean reconstructed PDS by +0.000385, with 2/5 line-average wins and 8/15 paired-seed wins. No deployment refit was triggered; the existing official b19 remains the retained model.

This evaluates the nested selection procedure. Its per-line choices differ; these outer scores are not a direct validation of one universal final configuration.

## Did b19 need more epochs?

| Schedule | Sampled updates | Weighted source latent MSE | Mean inner reconstructed PDS |
| --- | --- | --- | --- |
| constant | 2500 | 0.416567 | 0.772036 |
| constant | 10000 | 0.368568 | 0.770536 |
| constant | 25000 | 0.348465 | 0.767324 |
| cosine | 2500 | 0.416515 | 0.772073 |
| cosine | 10000 | 0.365971 | 0.770536 |
| cosine | 25000 | 0.344988 | 0.767051 |

Longer training substantially reduces the fitted source residual error but reduces transfer PDS on average. Therefore 2,500 updates were not full source-loss convergence, yet simply training longer is not a reliable improvement. Three outer selections chose 2,500 updates, two chose 10,000, and none chose 25,000.

An exhaustive coverage pass precedes weighted sampling. “Updates” are minibatches of 256, not uniform epochs. Historical full-data b19 saw 743,295 examples over 103,295 rows (7.196 dataset-equivalent exposures), with only the first pass guaranteeing every row. Cosine scheduling used a fixed 25,000-update horizon and learning rate 0.003 to 0.00003.

![Convergence diagnostic](convergence.png)

## Architecture and training changes

Prediction retains the empirical rank-32 source denoiser plus beta × global residual and gamma × context residual. Both branches already exist in b19. This study separates their output penalties and shrinkage, and compares joint training with global-first training followed by a frozen-global context phase. It does not introduce a generative model.

Every staged phase covers all permitted rows. Total optimizer updates match the joint budget; the second coverage pass replaces some weighted sampled updates. Global penalty is 0.5; separate context penalties are 0.5 or 2. Weight decay is 0.03. Source-inner native PDS shortlists branch strengths, then reconstructed-count PDS selects them. Zero residual, global-only, and context-only candidates are included.

| Outer line | Training | Updates | Schedule | beta | gamma | Rank | Hidden |
| --- | --- | --- | --- | --- | --- | --- | --- |
| K562 | joint_total | 2500 | cosine | 1.0 | 1.0 | 32 | 64 |
| RPE1 | staged | 2500 | cosine | 0.5 | 0.25 | 32 | 64 |
| HepG2 | staged | 10000 | constant | 1.0 | 0.0 | 32 | 64 |
| Jurkat | staged | 10000 | cosine | 0.5 | 0.5 | 32 | 64 |
| KOLF | joint_split | 2500 | constant | 0.5 | 0.0 | 32 | 64 |

All five selections retain rank 32 / hidden 64 (24,192 parameters). Rank 64 (30,336 parameters) was tested in 20 directed inner folds and passed none of the expansion gates, so rank 128 was not run. The targeted hidden-128 model (72,448 parameters) ran in four Jurkat-excluded inner folds and lost selection. Bigger capacity did not help in this bounded study.

## Three-seed outer comparison

| Candidate | Mean | K562 | RPE1 | HepG2 | Jurkat | KOLF |
| --- | --- | --- | --- | --- | --- | --- |
| no_change | 0.500000 | 0.500000 | 0.500000 | 0.500000 | 0.500000 | 0.500000 |
| mean | 0.736193 | 0.760364 | 0.807291 | 0.747179 | 0.811446 | 0.554683 |
| b16 | 0.754619 | 0.786061 | 0.741880 | 0.813597 | 0.835422 | 0.596134 |
| refit | 0.800538 | 0.882070 | 0.810684 | 0.884560 | 0.859588 | 0.565788 |
| b19_reference | 0.795777 | 0.883530 | 0.810312 | 0.856429 | 0.859666 | 0.568947 |
| duration | 0.794996 | 0.883501 | 0.810290 | 0.855284 | 0.856957 | 0.568947 |
| branches | 0.796162 | 0.883501 | 0.814902 | 0.855831 | 0.859718 | 0.566858 |
| selected | 0.796162 | 0.883501 | 0.814902 | 0.855831 | 0.859718 | 0.566858 |

All values above are local raw PDS averaged over seeds 0/1/2. `b19_reference` refits the fixed b19 recipe without each outer biology; it is not the all-source submitted checkpoint applied to its training data. `b16` uses the eligible archived operator on the same permitted-source mean. `refit` is the rank-32 baseline with zero new residual. Branch and final selected rows coincide because no larger model won.

| Candidate | PDS minus shuffled null | Relative excess vs fixed b19 |
| --- | --- | --- |
| mean | 0.238534 | 0.808247 |
| b16 | 0.255560 | 0.945192 |
| refit | 0.300235 | 1.006815 |
| b19_reference | 0.295204 | 1.000000 |
| duration | 0.294424 | 0.997811 |
| selected | 0.295534 | 0.997363 |

For each paired line/seed, the last column is (candidate PDS − its 32-derangement mean)/(fixed b19 PDS − that same candidate null). Values are then averaged; 1 means the fixed b19 discrimination level. See null_adjusted_comparison.csv for every paired calculation.

![Paired-seed outer gains](outer_gains.png)

The promotion gate requires mean gain ≥0.001, at least 3/5 line wins and 9/15 paired wins, with all predictions valid. Seeds are optimization/reconstruction repetitions, not additional biological replications.

## Does context contribute?

| Line | Correct | Trained no-context | Training context shuffle | Zero query | Wrong query |
| --- | --- | --- | --- | --- | --- |
| K562 | 0.888372 | 0.883333 | 0.766109 | 0.884582 | 0.837659 |
| RPE1 | 0.813679 | 0.813166 | 0.812843 | 0.813166 | 0.814247 |
| HepG2 | 0.858874 | 0.858874 | 0.858874 | 0.858874 | 0.858874 |
| Jurkat | 0.855897 | 0.853868 | 0.855987 | 0.853868 | 0.856566 |
| KOLF | 0.566515 | 0.567585 | 0.569016 | 0.566515 | 0.566515 |

These are matched seed-0 controls at the selected configuration, not independently retuned alternatives. Gamma zero explicitly removes context at inference. A benefit from changing training or denoising is not automatically evidence for useful context transfer. K562 is an exploratory positive: correct context beats all four controls at seed 0, including a +0.005039 PDS gain over the independently trained no-context model. RPE1 and Jurkat instead lose to wrong-query controls, and HepG2/KOLF select gamma zero. This is localized evidence, not proof that context is absent or a consistent cross-line context benefit.

## Installed full six-metric diagnostics

Fixed b19-reference and selected predictions were reconstructed and run through installed `cell-eval2 --preset vcc2026`, seed 0, 400 cells per target, scale 3. Query controls are used for prediction; disjoint reference controls replace real controls. Reconstructed pseudobulk and full-scorer PDS were checked against the primary path.

| Line | Candidate | PDS ↑ | MSE ↓ | LFC NMAE ↓ | Fidelity ↑ | Reach ↑ | Jaccard ↑ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| K562 | b19_reference | 0.888439 | 5.195094 | 0.803670 | 0.668053 | 0.615208 | 0.028327 |
| K562 | selected | 0.888372 | 5.197508 | 0.806897 | 0.668034 | 0.614927 | 0.028337 |
| RPE1 | b19_reference | 0.810981 | 1.475666 | 0.849752 | 0.614628 | 0.327468 | 0.086334 |
| RPE1 | selected | 0.813679 | 1.458953 | 0.846249 | 0.616802 | 0.312324 | 0.081953 |
| HepG2 | b19_reference | 0.851048 | 2.679304 | 0.820947 | 0.607954 | 0.628601 | 0.057864 |
| HepG2 | selected | 0.858874 | 2.213679 | 0.799660 | 0.620879 | 0.600178 | 0.054415 |
| Jurkat | b19_reference | 0.861661 | 3.950383 | 0.757775 | 0.650139 | 0.598686 | 0.034249 |
| Jurkat | selected | 0.855897 | 4.578977 | 0.769361 | 0.636455 | 0.596842 | 0.034023 |
| KOLF | b19_reference | 0.565934 | 5.154076 | 0.963143 | 0.498465 | 0.100879 | 0.012152 |
| KOLF | selected | 0.566515 | 5.179351 | 0.964781 | 0.498946 | 0.113731 | 0.012661 |

These are raw local metrics; do not average them into an official overall score. Error metrics have different directions and calibration from published leaderboard scores. Any missing metric is explicit in the corresponding receipt. Detailed count clipping, count limits and normalization checks are in `full_metrics/*/*/reconstruction.json` and primary score receipts.

## Data, scope and limitations

- All fitting uses the full permitted genetic population in each scope: the existing table has 103,295 genetic rows / 18,711 targets before exclusions. Replogle, Nadig, Jiang and existing Orion/CD4/Feng/KOLF/VCC2025 assets are reused. No Tahoe/drug rows are trained in this follow-up; earlier drug comparisons are unchanged.
- Five whole-line outer holdouts; ten excluded-line-pair source-only inner preparations. Training anchors additionally exclude each training example’s destination biology. Basis fitting, preprocessing statistics and model selection exclude outer responses.
- Inner scoring uses 64 targets per line; outer RN panels use 300 mock targets disjoint from the official 300, while KOLF uses 279 official targets. Native measured gene axes and local target-gene exclusion differ from the complete official 300-target / 18,533-gene panel, so neither raw PDS levels nor auxiliary metrics exactly reproduce official validation. Residuals occupy a source-fitted low-rank subspace of the 5,592 shared learned coordinates; the baseline is preserved elsewhere.
- Control inputs are disjoint from evaluation-reference controls. Existing gene mappings, log1p CP10K training effects, target-centered residual decoding, native-axis CP50K scoring and submission-time stochastic count reconstruction are reused.
- These previously studied local panels are not new independent biological validation. This is measured-source/seen-target transfer; wholly unseen-target and biological-repeat conclusions remain in prior reports.
- The deployment consensus rule was frozen before outer results: most frequent exact source-inner choice, then supporting-fold inner mean and lexical tie-break. Here every exact choice occurs once; the tie-break selects the KOLF-excluded choice. Different inner line difficulties limit this consensus rule.

## Deployment and official records

Consensus configuration: `joint_split_constant_s2500_r32_h64_c0.5_b0.5_g0`. See `final/decision.json` for the recorded decision not to refit/promote.

| Official entry | Submission ID | Published overall | Published PDS |
| --- | --- | --- | --- |
| b16 | preserved historical record | 0.12919304282260793 | 0.554670505212238 |
| b19 | pamfw9b1glTBvlXvEyQA | 0.13646991123621474 | 0.5431564640579106 |
| This follow-up | none | not submitted | not submitted |

B19 checkpoint/package remain in `../context_forcing_strength_2026-09-16/final/model.pt` and `../context_forcing_strength_2026-09-16/submission/b19.vcc`. Official receipt: `../context_forcing_strength_2026-09-16/submission/receipt.json`. No score is estimated or presented as official.

## Reproducibility and verification

169 unique neural cross-validation fits / 209 saved checkpoints; 120 convergence scores, 1,200 native branch-screen scores, 120 outer score rows and 20 context controls. See CSVs and immutable per-fold requests, model receipts and score JSONs. The final audit checked 1292 files and preserved all 20 frozen historical references.

Fourteen core pytest tests pass, plus the exact bulk/disjoint-control test and an actual-score binding/corruption check. Independent reviews are saved. A stale frozen-gradient clipping issue was corrected before any staged scientific fit; the original and amended source hashes and 30 then-existing joint checkpoints are preserved in `pre_staged_gradient_amendment/`. Joint training remains bit-exact with the original trainer in the regression test. Later receipt checks were tightened before full-metric scientific runs.

No fit, local score job, upload or public score poll is pending. See `CLOSEOUT_VERIFICATION.json`, `research_status.json`, `DELIVERABLES.md`, `protocol.json`, and `finish_sources_frozen.json`.
