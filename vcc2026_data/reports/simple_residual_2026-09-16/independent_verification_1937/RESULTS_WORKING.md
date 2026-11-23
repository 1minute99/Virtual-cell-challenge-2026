# Simple residual experiment — verified results during completion

**Status: official upload and score pending.** Implementation, sweeps, full-population
confirmation, all raw-six/null checks, KOLF and Jiang diagnostics, full-data
challenger training, and official packaging/validation are complete. The selected
predictor is preserved b16 with zero new context residual; the expanded-data
challenger is archived unused. No new official score is asserted here.

## Current primary result

The source-refitted rank32 denoiser with residual alpha zero has the highest
mean in the completed primary comparison. All candidates use identical outer
biological exclusions, evaluation panels, reconstruction and paired seeds.

| Candidate | Four-line, three-seed reconstructed raw PDS | Difference from refitted baseline |
|---|---:|---:|
| No change | 0.500000000 | -0.402890190 |
| Historical-source mean | 0.763470829 | -0.139419361 |
| Archived b16 | 0.854646971 | -0.048243218 |
| Source-refitted rank32 baseline, alpha0 | 0.902890190 | 0 |
| Genetic residual | 0.897102378 | -0.005787811 |
| Genetic extended-budget control | 0.897290041 | -0.005600149 |
| Genetic joint-budget control | 0.896550539 | -0.006339651 |
| Joint genetic + drug residual | 0.896723337 | -0.006166852 |
| Drug pretraining then genetic fine-tuning | 0.900244333 | -0.002645857 |

The refitted baseline exceeds archived b16 on all four line averages. Its mean
gain is 0.048243218 raw PDS. The denoiser, data pool and output-axis handling
change relative to archived b16; this is a pipeline comparison and does not
isolate the learned operator. The historical mean also uses six older source
stores. A separately frozen same-all-source mean diagnostic addresses
that data-pool difference under disjoint query/reference controls: mean transfer
scores 0.822826087 versus the refitted denoiser’s 0.859225195. The denoiser gains
0.036399108 on average, with 3/4 line and 9/12 paired wins, but loses 0.047711
on RPE1. Source expansion and denoising therefore have different contributions
by biology. This supplement is separate from the shared-control primary table. Every residual increment above uses the **same newly
refitted baseline** as its comparator, so those increments do isolate adding the
residual under the evaluated protocol.

The strongest residual arm is pretraining/fine-tuning, but it remains below
alpha zero. Architecture choices vary by outer line. The selected residual
architectures are linear, bilinear and a small MLP, with residual-based output
bases; no single context-dependent architecture establishes a robust gain.
The linear model's target-independent context contribution cancels under query
target centering and must not be described as learned context transfer.

## Dataset and training scope

The assembled source table contains 116,217 pseudobulk rows on 5,592 common
response coordinates: 103,295 genetic rows over 18,711 targets and 12,922 drug
rows. Genetic data include the existing Orion, CD4, Feng, H1, KOLF, Replogle,
Nadig and four usable Jiang-stimulus assets. The full eligible local Tahoe
acquisition contributes 4,431,233 cells from 157 shards, 48 lines and 185
drug-dose conditions. This is the local acquisition, not the entire remote
Tahoe corpus. No new download occurred.

Response labels are source-native log1p(CP10K) treated-minus-matched-control
effects. Full native library sums precede the common-gene intersection. RN
controls are batch matched for response construction; Tahoe uses matched DMSO;
Jiang retains stimulus and biological-replicate identities. Effective cell
counts are weighting proxies, not gene-wise variance estimates. Unmeasured
response values are excluded through the explicit common axis rather than
treated as observed zeros.

The genetic input uses existing ESM2 sequence embeddings projected to 64
dimensions. Drugs have a separate adapter on hashed SMILES features projected
to 64 dimensions. The small residual models use source-only control PCA followed
by a learned linear encoder, then linear/bilinear/MLP interactions and a rank16
or rank32 decoder. The original baseline remains outside the residual span.

The staged search includes 27 initial configurations and 162 native fits;
900 reconstructed inner scores select 20 outer-specific settings. Full
confirmation uses all eligible source rows, 1,200 main updates and an explicit
coverage pass, with three seeds on each of four outer lines. Drug arms add 400
auxiliary updates and drug coverage; matched controls use additional genetic
updates. Exact main genetic minibatch streams match, but full total FLOPs do
not: drug arms have an extra coverage pass. Every eligible row is visited;
that fact alone is not evidence of convergence.

The search covers context dimensions4/8, output ranks16/32,
effect/residual/control-plus-effect bases, weight decay0.001/0.03/0.1,
residual penalties0.1/0.5/1 and shrinkage0/0.1/0.25/0.5/1. Sharing ablations
compare context only, context plus interaction, and all shared downstream
modules, while gene/drug adapters remain separate. `tables/sharing_ablation.csv`
records24 source-inner comparisons at a matched leading structure and drug
weight. They show many alpha-zero ties and small inner differences, not a
robust downstream sharing advantage.

## Drug attribution and context checks

Architecture-matched drug counterfactuals produce a smaller effect than the
independently selected-arm differences:

| Intervention versus matched genetic control | Mean raw-PDS increment | Line signs |
|---|---:|---|
| Joint drug addition | +0.001067447 | 3 positive, 1 alpha0 tie |
| Drug pretraining | -0.000266630 | 2 positive, 1 negative, 1 alpha0 tie |

Neither primary drug arm exceeds the refitted alpha-zero baseline. The 36
no-context/training-context-shuffle scores also do not establish a consistent
benefit from correct context. These are bounded negative/weak findings for
this preprocessing and model family, not proof that drug transfer or context
learning is impossible.

The separate wholly-unseen100 experiment excludes those genes from every
genetic label fit, baseline, basis and selection input. Its global predictor
and residual variants score approximately0.5 raw PDS. It uses a pilot source
population of1,068 genetic targets rather than the full18,711; it must not be
presented as a full-data unseen-gene verdict. The final official300 setting
has empirical source measurements and is reported separately.

## Evaluation fidelity and supplements

Outer lines are K562, RPE1, HepG2 and Jurkat. Their treated responses are excluded
across genetic and drug modalities; source-only whole-line inner selection is
separate. Training anchors additionally exclude each example's own destination
biology. Actual response/control poisoning and source-axis checks pass.

The main mock uses300 targets, up to200 real cells per target and400 reconstructed
predicted cells, scale3 and the existing stochastic count builder. The installed
cell-eval2 vcc2026 preset uses CP50K pseudobulks and panel-target exclusions.
Count validity and reconstruction changes were checked for all1,620 initial
candidate/alpha/direction evaluations. All 21 physical / 24 logical full raw-six evaluations are complete and use the same installed
scorer with its CPU pdex backend; they cannot reproduce the hidden official
per-context calibration anchors and are not a local official overall score.

Primary mock reconstruction and scoring share control cells. The completed
baseline supplement instead partitions800 stored controls into400 query and400
reference cells, without expression-based selection. The refitted baseline
scores0.859225195 versus b16's0.794240059, with4/4 line and12/12 paired wins.
On the identical predicted counts, scoring against the query-control origin
gives0.902245448 and0.856211446. Thus control-origin choice materially changes
absolute scores. The expanded 108-row candidate comparison is complete and independently checked.
Genetic-only leads alpha zero by just 0.000131921, with two improving lines,
one losing line and one exact tie; joint and pretraining lose by 0.000634522
and 0.000226682. This small ranking reversal does not establish a stable
context benefit. The original primary decision is preserved.
For residual candidates, original untreated context features remain frozen;
that extension tests reconstruction/control-origin sensitivity rather than a
fully independent re-query of the context encoder.

The original target shuffle is an alphabetical cyclic shift. The new32-random-
derangement diagnostic permutes already reconstructed groups and verifies
identity PDS against each full-score receipt. All 21 checks find random means near 0.5 while the cyclic control is around
0.54–0.58; every cyclic value exceeds all 32 random draws for its candidate. Preserve both, label
the structured cyclic null correctly, and do not use it as an unbiased chance
floor. Neither diagnostic changes existing raw scores or model settings. Identity
parity error is at most 2.22e-16 and all counts are valid, with maximum cell
total 113,424. The refit’s seed-0 equal-line random-null-adjusted ratio versus
b16 is 1.141884, compared with 1.176228 using the biased cyclic floor. This
seed-0 audit is separate from the three-seed means above. Raw-six exports retain
all per-metric coverage counts; undefined reach rows are not filled with zero.

Jiang's138 separate-replicate arm folds,46 directions and414 postprocessed
variants are complete. The six-line equal-weight total contrast recovery is
-117.0654% genetic, -116.9946% joint and -116.9872% pretraining, versus
-120.7913% for the refitted alpha-zero auxiliary baseline. Residual increments
are therefore +3.7259/+3.7968/+3.8041 percentage points, positive in5/6 lines,
but total error remains worse than predicting no source-to-destination contrast.
Native PDS also falls:0.702419/0.701969/0.700728 versus0.705701 for alpha zero.
These numbers concern repeatable contrast recovery on a pilot population and
cannot replace count-reconstructed PDS or official scoring. They do not prove
that query context caused the small increment. All denominators are identified;
source-mean contrast recovery numerator is exactly zero.

The KOLF check is complete: archived b16 scores 0.596133846, refit 0.565788230,
historical mean 0.554683135, matched-source mean 0.553415332 and no-change 0.5.
Refit loses to b16 in all three paired seeds by a mean of 0.030345616. This
endpoint has 17,971 measured official genes and 279 actual official targets;
the four RN panels contain zero official targets. The archived operator and
source mean both exclude KOLF; their provenance and predictions were checked
independently. Following the requested priority on VCC-like target/gene coverage,
the frozen official choice is therefore archived b16 with zero new context
residual. This supplemental decision preserves the original RN choice. It does
not pool KOLF into the RN score or establish a universal biology ranking.

## Selection, full-data archive and preservation

The original finish02 promotion omitted the evaluated refitted alpha-zero
baseline. A separate final-decision companion repaired that user-requirement
omission after partial local results had been seen, before the full confirmation
finished. Its timing is recorded; it is not a fully prospective pre-outcome
selection rule. Original fits, results and promotion are preserved. The corrected
primary decision selects the refitted rank32 baseline at alpha zero. The later
closer-panel KOLF check supersedes that choice for the one official upload;
`submission_choice_2026-09-16/spec.json` freezes the b16 fallback and all evidence
hashes before packaging. It is not a new or improved predictor.

The final full-data challenger fit is complete and independently verified. It
visited all 103,295 eligible genetic rows (18,711 targets) and zero drug rows.
The archived residual is a 6,776-parameter MLP with context dimension 4, response
rank 16, hidden width 24, residual basis, weight decay 0.03, residual energy
penalty 0.5, learning rate 0.003, batch size 256, 1,200 main updates and a complete
coverage pass, seed 0. Its source-refitted baseline has rank 32. Alpha zero makes
all A/B/C deltas exactly equal that refitted baseline. All 300 official targets
have empirical source measurements. Saved checkpoint state, source arrays, axes
and prediction hashes passed independent verification. The existing exposure
counter overstates the last short coverage batch by 129 examples; exact unique-row
coverage is correct, and the actual example count is 410,495.

This new archive is not the selected official predictor. The selected historical
b16 operator has rank 128 over 17,742 correction coordinates, trained on 1,661
nonofficial targets from HCT116, HEK293T, CD4 and Feng iPSC sources. Its empirical
base model uses the preserved Orion, CD4, KOLF, Replogle, Feng and VCC2025 source
means. New context-residual alpha is zero; the historical signed covariance
correction remains active. No Tahoe, new Jiang or new RN residual weights enter
the selected prediction. Do not describe b16 as newly retrained on the expanded
103,295-row table.

The fresh official dry-run, full prep and archive validator passed: 360,000 cells,
18,533 genes, 120,000 cells per A/B/C, 400 per target, no controls, no dropped or
reordered coordinates, and counts preserved. Package b18.vcc is byte-identical to
the preserved b16 package (SHA256
010e77fee8e078b58bcc13cec5e4679f8bbc9b35195c5269faf55fa49c3e0abb). Live authorization
and allowance checks passed at 21:11 UTC. Exactly one opaque b18 upload is underway;
its published result is still pending. This is a baseline revalidation, not a
new model or an estimated score improvement.

Independent verification passed124 distinct automated tests, plus source/delta
identity, poisoning and reconstruction checks. The uploader's duplicate-retry
defect was reproduced and corrected, with the original source preserved and26
uploader tests passing. All13 original reference hashes, metadata for11,659
historical research files, the b16 package and all11 assembled table files match
their frozen records. The fresh end-of-training preservation audit also passes all 13 reference hashes,
21 key prior report hashes, metadata for 11,659 prior research files and the
original b16 package hash.

See `../METHODS.md`, `../tables/`, `../matched_drug/`,
`../independent_controls_1944/REPORT.md`,
`../independent_jiang_resume_1959/REPORT.md`, and the two independent review
reports for the exact scopes, scripts, manifests and qualifications.
