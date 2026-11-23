# 04 — What did not work

This file catalogues what the team tried during the Virtual Cell Challenge 2026 that did not raise the score, with
the measurement behind each verdict. The final system is in `01-challenge-and-final-system.md`, the upload sequence in
`02-submission-history.md`, what did raise the score in `03-what-raised-the-score.md`, and the reasons for the
failures in `05-why-it-does-not-work.md`.

## How to read the tables

Every result says where it was measured (column "By"):

- **O — official leaderboard.** Validation-phase composite on contexts A/B/C (mean of six normalised members). One
  upload per arm, no seed replicates.
- **R — real local scorer.** `cell-eval2 --preset vcc2026` run locally on held-out panels built from public screens,
  with the destination's own data removed from every source. Until 2026-09-22 this was the five-line harness (K562,
  RPE1, HepG2, Jurkat, KOLF); from 2026-09-23 the H1 panel, later joined by KOLF, OR_HCT116(X) and the four RN
  panels. "PDS" means the raw `pds_cosine` member only.
- **X — offline or proxy.** Exact offline sum-side evaluators, pseudobulk PDS on external benchmarks, residual cosine,
  and the screening proxy J = 2.23 ΔPDS − 1.05 ΔMSE (≈ 6 × Δcomposite).

Deltas are composite unless a member is named. On the five-line harness the paired-seed noise floor was ±0.010
composite; from 2026-09-27 no upload was spent on an expected official gain below +0.005, because local gains shrank on
the leaderboard (`05-why-it-does-not-work.md`, Part B). Submission labels follow
`output/data-analysis/results/official_record.csv`.

Sources: `AGENTS.md`; `vcc2026_data/reports/PROJECT_HISTORY.md`.

## (a) Learned architectures

| What was tried | Result | By |
|---|---|---|
| Seven neural nets on the KOLF harness (Sept 13–14): raw-delta and residual MLPs, cosine term, magnitude weighting, 58 contexts, learned gate, DeepSets set predictor | None transferred; the raw-delta net collapsed to ~0, the gate learned multipliers of 0.97–1.01. The set predictor's cosine gain was a shared all-target component (42.5 % of its output energy vs 2.6 % for mean transfer): PDS 0.0152 vs 0.1570, composite 0.1333 vs 0.1972. | R |
| Residual strength on the refit baseline (b19 line, Sept 16), α 0–4 | Five-line PDS fell with α on four lines (α 0: 0.800538, α 1: 0.794715, α 2: 0.785289); only KOLF rose. The shipped b19 (α 1) gained +0.0073 officially over b18, entirely through NMAE; its PDS was 0.0115 lower. | R, O |
| Training length, capacity and optimiser grids for that residual (169 fits; 30 LR configurations; weight decay × batch) | Longer training lowered source error and inner PDS (0.7721 at 2,500 updates, 0.7671 at 25,000). Rank 64 and hidden 128 lost. Spread over learning rates 0.0062, over schedules 0.0002, over weight decay × batch 0.0028. Promotion gates failed (+0.000385; −0.001689). | R |
| Seed ensembling of the residual | K = 3: +0.0012 PDS on 4/5 lines — the only promotion gate passed in the project, but smaller than the per-line seed SD on 4/5 lines and negative on KOLF; not deployed. K = 5: +0.00004. | R |
| Direct prediction with no empirical target mean in input, offset or output | 0.574450 vs 0.786234 PDS on identical gene support; lost 15/15 paired seeds. | R |
| Learned change models on the fusion pipeline (GBM/RF stacking, learned gain, learned family weights, GBM/MLP residuals; Sept 25) | Lost or nil on H1 and KOLF. | R |
| 17-architecture screen (Sept 22; 104 runs, K562 + HepG2, seed 0, α 1) | Every joint, fully conditional model (CPA-style, biolord-style, joint MLP, conditional flow, flow matching) lost 0.04–0.22, almost entirely through PDS. Residual models on the refit differed by a few thousandths; the 141,888-parameter hypernetwork, the 86,496-parameter token Transformer and cross-attention ranked below 2k–35k-parameter models. | R |
| Seed replication of six of those (seeds 0–2, Nov 2026) | No architecture beat the refit baseline on the two-line mean with a CI excluding 0. Best: additive shift +0.0075 [−0.0010, +0.0161]. Additive helped K562 (+0.0240 [+0.0048, +0.0433]) and hurt HepG2 (−0.0090 [−0.0132, −0.0047]). MoE −0.0099, bilinear −0.0028, cross-attention −0.0166 on the two-line mean. | R |

Sources: `AGENTS.md`; `vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md`;
`output/data-analysis/results/rq2_arch/rq2_arch_seeds_summary.json`; `output/data-analysis/results/official_record.csv`.

Notes. The CPA- and biolord-style models were re-implementations under the screen's common protocol
(`vcc2026_data/scripts/arch_zoo.py`; cf. Lotfollahi et al. 2023, `lotfollahi-2023-cpa`; Piran et al. 2024,
`piran-2024-biolord`), not the released packages. No conditional
model recovered the per-target mean that the refit supplies for free. Learned components were not useless as such:
the S4 part (refit plus a FiLM residual at α 0.25) stayed in the system, and replacing it by plain fold changes cost
−0.0046 officially (S10e). What failed was adding capacity on top of the transfer baseline.

## (b) Context-conditioned and context-forcing models

| What was tried | Result | By |
|---|---|---|
| Eight context architectures on the Jiang CRISPRi screen (Jiang et al. 2025, `jiang-2025-pathway-signatures`; Sept 15; DeltaNet, hypernetwork, DeepSets, MoE, tensor, CNP, learned metric, GP; 184 outer folds, 3,300 fits) | Best gain over the denoiser +0.000267 PDS, below its own seed SD (0.000470); every Holm-adjusted p = 1.0. The GP trained on shuffled labels (0.779511) beat the correctly trained GP; the metric model's zero-query control beat its correct query. DeltaNet had a confirmed objective mismatch and was never re-run. | X |
| Context-interaction ridge on `[z, z × h]` | 0.77835 vs denoiser 0.77902; lost to its zero-context ablation in 5/6 lines. | X |
| LINCS context models (Subramanian et al. 2017, `subramanian-2017-connectivity-map`; 8 lines × 200 shRNA targets, CCLE expression as context) | Every context model lost to the context-free denoiser on 8/8 lines (0.937 / 0.935 vs 0.953). The ridge's correct-minus-no-context output was exactly a target-independent vector. | X |
| Learned signed correction on KOLF (Sept 14) | Normalised PDS 0.1570 → 0.2359, but a training-label-shuffled operator scored higher (0.2449): the mechanism was source-disagreement covariance denoising, not context. | R |
| Context forcing (Sept 16; weaker penalties, input dropout, interaction terms, contrast losses; 158 fits) | Zero query context (0.799576) beat the correct context (0.795414); correct beat a separately trained no-context model in 2/5 lines. In the linear family the context term cancels under query target-centring (replay change 1–2e-7). | R |
| Gated two-encoder context residual on the refit (Sept 19) | +0.0003 (0.8008 vs 0.8005) at the best α; performance fell with α; a gene-shuffled context scored above the correct one. | R |
| Pooled basal-state context (Sept 20; 9 biologies, linear / bilinear) | On the perturbation-specific 94 % of the deviation R² −0.002 to +0.001 against a noise ceiling of 0.28; correct context beat zero context in 5/9 folds. | X |
| PCA of pooled controls as context on the refit (Sept 20) | Positive in cosine on the unseen-context benchmark (+0.0268, 22.3 % of the oracle headroom). On the harness the wrong context (+0.0003 PDS) tied the best correct arm (+0.0006); HepG2 −0.0727 at α 1; nested α 0.7650, worse than any fixed α. Its upload (S2, −0.0020) used scrambled training profiles and is not evidence either way. | X, R |
| More real contexts for the FiLM residual (Sept 21–22; 12 CD4 donor × state + 19 HipSci lines, 38 in total) | Scaling 6 → 12 → 24 → 35 contexts was flat: PDS moved < 0.01 on every line and correct-minus-wrong did not grow (it fell on K562). Letting the new contexts vote as biologies cost −0.057 (refit) and −0.122 (mean transfer). The best arm (+0.0247, 15/15 seeds) gained through NMAE, reach and fidelity with PDS −0.014, and shipped as S4 at +0.0013 official. | R, O |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `AGENTS.md`; `vcc2026_data/reports/basal_context_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/ext_contexts_2026-09-21/FINDINGS.md`; `vcc2026_data/reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`;
`output/data-analysis/results/official_record.csv`.

Notes. Context information was not absent. In the architecture screen the correct context beat a real wrong context
for 10/11 residual models on K562 and 11/11 on HepG2 (+0.01 to +0.06), and a same-context oracle beat the
context-free prediction by +0.1953 cosine without leakage (composite headroom +0.27 to +0.48). What failed was turning
it into a net gain from 6–9 training biologies (`05-why-it-does-not-work.md`, A2–A3). Test-time signal from the destination's own
control-cell covariance discriminated at chance (0.49–0.52). The covariance denoiser itself shipped as B16
(`03-what-raised-the-score.md`).

Sources: `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md`; `vcc2026_data/reports/PROJECT_HISTORY.md`; `AGENTS.md`.

## (c) Foundation models and embeddings

| What was tried | Result | By |
|---|---|---|
| scLong (1.01 B parameters, frozen), zero-shot in-silico knockout (Sept 5–7) | No perturbation signal: mean \|Δ\| 0.0017 over 27,875 genes, the same 3–4 genes moved for every target. The model has no perturbation head. | X |
| scLong encoder + a trained GEARS head (34.5 M; cf. Roohani et al. 2024, `roohani-2024-gears`) | Epoch-0 validation MSE 0.06336 against 0.02929 for predicting the control; abandoned before any prediction. | X |
| scLong, STATE SE and PCA as context encoders that predict transferability (6 Replogle–Nadig line pairs) | Spearman with transfer quality: PCA −0.943, STATE SE −0.886, scLong +0.029 (Pearson wrong-signed). Six pairs only, underpowered. | X |
| Released STATE ST checkpoint (Adduri et al. 2026, `adduri-2026-state`) | Not runnable on the task: one-hot perturbation vocabulary of 2,024 essential genes, 0/300 official targets. | — |
| STATE ST trained from scratch on four lines with ESM-2 features (prior B) | Official −0.1910, below no change and 0.168 below the training-free transfer (prior A, −0.0228). In a K562 → RPE1 pilot it lost to no change on 100/100 perturbations (macro MSE 0.0350 vs 0.0162). | O, X |
| STATE SE-600M and UCE as frozen context encoders on the refit (Sept 20) | Not used. Five of seven biologies' control profiles were scrambled (sorted gene names over source-order counts), so the comparison with PCA is invalid. STATE SE identified cell lines with 100 % accuracy, which says nothing about responses. | — |
| ESM-2 15B gene embeddings (and ESM-C 600M at widths 32–256) in the residual's gene branch | Deleting the gene features was the best arm (+0.0021 PDS over the ESM-2 reference); random Gaussian vectors +0.0010; wider projections −0.0012 to −0.0024. The refit cannot use embeddings at all: every evaluation target is seen, so its prediction changes by exactly 0.000. | R |
| ESM-2 in the discriminative (InfoNCE) residual; ESM-2, GenePT, GO-BP and STRING as neighbour sources | With ESM-2 the residual acted as a lookup key and hurt H1; without it J > 0 on 7/7 panels. Seven alternative neighbour sources were each negative on H1 and KOLF; DepMap co-essentiality stayed best. | R, X |
| PIE (Verma et al. 2026, `verma-2026-pie`) as a DE-side re-ranker | Rejected on a cross-dataset unseen evaluation: discrimination 0.53, direction 0.65. | X |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `RESULTS.md`; `vcc2026_data/reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`;
`AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`; `vcc2026_data/reports/arch_screen_2026-09-22/FINDINGS.md`;
`vcc2026_data/reports/final_2026-10-27/README.md`; `output/data-analysis/results/official_record.csv`.

## (d) Drug data as a signal

| What was tried | Result | By |
|---|---|---|
| Tahoe-100M context models (Zhang et al. 2026, `zhang-2026-tahoe`; Sept 15; 16 lines × 150 drugs, 8 model families, 128 outer folds / 3,600 fits) | Largest increment 0.0000364 PDS; all eight predefined gates failed. Native pseudobulk PDS was ~0.999, so headroom was small. A full-response GP later recovered 1.445 % [0.540, 2.495] of repeatable contrast; not enough to pass. | X |
| PerturbedVAE variants on Tahoe | A deterministic anchored residual (+2.56 % repeat recovery, 14/16 lines) beat every VAE; unshrunk it was −1.874 %, so the shrinkage was the method. Frozen drug-pretrained transfer lost to frozen-random on Jiang (−0.2627 pp). | X |
| Tahoe rows in a small genetic residual (Sept 16; four RN lines) | Refit α 0: 0.902890; joint genetic + drug 0.896723; drug pretraining then genetic fine-tuning 0.900244. No residual beat α 0. | R |
| Same-target drug → gene concordance (Sept 16; A549, HT-29; Jiang CRISPRi, LINCS shRNA) | Same-target Spearman 0.0059–0.0383 (Jiang, four stimuli) and 0.0177 (LINCS, vs 0.0069 for wrong targets). Jiang's own replicate correlation for the same knockdown was 0.069–0.115. Reliability-corrected bound (post-challenge): ρ* = 0.058 [−0.019, 0.103], i.e. at most ~0.3 % of reproducible CRISPRi variance; most permissive stimulus TNFA 0.153 [0.036, 0.253]. | X |
| Drug-context pretraining, drug side (Sept 20; Tahoe 48 lines + GSE306429; FiLM r32) | Context discrimination grew ~4× from 6 to 40 contexts (correct − wrong 0.025 → 0.098, 6/6 folds), but the total-response gain stayed negative (−0.025, −0.028, −0.008 at 6, 12, 24) and reached break-even only at 40 (+0.002). GSE contexts were net useful on 0/10 after a key-collision fix removed an artefactual signal. | X |
| Drug-context pretraining, gene side, and its upload S3 | **Invalid; not used.** The gene-side context vectors were scrambled for 5 of 7 biologies, so the drug → gene verdicts of that experiment are not established. Its +0.0180 local arm shipped as S3 at −0.0005 official; that contrast is invalid too. | — |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `AGENTS.md`; `vcc2026_data/reports/drug_gene_concordance_2026-09-16/RESULTS.md`;
`output/data-analysis/results/appx_drug_gene.json`; `vcc2026_data/reports/drug_context_pretrain_2026-09-20/FINDINGS.md`;
`vcc2026_data/reports/ext_contexts_2026-09-21/CONTROLS_COMMON_BUG.md`; `output/data-analysis/results/e4_official.json`.

Notes. The ρ* bound may be anti-conservative if replicates share control cells (`knowledge/synthesis/appx-drug-gene-bound.md`).
The pre-registered post-challenge learning curve on Tahoe, which replaces the in-challenge drug-side curve, is in
`05-why-it-does-not-work.md`, A3.

## (e) Other public sources and learned source weighting

| What was tried | Result | By |
|---|---|---|
| LINCS L1000 shRNA (GSE106127, 9 lines) as extra source contexts | Cosine with the truth on landmark genes H1 0.007–0.026 vs ours 0.035–0.113; on strong targets 0.005–0.065 vs 0.11–0.24. The LINCS HepG2 line was not even the best LINCS line on the HepG2 panel. | X |
| LINCS2020 CRISPR-KO consensus (Oct 29; 5,157 knocked-out genes; 75/300 final-panel targets) | Alone, 1/2 to 1/6 of our cosine (H1 0.018 vs 0.105). As an extra family mean J ≈ +0.002 (≈ +0.0003 composite), KOLF −0.0065, higher weights worse. Rejected. | X |
| VIPerturb-seq K562 (Flex chemistry) | Worse than Replogle K562 on the same targets (~47 cells per target). As an extra family (Oct 28): OR_X J +0.019, H1 −0.0025; weight ≥ 4 negative everywhere; replacing K562 cost 0.018–0.024 PDS. Rejected. | X |
| GSE249595 Jurkat (high MOI, ~12 guides per cell) as a same-line source for context A | Split-half per-target cosine +0.003 on official targets; a per-cell regression re-estimate agreed with the marginal one (cosine 0.98). Unusable. | X |
| Small public screens (melanoma, lung, HEK, iPSC neurons, THP-1, primary T; Sept 29) | Each within ±0.005 J of zero at weight 1; they overlap a 300-target panel at chance (~3.5 genes). Two new stem-cell screens covered 0–7 of 300 targets. | X |
| DepMap beyond the co-essentiality neighbour blend | Destination-line essentiality ranked the true per-target energy at Spearman 0.16–0.36, below our own delta norm (0.46–0.61). A larger blend (0.5) bundled with two knobs (S10f) scored −0.0015 officially against a harness +0.006. A fitness-from-expression ridge passed its gates (R² 0.045) but predicts fitness, not RNA responses. | X, O |
| Response-space kNN transfer of K562 GWPS (Oct 28; 9,729 perturbations) | MSE down on 6/6 panels, PDS down on 5/6; best J +0.0075; gated to the strongest 30 % J +0.011 (≈ +0.002 composite). | X |

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`; `vcc2026_data/reports/final_2026-10-27/README.md`; `AGENTS.md`;
`output/data-analysis/results/official_record.csv`; `vcc2026_data/reports/PROJECT_HISTORY.md`.

**Learned source weighting — six failures.** The per-target non-negative least-squares oracle over the sources was
large (H1 +0.048 PDS, Orion +0.130; full-recipe J on H1 +0.125), so the headroom was real. Six learners tried to
predict it from test-time features:

| # | Learner | Result (X) |
|---|---|---|
| 1 | MODZ mutual-agreement weights (powers 0.5 / 1 / 2) | Worse than fixed weights in the same construction: H1 PDS 0.8589 vs 0.8640, Orion 0.7007 vs 0.7117. |
| 2 | Leave-one-panel-out GBM on test-time features | Recovered −17 % of the oracle on average, indistinguishable from its shuffled control. |
| 3 | Similarity of destination and source control profiles | Pooled within-panel Spearman with family value +0.189 (−0.140 to +0.392 per panel). |
| 4 | GBM trained on the oracle weight vectors, full recipe rebuilt faithfully | H1 J +0.004 / +0.005 / −0.006 at λ 0.25 / 0.5 / 1 (PDS −0.003 to −0.014); negative on all four RN panels (down to −0.098). |
| 5 | Low-dimensional learned gate on family weights, trained on pseudo-destinations | PDS fell on 6/6 panels; J from −0.046 (RN_Jurkat) to +0.029 (OR_HCT116); stopped by its pre-registered kill rule. |
| 6 | Mixer trained on the official-pipeline H1 truth (MSE + own cosine; InfoNCE variant) | Out-of-fold PDS −0.0098 / −0.0085, MSE −0.0346 / −0.0299, J +0.0145 / +0.0125 against a kill bar of 0.03. |

Every learner moved along the PDS–MSE trade-off (MSE down, PDS down) rather than beyond it. A first full-recipe
measurement of learner 4 was withdrawn because of an injection defect; the faithful rerun is quoted.

Sources: `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§18.4, §19, §20.16, §22.5, §23.6);
`vcc2026_data/reports/sc_hybrid_2026-09-25/debate_loop/ROUNDS.md`; `AGENTS.md`.

## (f) Amplitude, budget and metric tricks

| What was tried | Result | By |
|---|---|---|
| Global delta scale beyond ×3 | ×4 (B10) 0.0868 vs ×3 (B8) 0.1217. A later five-line sweep put the optimum at 2.5, +0.004 over 3.0 (inside seed noise), with per-line optima 2.0–3.0. | O, R |
| Expression-gated per-context multipliers (B11, eight fitted numbers) | Harness +0.0085 (4/4 seeds); official 0.1156 vs 0.1217. | R, O |
| Shrinking emitted cells toward the control mean (λ 0.25–0.5) | −0.03 to −0.11 on K562 at every scale. | R |
| Per-target amplitude on the final-phase harness (Sept 24) | Per-target amplitude −0.091 (0/6 panels), amplitude shrink −0.047, expression multiplier −0.010, low-rank lr64 −0.150. | R |
| Per-target amplitude, shape half (`tzb0.25`, budget held) | H1 +0.0022, OR_X −0.0108 (NMAE +0.0373), RN_Jurkat −0.0013; gate composite −0.0042. Letting the DE side follow the shrink made NMAE worse (+0.0530). | R |
| Per-target amplitude, budget half (`tz0.25x0.9`) | Offline ΔMSE −0.0211, real ΔMSE +0.0055; composite +0.0001. | X, R |
| Better learned amplitude predictors (v2, v3) in M2's shipped rank form | +0.0003 / −0.0000; the v3 features raised mean Spearman only from 0.6722 to 0.6736. | X |
| MSE kink factor 1.17 → 1.35 | OR_X composite 0.3942 → 0.2130 (−0.181; MSE 0.0299 → 1.2108) where the offline probe had claimed +0.0768; RN_Jurkat +0.0075 for the same +7 % budget. | R |
| Official sum-budget probes from S16 | Up (S17u, margin 1.2): −0.0462, MSE member −0.2737. Down (S17d, margin 0.6, floor 0): −0.0073. | O |
| More emitted per-cell dispersion | Impossible to use: the scorer's cap already clipped our claim on 263/263 OR_X targets. | R |
| Strength-gated hubness (`disw`); continuous DE rank calibration (`drp`); DE-side `lcp0.25` | +0.0007; −0.0013; gate composite +0.0026 (kept dormant). | X, R |
| Ensembles without member diversity | Neighbour-cut ensemble `nbe` +0.0002 / +0.0005 (neighbour lists ≤ 20 long, so cuts 20/40/80 were identical); seed ensemble K = 5 +0.00004; rank-pool members ~0; averages over PC-removal λ or centring are provable no-ops; seven alternative neighbour sources were biased the same way, so their mean lost; drug-pretrained ensembles +0.0010 (plain) and 0.0000 (gated), but these used the scrambled gene-context vectors and are invalid (group (d)). | R, X |
| Small harness-tuned knobs | S10f (blend 0.5, concentration share 0.2, K 5000): harness H1 +0.006 → official −0.0015. S7a (fusion share 0.75): −0.0063. The two DE-side fixes in S11 moved official NMAE by −0.0000165 and reach by +0.0021. | O |
| Predicting the DE call set and its size | Count per target weakly predictable (Spearman H1 0.39, KOLF 0.16); a leave-one-panel-out count under-called badly (H1 fidelity 0.23–0.34 vs an oracle 0.590); over-calling stayed. | R |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `vcc2026_data/reports/recon_sweep_2026-09-22/FINDINGS.md`; `AGENTS.md`;
`vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md` (§4, §24.6, §25.1–25.12); `vcc2026_data/reports/drug_context_pretrain_2026-09-20/FINDINGS.md`;
`output/data-analysis/results/official_record.csv`.

## (g) The late "big-lever" search (Oct 28–29)

After the validation phase closed, every untested public lever was screened with the exact offline evaluator on the
S16 base (CPU only). OR_X offline MSE was known to be unreliable, so PDS was read there.

| Lever | Result (X) |
|---|---|
| Response-space kNN of K562 GWPS | Best J +0.0075; gated to strong targets +0.011 (≈ +0.002 composite). |
| Per-gene transfer weights learned from truth (LOPO, or H1 only) | H1-trained: OR_X PDS +0.014, KOLF +0.0067; H1 itself under LOPO +0.0007 PDS (nil). The weights mostly shrank genes (median 0.65). |
| VIPerturb-seq as a family | Rejected (group (e)). |
| Sum-side expression down-weighting (`ew`) re-tested with budget match | OR_X PDS +0.0139, KOLF +0.0066, H1 +0.0002 PDS / −0.0005 MSE (nil). Kept in F2 because it added to the deflation on every panel (OR_X +0.012 PDS over F1, H1 +0.0007); its official effect is not known. |
| Stacks on the S18 deflation, M3 (ΔPDS OR_X / H1) | Deflation alone +0.011 / +0.0008; + `ew` +0.0227 / +0.0015; + `ew` + gene weights +0.0267 / +0.0006; a stronger deflation + `ew` +0.040 / −0.0027 (alone on H1 −0.0043). |
| Strong-target direction levers (top 30 % by predicted amplitude, Oct 29) | Boosting along the destination's control PCs: H1 J +0.0112, but PDS fell on every panel. Down-weighting K562 on strong rows: J +0.007 to +0.016 on four panels, −0.016 on Jurkat. Best single lever ≈ +0.002 composite on H1. |
| LINCS CRISPR-KO consensus; PIE | Rejected (groups (e), (c)). |

The recorded verdict: no public mechanism gave a verified +0.03 to +0.05; the sum side was saturated at ≤ +0.005 per
lever; the stack was worth about +0.003 to +0.006 composite only if official PDS followed OR_X (H1 said ≈ 0).

Sources: `vcc2026_data/reports/final_2026-10-27/README.md`; `AGENTS.md`.

## (h) Source combination, preprocessing and estimator variants

| What was tried | Result | By |
|---|---|---|
| Sign-agreement weighting of sources (B7) | −0.0994 against its parent B6. | O |
| Precision / effective-cell-count source weighting (B14 at ×3.5; B15 at ×3) | B14 0.0888; B15 0.1067 vs B8 0.1217 (PDS −0.0455), after the KOLF harness had predicted a gain. | O |
| Per-source amplitude calibration; true log fold change instead of the log1p compression | −0.003 PDS; −0.025 PDS. | R |
| Robust pseudobulk aggregation and inverse-SE weighting within biologies | Every arm lost on all five lines locally (weighting −0.0006, robust −0.002, median aggregation −0.029); on the local official-preset composite (seed 0) −0.0007 (weighting) and −0.0004 (nested robust). | R |
| Gene subsets ranked by cross-biology reliability (500–5,000 genes) | Never beat all 5,592 genes (−0.029 to −0.001). | R |
| Wide gene axis (17,282 genes, masked refit) | +0.00105 PDS, but only 2/5 lines and 4/15 seeds; KOLF lost. Failed the gate. | R |
| Per-source z-shrinkage before fusion; rank-space fusion | Fusion-only shrinkage: PDS down on 5/6 panels. Rank space: H1 J −0.045, Orion −0.213. | X |
| Removing K562-specific (erythroid / myeloid) programs; per-gene source-expression weighting `xw` | RN −0.015 to −0.064; mean J −0.0001 / +0.0019 (4/7 positive), again an MSE-for-PDS move. | X |
| Cell-cycle composition correction; knockdown-escape (responder) weighting | ≈ 0 on H1 (composition is 0.1–0.5 % of H1/KOLF response energy); responder-weighted Orion sources null (0.2156 vs 0.2157). | R |

Sources: `vcc2026_data/reports/PROJECT_HISTORY.md`; `AGENTS.md`; `vcc2026_data/reports/sc_hybrid_2026-09-25/FINDINGS.md`;
`output/data-analysis/results/official_record.csv`.

## Pattern across the groups

Four patterns recur. (1) Local gains that were not carried by a structural change were either one-line effects (RPE1,
HepG2, K562) or shrank on the leaderboard. (2) Learned weighting, residuals and amplitude reshaping moved along a
PDS-for-MSE trade-off instead of improving both. (3) The destination's context carried real information, but no model
fitted on 6–9 training biologies converted it into a net gain. (4) External data either did not cover the 300-target
panel or did not transfer across cell types or perturbation modalities. `05-why-it-does-not-work.md` gives the
reasons.
