---
title: "Literaturguide — zero-shot limits & context sample complexity"
type: note
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
plan: zero-shot-limits-context-complexity-plan
---

> ## ⚠ STOP-FLAG (qualified — raise with the user before Task 8)
>
> **Ji et al. 2024, Prophet** (`ji-2024-universal-prediction`, bioRxiv 10.1101/2024.08.12.607533, unpublished
> preprint, 11 authors, Theis lab). It reports an **empirical learning curve over the share of cell lines (cell states)
> seen in training, against a context-free mean baseline**, and says the curve crosses that baseline at a
> dataset-dependent point:
> - Results text (v2 full text, read 2026-11-05): *"we trained Prophet-individual models to predict the effect of either
>   unseen treatments or unseen cell states, including 5-80% of the data along the respective axis. The amount of data
>   required to outperform the mean estimation varies per dataset … (fig. S4, S5)."*
> - Fig. S4A legend: *"Performance of Prophet along the cell line axis."* Fig. S5A defines the baseline: *"The mean
>   baseline estimates that the effect of an intervention in an unseen cell state is the average value of that
>   intervention in all previously seen states."* That is exactly our context-free model (μ + a_t).
> - *"Prophet encounters challenges in outperforming the mean estimator in certain scenarios, such as predicting the
>   effects of observed treatments in unobserved cell lines within SCORE"*. They attribute this to low replicate
>   correlation (Table S4), not to the number of contexts.
>
> **Why this is a qualified flag, not a pre-emption:**
> 1. It is a share-of-lines curve (5–80 %), without a reported break-even count C*.
> 2. There is no theory: nothing on identifiability rank ≤ C−1, descriptor dimension d, or a noise term.
> 3. The data are mostly viability/IC50 (GDSC, CTRP, PRISM), CRISPR viability (SCORE) and L1000. There is no CRISPRi
>    Perturb-seq and no unseen-line transcriptome break-even.
> 4. The figure image could not be retrieved (bioRxiv blocks automated access), so the per-dataset crossing points are
>    **unread**. **The user must look at fig. S4A manually** and record whether a transcriptomic dataset (LINCS) shows a
>    crossing and at what share of lines.
>
> **Implication:** H2a's empirical curve stays new for CRISPRi/Tahoe transcriptomes, but the paper must cite Prophet
> fig. S4 as the first "crosses the mean baseline only after enough lines" observation. The novelty claim moves to the
> explicit C* formula, its independent-parameter prediction (H2b), and the genetic-screen bound (H2c).
>
> **Near-misses (no flag):**
> - **Dibaeinia et al. 2026, Fig. 4A–F.** Observational stratification: DEG-F1 rises from < 0.1 to ≈ 0.19 between
>   ≤ 3 and 8 training contexts per perturbation. There is no context-free baseline on that curve and no crossing. The
>   text says: *"We do not address how much contextual diversity is “enough,” as this likely depends on the target application."*
> - **Molina & Zhang 2026, Fig. 3a.** 1 vs 3 source lines for a model without a context descriptor; no break-even.
> - **Shen 2026.** Tahoe leave-one-line-out at a single C = 46: no model beats the drug mean (best R² +0.004). This is a
>   single point, not a curve, and it is a *competing explanation* (reliability), which C* absorbs through
>   R²·Var(i).
> - **Wei et al. Nat Methods**: more contexts do not help, similarity matters (Suppl. Note 11).
> - **PerturBench Fig. 3**: nested states, no baseline crossing.
> - **Partin 2021**: learning curves over cell–drug pairs, not over lines.
>
> **Theory precedent (not a flag, but it bounds the novelty claim):** "you need more than d environments" is
> established in domain-generalization theory: Rosenfeld et al. 2021 (E > d_s); Chen et al. 2022 (lower bound at
> E ≤ d_s); Wang et al. 2022 (ISR-Mean needs d_s+1, the same C−1 span argument as Proposition 1). The random-X OLS risk
> σ²d/(n−d−1) behind Proposition 2 is classical (Breiman & Freedman 1983; Rosset & Tibshirani 2020; Hastie et al. 2022).
> Present Propositions 1 and 2 as **transfers of known results to cell-line contexts plus a break-even against a
> context-free baseline**, not as new theory.

# Literaturguide: zero-shot CRISPRi prediction in unseen cell lines — ceilings and context sample complexity

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (pre-registered), Task 2. Search date 2026-11-05.
Bibkeys follow the project schema `surname-year-shorttitle`. The year is the print-issue year where Crossref has one,
otherwise the first-publication year. All 102 entries in `output/bibtex/references.bib` were generated
programmatically from Crossref `works/{doi}` or arXiv API metadata. No author, DOI or venue was typed by hand. Venue
fields for conference papers that only arXiv could confirm carry a `note` where the proceedings record was not
re-verified.

## 1. Quellentabelle (gewichtet)

**Counts:** A = 33 · B = 37 · C = 32 (total 102).

The last column is an addition to the fixed five-column schema: bibkey plus one line of relevance to
H1 / H2a–c / H3 / H4. Prop. 1 = identifiability (rank ≤ C−1); Prop. 2 = break-even C*. "Appendix (drug ↔ genetic)"
refers to the drug–gene concordance appendix.

| Grade | Autor Jahr | Kurztitel | OA/Zugang | DOI/Link | Bibkey — relevance |
|-------|------------|-----------|-----------|----------|--------------------|
| A | Breiman & Freedman 1983 | how-many-variables | ◐ UB | [10.1080/01621459.1983.10477941](https://doi.org/10.1080/01621459.1983.10477941) | `breiman-1983-how-many-variables` — H2b: classical random-X break-even for adding regressors; antecedent of C*. |
| A | Maurer et al. 2016 | multitask-representation | ✓ OA | [arXiv:1505.06279](https://arxiv.org/abs/1505.06279) | `maurer-2016-multitask-representation` — H2: O(1/√T) bounds in number of tasks. |
| A | Rosset & Tibshirani 2020 | random-x | ◐ UB | [10.1080/01621459.2018.1424632](https://doi.org/10.1080/01621459.2018.1424632) | `rosset-2020-random-x` — H2b: exact random-X OLS excess risk σ²p/(n−p−1) — the C−d−1 term. |
| A | Tripuraneni et al. 2020 | task-diversity | ✓ OA | [arXiv:2006.11650](https://arxiv.org/abs/2006.11650) | `tripuraneni-2020-task-diversity` — Prop. 1: task diversity must span the representation. |
| A | Rosenfeld et al. 2021 | risks-irm | ✓ OA | [arXiv:2010.05761](https://arxiv.org/abs/2010.05761) | `rosenfeld-2021-risks-irm` — Prop. 1 / H2: linear-Gaussian IRM needs E > d_s environments. |
| A | Blanchard et al. 2021 | marginal-transfer | ✓ OA | [arXiv:1711.07910](https://arxiv.org/abs/1711.07910) | `blanchard-2021-marginal-transfer` — Prop. 2 framing: predictor takes a marginal (control) descriptor; risk bounds in number of training domains. |
| A | Tripuraneni et al. 2021 | linear-representations | ✓ OA | [arXiv:2002.11684](https://arxiv.org/abs/2002.11684) | `tripuraneni-2021-linear-representations` — Prop. 1/2: needs t ≥ r tasks; explicit negative transfer regime (break-even analog). |
| A | Du et al. 2021 | few-shot-representation | ✓ OA | [arXiv:2002.09434](https://arxiv.org/abs/2002.09434) | `du-2021-few-shot-representation` — Prop. 1: source tasks must have rank k (T ≥ k). |
| A | Rothenhäusler et al. 2021 | anchor-regression | ✓ OA | [10.1111/rssb.12398](https://doi.org/10.1111/rssb.12398) | `rothenhausler-2021-anchor-regression` — Prop. 1: robustness only along directions spanned by training anchor (environment) variation. |
| A | Replogle et al. 2022 | genome-scale | ◐ UB (preprint ✓ OA) | [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013) | `replogle-2022-genome-scale` — Data: K562 genome-wide + K562/RPE1 essential CRISPRi (harness panels RN K562, RN RPE1). |
| A | Hastie et al. 2022 | ridgeless-interpolation | ✓ OA | [10.1214/21-AOS2133](https://doi.org/10.1214/21-AOS2133) | `hastie-2022-ridgeless-interpolation` — H2b + Prop. 1: underparameterized risk σ²γ/(1−γ); min-norm solution confined to row span (C < d). |
| A | Chen et al. 2022 | iterative-feature-matching | ✓ OA | [arXiv:2106.09913](https://arxiv.org/abs/2106.09913) | `chen-2022-iterative-feature-matching` — Prop. 1: information-theoretic lower bound, no generalization with E ≤ d_s; O(log d_s) under extra structure (reviewer question). |
| A | Wang et al. 2022 | subspace-recovery | ✓ OA | [arXiv:2201.12919](https://arxiv.org/abs/2201.12919) | `wang-2022-subspace-recovery` — Prop. 1: ISR-Mean recovers invariant subspace from E class means, needs E ≥ d_s+1 — same C−1 span argument. |
| A | Ji et al. 2024 | universal-prediction | ✓ OA | [10.1101/2024.08.12.607533](https://doi.org/10.1101/2024.08.12.607533) | `ji-2024-universal-prediction` — STOP-FLAG (qualified), H2a: Prophet; fig. S4 performance vs share (5-80 %) of cell lines seen vs mean-over-seen-lines baseline; "amount of data required to outperform the mean estimation varies per dataset"; never beats mean on SCORE unseen lines. |
| A | Ahlmann-Eltze et al. 2025 | linear-baselines | ✓ OA | [10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6) | `ahlmanneltze-2025-linear-baselines` — H4: deep models do not beat linear/additive baselines. |
| A | Nadig et al. 2025 | transcriptome-wide | ◐ UB (preprint ✓ OA) | [10.1038/s41588-025-02169-3](https://doi.org/10.1038/s41588-025-02169-3) | `nadig-2025-transcriptome-wide` — Data + H1: HepG2/Jurkat essential CRISPRi; TRADE noise-aware effect sizes. |
| A | Huang et al. 2025 | x-atlas-orion | ✓ OA | [10.1101/2025.06.11.659105](https://doi.org/10.1101/2025.06.11.659105) | `huang-2025-x-atlas-orion` — Data: X-Atlas/Orion genome-wide Perturb-seq (OR_HCT116X panel); sgRNA dose/depth. |
| A | Molina & Zhang 2026 | response-decomposition | ✓ OA | [10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459) | `molina-2026-response-decomposition` — H1 H3 H2: source of the decomposition and split-half noise correction (template 27.8 / conserved 29.4 / interaction 23.5 / noise 19.3 %); projects the response on its own line's control PCs (random-basis null only); 1 vs 3 source lines. Overlap gate: RQ3 keep. |
| A | Huu Loc & Van 2026 | cell-quota | ✓ OA | [10.64898/2026.10.02.754788](https://doi.org/10.64898/2026.10.02.754788) | `huuloc-2026-cell-quota` — H1: closed-form cell quota per perturbation direction (cells, not contexts); Tahoe 55 % control-pool-limited. Complements the C-axis (contexts). |
| A | Li et al. 2026 | systematic-comparison | ✓ OA | [10.1126/sciadv.aed3414](https://doi.org/10.1126/sciadv.aed3414) | `li-2026-systematic-comparison` — H1 H2: 13 methods × 25 datasets; no method generalizes across cell types; split-half control noise floor. |
| A | Dibaeinia et al. 2026 | context-scale | ✓ OA | [10.64898/2026.02.04.703804](https://doi.org/10.64898/2026.02.04.703804) | `dibaeinia-2026-context-scale` — H2: position paper, context coverage is the bottleneck; Fig. 4 = DEG-F1 vs contexts per perturbation (≤3 → 8), no baseline crossing; "we do not address how much diversity is enough". |
| A | Piao et al. 2026 | adaptive-contexts | ✓ OA | [arXiv:2602.18885](https://arxiv.org/abs/2602.18885) | `piao-2026-adaptive-contexts` — H4 (H1-adjacent): AdaPert, ICML 2026; mean collapse and noise; unseen-perturbation focus. |
| A | Adduri et al. 2026 | state | ◐ UB (preprint ✓ OA) | [10.1016/j.cell.2026.07.052](https://doi.org/10.1016/j.cell.2026.07.052) | `adduri-2026-state` — H2 H4: STATE; Replogle-Nadig "leave-one-line-out" is few-shot (30 % of target perturbations); Fig. S9D scales perturbations, not lines. |
| A | Zhang et al. 2026 | tahoe | ◐ UB (preprint ✓ OA) | [10.1016/j.cell.2026.08.035](https://doi.org/10.1016/j.cell.2026.08.035) | `zhang-2026-tahoe` — H2a/H2b data: 50 lines × 1,100 drug-dose conditions; no sweep over number of lines. |
| A | Viñas Torné et al. 2026 | systema | ✓ OA | [10.1038/s41587-025-02777-8](https://doi.org/10.1038/s41587-025-02777-8) | `vinastorne-2026-systema` — H3b H4: systematic variation; Systema-style common-response removal is an H3b alternative. |
| A | Adduri et al. 2026 | virtual-cell-challenge | ◐ UB | [10.1016/j.cell.2026.08.004](https://doi.org/10.1016/j.cell.2026.08.004) | `adduri-2026-virtual-cell-challenge` — H2 framing: VCC 2026 = zero-shot CRISPRi prediction in unseen cell lines (exactly our setting). |
| A | Pan et al. 2026 | global-cell-state | ✓ OA | [10.64898/2026.05.16.725005](https://doi.org/10.64898/2026.05.16.725005) | `pan-2026-global-cell-state` — Data H4: GSE295214 genome-scale CRISPRi in hESC (H1GW source); conserved vs context-specific responses. |
| A | Shen 2026 | basal-expression | ✓ OA | [10.21203/rs.3.rs-10846736/v1](https://doi.org/10.21203/rs.3.rs-10846736/v1) | `shen-2026-basal-expression` — H1 H2 H3: Tahoe 47-line LOCO (C = 46); context residual vs drug-mean: best R² +0.004, basal-PC ridge −0.068; skill tracks repeat reliability (ρ = 0.76). Competing explanation to C < C*. |
| A | Wang et al. 2026 | reliable-perturbations | ✓ OA | [10.64898/2026.08.11.744177](https://doi.org/10.64898/2026.08.11.744177) | `wang-2026-reliable-perturbations` — H1: split-half reliability, ceiling |r| ≤ √ρ, "ceiling efficiency"; 65 % of 7,170 perturbations unreliable. |
| A | Huang 2026 | shared-clones | ✓ OA | [arXiv:2610.32962](https://arxiv.org/abs/2610.32962) | `huang-2026-shared-clones` — H1 caveat: random cell halves inflate split-half ceilings via shared clones; small in CRISPRi cell-line screens. |
| A | Verma et al. 2026 | pie | ✓ OA | [10.64898/2026.10.02.756297](https://doi.org/10.64898/2026.10.02.756297) | `verma-2026-pie` — H2 H3: PIE, current SOTA for unseen context + unseen perturbation on Replogle-Nadig; natural comparator. |
| A | Wei et al. 2026 | benchmarking-algorithms | ◐ UB | [10.1038/s41592-025-02980-0](https://doi.org/10.1038/s41592-025-02980-0) | `wei-2026-benchmarking-algorithms` — H2: 27 methods × 29 datasets; Suppl. Note 11 more contexts do not help; context similarity matters (Fig. 3e). |
| A | Wenkel et al. 2026 | txpert | ✓ OA | [10.1038/s41587-026-03113-4](https://doi.org/10.1038/s41587-026-03113-4) | `wenkel-2026-txpert` — H1 H2: split-half / resampled reproducibility ceiling; leave-one-line-out on 4 lines. |
| B | Copas 1983 | regression-shrinkage | ◐ UB | [10.1111/j.2517-6161.1983.tb01258.x](https://doi.org/10.1111/j.2517-6161.1983.tb01258.x) | `copas-1983-regression-shrinkage` — H2b: random-X OLS prediction risk and shrinkage. |
| B | Baxter 2000 | inductive-bias | ✓ OA | [10.1613/jair.731](https://doi.org/10.1613/jair.731) | `baxter-2000-inductive-bias` — H2: first number-of-tasks bound. |
| B | Maas & Hox 2005 | sample-sizes | ◐ UB | [10.1027/1614-2241.1.3.86](https://doi.org/10.1027/1614-2241.1.3.86) | `maas-2005-sample-sizes` — H2c: few groups bias group-level variance components (random-effects view of C). |
| B | Dicker 2013 | equivariant-prediction | ✓ OA | [10.1214/13-EJS826](https://doi.org/10.1214/13-EJS826) | `dicker-2013-equivariant-prediction` — H2b: exact Gaussian-design OLS predictive risk (inverse-Wishart moment). |
| B | McNeish & Stapleton 2016 | small-sample | ◐ UB | [10.1007/s10648-014-9287-x](https://doi.org/10.1007/s10648-014-9287-x) | `mcneish-2016-small-sample` — H2c: small number of clusters in two-level models. |
| B | Arjovsky et al. 2019 | irm | ✓ OA | [arXiv:1907.02893](https://arxiv.org/abs/1907.02893) | `arjovsky-2019-irm` — Prop. 1: environments in linear general position. |
| B | Kong et al. 2020 | mixed-linear-regression | ✓ OA | [arXiv:2002.08936](https://arxiv.org/abs/2002.08936) | `kong-2020-mixed-linear-regression` — H2: tasks needed for meta linear regression. |
| B | Partin et al. 2021 | learning-curves | ✓ OA | [10.1186/s12859-021-04163-y](https://doi.org/10.1186/s12859-021-04163-y) | `partin-2021-learning-curves` — H2: drug-response learning curves over cell–drug pairs (not over number of lines). |
| B | Kamath et al. 2021 | irm-invariance | ✓ OA | [arXiv:2101.01134](https://arxiv.org/abs/2101.01134) | `kamath-2021-irm-invariance` — Prop. 1: limits of IRM with few environments. |
| B | Lucas et al. 2021 | meta-learning-bounds | ✓ OA | [arXiv:2010.07140](https://arxiv.org/abs/2010.07140) | `lucas-2021-meta-learning-bounds` — H2: minimax lower bounds in number of tasks. |
| B | Christiansen et al. 2022 | causal-framework | ✓ OA | [10.1109/TPAMI.2021.3094760](https://doi.org/10.1109/TPAMI.2021.3094760) | `christiansen-2022-causal-framework` — Prop. 1: generalization only within span of training interventions. |
| B | Lotfollahi et al. 2023 | cpa | ✓ OA | [10.15252/msb.202211517](https://doi.org/10.15252/msb.202211517) | `lotfollahi-2023-cpa` — Baseline model with cell-line covariate embedding (CPA). |
| B | Müller et al. 2023 | context-aware-dg | ✓ OA | [arXiv:2312.10107](https://arxiv.org/abs/2312.10107) | `muller-2023-context-aware-dg` — H2b: when a set/context descriptor helps vs hurts (break-even analog). |
| B | Roohani et al. 2024 | gears | ✓ OA | [10.1038/s41587-023-01905-6](https://doi.org/10.1038/s41587-023-01905-6) | `roohani-2024-gears` — Baseline model (GEARS). |
| B | Cui et al. 2024 | scgpt | ◐ UB (preprint ✓ OA) | [10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0) | `cui-2024-scgpt` — Baseline foundation model (scGPT). |
| B | Wenteler et al. 2024 | perteval | ✓ OA | [10.1101/2024.10.02.616248](https://doi.org/10.1101/2024.10.02.616248) | `wenteler-2024-perteval` — H4: PertEval-scFM, limited gains under distribution shift. |
| B | Roohani et al. 2025 | virtual-cell-challenge | ◐ UB | [10.1016/j.cell.2025.06.008](https://doi.org/10.1016/j.cell.2025.06.008) | `roohani-2025-virtual-cell-challenge` — Framing: VCC 2025 (H1 hESC CRISPRi) design and metrics. |
| B | Mejia et al. 2025 | diversity-design | ✓ OA | [arXiv:2506.22641](https://arxiv.org/abs/2506.22641) | `mejia-2025-diversity-design` — H1: technical-duplicate (split-half) baseline as ceiling; calibrated metrics. |
| B | Yu et al. 2025 | perturbnet | ✓ OA | [10.1038/s44320-025-00131-3](https://doi.org/10.1038/s44320-025-00131-3) | `yu-2025-perturbnet` — H4 context: PerturbNet unseen perturbations. |
| B | Csendes et al. 2025 | foundation-models | ✓ OA | [10.1186/s12864-025-11600-2](https://doi.org/10.1186/s12864-025-11600-2) | `csendes-2025-foundation-models` — H4: train-mean baseline beats scGPT/scFoundation. |
| B | Wu et al. 2025 | perturbench | ✓ OA | [10.52202/085713-3225](https://doi.org/10.52202/085713-3225) | `wu-2025-perturbench` — H2: covariate transfer; nested-state data scaling (Fig. 3), no break-even. |
| B | Dwork et al. 2025 | how-many-domains | ✓ OA | [arXiv:2506.16704](https://arxiv.org/abs/2506.16704) | `dwork-2025-how-many-domains` — H2: domain sample complexity via domain shattering dimension. |
| B | Chen et al. 2026 | signature-recontextualization | ✓ OA | [10.64898/2026.08.14.744937](https://doi.org/10.64898/2026.08.14.744937) | `chen-2026-signature-recontextualization` — H3 H4: control-only regime never beats no-change; target coverage 10 % / 90 %. |
| B | Huang et al. 2026 | recoverable-resolution | ✓ OA | [10.64898/2026.09.30.755756](https://doi.org/10.64898/2026.09.30.755756) | `huang-2026-recoverable-resolution` — H1 H3: models collapse context-specific responses ("context geometry compression"). |
| B | Miller et al. 2026 | calibrated-metrics | ✓ OA | [10.1038/s41587-026-03307-w](https://doi.org/10.1038/s41587-026-03307-w) | `miller-2026-calibrated-metrics` — H1 H4: positive-control baselines and metric calibration. |
| B | Nicol et al. 2026 | spurious-correlation | ✓ OA | [10.64898/2026.05.07.723486](https://doi.org/10.64898/2026.05.07.723486) | `nicol-2026-spurious-correlation` — H1/H3 method: shared-control reuse inflates Pearson-Δ; use split controls (cf. RQ3 three-way control split). |
| B | Cheng 2026 | reliability-weighted | ✓ OA | [10.64898/2026.07.13.738312](https://doi.org/10.64898/2026.07.13.738312) | `cheng-2026-reliability-weighted` — H1: generalizability-theory variance components with NTC error floor. |
| B | Bianco 2026 | projection-basis | ✓ OA | [10.64898/2026.07.07.737004](https://doi.org/10.64898/2026.07.07.737004) | `bianco-2026-projection-basis` — H3: r² bounded by variance explained by the projection basis. |
| B | Feng et al. 2026 | crispri-map | ◐ UB (preprint ✓ OA) | [10.1016/j.xgen.2025.101076](https://doi.org/10.1016/j.xgen.2025.101076) | `feng-2026-crispri-map` — H2c data: HipSci CRISPRi across 34 iPSC lines / 26 donors — counts toward C_public. |
| B | Nourreddine et al. 2026 | crispri-atlas | ◐ UB (preprint ✓ OA) | [10.1038/s41587-026-03199-w](https://doi.org/10.1038/s41587-026-03199-w) | `nourreddine-2026-crispri-atlas` — H2c data: KOLF2.1J genome-scale CRISPRi (KOLF panel). |
| B | Zhu et al. 2026 | cd4-perturb-seq | ◐ UB (preprint ✓ OA) | [10.1016/j.cell.2026.08.002](https://doi.org/10.1016/j.cell.2026.08.002) | `zhu-2026-cd4-perturb-seq` — H2c data: CD4+ T-cell genome-scale CRISPRi (donor × stimulation contexts). |
| B | Luan et al. 2026 | centered-transfer | ✓ OA | [10.64898/2026.10.20.760514](https://doi.org/10.64898/2026.10.20.760514) | `luan-2026-centered-transfer` — H4: centered KOLF→H1 transfer matches VCC 2025 winner. |
| B | Dhaka et al. 2026 | hypermap | ✓ OA | [10.64898/2026.04.23.720505](https://doi.org/10.64898/2026.04.23.720505) | `dhaka-2026-hypermap` — H4: few-shot seed perturbations for new lines/donors. |
| B | Scholkemper & Mukherjee 2026 | c3tl | ✓ OA | [arXiv:2603.13051](https://arxiv.org/abs/2603.13051) | `scholkemper-2026-c3tl` — H2 H4: C3TL lightweight context transfer. |
| B | Chen et al. 2026 | scale | ✓ OA | [arXiv:2603.17380](https://arxiv.org/abs/2603.17380) | `chen-2026-scale` — H2: SCALE on Tahoe; v2 notes dominant cell-line effects. |
| B | Fu et al. 2026 | strand | ✓ OA | [arXiv:2602.10156](https://arxiv.org/abs/2602.10156) | `fu-2026-strand` — H2: leave-one-line-out; STATE ≈ identity. |
| B | Yeung et al. 2026 | multiply-perturbed | ✓ OA | [10.64898/2026.07.10.737863](https://doi.org/10.64898/2026.07.10.737863) | `yeung-2026-multiply-perturbed` — H2c data: DLD-1 CRISPRi Perturb-seq (GSE337988) — counts toward C_public. |
| C | Hocking 1976 | variable-selection | ◐ UB | [10.2307/2529336](https://doi.org/10.2307/2529336) | `hocking-1976-variable-selection` — H2b: Sp criterion (random-X variable selection). |
| C | Muandet et al. 2013 | invariant-feature | ✓ OA | [arXiv:1301.2115](https://arxiv.org/abs/1301.2115) | `muandet-2013-invariant-feature` — Background: domain generalization via invariant features. |
| C | Peters et al. 2016 | invariant-prediction | ✓ OA | [10.1111/rssb.12167](https://doi.org/10.1111/rssb.12167) | `peters-2016-invariant-prediction` — Background: invariant causal prediction. |
| C | Subramanian et al. 2017 | connectivity-map | ✓ OA | [10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049) | `subramanian-2017-connectivity-map` — Appendix (drug ↔ genetic) + LINCS shRNA source: CMap L1000. |
| C | Pabon et al. 2018 | protein-targets | ✓ OA | [10.1371/journal.pcbi.1006651](https://doi.org/10.1371/journal.pcbi.1006651) | `pabon-2018-protein-targets` — Appendix (drug ↔ genetic): L1000 drug vs knockdown correlation for target ID. |
| C | Lotfollahi et al. 2019 | scgen | ◐ UB | [10.1038/s41592-019-0494-8](https://doi.org/10.1038/s41592-019-0494-8) | `lotfollahi-2019-scgen` — Context: scGen cross-cell-type latent arithmetic. |
| C | Belkin et al. 2019 | bias-variance | ✓ OA | [10.1073/pnas.1903070116](https://doi.org/10.1073/pnas.1903070116) | `belkin-2019-bias-variance` — H2b: double descent, risk peak at n ≈ p (C ≈ d+1 singularity). |
| C | Gonçalves et al. 2020 | drug-moa | ✓ OA | [10.15252/msb.20199405](https://doi.org/10.15252/msb.20199405) | `goncalves-2020-drug-moa` — Appendix (drug ↔ genetic): drug–CRISPR viability associations. |
| C | Kurilov et al. 2020 | modelling-strategies | ✓ OA | [10.1038/s41598-020-59656-2](https://doi.org/10.1038/s41598-020-59656-2) | `kurilov-2020-modelling-strategies` — H2: drug-response R² 0.267 with 329 lines vs 0.102 with 26–68 lines. |
| C | Bühlmann 2020 | invariance-causality | ✓ OA | [10.1214/19-STS721](https://doi.org/10.1214/19-STS721) | `buhlmann-2020-invariance-causality` — Background: invariance and robustness. |
| C | Umarov et al. 2021 | deepcellstate | ✓ OA | [10.1371/journal.pcbi.1009465](https://doi.org/10.1371/journal.pcbi.1009465) | `umarov-2021-deepcellstate` — H2: L1000 cross-cell-type transfer improves with more cell types. |
| C | Ahuja et al. 2021 | information-bottleneck | ✓ OA | [arXiv:2106.06607](https://arxiv.org/abs/2106.06607) | `ahuja-2021-information-bottleneck` — Background: IRM + information bottleneck. |
| C | Hetzel et al. 2022 | chemcpa | ✓ OA | [10.52202/068431-1937](https://doi.org/10.52202/068431-1937) | `hetzel-2022-chemcpa` — Context: chemCPA. |
| C | Eastwood et al. 2022 | quantile-risk | ✓ OA | [arXiv:2207.09944](https://arxiv.org/abs/2207.09944) | `eastwood-2022-quantile-risk` — Background: bounds depending on number of training domains. |
| C | Bunne et al. 2023 | cellot | ✓ OA | [10.1038/s41592-023-01969-x](https://doi.org/10.1038/s41592-023-01969-x) | `bunne-2023-cellot` — Context: CellOT. |
| C | Engler Hart et al. 2023 | transcriptomic-response | ✓ OA | [10.1186/s12859-023-05337-6](https://doi.org/10.1186/s12859-023-05337-6) | `englerhart-2023-transcriptomic-response` — Appendix (drug ↔ genetic): weak drug–target transcriptomic concordance. |
| C | Piran et al. 2024 | biolord | ✓ OA | [10.1038/s41587-023-02079-x](https://doi.org/10.1038/s41587-023-02079-x) | `piran-2024-biolord` — Context: biolord disentanglement. |
| C | Bendidi et al. 2024 | transcriptomics-foundation | ✓ OA | [arXiv:2410.13956](https://arxiv.org/abs/2410.13956) | `bendidi-2024-transcriptomics-foundation` — H4: "one PCA still rules them all". |
| C | Chandrasekaran et al. 2024 | morphological-profiles | ✓ OA | [10.1038/s41592-024-02241-6](https://doi.org/10.1038/s41592-024-02241-6) | `chandrasekaran-2024-morphological-profiles` — Appendix (drug ↔ genetic): matched compound/CRISPR morphology (CPJUMP1). |
| C | Peidli et al. 2024 | scperturb | ◐ UB (preprint ✓ OA) | [10.1038/s41592-023-02144-y](https://doi.org/10.1038/s41592-023-02144-y) | `peidli-2024-scperturb` — Data resource: scPerturb (RN panels loaded from it). |
| C | Kernfeld et al. 2025 | expression-forecasting | ✓ OA | [10.1186/s13059-025-03840-y](https://doi.org/10.1186/s13059-025-03840-y) | `kernfeld-2025-expression-forecasting` — H4: GGRN/PEREGGRN expression forecasting. |
| C | Miladinovic et al. 2025 | perturbation-models | ✓ OA | [10.1038/s43588-025-00870-1](https://doi.org/10.1038/s43588-025-00870-1) | `miladinovic-2025-perturbation-models` — Appendix (drug ↔ genetic): inhibitors co-embed with target CRISPR. |
| C | Jiang et al. 2025 | pathway-signatures | ◐ UB (preprint ✓ OA) | [10.1038/s41556-025-01622-z](https://doi.org/10.1038/s41556-025-01622-z) | `jiang-2025-pathway-signatures` — Data ("Jiang CRISPRi" in plan): multi-line × cytokine Perturb-seq — CRISPR modality (KO vs CRISPRi) UNVERIFIED. |
| C | DeMeo et al. 2025 | active-learning | ◐ UB | [10.1126/science.adi8577](https://doi.org/10.1126/science.adi8577) | `demeo-2025-active-learning` — Data: GSE306429 (104 compounds × 10 lines). |
| C | Dip & Zhang 2026 | stable-shift | ✓ OA | [10.1145/3807503.3820871](https://doi.org/10.1145/3807503.3820871) | `dip-2026-stable-shift` — Context: Stable-Shift, unseen genes in K562 only. |
| C | Sun et al. 2026 | de-classification | ✓ OA | [10.64898/2026.07.20.739620](https://doi.org/10.64898/2026.07.20.739620) | `sun-2026-de-classification` — Context: cross-line DE classification; STATE runs few-shot. |
| C | Fullmer et al. 2026 | additive-input | ✓ OA | [10.64898/2026.09.30.755776](https://doi.org/10.64898/2026.09.30.755776) | `fullmer-2026-additive-input` — H1: held-out error stays at predict-zero. |
| C | Gallone et al. 2026 | confidence-calibrated | ✓ OA | [10.64898/2026.10.12.758706](https://doi.org/10.64898/2026.10.12.758706) | `gallone-2026-confidence-calibrated` — H1: identity-predictor echo ceiling leakage. |
| C | Wei et al. 2026 | cancer-heterogeneity | ✓ OA | [10.64898/2026.08.10.743942](https://doi.org/10.64898/2026.08.10.743942) | `wei-2026-cancer-heterogeneity` — H2: 162 lines; cancer-line heterogeneity as bottleneck. |
| C | Wang et al. 2026 | x-cell | ✓ OA | [10.64898/2026.03.18.712807](https://doi.org/10.64898/2026.03.18.712807) | `wang-2026-x-cell` — H2: X-Cell; "perturbation-context combinations are the binding constraint"; no C sweep. |
| C | Xing & Song 2026 | world-model | ◐ UB | [10.1016/j.cell.2026.08.042](https://doi.org/10.1016/j.cell.2026.08.042) | `xing-2026-world-model` — Framing: Cell commentary on virtual-cell world models. |
| C | Li et al. 2026 | unipert | ◐ UB | [10.1016/j.cell.2026.06.005](https://doi.org/10.1016/j.cell.2026.06.005) | `li-2026-unipert` — Appendix (drug ↔ genetic): UniPert-G2CP bridges genetic and chemical screens. |

**Legend.**
- **Grade:** A core / must-read · B supporting · C peripheral / context.
- **OA/Zugang:** ● Bibliothek (already in the library — none yet, the library is empty) · ✓ OA (auto-downloadable;
  bioRxiv/arXiv/CC-BY) · ◐ UB (licence needed; "preprint ✓ OA" means a legal preprint exists, but ingest needs the
  version of record) · ○ Fernleihe · ✗.
- **Plan-listed works** (the plan's A set) are graded A whatever their topical weight.

## 2. Research question and how the literature bears on it

**RQ.** Is zero-shot CRISPRi response prediction in unseen cell lines limited by data rather than model capacity? The
three candidate limits are (i) a measurement-noise ceiling, (ii) a non-transferable line × perturbation interaction,
and (iii) too few contexts to learn that interaction (Proposition 1: rank ≤ C−1; Proposition 2:
C* = d + 1 + dσ²/(R²Var(i))).

- **H1 (ceiling).**
  - Ingredients: Molina & Zhang 2026, the four-way split-half decomposition. A ceiling of ≈0.84 is implied by their
    Fig. 1c; see `input/notes/overlap-molina-zhang.md`.
  - Reliability-ceiling methods: Wang et al. 2026 (√ρ); TxPert (split-half / resampled); Mejía et al. 2025
    (technical duplicate); Li et al. 2026 (control noise floor).
  - Caveats: Huang 2026 (shared clones inflate random-half ceilings; small effect in CRISPRi lines) and Nicol et al.
    2026 (shared-control reuse).
  - **No paper** converts the decomposition into a per-destination zero-shot transfer ceiling and places real methods
    on it.
- **H2a/H2b (context break-even, Tahoe).**
  - Empirical precedent: Prophet fig. S4 (STOP-FLAG above); Dibaeinia Fig. 4; Molina Fig. 3a; Wei NM Suppl. Note 11;
    Shen 2026 (C = 46, no gain).
  - Theory: random-X OLS risk (Breiman & Freedman; Rosset & Tibshirani; Hastie et al.; Dicker; Copas); environment
    counts in domain generalization (Rosenfeld; Chen; Wang; Arjovsky; Blanchard); task diversity in meta-learning
    (Tripuraneni ×2; Du; Maurer; Baxter; Kong; Lucas).
- **H2c (genetic data below break-even).** C_public registry sources:
  - Replogle 2022, Nadig 2025, X-Atlas/Orion 2025, Pan 2026 (GSE295214, hESC).
  - KOLF2.1J (Nourreddine 2026); HipSci, 34 iPSC lines (Feng 2026) — decide how donor lines count toward C_public.
  - CD4 T cells (Zhu 2026); DLD-1 (Yeung 2026, GSE337988).
  - Random-effects view of few groups: Maas & Hox 2005; McNeish & Stapleton 2016.
- **H3 (residual on destination control PCs + correction comparison).**
  - Precedent: Molina & Zhang Supp. Note 3 projects the *response* onto the line's *own* control PCs against a random
    basis only. **Verdict: RQ3 keep as planned** (see overlap note).
  - Related: Bianco 2026 (basis ceiling); Systema (common-response removal, an H3b alternative); Chen et al. 2026
    (control-only regime fails); Huang et al. 2026 (context geometry compression).
- **H4 (same-line data vs modelling).**
  - Molina & Zhang: 30 % target perturbations recover γ; calibration saturates at 100–150 perturbations.
  - STATE: the 30 %-target protocol is few-shot. Also HyperMap, Luan 2026 (centered KOLF → H1 transfer ≈ VCC
    winner), and the baselines literature (Ahlmann-Eltze; Csendes; PertEval; Bendidi; Kernfeld; Miller 2026).

## 3. Primary sources / data

- Replogle 2022 (K562 GW + K562/RPE1 essential) and Nadig 2025 (HepG2, Jurkat).
- X-Atlas/Orion 2025 (HCT116 and others); Tahoe-100M (Cell 2026; 50 lines).
- Pan 2026: GSE295214, genome-scale CRISPRi in hESC. GEO brief says "hESCs"; confirm the line is H1 in the sample
  records.
- KOLF2.1J (Nourreddine 2026), HipSci (Feng 2026), CD4 (Zhu 2026), DLD-1 (Yeung 2026, GSE337988).
- GSE306429 = DeMeo et al. 2025 Science (104 compounds × 10 lines; chemical, not CRISPRi).
- LINCS shRNA = CMap L1000 (Subramanian 2017).
- "Jiang CRISPRi" was resolved provisionally to Jiang et al. 2025 NCB. **Its CRISPR modality is unverified (it may be
  Cas9 KO)**; check before counting it in C_public.
- Challenge framing: Roohani 2025 (VCC 2025) and Adduri 2026 (VCC 2026: zero-shot across unseen lines).

## 4. Debates

1. **Do learned models beat simple baselines?** No: Ahlmann-Eltze; Csendes; Wei NM; Li SciAdv; Bendidi; PertEval. Yes,
   on calibrated metrics: Miller 2026 NBT; Mejía 2025.
2. **Scale vs context diversity.**
   - For diversity: Dibaeinia 2026; Wang 2026 X-Cell ("binding constraint").
   - Similarity, not count: Wei NM 2025.
   - Our C* quantifies "how much diversity is enough", which Dibaeinia explicitly leave open.
3. **Is the cross-line residual signal or noise?**
   - Molina & Zhang: the interaction is reproducible (23.5 %) and learnable with target data.
   - Shen 2026 (Tahoe): the context residual has median reliability 0.067.
   - Prophet: SCORE fails because of low replicate correlation.
   - Our framework separates the two: σ²_ε/Var(i) enters C*_min directly.
4. **Does control expression carry context information?** Negative: Molina & Zhang (descriptor ablations, Supp.
   Tables 24/26); Chen et al. 2026 (control-only regime); Shen (basal PCs). Hopeful: PIE, STATE, CPA-type embeddings.
   H3 tests the residual side, a different question.

## 5. Methods worth borrowing

- Split-half cross-product estimator: Molina §4.8.4, 50 splits; adopt and cite.
- Reliability ceiling √ρ: Wang 2026.
- Guide-split rather than random-half splits: Huang 2026.
- Split controls: Nicol 2026.
- Perturbed-reference Pearson: Systema; Molina §4.2.3.
- Random-X risk algebra: Rosset & Tibshirani; Dicker.
- Rank / span arguments: Wang et al. 2022 (ISR-Mean); Chen et al. 2022 lower-bound construction.

## 6. Open access

- 79 of 102 entries are ✓ OA (bioRxiv, arXiv, CC-BY); 23 are ◐ UB.
- Paywalled versions of record with OA preprints: STATE, Tahoe, Replogle, Nadig, Zhu, Feng, KOLF, scGPT, scPerturb,
  Jiang.
- No OA copy located: Breiman & Freedman, Rosset & Tibshirani, Copas, Hocking, Maas & Hox, McNeish, Wei NM, VCC 2025
  and VCC 2026 commentaries, Xing 2026, DeMeo 2025, scGen, UniPert.
- `acquire-sources` (Task 3) handles all of these. Nothing was downloaded here.

## 7. Gaps (what this preprint fills)

1. A context-count break-even C* with a closed form, predicted from independently estimated parameters (H2b). The
   precedents are only an empirical share-of-lines curve (Prophet) and qualitative calls for diversity.
2. A model-free lower bound C*_min for public genetic screens, against a counted C_public (H2c).
3. A per-destination zero-shot transfer ceiling with real methods placed on it (H1).
4. Transfer-prediction residual vs other lines' / permuted control PCs under split-half control, plus a correction
   comparison (H3).
5. Same-line data vs modelling gains on one official scorer (H4).

## 8. Recommended reading order

1. Molina & Zhang 2026 (overlap note first).
2. Prophet fig. S4 / S5 (manual figure check).
3. Shen 2026.
4. Dibaeinia 2026.
5. Wang et al. 2022 ISR and Rosenfeld et al. 2021.
6. Rosset & Tibshirani 2020 and Hastie et al. 2022.
7. Wang et al. 2026 (reliability) and Huang 2026 (shared clones).
8. STATE (Cell) and PIE.
9. Li SciAdv 2026, Wei NM 2025.
10. Systema; Ahlmann-Eltze.
11. Data papers as needed for Task 4.

## 9. Follow-up searches / unresolved

- **Prophet fig. S4A** crossing points (manual; bioRxiv blocks scripted access). Also check whether Prophet has a
  journal version; Crossref shows only the preprint.
- **Not added to the bib (unverifiable):**
  - PertResolve (Li B. et al., "Distinguishing measurement and model limitations in perturbation prediction"): only
    a GitHub-hosted manuscript, no DOI.
  - Stein 1960, "Multiple regression" (Hotelling Festschrift): no DOI; cite from the book if needed.
  - Denevi, Pontil & Ciliberto, NeurIPS 2020 (conditional meta-learning): arXiv id unverified.
  - Hodos et al., PSB 2018: DOI unverified.
- Rosenfeld et al. 2021 / Chen et al. 2022 / Dwork et al. 2025 / Lucas et al. 2021: confirm proceedings records
  (OpenReview / NeurIPS) before submission.
- Huu Loc & Van 2026: name order as deposited in Crossref; confirm on the PDF.
- C_public counting rule: do HipSci's 34 donor iPSC lines count as 34 contexts or one "iPSC" context? Decide in Task 4
  and log it. This changes H2c materially.
- Semantic Scholar was mostly unreachable (bot check / timeouts) and OpenAlex was not used (shared budget). A second
  pass with an OpenAlex API key could widen the theory side (e.g. "environment complexity", "number of training
  domains").

## 10. BibTeX overview

`output/bibtex/references.bib`: 102 entries. The grade is in `keywords = {grade-X}`. Provisional plan keys were renamed
to the schema shape; renames are logged in `knowledge/_meta/log.md`:

| Plan key | Bibkey |
|---|---|
| molina-zhang-2026 | molina-2026-response-decomposition |
| cell-quota-2026 | huuloc-2026-cell-quota |
| sciadv-benchmark-2026 | li-2026-systematic-comparison |
| dibaeinia-2026 | dibaeinia-2026-context-scale |
| ahlmann-eltze-2025 | ahlmanneltze-2025-linear-baselines |
| adapert-2026 | piao-2026-adaptive-contexts |
| state-2025 | adduri-2026-state |
| tahoe-2025 | zhang-2026-tahoe |
| systema-2025 | vinastorne-2026-systema |
| vcc-2026 | adduri-2026-virtual-cell-challenge |
| replogle-2022 | replogle-2022-genome-scale |
| nadig-2025 | nadig-2025-transcriptome-wide |
| orion-2025 | huang-2025-x-atlas-orion |
| pan-2026 | pan-2026-global-cell-state |

Audit log: `input/bibliography/audit-log-2026-11-05.json`.
