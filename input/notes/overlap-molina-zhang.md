---
title: "Overlap gate — Molina & Zhang 2026 vs RQ1/RQ2/RQ3"
type: note
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
plan: zero-shot-limits-context-complexity-plan
bibkey: molina-2026-response-decomposition
---

# Overlap gate: Molina & Zhang (2026) vs this preprint

**Source read.** Molina, A. & Zhang, X. (2026). *Perturbation response decomposition enables biologically aligned
generalization to unseen perturbations and cellular contexts.* bioRxiv, posted 2026-07-27,
doi:10.64898/2026.07.24.740459 (CC BY 4.0; metadata verified against Crossref 2026-11-05). Code:
github.com/xinyizhanglab/perturbation-decomposition.

**What was read.** The bioRxiv `v1.full` HTML (fetched via Exa, 182 k characters, 2026-11-05): abstract, §1–3, Methods
§4.1–4.11 (incl. §4.5.1–4.5.4, §4.7, §4.8.1–4.8.4, §4.9.1–4.9.7, §4.10.1–4.10.4), figure legends Fig. 1–9, and
Supplementary Notes 1–16. Not read: §4.6 (Norman combinatorial benchmark — irrelevant here) in detail, and the
supplementary *tables* themselves (the HTML gives table titles only; their numbers are quoted in the Notes). Equations
are images in the HTML and are missing from the text; definitions below are taken from the surrounding prose. Direct
curl to biorxiv.org is blocked (Cloudflare 1015); the PDF goes to the library in Task 3.

**Their data.** Four CRISPRi Perturb-seq screens: K562, RPE1 (Replogle 2022), HepG2, Jurkat (Nadig 2025); 2,000 HVGs
from 6,642 shared genes; balanced tensor of **638 perturbations × 4 lines** (§4.8.2). Tahoe-100M: 5 lines × 268
compounds × 3 doses (§4.7.1). Notation map: their μ(g), α_c(g), β_p(g), γ(c,p,g) = our μ, b_c, a_t, i_{c,t}.

---

## (i) RQ1 — variance decomposition, noise correction, ceiling

**What they show (overlaps directly with the *ingredients* of H1):**

- Split-half, noise-corrected four-way decomposition. §4.8.4: *"Within each cell-line-by-perturbation condition, cells
  were randomly divided into two non-overlapping subsets … reproducible signal was estimated using the cross-half
  Frobenius inner product"*; 50 random splits (§4.11). This is the same estimator as our operational definition 4
  (cross-product ⟨i_A, i_B⟩).
- Numbers (§2.1, Fig. 1c, Supp. Note 1 / Supp. Table 1): *"reproducible variance was distributed across the
  cell-line-specific component (the cell-line template; 27.8%), the conserved perturbation-specific component (29.4%),
  and the cell-line-by-perturbation interaction (23.5%), with the remaining 19.3% attributable to measurement noise."*
  (Supp. Note 1 gives 27.9 / 29.4 / 23.5 / 19.2.) After projective template removal: conserved 35.0 %, interaction
  27.7 %, noise 37.4 % (Supp. Note 1).
- They state the interpretation that sits under our ceiling: §4.8.4 *"The reproducible portion of γ c,p,g represents
  target-specific interaction structure that is consistently detected across independent cell subsets and is
  therefore, in principle, available for model recovery."*
- Zero-shot attribution: Fig. 4b / Supp. Note 16 (Table 25): DepMap Ridge/MLP *"aligned positively with the conserved
  component but not with the target interaction"*; *"correlations with γ(c, p, g) were negative across all targets."*
  STATE/MORPH recover neither β nor γ zero-shot (Fig. 4c, Supp. Note 9: perturbed-reference r 0.046 / −0.013 K562→RPE1).
- Upper *references* exist, but they are empirical, not ceilings: Fig. 8e / Supp. Note 14 (*"full within cell line
  training provided an upper reference for matched context prediction"*), and target calibration saturating after
  ~100–150 target perturbations (Supp. Note 15).

**What they do NOT do (H1 stays new):**

- No zero-shot **transfer ceiling**: they never convert the decomposition into the best correlation achievable by an
  oracle that knows a_t and b_c but not i (our r_max = √[(Var(a)+Var(b))/Var_reproducible]). No method is placed as a
  fraction of such a bound, and no share of remaining error is assigned to transferable vs non-transferable components.
- Their decomposition is on 4 lines × 638 essential-gene perturbations; ours is per H1 panel (H1 hESC destination, 7
  harness panels), at the VCC scorer's gene space.
- **Caveat for the paper:** a reader can compute an implied ceiling from their Fig. 1c numbers:
  (27.8 + 29.4)/80.7 ≈ 0.71 of reproducible variance is transferable → r_max ≈ 0.84 on their 4-line tensor. We must
  cite this explicitly as "implied by Molina & Zhang's decomposition" and present our contribution as (a) making the
  ceiling operational and per-destination, (b) placing real methods (our frozen pipeline + 17 architectures) on it.

## (ii) RQ3 — what is projected onto control PCs, and against which null

**What they project: the RESPONSE (template and template-removed response residual), not a transfer-prediction
residual.**

- Supp. Note 3 (*"Projection of perturbation responses onto control expression subspaces"*) and Fig. 1i / Fig. 5e:
  *"we fitted PCA on centered non-targeting control cells from each cell line and projected the template and residual
  responses onto the resulting control subspace."* Here "residual" = the observed response minus its projection on
  the same line's template (§4.9.4: *"The template-removed residual was defined as …"*), i.e. a property of the truth,
  not of any model's error.
- Numbers (Supp. Note 3): first 100 control PCs capture 92.0 / 98.3 / 96.8 / 95.2 % of the template (K562 / RPE1 /
  HepG2 / Jurkat) but only 36.9 / 34.5 / 39.3 / 38.6 % of residual variance; at 5 dimensions control PCs capture 13.4 /
  7.1 / 11.2 / 11.9 % of residual variance vs 34.4 / 32.1 / 34.9 / 31.2 % for residual PCs.
- The control subspace is **the same line's own** NT controls throughout. In the genome-wide K562 analysis (§4.9.6)
  the "control" subspace is built from **non-significant perturbations** (energy-test p ≥ 0.5) as a control proxy, with
  a Marchenko–Pastur threshold for noise.

**Which null: random bases, not other lines' controls, not permuted controls.**

- Supp. Note 3: *"We then compared control PCs, perturbation PCs and random bases … Random bases captured less than
  0.5% in all cases"* (≈0.3 % after template removal). Their only null is a random orthonormal basis.
- **No comparison against other cell lines' control PCs** anywhere in the paper. The only cross-line subspace analysis is
  source **response** basis vs target response residual (§4.5.4, Fig. 3c, Supp. Note 13: source bases capture 20–30 %
  of target residual variance at rank 20) — response-on-response, not residual-on-control.
- **No gene-wise permuted-control null**, and **no split-half noise control on the projection fractions** (the
  projected residual energy includes measurement noise; noise control is applied only to the four-way decomposition).
- Related but different: control expression as a model **input** does not help (Supp. Note 16, Tables 24 & 26:
  *"Concatenating target or source control descriptors with DepMap left perturbed reference performance and component
  alignment essentially unchanged"*; target-control fine-tuning of MORPH does not recover γ, Fig. 4c). This is a
  statement about descriptors, not about where a transfer residual lives.
- No PC-deflation **correction** is proposed or compared (our H3b: M3 vs other-line / permuted / prediction-matrix PC
  deflation / Systema-style removal / none).

**Tension to address in the paper.** Their headline reading (§2.3) is *"the top PCs of the control cells do not
capture variance beyond cell-line template"*. Our H3a claims the **transfer residual** concentrates on destination
control PCs more than on other lines' PCs. These are compatible (template ≈ in control space; a transfer residual is
dominated by the mis-predicted b_c/template part), but the manuscript must say so explicitly and report our k = 8
fractions next to their 5- and 100-PC numbers.

## (iii) RQ2 — number of contexts

- The only context-count manipulation is **1 vs 3 source lines** (Fig. 3a; §4.5.1; Supp. Notes 10–11):
  *"Incorporation of more source cell lines improved generalization performance for both models and in all held-out
  cell lines"*; for MLP *"the improvement from one to three source cell lines was significant in every target cell line
  (… p < 10^−10 for all comparisons)."*
- This is a two-point curve for a **context-free** model: §4.5.1 *"Ridge and MLP received no explicit cell-line
  descriptor and consequently estimated a source-averaged response mapping."* More sources improve the estimate of the
  conserved component β, which is exactly what our Proposition 2 predicts for the context-free arm. There is **no**
  learned context correction, **no** gain-vs-C curve against a context-free baseline, **no** zero crossing, **no**
  break-even, **no** sample-complexity bound in C.
- Sample size in **perturbations** (not contexts) for target calibration: Supp. Note 15 — *"Target calibration
  saturates after approximately 100 to 150 perturbations."* Relevant as a contrast (few-shot within target vs zero-shot
  across contexts), and for H4 (same-line data vs modelling).
- Their claim *"the cell-line-by-perturbation interaction is learnable from currently available data"* (§2.4, Fig. 4d,
  Supp. Fig. 5f) refers to **target-calibrated** training with 30 % of target perturbations, not zero-shot. Cite with
  this qualification; it supports H4 (target data beat modelling) rather than contradicting H2.

**No STOP-FLAG from this paper.**

---

## Verdict: **RQ3 keep as planned**

Reasons:
1. Molina & Zhang project the observed **response** (template / template-removed residual) onto the **same line's**
   control PCs; they never project a **transfer-prediction residual** (truth − pre-M3 prediction).
2. Their null is a **random basis**; they compare against neither **other cell lines' control PCs** nor **gene-wise
   permuted controls** — the two nulls of H3a.
3. No split-half noise control is applied to the projection fractions; H3a's three-way control split (C_pc / C_a / C_b)
   and cross-half estimator F(U) are new.
4. No correction comparison exists (H3b is entirely new).

The plan's shrink condition (*"If it already compares transfer residuals against other lines' control axes, shrink RQ3
to the correction comparison (H3b)"*) is **not met**. Required citations/framing: H3a's methods paragraph must cite
Supp. Note 3 / Fig. 1i / Fig. 5e as the response-level precedent and state the three differences above; RQ1's ceiling
must cite Fig. 1c / Supp. Note 1 as the source of the variance components and acknowledge the implied ≈0.84 ceiling.
