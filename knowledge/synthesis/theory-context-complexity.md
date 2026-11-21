---
title: "Theory: identifiability and the break-even number of contexts for learned context corrections"
type: synthesis
created: 2026-11-05
updated: 2026-11-06
status: stable
author: mixed
---

# Theory: how many cell contexts does a learned context correction need?

Plan: `input/ideas/zero-shot-limits-context-complexity-plan.md` (Notation block, H2a–c, Amendments A1–A3).
Empirical counterparts: [[rq2_h2a-assumptions]] (Tahoe curve), [[rq2_params-assumptions]] (Tahoe parameters),
[[rq2_gene_params-assumptions]] (genetic parameters). Concepts: [[context-count-learning-curve]],
[[perturbation-response-decomposition]], [[zero-shot-transfer-ceiling]].

## 1. Setup

For cell context $c$ and perturbation target $t$ the pseudobulk response (log1p CP10K, gene vector in $\mathbb{R}^G$) is
decomposed as in [[molina-2026-response-decomposition]]:
$$\delta_{c,t} = \mu + a_t + b_c + i_{c,t} + \varepsilon_{c,t},$$
with a shared response $\mu$, a conserved target component $a_t$, a context component $b_c$, a context × target
interaction $i_{c,t}$ and measurement noise $\varepsilon_{c,t}$ (variance $\sigma^2_\varepsilon$ per element at full depth).
In zero-shot prediction for an unseen context $c^\ast$, the only information about $c^\ast$ is its control-cell profile,
summarised as a descriptor $x_{c^\ast}\in\mathbb{R}^p$. Everything context-specific ($b_{c^\ast}$, $i_{c^\ast,t}$) must
be predicted from $x_{c^\ast}$ through a map learned on $C$ source contexts. The number of independent samples for
that map is $C$, regardless of how many cells or perturbations each context contributes.

## 2. Proposition 1 (identifiability; cf. multi-environment results)

Let the learned correction be linear in the descriptor, $\hat i_{c,t} = W_t x_c$, and let the source descriptors be
centred on their training mean, $X=[x_1,\dots,x_C]$, so that $\operatorname{rank}(X)\le C-1$. Any estimator that sees
only $\{(x_c, i_{c,t})\}_{c\le C}$ observes $W_t$ only through $W_tX$: two maps with $(W_t-W'_t)X=0$ generate identical
training data but predict differently for any $x^\ast$ with a component outside $\operatorname{span}(X)$, the difference
being $(W_t-W'_t)P^\perp_X x^\ast$. Hence, over the class $\{\|W_t\|_F\le B\}$, the worst-case error of every estimator at
$x^\ast$ is at least of order $B^2\|P^\perp_X x^\ast\|^2$. At most $C-1$ descriptor directions are identifiable,
whatever the model capacity. For a nonlinear map the same statement holds pointwise (the function is constrained only at
$C$ points), and under an $s$-smoothness assumption the minimax rate is the nonparametric rate with sample size $C$,
$C^{-2s/(2s+p_{\text{eff}})}$ [@stone-1982-optimal-rates].

This is the context-count analogue of the environment-count results of domain-generalisation theory: invariant risk
minimisation is optimal only with more training environments than the dimension of the environmental feature, in a
linear, infinite-sample, binary setting without a noise term ([[rosenfeld-2021-risks-irm]], Thm 5.1, Cor. 5.2;
Thm 6.1 for failure in the nonlinear case). We therefore cite Proposition 1 as "cf.", not as a new theorem.

## 3. Proposition 2 (break-even, linear-Gaussian OLS) — tested and REJECTED quantitatively

Assume a Gaussian random design over contexts, a $d$-dimensional descriptor (effective dimension = participation ratio
of the centred descriptor spectrum), a per-element linear model with explained fraction $R^2$ of the noise-free
interaction variance $\operatorname{Var}(i)$, and residual variance $\sigma^2_{\text{res}}=(1-R^2)\operatorname{Var}(i)+\sigma^2_\varepsilon$.
With $C>d+1$ contexts the expected out-of-sample excess risk of OLS is $\sigma^2_{\text{res}}\,d/(C-d-1)$
[@breiman-1983-how-many-variables; @rosset-2020-random-x; @dicker-2013-equivariant-prediction]. Against the
context-free predictor (predict $i=0$) the expected gain per element is
$$\operatorname{gain}(C)=R^2\operatorname{Var}(i)-\sigma^2_{\text{res}}\frac{d}{C-d-1},\qquad
C^\ast = d+1+\frac{d\,\sigma^2_{\text{res}}}{R^2\operatorname{Var}(i)},\qquad
C^\ast_{\min}=d+1+\frac{d\,\sigma^2_\varepsilon}{\operatorname{Var}(i)}\ (R^2=1).$$

**Empirical test (pre-registered H2b).** On Tahoe (48 lines, 185 drugs), parameters estimated blind to the learning
curve gave $d=20.6$, noise/interaction ratio 2.87, $R^2_{\text{adj}}=0.35$, hence $C^\ast_{\text{theory}}=230$
[192, 287] and $C^\ast_{\min}=80.6$ [68.5, 97.6]. The measured break-even of the context pathway of a FiLM model is
$C^\ast_{\text{emp}}=37.4$ [33.3, > 40] (H2a confirmed). The ratio 6.2 falls outside the registered band [0.5, 2]:
**the quantitative form of Proposition 2 is rejected**, and even its "lower bound" $C^\ast_{\min}$ exceeds the measured
break-even (×2.2). Pre-specified sensitivities do not rescue it: using the model's actual target $b+i$ gives
$C^\ast_{\text{theory}}=192$, $C^\ast_{\min}=65$.

**Why the per-gene OLS bound is pessimistic (interpretation, not tested):** (i) the learner is regularised (weight
decay, early stopping at the first checkpoint in ≥ 94 % of runs), so it never pays the OLS variance blow-up near
$C\approx d$ [@hastie-2022-ridgeless-interpolation; @copas-1983-regression-shrinkage]; (ii) it predicts through a
rank-32 response basis shared across genes and drugs, so its effective noise per parameter is far below the per-gene
$\sigma^2_\varepsilon$; (iii) its target includes $b_c$, one vector per line shared by all drugs, which is cheap to
learn ($C^\ast_b$ = 6–22 depending on the noise model). The surviving, model-free content is Proposition 1 plus the
measured curve: the context pathway is harmful below ~24 contexts and breaks even near 37.

## 4. Consequence for public genetic data (H2c)

With 11 distinct public cell types (Count 1) the genetic analogue gives $C^\ast_{\min}=16.45$ [16.09, 16.82]
($d_{\text{gene}}=5.1$, capped at $C-1$; noise/interaction 2.04) — above 11 but below the 25 contexts obtained when
HipSci donor lines are counted separately (Count 2). Verdict: "holds for distinct cell types, not for donor lines". Because
the same estimator overestimated the Tahoe break-even ×2.2 and its line-level bootstrap spans [7.95, 17.07], the
genetic comparison is reported as weak support only; any drug-calibrated genetic break-even is exploratory.

## 5. What the paper claims and does not claim

- Claims: the effective sample size for cross-context generalisation is the number of contexts (Prop. 1, cited as cf.);
  a matched-ablation learning curve on 48 drug-treated lines shows harm below ~24 contexts and break-even at ≈ 37;
  public genetic Perturb-seq covers ~11 distinct cell types.
- Does not claim: a validated closed-form $C^\ast$; a valid noise lower bound; that drug-side break-even transfers
  quantitatively to genetic perturbations.

## Evidence base

- Tahoe curve and break-even: `output/data-analysis/results/rq2_h2a.json`, `rq2_tahoe_curve.parquet` ([[rq2_h2a-assumptions]]); 3 seeds × 6 folds × 7 C values; cluster bootstrap.
- Tahoe parameters (blind to the curve): `rq2_params_tahoe.json` ([[rq2_params-assumptions]]), independently reviewed (PASS-WITH-NOTES).
- Genetic parameters: `rq2_params_gene.json`, `rq2_h2c.json` ([[rq2_gene_params-assumptions]]).
- Literature: [[rosenfeld-2021-risks-irm]] (environment-count identifiability), [[ji-2024-universal-prediction]] (empirical curves on scalar readouts), regression-risk results [@breiman-1983-how-many-variables; @rosset-2020-random-x; @dicker-2013-equivariant-prediction; @hastie-2022-ridgeless-interpolation; @stone-1982-optimal-rates].
- Gap: one learner family (FiLM r32) and one drug dataset (a 4.43 M-cell Tahoe subset) carry the empirical break-even.

## What stays open

- Whether the ≈ 37-context break-even holds for other learners (stronger regularisation, larger models) and for genetic perturbations directly — no public genetic dataset has enough contexts to measure it.
- A learner-aware theory (ridge / low-rank shared basis, target b + i) that predicts the measured break-even; the per-gene OLS form is ×6 too pessimistic.
- How to count donor-line contexts (Count 2) in a way that reflects descriptor diversity rather than line labels.
