---
title: "Rosenfeld, Ravikumar & Risteski 2021 — The Risks of Invariant Risk Minimization"
type: source
created: 2026-11-05
updated: 2026-11-05
status: review
author: llm
bibkey: rosenfeld-2021-risks-irm
depth: standard
tags: [domain-generalization, invariant-risk-minimization, environment-count, identifiability, linear-gaussian, theory, context-count, proposition-1, h2b]
relations:
  - target: context-count-learning-curve
    type: supports
    confidence: extracted
    because: "Supplies the theoretical floor and a synthetic count curve against a context-free (ERM) comparator: «up until E = de , IRM essentially matches ERM in performance … As soon as we cross that threshold of observed environments, the predictor learned via IRM begins to perform drastically better under distribution shift» (p. 17; Fig. C.1, p. 18)."
---

# The Risks of Invariant Risk Minimization

## Bibliographic Details

Rosenfeld, Elan, Pradeep Ravikumar, and Andrej Risteski. 2021. "The Risks of Invariant Risk Minimization." In
*International Conference on Learning Representations (ICLR 2021)*. arXiv:2010.05761v2 (27 Mar 2021),
doi:10.48550/arXiv.2010.05761. Carnegie Mellon University, Machine Learning Department.
Read from `<library>/pdf/rosenfeld-2021-risks-irm.pdf` (arXiv v2, 27 pp.; printed page numbers = PDF page numbers;
main text pp. 1–9, references pp. 10–12, appendices A–F pp. 13–27). Citekey: @rosenfeld-2021-risks-irm.

A pure-theory paper (plus one synthetic experiment). It contains no data set, no regression of a continuous response
and no biology. Theorem / lemma numbers below are the paper's own.

<!-- ═══ ERSCHLIESSUNG ═══ -->

## Section map

| Covered | Section | pp. |
|---|---|---|
| ✓ F1 | Abstract; 1 Introduction | 1–2 |
| ✓ F1 | 2 Related work (domain-invariant vs conditionally invariant representations; concurrent results) | 2–3 |
| ✓ F1 | 3 Model and informal results (eqs. 1–3, Fig. 3.1, Definition 1, Theorems 3.1–3.3 informal) | 3–5 |
| ✓ F1 | 4 Background on IRM and its alternatives (eqs. 4–5, Proposition 4.1) | 5 |
| ✓ F1 | 5 The difficulties of IRM in the linear regime (Theorem 5.1, Corollary 5.2, Arjovsky comparison, Theorem 5.3) | 6–7 |
| ✓ F1 | 6 The failure of IRM in the non-linear regime (Definition 2, Assumption 1, Theorem 6.1) | 7–9 |
| ✓ F1 | 7 Conclusion | 9 |
| — | References | 10–12 |
| — | A Notation; B Proof of Proposition 4.1 | 13 |
| ✓ F1 | C.1 Proof of Theorem 5.1 (non-degeneracy conditions 9–10; Lemmas C.1–C.3) | 13–17 |
| ✓ F1 | C.2 Experiments for Theorem 5.1 (Fig. C.1) | 17–18 |
| ✓ F1 | C.3 Proof of Theorem 5.3; C.4 simulations of μ̃ (Fig. C.2) | 18–19 |
| ✓ F1 | D.1 Proof of Theorem 6.1 (Lemmas D.1–D.2, Theorem D.3 full version) | 19–23 |
| ✓ F1 | D.2 Discussion of conditions and Assumption 1 | 23 |
| ✓ F1 | E Extensions to alternative objectives (Corollaries E.1–E.2; Fig. E.1) | 23–25 |
| — | F Technical lemmas (F.1 strict-subset risk; F.2 optimal coefficient 2Σ̄⁺μ̄; F.3–F.5 Gaussian tails) | 25–27 |

## Kernthesen

### 1. ⭐⭐ The data model: a linear-Gaussian latent SEM with a fixed invariant block and an environment-varying block (pp. 3–4; eqs. 1–3; Fig. 3.1)

Data come from E training environments ℰ = {e₁, …, e_E}, and «we know from which environment each sample is drawn»
(p. 3). In every environment a binary label is drawn, y = 1 w.p. η, −1 otherwise (eq. 1); then
«zc ∼ N(y · µc, σc² I), ze ∼ N(y · µe, σe² I)» with µc ∈ R^{d_c}, µe ∈ R^{d_e} (eq. 2, p. 3), and the observation is
x = f(zc, ze) with f injective — «this our only assumption on f, even when it is non-linear» (p. 3). What is fixed
across environments: «the constancy of the invariant parameters η, µc, σc², f across environments—the dependence of
µe, σe on the environment allows for varying distributions» (p. 3). So environments differ in **both the
environmental mean µe (a d_e-vector) and a scalar isotropic variance σe²**; there is no prior over environmental
parameters (p. 3). The expected regime is stated: «typically, for complex, high-dimensional data we would expect
E < dc ≪ de» (p. 3). The label causes the features (anti-causal); the log-odds of y are linear in zc with invariant
βc = 2µc/σc² and β₀ = log(η/(1−η)) (p. 4). The opposite-sign class means are «merely for clarity»; «it is
straightforward to extend our results to arbitrary means» (p. 3).

### 2. ⭐⭐ Everything is population-level: infinite samples per environment (p. 4)

«We consider the setting where we are given infinite samples from each environment; this allows us to isolate the
behavior of the objectives themselves, rather than finite-sample effects» (p. 4). The goal is worst-case risk on an
unseen environment, not average risk: «we are not hoping to minimize risk in expectation over the environments; this
is already accomplished via ERM or distributionally robust optimization» (p. 4). No theorem in the paper contains a
sample size or a noise-to-signal term in the environment count; the only finite-sample remark is in the non-linear
discussion (Kernthese 12).

### 3. The target is the "optimal invariant predictor", which is not the Bayes classifier (Definition 1, p. 4)

Φ*(x) = [I 0; 0 0] ∘ f⁻¹(x) = [zc] with β̂* = [2µc/σc², log(η/(1−η))] (Definition 1, p. 4). It «must recover all
invariant features», unlike Arjovsky et al.'s Definition 3 (p. 4). «the Bayes classifier uses environmental features
which are informative of the label but non-invariant; the optimal invariant predictor explicitly ignores these
features» (p. 4). So any predictor that beats it on the training environments must use environmental features (data
processing inequality, p. 6).

### 4. IRM is a bilevel objective; its penalised form and the alternatives (REx etc.) are treated together (pp. 5, 23–25)

IRM (eq. 4): minimise the average risk over environments s.t. β̂ ∈ argmin_β R^e(Φ, β) for every e; the penalised
version (eq. 5) adds λ‖∇_β̂ R^e(Φ, β̂)‖²₂, equivalent to eq. 4 at λ = ∞ (p. 5). The optimal invariant predictor is a
stationary point of eq. 4 (Proposition 4.1, p. 5; proof App. B, p. 13). All results are extended by corollaries to
REx-style risk-variance penalties and similar objectives (App. E, pp. 23–25).

### 5. ⭐⭐ Theorem 5.1 (linear case): IRM recovers the invariant predictor iff E > d_e (p. 6; restated and proved pp. 16–17)

Verbatim (p. 6):

> «**Theorem 5.1** (Linear case). Assume f is linear. Suppose we observe E training environments. Then the following
> hold:
> 1. Suppose E > de. Consider any linear featurizer Φ which is feasible under the IRM objective (4), with invariant
>    optimal classifier β̂ ≠ 0, and write Φ(f(zc, ze)) = Azc + Bze. Then under mild non-degeneracy conditions, it
>    holds that B = 0. Consequently, β̂ is the optimal classifier for all possible environments.
> 2. If E ≤ de and the environmental means µe are linearly independent, then there exists a linear Φ—where
>    Φ(f(zc, ze)) = Azc + Bze with rank(B) = de + 1 − E—which is feasible under the IRM objective. Further, both the
>    logistic and 0-1 risks of this Φ and its corresponding optimal β̂ are strictly lower than those of the optimal
>    invariant predictor.»

The count is **strict: E > d_e, i.e. E ≥ d_e + 1**; E = d_e is in the failure case (part 2). The informal version
(Theorem 3.1, p. 4) adds that for E > d_e the solution «generalizes to all environments with minimax-optimal risk».
Both parts assume **linear f and linear Φ**; part 1 assumes β̂ ≠ 0 and the non-degeneracy conditions (Kernthese 7);
part 2 assumes linearly independent environmental means. The theorem is about **feasible points / the global
minimum of the IRM objective in a binary classification problem**, not about a learning algorithm or a sample size.

### 6. ⭐⭐ Corollary 5.2: the threshold is an "if and only if" (p. 6)

«For both logistic and 0-1 loss, the optimal invariant predictor is the global minimum of the IRM objective if and
only if E > de» (Corollary 5.2, p. 6). With E ≤ d_e «the global minimum necessarily uses these non-invariant features
and therefore will not universally generalize to unseen environments» (p. 6). The authors call E > d_e «the (perhaps
unlikely) case» (p. 6) and frame the result as a lower bound: «observing a large number of environments—linear in the
number of environmental features—is necessary for generalization in the linear regime» (p. 6).

### 7. ⭐ The conditions: general position of means, unequal variance ratios, "measure zero", "with probability 1" (pp. 6, 14, 16)

The "mild non-degeneracy conditions" of part 1 are eqs. 9–10 (p. 14). (9): for any way of writing one environmental
mean as a linear combination of the others, µe = Σᵢ αᵢᵉ µᵢ, the coefficients satisfy Σᵢ αᵢᵉ ≠ 1 — «none of the
environmental means are are an affine combination of the others; in other words, they lie in general linear position,
which is the same requirement as Arjovsky et al. (2019)» (p. 14; the doubled "are" is in the original). (10): the
ratio (σa² − Σᵢ αᵢᵃ σᵢ²)/(1 − Σᵢ αᵢᵃ) differs for at least two environments a, b — «a similarly lax non-degeneracy
requirement regarding the relative scale of the variances» (p. 14). Genericity: «the set of environments which do not
satisfy Theorem 5.1 has measure zero under any absolutely continuous density over environmental parameters» (p. 6),
and «if E ≤ de then the environmental means will be linearly independent almost surely» (p. 14). In the proof of
Lemma C.3: «We will show that for any β̂, A, with probability 1 only B = 0 can satisfy Equation 12 for every
environment» (p. 16). Gaussianity enters through the closed-form optimal coefficient 2Σ̄⁺µ̄ (Lemma F.2, pp. 25–26).

### 8. ⭐⭐ Each environment removes one degree of freedom: d_e + 1 − E environmental directions stay feasible (Lemmas C.1–C.2, pp. 14–15; p. 7)

Lemma C.1 (p. 14): for E ≤ d_e linearly independent means there is a unique unit vector p with pᵀµe = σe² µ̃ for every
e (µ̃ maximal), constructed as p_α = U⁻ᵀσ / ‖U⁻ᵀσ‖ with µ̃ = 1/‖U⁻ᵀσ‖. On the feature z̃e = pᵀze the label-conditional
law is N(y·σe² µ̃, σe²), whose optimal coefficient 2µ̃ is the same in every environment (pp. 6–7, 15). Lemma C.2
(p. 14): «there exists a linear Φ with rank(Φ) = dc + de + 1 − E whose output depends on the environmental features,
yet the optimal classifier on top of Φ is invariant». The construction adds one direction per environment removed
(«Working backwards, we can repeat this process for each additional mean, such that rank(Φ) = dc + 1 + (de − E)»,
p. 15), and «for E = 1 any Φ will be vacuously feasible» (p. 15). The authors' own reading: «each environment reduces a
"degree of freedom" of the solution, such that only the invariant features remain feasible if enough environments are
seen» (p. 7). In arithmetic: E environments exclude E − 1 of the d_e environmental directions; all d_e are excluded
only at E = d_e + 1. *Typo in the source:* p is declared «p ∈ R^{dc}» on p. 15; it must be R^{d_e}.

### 9. ⭐ Theorem 5.3: below the threshold a purely environmental predictor can beat the invariant one, until a signal-dependent count (p. 7; proof pp. 18–19; Fig. C.2, p. 19)

«Suppose we observe E ≤ de environments, such that all environmental means are linearly independent. Then there
exists a feasible Φ, β̂ which uses only environmental features and achieves lower 0-1 risk than the optimal invariant
predictor on every environment e such that σe µ̃ > σc⁻¹‖µc‖₂ and 2σe µ̃ σc⁻¹‖µc‖₂ ≥ |β0|» (Theorem 5.3, p. 7). The
first condition compares the strength of the best invariant-classifier-compatible environmental direction (σe µ̃)
with the invariant signal-to-noise ‖µc‖/σc. It is «more likely to be met when E ≪ de» (p. 7). Simulations (App. C.4,
µe ∼ N(b, 4I), 20 runs, Fig. C.2, p. 19) show σe µ̃ tracking √(de − E) when σe² = 1, with «the crossover point …
approximately de − σe² dc» and «the number of environments observed before the condition ceases to hold is quite
high—on the order of de − dc» (pp. 7, 19). This crossover is a **signal-strength–dependent environment count, distinct
from the identifiability threshold d_e + 1**, and it is established only by simulation, not proved.

### 10. ⭐ The synthetic experiment reproduces the threshold — and shows IRM is hard to optimise (App. C.2, pp. 17–18; Fig. C.1)

d_c = 3, d_e = 6, η = 0.5, environmental means from a Gaussian prior; invariant features alone give 71.9 % accuracy,
environmental features 99.8 % in distribution (p. 17). Five runs. «optimizing for the IRM objective was quite
unstable, frequently collapsing to the ERM solution unless λ and the optimizer learning rate were carefully tuned»
(p. 17). Test environments are drawn from the same prior, and shift is simulated by negating the environmental means.
Result: «up until E = de , IRM essentially matches ERM in performance both in-distribution and under distribution
shift. As soon as we cross that threshold of observed environments, the predictor learned via IRM begins to perform
drastically better under distribution shift» (p. 17). Exception: «occasionally the invariant solution would be found
after only E = de = 6 environments» — conjectured to be an optimiser effect because the feasible non-invariant
solution is then a single point (p. 17). Fig. C.1 (p. 18) plots accuracy against E with ERM as dashed lines and the
optimal invariant predictor as a dashed green line. Data were generated «precisely according to our model» (p. 17),
and sample sizes per environment are not reported.

### 11. ⭐ Arjovsky et al. (2019) gave only a sufficient condition, needing more environments (p. 6)

«Roughly, their theorem says that for a learned Φ of rank r with invariant optimal coefficient β̂, if the training set
contains d − r + d/r "non-degenerate" environments, then β̂ will be optimal for all environments» (p. 6), with d the
observation dimension. Criticisms: no result ties the rank of Φ to the objective; for high-dimensional x it needs
Ω(d) environments; even for d = dc + de and r = dc it «would require E > de + de/dc , a stronger requirement than
ours»; and «they give no lower bound on the number of required environments» (p. 6). Rosenfeld et al. call their own
result a «substantially stronger (and simplified) upper bound in the classification setting, along with a matching
lower bound» (p. 6). Arjovsky's result is for regression (p. 6: «invariant optimal regression vector»).

### 12. ⭐⭐ Non-linear regime (Theorem 6.1; full Theorem D.3): no environment count rescues IRM; it generalises only near the training environments (pp. 7–9, 21–23)

«In the non-linear case, it's not clear how to capture this idea of restricting the "degrees of freedom"—and in fact
our results imply that this intuition is simply wrong. Instead, we show that the solution generalizes only to test
environments that are sufficiently similar to the training environments» (p. 7). Theorem 6.1 (simplified, σe² = 1,
p. 8): for any ε > 1 there is a featurizer Φ_ε which, with the ERM-optimal classifier, (1) has IRM penalty
O(p_ε²(c_ε de + ‖µ‖²)) with p_ε = exp{−de·min(ε − 1, (ε − 1)²)/8}, exponentially small in d_e; (2) equals the optimal
invariant predictor on ≥ 1 − p_ε of the training set and the ERM solution elsewhere; (3) on any test environment
whose mean satisfies min_y ‖µ_{E+1} − y·µe‖₂ ≥ (√ε + δ)√de for all training e (eq. 7), equals the ERM predictor on
≥ 1 − q of test points, q = (2E/(√π δ)) exp{−δ²}; (4) under Assumption 1 and reversed correlations
(µ_{E+1} = −Σ αe µe with eq. 8), has 0-1 risk > 0.975 − q. The theorem holds for **every E**; E enters only through
the union bound in q. The construction: Φ_ε = [zc] on balls of radius √(ε σe² de) around ±µe and [zc, ze] outside
(eqs. 18–19, p. 22). Assumption 1 (p. 8) requires the pooled ERM classifier on the environmental features to be
γ-close to optimal in every training environment (Definition 2, p. 7). Full version with arbitrary σe²: Theorem D.3,
p. 21.

### 13. "Covering" the environments would work — but then ERM would too (p. 9; p. 23)

«We can of course guarantee generalization if the training distributions "cover" (or approximately cover) the full
space of environments in order to tie down the performance on future distributions. But in such a scenario, there
would no longer be a need for ICP; we could expect ERM or DRO to perform just as well» (p. 9). Under a Gaussian prior
on the environmental means, training means are O(√de) apart; a test mean Ω(√(de log E)) away already makes q ∈ O(1/E)
(p. 9; App. D.2, p. 23).

### 14. Finite samples make the problem worse, not better (p. 9)

On part 2 of Theorem 6.1: «for large enough de, the "fake" predictor will often be a preferred solution. In the
finite-sample setting, we would need exponentially many samples to even distinguish between the two!» (p. 9). This
is the paper's only sample-size statement, and it is qualitative.

### 15. Alternatives to IRM fail in the same settings (App. E, pp. 23–25)

Corollary E.1 (p. 24): with E ≤ de, linearly independent means **and identical variance σe²**, the Lemma C.2
featurizer beats the optimal invariant predictor under the REx risk-variance penalty «for any choice of
regularization parameter λ ∈ [0, ∞]». Corollary E.2 (p. 24) extends Theorem 5.3 to REx; Fig. E.1 (p. 25) repeats
Fig. C.2. The non-linear construction carries over unchanged (App. E.2, p. 24).

### 16. ⭐ Bottom line: IRM and its alternatives do not improve over ERM (pp. 1, 3, 9)

«Thus, in this setting we find that IRM and its alternatives fundamentally do not improve over standard Empirical
Risk Minimization» (abstract, p. 1). The authors' concurrent paper (Rosenfeld, Ravikumar & Risteski 2021, "An online learning approach to interpolation
and extrapolation in domain generalization", arXiv:2102.13128; not in the project bib) is cited for «ERM remains optimal for both interpolation and extrapolation in the linear setting and
that the latter is exponentially harder than the former» (p. 3).

## Figures and tables

- Fig. 3.1 (p. 4): Bayesian network of the model (e → ze; y → ze, zc; ze, zc → x).
- **Fig. C.1 (p. 18): IRM vs ERM accuracy against the number of training environments E (d_c = 3, d_e = 6), with and
  without reversed correlation; the switch occurs after E = d_e.**
- Fig. C.2 (p. 19): σe µ̃ against environments observed for (d_c, d_e, σe²) ∈ {(40, 200, 2), (100, 200, 1),
  (50, 300, 2), (60, 300, 1)}, with √(de − E) and √dc reference lines.
- Fig. E.1 (p. 25): the same for the REx variant (σe² ∈ {1, 3}).
- No tables.

## What this source does not address

- **No interaction term and no context-conditioned predictor.** The environment's identity and parameters are never
  an input to Φ or β̂; the goal is to *ignore* environment-specific structure, not to model it (pp. 3–4).
- **No rank ≤ E − 1 statement about a learned map or a response matrix.** The nearest statement is the residual
  feasible rank d_e + 1 − E of the environmental block of a featurizer (Theorem 5.1(2), Lemma C.2, pp. 6, 14–15).
- **No regression.** All results are for binary classification with logistic and 0-1 loss (p. 4); Arjovsky et al.'s
  regression result is discussed but not re-proved (p. 6).
- **No finite-sample or noise-dependent environment count.** Population-level throughout (p. 4); the only noise-like
  dependence is the signal-strength crossover of Theorem 5.3, found by simulation (p. 19).
- **No effective dimension.** d_e is the exact dimension of an isotropic latent block whose environmental means must
  span R^{d_e} (pp. 14, 16); spectra, participation ratios or low-rank µe are not considered.
- **No real data.** One synthetic experiment with d_c = 3, d_e = 6 (p. 17).

<!-- ═══ FOKUS ═══ -->

## Focus: Limits / ceilings / context count / what is already shown — and exactly what overlaps with the pre-registered hypotheses H1, H2a, H2b, H2c, H3a, H3b, H4 — 2026-11-05

**Carries Kernthesen 1, 2, 5, 6, 7, 8, 12** (core), and 9, 10, 11, 13, 14 (context). The hypotheses are those of the
plan at `input/ideas/zero-shot-limits-context-complexity-plan.md` (H1–H4 at the top; Amendments A1–A3 at the end).
Amendment A3 states that «Proposition 1 is cited from domain-generalisation theory (rosenfeld-2021-risks-irm)»; this
block checks whether that citation is safe.

### What this focus needs beyond the Kernthesen

1. **The exact threshold.** Linear f, linear Φ, binary labels, infinite samples per environment: the optimal
   invariant predictor is the IRM global minimum **iff E > d_e** (Corollary 5.2, p. 6), i.e. iff E ≥ d_e + 1. Part 1
   (success) needs β̂ ≠ 0 and the non-degeneracy conditions 9–10, which fail only on a measure-zero set (pp. 6, 14).
   Part 2 (failure) needs only linearly independent environmental means, which hold almost surely for E ≤ d_e
   (p. 14). With E ≤ d_e the surviving non-invariant block has rank d_e + 1 − E (Theorem 5.1(2), p. 6).
2. **Mapping to the project's notation** (δ_{c,t} = μ + a_t + b_c + i_{c,t} + ε_{c,t}; context descriptor x_c of
   effective dimension d; the plan's Notation block). Reader's mapping, not in the source:

   | Rosenfeld et al. 2021 | Project | Fit | Anchor |
   |---|---|---|---|
   | E training environments | C training contexts (cell lines) | good: both are discrete, labelled groups | p. 3 |
   | environment parameters (µe ∈ R^{d_e}, σe²) | context descriptor x_c ∈ R^d, and the context terms b_c, i_{c,t} | partial: µe is a latent, never an input; x_c is an observed model input | pp. 3–4 |
   | d_e (exact dimension of ze; means must span R^{d_e}) | d (participation ratio of the descriptor matrix, ≈ 20.6 on Tahoe) | weak: exact rank vs soft effective dimension | pp. 14, 16 |
   | invariant features zc, joint law with y fixed (η, µc, σc², f) | conserved component a_t (with μ) | good in spirit | p. 3 |
   | environmental features ze | b_c + i_{c,t} | partial | p. 3 |
   | binary label y, logistic / 0-1 loss | continuous response over genes, residual MSE | none: classification only | p. 4 |
   | IRM: learn Φ that **ignores** ze | context model: **use** x_c to predict i_{c,t} | inverted goal | pp. 4–5 |
   | ERM (pooled over environments) | context-free model (μ + â_t, pooled) | good | pp. 1, 8–9 |
   | unseen test environment | held-out cell line (zero-shot) | good | p. 4 |
   | infinite samples per environment | finite cells; σ²_ε, Var(i), R² | none | p. 4 |
   | "degrees of freedom" removed: E − 1 of d_e | rank ≤ C − 1 of the identifiable interaction (Prop. 1) | same arithmetic, different object | pp. 7, 15 |

3. **What the analogy supports.** (i) A published, population-level, "dimension + 1" floor on the number of
   environments for linear identification (Theorem 5.1, Corollary 5.2, p. 6), with a matching lower bound — the
   qualitative claim "you need more contexts than descriptor dimensions". (ii) The same counting as Proposition 1:
   E environments eliminate E − 1 directions, leaving d_e + 1 − E (Lemma C.2, pp. 14–15; «each environment reduces a
   "degree of freedom"», p. 7). (iii) Below the floor, the richer objective is no better than the pooled one (ERM)
   (p. 1; Fig. C.1, p. 18) — the theoretical cousin of a context-free baseline winning at small C. (iv) A synthetic
   learning curve over E with a pooled comparator and a switch at E = d_e + 1 (Fig. C.1).
4. **What it does not support.** (i) **No rank ≤ C − 1 statement for an interaction term, a response matrix or a
   context-conditioned predictor exists in the paper** (checked pp. 1–27). (ii) No regression, no continuous
   response (p. 4). (iii) No finite-sample or noise term, so nothing for the noise part of the rejected
   C* = d + 1 + d·σ²_res/(R²·Var(i)) (p. 4); the paper speaks only to the **d + 1 term**. (iv) No effective
   dimension: d_e is a nominal rank (p. 14). (v) The project's context model (FiLM r32, plan definition 1) is
   non-linear; for non-linear models the paper says the degrees-of-freedom intuition «is simply wrong» (p. 7) and
   proves a failure for every E (Theorem 6.1, p. 8). The analogy therefore covers a linear-in-descriptor context
   model only.
5. **Verdict on Amendment A3.** Citing this paper *for* Proposition 1 is **safe only as an analogy** ("cf."): it is
   the published precedent for an environment-count threshold E > d_e and for the E − 1 degrees-of-freedom count in
   linear invariant prediction (Theorem 5.1, Corollary 5.2, Lemma C.2). It is **not safe** to write that Proposition 1
   "is proven in" or "follows from" it: the paper has no interaction, no regression, no context-conditioned
   predictor and no rank statement about a learned map. Suggested wording (reader's proposal): *"Proposition 1 is the
   cell-line analogue of the environment-count threshold in linear invariant prediction (Rosenfeld et al. 2021,
   Thm. 5.1 and Cor. 5.2; cf. Arjovsky et al. 2019); its proof is elementary — the centred C × d descriptor matrix
   has rank ≤ min(C − 1, d)."* Closer formal matches named in `input/bibliography/literaturguide.md`
   (@wang-2022-subspace-recovery, ISR-Mean, E ≥ d_s + 1; @tripuraneni-2021-linear-representations;
   @du-2021-few-shot-representation; @rothenhausler-2021-anchor-regression) are not yet ingested and should be
   checked before the Theory section is drafted.
6. **Reader's derivation: the threshold is an affine-span condition.** Condition 9 is about *affine* combinations
   (Σαᵢ ≠ 1, p. 14). In eq. 12 (p. 16), if all environments share one variance σ², the term σ²BBᵀβ̂ is constant, so
   Bµe must be the same for every e; hence B(µe − µ₁) = 0 for all e, and B = 0 as soon as the E − 1 differences span
   R^{d_e}, i.e. E ≥ d_e + 1. This is the same "C − 1 centred directions" arithmetic as Proposition 1. It also shows
   that the equal-variance case, which condition 10 *as written* excludes (the ratio is then σ² for every
   environment), still gives B = 0. This derivation is not in the paper; the paper's own proof route uses
   condition 10.
7. **Numbers against the project (from `knowledge/_meta/log.md`).** With d ≈ 20.6 the floor is d + 1 ≈ 21.6. The
   empirical break-even C*_emp = 37.4 [33.3, > 40] (H2a) lies above it, so the measured curve does not contradict an
   E > d-type floor. The floor is necessary, not sufficient, and its transfer assumes d_e ↔ participation-ratio d,
   which the paper does not license. The signal-dependent crossover of Theorem 5.3 (≈ d_e − σe² d_c, Fig. C.2, p. 19)
   is the paper's only count that depends on signal strength; it is simulated, concerns a different comparison
   (environmental-only vs invariant predictor on the training environments), and cannot be used to rescue the
   rejected Proposition 2.
8. **Cited related counts.** Arjovsky et al. 2019: d − r + d/r non-degenerate environments for a Φ of rank r
   (observation dimension d), i.e. E > d_e + d_e/d_c in the ideal case, sufficient only (p. 6;
   @arjovsky-2019-irm, in the project bib, not yet ingested). Non-linear: any E fails for
   test means (√ε + δ)√d_e away (eq. 7, p. 8).

### Overlap with H1–H4

| Hyp. | Verdict | What the source shows exactly | Anchors |
|---|---|---|---|
| **H1** (frozen pipeline ≥ 0.7 of r_max; < 50 % of remaining error in a_t + b_c) | **none** | No ceiling, no variance components, no noise. The "optimal invariant predictor" (uses only the invariant block) is only loosely like the r_max oracle that knows a_t and b_c. | Def. 1, p. 4 |
| **H2a** (Tahoe: learned context correction negative at C = 12; C*_emp ≥ 24) | **partial (theoretical analogue only)** | Below E = d_e + 1 the richer objective does no better than pooled ERM, and a synthetic curve over E switches only after E = d_e (Fig. C.1). It is about IRM vs ERM on binary Gaussian data, not a context correction vs a context-free regression, so it predicts no C*_emp for Tahoe. | Thm. 5.1, Cor. 5.2, p. 6; pp. 17–18 |
| **H2b** (theory predicts the break-even within ×2) | **partial (floor term only)** | Supports only the noise-free d + 1 term of C*, at population level, for linear models. Nothing on σ²_res, R² or Var(i); infinite samples by assumption. Consistent with the H2b rejection: the paper offers no noise-dependent count to fall back on (the Thm. 5.3 crossover is simulated and answers a different question). | p. 4; p. 6; p. 19 |
| **H2c** (genetic C*_min > C_public) | **partial (floor term only, inferred)** | If the analogy is accepted, C_public = 11 (Count 1) or 25 (Count 2) must at least exceed d for any linear context model to be identifiable; the noise term of C*_min has no counterpart. The genetic d is not yet estimated, so no verdict follows from this paper. | Cor. 5.2, p. 6 |
| **H3a** (destination control-PC residual energy) | **none** | No residual or subspace analysis of prediction errors. | — |
| **H3b** (M3 vs five corrections) | **none** | No correction methods compared. | — |
| **H4** (same-line data vs modelling gain) | **none (thematic echo only)** | «IRM and its alternatives fundamentally do not improve over standard Empirical Risk Minimization» (p. 1) echoes "data beat modelling", but no same-environment data value is analysed. | p. 1 |
| *Prop. 1 (A3 citation)* | **analogy, not source** | E > d_e threshold and d_e + 1 − E residual rank are in the paper; rank ≤ C − 1 of an interaction is not (focus item 5). | pp. 6, 14–15 |

### Direct quotes (for the draft)

> "If E > de , then Φ, β̂ is precisely the optimal invariant predictor; it uses only invariant features and
> generalizes to all environments with minimax-optimal risk. If E ≤ de , then Φ, β̂ relies upon non-invariant
> features." (Theorem 3.1, informal, p. 4)

> "For both logistic and 0-1 loss, the optimal invariant predictor is the global minimum of the IRM objective if and
> only if E > de ." (Corollary 5.2, p. 6)

> "we demonstrate that observing a large number of environments—linear in the number of environmental features—is
> necessary for generalization in the linear regime." (p. 6)

> "We consider the setting where we are given infinite samples from each environment; this allows us to isolate the
> behavior of the objectives themselves, rather than finite-sample effects." (p. 4)

> "Our results—and those of Arjovsky et al. (2019)—intuitively proceed by observing that each environment reduces a
> "degree of freedom" of the solution, such that only the invariant features remain feasible if enough environments
> are seen. In the non-linear case, it's not clear how to capture this idea of restricting the "degrees of
> freedom"—and in fact our results imply that this intuition is simply wrong." (p. 7)

> "Similar to Arjovsky et al. (2019), the set of environments which do not satisfy Theorem 5.1 has measure zero under
> any absolutely continuous density over environmental parameters." (p. 6)

### Examples & illustrations (for later drafting)

- Fig. C.1 (p. 18): d_c = 3, d_e = 6. IRM tracks ERM for E ≤ 6 and approaches the optimal invariant predictor for
  E ≥ 7 under reversed correlations — a ready-made picture of a "d + 1" switch with a pooled comparator.
- Arjovsky et al.'s condition with d = dc + de, r = dc gives E > de + de/dc (p. 6): for d_c = 3, d_e = 6 that is
  E > 8, against E > 6 here (reader's arithmetic).
- Fig. C.2 (p. 19): for d_e = 200, d_c = 100, σe² = 1 the environmental-only predictor keeps beating the invariant one
  until roughly d_e − σe² d_c ≈ 100 environments (reader's reading of the stated rule). The count depends on the
  relative strength of the two signals, yet stays inside the E ≤ d_e failure region of Theorem 5.1.

### Boundary of this focus

This focus does not use the IRM-vs-alternatives material (App. E) or the proof machinery (App. F) beyond the
conditions they impose. A reader looking for a noise-dependent environment count, a regression result, or a rank
bound on a learned context-conditioned map will not find one here.

## Mentioned entities

- Concepts: [[context-count-learning-curve]]
- Related project sources: [[ji-2024-universal-prediction]] (empirical count curve with a context-free baseline;
  this paper supplies the theoretical floor that A3 pairs with it).

## Connections

- Supports: [[context-count-learning-curve]] (theoretical floor E > d_e and synthetic curve, Fig. C.1, p. 18).
- Cites (not yet ingested): @arjovsky-2019-irm (the sufficient d − r + d/r condition, p. 6).
- Named as closer formal matches for Proposition 1 in the literaturguide (not cited by this paper, not yet ingested):
  @wang-2022-subspace-recovery, @tripuraneni-2021-linear-representations, @du-2021-few-shot-representation,
  @rothenhausler-2021-anchor-regression, @chen-2022-iterative-feature-matching.
- Earlier reading: `input/bibliography/literaturguide.md` pre-read note («Rosenfeld et al. 2021 (E > d_s)»; «linear-
  Gaussian IRM needs E > d_s environments») — verified: correct as to the count, but the guide's d_s is the paper's
  d_e (environmental, "spurious" dimension; the paper avoids the word "spurious", p. 3, fn. 1), and the guide
  attributes the "same C−1 span argument as Proposition 1" to Wang et al. 2022 (unverified here), not to this paper —
  consistent with focus item 4(i).
