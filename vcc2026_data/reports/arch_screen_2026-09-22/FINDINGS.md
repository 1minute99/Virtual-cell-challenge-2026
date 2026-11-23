# Architecture screen: B(g) + context residual versus fully conditional models — 2026-09-22

One fixed protocol for 17 architectures (`scripts/arch_zoo.py`, driver `arch_screen.py`): seed 0,
extended gene contexts (35 training; JURKAT held for early stopping), identical inputs (context =
training-only PCA of the pooled control profile, 36-d, or the same profile cut into 32 gene-block
tokens; perturbation = rank-32 PCA of the refit programme, or the 64-d frozen ESM features for the
joint models that must not see B(g)), identical bases and evaluation. Early stopping: ≥ 2,000
steps, patience 30 evals (3,000 steps), max 20,000. Held-out K562 (no measurable context
specificity) and, as a replicate, HepG2 (real specificity). Official `cell-eval2 vcc2026` preset
for correct / wrong / zero context at alpha 1 (104 runs), harness PDS also at alpha 0.25.
Base = extended refit (refit_v2): K562 0.3481, HepG2 0.4270 official composite.

## Ranking (official composite Δ vs the base, alpha 1, correct context)

| architecture | family | K562 | HepG2 | mean | correct − wrong K / H | harness PDS Δ at α 0.25 K / H | params |
|---|---|---|---|---|---|---|---|
| additive shift (context only) | residual | +0.0154 | −0.0083 | **+0.0036** | +0.026 / +0.001 | +0.0048 / −0.0303 | 2,304 |
| mixture of experts | residual | −0.0058 | +0.0127 | +0.0034 | +0.043 / +0.046 | +0.0029 / −0.0122 | 24,552 |
| bilinear | residual | −0.0132 | +0.0152 | +0.0010 | +0.017 / +0.058 | +0.0007 / −0.0034 | 35,873 |
| full FiLM (γ and β, 2 blocks) | residual | +0.0008 | +0.0002 | +0.0005 | +0.027 / +0.055 | +0.0021 / −0.0202 | 31,648 |
| full hypernetwork | residual | +0.0003 | −0.0025 | −0.0011 | +0.009 / +0.006 | +0.0045 / −0.0221 | 141,888 |
| LoRA hypernetwork | residual | −0.0211 | +0.0169 | −0.0021 | +0.010 / +0.024 | −0.0011 / −0.0120 | 18,080 |
| token Transformer (STATE proxy) | residual | −0.0156 | +0.0088 | −0.0034 | +0.033 / +0.058 | +0.0004 / −0.0145 | 86,496 |
| concat MLP F([h_g, c]) | residual | +0.0029 | −0.0100 | −0.0035 | +0.030 / +0.030 | +0.0012 / −0.0371 | 10,912 |
| FiLM (S4's model) | residual | −0.0118 | −0.0061 | −0.0089 | +0.024 / +0.024 | +0.0035 / −0.0181 | 4,416 |
| gated / multiplicative | residual | −0.0184 | −0.0003 | −0.0094 | +0.026 / +0.033 | +0.0036 / −0.0191 | 15,136 |
| cross-attention | residual | −0.0584 | +0.0075 | −0.0254 | −0.013 / +0.035 | −0.0028 / −0.0130 | 32,800 |
| biolord-style | joint | −0.0438 | −0.0639 | −0.0539 | +0.019 / +0.036 | – | 25,728 |
| CPA-style | joint | −0.0454 | −0.1339 | −0.0897 | +0.014 / −0.036 | – | 29,888 |
| joint MLP + programme | joint | −0.1464 | −0.0902 | −0.1183 | −0.013 / +0.005 | – | 34,112 |
| conditional flow | joint | −0.1624 | −0.1264 | −0.1444 | +0.030 / +0.032 | – | 50,688 |
| flow matching | joint | −0.2165 | −0.1015 | −0.1590 | −0.130 / −0.016 | – | 42,432 |
| joint MLP (ESM + context) | joint | −0.1933 | −0.2117 | −0.2025 | +0.022 / +0.061 | – | 109,120 |

Per-line official tables: `official_K562_seed0.json`, `../arch_screen_2026-09-22_HepG2/`; harness
and cheap diagnostics: `results.json` in each; `final_ranking.json`.

## Answers

1. **Residual on B(g) or joint?** Residual, decisively. Every joint model loses 0.04–0.22 on both
   lines, almost entirely through PDS (−0.24 to −0.90 normalised): no conditional model recovered
   the per-target mean that the refit supplies for free. CPA/biolord-style latent composition is
   the least bad joint family and still −0.05 to −0.13.
2. **Among residual models the differences are within a few thousandths** and line-dependent:
   K562 favours the perturbation-blind additive shift (+0.015, a calibration effect, no context
   specificity), HepG2 favours bilinear / LoRA / MoE (+0.013–0.017). Only full FiLM is non-negative
   on both lines, at +0.0005. At alpha 1 every residual model destroys HepG2's harness PDS
   (−0.05 to −0.53 normalised); at alpha 0.25 the damage is −0.003 (bilinear) to −0.037.
3. **Correct beats wrong context** for 10/11 residual models on K562 and 11/11 on HepG2 on the
   official composite, by +0.01 to +0.06, even where the residual-cosine diagnostic is negative
   (K562). Cross-attention is the exception on K562.
4. **Capacity and attention buy nothing.** The 141k-parameter hypernetwork, the 86k token
   Transformer and cross-attention all rank below 2k–35k-parameter models; cross-attention is the
   worst residual model on both lines.
5. **Generous early stopping changed nothing material**: best steps ranged 100–18,100 with no
   relation to rank.

## Candidates for multi-seed / multi-line validation

1. **Bilinear** (c^T A_k p, rank-16 factors, 36k params): the only model that is positive on the
   informative line (+0.015 HepG2) with the strongest correct−wrong margin there (+0.058), the
   smallest alpha-0.25 harness loss (−0.003), and the largest residual gain at alpha 0.25 among
   those; K562 −0.013 at alpha 1 needs the nested alpha.
2. **Mixture of experts** (context router over 8 linear experts): best mean of the
   perturbation-aware models (+0.0034) and the largest correct−wrong margin on both lines
   (+0.043 / +0.046).
Full FiLM is the safe third (never negative). The additive shift's K562 lead is a calibration
artefact worth understanding but not a context model. S4's FiLM ranks 9th here at alpha 1; its
leaderboard win came at alpha 0.25 on the nested refit, which this screen's alpha-1 protocol
penalises.
