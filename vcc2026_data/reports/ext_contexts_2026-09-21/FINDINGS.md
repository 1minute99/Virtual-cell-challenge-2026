# Do more real gene-perturbation contexts (Feng HipSci lines + primary CD4 donors) help VCC 2026?

**Dates:** data 2026-09-21, models and official scoring 2026-09-22. **Status:** COMPLETE. 222 official
preset runs (74 arms x 3 seeds over 5 lines), all PDS-parity checked.

Everything here reuses the current pipeline: the frozen training table and its 5,592-gene axis,
`LeanSourceBaseline` rank 32 for the refit, `simple_residual_baselines.empirical` for mean
transfer, the FiLM r32 residual of the drug-context experiment, the five-line harness for
development PDS and `cell-eval2 --preset vcc2026` for every final number. Strict
leave-one-scored-biology-out throughout: the held-out line contributes only its control profile;
on the KOLF fold the two KOLF-background HipSci lines are held out as well.

## 1. What was added (Part 1–2)

| dataset | what exists | usable per-context data | decision |
|---|---|---|---|
| Feng broad screens (34 lines, 4,959 + 2,259 targets) | pooled counts + per-cell metadata on disk | **1–2 cells per (line, target)**; 0 and 8 pairs at the table's 25-cell floor; 12/36 confirmed NTC cells in total | kept **pooled** as biology FENG (already in training); the fitness arm's 1,726 pooled pairs were added to FENG |
| Feng targeted screen | 1.16 M cells, 19 lines, 445 targets | 23 cells median per pair, ~490 NTC cells per line; 8,226 pairs, **7,435 at ≥25 cells** | **19 line contexts** (17 with ≥200 targets) |
| CD4 Perturb-seq (Marson) | author pseudobulk, 42 GB | re-aggregated to donor x state over ALL targets: 12 contexts, ~45,000 pairs per condition, 68–88k NTC cells per context | **12 donor x state contexts**, 123,543 rows at ≥25 cells, replacing the 35,832 pooled CD4 rows |

Coverage (`coverage.json`): FiLM contexts 9 old + 12 CD4 + 17 Feng = **38**; residual panel
13,070 targets (10,580 measured in ≥12 contexts). Genes on the official axis measured in
≥5 / ≥10 / ≥20 contexts: **18,527 / 18,527 / 18,513** (before: 18,011 / – / – with 9 contexts).
Extended table: 213,085 rows (191k genetic) versus 116,217.

Two labelings, because the refit averages per BIOLOGY: **v1** every new context is its own
biology (contexts vote), **v2** the new contexts nest inside CD4 / FENG (finer data, same votes).

## 2. Baselines (Part 3) — official vcc2026 composite, refit_cur = 0.2755

| model | K562 | RPE1 | HepG2 | Jurkat | KOLF | 5-line Δ | paired seeds |
|---|---|---|---|---|---|---|---|
| refit, current | 0.3351 | 0.0553 | 0.4184 | 0.3638 | 0.2048 | 0 | – |
| **refit + Feng + CD4, nested (v2)** | +0.0130 | −0.0098 | +0.0080 | +0.0133 | +0.0056 | **+0.0060** | **12/15**, 4/5 lines |
| refit + Feng + CD4, contexts as biologies (v1) | −0.0299 | −0.1766 | −0.0442 | −0.0348 | +0.0008 | −0.0569 | 2/15 |
| mean transfer, current | −0.0996 | −0.0115 | −0.0914 | +0.0009 | −0.0093 | −0.0422 | 2/15 |
| mean transfer + Feng + CD4, v2 | −0.0777 | −0.0330 | −0.0837 | +0.0099 | −0.0052 | −0.0380 | 3/15 |
| mean transfer + Feng + CD4, v1 (36 voters) | −0.1330 | −0.2350 | −0.1557 | −0.0806 | −0.0049 | −0.1218 | 0/15 |

Development PDS agreed in sign on every row. Reading: the new rows help the refit **only when
they do not change the cross-biology weighting** (+0.006, small but 12/15 seeds and four lines);
letting 31 CD4/iPSC contexts vote as biologies costs −0.06 on the refit and −0.12 on mean
transfer, with KOLF (itself iPSC) the one line that does not lose. Mean transfer gains +0.004
from finer CD4/Feng data and is far below the refit either way.

## 3. FiLM on gene contexts (Part 4)

Same FiLM r32 recipe as the drug-context `scratch` arm, contexts from the extended table, context
profiles = the table's own pooled control profiles (the earlier profiles were scrambled, see
`CONTROLS_COMMON_BUG.md`), context PCA on training gene contexts only, inner early-stopping
biology as before, alpha by nested selection on the other four lines.

**With correct profiles the model finally trains.** Best inner-validation step 500–1,800 (was 100
in 20/20 scrambled runs) and it beats the zero residual on the inner biology in 4/5 folds (was 0/5).

Residual cosine on the held-out line (`film_diagnostics.json`), 35 training contexts:

| line | correct | wrong ctx | zero ctx | correct − wrong |
|---|---|---|---|---|
| K562 | +0.036 | +0.190 | +0.283 | −0.154 |
| RPE1 | +0.274 | +0.059 | −0.258 | +0.215 |
| HepG2 | +0.195 | +0.046 | −0.071 | +0.149 |
| Jurkat | +0.242 | +0.161 | +0.240 | +0.080 |
| KOLF | +0.174 | +0.184 | +0.198 | −0.010 |

### 3a. Context-count scaling (6 → 12 → 24 → 35, fixed order: old biologies, CD4, HipSci)

| line | correct − wrong at 6 / 12 / 24 / 35 | harness PDS (α 0.25) at 6 / 12 / 24 / 35 |
|---|---|---|
| K562 | +0.04 / +0.02 / −0.09 / −0.15 | 0.845 / 0.842 / 0.854 / 0.851 |
| RPE1 | +0.24 / +0.16 / +0.30 / +0.22 | 0.796 / 0.797 / 0.789 / 0.783 |
| HepG2 | +0.05 / +0.14 / +0.07 / +0.15 | 0.821 / 0.823 / 0.831 / 0.829 |
| Jurkat | +0.18 / −0.02 / +0.04 / +0.08 | 0.827 / 0.827 / 0.828 / 0.833 |
| KOLF | −0.11 / −0.01 / 0.00 / −0.01 | 0.573 / 0.569 / 0.576 / 0.573 |

No monotone trend in either quantity; PDS moves by < 0.01 between 6 and 35 contexts on every
line. (The base for these is refit_v1, which is 0.04 below the current refit on four lines.)

### 3b. Official composite of the FiLM arms

Base refit_v1 (`film_ext`), alpha 1 so correct/wrong/zero are matched:

| line | correct | wrong | zero | refit_v1 | refit_cur |
|---|---|---|---|---|---|
| K562 | 0.3318 | 0.2751 | 0.2589 | 0.3052 | 0.3351 |
| RPE1 | 0.0695 | −0.1036 | −0.2120 | −0.1214 | 0.0553 |
| HepG2 | 0.3869 | 0.3669 | 0.3114 | 0.3742 | 0.4184 |
| Jurkat | 0.3614 | 0.3307 | 0.2865 | 0.3289 | 0.3638 |
| KOLF | 0.1997 | 0.2013 | 0.2068 | 0.2057 | 0.2048 |

Correct > wrong > zero on 4/5 lines (KOLF flat); correct − wrong averages **+0.056** on the
official aggregate, and correct beats its own base refit_v1 on 4/5 lines by +0.01 to +0.19. It
does **not** beat the current refit: −0.0056 at alpha 1 (3/15 seeds), −0.038 at the nested alphas.
The comparator `film_old` (6 old contexts, scrambled profiles, base refit_cur) is +0.018 locally
and lost officially on the leaderboard (S3).

### 3c. film_v2 — the same trained residuals on the stronger refit_v2 base (the best model)

Nested alpha 0.25 on K562/RPE1/Jurkat/KOLF, 0.5 on HepG2. Official composite vs refit_cur:

| line | refit_cur | refit_v2 | **film_v2 (nested α)** | Δ vs refit_cur | seeds | film_v2 α1 | wrong α1 | zero α1 |
|---|---|---|---|---|---|---|---|---|
| K562 | 0.3351 | 0.3481 | **0.3565** | **+0.0214** | 3/3 | 0.3416 | 0.3077 | 0.2997 |
| RPE1 | 0.0553 | 0.0455 | **0.0914** | **+0.0361** | 3/3 | 0.1846 | 0.0069 | −0.0952 |
| HepG2 | 0.4184 | 0.4265 | **0.4570** | **+0.0386** | 3/3 | 0.4213 | 0.4001 | 0.3431 |
| Jurkat | 0.3638 | 0.3771 | **0.3862** | **+0.0224** | 3/3 | 0.3935 | 0.3657 | 0.3299 |
| KOLF | 0.2048 | 0.2104 | **0.2096** | **+0.0048** | 3/3 | 0.2028 | 0.2026 | 0.2169 |
| five-line | 0.2755 | 0.2815 | **0.3001** | **+0.0247** | **15/15** | 0.3088 (+0.0333, 12/15) | 0.2566 (−0.0189) | 0.2189 (−0.0566) |

Against its own base refit_v2 the residual adds +0.0186 (12/15 seeds). Correct − wrong at alpha 1:
K562 +0.034, RPE1 +0.178, HepG2 +0.021, Jurkat +0.028, KOLF +0.000 (mean **+0.052**); correct −
zero +0.042 / +0.280 / +0.078 / +0.064 / −0.014. Per metric (normalised, five-line mean vs
refit_cur): NMAE +0.097, reach +0.045, fidelity +0.014, Jaccard +0.006, **PDS −0.014**.

Unlike every previous local gain in this project, this one is **not carried by one line**: four
lines gain +0.021 to +0.039 and KOLF is flat-positive (RPE1 supplies 29 % of the mean, not 80 %).
It is still the local composite with locally built anchors, the quantity that predicted +0.021 for
the PCA interaction (leaderboard −0.002) and +0.018 for the old scratch FiLM (−0.0005); PDS, the
only target-discriminating metric, is slightly negative again. The difference from those two is
consistency (15/15, 5/5) and the correct-vs-wrong margin (+0.052 vs +0.026 before).

## 4. Answers (Part 6) — items 7 and 8 finalised when film_v2 completes

1. **Contexts gained:** 12 CD4 donor x state (all targets) + 19 HipSci lines (445 targets, 17 with
   ≥200 usable) = 31; 38 FiLM contexts in total. The 34-line broad screen yields none (1–2 cells
   per pair) and stays pooled.
2. **Mean transfer:** no. Finer data +0.004, contexts as voters −0.12.
3. **Refit:** yes, slightly, only when the new rows stay inside their biologies: **+0.0060**,
   12/15 seeds, 4/5 lines (RPE1 −0.010). As separate biologies −0.057.
4. **Multi-context FiLM vs the same FiLM without Feng/CD4:** on the matched base the 35-context
   model beats the 6-context model on harness PDS by < 0.01 on three lines and loses on two (§3a);
   officially, film_v2 (35 contexts, corrected profiles, extended refit base) is +0.0247 where the
   old 6-context film_old was +0.0180 locally and −0.0005 on the leaderboard. The improvement over
   film_old comes with three confounded changes at once (contexts, corrected control profiles,
   base), and the scaling curve attributes none of it to the context count.
5. **Performance vs context count:** flat. Discrimination (correct − wrong) does not grow from 6
   to 35 contexts; it falls on K562.
6. **Correct vs wrong context:** yes, clearly, on 4/5 lines both in residual cosine and on the
   official aggregate (+0.056 mean). K562 and KOLF show no context specificity.
7. **Best model vs refit on the official composite:** yes. **film_v2 (FiLM residual on the
   extended refit, nested alpha) +0.0247, 15/15 paired seeds, 5/5 lines**; refit_v2 alone +0.0060
   (12/15). Gains are NMAE, reach and fidelity; PDS −0.014. Spread across four lines, not RPE1-only.
8. **Add more gene datasets (Orion etc.)?** Not for context count: 31 real contexts moved
   neither discrimination nor generalisation, so the scaling answer is flat. What did move the
   composite was (a) feeding the refit finer data without changing its biology weighting (+0.006)
   and (b) a context residual whose inputs are finally correct (+0.025, 15/15). The next step
   should therefore be to test film_v2 where it counts — one leaderboard slot — rather than to
   ingest another dataset; if it holds officially, extend the refit with more biologies in the
   nested labeling and retrain the same residual; if it does not, the composite gain is the third
   local-only artefact and the residual direction closes for good.

## Final table (official vcc2026 composite, five-line means; `final_table.json`)

| model | contexts | composite | Δ vs refit | seeds | lines | PDS | NMAE | fid | reach | jac |
|---|---|---|---|---|---|---|---|---|---|---|
| refit, current | 26 bios | 0.2755 | 0 | – | – | 0 | 0 | 0 | 0 | 0 |
| mean transfer, current | 6 stores | 0.2333 | −0.0422 | 2/15 | 1/5 | −0.176 | +0.003 | −0.012 | −0.059 | −0.009 |
| mean transfer + Feng + CD4 (v2) | 8 voters | 0.2375 | −0.0380 | 3/15 | 1/5 | −0.180 | +0.021 | −0.005 | −0.056 | −0.008 |
| mean transfer + Feng + CD4 (v1) | 36 voters | 0.1537 | −0.1218 | 0/15 | 0/5 | −0.289 | −0.248 | −0.035 | −0.148 | −0.012 |
| refit + Feng + CD4 (v2) | 26 bios | 0.2815 | +0.0060 | 12/15 | 4/5 | −0.001 | +0.022 | +0.006 | +0.007 | +0.003 |
| refit + Feng + CD4 (v1) | 56 bios | 0.2185 | −0.0569 | 2/15 | 1/5 | −0.080 | −0.157 | −0.022 | −0.079 | −0.004 |
| FiLM, old 6 contexts (scrambled profiles, base refit) | 6 | 0.2935 | +0.0180 | 10/15 | 3/5 | −0.003 | +0.065 | +0.008 | +0.034 | +0.004 |
| FiLM, 35 contexts, base refit_v1, nested α | 35 | 0.2374 | −0.0381 | 0/15 | 0/5 | −0.108 | −0.083 | −0.010 | −0.031 | +0.002 |
| FiLM, 35 contexts, base refit_v1, α 1 | 35 | 0.2699 | −0.0056 | 3/15 | 1/5 | −0.091 | +0.051 | +0.003 | −0.003 | +0.006 |
| **FiLM, 35 contexts, base refit_v2, nested α** | 35 | **0.3001** | **+0.0247** | **15/15** | **5/5** | −0.014 | +0.097 | +0.014 | +0.045 | +0.006 |
| FiLM, 35 contexts, base refit_v2, α 1 | 35 | 0.3088 | +0.0333 | 12/15 | 4/5 | −0.044 | +0.143 | +0.021 | +0.071 | +0.009 |
| same, wrong context, α 1 | 35 | 0.2566 | −0.0189 | 3/15 | 1/5 | −0.038 | −0.024 | −0.008 | −0.045 | +0.002 |
| same, zero context, α 1 | 35 | 0.2189 | −0.0566 | 3/15 | 1/5 | −0.040 | −0.163 | −0.021 | −0.113 | −0.003 |

Artifacts: tables `table_v1/`, `table_v2/`; assets `ext_effects.npy`, `ext_gene_contexts.npz`;
deltas `<line>/<arm>.npy`; harness scores `scores.json`; diagnostics `baseline_diagnostics.json`,
`film_diagnostics.json`, `curve_diagnostics.json`, `film_curves.json`; official runs
`preset/runs/<line>/<arm>/seed*/{receipt.json,scored.csv,cell_eval2/}`; `coverage.json`,
`final_table.json`, `CONTROLS_COMMON_BUG.md`, `CRASH_RECOVERY_NOTES.md`. No official upload was
made in this experiment.
