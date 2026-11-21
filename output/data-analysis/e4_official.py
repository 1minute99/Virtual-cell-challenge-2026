"""Task 19 of plan zero-shot-limits-context-complexity (input/ideas/zero-shot-limits-context-complexity-plan.md; A3).
Official VCC 2026 record: score trajectory, final validation rank, gains attributable to single-change submission
pairs, local -> official shrinkage, leaderboard snapshots in aggregate only.

Writes results/e4_official.csv (one row per official submission, with official parent and attribution class),
results/e4_official.json (summary), results/e4_official_trajectory.png.

VCC RULES (input/notes/external-gates.md section 1): only OUR scores and derived statistics are reported. Leaderboard
snapshots are reduced to aggregate quantiles + our rank; no other team's name, id, organisation or row is written.

ASSUMPTIONS (stated; checked below, outcomes in e4_official.json["checks"]):
  E1  official_record.csv (Task 4) is the complete list: 38 entries = 36 validation + 2 final (A3). B10/B11/B14 have
      no entry id on disk; their score/rank come from PROJECT_HISTORY.md only.
  E2  The composite score_avg = mean of the six component scores (pds, mse, nmae, fid, reach, jac). Checked per row.
  E3  Submission order: the B-series/b18/b19 in the PROJECT_HISTORY.md "Official submissions, all 14" order
      (no timestamps on disk), then S/F entries by submission_date.
  E4  A leaderboard snapshot shows one row per team = that team's LAST submission (not its best). Our row is found
      by our entry ids; rank = 1 + #(teams with strictly higher score_avg). Checked: rows == distinct teams, and our
      displayed entry == our latest validation entry at snapshot time.
  E5  "rank at submission" in official_record.csv is the rank on the day the score was retrieved; the board grew
      from 1,036 to 1,635 teams, so ranks across submissions are not comparable; score is.
  E6  Attribution: a gain is attributed to a component only if the child differs from its OFFICIAL parent (the
      nearest ancestor that was itself uploaded and scored) by ONE lever. Lever counts come from a recursive diff
      of the two build specs (`<label>.pred_build.json`, key `spec`, notes excluded) where both exist, otherwise
      from the submission notes / AGENTS.md / PROJECT_HISTORY.md (marked `evidence = notes`).
  E7  Local -> official shrinkage ratio = official delta / local delta for the SAME contrast vs the refit S1.
No randomness (deterministic)."""
import ast, glob, json, pathlib, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); RP = ROOT / "vcc2026_data/reports"
OUT = ROOT / "output/data-analysis/results"; SUB = RP / "submission_2026-09-23"
COMP = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]
checks = {}

rec = pd.read_csv(OUT / "official_record.csv")
rec["name"] = rec["s_label"].fillna(rec["label"])
checks["E1_n_total"] = int(len(rec)); checks["E1_n_validation"] = int((rec.phase == "validation").sum())
checks["E1_n_final"] = int((rec.phase == "final").sum()); checks["E1_no_id"] = sorted(rec.loc[rec.entry_id.str.startswith("NO_ID"), "name"])
full = rec.dropna(subset=COMP)
checks["E2_max_abs_avg_minus_mean6"] = float((full[COMP].mean(axis=1) - full.score_avg).abs().max())

# ---- E3 order ------------------------------------------------------------------------------------------------
HIST = ["이전A", "이전B", "B1", "B2", "B6", "B7", "B8", "B10", "B11", "B14", "B15", "B16", "b18", "b19"]
s = rec[~rec.name.isin(HIST)].sort_values("submission_date")
rec = pd.concat([rec.set_index("name").loc[HIST].reset_index(), s]).reset_index(drop=True)
rec["order"] = np.arange(1, len(rec) + 1)
rec["display"] = rec["name"].replace({"이전A": "priorA", "이전B": "priorB"})

# ---- E6 lineage: child -> (official parent, evidence, levers, description) ------------------------------------
# Official parent = nearest scored ancestor. Unscored intermediates (S9, S10i, S17) are named in `via`.
L = {
 "B2":   ("B1",  "notes", 1, "CD4 alone -> 3-source mean"),
 "B6":   ("B2",  "notes", 2, "3 -> 6 sources AND deterministic -> stochastic rounding"),
 "B7":   ("B6",  "notes", 1, "6-source mean -> sign-agreement weighting"),
 "B8":   ("B6",  "notes", 1, "delta x3 amplitude scale"),
 "B10":  ("B8",  "notes", 2, "plain mean transfer AND scale x3 -> x4 (stochastic rounding status not recorded)"),
 "B11":  ("B8",  "notes", 2, "expression-gated per-context AND scale x2.25"),
 "B14":  ("B8",  "notes", 2, "precision-weighted sources AND scale x3.5"),
 "B15":  ("B8",  "notes", 1, "effective-cell-count source weighting (scale 3 kept)"),
 "B16":  ("B8",  "notes", 1, "frozen rank-128 covariance denoiser (scale 3 kept)"),
 "b18":  ("B16", "notes", 0, "byte-identical package to B16 (scorer determinism check)"),
 "b19":  ("b18", "notes", 1, "context-forcing residual (rank 32, alpha 1, x3)"),
 "S2":   ("S1",  "notes", 1, "refit + PCA interaction; training context profiles SCRAMBLED (CONTROLS_COMMON_BUG.md)"),
 "S3":   ("S1",  "notes", 1, "refit + 0.25 x scratch FiLM residual; training context vectors SCRAMBLED"),
 "S4":   ("S1",  "notes", 2, "extended-table refit_v2 AND 0.25 x FiLM residual on 38 gene contexts (submission_2026-09-22/FINDINGS.md)"),
 "S4gp": ("S4",  "notes", 3, "deltas x3 AND gated cellshift AND promoter prior"),
 "S4mix":("S4gp","spec", 2, "0.5 x S4 + 0.5 x transfer_v2 fusion x1.6 (new part + weight)"),
 "S7a":  ("S4mix","spec", 1, "fusion share 0.5 -> 0.75 (part weights)"),
 "S8":   ("S4mix","spec", 1, "DepMap co-essentiality neighbour blend w 0.2"),
 "S10b": ("S8",  "spec", 2, "blend w 0.2 -> 0.3 (via unsubmitted S9) AND dual-moment emission"),
 "S10c": ("S10b","spec", 1, "cpc composition pin"),
 "S10e": ("S10c","spec", 1, "S4 half replaced by equal-weight fold-change part (one part swapped)"),
 "S10f": ("S10c","spec", 3, "blend w 0.5 AND re-allocation K 5000 AND down-gene share 0.2"),
 "S10k": ("S10f","spec", 2, "via unsubmitted S10i (kink budget, sum-side centre/ew, blend/K reverted) AND no-KOLF sum side + learned residual"),
 "S11":  ("S10k","spec", 3, "min_cells ensemble AND ci5 fixed point AND depth-balanced block"),
 "S12":  ("S11", "spec", 1, "DE side from the sum-side recipe override"),
 "S14fk":("S12", "spec", 2, "S4 part without FENG AND K562-denoised DE source"),
 "S14fk2":("S14fk","spec",1, "DE source symlink fix (intended all-family K562-denoised DE side)"),
 "S15":  ("S14fk2","spec",1, "M2: agreement amplitude (tamp beta 0.5)"),
 "S16":  ("S15", "spec", 1, "learned per-target amplitude (tamp file swap, beta 0.75)"),
 "S17u": ("S16", "spec", 1, "sum-budget probe UP (kink margin 1.2, floor 1; ew key removed but ew already inactive in S16 build log)"),
 "S17d": ("S16", "spec", 1, "sum-budget probe DOWN (kink margin 0.6, floor 0; same ew note)"),
 "S18":  ("S16", "spec", 4, "via unsubmitted S17: +4 source families AND LFC-pseudocount ensemble AND DE calibration drc, THEN M3 control-PC deflation (k8, eta 0.75)"),
}

def spec(n):
    p = SUB / f"{n}.pred_build.json"
    if not p.exists(): return None
    sp = json.loads(p.read_text())["spec"]
    return ast.literal_eval(sp) if isinstance(sp, str) else sp

def flat(x, p=""):
    o = {}
    if isinstance(x, dict):
        for k, v in x.items(): o.update(flat(v, f"{p}.{k}"))
    elif isinstance(x, list):
        for i, v in enumerate(x): o.update(flat(v, f"{p}[{i}]"))
    else: o[p] = x
    return o

def specdiff(a, b):
    A, B = spec(a), spec(b)
    if A is None or B is None: return None
    A, B = flat(A), flat(B)
    keys = [k for k in set(A) | set(B) if A.get(k) != B.get(k) and k != ".note"]
    return {"n_keys": len(keys), "blocks": sorted({re.split(r"[.\[]", k)[1] for k in keys})}

# spec-level checks of the pairs named in the task, including the unscored S17 step
for a, b in [("S14fk2", "S15"), ("S15", "S16"), ("S16", "S17"), ("S17", "S18"), ("S16", "S18"), ("S16", "S17u"), ("S16", "S17d")]:
    checks[f"specdiff_{a}_{b}"] = specdiff(a, b)
checks["S17_officially_scored"] = bool((rec.name == "S17").any())
checks["S17_note"] = flat(spec("S17")).get(".note", "")[:400]

by = rec.set_index("name")
rows = []
for _, r in rec.iterrows():
    d = {"order": int(r.order), "submission": r.display, "entry_id": r.entry_id, "phase": r.phase,
         "submission_date": r.submission_date, "score_avg": r.score_avg, "rank_at_retrieval": r["rank"],
         "best_so_far": np.nan, "note": r.note}
    for c in COMP: d[c] = r[c]
    if r["name"] in L:
        par, ev, nlev, desc = L[r["name"]]
        d.update(official_parent=by.loc[par, "display"], evidence=ev, n_levers=nlev, change=desc,
                 delta_avg=r.score_avg - by.loc[par, "score_avg"])
        for c in COMP: d["d_" + c] = r[c] - by.loc[par, c]
        sd = specdiff(par, r["name"]) if ev == "spec" else None
        d["spec_keys_changed"] = sd["n_keys"] if sd else np.nan
        d["spec_blocks_changed"] = ";".join(sd["blocks"]) if sd else ""
        d["lineage_confidence"] = "high (build-spec diff)" if ev == "spec" else ("low (parent inferred from 1-line notes)" if r["name"] in ("B15", "B16", "B10", "B11", "B14") else "medium (notes/FINDINGS)")
        d["attribution"] = ("determinism_check" if nlev == 0 else "single_change" if nlev == 1 else "MULTI_CHANGE_not_attributable")
        if r["name"] in ("S2", "S3"): d["attribution"] += "|INVALID_scrambled_context_vectors"
    else:
        d.update(official_parent="", evidence="", n_levers=np.nan, change="lineage root / pipeline restart (not attributable)",
                 attribution="root", delta_avg=np.nan)
    rows.append(d)
T = pd.DataFrame(rows)
T["best_so_far"] = T.score_avg.where(T.phase == "validation").cummax()
T.to_csv(OUT / "e4_official.csv", index=False)

# ---- leaderboard snapshots (aggregate only) -------------------------------------------------------------------
ours_ids = set(rec.entry_id); snaps = []
for f in sorted(glob.glob(str(RP / "leaderboard_*.json"))):
    e = json.load(open(f))["entries"]
    sc = np.array([x["score_avg"] for x in e if x["score_avg"] is not None], float)
    mine = [x for x in e if x["id"] in ours_ids]
    assert len(mine) == 1, f
    m = mine[0]; ours = float(m["score_avg"]); rank = int((sc > ours).sum() + 1)
    date = re.search(r"(\d{4}-\d{2}-\d{2}b?)", f).group(1)
    snaps.append({"snapshot": date, "n_teams": len(e), "n_distinct_teams": len({x["team_id"] for x in e}),
                  "our_displayed_entry": by.reset_index().set_index("entry_id").loc[m["id"], "name"],
                  "our_score": ours, "our_rank": rank, "our_rank_fraction": rank / len(sc),
                  "our_percentile": 100 * (sc < ours).mean(), "our_submission_count": m.get("submission_count"),
                  "top1": float(sc.max()), "top10_threshold": float(np.sort(sc)[::-1][9]),
                  "q99": float(np.quantile(sc, .99)), "q90": float(np.quantile(sc, .90)),
                  "q75": float(np.quantile(sc, .75)), "median": float(np.median(sc)), "q25": float(np.quantile(sc, .25))})
snaps = sorted(snaps, key=lambda x: x["snapshot"])
checks["E4_rows_equal_distinct_teams"] = all(s_["n_teams"] == s_["n_distinct_teams"] for s_ in snaps)
fin = snaps[-1]
checks["final_rank_verified"] = {"rank": fin["our_rank"], "n": fin["n_teams"], "agents_md": "S18 0.2545 ended rank 129; top-10 >= 0.319",
                                 "match": fin["our_rank"] == 129 and fin["n_teams"] == 1635 and fin["our_displayed_entry"] == "S18"}
checks["final_submission_count_on_board"] = fin["our_submission_count"]

# ---- local -> official shrinkage (E7) -------------------------------------------------------------------------
S1 = by.loc["S1", "score_avg"]
shrink = {
  "S4": {"local_delta": 0.0247, "local_source": "ext_contexts_2026-09-21/FINDINGS.md via AGENTS.md: film_v2 on refit_v2 vs refit, five-line local composite, 15/15 paired seeds",
         "official_delta": float(by.loc["S4", "score_avg"] - S1), "contrast": "S4 vs S1 refit (bundle: refit_v2 +0.0060 locally on its own, plus FiLM)",
         "valid": True},
  "S3": {"local_delta": 0.0180, "local_source": "drug_context_pretrain_2026-09-20/FINDINGS.md via AGENTS.md: scratch FiLM residual vs refit, five-line local composite",
         "official_delta": float(by.loc["S3", "score_avg"] - S1), "contrast": "S3 vs S1 refit",
         "valid": False, "why_invalid": "both the local +0.0180 and the S3 model used the scrambled gene-context vectors (CONTROLS_COMMON_BUG.md; FINDINGS caveat); the local number comes from the invalid drug-pretraining experiment"},
}
for v in shrink.values(): v["ratio_official_over_local"] = v["official_delta"] / v["local_delta"]
checks["S4_official_delta_matches_FINDINGS_+0.001310"] = abs(shrink["S4"]["official_delta"] - 0.001310) < 1e-6
checks["S3_official_delta_matches_-0.000501"] = abs(shrink["S3"]["official_delta"] + 0.000501) < 1e-6

single = T[T.attribution.str.startswith("single_change")][["submission", "official_parent", "change", "lineage_confidence", "attribution", "delta_avg"] + ["d_" + c for c in COMP]]
summary = {
  "plan": "input/ideas/zero-shot-limits-context-complexity-plan.md", "task": 19,
  "n_submissions": checks["E1_n_total"], "n_validation": checks["E1_n_validation"], "n_final": checks["E1_n_final"],
  "first_score": float(T.score_avg.iloc[0]), "best_validation": {"submission": "S18", "score": float(by.loc["S18", "score_avg"])},
  "final_validation_standing": {k: fin[k] for k in ["snapshot", "n_teams", "our_rank", "our_rank_fraction", "our_percentile", "top10_threshold", "median", "q90", "q99"]},
  "gap_to_top10": fin["top10_threshold"] - fin["our_score"],
  "M2": {"pair": "S14fk2 -> S15", "single_change": True, "delta_avg": float(by.loc["S15", "score_avg"] - by.loc["S14fk2", "score_avg"]),
         "delta_components": {c: float(by.loc["S15", c] - by.loc["S14fk2", c]) for c in COMP}},
  "M3_bundle": {"pair_requested": "S17 -> S18", "S17_officially_scored": False,
                "official_pair": "S16 -> S18", "delta_avg": float(by.loc["S18", "score_avg"] - by.loc["S16", "score_avg"]),
                "delta_components": {c: float(by.loc["S18", c] - by.loc["S16", c]) for c in COMP},
                "changes": ["+HEPG2/RNJURKAT/RPE1/DLD1 source families (weight 1)", "LFC-pseudocount ensemble members (c 5/10/20/100)",
                            "DE magnitude calibration drc 0.8/1.1", "M3 control-PC deflation of weaker half (k 8, eta 0.75)"],
                "S17_local_harness_bundle_R21": 0.0055,
                "verdict": "MULTI-CHANGE: the official +gain cannot be attributed to M3 alone; S17->S18 is a single change at spec level but S17 has no official score"},
  "single_change_pairs": single.to_dict("records"),
  "shrinkage": shrink,
  "final_phase": {"F1": "uploaded 2026-10-27, superseded", "F2": "uploaded 2026-10-28, last (counting) final entry",
                  "scores": "not yet published (rules: final leaderboard released when winners are announced, late November 2026)"},
  "leaderboard_snapshots_aggregate": snaps, "checks": checks,
}
(OUT / "e4_official.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=float))

# ---- figure: two panels, one axis each (score trajectory; rank fraction over snapshots) -----------------------
V = T[T.phase == "validation"]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={"width_ratios": [2.2, 1]})
a = ax[0]
a.plot(V.order, V.score_avg, "o", ms=4, color="#2a78d6", label="submission")
a.step(V.order, V.best_so_far, where="post", lw=2, color="#eb6834", label="best so far")
for nm in ["B8", "b19", "S4", "S10c", "S15", "S18", "S17u"]:
    rr = V[V.submission == nm].iloc[0]; a.annotate(nm, (rr.order, rr.score_avg), textcoords="offset points", xytext=(3, 5), fontsize=7, color="#52514e")
a.axhline(fin["top10_threshold"], ls=":", lw=1, color="#52514e"); a.text(1, fin["top10_threshold"] + .004, "top-10 threshold (final validation board)", fontsize=7, color="#52514e")
a.axhline(fin["median"], ls=":", lw=1, color="#52514e"); a.text(1, fin["median"] + .004, "median team (final validation board)", fontsize=7, color="#52514e")
a.set_xlabel("official validation submission (order)"); a.set_ylabel("official composite score"); a.legend(frameon=False, fontsize=8, loc="lower right")
a.spines[["top", "right"]].set_visible(False); a.grid(axis="y", lw=.3, alpha=.5)
b = ax[1]; x = np.arange(len(snaps))
b.plot(x, [100 * s_["our_rank_fraction"] for s_ in snaps], "o-", lw=2, ms=5, color="#2a78d6")
b.set_xticks(x); b.set_xticklabels([s_["snapshot"][5:] for s_ in snaps], rotation=60, fontsize=7)
b.set_ylabel("our rank / teams (%)  (lower = better)"); b.invert_yaxis(); b.spines[["top", "right"]].set_visible(False); b.grid(axis="y", lw=.3, alpha=.5)
b.set_title("displayed entry = latest, not best", fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "e4_official_trajectory.png", dpi=150)
print(json.dumps({k: summary[k] for k in ["n_submissions", "best_validation", "final_validation_standing", "gap_to_top10", "M2"]}, indent=1, default=float))
print("M3 bundle S16->S18:", summary["M3_bundle"]["delta_avg"]); print("shrink:", {k: round(v["ratio_official_over_local"], 4) for k, v in shrink.items()})
print(single[["submission", "official_parent", "delta_avg"]].to_string()); print(json.dumps(checks, indent=1, default=str, ensure_ascii=False)[:3000])
