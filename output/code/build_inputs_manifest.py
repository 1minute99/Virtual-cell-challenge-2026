"""Task 4 (plan zero-shot-limits-context-complexity): enumerate every input file of plan Tasks 6-21, grouped by task,
hash them (sha256) and write output/code/inputs_filelist.tsv + output/code/inputs.sha256.
Run:  ionice -c3 nice -n 10 python output/code/build_inputs_manifest.py   (from /home/wonmin-kim/AIVC)
Deterministic: file order = group order, then sorted path. No randomness."""
import glob, hashlib, os, pathlib, subprocess, sys
ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); V = "vcc2026_data"; RP = f"{V}/reports"; SC = f"{V}/scripts"
PY = os.path.expanduser("~/miniconda3/envs/vcc2026/bin/python")

def closure(*roots):
    out = subprocess.run([PY, str(ROOT / "output/code/import_closure.py"), *[str(ROOT / SC / r) for r in roots]],
                         capture_output=True, text=True, check=True).stdout.split()
    return [str(pathlib.Path(p).relative_to(ROOT)) for p in out]

def g(*pats):
    res = []
    for p in pats:
        hits = sorted(glob.glob(str(ROOT / p)))
        res += [str(pathlib.Path(h).relative_to(ROOT)) for h in hits if os.path.isfile(h)] if hits else [f"MISSING::{p}"]
    return res

PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
GROUPS = [
 ("T14,T17,T18 frozen pipeline S18/F2 (tv2_final.py + builder modules, specs, receipts)",
  closure("tv2_final.py", "tv2_submit_build.py", "tv2_submit_chain.sh", "tv2_tlrn_prep.py") +
  g(f"{RP}/final_2026-10-27/README.md", f"{RP}/final_2026-10-27/run_F1_spec.sh", f"{RP}/final_2026-10-27/run_F1_build.sh",
    f"{RP}/final_2026-10-27/run_F2_build.sh", f"{RP}/final_2026-10-27/F1/F1_spec.json", f"{RP}/final_2026-10-27/F1/F2_spec.json",
    f"{RP}/final_2026-10-27/F1/identity.json", f"{RP}/final_2026-10-27/F1/report.json", f"{RP}/final_2026-10-27/F1/targets.txt",
    f"{RP}/submission_2026-09-23/F2.sha256", f"{RP}/submission_2026-09-23/F2.upload_receipt.json", f"{RP}/submission_2026-09-23/F2.status.json",
    f"{RP}/submission_2026-09-23/S18.pred_build.json", f"{RP}/submission_2026-09-23/S18.status.json", f"{RP}/submission_2026-09-23/S18.sha256",
    f"{RP}/submission_2026-10-11/upload_s17d_s18.sh")),
 ("T9,T12-T18 harness truth (cell level) + harness config + real moments",
  g(f"{RP}/transfer_v2_2026-09-23/h1/real_thin.h5ad", f"{RP}/sc_hybrid_2026-09-25/orion_hct116/real_panel.h5ad",
    f"{RP}/sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad", f"{RP}/b19_sweep2_2026-09-18/full_metrics/KOLF/real_disjoint.h5ad",
    f"{V}/harness/kolf_real.h5ad", f"{V}/harness/jurkat_real.h5ad",
    *[f"{RP}/simple_residual_2026-09-16/mock/{l}/panel300/*" for l in ("K562", "HepG2", "Jurkat", "RPE1")],
    f"{V}/raw/scperturb/ReplogleWeissman2022_K562_gwps.h5ad", f"{V}/raw/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
    f"{V}/raw/scperturb/ReplogleWeissman2022_rpe1.h5ad", f"{V}/raw/scperturb/NadigOConner2024_hepg2.h5ad",
    f"{V}/raw/scperturb/NadigOConner2024_jurkat.h5ad",
    f"{RP}/final_harness_2026-09-24/panels.json", f"{RP}/final_harness_2026-09-24/neighbours.json",
    f"{RP}/final_harness_2026-09-24/targets.txt",
    *[f"{RP}/final_harness_2026-09-24/cache_real_{p}/pseudobulk_moments_bulk_lognorm-*.moments.npz" for p in PANELS7 if p != "H1"],
    f"{RP}/transfer_v2_2026-09-23/cache_real_H1/pseudobulk_moments_bulk_lognorm-*.moments.npz")),
 ("T12-T17 ERA offline packs + exact evaluator (sch_sumeval, ERA evaluator)",
  g(f"{RP}/era_2026-10-07/packs/*") + closure("sch_sumeval.py") +
  g("external/era-aivc/aivc/tasks/vcc_sum/evaluator.py", "external/era-aivc/aivc/tasks/vcc_sum/prep.py",
    "external/era-aivc/aivc/tasks/vcc_sum/task.py", "external/era-aivc/aivc/tools/ctrl_pcs.py",
    "external/era-aivc/aivc/tools/ctrl_var_prep.py", "external/era-aivc/aivc/tools/stack_screen.py")),
 ("T9,T11,T14,T18 M1 source registry (source/family npz) + pack-builder inputs + S4 parts",
  g(f"{RP}/transfer_v2_2026-09-23/sources/*", f"{RP}/transfer_v2_2026-09-23/sources_new/*",
    f"{RP}/final_harness_2026-09-24/sources/*", f"{RP}/sc_hybrid_2026-09-25/h1fuse/*_data.npz",
    f"{RP}/s4_improve_2026-09-27/K562/film_noKOLFFENG.npy",
    f"{RP}/s4_improve_2026-09-27/HepG2/film_noKOLFFENG.npy", f"{RP}/s4_improve_2026-09-27/Jurkat/film_noKOLFFENG.npy",
    f"{RP}/s4_improve_2026-09-27/RPE1/film_noKOLFFENG.npy", f"{RP}/s4_improve_2026-09-27/H1/film_orig.npy",
    f"{RP}/s4_improve_2026-09-27/HCT116/film_orXnoKOLFFENG.npy", f"{RP}/s4_improve_2026-09-27/HCT116/film_orRnoKOLFFENG.npy",
    f"{RP}/ext_contexts_2026-09-21/KOLF/film_v2_a0.25.npy",
    f"{V}/raw/feng/*Cell-Metadata.tsv.gz") + closure("ml_h1fuse.py", "fh_harness.py", "sch_emit.py", "s4i_refit_excl.py", "s4i_refit_h1.py")),
 ("T6-T8 Tahoe context curve (dcp protocol, seed-0 results)",
  g(f"{RP}/simple_residual_2026-09-16/data/tahoe/table.npz", f"{V}/data/processed/drug_context_pretrain/*",
    f"{RP}/drug_context_pretrain_2026-09-20/context_splits.json", f"{RP}/drug_context_pretrain_2026-09-20/scaling_table.json",
    f"{RP}/drug_context_pretrain_2026-09-20/FINDINGS.md", f"{RP}/drug_context_pretrain_2026-09-20/curves.json",
    f"{RP}/drug_context_pretrain_2026-09-20/final_table.json") +
  closure("dcp_data.py", "dcp_model.py", "dcp_prepare.py", "dcp_select.py")),
 ("T10 learned-architecture replication (arch screen)",
  closure("arch_zoo.py") + g(f"{RP}/arch_screen_2026-09-22/*.json", f"{RP}/arch_screen_2026-09-22/FINDINGS.md",
    f"{RP}/arch_screen_2026-09-22_HepG2/*.json", f"{RP}/ext_contexts_2026-09-21/FINDINGS.md")),
 ("T11,T18 source-count curve + H1GW value",
  closure("ml_ctx_curve.py", "ml_h1gw_value.py") +
  g(f"{RP}/transfer_v2_2026-09-23/sources_new/H1GW.npz", f"{RP}/transfer_v2_2026-09-23/sources_new/H1GW.json",
    f"{RP}/sc_hybrid_2026-09-25/ctx_curve.json", f"{RP}/sc_hybrid_2026-09-25/err_decomp/h1gw_value.json")),
 ("T13,T19,T20,T22+ reports named in the plan",
  g(f"{RP}/arch_screen_2026-09-22/FINDINGS.md", f"{RP}/ext_contexts_2026-09-21/FINDINGS.md", f"{RP}/sc_hybrid_2026-09-25/FINDINGS.md",
    f"{RP}/PROJECT_HISTORY.md", f"{RP}/drug_gene_concordance_2026-09-16/*")),
 ("T19 official record (scores, receipts, leaderboard snapshots for aggregate quantiles)",
  g(f"{RP}/official_scores.csv", f"{RP}/simple_residual_2026-09-16/submission/official_scores_new.csv",
    f"{RP}/context_forcing_strength_2026-09-16/submission/official_scores_new.csv",
    f"{RP}/submission_2026-09-20/S1_status.log", f"{RP}/submission_2026-09-20/S2_status.log",
    f"{RP}/submission_2026-09-21/status.json", f"{RP}/submission_2026-09-22/status.json",
    f"{RP}/submission_2026-09-23/*.status.json", f"{RP}/submission_2026-09-23/*.upload_receipt.json",
    f"{RP}/submission_receipts/*.json", f"{RP}/leaderboard_*.json")),
]

seen, rows, missing = set(), [], []
for gname, files in GROUPS:
    for f in files:
        if f.startswith("MISSING::"):
            missing.append((gname, f[9:])); continue
        if "/__pycache__/" in f: continue
        rows.append((gname, f, f in seen)); seen.add(f)

def sha(p, bs=1 << 23):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(bs)
            if not b: break
            h.update(b)
    return h.hexdigest()

hashes, prev = {}, {}
SHAF = ROOT / "output/code/inputs.sha256"
if SHAF.exists():                       # reuse a hash only if the file is older than the previous sha file
    t0 = SHAF.stat().st_mtime
    for line in SHAF.read_text().splitlines():
        h, f = line.split("  ", 1)
        if (ROOT / f).exists() and (ROOT / f).stat().st_mtime < t0: prev[f] = h
for i, f in enumerate(sorted(seen)):
    hashes[f] = prev.get(f) or sha(ROOT / f)
    if i % 50 == 0: print(f"[{i}/{len(seen)}] {f}", flush=True)
with open(ROOT / "output/code/inputs.sha256", "w") as fh:
    for f in sorted(seen): fh.write(f"{hashes[f]}  {f}\n")
with open(ROOT / "output/code/inputs_filelist.tsv", "w") as fh:
    fh.write("group\tpath\tbytes\tsha256\n")
    for gname, f, dup in rows:
        fh.write(f"{gname}\t{f}\t{os.path.getsize(ROOT / f)}\t{hashes[f]}\n")
    for gname, f in missing:
        fh.write(f"{gname}\t{f}\tMISSING\tMISSING\n")
print("FILES", len(seen), "BYTES", sum(os.path.getsize(ROOT / f) for f in seen), "MISSING", missing)
