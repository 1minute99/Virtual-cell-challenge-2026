"""Final phase: DepMap-neighbour estimates for a NEW target panel (feeds the builder's `blend` key).

  tv2_final_knn.py TARGETS_TXT OUT_DIR

1. DepMap 23Q4 co-essentiality neighbours (top 20, corr > 0) of every panel target among genes some public screen
   perturbed (as fh_prep.py).
2. Source library (tv2 format, RN HepG2/Jurkat/RPE1 included) for panel targets + their top-3 neighbours
   (as fh_sources.py) in OUT_DIR/sources; existing files are reused.
3. Recipe delta of every library target (fh_harness base recipe, all families), own knockdown column zeroed;
   neighbour estimate = corr^2-weighted mean of the top-10 covered neighbours (t itself excluded).
Writes OUT_DIR/knn.npz: targets, knn (n x 18533), has (neighbour estimate exists), covered (t has >= 20 cells in
some source). UNCOVERED targets (the builder gives them the neighbour estimate only) use CONCATENATED-prior neighbours,
similarity 0.75 * DepMap cosine + 0.25 * STRING v12 network-embedding cosine (2026-09-27: harness unseen fallback
H1 0.136 -> 0.156, KOLF 0.028 -> 0.041; neighbour-estimate quality +9-18 % on K562 / HCT116). Covered targets keep DepMap. Validated on the six-panel harness: blend 0.3 (energy-matched) 18/18 vs base; unseen-target fallback
PDS 0.59-0.70 vs 0.51.
"""
import json, pathlib, sys
import numpy as np, pandas as pd
from cell_eval2.moments import jackknife_correction  # noqa: F401
_argv = list(sys.argv); sys.argv = ["x"]
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_sources as TS

R = TS.R; canon = TS.canon
tg_file, out = pathlib.Path(_argv[1]), pathlib.Path(_argv[2]); out.mkdir(parents=True, exist_ok=True)
EXCL = set(x for x in (_argv[3].split(",") if len(_argv) > 3 else []) if x)   # optional families left out of the library
                                                                             # (S10k: somatic destinations drop KOLF; 2026-09-29)
targets = [t for t in tg_file.read_text().split() if t]

# 1. neighbours
pert = set()
for store in ["orion_target_cellline_all", "cd4_target_condition_all", "kolf_target"]:
    o = pd.read_csv(R / "pseudobulk" / store / "obs.csv", usecols=["target_symbol"], low_memory=False)
    pert |= {canon(x) for x in o.target_symbol.astype(str).unique()}
z = np.load(R / "reports/simple_residual_2026-09-16/data/replogle_nadig/K562/pseudobulk.npz", allow_pickle=True)
pert |= {canon(str(x)) for x in z["targets"]}
dm = np.load(R / "processed/depmap_23q4/gene_effects.npz", allow_pickle=True)
g = np.array([canon(str(x)) for x in dm["genes"]]); V = dm["values"].astype(np.float64)
V = V - np.nanmean(V, 1, keepdims=True); V = np.nan_to_num(V); V /= np.linalg.norm(V, axis=1, keepdims=True) + 1e-12
gi = {x: i for i, x in enumerate(g)}; cand = np.array([i for i, x in enumerate(g) if x in pert])
nb = {}
for t in targets:
    i = gi.get(canon(t))
    if i is None:
        continue
    c = V[cand] @ V[i]; o = np.argsort(-c)
    nb[canon(t)] = [(g[cand[j]], float(c[j])) for j in o[:30] if g[cand[j]] != canon(t)][:20]
(out / "neighbours.json").write_text(json.dumps(nb))
# concatenated DepMap + STRING neighbours (for uncovered targets)
import h5py
SA = 0.25; S_DIR = R / "raw/string"
_info = pd.read_csv(S_DIR / "9606.protein.info.v12.0.txt.gz", sep="\t"); _nm = dict(zip(_info["#string_protein_id"], _info["preferred_name"]))
_f = h5py.File(S_DIR / "9606.protein.network.embeddings.v12.0.h5", "r")
_sg = np.array([canon(str(_nm.get(p.decode() if isinstance(p, bytes) else str(p), ""))) for p in _f["proteins"][:]])
_S = _f["embeddings"][:].astype(np.float64); _S -= _S.mean(0); _S /= np.linalg.norm(_S, axis=1, keepdims=True) + 1e-12
_si = {}
for i, x in enumerate(_sg):
    _si.setdefault(x, i)
cand_g = g[cand]; cand_s = np.array([_si.get(x, -1) for x in cand_g])
nbc = {}
for t in targets:
    i = gi.get(canon(t)); j = _si.get(canon(t))
    if i is None and j is None:
        continue
    cD = V[cand] @ V[i] if i is not None else np.full(len(cand), np.nan)
    cS = np.where(cand_s >= 0, _S[np.maximum(cand_s, 0)] @ _S[j], np.nan) if j is not None else np.full(len(cand), np.nan)
    c = np.where(np.isfinite(cD) & np.isfinite(cS), (1 - SA) * cD + SA * cS, np.where(np.isfinite(cD), cD, cS))
    c = np.nan_to_num(c, nan=-1.0); o = np.argsort(-c)
    nbc[canon(t)] = [(cand_g[k], float(c[k])) for k in o[:30] if cand_g[k] != canon(t)][:20]
(out / "neighbours_concat.json").write_text(json.dumps(nbc))
# top-3 per target (+ incidental members through other targets) on purpose: the harness library that validated blend 0.3
# held on average 5.45 of each official target's top-10 neighbours, this construction 5.02 on A/B/C, a top-10 library
# ~10 -> a different, unvalidated estimator (Codex R15 flagged the 3 vs 10 asymmetry; measured 2026-09-29)
need = sorted({canon(t) for t in targets} | {x for t in nb for x, _ in nb[t][:3]} | {x for t in nbc for x, _ in nbc[t][:3]})
print("targets", len(targets), "with DepMap neighbours", len(nb), "library targets", len(need), flush=True)

# 2. source library
TS.OUT = out / "sources"; TS.OUT.mkdir(exist_ok=True); nset = set(need)
jobs = [("K562", lambda: TS.from_rn("K562", "K562", nset)), ("RN_HepG2", lambda: TS.from_rn("RN_HepG2", "HepG2", nset)),
        ("RN_Jurkat", lambda: TS.from_rn("RN_Jurkat", "Jurkat", nset)), ("RN_RPE1", lambda: TS.from_rn("RN_RPE1", "RPE1", nset)),
        ("HCT116", lambda: TS.from_store("HCT116", "orion_target_cellline_all", nset, ("cell_context", "HCT116"))),
        ("HEK293T", lambda: TS.from_store("HEK293T", "orion_target_cellline_all", nset, ("cell_context", "HEK293T"))),
        *[(f"CD4_{s}", (lambda s=s: TS.from_store(f"CD4_{s}", "cd4_target_condition_all", nset, ("stimulation", s))))
          for s in ["Rest", "Stim8hr", "Stim48hr"]],
        ("KOLF", lambda: TS.from_store("KOLF", "kolf_target", nset)),
        ("FENG_T", lambda: TS.from_store("FENG_T", "feng_targetedscreen_target", nset)),
        ("FENG_GWF", lambda: TS.from_store("FENG_GWF", "feng_genomewidescreen_fitnessgenes_target", nset)),
        ("FENG_GWN", lambda: TS.from_store("FENG_GWN", "feng_genomewidescreen_nonfitnessgenes_target", nset)),
        ("H1", lambda: TS.from_store("H1", "vcc2025_target_batch", nset, batch_matched=True))]
for name, f in jobs:
    if not (TS.OUT / f"{name}.npz").exists():
        f()

# 2b. optional (2026-09-30, S14fk DE side): program-denoise one source of the library before the neighbour estimates.
#     TV2_KNN_DENOISE="K562:100:0.25" -> ml_denoise_apply.py on this library (basis from the FULL screen, panel-agnostic).
import os, subprocess
_DN = os.environ.get("TV2_KNN_DENOISE", ""); LIBSRC = TS.OUT
if _DN:
    _src, _k, _rho = _DN.split(":"); LIBSRC = out / f"sources_dn{_src}"
    if not (LIBSRC / f"{_src}.npz").exists():
        subprocess.run([sys.executable, "-P", str(R / "scripts/ml_denoise_apply.py"), str(TS.OUT), str(LIBSRC), _k, _rho, _src], check=True)
    print("neighbour library from the denoised copy", LIBSRC, flush=True)

# 3. neighbour estimates with the harness recipe, repointed at this library
import fh_harness as H
H.TM.SRC = LIBSRC; H.TM._cache.clear(); H.CACHE = out / ("cache_dn" if _DN else "cache"); H.CACHE.mkdir(exist_ok=True)
H.LIB = sorted({str(t) for s in H.TM.FAMILY if (TS.OUT / f"{s}.npz").exists() for t in np.load(TS.OUT / f"{s}.npz")["targets"]})
dlib, _ = H.base(H.LIB, EXCL, "lib_all" + ("_ex" + "-".join(sorted(EXCL)) if EXCL else ""))
for i, t in enumerate(H.LIB):
    j = H.GIC.get(canon(t))
    if j is not None:
        dlib[i, j] = 0.0
cov = np.abs(dlib).sum(1) > 0; pos = {canon(t): i for i, t in enumerate(H.LIB)}
covered = np.zeros(len(targets), bool)
for s in H.TM.FAMILY:
    if (TS.OUT / f"{s}.npz").exists() and H.TM.FAMILY[s] not in EXCL:
        zz = H.TM.load(s)
        covered |= np.array([(zz["pos"].get(canon(t)) is not None) and zz["ncell"][zz["pos"][canon(t)]] >= 20 for t in targets])
knn = np.zeros((len(targets), H.NG), np.float32); has = np.zeros(len(targets), bool)
for k, t in enumerate(targets):
    src = nb if covered[k] else nbc                  # uncovered targets: concatenated DepMap + STRING neighbours
    cands = [(pos[x], c) for x, c in src.get(canon(t), []) if x in pos and x != canon(t) and cov[pos[x]] and c > 0][:10]
    if cands:
        w = np.array([c for _, c in cands]) ** 2
        knn[k] = (w[:, None] * dlib[[i for i, _ in cands]]).sum(0) / w.sum(); has[k] = True
np.savez(out / "knn.npz", targets=np.array(targets), knn=knn, has=has, covered=covered, excluded=np.array(sorted(EXCL)))
print(json.dumps(dict(targets=len(targets), with_neighbours=int(has.sum()), covered=int(covered.sum()),
                      uncovered=[t for t, c in zip(targets, covered) if not c])), flush=True)
