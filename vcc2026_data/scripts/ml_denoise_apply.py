"""A-2 official inputs (2026-09-30): apply the SAME program-space denoising operator used for the harness library
(ml_denoise_sources.py: basis = top-K right singular vectors of each source's FULL screen, whitened; output = projection +
rho x residual) to another source directory's files (the official transfer_v2 sources), writing OUT_DIR with every
other file symlinked. Usage: ml_denoise_apply.py ORIG_DIR OUT_DIR K RHO SRC [SRC ...]"""
import sys, os, shutil, pathlib, json
import numpy as np
from sklearn.utils.extmath import randomized_svd
ORIG, OUT, K, RHO, SRCS = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve(), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5:]
# resolve(): relative ORIG made RELATIVE symlinks that do not resolve from OUT (2026-09-30 bug: official_dnK562 / official_dnALL
# had all 26 links broken, so the builder's de.src fusion saw K562 ONLY; see FINDINGS 22.6)
FULL = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data/reports/sc_hybrid_2026-09-25/denoise_src/full"); OUT.mkdir(parents=True, exist_ok=True)
rep = {}
for src in SRCS:
    zf = np.load(FULL / f"{src}.npz", allow_pickle=True); n = zf["ncell"].astype(float); keep = n >= 20
    fp = zf["fp"][keep].astype(np.float64); fc = zf["fc"][keep].astype(np.float64); n = n[keep]
    Y = np.log((fp * 1e6 + 50) / (fc * 1e6 + 50)); meas = np.isfinite(Y).mean(0) >= 0.9
    g = np.flatnonzero(meas); Y = Y[:, g]; mu = np.nanmean(Y, 0); Y = np.where(np.isfinite(Y), Y, mu[None, :]) - mu[None, :]
    sg = np.median(np.abs(np.sqrt(n)[:, None] * Y), 0) / 0.6745 + 1e-6; W = np.sqrt(n)[:, None] * Y / sg[None, :]
    _, _, Vt = randomized_svd(W, n_components=K, n_iter=5, random_state=0); del W, Y, fp, fc
    zl = dict(np.load(ORIG / f"{src}.npz", allow_pickle=True)); fpl = zl["fp"].astype(np.float64); fcl = zl["fc"].astype(np.float64); nl = zl["ncell"].astype(float)
    y = np.log((fpl[:, g] * 1e6 + 50) / (fcl[:, g] * 1e6 + 50)); fin = np.isfinite(y); yc = np.where(fin, y, mu[None, :]) - mu[None, :]
    w = np.sqrt(nl)[:, None] * yc / sg[None, :]; c = w @ Vt.T; wh = c @ Vt + RHO * (w - c @ Vt)
    yh = wh * sg[None, :] / np.sqrt(nl)[:, None] + mu[None, :]
    fpn = ((fcl[:, g] * 1e6 + 50) * np.exp(yh) - 50) / 1e6; clip = float((fpn < 0).mean()); fpn = np.where(fin, np.maximum(fpn, 0.0), fpl[:, g])
    fpo = fpl.copy(); fpo[:, g] = fpn; zl["fp"] = fpo.astype(np.float32); np.savez_compressed(OUT / f"{src}.npz", **zl)
    if (ORIG / f"{src}.json").exists(): shutil.copy(ORIG / f"{src}.json", OUT / f"{src}.json")
    rep[src] = dict(targets=int(len(nl)), genes=int(len(g)), clip=round(clip, 5)); print(OUT.name, src, rep[src], flush=True)
for f in ORIG.iterdir():
    if not (OUT / f.name).exists(): os.symlink(f.resolve(), OUT / f.name)
bad = [f.name for f in OUT.iterdir() if f.is_symlink() and not f.exists()]; assert not bad, f"broken links: {bad[:5]}"
json.dump(rep, open(OUT / "denoise_report.json", "w"), indent=1); print("DENOISE_APPLY_DONE", OUT, flush=True)
