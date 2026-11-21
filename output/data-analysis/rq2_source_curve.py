"""Task 11 (plan `input/ideas/zero-shot-limits-context-complexity-plan.md`, RQ2 contrast + RQ4 input; DESCRIPTIVE, no
hypothesis): transfer-pipeline source-count curve PDS(k) / MSE(k) over the number k of source FAMILIES, extending
`vcc2026_data/scripts/ml_ctx_curve.py` (AGENTS.md 2026-09-29: PDS ~ 0.796 + 0.048 ln k on H1, older S10 two-part recipe).

Subcommands (vcc2026 env, `python -P`; nice 10; OMP/MKL/OpenBLAS 4 threads; one process <= 8 GB):
  libfam        per-family, per-recipe-part library fusion arrays on the 2,581-target neighbour library
                (results/rq2_source_curve/libfam/<FAM>_p<k>.{fm,has}.npy; FENG also kolf-free, for the KOLF panel, A3)
  curve P ...   per destination panel: every subset build + offline score -> results/rq2_source_curve/sub/<P>/*.npz
  collect       -> results/rq2_source_curve.parquet (one row per panel x subset x arm) + results/rq2_source_curve.json
                (fits, AIC, LOFO marginal value with bootstrap CI, composite on H1, budget-invariance flags)
  spot_prep     rebuild the 3 chosen H1 subsets, save sums/overrides
  spot_build    H1 real-scorer packages at k = 1, ceil(K/2), K (one subset each: the median-PDS subset of the primary arm)
  spot_real     real scorer (`cell-eval2 run --preset vcc2026 --set de.backend=pdex`), one package after the other; the
                queue wraps every call in `flock results/.celleval.lock`
  (collect also folds the spot-check results into the JSON once they exist)

DESIGN. Per destination panel P the eligible families are the panel's M1 fusion families after the destination exclusions
(`results/panels.json` pack_families: own family / same dataset out, KOLF out for somatic destinations). K = their number.
For k = 1..K: all C(K, k) subsets if <= 20, else R = 20 distinct subsets drawn uniformly without replacement
(rng = default_rng([0, panel index, k])); k = K-1 is therefore exactly the leave-one-family-out (LOFO) set.

ASSUMPTIONS (stated; checks in code and in the JSON `parity` block; note knowledge/synthesis/rq2_source_curve-assumptions.md)
  A1  "M1 + M2 (no M3)" = Task 14's pre-M3 recipe (rq3_build_preds.py), rebuilt on the subset S:
      harness M1 = 0.5 x S4 + 0.5 x fused(S), 0.3 DepMap-neighbour blend from the library(S) where a covered neighbour
      exists, energy-matched to (0.5 S4 + 0.5 fused(S)) (ml_h1fuse.build); fused(S) = ml_h1fuse.fixed_fusion over the 8
      min_cells members with the frozen family weights zeroed outside S (member data = h1fuse/<P>_data.npz, the exact
      arrays the ERA pack was built from); library(S) = both recipe parts (part 0: c 100 uncentred, equal weights, x3;
      part 1: c 50 centred, W_NEW, x1.6; 0.5/0.5) fused over S only -- linear in the per-family means, so it is
      assembled exactly from per-family arrays (A4); then F2's sum side: ew (1+cpm/50)^-0.25 -> centring 0.5 ->
      guarded top-3 library-PC removal with V3 = PCs of library(S) -> M2 = rank(cross-family agreement over S)^0.5 ->
      budget match.
  A2  S4 (the learned refit + FiLM residual model, 50 % of M1) is HELD FIXED at its frozen leave-destination-out fit for
      every subset (KOLF: the A3 kolf-free S4 of Task 14). Refitting it 841 times is out of budget and its training
      table mixes sources outside the family registry. So the primary curve prices a family through the fusion part,
      the neighbour library, the library PCs and M2, CONDITIONAL on S4. Sensitivity arm `noS4`: M1 without S4
      (total = fused(S) with the same neighbour blend, energy-matched to fused(S)), i.e. a pure source-count curve.
  A3  M2 is undefined when < 2 families of S measure a target (verbatim agreement -> NaN -> median rank). If no target
      has an agreement (always at k = 1), M2 is off (constant amplitude) -- the only change to the verbatim code.
  A4  Library per family (exact): fam_mean_f(t, g) = mean over the family's sources of log((fp 1e6 + c)/(fc 1e6 + c)) - centre
      (sources with >= 20 cells for t); part L = clip(sum_f w_f fam_mean_f / sum_f w_f has_f, +-3). Stored float32.
  A5  Budget. PRIMARY arm `own`: the recipe's own budget, Bt(S) = B(max(1, kink(S))) from the subset's own M1 (kink factor
      1.17, margin 0.8; C_real(n) jackknife is control-only and computed once per panel, rng 0) -- what the builder would
      do with fewer sources. Because Bt(S) != Bt(K) in general, offline numbers of this arm are NOT guaranteed exact
      (AGENTS.md: the offline evaluator is exact iff budget-invariant) -> the `bt_rel` column flags it and the H1 real
      scorer spot check (k = 1, ceil(K/2), K) decides whether the curve is reported in absolute or relative terms
      (threshold |PDS_real - PDS_offline| > 0.002). Arm `fixBt`: same direction, budget matched to the full-set Bt(K)
      (budget-invariant relative to the K build -> offline-exact for the comparison with k = K).
  A6  Offline scorer = the ERA pack evaluator (`score_<P>.npz` real moments), with the subset's own override Od (the
      promoter-floor override depends on the fused delta). KOLF: Task 14's kolf-free FENG stores + S4 (Amendment A3).
  A7  Fits (descriptive; subsets overlap, so points are not independent and AIC is a descriptive comparison only):
      log-linear PDS = a + b ln k (OLS) vs saturating PDS = a - c exp(-k / tau) (least squares, tau in [0.05, 200]);
      Gaussian AIC = n ln(RSS / n) + 2 p (p = 3 / 4 incl. the variance), on subset-level points (primary) and on the K
      k-means (sensitivity). Delta AIC < 2 -> "indistinguishable".
  A8  LOFO marginal value = mean over families f and targets t of pds_t(K) - pds_t(K \\ f) (positive = the family helps);
      95 % percentile CI from B = 2,000 bootstrap draws resampling families (the K LOFO subsets) and targets
      independently (rng 0). Composite on H1: sum-side composite change from the vcc2026 preset scaling constants
      (tv2_score.CONST: PDS (x - 0.5)/0.4515, MSE clip((x - 0.984)/(0.043 - 0.984), 0, 1), each / 6); the DE side is
      not rebuilt (NMAE/reach/fid/jaccard fixed), so this is the PDS + MSE share only.
Seeds: subset draws default_rng([0, panel index, k]); kink jackknife rng 0 (4 draws); randomized_svd random_state 0;
bootstrap rng 0; emission rng 0 (spot check)."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import itertools, json, math, pathlib, shutil, subprocess, sys, time, zipfile
import numpy as np

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); V = ROOT / "vcc2026_data"; RP = V / "reports"; SCR = V / "scripts"
OUTD = ROOT / "output/data-analysis"; RES = OUTD / "results"; W = RES / "rq2_source_curve"
LIBF = W / "libfam"; SUB = W / "sub"; SPOT = W / "spot"
PACKS = RP / "era_2026-10-07/packs"; H1FUSE = RP / "sc_hybrid_2026-09-25/h1fuse"; HSRC = RP / "final_harness_2026-09-24/sources"
PRED = RES / "preds"; FIX = PRED / "kolf_fix"; KF_SRC = FIX / "sources_kolffree"
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
SOMATIC = {"OR_HCT116X", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"}
R_SUB = 20; B_BOOT = 2000; SPOT_THR = 0.002
ARMS = ("own", "fixBt", "noS4")
CE = pathlib.Path("/home/wonmin-kim/miniconda3/envs/vcc2026/bin/cell-eval2")
CONST_PDS = (0.5, 0.9515); CONST_MSE = (0.984, 0.043)
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.0f}s]", *a, flush=True)


def jdump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)))


def imports():
    from cell_eval2.moments import jackknife_correction  # noqa: F401  (import order: BEFORE scripts/ joins sys.path)
    from cell_eval2.prep import bulk_lognorm_means  # noqa: F401
    sys.argv = ["x"]; sys.path.append(str(SCR)); sys.path.insert(0, str(OUTD))
    sys.path.insert(0, str(ROOT / "external/era-aivc/aivc"))
    import fh_harness as H, sch_emit as EMT, tv2_model as TM
    import ml_h1fuse as F
    import common as CM
    import rq3_build_preds as R3, rq4_sameline as R4
    from tasks.vcc_sum import evaluator as EV
    return H, EMT, TM, F, CM, R3, R4, EV


def npz_member(path, key):
    """Read-only memmap of an UNCOMPRESSED .npz member (h1fuse data files are np.savez, stored)."""
    zf = zipfile.ZipFile(path); info = zf.getinfo(key + ".npy"); assert info.compress_type == 0, (path, key)
    with open(path, "rb") as f:
        f.seek(info.header_offset); h = f.read(30)
        off = info.header_offset + 30 + int.from_bytes(h[26:28], "little") + int.from_bytes(h[28:30], "little")
        f.seek(off); ver = np.lib.format.read_magic(f)
        shape, fortran, dtype = (np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0)(f)
        doff = f.tell()
    return np.memmap(path, dtype=dtype, mode="r", shape=shape, offset=doff, order="F" if fortran else "C")


def composite_sum(pds, mse):
    s_p = (pds - CONST_PDS[0]) / (CONST_PDS[1] - CONST_PDS[0])
    s_m = min(max((mse - CONST_MSE[0]) / (CONST_MSE[1] - CONST_MSE[0]), 0.0), 1.0)
    return (s_p + s_m) / 6.0


def panel_families():
    pj = json.loads((RES / "panels.json").read_text())["panels"]
    return {p: list(pj[p]["pack_families"]) for p in PANELS7}


# ======================================================================================= libfam
PART_CFG = [dict(c=100.0, center=False, alpha=3.0), dict(c=50.0, center=True, alpha=1.6)]


def cmd_libfam():
    H, EMT, TM, F, CM, R3, R4, EV = imports()
    LIBF.mkdir(parents=True, exist_ok=True)
    LIB = list(H.LIB)
    (LIBF / "lib_targets.json").write_text(json.dumps(LIB))
    jobs = [(f, HSRC, f) for f in H.FAMS] + [("FENG", KF_SRC, "FENG_kf")]
    for fam, src, tag in jobs:
        if all((LIBF / f"{tag}_p{k}.has.npy").exists() for k in range(2)):
            log(tag, "exists"); continue
        S0 = TM.SRC
        try:
            TM.SRC = pathlib.Path(src); TM._cache.clear()
            srcs = [s for s, f_ in TM.FAMILY.items() if f_ == fam and (TM.SRC / f"{s}.npz").exists()]
            assert srcs, (fam, src)
            for k, pc in enumerate(PART_CFG):
                fnum = np.zeros((len(LIB), H.NG)); fden = np.zeros((len(LIB), H.NG), np.float32)
                for s in srcs:
                    z = TM.load(s)
                    cen = np.where(np.isfinite(z["center_c50"]), z["center_c50"], 0.0) if pc["center"] else np.zeros(H.NG)
                    seen = {}                                   # predict_lfc verbatim: row = z["pos"].get(canon(LIB target))
                    for li, t in enumerate(LIB):
                        i = z["pos"].get(H.canon(t))
                        if i is None or z["ncell"][i] < 20 or H.canon(t) in TM.DROP_TARGETS:
                            continue
                        seen[li] = i
                    for li, i in seen.items():
                        v = np.log((z["fp"][i].astype(np.float64) * 1e6 + pc["c"]) / (z["fc"][i].astype(np.float64) * 1e6 + pc["c"])) - cen
                        m = np.isfinite(v); fnum[li, m] += v[m]; fden[li, m] += 1.0
                fm = np.where(fden > 0, fnum / np.maximum(fden, 1e-12), 0.0).astype(np.float32)
                np.save(LIBF / f"{tag}_p{k}.fm.npy", fm); np.save(LIBF / f"{tag}_p{k}.has.npy", fden > 0)
                log(tag, f"part {k}", "sources", srcs, "rows", int((fden > 0).any(1).sum())); del fnum, fden, fm
        finally:
            TM.SRC = S0; TM._cache.clear()
    log("LIBFAM_DONE")


# ======================================================================================= curve
def subsets_for(fams, pidx):
    K = len(fams); out = {}
    for k in range(1, K + 1):
        allc = list(itertools.combinations(fams, k))
        if len(allc) <= R_SUB:
            out[k] = allc
        else:
            rng = np.random.default_rng([0, pidx, k]); idx = sorted(rng.choice(len(allc), R_SUB, replace=False))
            out[k] = [allc[i] for i in idx]
    return out


class Panel:
    """Everything subset-independent for one destination panel, plus build(S)."""

    def __init__(self, p):
        H, EMT, TM, F, CM, R3, R4, EV = self.mods = imports()
        self.p = p; self.H, self.EMT, self.TM, self.F, self.R3, self.R4 = H, EMT, TM, F, R3, R4
        self.fams = panel_families()[p]
        dpath = H1FUSE / f"{p}_data.npz"; zs = np.load(dpath, allow_pickle=False)
        self.fin = [str(f) for f in zs["fin"]]; self.wb = zs["wb"].astype(np.float64); tg = [str(t) for t in zs["tg"]]
        assert self.fin == self.fams, (p, self.fin, self.fams)
        zi = np.load(PACKS / f"inputs_{p}.npz", allow_pickle=False)
        assert [str(t) for t in zi["targets"]] == tg and [str(f) for f in zi["families"]] == self.fin
        self.pack_fused = zi["fused"]; self.pack_knn = zi["knn"]; self.pack_has = zi["knn_has"]; self.pack_total = zi["base_total"]
        del zi
        self.tg = tg; self.T = len(tg)
        self.X = npz_member(dpath, "x"); self.M = npz_member(dpath, "m")       # (8, T, NG, F) float32 / bool, memmapped
        self.iF = {f: i for i, f in enumerate(self.fin)}
        # KOLF (A3): kolf-free FENG family column for the 8 members
        self.feng_kf = None
        if p == "KOLF":
            fk = W / "kolf_feng_members.npy"
            if not fk.exists():
                S0 = TM.SRC
                try:
                    TM.SRC = KF_SRC; TM._cache.clear()
                    arr = np.stack([F.fam_means(tg, ["FENG"], 20 if mc is None else mc)["FENG"].astype(np.float32) for mc in F.MC])
                finally:
                    TM.SRC = S0; TM._cache.clear()
                np.save(fk, arr)
            self.feng_kf = np.load(fk, mmap_mode="r")
        # cfg / promoter / target index of the fusion part (part 1 cfg, alpha 1.6) -- source-independent
        c1 = dict(H.PARTS[1][1]); self.a1 = H.PARTS[1][2]; cfg = TM.Cfg()
        for k_, v_ in c1.items(): setattr(cfg, k_, v_)
        self.cfg = cfg
        _, self.prom, self.tidx = TM.predict_lfc(tg, cfg, exclude=set(TM.FAMILY.values()) | set(TM.FAMILY))
        # S4 (fixed, A2)
        if p == "KOLF":
            self.S3 = 3.0 * np.load(FIX / "s4_KOLF_film_nopooledFENG.npy").astype(np.float64)
        else:
            key = F.S4KEY[p]; line = EMT.S4_LINE[p]
            fS = {"F": EMT.S4_EXT / line / "film_v2_a0.25.npy"}.get(key, EMT.S4_NEW / line / f"film_{key[1:]}.npy")
            self.S3 = 3.0 * np.load(fS).astype(np.float64)
        # library per family (A4)
        self.LIB = json.loads((LIBF / "lib_targets.json").read_text()); assert self.LIB == list(H.LIB)
        self.libfm, self.libhas = {}, {}
        for f in self.fin:
            tag = "FENG_kf" if (p == "KOLF" and f == "FENG") else f
            for k in range(2):
                self.libfm[(f, k)] = np.load(LIBF / f"{tag}_p{k}.fm.npy", mmap_mode="r")
                self.libhas[(f, k)] = np.load(LIBF / f"{tag}_p{k}.has.npy", mmap_mode="r")
        self.wpart = [{f: (0.0 if (p in SOMATIC and f == "KOLF") else 1.0) for f in self.fin},
                      {f: (0.0 if (p in SOMATIC and f == "KOLF") else F.W_NEW[f]) for f in self.fin}]
        self.own_lib = [(i, H.GIC.get(H.canon(t))) for i, t in enumerate(self.LIB)]
        self.own_lib = [(i, j) for i, j in self.own_lib if j is not None]
        self.lpos = {H.canon(t): i for i, t in enumerate(self.LIB)}
        self.nb = [[(self.lpos[x], c) for x, c in H.NB.get(H.canon(t), []) if x in self.lpos and x != H.canon(t) and c > 0] for t in tg]
        # destination side (controls, scorer)
        sp = np.load(PACKS / f"score_{p}.npz"); self.n = int(sp["n"]); self.gmap = sp["gmap"]; self.ok = self.gmap >= 0
        self.pack_Od = sp["Od"]; self.pack_V3 = sp["V3"].astype(np.float64); self.pack_Bt = float(sp["Bt"])
        ctrl, genes = CM.controls(p); self.genes = np.asarray(genes); self.ctrl = ctrl
        assert np.array_equal(CM.to_official(self.genes), self.gmap)
        self.fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; assert np.allclose(self.fm, sp["fm"], rtol=1e-6, atol=1e-9)
        self.cpm = R4.cpm_of(ctrl, self.ok)
        self.pk = EV.Pack(PACKS / f"score_{p}.npz")
        # pairwise family cosines per target for the verbatim M2 agreement (A3): member 0 = min_cells 20
        Fm = {}
        x0 = np.asarray(self.X[0]); m0 = np.asarray(self.M[0])
        for f in self.fin:
            i = self.iF[f]
            if self.feng_kf is not None and f == "FENG":
                v = np.asarray(self.feng_kf[0], dtype=np.float64)
            else:
                v = np.where(m0[:, :, i], x0[:, :, i].astype(np.float64), np.nan)
            Fm[f] = v
        del x0, m0
        nF = len(self.fin); self.big = np.zeros((self.T, nF), bool); self.cos = np.full((self.T, nF, nF), np.nan)
        for k in range(self.T):
            for i, f in enumerate(self.fin):
                self.big[k, i] = np.isfinite(Fm[f][k]).sum() > 500
            for i in range(nF):
                for j in range(i + 1, nF):
                    a, b = Fm[self.fin[i]][k], Fm[self.fin[j]][k]; mm = np.isfinite(a) & np.isfinite(b)
                    self.cos[k, i, j] = float(a[mm] @ b[mm] / (np.linalg.norm(a[mm]) * np.linalg.norm(b[mm]) + 1e-12))
        self.Fm_full_agree = R3.agreement(Fm, self.fin, self.T)       # verbatim, for the parity check
        del Fm
        self.cr = None

    # ------------------------------------------------------------------ pieces of build(S)
    def agreement(self, S):
        idx = [self.iF[f] for f in self.fin if f in S]; ag = np.full(self.T, np.nan)
        for k in range(self.T):
            cols = [i for i in idx if self.big[k, i]]
            cs = [self.cos[k, cols[a], cols[b]] for a in range(len(cols)) for b in range(a + 1, len(cols))]
            if cs: ag[k] = float(np.mean(cs))
        return ag

    def fused(self, S):
        wS = np.array([self.wb[i] if f in S else 0.0 for i, f in enumerate(self.fin)])
        keep = np.flatnonzero(wS > 0); wk = wS[keep]; d1s = []; ov = None
        for mi in range(self.X.shape[0]):
            xs = np.asarray(self.X[mi][:, :, keep], dtype=np.float32); ms = np.asarray(self.M[mi][:, :, keep])
            if self.feng_kf is not None and "FENG" in S:
                j = list(keep).index(self.iF["FENG"]); v = np.asarray(self.feng_kf[mi]); ms[:, :, j] = np.isfinite(v)
                xs[:, :, j] = np.nan_to_num(v)
            a = wk * ms; A = a.sum(-1)
            yh = np.where(A > 0, (a * xs).sum(-1) / np.maximum(A, 1e-300), 0.0); del a, xs, ms
            L = np.clip(yh, -self.cfg.clip, self.cfg.clip)
            dd, oo = self.TM.split(L, self.prom, self.tidx, self.a1, self.cfg); d1s.append(dd); ov = oo if ov is None else ov
            del yh, L
        return np.mean(d1s, 0), ov

    def library(self, S):
        dl = np.zeros((len(self.LIB), self.H.NG))
        for k, pc in enumerate(PART_CFG):
            num = np.zeros_like(dl); den = np.zeros_like(dl)
            for f in S:
                w = self.wpart[k][f]
                if w <= 0: continue
                num += w * self.libfm[(f, k)]; den += w * self.libhas[(f, k)]
            L = np.clip(np.where(den > 0, num / np.maximum(den, 1e-12), 0.0), -3.0, 3.0); del num, den
            dl += 0.5 * pc["alpha"] * L; del L
        for i, j in self.own_lib:
            dl[i, j] = 0.0
        return dl

    def knn_v3(self, dl):
        cov = np.abs(dl).sum(1) > 0; Kn = np.zeros((self.T, self.H.NG)); has = np.zeros(self.T, bool)
        for k_ in range(self.T):
            cand = [(i, c) for i, c in self.nb[k_] if cov[i]][:self.EMT.NB_TOPK]
            if cand:
                wv = np.array([c for _, c in cand]) ** 2; Kn[k_] = (wv[:, None] * dl[[i for i, _ in cand]]).sum(0) / wv.sum(); has[k_] = True
        from sklearn.utils.extmath import randomized_svd
        Dl = np.zeros((dl.shape[0], len(self.genes))); Dl[:, self.ok] = dl[:, self.gmap[self.ok]]; Dl = Dl[np.abs(Dl).sum(1) > 0]
        V3 = randomized_svd(Dl - Dl.mean(0), n_components=3, random_state=0)[2]
        return Kn, has, V3, int(cov.sum())

    def m1(self, d1, Kn, has, s4=True):
        total = (0.5 * self.S3 + 0.5 * d1) if s4 else d1.copy()
        e0 = float((total ** 2).sum()); w = 0.3 * has
        total = (1 - w)[:, None] * total + w[:, None] * Kn
        return total * np.sqrt(e0 / max(float((total ** 2).sum()), 1e-12))

    def budget(self, L0, Od):
        """Bt = B(max(1, kink)) exactly as rq3_build_preds.budget / sch_emit.kink_gamma (C_real cached per panel)."""
        if self.cr is None:
            _, self.cr = self.EMT.kink_gamma(self.ctrl, L0, Od, self.n, factor=1.17, margin=0.8)
        B = self.R4.Bfun(self.fm, L0, Od, self.n); claim = L0.shape[0] * 1.17 * 0.8 * self.cr
        lo, hi = 0.3, 6.0
        if B(hi) < claim:
            gk = hi
        else:
            for _ in range(40):
                mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(mid) < claim else (lo, mid)
            gk = 0.5 * (lo + hi)
        return B(max(1.0, gk)), float(gk)

    def bmatch(self, X, Od, bt):
        """Task 14 bmatch ([0.05, 8], 40 steps). DEVIATION (logged): if the scale lands at a bisection edge (a sparse
        single-family M1 can need s < 0.05), redo the bisection in log-space on [1e-4, 100]; flagged per row."""
        try:
            return self.R4.bmatch(self.fm, X, Od, self.n, bt), 0
        except AssertionError:
            B = self.R4.Bfun(self.fm, X, Od, self.n); lo, hi = math.log(1e-4), math.log(100.0)
            for _ in range(60):
                mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if B(math.exp(mid)) < bt else (lo, mid)
            s = math.exp(0.5 * (lo + hi)); assert 1.1e-4 < s < 99.0, f"budget unreachable even on [1e-4, 100]: {s}"
            return s, 1

    def build(self, S, Bt_ref=None, want_L=False):
        S = set(S); d1, ov = self.fused(S); dl = self.library(S); Kn, has, V3, ncov = self.knn_v3(dl); del dl
        ag = self.agreement(S); rk = None if np.all(np.isnan(ag)) else self.R3.rank01(ag)
        info = dict(knn_has=int(has.sum()), lib_covered=ncov, m2_defined=int(np.isfinite(ag).sum()),
                    targets_with_source=int((np.abs(d1).sum(1) > 0).sum()))
        Od = np.full((self.T, len(self.genes)), np.nan); Od[:, self.ok] = ov[:, self.gmap[self.ok]]
        self.pk.Od, self.pk.ovj = Od, np.isfinite(Od)
        out, Ls = {}, {}
        for arm_s4, s4 in (("frozen", True), ("noS4", False)):
            tot = self.m1(d1, Kn, has, s4=s4); L0 = self.R4.on_panel(tot, self.gmap)
            if not s4 and not np.any(L0):
                # the subset's families cover NO panel target: fused(S) = 0, the noS4 energy reference is 0 and the recipe
                # returns the zero delta (overrides only); scored as such, flagged (scale / budget undefined)
                e = self.pk.evaluate(np.zeros_like(L0))
                out["noS4"] = dict(pds=e["pds"], mse=e["mse"], scale=float("nan"), bmatch_edge=2.0, Bt=float("nan"), kink=float("nan"),
                                   Bt_own=float("nan"), pds_t=e["pds_t"], num=e["num"], dist=e["dist"], pcr_energy=float("nan"),
                                   pcr_lambda=float("nan"))
                if want_L: Ls["noS4"] = np.zeros_like(L0)
                del L0; continue
            Bt, gk = self.budget(L0, Od)
            X, tinfo = self.R4.sum_side(L0, Od, V3, self.cpm, rk)
            budgets = [("own", Bt)] if not s4 else [("own", Bt)] + ([("fixBt", Bt_ref)] if Bt_ref is not None else [])
            for nm, bt in budgets:
                arm = nm if s4 else "noS4"
                s, edge = self.bmatch(X, Od, bt); L = X * s; e = self.pk.evaluate(L)
                out[arm] = dict(pds=e["pds"], mse=e["mse"], scale=s, bmatch_edge=float(edge), Bt=bt, kink=gk, Bt_own=Bt, pds_t=e["pds_t"], num=e["num"],
                                dist=e["dist"], **tinfo)
                if want_L: Ls[arm] = L
            if s4: tot_frozen, V3_, Bt_frozen = tot, V3, Bt
            del L0, X
        return out, info, dict(total=tot_frozen, d1=d1, Kn=Kn, has=has, V3=V3_, Od=Od, Bt=Bt_frozen, rk=rk, ag=ag, Ls=Ls)


def cmd_curve(panels):
    for p in panels:
        pidx = PANELS7.index(p); dd = SUB / p; dd.mkdir(parents=True, exist_ok=True)
        tp = time.time(); P = Panel(p); fams = P.fams; K = len(fams); subs = subsets_for(fams, pidx)
        log("=====", p, "K", K, "families", fams, "subsets", {k: len(v) for k, v in subs.items()}, "setup secs", round(time.time() - tp, 1))
        # full set first: parity + reference budget
        full = tuple(fams); out, info, aux = P.build(full)
        Bt_ref = aux["Bt"]; par = {}
        if p != "KOLF":
            par["fused_maxabs"] = float(np.abs(aux["d1"] - P.pack_fused).max())
            par["knn_maxabs"] = float(np.abs(aux["Kn"] - P.pack_knn).max()); par["knn_has_equal"] = bool(np.array_equal(aux["has"], P.pack_has))
            par["base_total_maxabs"] = float(np.abs(aux["total"] - P.pack_total).max())
            par["V3_absdot_min"] = float(np.abs((aux["V3"] * P.pack_V3).sum(1)).min())
            par["Od_equal"] = bool(np.array_equal(np.isfinite(aux["Od"]), np.isfinite(P.pack_Od)) and
                                   np.allclose(np.nan_to_num(aux["Od"]), np.nan_to_num(P.pack_Od), atol=1e-6))
            par["Bt_rel"] = float(abs(Bt_ref - P.pack_Bt) / P.pack_Bt)
        zt = np.load(PRED / f"{p}_prem3.npz", allow_pickle=True)
        par["m2_rank_equal_task14"] = bool(np.allclose(aux["rk"], zt["m2_rank"]))
        par["agree_equal_verbatim"] = bool(np.allclose(np.nan_to_num(aux["ag"], nan=9), np.nan_to_num(P.Fm_full_agree, nan=9)))
        out2, _, aux2 = P.build(full, Bt_ref=Bt_ref, want_L=True)
        Lk = aux2["Ls"]["own"]; par["prem3_lfc_maxrel_task14"] = float(np.abs(Lk - zt["lfc"]).max() / np.abs(zt["lfc"]).max())
        par["prem3_pds_task14"] = float(zt["pds_offline"]); par["prem3_pds_here"] = out2["own"]["pds"]
        par["prem3_mse_task14"] = float(zt["mse_offline"]); par["prem3_mse_here"] = out2["own"]["mse"]
        del Lk, aux2, zt
        log(p, "PARITY", par)
        ok = (par["prem3_lfc_maxrel_task14"] < 1e-4 and abs(par["prem3_pds_task14"] - par["prem3_pds_here"]) < 1e-6)
        if p != "KOLF":
            ok = ok and par["base_total_maxabs"] < 1e-4 and par["V3_absdot_min"] > 0.999 and par["Bt_rel"] < 1e-5 and par["Od_equal"]
        par["pass"] = bool(ok); jdump(dd / "parity.json", dict(panel=p, K=K, families=fams, Bt_ref=Bt_ref, parity=par,
                                                                subsets={k: ["+".join(s) for s in v] for k, v in subs.items()}))
        assert ok, par
        order = [K, K - 1, 1] + [k for k in range(2, K - 1)]
        order = list(dict.fromkeys(k for k in order if k >= 1))
        for k in order:
            for si, S in enumerate(subs[k]):
                f = dd / f"k{k:02d}_s{si:02d}.npz"
                if f.exists(): continue
                ts = time.time(); o, info, _ = P.build(S, Bt_ref=Bt_ref)
                row = dict(panel=p, k=k, sub_idx=si, subset="+".join(S), K=K, **info)
                arrs = {}
                for arm, r in o.items():
                    for key in ("pds", "mse", "scale", "bmatch_edge", "Bt", "kink", "Bt_own", "pcr_energy", "pcr_lambda"):
                        row[f"{arm}_{key}"] = float(r[key])
                    arrs[f"{arm}_pds_t"] = r["pds_t"].astype(np.float32); arrs[f"{arm}_num"] = r["num"]; arrs[f"{arm}_dist"] = r["dist"]
                row["bt_rel"] = row["own_Bt"] / Bt_ref - 1.0; row["secs"] = round(time.time() - ts, 1)
                np.savez(f, row=np.array(json.dumps(row)), **arrs)
                log(p, f"k={k} s{si:02d} {row['subset']}: PDS own {row['own_pds']:.4f} fixBt {row['fixBt_pds']:.4f} noS4 {row['noS4_pds']:.4f}"
                       f" | MSE own {row['own_mse']:.4f} | bt_rel {row['bt_rel']:+.4f} | {row['secs']}s")
        log(p, "PANEL_DONE secs", round(time.time() - tp, 1)); del P
        import gc; gc.collect()
    log("CURVE_DONE")


# ======================================================================================= collect / fits
def fit_curves(k, y):
    from scipy.optimize import curve_fit
    k = np.asarray(k, float); y = np.asarray(y, float); n = len(y); res = {}
    A = np.c_[np.ones(n), np.log(k)]; coef, *_ = np.linalg.lstsq(A, y, rcond=None); rss = float(((y - A @ coef) ** 2).sum())
    res["loglin"] = dict(a=float(coef[0]), b=float(coef[1]), rss=rss, aic=n * math.log(max(rss, 1e-300) / n) + 2 * 3)
    f = lambda kk, a, c, tau: a - c * np.exp(-kk / tau)
    best = None
    for tau0 in (0.5, 1, 2, 4, 8, 16, 50):
        try:
            p0 = (float(y.max()), float(y.max() - y.min()) * math.exp(1 / tau0), tau0)
            pp, _ = curve_fit(f, k, y, p0=p0, bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 200.0]), maxfev=20000)
            r = float(((y - f(k, *pp)) ** 2).sum())
            if best is None or r < best[1]: best = (pp, r)
        except Exception:  # noqa
            pass
    if best is not None:
        pp, r = best
        res["saturating"] = dict(a=float(pp[0]), c=float(pp[1]), tau=float(pp[2]), rss=r, aic=n * math.log(max(r, 1e-300) / n) + 2 * 4,
                                 tau_at_bound=bool(pp[2] > 199.0 or pp[2] < 0.051),
                                 degenerate_step=bool(pp[2] < 0.5))      # tau < 0.5: the "saturation" is a k=1 -> k>=2 step
        d = res["saturating"]["aic"] - res["loglin"]["aic"]
        res["dAIC_sat_minus_loglin"] = float(d)
        res["verdict"] = "indistinguishable" if abs(d) < 2 else ("saturating" if d < 0 else "log-linear")
    res["n"] = n
    return res


def lofo(rows_full, rows_lofo, arrs, panel, arm="own", rng=None):
    """A8: Delta_{f,t} = pds_t(K) - pds_t(K \\ f); composite share uses MSE = sum num / sum dist over resampled targets."""
    pf = arrs[(panel, rows_full["k"], rows_full["sub_idx"])]; P_K = pf[f"{arm}_pds_t"].astype(np.float64)
    numK, distK = pf[f"{arm}_num"], pf[f"{arm}_dist"]
    fams = rows_full["subset"].split("+"); D, NUM, DIST, drop = [], [], [], []
    for r in rows_lofo:
        a = arrs[(panel, r["k"], r["sub_idx"])]; miss = sorted(set(fams) - set(r["subset"].split("+")))
        assert len(miss) == 1; drop.append(miss[0])
        D.append(P_K - a[f"{arm}_pds_t"].astype(np.float64)); NUM.append(a[f"{arm}_num"]); DIST.append(a[f"{arm}_dist"])
    D = np.array(D); NUM = np.array(NUM); DIST = np.array(DIST); nF, T = D.shape
    rng = rng or np.random.default_rng(0)
    comp = lambda ip, im: composite_sum(ip, im)
    cK = comp(float(P_K.mean()), float(numK.sum() / distK.sum()))
    c_f = np.array([cK - comp(float(P_K.mean() - D[i].mean()), float(NUM[i].sum() / DIST[i].sum())) for i in range(nF)])
    bs, bc = np.empty(B_BOOT), np.empty(B_BOOT)
    for b in range(B_BOOT):
        fi = rng.integers(0, nF, nF); ti = rng.integers(0, T, T)
        bs[b] = D[np.ix_(fi, ti)].mean()
        pK = P_K[ti].mean(); mK = numK[ti].sum() / distK[ti].sum(); cKb = comp(pK, mK)
        bc[b] = np.mean([cKb - comp(pK - D[i, ti].mean(), NUM[i, ti].sum() / DIST[i, ti].sum()) for i in fi])
    per = {}
    for i, f in enumerate(drop):
        bt = np.array([D[i, rng.integers(0, T, T)].mean() for _ in range(B_BOOT)])
        per[f] = dict(dPDS=float(D[i].mean()), ci95=[float(np.quantile(bt, 0.025)), float(np.quantile(bt, 0.975))],
                      dcomposite_sum=float(c_f[i]))
    return dict(mean_dPDS=float(D.mean()), ci95=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))],
                mean_dcomposite_sum=float(c_f.mean()), ci95_composite_sum=[float(np.quantile(bc, 0.025)), float(np.quantile(bc, 0.975))],
                mean_dcomposite_pds_only=float(D.mean() / (6 * (CONST_PDS[1] - CONST_PDS[0]))),
                per_family=per, n_families=nF, n_targets=T)


def cmd_collect():
    import pandas as pd
    rows, arrs = [], {}
    for p in PANELS7:
        for f in sorted((SUB / p).glob("k*_s*.npz")):
            z = np.load(f, allow_pickle=False); r = json.loads(str(z["row"])); rows.append(r)
            arrs[(p, r["k"], r["sub_idx"])] = {k: z[k] for k in z.files if k != "row"}
    df = pd.DataFrame(rows)
    long = []
    for r in rows:
        for arm in ARMS:
            if f"{arm}_pds" not in r: continue
            long.append(dict(panel=r["panel"], K=r["K"], k=r["k"], sub_idx=r["sub_idx"], subset=r["subset"], arm=arm,
                             pds=r[f"{arm}_pds"], mse=r[f"{arm}_mse"], budget_scale=r[f"{arm}_scale"],
                             bmatch_wide=bool(r.get(f"{arm}_bmatch_edge", 0) == 1), zero_delta=bool(r.get(f"{arm}_bmatch_edge", 0) == 2), Bt=r[f"{arm}_Bt"],
                             Bt_own=r[f"{arm}_Bt_own"], kink=r[f"{arm}_kink"], pcr_lambda=r[f"{arm}_pcr_lambda"],
                             bt_rel_own_vs_K=r["bt_rel"] if arm != "noS4" else None, knn_has=r["knn_has"], m2_defined=r["m2_defined"],
                             targets_with_source=r["targets_with_source"], composite_sum=composite_sum(r[f"{arm}_pds"], r[f"{arm}_mse"]),
                             secs=r["secs"]))
    L = pd.DataFrame(long); L.to_parquet(RES / "rq2_source_curve.parquet", index=False)
    out = dict(plan="input/ideas/zero-shot-limits-context-complexity-plan.md", task=11, hypothesis=None,
               note="descriptive (no hypothesis attached)", recipe=__doc__.split("ASSUMPTIONS")[1].split("Seeds:")[0].strip(),
               arms=dict(own="PRIMARY: frozen M1 (S4 fixed) + M2, own budget (recipe-faithful)",
                         fixBt="same direction, budget matched to the full-set Bt(K) (offline-exact relative to k = K)",
                         noS4="sensitivity: M1 without S4 (pure source-count curve), own budget"),
               R=R_SUB, B=B_BOOT, panels={})
    rng = np.random.default_rng(0)
    for p in PANELS7:
        d = L[L.panel == p]
        if d.empty: continue
        pj = json.loads((SUB / p / "parity.json").read_text()); K = pj["K"]; po = dict(K=K, families=pj["families"], parity=pj["parity"])
        po["n_subsets_done"] = int(d[d.arm == "own"].shape[0]); po["complete"] = bool(po["n_subsets_done"] == sum(len(v) for v in pj["subsets"].values()))
        for arm in ARMS:
            a = d[d.arm == arm]
            if a.empty: continue
            km = a.groupby("k").agg(pds=("pds", "mean"), pds_sd=("pds", "std"), mse=("mse", "mean"), n=("pds", "size"),
                                    composite_sum=("composite_sum", "mean")).reset_index()
            ao = dict(curve={int(r.k): dict(pds=r.pds, pds_sd=(None if np.isnan(r.pds_sd) else r.pds_sd), mse=r.mse, n=int(r.n),
                                            composite_sum=r.composite_sum) for r in km.itertuples()})
            if a.k.nunique() >= 3:
                ao["fit_subsets"] = fit_curves(a.k.values, a.pds.values)
                ao["fit_kmeans"] = fit_curves(km.k.values, km.pds.values)
                A_ = np.c_[np.ones(len(a)), np.log(a.k.values)]; ao["mse_loglin_slope"] = float(np.linalg.lstsq(A_, a.mse.values, rcond=None)[0][1])
            if K in set(a.k) and (K - 1) in set(a.k) and (a.k == K - 1).sum() == K:
                rf = df[(df.panel == p) & (df.k == K)].iloc[0].to_dict(); rl = df[(df.panel == p) & (df.k == K - 1)].to_dict("records")
                ao["lofo"] = lofo(rf, rl, arrs, p, arm=arm, rng=rng)
                ao["last_family_marginal_kmeans"] = float(km[km.k == K].pds.iloc[0] - km[km.k == K - 1].pds.iloc[0])
            po[arm] = ao
        if "bt_rel_own_vs_K" in d:
            b = d[d.arm == "own"].bt_rel_own_vs_K.astype(float)
            po["budget"] = dict(max_abs_bt_rel=float(b.abs().max()), median_abs_bt_rel=float(b.abs().median()),
                                budget_changes=bool(b.abs().max() > 1e-5),
                                own_minus_fixBt_pds_maxabs=float((d[d.arm == "own"].pds.values - d[d.arm == "fixBt"].pds.values).__abs__().max())
                                if (d.arm == "fixBt").sum() == (d.arm == "own").sum() else None)
        out["panels"][p] = po
    # spot check
    spot = {}
    for f in sorted(SPOT.glob("*/result.json")):
        r = json.loads(f.read_text()); spot[r["tag"]] = r
    if spot:
        dev = {t: r["raw"]["pds_cosine"] - r["meta"]["pds_offline"] for t, r in spot.items()}
        out["spot_check"] = dict(panel="H1", threshold=SPOT_THR, runs={t: dict(k=r["meta"]["k"], subset=r["meta"]["subset"],
                                 pds_offline=r["meta"]["pds_offline"], pds_real=r["raw"]["pds_cosine"], dPDS=dev[t],
                                 mse_offline=r["meta"]["mse_offline"], mse_real=r["raw"]["expr_mse_unbiased_capped_norm"],
                                 composite_real=r["composite"]) for t, r in spot.items()},
                                 max_abs_dPDS=float(max(abs(v) for v in dev.values())), n_runs=len(spot), complete=len(spot) == 3)
        out["spot_check"]["curves_reported_as"] = ("relative" if out["spot_check"]["max_abs_dPDS"] > SPOT_THR else
                                                   ("absolute" if len(spot) == 3 else "pending (incomplete spot check)"))
    else:
        out["spot_check"] = dict(status="queued (real scorer, unit aivc-paper-t11-spot)", curves_reported_as="pending")
    # k = K proxy (available before the queued spot runs): Task 18 arm noM3 = Task 14 pre-M3 sums on H1 = this script's
    # k = K build (parity: lfc max rel < 1e-7, same Od, same budget), scored with the real scorer (DE side = F2 mirror)
    t18 = RES / "rq4/pkg/noM3/result.json"
    if t18.exists():
        r = json.loads(t18.read_text()); par = out["panels"].get("H1", {}).get("parity", {})
        out["spot_check"]["k_eq_K_proxy_task18_noM3"] = dict(pds_offline=r["meta"]["pds_offline"], pds_real=r["raw"]["pds_cosine"],
                                                             dPDS=r["raw"]["pds_cosine"] - r["meta"]["pds_offline"],
                                                             mse_offline=r["meta"]["mse_offline"], mse_real=r["raw"]["expr_mse_unbiased_capped_norm"],
                                                             sums_identity_maxrel=par.get("prem3_lfc_maxrel_task14"))
    try:
        plot_curves(L, out)
    except Exception as e:  # noqa
        log("plot failed", repr(e))
    h = out["panels"].get("H1", {}).get("own", {}).get("lofo")
    out["per_context_value_H1"] = (dict(source="Task 11 LOFO, primary arm (own budget), H1",
                                        dPDS=h["mean_dPDS"], dPDS_ci95=h["ci95"], dcomposite_sum=h["mean_dcomposite_sum"],
                                        dcomposite_sum_ci95=h["ci95_composite_sum"], dcomposite_pds_only=h["mean_dcomposite_pds_only"],
                                        note="sum-side composite (PDS + MSE legs, vcc2026 preset constants); DE side not rebuilt; "
                                             "rq4_sameline.py (Task 18) hard-codes per_context_value=None and rewrites rq4_h4.json on "
                                             "every verdict call, so the value is kept HERE and not written into rq4_h4.json")
                                   if h else None)
    jdump(RES / "rq2_source_curve.json", out); log("COLLECT_DONE", len(L), "rows")


def plot_curves(L, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(16, 7.5), sharex=False); axs = axs.ravel()
    for ax, p in zip(axs, PANELS7):
        for arm, col in (("own", "#1f5fa8"), ("noS4", "#c0582b")):
            a = L[(L.panel == p) & (L.arm == arm)]
            if a.empty: continue
            ax.scatter(a.k, a.pds, s=8, alpha=0.35, color=col)
            km = a.groupby("k").pds.mean(); ax.plot(km.index, km.values, "-o", color=col, ms=3, label=f"{arm} (mean)")
            fl = out["panels"][p][arm].get("fit_subsets", {}).get("loglin")
            if fl:
                kk = np.linspace(1, a.k.max(), 50); ax.plot(kk, fl["a"] + fl["b"] * np.log(kk), ":", color=col, lw=1)
        ax.set_title(f"{p} (K={out['panels'][p]['K']})", fontsize=10); ax.set_xlabel("source families k"); ax.set_ylabel("PDS (offline)")
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8); axs[-1].axis("off")
    axs[-1].text(0.0, 0.5, "own = frozen M1 (S4 fixed) + M2, own budget\nnoS4 = M1 without S4\ndots = subsets; dotted = a + b ln k",
                 fontsize=9, va="center")
    fig.tight_layout(); fig.savefig(RES / "rq2_source_curve.png", dpi=130); plt.close(fig)


# ======================================================================================= spot check (H1, real scorer)
def spot_choice():
    pj = json.loads((SUB / "H1" / "parity.json").read_text()); K = pj["K"]; ks = [1, math.ceil(K / 2), K]; ch = {}
    for k in ks:
        rs = []
        for f in sorted((SUB / "H1").glob(f"k{k:02d}_s*.npz")):
            r = json.loads(str(np.load(f)["row"])); rs.append(r)
        assert rs, f"H1 k={k} not computed yet"
        rs.sort(key=lambda r: (r["own_pds"], r["sub_idx"])); ch[k] = rs[(len(rs) - 1) // 2]     # median (lower) PDS subset
    return ch


def cmd_spot_prep():
    """Rebuild the three chosen H1 subsets (deterministic) and save their sums / overrides; no emission (memory)."""
    ch = spot_choice(); P = Panel("H1"); pj = json.loads((SUB / "H1" / "parity.json").read_text()); Bt_ref = pj["Bt_ref"]
    for k, row in ch.items():
        tag = f"k{k}"; out = SPOT / tag; out.mkdir(parents=True, exist_ok=True)
        if (out / "inputs.npz").exists(): log(tag, "inputs exist"); continue
        S = row["subset"].split("+"); o_, info, aux = P.build(S, Bt_ref=Bt_ref, want_L=True)
        assert abs(o_["own"]["pds"] - row["own_pds"]) < 1e-9, "rebuild != stored subset result"
        np.savez_compressed(out / "inputs.npz", Lsum=aux["Ls"]["own"], Od=aux["Od"])
        meta = dict(tag=tag, k=k, subset=row["subset"], sub_idx=row["sub_idx"], pds_offline=o_["own"]["pds"], mse_offline=o_["own"]["mse"],
                    budget_scale=o_["own"]["scale"], bt_rel=o_["own"]["Bt"] / Bt_ref - 1.0, choice="median offline PDS (primary arm) at this k",
                    de_side="Task 18 de_delta_H1.npz (full-set denoised DE delta, drc 0.8/1.1) -- NOT rebuilt per subset: only PDS/MSE "
                            "(pinned pseudobulk sums) are compared; DE metrics of these runs are not interpretable")
        jdump(out / "meta_prep.json", meta); log(tag, meta)
    log("SPOT_PREP_DONE")


def cmd_spot_build():
    import anndata as ad, pandas as pd, gc
    from scipy import sparse
    H, EMT, TM, F, CM, R3, R4, EV = imports(); SE = EMT
    zd = np.load(RES / "rq4/de_delta_H1.npz"); d_de = zd["d_de"].astype(np.float64)          # Task 18 DE side (full set), fixed
    tg = [str(t) for t in zd["targets"]]; assert tg == [str(t) for t in H.PANELS["H1"]["targets"]]
    sp = np.load(PACKS / "score_H1.npz"); gmap = sp["gmap"]; ok = gmap >= 0; n = int(sp["n"])
    r = H.real("H1"); r.pop("R", None); ctrl = r["ctrl"]; genes = np.asarray(r["genes"]).astype(str)
    fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; assert np.allclose(fm, sp["fm"], rtol=1e-6, atol=1e-9)
    fl = H.fam_lfc(tg, H.EXCL["H1"]); g_ = gmap[ok]
    S_ = np.stack([np.nan_to_num(v[:, g_]) for v in fl.values()]); M_ = np.stack([np.isfinite(v[:, g_]) for v in fl.values()])
    mean_ = np.where(M_.sum(0) > 0, S_.sum(0) / np.maximum(M_.sum(0), 1), 0.0)
    fr_ = np.where(M_.sum(0) > 0, (M_ & (np.sign(S_) == np.sign(mean_)[None])).sum(0) / np.maximum(M_.sum(0), 1), 0.5)
    AG = np.zeros((len(tg), len(genes))); AG[:, ok] = np.clip(2 * fr_ - 1, 0, 1); del S_, M_, fl; gc.collect()
    o = SE.parse(R4.EMIT_ARM); gate = r["cpm"] >= o.get("gt", 5.0); low = ~gate
    for out in sorted(SPOT.glob("k*")):
        tag = out.name
        if (out / "pred.h5ad").exists() and (out / "meta.json").exists() or (out / "result.json").exists():
            log(tag, "package exists / scored"); continue
        zi = np.load(out / "inputs.npz"); Lsum = zi["Lsum"]; Od = zi["Od"]; meta = json.loads((out / "meta_prep.json").read_text())
        LdDE = R4.drc(R4.on_panel(d_de, gmap), Od, *R4.DRC); rs = np.abs(LdDE) * AG
        # 2026-11-06 execution fix (memory only, content unchanged): the first run materialised all emitted blocks AND their
        # vstack copy (~2 x 5 GB) and was OOM-killed at the 8 GB scope cap. The SAME matrix (blocks in target order, then the
        # float32 controls; same obs/var) is now streamed block by block into the h5ad, and the written file is verified
        # against running hashes of the emitted data / indices / row lengths. Pre-fix copy: rq2_source_curve.pre_stream_2026-11-06.py
        import h5py, hashlib
        from anndata.io import write_elem, sparse_dataset
        nrow = len(tg) * n + ctrl.shape[0]; chk = {0, len(tg) // 2, len(tg) - 1}
        tmp = out / "pred.h5ad.part"; tmp.unlink(missing_ok=True)
        ad.AnnData(obs=pd.DataFrame({"target": pd.Categorical([*np.repeat(tg, n), *["non-targeting"] * ctrl.shape[0]])},
                                    index=[f"p{i}" for i in range(nrow)]), var=pd.DataFrame(index=list(genes))).write_h5ad(tmp)
        hd, hi, hr = hashlib.sha256(), hashlib.sha256(), hashlib.sha256(); nnz = 0; nt = 0
        def put(b, ds):
            nonlocal_nnz = b.nnz
            hd.update(np.ascontiguousarray(b.data, np.float32).tobytes()); hi.update(np.ascontiguousarray(b.indices, np.int64).tobytes())
            hr.update(np.diff(b.indptr).astype(np.int64).tobytes())
            if ds is None:
                b0 = b.copy(); b0.indptr = b0.indptr.astype(np.int64); write_elem(f, "X", b0); return sparse_dataset(f["X"]), nonlocal_nnz
            ds.append(b); return ds, nonlocal_nnz
        with h5py.File(tmp, "a") as f:
            ds = None
            for t, b in enumerate(SE.emit(ctrl, Lsum, Od, n, o, gate, low, rs, "H1", LdDE)):
                assert b.format == "csr" and b.dtype == np.float32 and b.shape == (n, len(genes)), f"{tag}: block {t} {b.format} {b.dtype} {b.shape}"
                if t in chk:
                    e = np.exp(Lsum[t]); j = np.isfinite(Od[t]); e[j] = np.exp(Od[t, j])
                    assert np.array_equal(np.asarray(b.sum(0)).ravel(), np.round(fm * e * n)), f"{tag}: pinned sums broken"
                ds, k_ = put(b, ds); nnz += k_; nt += 1; del b
            assert nt == len(tg), f"{tag}: emitted {nt} blocks for {len(tg)} targets"
            del LdDE, rs; gc.collect()
            ds, k_ = put(ctrl.astype(np.float32).tocsr(), ds); nnz += k_
        with h5py.File(tmp, "r") as f:                       # read-back verification of the written matrix
            g = f["X"]; assert g.attrs["encoding-type"] == "csr_matrix" and tuple(g.attrs["shape"]) == (nrow, len(genes))
            ip = g["indptr"][:]; assert int(ip[-1]) == nnz and len(ip) == nrow + 1
            assert hashlib.sha256(np.diff(ip).astype(np.int64).tobytes()).hexdigest() == hr.hexdigest(), f"{tag}: row lengths differ"
            rd, ri = hashlib.sha256(), hashlib.sha256(); CH = 50_000_000
            for s0 in range(0, nnz, CH):
                rd.update(np.ascontiguousarray(g["data"][s0:s0 + CH], np.float32).tobytes())
                ri.update(np.ascontiguousarray(g["indices"][s0:s0 + CH], np.int64).tobytes())
            assert rd.hexdigest() == hd.hexdigest() and ri.hexdigest() == hi.hexdigest(), f"{tag}: written X differs from emitted blocks"
        tmp.rename(out / "pred.h5ad"); gc.collect()
        log(tag, "STREAMED", dict(rows=nrow, nnz=nnz, data_sha256=hd.hexdigest()))
        meta.update(emit_arm=R4.EMIT_ARM, emit_seed=o["seed"]); jdump(out / "meta.json", meta); log(tag, "PACKAGE", meta)
    log("SPOT_BUILD_DONE")


def cmd_spot_real():
    import pandas as pd
    from cell_eval2.moments import jackknife_correction  # noqa: F401
    sys.path.append(str(SCR)); import tv2_score as S
    truth = RP / "transfer_v2_2026-09-23/h1/real_thin.h5ad"; cache_real = RES / "rq3_cache/cache_real_H1"
    for d in sorted(SPOT.glob("k*")):
        res_f = d / "result.json"
        if res_f.exists(): log(d.name, "already scored"); continue
        ta = time.time()
        cmd = [str(CE), "run", "-ap", str(d / "pred.h5ad"), "-ar", str(truth), "--preset", "vcc2026", "--pert-col", "target",
               "--control", "non-targeting", "--input-type", "counts", "-o", str(d / "cell_eval2"), "--cache-real", str(cache_real),
               "--set", "de.backend=pdex", "--set", "num_threads=4", "--set", "gather_threads=4"]
        env = dict(os.environ, NUMBA_NUM_THREADS="4", OMP_NUM_THREADS="4")
        with (d / "scorer.log").open("w") as g:
            subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
        agg = pd.read_csv(d / "cell_eval2/agg_results.csv", index_col=0)
        raw = {k: float(agg.loc["mean", k]) for k in S.CONST}; comp, scaled = S.composite(raw)
        shutil.copy(d / "cell_eval2/agg_results.csv", d / "agg_results.csv")
        meta = json.loads((d / "meta.json").read_text())
        jdump(res_f, dict(tag=d.name, raw=raw, scaled=scaled, composite=comp, seconds=round(time.time() - ta, 1), meta=meta, scorer=" ".join(cmd[1:])))
        shutil.rmtree(d / "cell_eval2", ignore_errors=True); (d / "pred.h5ad").unlink(missing_ok=True)
        log(d.name, "REAL pds", raw["pds_cosine"], "offline", meta["pds_offline"], "d", raw["pds_cosine"] - meta["pds_offline"])
        return                                            # one package per call (the queue takes the lock per call)
    log("SPOT_REAL_DONE")


if __name__ == "__main__":
    cmd = sys.argv[1]; args = sys.argv[2:]
    if cmd == "libfam":
        cmd_libfam()
    elif cmd == "curve":
        cmd_curve(args or PANELS7)
    elif cmd == "collect":
        cmd_collect()
    elif cmd == "spot_prep":
        cmd_spot_prep()
    elif cmd == "spot_build":
        cmd_spot_build()
    elif cmd == "spot_real":
        cmd_spot_real()
    else:
        raise SystemExit(f"unknown subcommand {cmd}")
