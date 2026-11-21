"""Build a transfer_v2 official submission h5ad (A/B/C x 300 targets x 400 cells, raw integer counts).

SPEC (json):
  {"mode": "shared_cellshift_gated" | "shared_mult" | "iid_mult" | "shared_cellshift" | "dual_mult",
   "dual": {"gamma": 1.0, "K": 3000, "depth_up": 0.5, "agree": npz, "cpc": false},  # dual_mult only (cpc: DE-functional pin)
   "seed": int,
   "parts": [ {"w": 1.0, "alpha": 3.0, "cfg": {...tv2_model.Cfg overrides...}},          # fusion part
              {"w": 1.0, "alpha": 3.0, "delta": {"A": path, "B": path, "C": path}, "promoter": {...cfg}} ],  # file part
   "cfg_ctx": {"A": {...overrides applied to every fusion part for that context}},
   "blend": {"knn": npz from tv2_blend_prep.py, "w": 0.2, "energy": true}}   # optional DepMap-neighbour blend
The per-context delta = sum_i w_i * alpha_i * part_i, then (for fusion parts) the promoter floor.
Emission: tv2_emit (one draw of 400 of that context's own control cells for shared_* modes; genes under
5 CPM in the context's controls are always multiplicative in *_gated mode, so no counts are invented for
genes the context does not express). Streams rows to an on-disk CSR.
Usage: tv2_submit_build.py SPEC_JSON OUT_H5AD
"""
import json, pathlib, sys, time
import h5py, numpy as np, pandas as pd
from scipy import sparse
from cell_eval2.moments import jackknife_correction  # noqa: F401  (before scripts/ joins sys.path)
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_model as TM, tv2_emit as EM

OFF = TM.OFF
import os as _os
CELLS = int(_os.environ.get("TV2_CELLS", 400))   # 400 = the official task; holdout rehearsals set it to the
                                                 # panel's own cells/target so the scorer's n_pred matches n_real
MAX_COUNTS_PER_CELL = 1_000_000
MAX_STORED = 4.75e9


def load_controls(ctx, base=None):
    with h5py.File((pathlib.Path(base) if base else OFF) / f"context_{ctx}.h5ad", "r") as h:
        g = h["X"]
        M = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(int(v) for v in g.attrs["shape"]))
    return M


def context_delta(spec, ctx, targets):
    total = np.zeros((len(targets), TM.NG)); info = []
    fusion_i = -1
    for part in spec["parts"]:
        if "delta" in part:
            d = np.load(part["delta"][ctx]).astype(np.float64)
            assert d.shape == (len(targets), TM.NG)
            total += part["w"] * part["alpha"] * d
            info.append(dict(kind="file", path=part["delta"][ctx], w=part["w"], alpha=part["alpha"]))
        else:
            fusion_i += 1
            cfg = TM.Cfg()
            per_part = spec.get("cfg_ctx_parts", {}).get(ctx, [])
            extra = per_part[fusion_i] if fusion_i < len(per_part) else {}
            merged = {**part.get("cfg", {}), **spec.get("cfg_ctx", {}).get(ctx, {}), **extra}
            for k, v in merged.items():
                setattr(cfg, k, v)
            L, prom, tidx, covered = TM.predict_lfc(targets, cfg, return_parts=True)
            if part.get("cfg_ens"):
                # HYPERPARAMETER ENSEMBLE: this part's L becomes the MEAN over the shipped cfg plus each override in
                # `cfg_ens`. Measured with min_cells 5/10/15/30/50/100/200: PDS and MSE improve TOGETHER on H1 / KOLF /
                # Orion, which no re-weighting of a single prediction does (FINDINGS sections 16.4 and 17).
                Ls = [L]
                for ov_ in part["cfg_ens"]:
                    c2 = TM.Cfg()
                    for k, v in {**merged, **ov_}.items():
                        setattr(c2, k, v)
                    Ls.append(TM.predict_lfc(targets, c2)[0])
                L = np.mean(Ls, axis=0)
                info.append(dict(kind="cfg_ens", members=len(Ls), overrides=part["cfg_ens"]))
            total += part["w"] * part["alpha"] * L
            info.append(dict(kind="fusion", w=part["w"], alpha=part["alpha"], covered=int(covered.sum()),
                             cfg=dict(c=cfg.c, center=cfg.center, center_frac=cfg.center_frac, weights=cfg.weights, promoter=cfg.promoter)))
    bl = spec.get("blend")
    if bl:   # DepMap-neighbour blend (fh_harness blend20_e): (1-w) * delta + w * neighbour estimate, same total energy
        z = np.load(bl["knn"]); pos = {str(t): i for i, t in enumerate(z["targets"])}
        K = np.zeros_like(total); has = np.zeros(len(targets), bool); unc = np.zeros(len(targets), bool)
        for k, t in enumerate(targets):
            i = pos.get(t)
            if i is not None and z["has"][i]:
                K[k] = z["knn"][i]; has[k] = True
                if "covered" in z.files and not z["covered"][i]:
                    unc[k] = True
        e0 = float((total ** 2).sum()); w = bl.get("w", 0.2) * has
        if bl.get("uncovered_full", True):       # a target no screen measured: neighbour estimate only (fh unseen_knn)
            w[unc] = 1.0
        total = (1 - w)[:, None] * total + w[:, None] * K
        if bl.get("energy", True):
            total *= np.sqrt(e0 / max(float((total ** 2).sum()), 1e-12))
        info.append(dict(kind="blend", knn=bl["knn"], w=bl.get("w", 0.2), energy=bl.get("energy", True), with_neighbours=int(has.sum()),
                         uncovered_full=int(unc.sum())))
    # promoter floor (CRISPRi neighbour knockdown) and on-target knockdown, applied once to the blend
    pcfg = TM.Cfg(promoter=spec.get("promoter", True))
    _, prom, tidx = TM.predict_lfc(targets, TM.Cfg(weights={k: 0 for k in TM.Cfg().weights}, promoter=pcfg.promoter))
    delta, ov = TM.split(total, prom, tidx, 1.0, pcfg)
    return delta, ov, info


def main():
    spec = json.loads(pathlib.Path(sys.argv[1]).read_text()); out = pathlib.Path(sys.argv[2])
    base = pathlib.Path(spec.get("controls_dir", str(OFF)))
    ctxs = spec.get("contexts", ["A", "B", "C"])
    targets = pd.read_csv(base / "pert_counts.csv").iloc[:, 0].astype(str).tolist()
    genes = pd.read_csv(base / "gene_names.csv").iloc[:, 0].astype(str).tolist()
    assert list(TM.GENES) == genes, "gene axis differs from the official 18,533 axis"
    nT, nG = len(targets), len(genes); n_rows = len(ctxs) * nT * CELLS
    t0 = time.time(); stats = {}
    obs_pert = np.empty(n_rows, np.int16); obs_ctx = np.empty(n_rows, np.int8)
    audit = bool(spec.get("audit_only"))     # audits never touch `out`: in-memory HDF5, results in *_audit.json (Codex R16)
    with (h5py.File(str(out) + ".audit", "w", driver="core", backing_store=False) if audit else h5py.File(out, "w")) as f:
        f.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        Xg = f.create_group("X")
        Xg.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": np.array([n_rows, nG], np.int64)})
        dd = Xg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(2 ** 20,))
        di = Xg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(2 ** 20,))
        dp = Xg.create_dataset("indptr", shape=(n_rows + 1,), dtype=np.int64); dp[0] = 0
        nnz = 0; row = 0
        for ci, ctx in enumerate(ctxs):
            D, OV, info = context_delta(spec, ctx, targets)
            np.save(out.with_name(out.stem + f"_override_{ctx}.npy"), OV.astype(np.float32))
            np.save(out.with_name(out.stem + f"_delta_{ctx}.npy"), D.astype(np.float32))
            M = load_controls(ctx, base); assert M.shape[1] == nG
            # stream target by target (flat memory; one RNG stream per context, identical to tv2_emit.emit)
            tots = []
            if spec["mode"] == "dual_mult":        # dual-moment emission (fh_dual; pinned sums x gamma, top-K DE re-allocation)
                dual = spec.get("dual", {}); AG = None
                if dual.get("agree"):
                    za = np.load(dual["agree"]); pa = {str(t): i for i, t in enumerate(za["targets"])}
                    AG = np.stack([za["agree"][pa[t]] if t in pa else np.full(TM.NG, 0.5, np.float32) for t in targets])
                g0 = dual.get("gamma", 1.0); gsum = g0; kinfo = None
                sm = spec.get("sum", {})               # SUM-side transforms (DE side keeps g0 * D)
                Dbase = D.copy()
                if sm.get("parts_w"):                 # the SUMS from a different part mix (e.g. more S4); DE side keeps D
                    import copy as _copy
                    spec2 = _copy.deepcopy(spec); spec2.pop("sum", None)
                    for pt, w_ in zip(spec2["parts"], sm["parts_w"]):
                        pt["w"] = float(w_)
                    Dbase, _, info2 = context_delta(spec2, ctx, targets)
                    print(ctx, "sum-side part weights", sm["parts_w"], flush=True)
                Dover = None                          # the override recipe BEFORE sum.resid (for de.from_sum_override)
                Dsl = None                            # main recipe with the same-line replacement (DE side without from_sum_override)
                ov_spec = sm.get("spec_override_ctx", {}).get(ctx, sm.get("spec_override"))   # per-context override (final
                if ov_spec:                           # phase: somatic lines drop KOLF, stem lines keep it) else the global one
                    import copy as _copy              # the SUMS from a different recipe (parts / blend replaced; S10k: S4 and
                                                      # fusion / neighbour library without KOLF); the DE side keeps D.
                    # NB: an override that should NOT inherit the main spec's per-context weights must carry its own
                    # "cfg_ctx_parts" -- spec3.update() only replaces the keys the override names.
                    spec3 = _copy.deepcopy(spec); spec3.pop("sum", None); spec3.update(_copy.deepcopy(ov_spec))
                    Dbase, _, info3 = context_delta(spec3, ctx, targets)
                    np.save(out.with_name(out.stem + f"_dbase_{ctx}.npy"), Dbase.astype(np.float32))   # 2026-10-02: sum override BEFORE resid (tlrn feature)
                    Dover = Dbase.copy()
                    dsrc_ = spec.get("de", {}).get("src")   # A-2 (2026-09-30): the DE-side snapshot from a DENOISED source library
                    if dsrc_:                                # (program-space projection, ml_denoise_sources.py); sums, budget
                        S0_ = TM.SRC                         # reference, overrides and the agree file stay on the original sources
                        bad_ = [f_.name for f_ in pathlib.Path(dsrc_["dir"]).iterdir() if f_.is_symlink() and not f_.exists()]
                        assert not bad_, f"de.src dir has BROKEN symlinks (2026-09-30 S14fk bug): {bad_[:5]}"
                        try:
                            TM.SRC = pathlib.Path(dsrc_["dir"]); TM._cache.clear()
                            spec_de = _copy.deepcopy(spec3)
                            if dsrc_.get("knn") and spec_de.get("blend"):   # knn: one path, or {context: path} (final phase:
                                kn_ = dsrc_["knn"][ctx] if isinstance(dsrc_["knn"], dict) else dsrc_["knn"]   # stem vs somatic)
                                spec_de["blend"] = dict(spec_de["blend"], knn=kn_)
                            Dover, _, _ = context_delta(spec_de, ctx, targets)
                        finally:
                            TM.SRC = S0_; TM._cache.clear()
                        cz_ = float((Dover * Dbase).sum() / (np.linalg.norm(Dover) * np.linalg.norm(Dbase) + 1e-12))
                        stats.setdefault("de_src", {})[ctx] = dict(dir=dsrc_["dir"], knn=dsrc_.get("knn"), cos_to_sum=round(cz_, 5))
                        print(ctx, "DE-side snapshot from the denoised source library", dsrc_, "cosine to the sum side", round(cz_, 5), flush=True)
                    print(ctx, "sum-side recipe override", sorted(ov_spec), [i.get("path", i.get("knn", i.get("kind"))) for i in info3], flush=True)
                sl = spec.get("same_line", {}).get(ctx)
                if sl:                                # SAME-LINE REPLACEMENT (FINDINGS 20.11/20.13): for targets a same-line LOW-MOI
                    zs = TM.load(sl["source"]); cen_ = np.nan_to_num(zs["center_c50"])
                    w_ = float(np.clip(float(sl.get("w", 1.0)), 0.0, 1.0))   # screen covers (>= 20 cells): sums and DE delta
                    Dsl = D.copy()                    # <- (1-w) D + w S |D| / |S| on the genes the SOURCE MEASURES (baseline kept
                    arrs_ = [Dbase, Dsl] + ([Dover] if Dover is not None else [])   # elsewhere), norms over genes neither
                    ncov = 0                          # overridden nor the target's own; the budget reference D and OV are NOT
                                                      # touched (Codex R15); Dsl = main recipe + replacement (DE side unless
                                                      # de.from_sum_override, which uses Dover); uncovered rows stay bit-equal
                    for k_, t_ in enumerate(targets):
                        i_ = zs["pos"].get(TM.canon(t_))
                        if i_ is None or zs["ncell"][i_] < 20:
                            continue
                        v_ = np.log((zs["fp"][i_].astype(np.float64) * 1e6 + 50) / (zs["fc"][i_].astype(np.float64) * 1e6 + 50)) - cen_
                        m_ = np.isfinite(v_) & ~np.isfinite(OV[k_])
                        jt_ = TM.GIC.get(TM.canon(t_))
                        if jt_ is not None: m_[jt_] = False
                        v_ = np.where(m_, v_, 0.0); nv_ = np.linalg.norm(v_)
                        if nv_ <= 0:
                            continue
                        for arr_ in arrs_:
                            na_ = np.linalg.norm(np.where(m_, arr_[k_], 0.0))
                            arr_[k_] = np.where(m_, (1 - w_) * arr_[k_] + w_ * v_ * (na_ / nv_), arr_[k_])
                        ncov += 1
                    stats.setdefault("same_line", {})[ctx] = dict(source=sl["source"], w=w_, covered=ncov)
                    print(ctx, f"same-line replacement from {sl['source']} (w {w_}) on {ncov} targets", flush=True)
                if sm.get("resid"):                   # learned no-embedding residual on the sums (ml_dr_train.py noesm_strict, seed mean):
                    zr = np.load(sm["resid"]["file"]) # Dbase[:, G] += alpha * r * |Dbase_t,G| / |X_t|   (sch_s10k_eval.py k2)
                    assert [str(x) for x in zr["targets"]] == list(targets)
                    gR = zr[f"gG_{ctx}"]; xn = zr[f"xnorm_{ctx}"].astype(np.float64)
                    sc_ = np.where(xn > 1e-6, np.linalg.norm(Dbase[:, gR], axis=1) / np.maximum(xn, 1e-6), 0.0)
                    Dbase = Dbase.copy(); Dbase[:, gR] += float(sm["resid"]["alpha"]) * zr[f"r_{ctx}"].astype(np.float64) * sc_[:, None]
                    print(ctx, "sum-side residual", sm["resid"], "genes", len(gR), "median scale", round(float(np.median(sc_)), 3), flush=True)
                tot_c = np.asarray(M.sum(0)).ravel(); cpm_c = tot_c / tot_c.sum() * 1e6
                Vdf = wdf = None
                if sm.get("ctrl_defl"):              # 2026-10-11 CONTROL-CELL PC deflation of the weaker targets' sums (public hint
                    cd_ = sm["ctrl_defl"]            # William v37a; harness `edefl`, ERA tools/ctrl_deflate_eval*.py): V = top-k
                    from sklearn.utils.extmath import randomized_svd   # axes of THIS context's controls (log1p CP10K, gene-centred,
                    Xc = M.tocsr().astype(np.float64); Lc_ = np.asarray(Xc.sum(1)).ravel()   # 6,000 cells, as tools/ctrl_pcs.py);
                    Xc = Xc.multiply(1e4 / np.maximum(Lc_, 1)[:, None]).tocsr(); Xc.data = np.log1p(Xc.data)
                    ix_ = np.random.default_rng(int(cd_.get("seed", 0))).choice(Xc.shape[0], min(int(cd_.get("ncells", 6000)), Xc.shape[0]), replace=False)
                    Dc_ = np.asarray(Xc[ix_].todense()); Dc_ -= Dc_.mean(0); del Xc
                    Vdf = randomized_svd(Dc_, n_components=int(cd_["k"]), n_iter=5, random_state=0)[2].astype(np.float64); del Dc_
                    za_ = np.load(cd_.get("file", sm.get("tamp", {}).get("file")))   # weak = the lower `weak` quantile of the
                    assert [str(x) for x in za_["targets"]] == [str(x) for x in targets], "ctrl_defl: target order"   # tamp amplitude
                    ag_ = za_[f"agree_{ctx}"].astype(np.float64); ag_ = np.nan_to_num(ag_, nan=float(np.nanmedian(ag_)))
                    amp_ = ((np.argsort(np.argsort(ag_, kind="stable"), kind="stable") + 1) / len(ag_)) ** float(cd_.get("beta", 0.75))
                    wdf = np.where(amp_ <= np.quantile(amp_, float(cd_["weak"])), float(cd_["eta"]), 0.0)
                    stats.setdefault("ctrl_defl", {})[ctx] = dict(k=int(cd_["k"]), eta=float(cd_["eta"]), weak=float(cd_["weak"]), n=int((wdf > 0).sum()))
                    print(ctx, f"ctrl_defl k {int(cd_['k'])} eta {cd_['eta']} on {int((wdf > 0).sum())}/{len(wdf)} targets", flush=True)
                def sum_side(use_ew, use_tamp=True):
                    X = Dbase.copy()
                    if sm.get("snr"):                 # source-family SNR weighting (tv2_snr_prep.py)
                        zs = np.load(sm["snr"]); ps = {str(x): i for i, x in enumerate(zs["targets"])}
                        X = X * np.stack([zs["snr"][ps[x]] if x in ps else np.full(TM.NG, 0.5, np.float32) for x in targets]).astype(np.float64)
                    if use_ew and sm.get("ew"):       # expression down-weighting (1 + cpm/c0)^-k, cpm of THIS context's controls
                        c0, kap = sm["ew"]; wew = (1 + cpm_c / float(c0)) ** (-float(kap))
                        X = X * (wew / np.average(wew, weights=np.maximum(cpm_c, 1e-3)))[None, :]
                    if sm.get("center"):              # partial centring over targets, override entries excluded
                        mt = np.nanmean(np.where(np.isfinite(OV), np.nan, X), 0)
                        X = X - float(sm["center"]) * np.nan_to_num(mt)[None, :]
                    if sm.get("pcr"):                 # partial removal of the top-k source-library principal axes (tv2_pc_prep.py)
                        pc = sm["pcr"]; V = np.load(pc["file"])["Vt"][: int(pc["k"])]; P_ = X @ V.T; lam = float(pc["lambda"])
                        if pc.get("guard"):           # panels whose OWN predicted energy sits on those axes (essential-gene-like
                            lo_, hi_ = pc["guard"]    # responses: harness RN 0.32-0.40 vs H1/KOLF 0.07-0.10, official 0.037) keep them
                            f_ = float((P_ ** 2).sum() / max(float((X ** 2).sum()), 1e-12))
                            lam = lam * float(np.clip((hi_ - f_) / (hi_ - lo_), 0.0, 1.0))
                            print(ctx, f"pcr guard: energy on top-{int(pc['k'])} library axes {f_:.4f} -> lambda {lam:.3f}", flush=True)
                        X = X - lam * P_ @ V
                    if use_tamp and sm.get("tamp"):   # 2026-10-01 per-target amplitude a_t = rank(cross-family agreement_t)^beta (PDS-invariant,
                        za = np.load(sm["tamp"]["file"]); ta_ = [str(x) for x in za["targets"]]   # harness `ta`, FINDINGS §23.8);
                        assert ta_ == [str(x) for x in targets], "tamp: target order"           # the budget match below restores B_pred
                        ag_ = za[f"agree_{ctx}"].astype(np.float64); ag_ = np.nan_to_num(ag_, nan=float(np.nanmedian(ag_)))
                        rk_ = (np.argsort(np.argsort(ag_, kind="stable"), kind="stable") + 1) / len(ag_)
                        X = X * (rk_ ** float(sm["tamp"]["beta"]))[:, None]
                        print(ctx, f"tamp beta {sm['tamp']['beta']}: amplitude median {float(np.median(rk_ ** float(sm['tamp']['beta']))):.3f}", flush=True)
                    if Vdf is not None:               # row-linear, so its place relative to the tamp amplitude is immaterial; the
                        X = X - wdf[:, None] * ((X @ Vdf.T) @ Vdf)   # budget match below restores B_pred
                    return X
                Dsum = sum_side(True)
                if sm.get("ew") and sm.get("ew_adaptive") and dual.get("kink"):
                    kk0 = dual["kink"]                # the ew decision is taken on the sums BEFORE the tamp amplitude (S15 audit 2026-10-02:
                    Dchk = sum_side(True, False) if sm.get("tamp") else Dsum   # tamp shrinks targets, which flipped ew back on)
                    gk_ew, _ = EM.kink_gamma(M, g0 * Dchk, OV, CELLS, factor=float(kk0.get("factor", 1.17)), margin=float(kk0.get("margin", 0.8)))
                    if gk_ew < 1.0:                   # budget already binding at gamma 1 -> expression weighting would overshoot it
                        Dsum = sum_side(False); print(ctx, f"ew dropped (kink with ew {gk_ew:.3f} < 1)", flush=True)
                    else:
                        print(ctx, f"ew kept (kink with ew {gk_ew:.3f})", flush=True)
                if sm:
                    print(ctx, "sum-side", {k: sm[k] for k in sm if k not in ("snr", "spec_override")}, "snr" in sm, flush=True)
                    np.save(out.with_name(out.stem + f"_dsum_{ctx}.npy"), Dsum.astype(np.float32))   # audit: sums before gamma
                if dual.get("kink"):                  # sums at max(floor, #348 budget kink); DE side stays at g0 * D
                    kk = dual["kink"]; fac_, mar_ = float(kk.get("factor", 1.17)), float(kk.get("margin", 0.8))
                    if kk.get("budget_match") and sm:  # transformed sums keep the UNtransformed recipe's budget position
                        gkb, cr = EM.kink_gamma(M, g0 * D, OV, CELLS, factor=fac_, margin=mar_); gb = max(float(kk.get("floor", 1.0)), gkb)
                        Bt = EM.budget_B(M, g0 * D, OV, CELLS)(gb); s_ = EM.budget_match(M, g0 * Dsum, OV, CELLS, Bt)
                        gk = gkb; gsum = g0 * s_
                        kinfo = dict(kink_base=round(gkb, 3), gamma_base=round(gb, 3), gamma_sum=round(gsum, 3), c_real=round(cr, 2), budget_match=True)
                    else:
                        gk, cr = EM.kink_gamma(M, g0 * Dsum, OV, CELLS, factor=fac_, margin=mar_)
                        gsum = g0 * max(float(kk.get("floor", 1.0)), gk); kinfo = dict(kink=round(gk, 3), gamma_sum=round(gsum, 3), c_real=round(cr, 2))
                    print(ctx, "kink", kinfo, flush=True)
                Dde = D; dcfg = spec.get("de", {})
                if Dsl is not None and not dcfg.get("from_sum_override"):
                    Dde = Dsl                         # same-line replacement reaches the DE side (covered rows only; Codex R16)
                if dcfg.get("from_sum_override"):     # DE side from the SUM-side recipe override (S11: no-KOLF S4 / fusion / library +
                    assert Dover is not None, "de.from_sum_override needs sum.spec_override"   # min_cells ensemble), snapshot
                    Dde = Dover                       # BEFORE sum.resid. The budget reference (kink / budget_match on D), the
                    stats.setdefault("de_side", {})[ctx] = "sum_override"   # overrides OV and the agree file stay the main
                    np.save(out.with_name(out.stem + f"_dde_{ctx}.npy"), Dde.astype(np.float32))   # recipe's (harness bro/brs arms,
                    print(ctx, "DE side from the sum-side recipe override; budget reference + overrides stay the main recipe", flush=True)
                if dcfg.get("pcr"):                   # DE-side removal of the top-k library axes (harness dpr<k>x<l>[g]); sums unchanged
                    pc = dcfg["pcr"]; V = np.load(pc["file"])["Vt"][: int(pc["k"])]; P_ = Dde @ V.T; lam = float(pc["lambda"])
                    if pc.get("guard"):
                        lo_, hi_ = pc["guard"]; f_ = float((P_ ** 2).sum() / max(float((Dde ** 2).sum()), 1e-12))
                        lam = lam * float(np.clip((hi_ - f_) / (hi_ - lo_), 0.0, 1.0))
                        print(ctx, f"de pcr guard: energy {f_:.4f} -> lambda {lam:.3f}", flush=True)
                    Dde = Dde - lam * P_ @ V
                if dcfg.get("drc"):               # within each target the top 5 % of genes by |DE lfc| x a and the bottom
                    a_, b_ = [float(v) for v in dcfg["drc"]]   # 70 % x b. The SUM target is untouched, so PDS / MSE are
                    Dde = Dde.copy()                           # EXACTLY invariant and only NMAE / reach move (harness `drc`,
                    for t_ in range(Dde.shape[0]):             # FINDINGS section 24; 3/3 panels positive, ride-along only)
                        idx_ = np.flatnonzero(~np.isfinite(OV[t_]))
                        o_ = idx_[np.argsort(-np.abs(Dde[t_, idx_]), kind="stable")]
                        nt_ = int(round(0.05 * len(o_))); nb_ = int(round(0.70 * len(o_)))
                        Dde[t_, o_[:nt_]] *= a_
                        if nb_: Dde[t_, o_[-nb_:]] *= b_
                    print(ctx, f"de drc {a_:g} / {b_:g}", flush=True)
                it = EM.emit_dual_iter(M, gsum * Dsum, CELLS, spec.get("seed", 0) + 1000 * ci,
                                       K=dual.get("K", 3000), override=OV, depth_up=dual.get("depth_up", 0.0), agree=AG,
                                       cpc=bool(dual.get("cpc", False)), top_frac=float(dual.get("nt", 0.1)),
                                       ci=int(dual.get("ci", 0)), lsb=bool(dual.get("lsb", False)),
                                       lcp=float(dual.get("lcp", 0.0)),
                                       lfc_de=(g0 * Dde if (kinfo or sm or dcfg.get("pcr") or dcfg.get("drc")) else None))
            else:
                it = EM.emit_iter(M, D, CELLS, spec.get("seed", 0) + 1000 * ci, spec["mode"], override=OV)
            if spec.get("audit_only"):            # parity audits: stop after the delta / sum / DE dumps and the kink (no emission)
                stats[ctx] = dict(parts=info, audit_only=True, kink=(kinfo if spec["mode"] == "dual_mult" else None),
                                  gamma_sum=(float(gsum) if spec["mode"] == "dual_mult" else None))
                print(ctx, "audit only", json.dumps(stats[ctx].get("kink")), f"{time.time()-t0:.0f}s", flush=True)
                del M
                continue
            for ti, b in enumerate(it):
                tot = np.asarray(b.sum(1)).ravel(); tots.append(tot)
                if tot.max() > MAX_COUNTS_PER_CELL:
                    raise SystemExit(f"cell total {tot.max()} exceeds cap")
                b.sort_indices(); k = b.nnz
                dd.resize((nnz + k,)); dd[nnz:nnz + k] = b.data
                di.resize((nnz + k,)); di[nnz:nnz + k] = b.indices
                dp[row + 1:row + 1 + CELLS] = nnz + b.indptr[1:].astype(np.int64)
                obs_pert[row:row + CELLS] = ti; obs_ctx[row:row + CELLS] = ci
                nnz += k; row += CELLS
            tot = np.concatenate(tots)
            stats[ctx] = dict(parts=info, delta_rms=float(np.sqrt((D ** 2).mean())), median_cell_total=float(np.median(tot)))
            print(ctx, json.dumps(stats[ctx]), f"nnz={nnz:,} {time.time()-t0:.0f}s", flush=True)
            del M
        if spec.get("audit_only"):
            out.with_name(out.stem + "_audit.json").write_text(json.dumps(dict(spec=spec, stats=stats, audit_only=True), indent=1))
            print("AUDIT_ONLY_DONE", flush=True)
            return
        if row != n_rows or nnz > MAX_STORED:
            raise SystemExit(f"rows {row}/{n_rows} nnz {nnz:,}")
        og = f.create_group("obs")
        og.attrs.update({"encoding-type": "dataframe", "encoding-version": "0.2.0", "_index": "_index",
                         "column-order": np.array(["target_gene", "context"], dtype=h5py.string_dtype("utf-8"))})
        ix = og.create_dataset("_index", data=np.array([f"c{i}" for i in range(n_rows)], dtype=object), dtype=h5py.string_dtype("utf-8"))
        ix.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
        for name, codes, cats in [("target_gene", obs_pert, targets), ("context", obs_ctx, list(ctxs))]:
            g = og.create_group(name)
            g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
            dc = g.create_dataset("categories", data=np.array(cats, dtype=object), dtype=h5py.string_dtype("utf-8"))
            dc.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
            dcod = g.create_dataset("codes", data=codes); dcod.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
        vg = f.create_group("var")
        vg.attrs.update({"encoding-type": "dataframe", "encoding-version": "0.2.0", "_index": "_index",
                         "column-order": np.array([], dtype=h5py.string_dtype("utf-8"))})
        dv = vg.create_dataset("_index", data=np.array(genes, dtype=object), dtype=h5py.string_dtype("utf-8"))
        dv.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
        for k in ["layers", "obsm", "obsp", "varm", "varp", "uns"]:
            f.create_group(k).attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
    out.with_name(out.stem + "_build.json").write_text(json.dumps(dict(spec=spec, stats=stats, nnz=int(nnz), n_rows=n_rows,
                                                                        seconds=round(time.time() - t0, 1)), indent=1))
    print(f"wrote {out} rows={n_rows} nnz={nnz:,} {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
