"""Final-phase pipeline (D/E/F on 2026-10-22; dry-run on A/B/C): identity -> source weights -> sources -> build.

  tv2_final.py --controls-dir DIR --contexts D,E,F --out OUT_DIR [--rebuild-sources] [--build]

DIR must hold context_<X>.h5ad, gene_names.csv, pert_counts.csv (the official controls package layout).
Steps
 1. identity: Spearman against HPA cell lines (tv2_identity) -> identity.json
 2. source weights per context: base recipe weights (S6 recipe); a source family whose own cell line is the
    context's identity (HPA rho >= SAME_LINE_RHO against that source line) is boosted xSAME_LINE_BOOST, but ONLY
    if the source passes the knockdown-quality gate (median on-target LFC <= KD_GATE on covered targets);
    CD4 (primary T) x2 when the identity is T-lymphoid.
 3. sources: tv2_sources + tv2_jurkat_source rebuilt with the new panel ($TV2_TARGETS) when --rebuild-sources
 4. spec per context -> tv2_submit_chain.sh (build, vcc prep, sha); upload stays a separate deliberate step.
Recipe (S6-like, without the S4 refit which is official-panel specific): 0.5 x equal-weight uncentred fold
changes (c=100) x3 + 0.5 x transfer_v2 fusion (K562 8, centred, c=50) x1.6; gated cellshift; promoter prior.
--blend W (default 0.3): DepMap-neighbour blend (tv2_final_knn.py; S8 = blend 0.2 scored +0.004 officially); targets
with no source get the neighbour estimate only. --dual gamma,K: dual-moment emission (fh_dual; +0.03 on H1/KOLF
harness), rule status unconfirmed.
"""
import argparse, json, os, pathlib, subprocess, sys
import numpy as np, pandas as pd
from cell_eval2.moments import jackknife_correction  # noqa: F401  (before scripts/ joins sys.path)
sys.path.append("/home/wonmin-kim/AIVC/vcc2026_data/scripts")
import tv2_identity as ID
import tv2_model as TM

R = TM.R
SAME_LINE_RHO, SAME_LINE_BOOST, KD_GATE = 0.85, 3.0, -0.8
# stem-vs-somatic control score (ml_kolf_policy.py markers): H1 +2.26, KOLF +2.34 | official A/B/C -1.33/-5.11/-2.43, Orion
# -2.23, RN K562/HepG2/RPE1 -1.65/-1.88/-2.80, RN Jurkat +0.71 (lymphoid lines lack the mesenchymal/epithelial somatic
# markers) -> threshold +1.5 separates the two stem panels from every somatic one (2026-09-29).
STEM_MARK = ["POU5F1", "NANOG", "SOX2", "LIN28A", "DNMT3B", "TDGF1", "SALL4", "ZFP42", "L1TD1", "USP44", "TERF1", "DPPA4"]
SOMA_MARK = ["VIM", "CDH2", "FN1", "COL1A1", "SPARC", "ANXA2", "S100A4", "TAGLN", "CAV1", "LGALS1", "THBS1", "SERPINE1"]
STEM_THR = 1.5
SAME_LINE_MARGIN = 0.03          # identity margin required for same-line replacement unless the runner-up is the same line
CROSS_LABEL = {"HEK293", "HCT 15"} # identity -> related source line: half weight
# HPA line name -> same-line source store (LOW-MOI only: GSE249595 Jurkat is high-MOI noise and deliberately absent)
SAME_LINE_SRC = {"K-562": "K562", "HCT 116": "HCT116", "HEK293": "HEK293T", "HCT 15": "DLD1", "Hep-G2": "RN_HepG2",
                 "hTERT-RPE1": "RN_RPE1"}          # exact tv2_identity/HPA names; HCT-15 and DLD-1 come from the same patient


KOLF_Y_ON, KOLF_Y_OFF = ["DDX3Y", "UTY", "USP9Y"], ["KDM5D", "EIF1AY"]


def kolf_signature(ctrl_h5ad, genes):
    """2026-10-03 (FINDINGS 23.16): KOLF2.1J carries Yq AZFa genes (DDX3Y / UTY / USP9Y 180-550 CPM) but NOT KDM5D / EIF1AY (0 CPM),
    unlike H1 / H1GW (KDM5D 8-40, EIF1AY 60-186), the official male lines and every female line -- a platform-independent
    presence / absence fingerprint for "this stem destination is KOLF2.1J" (our genome-scale KOLF screen is then same-line)."""
    import h5py
    from scipy import sparse
    with h5py.File(ctrl_h5ad, "r") as h:
        g = h["X"]; X = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(int(v) for v in g.attrs["shape"]))
    cs = np.asarray(X.sum(0)).ravel(); cpm = cs / cs.sum() * 1e6; gi = {x: i for i, x in enumerate(genes)}
    on = [cpm[gi[x]] for x in KOLF_Y_ON if x in gi]; off = [cpm[gi[x]] for x in KOLF_Y_OFF if x in gi]
    ok = len(on) == 3 and len(off) == 2 and float(np.mean(on)) >= 20.0 and max(off) < 1.0
    return bool(ok), dict(on=[round(float(v), 1) for v in on], off=[round(float(v), 2) for v in off])


def stem_score(ctrl_h5ad, genes):
    import h5py
    from scipy import sparse
    with h5py.File(ctrl_h5ad, "r") as h:
        g = h["X"]; X = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(int(v) for v in g.attrs["shape"]))
    cs = np.asarray(X.sum(0)).ravel(); cpm = cs / cs.sum() * 1e6; gi = {x: i for i, x in enumerate(genes)}
    ns_, no_ = sum(x in gi for x in STEM_MARK), sum(x in gi for x in SOMA_MARK)
    # Without this the score silently averages whatever markers survive the gene axis: the RN panels keep only
    # 1-3 of the 12 stem markers and score +0.71 for Jurkat against -1.33 for the SAME line on the official axis
    # (2026-10-04 audit). The official axis carries all 24, so a shortfall here means the controls are wrong.
    assert ns_ >= 10 and no_ >= 10, f"stem_score: only {ns_}/12 stem and {no_}/12 somatic markers on this gene axis"
    pm = lambda ns: float(np.mean([np.log1p(cpm[gi[x]]) for x in ns if x in gi]))
    return pm(STEM_MARK) - pm(SOMA_MARK)
SOURCE_LINES = {"K562": ["K-562"], "HCT116": ["HCT 116"], "HEK293T": ["HEK293", "HEK 293"], "JURKAT": ["Jurkat E6.1", "JURKAT"],
                "VIP": ["K-562"], "H1": [], "KOLF": [], "FENG": [], "CD4": []}
BASE_OLD = {"K562": 1.0, "HCT116": 1.0, "HEK293T": 1.0, "CD4": 1.0, "KOLF": 1.0, "FENG": 1.0, "H1": 1.0, "JURKAT": 0.0, "VIP": 0.0}
BASE_NEW = {"K562": 8.0, "HCT116": 1.0, "HEK293T": 0.5, "CD4": 1.0, "KOLF": 0.25, "FENG": 0.25, "H1": 2.0, "JURKAT": 0.0, "VIP": 0.0}


def source_kd_quality(fam):
    vals = []
    for s, f in TM.FAMILY.items():
        if f != fam or not TM._src_path(s).exists():          # sources_new too (DLD1, H1GW; Codex R15)
            continue
        z = TM.load(s)
        for t, i in z["pos"].items():
            j = TM.GIC.get(t)
            if j is None:
                continue
            fp, fc = z["fp"][i][j], z["fc"][i][j]
            if np.isfinite(fp) and np.isfinite(fc) and fc * 1e6 >= 10:
                vals.append(np.log((fp * 1e6 + 1) / (fc * 1e6 + 1)))
    return float(np.median(vals)) if vals else None, len(vals)


def context_weights(ident, M, lines):
    """Boosts for one context from its identity result."""
    q_top = {r["line"]: r["spearman"] for r in ident["top"]}
    boosts, notes = {}, []
    for fam, names in SOURCE_LINES.items():
        rho = max([q_top.get(n, -1) for n in names] + [-1])
        if rho >= SAME_LINE_RHO:
            kd, n = source_kd_quality(fam)
            if kd is not None and kd <= KD_GATE:
                boosts[fam] = SAME_LINE_BOOST; notes.append(f"{fam}: same line (rho {rho:.3f}), KD {kd:.2f} passes -> x{SAME_LINE_BOOST}")
            else:
                notes.append(f"{fam}: same line (rho {rho:.3f}) but KD gate failed ({kd}) -> no boost")
    # (the CD4 x2 rule for T-lymphoid identities was REMOVED 2026-09-29: on RN_Jurkat dropping CD4 is +0.0068 PDS, and S12
    #  -- the validated structure -- carries no boosts; Codex R15)
    return boosts, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls-dir", required=True); ap.add_argument("--contexts", default="A,B,C")
    ap.add_argument("--out", required=True); ap.add_argument("--rebuild-sources", action="store_true")
    ap.add_argument("--build", action="store_true"); ap.add_argument("--tag", default="final")
    ap.add_argument("--rn-sources", action="store_true",
                    help="also build the Replogle/Nadig HepG2 / Jurkat / RPE1 sources for the new panel and give them "
                         "--rn-weight in the fusion (they cover 0/300 A/B/C targets but 896 genes overall)")
    ap.add_argument("--rn-weight", type=float, default=1.0, help="fusion weight for each RN family when --rn-sources")
    ap.add_argument("--s4", action="store_true", help="rebuild the S4 refit+FiLM part for these contexts/targets (S6 recipe)")
    ap.add_argument("--fusion-share", type=float, default=0.5)
    ap.add_argument("--identity-boost", action="store_true",
                    help="apply the legacy x3 fusion boost to a source family whose line matches the context identity (off: S12)")
    ap.add_argument("--sum-tamp", type=float, default=0.0,
                    help="2026-10-01: per-target SUM amplitude rank(cross-family agreement)^BETA before the budget match "
                         "(harness token ta<beta>, FINDINGS 23.8; PDS-invariant MSE lever); 0 = off; needs --sum-policy")
    ap.add_argument("--sum-resid", type=float, default=0.0,
                    help="alpha of the no-embedding learned residual on the sums (S10k/S11/S12: 0.25); needs --sum-policy + --blend")
    ap.add_argument("--blend", type=float, default=0.3,
                    help="DepMap-neighbour blend weight (0 = off); targets no screen measured get the neighbour estimate only")
    ap.add_argument("--dual", default="", help="dual-moment emission 'gamma,K[,depth_up[,ag[,cpc[,nt<share>]]]]' (official best S10c: "
                                               "1.0,3000,0.5,ag,cpc; harness best S10f: 1.0,5000,0.5,ag,cpc,nt0.2 with --blend 0.5); empty = gated cellshift. "
                                               "User 2026-09-25: rules do not constrain emission.")
    ap.add_argument("--sum-center", type=float, default=0.0,
                    help="SUM-side partial centring over targets (DE side unchanged); harness 2026-09-27 PDS up on 6/6 panels at 0.5")
    ap.add_argument("--sum-ew", default="", help="'c0,kappa' (e.g. 50,0.25): SUM-side expression down-weighting (1+cpm/c0)^-kappa")
    ap.add_argument("--sum-pcr", default="", help="'k,lambda[g]' (e.g. 3,0.75g): SUM-side removal of lambda x the projection onto "
                    "the top-k source-library principal axes (tv2_pc_prep.py; harness pr<k>x<l>; a hubness correction for PDS)")
    ap.add_argument("--de-pcr", default="", help="'k,lambda[g]': DE-side removal of the top-k source-library axes (harness "
                    "dpr<k>x<l>[g]; sums unchanged; g = same guard as --sum-pcr)")
    ap.add_argument("--sum-snr", action="store_true",
                    help="SUM-side source-family SNR weighting (tv2_snr_prep.py on the final targets; DE side unchanged)")
    ap.add_argument("--kink", default="", help="'factor,margin' (e.g. 1.17,0.8): sums at max(1, margin x #348 budget kink), "
                                            "estimated from the context's controls; needs --dual")
    ap.add_argument("--sum-policy", default="off", choices=["off", "auto", "nokolf"],
                    help="S10k/S11 sum side per context: 'auto' = no-KOLF S4 / fusion / neighbour library for contexts whose "
                         "control stem-vs-somatic score is < STEM_THR (somatic; S10k +0.0061 officially), KOLF kept for stem "
                         "lines; 'nokolf' = every context. Needs --s4 --dual.")
    ap.add_argument("--de-denoise", default="off", choices=["off", "intended", "shipped"],
                    help="S14fk DE side (official 0.23273, 2026-09-30): K562 program-denoised (ml_denoise_apply.py K 100 rho 0.25) "
                         "DE-side snapshot + knn from the denoised library. 'intended' = all families with K562 denoised (harness "
                         "6/6 positive); 'shipped' = what S14fk actually uploaded (the DE-side fusion saw denoised K562 ONLY, "
                         "FINDINGS 22.6). Needs --de-follow-sum. Default off = the S12 path, bit-exact.")
    ap.add_argument("--s4-excl", default="KOLF",
                    help="biologies removed from the SOMATIC contexts' S4 refit (rows + empirical means; saved FiLM residual kept). "
                         "Default KOLF = S4noK via tv2_s4_nokolf_final.py (bit-exact S10k-S12 path); e.g. KOLF,FENG builds "
                         "S4noKF via tv2_s4_excl_final.py (2026-09-30 candidate, FINDINGS 21.9; adopt only if its gate passed).")
    ap.add_argument("--sum-cens", default="", help="LFC pseudocount values added to the sum-side fusion ensemble, e.g. "
                                                   "5,10,20,100 (only 1/5/10/20/50/100 have precomputed centring vectors). "
                                                   "RIDE-ALONG: R21 +0.0024, 3/3 panels (FINDINGS 25.13)")
    ap.add_argument("--de-drc", default="", help="DE-side magnitude calibration 'a,b': per target the top 5 %% of genes by "
                                                 "|DE lfc| x a and the bottom 70 %% x b; sums untouched so PDS/MSE are exactly "
                                                 "invariant. RIDE-ALONG: H1 +0.0018, 3/3 panels (FINDINGS 24)")
    ap.add_argument("--sum-defl", default="", help="'k,eta,weak' (S18 official: 8,0.75,0.5): CONTROL-cell PC deflation of the weaker "
                                                   "`weak` share of targets' sums (top-k axes of each context's own controls; weak = "
                                                   "lower quantile of the --sum-tamp agreement); S18 +0.0070 officially (2026-10-11)")
    ap.add_argument("--fam-weights", default="", help="extra fusion families 'FAM=W,...' (S17: DLD1=1 with --rn-sources); "
                                                      "joins both parts' weights like --rn-sources")
    ap.add_argument("--sum-ens", default="", help="min_cells values for the sum-side fusion ensemble, e.g. 5,10,15,30,50,100,200 "
                                                  "(S11 +0.00296 officially); needs --sum-policy")
    ap.add_argument("--same-line", type=float, default=0.0,
                    help="W > 0: same-line replacement (sums AND DE side) for a context whose HPA identity matches a line we hold a "
                         "LOW-MOI screen for (SAME_LINE_SRC; rho >= 0.85, margin >= 0.03, measured passing KD). RECOMMENDED W = 0.5 "
                         "(H1<-H1GW +0.077 at 0.5 vs +0.074 at 1); non-exact sources are capped: stem fallback H1GW min(W, 0.5), "
                         "cross-label maps min(W/2, 0.5) (FINDINGS 20.20-20.23). Audit each source's provenance / state first.")
    ap.add_argument("--de-follow-sum", action="store_true",
                    help="S12: the DE side follows the sum-side recipe (de.from_sum_override); needs --sum-policy")
    a = ap.parse_args()
    D = pathlib.Path(a.controls_dir); O = pathlib.Path(a.out); O.mkdir(parents=True, exist_ok=True)
    ctxs = a.contexts.split(",")
    targets = pd.read_csv(D / "pert_counts.csv").iloc[:, 0].astype(str).tolist()
    (O / "targets.txt").write_text("\n".join(targets))
    # 1 identity
    idents = {c: ID.infer(D / f"context_{c}.h5ad") for c in ctxs}
    (O / "identity.json").write_text(json.dumps(idents, indent=1))
    for c in ctxs:
        print(c, "->", idents[c]["best"], idents[c]["rho"], "margin", idents[c]["margin"], flush=True)
    # 3 sources for the new panel
    if a.rebuild_sources:
        env = dict(os.environ, TV2_TARGETS=str(O / "targets.txt"))
        py = sys.executable
        subprocess.run([py, "-P", str(R / "scripts/tv2_sources.py")], env=env, check=True)
        if a.rn_sources:                          # Replogle/Nadig lines: 896 extra targets, highest on-target KD (0.85+)
            subprocess.run([py, "-P", str(R / "scripts/tv2_sources.py"), "RN_HepG2", "RN_Jurkat", "RN_RPE1"], env=env, check=True)
        if (R / "raw/jurkat_gse249595").exists():
            subprocess.run([py, "-P", str(R / "scripts/tv2_jurkat_source.py")], env=env, check=True)
    cov = {}
    for fam in ["K562", "HCT116", "HEK293T", "CD4", "KOLF", "FENG", "H1", "JURKAT", "VIP"]:
        n = 0
        for s, f in TM.FAMILY.items():
            if f == fam and (TM.SRC / f"{s}.npz").exists():
                n = max(n, sum(TM.canon(t) in TM.load(s)["pos"] for t in targets))
        cov[fam] = n
    print("target coverage by family:", cov, flush=True)
    # S4 refit+FiLM part (bit-exact rebuild of the S4 code path, verified on A/B/C)
    s4_dir = O / "s4"
    if a.s4 and not all((s4_dir / f"S4_ctx{c}.delta.npy").exists() for c in ctxs):
        subprocess.run(["/usr/bin/bash", str(R / "scripts/tv2_s4_rebuild_run.sh"), str(s4_dir), str(D), ",".join(ctxs)], check=True)
    if a.rn_sources:                              # RN families join both parts' weights (inert if no source file exists)
        for _w in (BASE_OLD, BASE_NEW):
            for _f in ("HEPG2", "RNJURKAT", "RPE1"):
                _w.setdefault(_f, a.rn_weight)
        print(f"RN sources enabled, fusion weight {a.rn_weight} each", flush=True)
    for _kv in filter(None, a.fam_weights.split(",")):   # S17 families (e.g. DLD1=1)
        _f, _v = _kv.split("=")
        for _w in (BASE_OLD, BASE_NEW):
            _w[_f] = float(_v)
        print(f"extra fusion family {_f} weight {_v}", flush=True)
    # 2 weights + 4 spec
    M, _, lines = ID.reference()
    fs = a.fusion_share
    first = ({"w": 1 - fs, "alpha": 3.0, "delta": {c: str(s4_dir / f"S4_ctx{c}.delta.npy") for c in ctxs}} if a.s4 else
             {"w": 1 - fs, "alpha": 3.0, "cfg": {"c": 100.0, "center": False, "weights": BASE_OLD, "promoter": False}})
    spec = {"mode": "shared_cellshift_gated", "seed": 0, "promoter": True, "controls_dir": str(D), "contexts": ctxs,
            "parts": [first, {"w": fs, "alpha": 1.6, "cfg": {"c": 50.0, "center": True, "weights": BASE_NEW, "promoter": False}}],
            "cfg_ctx_parts": {}}
    if a.blend > 0:
        kdir = O / "knn"
        if not (kdir / "knn.npz").exists():
            subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_final_knn.py"), str(O / "targets.txt"), str(kdir)], check=True)
        spec["blend"] = {"knn": str(kdir / "knn.npz"), "w": a.blend, "energy": True, "uncovered_full": True}
    if a.dual:                                   # "gamma,K[,depth_up[,ag[,cpc]]]", e.g. 1.0,3000,0.5,ag,cpc (S10c, official best)
        f = a.dual.split(",")
        spec["mode"] = "dual_mult"; spec["dual"] = {"gamma": float(f[0]), "K": int(f[1])}
        if len(f) > 2:
            spec["dual"]["depth_up"] = float(f[2])
        if len(f) > 3 and f[3] == "ag":
            agf = O / "agree.npz"
            if not agf.exists():
                lib = O / "knn" / "sources"
                subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_agree_prep.py"), str(O / "targets.txt"), str(agf)]
                               + ([str(lib)] if lib.exists() else []), check=True)
            spec["dual"]["agree"] = str(agf)
        if "cpc" in f[2:]:                           # DE-functional pin (S10c official 0.21687, +0.025 over S10b)
            spec["dual"]["cpc"] = True
        for tok in f[2:]:                            # nt<share>: down-gene concentration share (S10d/S10f use nt0.2)
            if tok.startswith("nt"):
                spec["dual"]["nt"] = float(tok[2:])
            if tok.startswith("ci") and tok[2:].isdigit():   # S11/S12 parity (both measured nil officially)
                spec["dual"]["ci"] = int(tok[2:])
            if tok == "lsb":
                spec["dual"]["lsb"] = True
            if tok.startswith("lcp"):                    # 2026-10-01: weaken the DE-side re-allocation on high-expression genes
                spec["dual"]["lcp"] = float(tok[3:])     # (real scorer OR_HCT116X +0.0108 NMAE, H1 +0.0015 reach; PDS/MSE invariant)
    if a.kink:
        assert a.dual, "--kink needs --dual"
        kf = a.kink.split(","); spec["dual"]["kink"] = {"factor": float(kf[0]), "margin": float(kf[1]) if len(kf) > 1 else 0.8}
    sm = {}
    if a.sum_center:
        sm["center"] = a.sum_center
    if a.sum_ew:                                 # applied only where the kink leaves budget (kink with ew >= 1), else dropped
        sm["ew"] = [float(x) for x in a.sum_ew.split(",")]; sm["ew_adaptive"] = True
    if a.sum_pcr:
        pcf = O / "lib_pcs.npz"
        if not pcf.exists():
            subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_pc_prep.py"), str(pcf), "10"], check=True)
        kp_, lp_ = a.sum_pcr.split(","); sm["pcr"] = {"file": str(pcf), "k": int(kp_), "lambda": float(lp_.rstrip("g"))}
        if lp_.endswith("g"):                    # guard: lambda -> 0 when our own predicted energy sits on those axes
            sm["pcr"]["guard"] = [0.15, 0.25]
    if a.sum_snr:
        snf = O / "snr.npz"
        if not snf.exists():
            lib = O / "knn" / "sources"
            subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_snr_prep.py"), str(O / "targets.txt"), str(snf)]
                           + ([str(lib)] if lib.exists() else []), check=True)
        sm["snr"] = str(snf)
    if a.de_pcr:
        assert a.dual, "--de-pcr needs --dual"
        pcf = O / "lib_pcs.npz"
        if not pcf.exists():
            subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_pc_prep.py"), str(pcf), "10"], check=True)
        kd_, ld_ = a.de_pcr.split(","); spec["de"] = {"pcr": {"file": str(pcf), "k": int(kd_), "lambda": float(ld_.rstrip("g"))}}
        if ld_.endswith("g"):
            spec["de"]["pcr"]["guard"] = [0.15, 0.25]
    if sm:
        assert a.dual, "sum-side options need --dual (they act on the pinned sums only)"
        spec["sum"] = sm
        if spec["dual"].get("kink"):              # transformed sums keep the untransformed recipe's budget position
            spec["dual"]["kink"]["budget_match"] = True
    stem_scores = {}
    if a.same_line > 0 and a.sum_policy == "off":   # the stem fallback of --same-line needs the scores too (Codex R15)
        genes_off = pd.read_csv(D / "gene_names.csv").iloc[:, 0].astype(str).tolist()
        stem_scores = {c: stem_score(D / f"context_{c}.h5ad", genes_off) for c in ctxs}
    S4EXCL = sorted({b.strip() for b in a.s4_excl.split(",") if b.strip()})
    S4TAG = "S4noK" if S4EXCL == ["KOLF"] else ("S4no" + "".join(b[0] for b in S4EXCL) + "_" + "-".join(S4EXCL))
    if a.sum_policy != "off":
        assert a.s4 and a.dual, "--sum-policy needs --s4 and --dual"
        assert "sum" in spec, "--sum-policy needs at least one sum-side option (e.g. --sum-center)"
        genes_off = pd.read_csv(D / "gene_names.csv").iloc[:, 0].astype(str).tolist()
        stem_scores = {c: stem_score(D / f"context_{c}.h5ad", genes_off) for c in ctxs}
        if a.sum_policy == "nokolf" or any(v < STEM_THR for v in stem_scores.values()):
            if not all((s4_dir / f"{S4TAG}_ctx{c}.delta.npy").exists() for c in ctxs):   # panel-agnostic, bit-exact on A/B/C
                env = dict(os.environ, TV2_S4_OUT=str(s4_dir), VCC_OFFICIAL_DIR=str(D), TV2_CTX=",".join(ctxs))
                if S4TAG == "S4noK":
                    subprocess.run([str(pathlib.Path.home() / "miniconda3/envs/ct_transfer/bin/python"), str(R / "scripts/tv2_s4_nokolf_final.py")],
                                   env=env, check=True)
                else:
                    env.update(TV2_S4_EXCL=",".join(S4EXCL), TV2_S4_TAG=S4TAG)
                    subprocess.run([str(pathlib.Path.home() / "miniconda3/envs/ct_transfer/bin/python"), str(R / "scripts/tv2_s4_excl_final.py")],
                                   env=env, check=True)
            if a.blend > 0 and not (O / "knn_noK" / "knn.npz").exists():
                subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_final_knn.py"), str(O / "targets.txt"), str(O / "knn_noK"), "KOLF"], check=True)
    if a.de_follow_sum:
        assert a.sum_policy != "off", "--de-follow-sum needs --sum-policy"
        spec.setdefault("de", {})["from_sum_override"] = True
    report = {"identity": {c: idents[c]["best"] for c in ctxs}, "coverage": cov, "boosts": {},
              "blend": spec.get("blend"), "mode": spec["mode"], "dual": spec.get("dual"), "sum": spec.get("sum")}
    for c in ctxs:
        boosts, notes = context_weights(idents[c], M, lines)
        if not a.identity_boost and boosts:      # the x3 same-line fusion boost predates the S12 structure and was never validated
            notes = notes + [f"boost {boosts} NOT applied (opt-in --identity-boost)"]; boosts = {}   # in it (Codex R15/R16)
        report["boosts"][c] = notes
        new_w = {"weights": {k: v * boosts.get(k, 1.0) for k, v in BASE_NEW.items()}}
        spec["cfg_ctx_parts"][c] = ([new_w] if a.s4 else [{"weights": {k: v * boosts.get(k, 1.0) for k, v in BASE_OLD.items()}}, new_w])
        print(c, notes, flush=True)
        if a.same_line > 0:                      # same-line replacement (builder spec key "same_line")
            best, rho = idents[c]["best"], float(idents[c]["rho"]); src, w_sl = None, a.same_line
            top_ = sorted(idents[c]["top"], key=lambda r_: -r_["spearman"]); second = top_[1]["line"] if len(top_) > 1 else None
            same_alias = second is not None and SAME_LINE_SRC.get(second) == SAME_LINE_SRC.get(best)
            if rho >= SAME_LINE_RHO and best in SAME_LINE_SRC and (float(idents[c]["margin"]) >= SAME_LINE_MARGIN or same_alias):
                src = SAME_LINE_SRC[best]
                if best in CROSS_LABEL:          # related, not identical, line (HEK293 -> HEK293T, HCT 15 -> DLD1): half weight, capped at
                    w_sl = min(a.same_line / 2, 0.5)   # 0.5 -- a different-donor related line (DLD1 -> OR_HCT116) is +0.005 at w 0.5 but
                                                       # -0.035 at w 1 (FINDINGS 20.23)
            elif stem_scores.get(c, -9) >= STEM_THR:   # stem fallback: H1GW -> KOLF (another pluripotent line) is +0.021 at w 0.5
                ksig_, kinfo_ = kolf_signature(D / f"context_{c}.h5ad", genes_off)
                report.setdefault("kolf_signature", {})[c] = dict(match=ksig_, **kinfo_); print(c, "KOLF2.1J Y fingerprint", ksig_, kinfo_, flush=True)
                if ksig_:                           # exact line: our genome-scale KOLF2.1J screen at the full same-line weight
                    src, w_sl = "KOLF", a.same_line
                else:
                    src, w_sl = "H1GW", min(a.same_line, 0.5)   # and +0.024 at w 1 (one seed): every non-exact source capped at 0.5 (R17)
            if src and TM._src_path(src).exists():
                fam_ = TM.FAMILY.get(src, src); kd, _ = source_kd_quality(fam_)
                if kd is not None and kd <= KD_GATE:                  # a measured, passing knockdown is required (Codex R15)
                    spec.setdefault("same_line", {})[c] = {"source": src, "w": w_sl}
                    report.setdefault("same_line", {})[c] = dict(source=src, w=w_sl, identity=best, rho=round(rho, 3))
                    print(c, f"same-line source {src} (identity {best}, rho {rho:.3f}) at w {w_sl}", flush=True)
                else:
                    print(c, f"same-line source {src} fails the KD gate ({kd}) -> none", flush=True)
            else:
                print(c, f"no same-line low-MOI source for identity {best} (rho {rho:.3f})", flush=True)
        if a.sum_policy != "off":                # S10k / S11 / S12 structure, per context
            kind = "somatic" if (a.sum_policy == "nokolf" or stem_scores[c] < STEM_THR) else "stem"
            w_ov = {k: v for k, v in new_w["weights"].items() if not (kind == "somatic" and k == "KOLF")}
            fus = {"w": fs, "alpha": 1.6, "cfg": {"c": 50.0, "center": True, "weights": w_ov, "promoter": False}}
            if a.sum_ens:
                fus["cfg_ens"] = [{"min_cells": int(m)} for m in a.sum_ens.split(",")]
            if a.sum_cens:                        # the pseudocount members join the SAME cfg_ens union (harness `cens`)
                ok_ = (1.0, 5.0, 10.0, 20.0, 50.0, 100.0)
                cs_ = [float(v) for v in a.sum_cens.split(",")]
                bad_ = [v for v in cs_ if v not in ok_]
                assert not bad_, f"--sum-cens: c must be one of {ok_} (precomputed center_c<c>), got {bad_}"
                fus["cfg_ens"] = (fus.get("cfg_ens") or []) + [{"c": v} for v in cs_]
            s4f = s4_dir / (f"{S4TAG}_ctx{c}.delta.npy" if kind == "somatic" else f"S4_ctx{c}.delta.npy")
            ov = {"parts": [{"w": 1 - fs, "alpha": 3.0, "delta": {c: str(s4f)}}, fus],
                  "cfg_ctx_parts": {}}             # own weights: the main per-context boosts must NOT re-add KOLF
            if spec.get("blend"):
                ov["blend"] = dict(spec["blend"], knn=str((O / "knn_noK" / "knn.npz") if kind == "somatic" else spec["blend"]["knn"]))
            spec["sum"].setdefault("spec_override_ctx", {})[c] = ov
            report.setdefault("sum_policy", {})[c] = dict(kind=kind, stem_score=round(stem_scores[c], 3), s4=str(s4f), ens=a.sum_ens)
            print(c, f"sum policy: {kind} (stem score {stem_scores[c]:+.3f}); S4 {s4f.name}; ensemble {a.sum_ens or 'none'}", flush=True)
    if a.sum_resid > 0:                          # S10k/S11/S12 no-embedding learned residual on the sums (tv2_final_resid.py,
        assert a.sum_policy != "off" and a.blend > 0, "--sum-resid needs --sum-policy and the knn library (--blend > 0)"
        rf = O / "resid_final.npz"               # bit-exact to resid_S10k.npz on A/B/C); KOLF left out where the policy is somatic
        want = {c: (["KOLF"] if report["sum_policy"][c]["kind"] == "somatic" else []) for c in ctxs}
        tg_now = [t for t in (O / "targets.txt").read_text().split() if t]
        if rf.exists():                          # reuse only if targets AND per-context exclusions match (Codex R16)
            zr_ = np.load(rf); ok_ = list(map(str, zr_["targets"])) == tg_now and all(
                f"excl_{c}" in zr_.files and sorted(map(str, zr_[f"excl_{c}"])) == want[c] for c in ctxs)
            if not ok_:
                print("resid_final.npz does not match this panel / policy -> rebuilding", flush=True); rf.unlink()
        if not rf.exists():
            ctx_arg = ",".join(f"{c}:KOLF" if want[c] else c for c in ctxs)
            subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_final_resid.py"), str(O / "targets.txt"), str(D), ctx_arg,
                            str(O / "knn" / "sources"), str(rf), f"FIN{a.tag}"], check=True)
        spec["sum"]["resid"] = {"file": str(rf), "alpha": a.sum_resid}
        report["sum_resid"] = dict(file=str(rf), alpha=a.sum_resid)
    if a.sum_tamp > 0:                            # per-target agreement amplitude (tv2_tamp_prep.py; somatic contexts drop KOLF)
        assert a.sum_policy != "off", "--sum-tamp needs --sum-policy"
        tf = O / "tamp_final.npz"
        subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_tamp_prep.py"), str(tf),
                        *[f"{c}:{report['sum_policy'][c]['kind']}" for c in ctxs], "--targets", str(O / "targets.txt")], check=True)
        spec["sum"]["tamp"] = {"file": str(tf), "beta": a.sum_tamp}
        report["sum_tamp"] = dict(file=str(tf), beta=a.sum_tamp)
    if a.sum_defl:                                # S18 control-PC deflation (builder sum.ctrl_defl; weak split from the tamp file)
        assert a.sum_tamp > 0, "--sum-defl needs --sum-tamp (its weak/strong split reads the agreement file)"
        kd_, ed_, wd_ = a.sum_defl.split(",")
        spec["sum"]["ctrl_defl"] = {"k": int(kd_), "eta": float(ed_), "weak": float(wd_), "beta": 0.75, "ncells": 6000, "seed": 0}
        report["sum_defl"] = spec["sum"]["ctrl_defl"]
    if a.de_denoise != "off":                     # S14fk DE side (see --de-denoise)
        assert a.de_follow_sum and a.sum_policy != "off", "--de-denoise needs --de-follow-sum and --sum-policy"
        dn_all = O / "de_dnK562"
        if not (dn_all / "K562.npz").exists():
            subprocess.run([sys.executable, "-P", str(R / "scripts/ml_denoise_apply.py"), str(TM.SRC), str(dn_all), "100", "0.25", "K562"],
                           check=True)
        if a.de_denoise == "shipped":
            dn_dir = O / "de_dnK562only"; dn_dir.mkdir(exist_ok=True)
            for n_ in ("K562.npz", "K562.json"):
                if (dn_all / n_).exists() and not (dn_dir / n_).exists():
                    os.symlink((dn_all / n_).resolve(), dn_dir / n_)
        else:
            dn_dir = dn_all
        knn_ctx = {}
        for c in ctxs:
            kind = report["sum_policy"][c]["kind"]
            kd_ = O / ("knn_noK_dnK562" if kind == "somatic" else "knn_dnK562")
            if a.blend > 0:
                if not (kd_ / "knn.npz").exists():
                    src_lib = O / "knn" / "sources"   # reuse the library files already built (the exclusion only acts in the recipe step)
                    if src_lib.exists():
                        (kd_ / "sources").mkdir(parents=True, exist_ok=True)
                        for f_ in src_lib.iterdir():
                            if not (kd_ / "sources" / f_.name).exists():
                                os.symlink(f_.resolve(), kd_ / "sources" / f_.name)
                    env = dict(os.environ, TV2_KNN_DENOISE="K562:100:0.25")
                    subprocess.run([sys.executable, "-P", str(R / "scripts/tv2_final_knn.py"), str(O / "targets.txt"), str(kd_)]
                                   + (["KOLF"] if kind == "somatic" else []), env=env, check=True)
                knn_ctx[c] = str(kd_ / "knn.npz")
        spec.setdefault("de", {})["src"] = {"dir": str(dn_dir), **({"knn": knn_ctx} if knn_ctx else {})}
        report["de_denoise"] = dict(variant=a.de_denoise, dir=str(dn_dir), knn=knn_ctx)
        print("DE side denoised:", a.de_denoise, dn_dir, knn_ctx, flush=True)
    if a.de_drc:                                  # DE-side magnitude calibration; sums untouched => PDS/MSE exactly invariant
        ad_, bd_ = [float(v) for v in a.de_drc.split(",")]
        spec.setdefault("de", {})["drc"] = [ad_, bd_]
        report["de_drc"] = [ad_, bd_]; print(f"DE side drc {ad_:g} / {bd_:g} (ride-along; PDS/MSE invariant)", flush=True)
    (O / f"{a.tag}_spec.json").write_text(json.dumps(spec, indent=1))
    (O / "report.json").write_text(json.dumps(report, indent=1))
    if a.build:
        subprocess.run(["/usr/bin/bash", str(R / "scripts/tv2_submit_chain.sh"), str(O / f"{a.tag}_spec.json"), a.tag], check=True)


if __name__ == "__main__":
    main()
