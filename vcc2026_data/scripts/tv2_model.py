"""transfer_v2 fusion: fold-change transfer with source centering, source weights, promoter prior.

predict_lfc(targets, cfg, exclude=...) -> (n_targets, 18533) natural-log fold changes on the official axis,
to be applied MULTIPLICATIVELY to a destination's own control expression (never as an additive shift).
"""
import gzip, json, pathlib, re, sys
from dataclasses import dataclass, field
import numpy as np, pandas as pd

R = pathlib.Path("/home/wonmin-kim/AIVC/vcc2026_data")
SRC = R / "reports/transfer_v2_2026-09-23/sources"
SRC_EXTRA = R / "reports/transfer_v2_2026-09-23/sources_new"   # contexts added after the original library was built
import os as _os
if _os.environ.get("TV2_SRC_DIR"):            # 2026-10-03: holdout rehearsals point every consumer at a source dir WITHOUT the
    SRC = SRC_EXTRA = pathlib.Path(_os.environ["TV2_SRC_DIR"])   # held-out families (unset = the production library, unchanged)
sys.path.append(str(R / "scripts"))
import vcclib as V

OFF = pathlib.Path("/home/wonmin-kim/AIVC/data/2026_AIVC _validation_controls")
GENES = pd.read_csv(OFF / "gene_names.csv").iloc[:, 0].astype(str).values
NG = len(GENES)
GM = V.GeneMap(); canon = GM.canon
GIC = {}
for _i, _g in enumerate(GENES):
    GIC.setdefault(canon(_g), _i)
DROP_TARGETS = set()   # rehearsal: targets whose measurements are deleted from EVERY source
FAMILY = {"K562": "K562", "HCT116": "HCT116", "HEK293T": "HEK293T", "CD4_Rest": "CD4", "CD4_Stim8hr": "CD4",
          "CD4_Stim48hr": "CD4", "KOLF": "KOLF", "FENG_T": "FENG", "FENG_GWF": "FENG", "FENG_GWN": "FENG", "H1": "H1",
          "JURKAT_UNSTIM": "JURKAT", "JURKAT_STIM": "JURKAT", "VIP_K562": "VIP",
          # Replogle/Nadig lines: no source file exists for the A/B/C panel (they cover 0/300 official targets), so these
          # are inert there; the final-phase rebuild can create them for a NEW panel (2026-09-27, coverage insurance).
          "RN_HepG2": "HEPG2", "RN_Jurkat": "RNJURKAT", "RN_RPE1": "RPE1",
          # NEW CONTEXTS 2026-09-29 (`ml_newsrc_build.py`, `ml_dld1_pseudobulk.py`; files in `sources_new/`, found via
          # SRC_EXTRA). A context is only worth the ~+0.003 that `ml_ctx_curve.py` measured if it COVERS the panel:
          # DLD1 has 4,714 targets / 67 official, the scPerturb screens have 1-6, so their value is measured not assumed.
          "DLD1": "DLD1", "H1GW": "H1GW", "FrangiehIzar2021_RNA": "MELANOMA", "TianKampmann2019_iPSC": "TKIPSC",
          "TianKampmann2019_day7neuron": "TKNEURON", "TianKampmann2021_CRISPRi": "TKNEURONI",
          "WesselsSatija2023": "THP1", "PapalexiSatija2021_eccite_RNA": "THP1",
          "ShifrutMarson2018": "PRIMARYT", "SunshineHein2023": "SUNSHINE"}
GTF = pathlib.Path("/tmp/claude-1000/-home-wonmin-kim-AIVC/e6a8589e-1f1c-4fe1-add0-3faf80af9be9/scratchpad/gencode47.gtf.gz")
GTF_KEEP = R / "reports/transfer_v2_2026-09-23/gencode47_tss_official.json"


def _src_path(s):
    """A source npz lives in the original library or in `sources_new/` (contexts added later)."""
    p1 = SRC / f"{s}.npz"
    return p1 if p1.exists() else SRC_EXTRA / f"{s}.npz"


@dataclass
class Cfg:
    c: float = 10.0                    # CPM pseudocount of the source fold change
    center: bool = True                # subtract each source's genome-wide mean response
    center_frac: float = 1.0           # fraction of that mean to subtract
    weights: dict = field(default_factory=lambda: {"K562": 2.0, "HCT116": 1.0, "HEK293T": 1.0, "CD4": 0.5,
                                                   "KOLF": 1.0, "FENG": 0.5, "H1": 2.0})
    min_cells: int = 20
    ncw: float = 0.0        # >0: weight each (source, target) row by n/(n+ncw) instead of equally. The fusion combined
                            # rows with EQUAL weight and only a hard min_cells cut; averaging over min_cells thresholds
                            # (ml_hp_ensemble*.py `ens mc`) improved PDS and MSE together on H1/Orion/KOLF, and this is
                            # the smooth version of that step-function weighting.
    promoter: bool = True
    promoter_f0: float = 0.15          # remaining fraction at <=500 bp, ramping to 1 at 5 kb (log distance)
    clip: float = 3.0                  # |LFC| clip before amplitude
    kd_target: float = np.log(0.2)     # on-target knockdown (excluded from every metric; realism only)
    center_robust: float = 0.0   # >0: replace the stored per-source MEAN center by a TRIMMED mean over the source's
                                 # own rows, trimming this fraction at each tail (0.5 = median). The stored center is an
                                 # arithmetic mean, so a handful of huge-effect targets can drag it; |center|/|row| runs
                                 # 0.07 (HCT116, 9.8k rows) to 0.70 (PapalexiSatija, 25 rows). 0 = the shipped center.
    fwins: float = 0.0      # >0: before the weighted family average, winsorize each family's value at (t,g) to the
                            # [fwins, 1-fwins] quantiles ACROSS the families that measured it. The weighted mean is
                            # sensitive to one aberrant family (K562 carries weight 8 of ~12.5); this clips it toward
                            # the others per GENE instead of globally. 0 = bit-identical to the shipped combination.
    agree_shrink: float = 0.0          # >0: shrink each (t,g) by m^2/(m^2 + agree_shrink*var_between_families/n_fam)
    min_kd: float = 0.0                # >0: drop a source row whose own on-target knockdown 1-(fp+1)/(fc+1) is below this
                                       #     (only where the target is >= 5 CPM in that source's controls; 2026-09-27)
    zshrink: float = 0.0               # >0: per-SOURCE per-coordinate shrinkage v * z^2/(z^2 + zshrink) before fusion, z = v / sd,
                                       #     sd = robust across-target spread of that gene's FC in that source x sqrt(n_med/n_cells)
    xw: float = 0.0                    # >0: per-GENE source weighting by the gene's expression in that source row's matched
                                       #     controls, e = (cpm / (cpm + c))^xw (floor 1e-3): the c-pseudocount FC of a gene the
                                       #     source barely expresses is ~0 and used to DILUTE the other families at its full
                                       #     family weight (K562 = 63 % of the official fusion). Within a family rows are
                                       #     e-weighted; the family enters with w x its mean e. 0 = bit-identical (2026-09-30).
                                       #     MEASURED (ml_xw_diag.py, fusion only, offline): per-target cosine with the truth
                                       #     H1 +0.003 / OR +0.001 but discrimination H1 +0.001/-0.003, OR -0.008/-0.013
                                       #     (xw 0.5 / 1) -> NOT adopted, not worth a real-scorer run.


_cache = {}


def load(name):
    if name not in _cache:
        z = np.load(_src_path(name), allow_pickle=True)
        _cache[name] = {k: z[k] for k in z.files}
        _cache[name]["pos"] = {str(t): i for i, t in enumerate(z["targets"])}
    return _cache[name]


def tss_table():
    if GTF_KEEP.exists():
        return {int(k): v for k, v in json.loads(GTF_KEEP.read_text()).items()}
    feat = pd.read_csv(R / "metadata/vcc_official/official_features.csv")
    ens2idx = {e: i for i, e in zip(feat.feature_index, feat.ensembl_id) if isinstance(e, str)}
    tss = {}
    for line in gzip.open(GTF, "rt"):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "transcript":
            continue
        gid = re.search(r'gene_id "([^".]+)', f[8]).group(1)
        i = ens2idx.get(gid)
        if i is None:
            continue
        tss.setdefault(int(i), []).append([f[0], int(f[3]) if f[6] == "+" else int(f[4])])
    GTF_KEEP.write_text(json.dumps(tss))
    return tss


def promoter_neighbours(targets, max_bp=5000):
    tss = tss_table()
    bychrom = {}
    for i, L in tss.items():
        for c, t in L:
            bychrom.setdefault(c, []).append((t, i))
    for c in bychrom:
        bychrom[c].sort()
    starts = {c: np.array([a for a, _ in v]) for c, v in bychrom.items()}
    out = {}
    for tg in targets:
        i = GIC.get(canon(tg))
        best = {}
        for c, t in tss.get(i, []) if i is not None else []:
            arr = bychrom[c]; p = np.searchsorted(starts[c], t - max_bp)
            for a, j in arr[p:]:
                if a > t + max_bp:
                    break
                if j != i:
                    best[j] = min(best.get(j, 10 ** 9), abs(a - t))
        out[tg] = best
    return out


_noise = {}


def noise_sd(s, z, cen, c):
    """Robust per-gene noise SD of source s's fold changes (MAD over all its targets, rows scaled to the median cell
    count). Most targets leave most genes unchanged, so the across-target spread is dominated by sampling noise."""
    key = (s, float(c), float(np.nansum(cen)))
    if key not in _noise:
        V = np.log((z["fp"].astype(np.float64) * 1e6 + c) / (z["fc"].astype(np.float64) * 1e6 + c)) - cen[None, :]
        nc = np.maximum(z["ncell"].astype(np.float64), 1); nmed = float(np.median(nc))
        V = V * np.sqrt(nc / nmed)[:, None]
        med = np.nanmedian(V, 0); mad = np.nanmedian(np.abs(V - med[None, :]), 0) * 1.4826
        _noise[key] = (np.maximum(np.nan_to_num(mad, nan=0.0), 1e-3), nmed)
    return _noise[key]


_RCEN = {}
def robust_center(name, z, c, q, chunk=2000):
    """Trimmed-mean center over the source's OWN rows (q per tail; q>=0.5 -> median). Chunked over genes."""
    key = (name, float(c), float(q))
    if key in _RCEN: return _RCEN[key]
    fp = z["fp"]; fc = z["fc"]; G = fp.shape[1]; out = np.full(G, np.nan)
    for a0 in range(0, G, chunk):
        b0 = min(a0 + chunk, G)
        V = np.log((fp[:, a0:b0].astype(np.float64) * 1e6 + c) / (fc[:, a0:b0].astype(np.float64) * 1e6 + c))
        V = np.where(np.isfinite(V), V, np.nan)
        if q >= 0.5:
            out[a0:b0] = np.nanmedian(V, 0)
        else:
            lo = np.nanquantile(V, q, axis=0); hi = np.nanquantile(V, 1.0 - q, axis=0)
            out[a0:b0] = np.nanmean(np.where((V >= lo[None]) & (V <= hi[None]), V, np.nan), 0)
    _RCEN[key] = out
    return out


def predict_lfc(targets, cfg: Cfg, exclude=(), return_parts=False):
    """targets: list of symbols. exclude: families (or source names) to leave out (held-out context)."""
    _fw_m, _fw_h, _fw_w = [], [], []
    num = np.zeros((len(targets), NG)); den = np.zeros((len(targets), NG))
    fam_src = {}
    sq = np.zeros((len(targets), NG)); nfam = np.zeros((len(targets), NG))
    for s, fam in FAMILY.items():
        if fam in exclude or s in exclude or np.all(np.asarray(cfg.weights.get(fam, 0), dtype=np.float64) <= 0) or not _src_path(s).exists():
            continue
        fam_src.setdefault(fam, []).append(s)
    for fam, srcs in fam_src.items():
        # family = mean over its sub-sources that measure (t, g); family weight split evenly
        fnum = np.zeros((len(targets), NG)); fden = np.zeros((len(targets), NG))
        fcnt = np.zeros((len(targets), NG)) if cfg.xw > 0 else None
        for s in srcs:
            z = load(s)
            cen = (robust_center(s, z, cfg.c, cfg.center_robust) if cfg.center_robust > 0
                   else z[f"center_c{cfg.c:g}"]) if cfg.center else np.zeros(NG)
            cen = np.where(np.isfinite(cen), cen, 0.0) * cfg.center_frac
            if cfg.zshrink > 0:
                sd0, nmed = noise_sd(s, z, cen, cfg.c)
            for k, t in enumerate(targets):
                i = z["pos"].get(canon(t))
                if i is None or z["ncell"][i] < cfg.min_cells or canon(t) in DROP_TARGETS:
                    continue
                if cfg.min_kd > 0:
                    jt = GIC.get(canon(t))
                    if jt is not None:
                        fc_, fp_ = float(z["fc"][i][jt]), float(z["fp"][i][jt])
                        if np.isfinite(fc_) and np.isfinite(fp_) and fc_ * 1e6 >= 5 and 1 - (fp_ * 1e6 + 1) / (fc_ * 1e6 + 1) < cfg.min_kd:
                            continue
                v = np.log((z["fp"][i].astype(np.float64) * 1e6 + cfg.c) / (z["fc"][i].astype(np.float64) * 1e6 + cfg.c)) - cen
                if cfg.zshrink > 0:
                    z2 = (v / (sd0 * np.sqrt(nmed / max(float(z["ncell"][i]), 1.0)))) ** 2
                    v = v * z2 / (z2 + cfg.zshrink)
                m = np.isfinite(v)
                rw = 1.0 if cfg.ncw <= 0 else float(z["ncell"][i]) / (float(z["ncell"][i]) + cfg.ncw)
                if cfg.xw > 0:
                    cp_ = z["fc"][i].astype(np.float64) * 1e6
                    e_ = np.clip(np.where(np.isfinite(cp_), cp_, 0.0) / (np.where(np.isfinite(cp_), cp_, 0.0) + cfg.c), 1e-3, 1.0) ** cfg.xw
                    fnum[k, m] += rw * e_[m] * v[m]; fden[k, m] += rw * e_[m]; fcnt[k, m] += rw
                else:
                    fnum[k, m] += rw * v[m]; fden[k, m] += rw
        has = fden > 0
        # `maximum(fden, 1)` was safe while fden was an integer ROW COUNT (>= 1 wherever has); with fractional
        # `ncw` row weights it clamps the denominator to 1 and turns reliability weighting into amplitude
        # SHRINKAGE (a 20-cell single-source family at ncw=200 returned 0.091 x its response). Bit-identical at ncw=0.
        fam_mean = np.where(has, fnum / np.maximum(fden, 1e-12), 0.0)
        w = np.asarray(cfg.weights[fam], dtype=np.float64)
        if w.ndim == 1:                  # PER-TARGET family weights (learned weighting, ml_learn_weights2.py); a scalar = shipped
            w = w[:, None]
        if cfg.xw > 0:                   # the family enters with its mean expression weight at (t, g)
            w = w * np.where(fcnt > 0, fden / np.maximum(fcnt, 1e-12), 0.0)
        if cfg.fwins > 0:
            _fw_m.append(fam_mean); _fw_h.append(has); _fw_w.append(w)
        else:
            num += w * fam_mean; den += w * has; sq += w * fam_mean ** 2; nfam += has
    if cfg.fwins > 0:                        # winsorize across families at each (t, g), then combine as before
        M = np.stack(_fw_m); Hh = np.stack(_fw_h)
        Mn = np.where(Hh, M, np.nan)
        lo = np.nanquantile(Mn, cfg.fwins, axis=0); hi = np.nanquantile(Mn, 1.0 - cfg.fwins, axis=0)
        ok = np.isfinite(lo) & np.isfinite(hi)
        M = np.where(Hh & ok[None], np.clip(M, lo[None], hi[None]), M)
        for fm, hs, w in zip(M, Hh, _fw_w):
            num += w * fm; den += w * hs; sq += w * fm ** 2; nfam += hs
    L = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
    if cfg.agree_shrink > 0:
        var = np.where(den > 0, sq / np.maximum(den, 1e-12) - L ** 2, 0.0).clip(0)
        L = L * L ** 2 / (L ** 2 + cfg.agree_shrink * var / np.maximum(nfam, 1) + 1e-12)
    L = np.clip(L, -cfg.clip, cfg.clip)
    covered = (den > 0).any(1)
    prom = np.zeros_like(L)          # promoter prior kept separate so the amplitude does not scale it
    if cfg.promoter:
        nb = promoter_neighbours(targets)
        for k, t in enumerate(targets):
            for j, d in nb[t].items():
                f = cfg.promoter_f0 + (1 - cfg.promoter_f0) * np.clip((np.log10(max(d, 1)) - np.log10(500)) / 1.0, 0, 1)
                if f < 1:
                    prom[k, j] = np.log(f)
    tidx = np.array([GIC.get(canon(t), -1) for t in targets])
    if return_parts:
        return L, prom, tidx, covered
    return L, prom, tidx


def split(L, prom, tidx, alpha, cfg: Cfg):
    """(delta, override): delta = alpha * L for the emission transform; override = natural-log fold changes that
    must be applied MULTIPLICATIVELY whatever the transform (promoter-neighbour floor, on-target knockdown),
    NaN elsewhere. Never pass a knockdown through the log1p cell shift: it zeroes the gene and blows up log2FC."""
    delta = alpha * L
    ov = np.full(L.shape, np.nan)
    if cfg.promoter:
        m = (prom < 0) & (prom < delta)
        ov[m] = prom[m]
    rows = np.flatnonzero(tidx >= 0)
    ov[rows, tidx[rows]] = cfg.kd_target
    return delta, ov


def combine(L, prom, tidx, alpha, cfg: Cfg):
    """Final LFC = alpha * transferred LFC, floored by the promoter prior, on-target knockdown."""
    out = alpha * L
    out = np.minimum(out, np.where(prom < 0, prom, np.inf)) if cfg.promoter else out
    rows = np.flatnonzero(tidx >= 0)
    out[rows, tidx[rows]] = cfg.kd_target
    return out
