"""Shared loaders for plan `zero-shot-limits-context-complexity` (Tasks 6-21). Frozen in Task 4 (2026-11-05).

ASSUMPTIONS (stated, checked in task04_freeze_inputs.py):
  A1  The 7 harness panels' cell-level truth is the file the ERA packs were built from, i.e. the `path` in
      vcc2026_data/reports/final_harness_2026-09-24/panels.json (fh_harness.PANELS). That is the default (`which="truth"`).
      `which="full"` gives the larger raw file the panel was subsampled from (RN: the scPerturb h5ad; KOLF: harness/kolf_real;
      OR_HCT116X: real_thin, all 300 targets) -- use it only where more cells are needed (split halves) and say so.
      NOTE: RN_K562's truth is subsampled from Replogle K562 **GWPS**, not K562 *essential* as the plan text says.
  A2  Control cells are obs[label] == "non-targeting" (RN raw files: obs["gene"] == "non-targeting").
  A3  Gene identity across panels is resolved on the official 18,533-gene VCC axis with vcclib.GeneMap.canon
      (the harness's own mapping, fh_harness.GIC). G_cap ("G intersection") = official-axis genes that are on all 7 panels'
      var axes AND have >= 1 UMI in that panel's control cells.
  A4  Files are opened backed ("r"); only the requested rows are materialised. Never call .to_memory() on a whole panel.
No randomness here (split helpers take an explicit seed)."""
from __future__ import annotations
import functools, json, pathlib, sys
import numpy as np

ROOT = pathlib.Path("/home/wonmin-kim/AIVC"); V = ROOT / "vcc2026_data"; RP = V / "reports"
RESULTS = ROOT / "output/data-analysis/results"
PACKS = RP / "era_2026-10-07/packs"
HARNESS_PANELS_JSON = RP / "final_harness_2026-09-24/panels.json"
OFFICIAL_GENES_CSV = ROOT / "data/2026_AIVC _validation_controls/gene_names.csv"
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]

# which="full": larger parent file (label column, control label). None = the truth file is already the full panel.
FULL = {
    "H1": None,
    "OR_HCT116X": (RP / "sc_hybrid_2026-09-25/orion_hct116/real_thin.h5ad", "target", "non-targeting"),
    "KOLF": (V / "harness/kolf_real.h5ad", "target", "non-targeting"),
    "RN_K562": (V / "raw/scperturb/ReplogleWeissman2022_K562_gwps.h5ad", "gene", "non-targeting"),
    "RN_HepG2": (V / "raw/scperturb/NadigOConner2024_hepg2.h5ad", "gene", "non-targeting"),
    "RN_Jurkat": (V / "raw/scperturb/NadigOConner2024_jurkat.h5ad", "gene", "non-targeting"),
    "RN_RPE1": (V / "raw/scperturb/ReplogleWeissman2022_rpe1.h5ad", "gene", "non-targeting"),
}
# Paths the plan text names that do NOT exist / differ (logged as a Task-4 deviation):
PLAN_PATH_NOTES = {
    "KOLF": "plan: checkpoints/state/kolf_real.h5ad (absent); harness truth = reports/b19_sweep2_2026-09-18/full_metrics/KOLF/real_disjoint.h5ad; parent = vcc2026_data/harness/kolf_real.h5ad",
    "RN_Jurkat": "plan: checkpoints/state/jurkat_real.h5ad (absent); harness truth = mock/Jurkat/panel300/real.h5ad (from NadigOConner2024_jurkat); vcc2026_data/harness/jurkat_real.h5ad is an older 12k-NTC Jurkat panel",
    "RN_K562": "plan: raw/scperturb/ReplogleWeissman2022_K562_essential.h5ad; harness truth is subsampled from ReplogleWeissman2022_K562_gwps.h5ad (mock manifest `source`)",
    "OR_HCT116X": "plan: orion_hct116/real_thin.h5ad (300 targets); harness truth = orion_hct116/real_panel.h5ad (263 targets)",
}


@functools.lru_cache(None)
def harness_panels() -> dict:
    return json.loads(HARNESS_PANELS_JSON.read_text())


def truth_path(panel: str, which: str = "truth") -> pathlib.Path:
    if which == "truth" or FULL.get(panel) is None:
        return pathlib.Path(harness_panels()[panel]["path"])
    return pathlib.Path(FULL[panel][0])


def _label_col(panel, which):
    return ("target", "non-targeting") if which == "truth" or FULL.get(panel) is None else FULL[panel][1:]


def open_panel(panel: str, which: str = "truth"):
    """Backed AnnData (mode 'r'). Caller closes with `.file.close()`."""
    import anndata as ad
    return ad.read_h5ad(truth_path(panel, which), backed="r")


def labels(panel: str, which: str = "truth") -> np.ndarray:
    a = open_panel(panel, which)
    try:
        return a.obs[_label_col(panel, which)[0]].astype(str).to_numpy()
    finally:
        a.file.close()


def _rows(a, idx):
    """Materialise rows idx (sorted ints) of a backed AnnData as CSR float32."""
    from scipy import sparse
    idx = np.sort(np.asarray(idx))
    X = a.X[idx]          # backed CSR dataset -> csr; dense h5py dataset -> ndarray (needs increasing idx)
    return sparse.csr_matrix(X, dtype=np.float32)


def cells(panel: str, target: str | None = None, which: str = "truth", idx=None):
    """CSR counts for one target (or `idx` rows). target=None and idx=None -> controls."""
    a = open_panel(panel, which)
    try:
        col, ctrl = _label_col(panel, which)
        lab = a.obs[col].astype(str).to_numpy()
        if idx is None:
            idx = np.flatnonzero(lab == (ctrl if target is None else target))
        return _rows(a, idx), a.var_names.astype(str).to_numpy()
    finally:
        a.file.close()


def controls(panel: str, which: str = "truth"):
    """(CSR control counts, gene names)."""
    return cells(panel, None, which)


def split_halves(n: int, seed: int, k: int = 2) -> list[np.ndarray]:
    """Random k-way partition of range(n) (cell-level split halves / three-way control split)."""
    p = np.random.default_rng(seed).permutation(n)
    return [np.sort(x) for x in np.array_split(p, k)]


def pseudobulk(panel: str, which: str = "truth", targets=None, cp10k=True):
    """Per-target summed counts (T x G) and the control sum, streamed target by target. Returns dict."""
    a = open_panel(panel, which)
    try:
        col, ctrl = _label_col(panel, which); lab = a.obs[col].astype(str).to_numpy()
        tg = [t for t in (targets or sorted(set(lab) - {ctrl}))]
        G = a.n_vars; S = np.zeros((len(tg), G)); N = np.zeros(len(tg), int)
        for k, t in enumerate(tg):
            ii = np.flatnonzero(lab == t)
            if len(ii):
                S[k] = np.asarray(_rows(a, ii).sum(0)).ravel(); N[k] = len(ii)
        ci = np.flatnonzero(lab == ctrl); C = np.asarray(_rows(a, ci).sum(0)).ravel()
        out = dict(targets=np.array(tg), genes=a.var_names.astype(str).to_numpy(), sums=S, n=N, ctrl_sum=C, n_ctrl=len(ci))
    finally:
        a.file.close()
    if cp10k:
        out["lcp10k"] = np.log1p(S / np.maximum(S.sum(1, keepdims=True), 1) * 1e4)
        out["ctrl_lcp10k"] = np.log1p(C / max(C.sum(), 1) * 1e4)
    return out


# ---------------- ERA packs -----------------
PACK_KINDS = ("inputs", "score", "ctrlpc", "ctrlvar", "aux", "post", "sunny_cos0", "sunny_nce0")


def pack(panel: str, kind: str = "inputs"):
    """Lazy NpzFile of reports/era_2026-10-07/packs/<kind>_<panel>.npz (arrays load on key access)."""
    assert kind in PACK_KINDS, kind
    return np.load(PACKS / f"{kind}_{panel}.npz", allow_pickle=False)


# ---------------- Tahoe -----------------
TAHOE_PROCESSED = V / "data/processed/drug_context_pretrain/tahoe.npz"
TAHOE_TABLE = RP / "simple_residual_2026-09-16/data/tahoe/table.npz"     # DELETED by the 2026-10-11 cleanup (not recoverable on disk)


def tahoe(which: str = "processed"):
    """Lazy NpzFile. 'processed' = condition-level responses (13,011 conditions x 18,150 genes, 48 lines; dcp input).
    'table' = the original simple_residual Tahoe table -- MISSING since 2026-10-11 (raises)."""
    p = TAHOE_PROCESSED if which == "processed" else TAHOE_TABLE
    if not p.exists():
        raise FileNotFoundError(f"{p} is missing (deleted by reports/cleanup_2026-10-11; see output/code/INPUTS.md)")
    return np.load(p, allow_pickle=False)


# ---------------- gene axis / G_cap -----------------
@functools.lru_cache(None)
def official_genes() -> np.ndarray:
    import pandas as pd
    return pd.read_csv(OFFICIAL_GENES_CSV).iloc[:, 0].astype(str).to_numpy()


@functools.lru_cache(None)
def _canon():
    sys.path.append(str(V / "scripts"))
    import vcclib
    gm = vcclib.GeneMap(); gic = {}
    for i, g in enumerate(official_genes()):
        gic.setdefault(gm.canon(g), i)
    return gm.canon, gic


def to_official(genes) -> np.ndarray:
    """Panel gene symbols -> official-axis index (-1 = not on the axis); the harness's mapping (fh_harness.GIC)."""
    canon, gic = _canon()
    return np.array([gic.get(canon(str(g)), -1) for g in genes])


def g_cap(recompute: bool = False) -> np.ndarray:
    """Official-axis gene symbols measured in all 7 panels' controls (definition: on the var axis and >= 1 control UMI)."""
    f = RESULTS / "g_cap.json"
    if f.exists() and not recompute:
        return np.array(json.loads(f.read_text())["genes"])
    sets = {}
    for p in PANELS7:
        C, genes = controls(p); tot = np.asarray(C.sum(0)).ravel(); om = to_official(genes)
        sets[p] = set(om[(om >= 0) & (tot > 0)].tolist())
    keep = sorted(set.intersection(*sets.values()))
    g = official_genes()[keep]
    f.write_text(json.dumps(dict(definition="official-axis genes on all 7 panels' var axes with >=1 UMI in each panel's controls",
                                 n=len(g), per_panel={p: len(s) for p, s in sets.items()}, genes=g.tolist()), indent=0))
    return g
