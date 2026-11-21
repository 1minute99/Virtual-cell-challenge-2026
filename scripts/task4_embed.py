#!/usr/bin/env python
"""Task 4b - embed control cells of every context three ways (forward passes only).

Representations, one vector per context:
  sclong   : frozen scLong (1.01B params) encoder, mean-pooled over genes then cells
  state_se : frozen STATE SE-600M cell embedding, mean-pooled over cells
  pca      : PCA fit on pooled control cells of all contexts, projected + mean-pooled

Contexts: VCC validation A/B/C + Replogle k562/rpe1/jurkat/hepg2.
No fine-tuning anywhere. Logs VRAM, throughput and max batch size.

Run: ~/miniconda3/envs/AIVC/bin/python scripts/task4_embed.py --rep pca
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd, h5py
from scipy import sparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "work/cellline_cache")
VCC = os.path.join(ROOT, "data/2026_AIVC _validation_controls")
OUT = os.path.join(ROOT, "results"); os.makedirs(OUT, exist_ok=True)
CONTROL = "non-targeting"
LINES = ["k562", "rpe1", "jurkat", "hepg2"]
CTXS = ["A", "B", "C"]


def norm_counts(counts):
    tot = counts.sum(1, keepdims=True)
    return np.log1p(counts / np.maximum(tot, 1e-9) * 1e4).astype(np.float32)


def gather_controls(n_per_ctx, seed=0):
    """Return {context: (cells, |R|) log1p-CP10K matrix on the shared gene space R}."""
    G = np.load(os.path.join(ROOT, "work/gene_spaces.npz"), allow_pickle=True)
    R_rep, R_vcc, R_scl = G["R_in_rep"], G["R_in_aivc"], G["R_in_sclong"]
    rng = np.random.default_rng(seed)
    out = {}
    for L in LINES:
        X = np.load(os.path.join(CACHE, f"X_{L}.npy"), mmap_mode="r")
        p = np.asarray(np.load(os.path.join(CACHE, f"pert_{L}.npy"), allow_pickle=True), dtype=str)
        idx = np.where(p == CONTROL)[0]
        idx = np.sort(rng.choice(idx, min(n_per_ctx, len(idx)), replace=False))
        # cached X is log1p-normalised on the full 6546 space; undo then renormalise on R
        out[f"replogle:{L}"] = norm_counts(np.expm1(np.asarray(X[idx])[:, R_rep].astype(np.float32)))
        del X
    for c in CTXS:
        with h5py.File(os.path.join(VCC, f"context_{c}.h5ad"), "r") as h:
            g = h["X"]; shape = tuple(g.attrs["shape"])
            M = sparse.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        idx = np.sort(rng.choice(shape[0], min(n_per_ctx, shape[0]), replace=False))
        out[f"vcc:{c}"] = norm_counts(np.asarray(M[idx].todense())[:, R_vcc].astype(np.float32))
    return out, R_scl, len(G["sclong_genes"])


def rep_pca(ctrl, n_comp=50, seed=0):
    from sklearn.decomposition import PCA
    keys = list(ctrl)
    Xall = np.concatenate([ctrl[k] for k in keys])
    pca = PCA(n_components=min(n_comp, Xall.shape[0], Xall.shape[1]), random_state=seed)
    Z = pca.fit_transform(Xall)
    emb, i = {}, 0
    for k in keys:
        n = len(ctrl[k]); emb[k] = Z[i:i + n].mean(0); i += n
    return emb, dict(explained_var=float(pca.explained_variance_ratio_.sum()),
                     n_components=int(pca.n_components_))


def rep_sclong(ctrl, R_scl, nS, batch=4):
    import torch, pickle
    # torch>=2.13 removed torch.qr; scLong's vendored performer still calls it.
    def _qr(A, some=True, out=None):
        return torch.linalg.qr(A, mode="reduced" if some else "complete")
    torch.qr = _qr  # unconditional: torch>=2.13 keeps the attr but raises
    sys.path.insert(0, os.path.join(ROOT, "external/scLong"))
    from performer_pytorch_cont.ding_models import DualEncoderSCFM
    CK = os.path.join(ROOT, "checkpoints/sclong/")
    hp = pickle.load(open(CK + "gocont_4096_48m_pretrain_1b_mix.pkl", "rb"))
    hp["gene2vec_file"] = CK + "selected_gene2vec_27k.npy"
    model = DualEncoderSCFM(**hp).cuda().eval()
    sd = torch.load(CK + "gocont_4096_48m_pretrain_1b_mix_2024-02-05_16-23-37.pth",
                    map_location="cpu")["model_state_dict"]
    model.load_state_dict(sd)
    for p in model.parameters():
        p.requires_grad = False
    nparam = sum(p.numel() for p in model.parameters())
    col = torch.tensor(np.asarray(R_scl), dtype=torch.long, device="cuda")
    emb, stats = {}, {}
    torch.cuda.reset_peak_memory_stats()
    t0 = time.time(); ncell = 0
    for k, X in ctrl.items():
        acc = []
        for i in range(0, len(X), batch):
            xb = X[i:i + batch]
            x = torch.zeros(len(xb), nS, device="cuda")
            x[:, col] = torch.tensor(xb, device="cuda")
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                enc = model(x, return_encodings=True)["merged_decodings"]
            acc.append(enc[:, col, :].float().mean(1).cpu().numpy())
            ncell += len(xb)
        emb[k] = np.concatenate(acc).mean(0)
    el = time.time() - t0
    stats = dict(params=nparam, peak_vram_GB=torch.cuda.max_memory_allocated() / 1e9,
                 cells=ncell, seconds=el, cells_per_s=ncell / el, batch_size=batch,
                 gpu=torch.cuda.get_device_name(0))
    return emb, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep", default="pca", choices=["pca", "sclong", "state_se"])
    ap.add_argument("--n-per-ctx", type=int, default=256)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    ctrl, R_scl, nS = gather_controls(a.n_per_ctx, a.seed)
    print({k: v.shape for k, v in ctrl.items()})
    if a.rep == "pca":
        emb, stats = rep_pca(ctrl, seed=a.seed)
    elif a.rep == "sclong":
        emb, stats = rep_sclong(ctrl, R_scl, nS, a.batch)
    else:
        from state_se_embed import rep_state_se
        emb, stats = rep_state_se(ctrl, R_scl, batch=a.batch)
    keys = list(emb)
    np.save(os.path.join(OUT, f"task4_emb_{a.rep}.npy"), np.stack([emb[k] for k in keys]))
    json.dump(dict(contexts=keys, stats=stats),
              open(os.path.join(OUT, f"task4_emb_{a.rep}.json"), "w"), indent=1, default=str)
    print(a.rep, "stats:", stats)
    print("dim:", np.stack([emb[k] for k in keys]).shape)


if __name__ == "__main__":
    main()
