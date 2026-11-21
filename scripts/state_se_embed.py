"""STATE SE-600M cell embeddings (forward pass only, no fine-tuning).

Writes one h5ad of control cells per context on the shared gene space (mapped
ENSG -> HGNC symbol, which is what SE's ESM2 protein-embedding vocabulary is
keyed on), runs `state emb transform`, then mean-pools to one vector per context.
"""
import os, subprocess, sys
import numpy as np, anndata as ad, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SE = os.path.join(ROOT, "checkpoints/state/SE-600M")
TMP = os.path.join(ROOT, "work/state_se"); os.makedirs(TMP, exist_ok=True)
PY = sys.executable


def rep_state_se(ctrl, R_scl, batch=32):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from genemap import GeneMap
    gm = GeneMap()
    G = np.load(os.path.join(ROOT, "work/gene_spaces.npz"), allow_pickle=True)
    ens = [str(x) for x in G["R"]]
    syms = [gm.ens2sym.get(e) for e in ens]
    keep = np.array([i for i, s in enumerate(syms) if s])
    sym_keep = [syms[i] for i in keep]
    # de-duplicate symbols (SE keys on symbol)
    seen, uniq = set(), []
    for i, s in enumerate(sym_keep):
        if s not in seen:
            seen.add(s); uniq.append(i)
    keep, sym_keep = keep[uniq], [sym_keep[i] for i in uniq]
    print(f"SE gene space: {len(sym_keep)} symbols from {len(ens)} ENSG")

    emb, stats = {}, {}
    import torch
    for k, X in ctrl.items():
        tag = k.replace(":", "_")
        fin = os.path.join(TMP, f"{tag}.h5ad")
        fout = os.path.join(TMP, f"{tag}_emb.npy")
        # SE expects counts-like input; our matrices are log1p-CP10K -> invert
        M = np.expm1(X[:, keep]).astype(np.float32)
        A = ad.AnnData(M, obs=pd.DataFrame(index=[f"c{i}" for i in range(len(M))]),
                       var=pd.DataFrame(index=sym_keep))
        A.write_h5ad(fin)
        if not os.path.exists(fout):
            cmd = [PY, "-m", "state", "emb", "transform",
                   "--checkpoint", os.path.join(SE, "se600m_epoch16.ckpt"),
                   "--input", fin, "--output", fout,
                   "--protein-embeddings", os.path.join(SE, "protein_embeddings.pt"),
                   "--batch-size", str(batch)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stdout[-3000:]); print(r.stderr[-3000:])
                raise RuntimeError(f"state emb transform failed for {k}")
        Z = np.load(fout)
        emb[k] = Z.mean(0)
        print(f"  {k}: {Z.shape}")
    stats = dict(model="arcinstitute/SE-600M se600m_epoch16.ckpt",
                 n_genes=len(sym_keep), batch_size=batch,
                 gpu=torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu")
    return emb, stats
