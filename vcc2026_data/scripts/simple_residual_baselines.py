"""Original empirical source-store aggregation, with complete biological exclusions."""
from functools import lru_cache
import numpy as np
import pandas as pd
import vcclib as V
from simple_residual_data import ROOT, canonical_bio

STORES = ('orion_target_cellline_all','cd4_target_condition_all','kolf_target_batch',
          'replogle_k562_gwps','feng_genomewidescreen_nonfitnessgenes_target','vcc2025_target_batch')
BIO = {'cd4_target_condition_all':'CD4','kolf_target_batch':'KOLF',
       'replogle_k562_gwps':'K562','feng_genomewidescreen_nonfitnessgenes_target':'FENG',
       'vcc2025_target_batch':'H1'}


@lru_cache(None)
def metadata(store):
    p = pd.read_csv(ROOT/'effects'/store/'pairs.csv')
    gm = V.GeneMap()
    p['canonical'] = [gm.canon(x) for x in p.target_symbol.astype(str)]
    if store=='orion_target_cellline_all':
        p['bio'] = [canonical_bio(x) for x in p.cell_context]
    else:
        p['bio'] = BIO[store]
    return p


def empirical(targets,excluded_bios=(),held_targets=()):
    banned = {canonical_bio(x) for x in excluded_bios}
    gm = V.GeneMap()
    held = {gm.canon(t) for t in held_targets}
    canonical = [gm.canon(t) for t in targets]
    stack = []
    for store in STORES:
        p = metadata(store)
        ids = np.flatnonzero(~p.bio.isin(banned) & p.canonical.isin(canonical)
                             & ~p.canonical.isin(held))
        z = np.load(ROOT/'effects'/store/'delta.npy',mmap_mode='r')
        dest = np.full((len(targets),z.shape[1]),np.nan,dtype=np.float32)
        if len(ids):
            sub = p.iloc[ids]
            val = z[ids]
            np_ = sub.n_cells_pert.to_numpy(float)
            nc = sub.n_cells_control.to_numpy(float)
            w = 1/(1/np.maximum(np_,1e-9)+1/np.maximum(nc,1e-9))
            for j,t in enumerate(canonical):
                own = sub.canonical.to_numpy()==t
                if not own.any():
                    continue
                y = val[own]
                ok = np.isfinite(y)
                weights = ok*w[own,None]
                den = weights.sum(0)
                num = (np.nan_to_num(y)*weights).sum(0)
                dest[j] = np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
        stack.append(dest)
    a = np.stack(stack)
    den = np.isfinite(a).sum(0)
    return np.divide(np.nansum(a,axis=0),den,out=np.zeros(a.shape[1:],dtype=np.float32),
                     where=den>0).astype(np.float32)


def legacy_b16(mean, *, excluded_bios, held_targets):
    """Eligible only if no baseline-training biology or gene is held out."""
    banned = {canonical_bio(x) for x in excluded_bios}
    if held_targets or banned & {'HCT116','HEK293T','CD4','FENG'}:
        raise ValueError('The archived b16 operator is ineligible for this fitting scope')
    from apply_cross_source_residual import transform
    with np.load(ROOT/'processed/cross_source_residual/crsres128/model.npz') as model:
        return transform(mean,model)[0].astype(np.float32)
