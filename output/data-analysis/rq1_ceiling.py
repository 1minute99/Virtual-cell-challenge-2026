"""RQ1 Task 13 -- methods vs ceilings and the H1 verdict (library + CLI).

Plan: input/ideas/zero-shot-limits-context-complexity-plan.md -- H1 + its falsification criterion, operational
definition 6, Task 13; log decision 2026-11-05: H1 is evaluated AS REGISTERED on the H1 panel with the registered
r_max = sqrt[(Var(a) + Var(b)) / Var_rep] (Task 12: 0.423 [0.395, 0.462]) together with the plan's validity check;
r_proj (best-rescaled oracle) and the per-panel ill-posedness are EXPLORATORY.

ESTIMATOR (same family as Task 12, rq1_components.py): primary space = mu-centred response (delta - mu) over all
target x gene elements, G = 6,009 genes, the Task 12 'pack' target set of each panel, the 5 Task 12 split seeds of
the destination truth (dA, dB) and of the other lines' additive profile q (qA, qB). mu = mean_t q. For a predictor
P with halves P^A, P^B (a deterministic prediction has P^A = P^B = P):
   S(x, y) = sum_t w_t <x_t^A - mu^A, y_t^B - mu^B>   (symmetrised over A<->B, Gram matrices averaged over seeds)
   Var_rep = S(d, d)/N,  C_P = S(P, d)/N,  V_P = S(P, P)/N,  r_P = C_P / sqrt(V_P Var_rep)      (N = T G)
   -> the truth's measurement noise is removed (cross-half products); a deterministic prediction keeps its own
      noise (it is part of the method's error), so V_P = ||P - mu||^2 (the mu^A x mu^B term removes mu's noise only).
 Oracle-transfer predictor o = mu + a + b (Task 12; b = mean_t delta - mu):  ||o - mu||^2 = S(q,q) + Q(d,d),
   <o - mu, y - mu> = S(q, y) + Q(d, y), with Q(x, y) = T <xbar^A - mu^A, ybar^B - mu^B>.
 Definition 6: attained fraction = r_ours / r_max (r_ours = frozen pipeline M1 + M2 + M3, Task 14 `<panel>_m3.npz`);
   remaining-error share in the transferable components
   share = E_tr / (E_tr + Var(i)),  E_tr = ||o - P||^2 / N (noise-corrected: ||o||^2 from cross halves), Var(i) = Task 12.
 Uncertainty: target bootstrap, B = 2,000, the SAME multinomial draws as Task 12 (seed 20261113, chunked identically),
   every centring recomputed per draw, other-line fit held fixed; all ceilings and method r share the draws (paired).

PREDICTORS (rows of the CSV):
   m3        frozen pipeline (Task 14 `<panel>_m3.npz`, `delta` = log1p CP10K(fm e^L) - log1p CP10K(fm)) -- r_ours
   prem3     the same before M3 (`<panel>_prem3.npz`)
   xmean     plain cross-line mean transfer: unweighted mean over the pack lines that measure target t (>= 20 cells) of
             their full-depth truth response (Task 9 / Task 12 caches, same readers) -- deterministic
   xmean_nf  diagnostic: the same with its own source noise removed (its halves A/B enter the cross products)
   m3_fmsplit diagnostic: the pipeline's delta re-emitted with the control mean of the SAME truth control half
             (fm^A for P^A, fm^B for P^B) -> removes the shared-control dependence of the prediction (Task 14 note)
   oracle    diagnostic: mu + a + b (= Task 12 r_oracle; not a zero-shot method, it uses the destination's b)
   architectures: the arch_screen_2026-09-22 prediction arrays (K562/HepG2) are NOT on disk (only score JSONs remain);
             even if restored they are on another normalisation (CP50000, official 18,533 axis) against disjoint
             replacement controls -> reported as "not available / not comparable".

ASSUMPTIONS (stated before any method r was computed; checked in code where possible):
  C1 Truth, genes, targets, halves, mu, q, b and Var(i) are exactly Task 12's (pack set); parity: my Var_rep and
     r_oracle reproduce Task 12's point estimates (rel 1e-6) and the bootstrap reproduces its r_max CI.
  C2 Predictions are aligned by panel label (pack `labels` = prediction `targets`) and by gene position
     (`panel_cols`; asserted gmap[panel_cols] == genes_official_idx). KOLF uses the kolf-fixed main build (A3).
  C3 The predictions use ALL destination controls, the truth halves use half of them. The resulting shared control
     noise is a second-order effect (prediction depends on fm through the log1p normalisation); quantified by the
     m3_fmsplit diagnostic on every panel. Budget / ew cpm control dependence (<= 2 % of ||L||, Task 14) is ignored.
  C4 Validity check ("any method's noise-corrected r exceeds r_max with a 95 % CI excluding r_max"): PRIMARY reading =
     the method's own 95 % bootstrap CI lies entirely above the r_max point estimate; SENSITIVITY = the paired
     bootstrap CI of r_method - r_max lies entirely above 0. "Method" = the zero-shot predictors m3, prem3, xmean
     (and the architectures, unavailable). The diagnostics (xmean_nf, m3_fmsplit, oracle) are reported next to it;
     the oracle is the ceiling's own predictor, not a method -- if the main session counts it, H1 becomes untestable
     whenever r_oracle's CI excludes r_max (flagged in the JSON).
  C5 The H1 verdict uses POINT estimates of the attained fraction (threshold 0.7) and share (threshold 0.5), as the
     rule is written without a CI; CIs are reported. Sensitivity: the share with P rescaled by its optimal scalar
     alpha = C_P / V_P (the correlation is scale-free, the squared error is not).
  C6 The real-scorer runs (H1, PDS only) use the Task 17 emission (rq3_corrections.emit_share, pinned sums
     round(fm e^L n), n = 300 control cells rng 0) and command line. The oracle package: on the 6,009 Task 12 genes the
     normalised profile is set to expm1(o + log1p CP10K(fm)) (negatives clipped to 0); every other gene keeps the
     pipeline's profile; renormalised to 10^4 and given the pipeline's per-target library; Od dropped (all via L).
     They do not change the H1 verdict.
  C7 (added after the first run showed a NEGATIVE xmean share; before any verdict was written) Source-noise overlap:
     q (the oracle's a) and every source-built prediction use the SAME source cells, so a deterministic prediction's
     cross products with q^A / q^B contain shared source noise. This biases E_tr (hence the share) DOWNWARD and the
     mu-centring of r slightly (1/T). For xmean it is removed exactly: its cross terms with q use its own halves
     (twin 'xmean_nf'); 'xmean_naive' shows the uncorrected value. For m3 / prem3 no half-built prediction exists:
     the registered share is reported as computed, with a bound share_adj = E_tr + kappa x bias(xmean), kappa = the
     LS loading of m3 on the noisy xmean (kappa_ls) or the norm ratio (kappa_amp).
Seeds: bootstrap 20261113 (= Task 12); truth halves Task 12; emission rng 0 (Task 17); real-run paired bootstrap 20261115.
CLI (nice -n 10; OMP/MKL 4 threads):
   ct_transfer: python rq1_ceiling.py analyse [--panels ...] [--boot 2000]  -> results/rq1_ceiling.csv, rq1_h1.json, png
   vcc2026:     python -P rq1_ceiling.py build | real ARM | collect          (queued, under results/.celleval.lock)
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
from pathlib import Path
import argparse, json, sys, time
import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC')
DA = ROOT / 'output/data-analysis'
RESULTS = DA / 'results'
PRED = RESULTS / 'preds'
T12 = RESULTS / 'rq1_components'
OUT = RESULTS / 'rq1_ceiling'
BOOT_SEED = 20261113
REAL_BOOT_SEED = 20261115
N_BOOT = 2000
PANELS7 = ["H1", "OR_HCT116X", "KOLF", "RN_K562", "RN_HepG2", "RN_Jurkat", "RN_RPE1"]
METHODS = ['m3', 'prem3', 'xmean']                 # zero-shot methods (validity check, C4)
DIAG = ['xmean_nf', 'xmean_naive', 'm3_fmsplit']   # diagnostics (xmean_naive = xmean without the C7 twin correction)
ARCH_NOTE = ('not available / not comparable: vcc2026_data/reports/arch_screen_2026-09-22{,_HepG2} keep only score '
             'JSONs (the referenced <arm>.npy arrays are deleted); targets match the RN_K562 / RN_HepG2 panels (300/300) '
             'but the arrays were on the official 18,533-gene axis in log1p CP50000 against disjoint replacement '
             'controls (simple_residual independent_controls_1944), not the harness truth controls')
sys.path.insert(0, str(DA))


def log(**k):
    print(json.dumps(k, default=str), flush=True)


def delta_lcp(E, fm):
    """Task 14 deliverable delta for expected counts fm * E (rows = targets)."""
    S = fm[None, :] * E
    return np.log1p(S / S.sum(1, keepdims=True) * 1e4) - np.log1p(fm / fm.sum() * 1e4)[None, :]


def expL(z):
    E = np.exp(z['lfc'].astype(np.float64)); Od = z['Od']; j = np.isfinite(Od); E[j] = np.exp(Od[j].astype(np.float64))
    return E


# =========================================================================== inputs
def load_inputs(panel, fmsplit=True):
    import rq1_components as RC
    pk = dict(np.load(T12 / f'{panel}_pack.npz', allow_pickle=False))
    labels = [str(x) for x in pk['labels']]; cols = pk['panel_cols']; T = len(labels)
    preds, checks = {}, {}
    for arm in ('m3', 'prem3'):
        z = np.load(PRED / f'{panel}_{arm}.npz', allow_pickle=True)
        if panel == 'KOLF':
            assert 'rebuild' in json.loads(str(z['meta']))['kolf_fix'], 'KOLF must be the kolf-fixed build (A3)'
        tix = {str(t): i for i, t in enumerate(z['targets'])}
        rows = np.array([tix[l] for l in labels])
        assert np.array_equal(z['gmap'][cols], pk['genes_official_idx']), 'gene alignment'
        preds[arm] = z['delta'][rows][:, cols].astype(np.float64)
        if arm == 'm3':
            E = expL(z)
            fm = z['fm'].astype(np.float64)
            rec = delta_lcp(E, fm)
            checks['m3_delta_recompute_maxabs'] = float(np.abs(rec - z['delta']).max())
            genes_pred = np.asarray(z['genes']).astype(str); m3rows = rows
    # ---- cross-line mean transfer (pack lines), full depth and halves
    G_off, ccols = RC.gene_set()
    assert np.array_equal(G_off, pk['genes_official_idx'])
    srcs = RC.other_sets(panel)['pack']
    targets = [str(t) for t in pk['targets']]
    G = len(G_off); S = pk['dA'].shape[0]
    xf = np.zeros((T, G)); xA = np.zeros((S, T, G)); xB = np.zeros((S, T, G)); cnt = np.zeros(T)
    per_line = {}
    for s in srcs:
        resp, _ = RC.source_responses(s, targets, ccols,
                                      drop_donors=RC.FENG_DROP.get(panel, ()) if s == 'FENG' else ())
        per_line[s] = resp
        for j, t in enumerate(targets):
            if t in resp:
                xf[j] += resp[t]['full']; xA[:, j] += resp[t]['A']; xB[:, j] += resp[t]['B']; cnt[j] += 1
    assert cnt.min() >= 2, 'every pack target needs >= 2 other lines'
    xf /= cnt[:, None]; xA /= cnt[None, :, None]; xB /= cnt[None, :, None]
    # parity: rebuild q (additive fit) from the same responses and compare with Task 12's q
    used = sorted({s for s in srcs if any(t in per_line[s] for t in targets)})
    rws = [(used.index(s), j) for j, t in enumerate(targets) for s in srcs if s in used and t in per_line[s]]
    li = np.array([r[0] for r in rws]); ti = np.array([r[1] for r in rws])
    Pinv = RC.additive_design(li, ti, T, len(used))
    Yf = np.stack([per_line[used[k]][targets[j]]['full'] for k, j in rws])
    YA0 = np.stack([per_line[used[k]][targets[j]]['A'][0] for k, j in rws])
    checks['q_full_parity_maxabs'] = float(np.abs(RC.target_profile(Yf, Pinv, T, len(used)) - pk['q_full']).max())
    checks['qA0_parity_maxabs'] = float(np.abs(RC.target_profile(YA0, Pinv, T, len(used)) - pk['qA'][0]).max())
    checks['xmean_lines_per_target_median'] = float(np.median(cnt))
    del per_line, Yf, YA0
    halves = {'m3': (preds['m3'], preds['m3']), 'prem3': (preds['prem3'], preds['prem3']),
              'xmean': (xf, xf), 'xmean_nf': (xA, xB)}
    # ---- fm-split diagnostic: pipeline delta re-emitted with the truth's own control halves
    if fmsplit:
        import common as CM
        a = CM.open_panel(panel)
        try:
            lab = a.obs['target'].astype(str).to_numpy()
            assert list(a.var_names.astype(str)) == list(genes_pred), 'prediction gene axis = truth var_names'
            tg = sorted(set(lab) - {'non-targeting'}); g = len(tg)
            ii = np.flatnonzero(lab == 'non-targeting')
            X = CM._rows(a, ii).astype(np.float64)
        finally:
            a.file.close()
        checks['ctrl_mean_vs_fm_maxrel'] = float(np.abs(np.asarray(X.mean(0)).ravel() - fm).max() / fm.max())
        pid = RC.PANELS7.index(panel)
        pA = np.zeros((S, T, G)); pB = np.zeros((S, T, G))
        for s in range(S):
            perm = np.random.default_rng([RC.DEST_BASE, s, pid, g]).permutation(len(ii))
            inA = perm < len(ii) // 2
            fA = np.asarray(X[inA].mean(0)).ravel(); fB = np.asarray(X[~inA].mean(0)).ravel()
            pA[s] = delta_lcp(E, fA)[m3rows][:, cols]; pB[s] = delta_lcp(E, fB)[m3rows][:, cols]
        halves['m3_fmsplit'] = (pA, pB)
        del X
    return pk, halves, checks


# =========================================================================== Gram algebra
class Grams:
    """Symmetrised cross-half Gram of stacked rows [d; q; P_1; ...] averaged over split seeds."""

    def __init__(self, pk, halves):
        self.names = ['d', 'q'] + list(halves)
        S, T, G = pk['dA'].shape
        self.T, self.G, self.N = T, G, T * G
        k = len(self.names)
        M = np.zeros((k * T, k * T))
        for s in range(S):
            RA = [pk['dA'][s], pk['qA'][s]]; RB = [pk['dB'][s], pk['qB'][s]]
            for nm in halves:
                a, b = halves[nm]
                RA.append(a[s] if a.ndim == 3 else a); RB.append(b[s] if b.ndim == 3 else b)
            RA = np.vstack(RA).astype(np.float64); RB = np.vstack(RB).astype(np.float64)
            AB = RA @ RB.T
            M += 0.5 * (AB + AB.T)
        self.M = M / S
        self.ix = {n: slice(i * T, (i + 1) * T) for i, n in enumerate(self.names)}

    TWIN = {'xmean': 'xmean_nf'}

    def predictors(self):
        out = []
        for p in self.names[2:]:
            out.append((p, p, self.TWIN.get(p, p)))
            if p in self.TWIN:
                out.append((p + '_naive', p, p))
        return out

    def blk(self, x, y):
        return self.M[self.ix[x], self.ix[y]]

    def stats(self, W):
        W = np.atleast_2d(W); T, N = self.T, self.N
        dg, qf = {}, {}

        def D(x, y):
            if (x, y) not in dg:
                dg[(x, y)] = W @ np.diag(self.blk(x, y))
            return dg[(x, y)]

        def F(x, y):
            key = tuple(sorted((x, y)))
            if key not in qf:
                qf[key] = np.einsum('bt,ts,bs->b', W, self.blk(*key), W, optimize=True)
            return qf[key]

        def Sxy(x, y):
            return D(x, y) - (F(x, 'q') + F('q', y) - F('q', 'q')) / T

        def Qxy(x, y):
            return (F(x, y) - F(x, 'q') - F('q', y) + F('q', 'q')) / T

        Vrep = Sxy('d', 'd') / N
        Vo = (Sxy('q', 'q') + Qxy('d', 'd')) / N
        Co = (Sxy('q', 'd') + Qxy('d', 'd')) / N
        out = dict(Var_rep=Vrep, oracle_V=Vo, oracle_C=Co, oracle_r=Co / np.sqrt(Vo * Vrep))
        for outname, p, tw in self.predictors():
            # cross terms with q / mu use the source-noise-free twin (if any): a deterministic prediction built from
            # the SAME source cells as q would otherwise share their noise (C7)
            def Sx(x, y, tw=tw, p=p):
                xx = tw if (x == p and y == 'q') else x
                yy = tw if (y == p and x == 'q') else y
                if x == p and y == p:            # S(P,P) = D(P,P) - (2 F(P,q) - F(q,q))/T
                    return D(p, p) - (2 * F(tw, 'q') - F('q', 'q')) / T
                Dx = D(xx, yy)
                Fxq = F(tw if x == p else x, 'q'); Fqy = F('q', tw if y == p else y)
                return Dx - (Fxq + Fqy - F('q', 'q')) / T

            def Qx(x, y, tw=tw, p=p):
                Fxq = F(tw if x == p else x, 'q'); Fqy = F('q', tw if y == p else y)
                return (F(x, y) - Fxq - Fqy + F('q', 'q')) / T
            C = Sx(p, 'd') / N; V = Sx(p, p) / N; Op = (Sx('q', p) + Qx('d', p)) / N
            with np.errstate(invalid='ignore', divide='ignore'):
                r = C / np.sqrt(V * Vrep)
                al = C / V
            p = outname
            out.update({f'{p}|C': C, f'{p}|V': V, f'{p}|r': r, f'{p}|E_tr': Vo - 2 * Op + V,
                        f'{p}|E_tr_scaled': Vo - 2 * al * Op + al ** 2 * V, f'{p}|alpha': al,
                        f'{p}|MSE_nf': (Sxy('d', 'd') - 2 * C * N + V * N) / N})
        # C7 sensitivity: the pipeline is built from the same source cells as q, so its cross terms with q carry shared
        # source noise that cannot be removed without half-built predictions. Bound: m3's bias = kappa x xmean's bias
        # (xmean's is estimable from its halves); kappa_ls = LS loading of m3 on the noisy xmean, kappa_amp = norm ratio.
        if 'xmean' in self.names and 'm3' in self.names:
            Sdet = lambda x, y: D(x, y) - (F(x, 'q') + F('q', y) - F('q', 'q')) / T
            bias_x = out['xmean|E_tr'] - out['xmean_naive|E_tr']
            Vx = Sdet('xmean', 'xmean')
            for p in ('m3', 'prem3'):
                k_ls = Sdet(p, 'xmean') / Vx
                k_amp = np.sqrt(Sdet(p, p) / Vx)
                out[f'{p}|kappa_ls'] = k_ls; out[f'{p}|kappa_amp'] = k_amp
                out[f'{p}|E_tr_adj_ls'] = out[f'{p}|E_tr'] + k_ls * bias_x
                out[f'{p}|E_tr_adj_amp'] = out[f'{p}|E_tr'] + k_amp * bias_x
            out['xmean|E_tr_bias'] = bias_x
        return out


def derive(st, comp):
    """Add ceiling-relative quantities; st and comp are dicts of (B,) arrays for the same draws."""
    o = {}
    rmax, rproj, Vi = comp['r_max'], comp['r_proj'], comp['Var_i']
    o.update({f'ceil|{k}': comp[k] for k in ('r_max', 'r_proj', 'r_oracle', 'r_max_orth', 'ceiling_meas', 'Var_a',
                                             'Var_b', 'Var_i', 'C_ai', 'Var_rep')})
    o['ceil|oracle_r_recomputed'] = st['oracle_r']
    o['ceil|Var_rep_recomputed'] = st['Var_rep']
    preds = sorted({k.split('|')[0] for k in st if '|' in k})
    rows = preds + ['oracle']
    for p in rows:
        r = st['oracle_r'] if p == 'oracle' else st[f'{p}|r']
        with np.errstate(invalid='ignore', divide='ignore'):
            o[f'{p}|r'] = r
            o[f'{p}|frac_rmax'] = r / rmax
            o[f'{p}|diff_rmax'] = r - rmax
            o[f'{p}|frac_rproj'] = r / rproj
            o[f'{p}|diff_rproj'] = r - rproj
            o[f'{p}|r_vs_noisy_truth'] = r * comp['ceiling_meas']
            if p != 'oracle':
                E, Es = st[f'{p}|E_tr'], st[f'{p}|E_tr_scaled']
                o[f'{p}|E_tr'] = E
                o[f'{p}|share'] = E / (E + Vi)
                o[f'{p}|share_scaled'] = Es / (Es + Vi)
                o[f'{p}|alpha'] = st[f'{p}|alpha']
                o[f'{p}|MSE_nf'] = st[f'{p}|MSE_nf']
                o[f'{p}|cross_term'] = st[f'{p}|MSE_nf'] - E - Vi
                for tag in ('adj_ls', 'adj_amp'):
                    if f'{p}|E_tr_{tag}' in st:
                        Ea = st[f'{p}|E_tr_{tag}']
                        o[f'{p}|share_{tag}'] = Ea / (Ea + Vi)
                        o[f'{p}|kappa_{tag[4:]}'] = st[f'{p}|kappa_{tag[4:]}']
    return o


def boot_weights(T, n_boot, seed=BOOT_SEED, chunk=250):
    rng = np.random.default_rng(seed); out = []
    for lo in range(0, n_boot, chunk):
        b = min(chunk, n_boot - lo)
        out.append(np.stack([np.bincount(rng.integers(0, T, T), minlength=T) for _ in range(b)]).astype(np.float64))
    return np.vstack(out)


def analyse_panel(panel, n_boot):
    import rq1_components as RC
    t0 = time.time()
    pk, halves, checks = load_inputs(panel)
    gr = RC.mean_grams([RC.grams(pk['dA'][s], pk['dB'][s], pk['qA'][s], pk['qB'][s]) for s in range(pk['dA'].shape[0])])
    GM = Grams(pk, halves)
    del halves
    pt = derive(GM.stats(np.ones(GM.T)), RC.components(gr))
    pt = {k: float(v[0]) for k, v in pt.items()}
    W = boot_weights(GM.T, n_boot)
    acc = {}
    for lo in range(0, n_boot, 250):
        Wc = W[lo:lo + 250]
        d = derive(GM.stats(Wc), RC.components(gr, Wc))
        for k, v in d.items():
            acc.setdefault(k, []).append(v)
    acc = {k: np.concatenate(v) for k, v in acc.items()}
    ci = {}
    for k, v in acc.items():
        f = v[np.isfinite(v)]
        ci[k] = (float(np.quantile(f, 0.025)), float(np.quantile(f, 0.975)), int(len(f))) if len(f) else (np.nan, np.nan, 0)
    # parity with Task 12
    t12 = json.loads((RESULTS / 'rq1_components.json').read_text())['panels'][panel]['pack']
    checks['Var_rep_vs_T12_rel'] = abs(pt['ceil|Var_rep_recomputed'] / t12['point']['Var_rep'] - 1)
    ro = t12['point']['r_oracle']
    checks['r_oracle_vs_T12_abs'] = abs(pt['ceil|oracle_r_recomputed'] - ro) if np.isfinite(ro) else None
    lo12, hi12 = t12['ci']['r_max']
    checks['r_max_ci_vs_T12_abs'] = (max(abs(ci['ceil|r_max'][0] - lo12), abs(ci['ceil|r_max'][1] - hi12))
                                     if np.isfinite(lo12) else None)
    checks['n_targets'] = GM.T; checks['n_genes'] = GM.G; checks['secs'] = round(time.time() - t0, 1)
    log(event='panel_done', panel=panel, **{k: v for k, v in checks.items()},
        r_m3=pt['m3|r'], r_max=pt['ceil|r_max'], share_m3=pt['m3|share'])
    return dict(point=pt, ci=ci, checks=checks, n_targets=GM.T)


def validity(res, ceiling):
    """C4 check for one panel against 'r_max' or 'r_proj'."""
    pt, ci = res['point'], res['ci']
    c = pt[f'ceil|{ceiling}']
    out = {}
    for p in METHODS + DIAG + ['oracle']:
        lo, hi, _ = ci[f'{p}|r']; dlo, dhi, _ = ci[f'{p}|diff_{"rmax" if ceiling == "r_max" else "rproj"}']
        out[p] = dict(r=pt[f'{p}|r'], ci=[lo, hi], ceiling=c,
                      exceeds_primary=bool(np.isfinite(c) and lo > c),
                      paired_diff_ci=[dlo, dhi], exceeds_paired=bool(np.isfinite(dlo) and dlo > 0),
                      role='method' if p in METHODS else 'diagnostic')
    meth = [p for p in METHODS if out[p]['exceeds_primary']]
    return dict(per_predictor=out, ceiling_estimable=bool(np.isfinite(c)),
                any_method_exceeds_primary=bool(meth), methods_exceeding_primary=meth,
                any_method_exceeds_paired=any(out[p]['exceeds_paired'] for p in METHODS),
                oracle_exceeds_primary=out['oracle']['exceeds_primary'])


def verdict(res, ceiling='r_max'):
    pt = res['point']; key = 'frac_rmax' if ceiling == 'r_max' else 'frac_rproj'
    frac, share = pt[f'm3|{key}'], pt['m3|share']
    v = validity(res, ceiling)
    if not v['ceiling_estimable'] or not np.isfinite(frac) or not np.isfinite(share) or pt['ceil|Var_i'] <= 0:
        return dict(verdict='untestable (ceiling not estimable)', frac=frac, share=share, validity=v)
    if v['any_method_exceeds_primary']:
        vd = 'untestable (ceiling estimator invalid: a method r exceeds it)'
    elif frac < 0.7 and share >= 0.5:
        vd = 'rejected ("model-limited")'
    elif frac >= 0.7 and share < 0.5:
        vd = 'confirmed'
    else:
        vd = 'inconclusive'
    alt = None
    if v['oracle_exceeds_primary'] and not v['any_method_exceeds_primary']:
        alt = 'untestable if the oracle predictor itself is counted as a "method" (its CI excludes the ceiling)'
    return dict(verdict=vd, frac=frac, share=share, frac_ge_0p7=bool(frac >= 0.7), share_lt_0p5=bool(share < 0.5),
                share_scaled=pt['m3|share_scaled'], validity=v, if_oracle_counted=alt)


def cmd_analyse(panels, n_boot):
    import pandas as pd
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    allres = {}
    for p in panels:
        allres[p] = analyse_panel(p, n_boot)
        (OUT / f'{p}_point_ci.json').write_text(json.dumps(allres[p], indent=1, default=float))
    rows = []
    for p, r in allres.items():
        for k, v in r['point'].items():
            pred, q = k.split('|', 1)
            lo, hi, nb = r['ci'][k]
            rows.append(dict(panel=p, predictor=pred, quantity=q, estimate=v, ci_lo=lo, ci_hi=hi, n_boot_finite=nb,
                             n_targets=r['n_targets'], role=('method' if pred in METHODS else
                                                             'ceiling' if pred == 'ceil' else 'diagnostic'),
                             registered=(p == 'H1')))
        for arch in ('architectures (17, arch_screen_2026-09-22)',):
            rows.append(dict(panel=p, predictor=arch, quantity='r', estimate=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             n_boot_finite=0, n_targets=r['n_targets'], role='method', registered=(p == 'H1'),
                             note=ARCH_NOTE if p in ('RN_K562', 'RN_HepG2') else 'no predictions for this panel'))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / 'rq1_ceiling.csv', index=False)
    summ = {}
    for p, r in allres.items():
        pt, ci = r['point'], r['ci']
        summ[p] = dict(
            n_targets=r['n_targets'], checks=r['checks'],
            ceilings={k: dict(est=pt[f'ceil|{k}'], ci=ci[f'ceil|{k}'][:2]) for k in
                      ('r_max', 'r_proj', 'r_oracle', 'ceiling_meas', 'Var_i')},
            predictors={pr: {q: dict(est=pt[f'{pr}|{q}'], ci=ci[f'{pr}|{q}'][:2]) for q in
                             ('r', 'frac_rmax', 'frac_rproj', 'share', 'share_scaled', 'share_adj_ls', 'share_adj_amp', 'kappa_ls', 'kappa_amp', 'alpha', 'E_tr', 'MSE_nf')
                             if f'{pr}|{q}' in pt} for pr in METHODS + DIAG + ['oracle']},
            registered_rule_r_max=verdict(r, 'r_max'),
            exploratory_r_proj=verdict(r, 'r_proj'),
            note=('H1 panel = registered test' if p == 'H1' else 'secondary panel (descriptive; H1 is registered on H1 only)'))
    h1 = summ.get('H1')
    out = dict(created=time.strftime('%Y-%m-%d %H:%M'), task=13,
               plan='input/ideas/zero-shot-limits-context-complexity-plan.md',
               decision='log 2026-11-05: H1 evaluated as registered on the H1 panel with r_max; r_proj exploratory',
               n_boot=n_boot, boot_seed=BOOT_SEED,
               H1=None if h1 is None else dict(
                   r_max=h1['ceilings']['r_max'], r_ours_m3=h1['predictors']['m3']['r'],
                   attained_fraction=h1['predictors']['m3']['frac_rmax'],
                   remaining_error_share=h1['predictors']['m3']['share'],
                   remaining_error_share_scaled_sensitivity=h1['predictors']['m3']['share_scaled'],
                   remaining_error_share_source_noise_bound=dict(
                       ls=h1['predictors']['m3']['share_adj_ls'], amp=h1['predictors']['m3']['share_adj_amp']),
                   verdict=h1['registered_rule_r_max']['verdict'],
                   validity=h1['registered_rule_r_max']['validity'],
                   if_oracle_counted=h1['registered_rule_r_max']['if_oracle_counted'],
                   exploratory_r_proj=dict(r_proj=h1['ceilings']['r_proj'],
                                           attained_fraction=h1['predictors']['m3']['frac_rproj'],
                                           verdict=h1['exploratory_r_proj']['verdict'],
                                           validity=h1['exploratory_r_proj']['validity'])),
               architectures=ARCH_NOTE, panels=summ, real_scorer='queued (see collect)',
               runtime_s=round(time.time() - t0, 1))
    prev = RESULTS / 'rq1_h1.json'
    if prev.exists():
        old = json.loads(prev.read_text())
        if isinstance(old.get('real_scorer'), dict):
            out['real_scorer'] = old['real_scorer']
    prev.write_text(json.dumps(out, indent=1, default=float))
    plot(summ)
    log(event='done', secs=round(time.time() - t0, 1))


def plot(summ):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    panels = [p for p in PANELS7 if p in summ]
    ser = [('m3', 'pipeline M1+M2+M3', '#2a78d6'), ('prem3', 'pre-M3', '#eb6834'),
           ('xmean', 'cross-line mean transfer', '#1baf7a'), ('oracle', 'oracle mu+a+b (diagnostic)', '#eda100')]
    fig, ax = plt.subplots(figsize=(9, 4.4), dpi=150)
    off = np.linspace(-0.27, 0.27, len(ser))
    for i, p in enumerate(panels):
        c = summ[p]['ceilings']
        for key, ls, col in (('r_max', '-', '#52514e'), ('r_proj', (0, (3, 2)), '#e87ba4')):
            v = c[key]['est']
            if np.isfinite(v):
                ax.hlines(v, i - 0.42, i + 0.42, colors=col, lw=2, linestyles=ls, zorder=1)
    for (k, lab, col), o in zip(ser, off):
        y = np.array([summ[p]['predictors'][k]['r']['est'] for p in panels])
        lo = np.array([summ[p]['predictors'][k]['r']['ci'][0] for p in panels])
        hi = np.array([summ[p]['predictors'][k]['r']['ci'][1] for p in panels])
        x = np.arange(len(panels)) + o
        ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt='o', ms=6, color=col, elinewidth=1.5, capsize=0, label=lab,
                    mec='white', mew=1, zorder=3)
    ax.plot([], [], color='#52514e', lw=2, label='r_max (registered)')
    ax.plot([], [], color='#e87ba4', lw=2, ls=(0, (3, 2)), label='r_proj (exploratory)')
    ax.set_xticks(np.arange(len(panels))); ax.set_xticklabels(panels, fontsize=9); ax.set_xlim(-0.5, len(panels) - 0.5)
    ax.set_ylabel('noise-corrected r (mu-centred)', fontsize=9); ax.set_ylim(0, 1.0)
    ax.grid(axis='y', color='#e6e5e1', lw=0.8); ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=7.5, frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.09))
    ax.set_title('Task 13: methods vs transfer ceilings (95 % target-bootstrap CI; OR r_max not estimable)',
                 fontsize=10, loc='left')
    fig.tight_layout()
    fig.savefig(RESULTS / 'rq1_ceiling.png'); fig.savefig(RESULTS / 'rq1_ceiling.svg')


# =========================================================================== real scorer (H1; vcc2026 env)
def cmd_build():
    z = np.load(PRED / 'H1_m3.npz', allow_pickle=True)
    pk = np.load(T12 / 'H1_pack.npz', allow_pickle=False)
    tg = [str(t) for t in z['targets']]
    assert [str(x) for x in pk['labels']] == tg
    cols = pk['panel_cols']; assert np.array_equal(z['gmap'][cols], pk['genes_official_idx'])
    fm = z['fm'].astype(np.float64)
    E = expL(z); F = fm[None, :] * E; tot = F.sum(1, keepdims=True)
    y = F / tot * 1e4
    c = fm / fm.sum() * 1e4
    o = pk['q_full'].astype(np.float64) + pk['b'].astype(np.float64)[None, :]
    yo = np.expm1(o + np.log1p(c[cols])[None, :])
    n_clip = int((yo < 0).sum()); yo = np.maximum(yo, 0)
    y2 = y.copy(); y2[:, cols] = yo
    y2 *= 1e4 / y2.sum(1, keepdims=True)
    F2 = y2 / 1e4 * tot
    with np.errstate(divide='ignore', invalid='ignore'):
        L2 = np.where(fm[None, :] > 0, np.log(np.maximum(F2, 1e-300) / np.where(fm > 0, fm, 1)[None, :]), 0.0)
    Od2 = np.full_like(z['Od'], np.nan)
    d2 = delta_lcp(np.exp(L2), fm)[:, cols]
    cos = float((d2 * o).sum() / np.sqrt((d2 ** 2).sum() * (o ** 2).sum()))
    # same-gene check against the pipeline's delta outside the oracle genes
    other = np.setdiff1d(np.arange(len(fm)), cols)
    dp = delta_lcp(E, fm)
    d2all = delta_lcp(np.exp(L2), fm)
    chk = dict(oracle_genes=int(len(cols)), clipped_negative_profile_entries=n_clip,
               cos_emitted_delta_vs_oracle=cos, maxabs_emitted_vs_oracle=float(np.abs(d2 - o).max()),
               median_abs_shift_non_oracle_genes=float(np.median(np.abs(d2all[:, other] - dp[:, other]))),
               lfc_pipeline_vs_T17_M3_dest_maxabs=float(np.abs(
                   np.load(RESULTS / 'rq3_corr/H1_real_arms.npz')['lfc_M3_dest'] - z['lfc']).max()))
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / 'H1_real_arms.npz', targets=np.array(tg), genes=np.asarray(z['genes']).astype(str), fm=z['fm'],
             n_cells=int(z['n_cells']), lfc_pipeline=z['lfc'], Od_pipeline=z['Od'],
             lfc_oracle=L2.astype(np.float32), Od_oracle=Od2)
    (OUT / 'build.json').write_text(json.dumps(chk, indent=1))
    log(event='build', **chk)


def cmd_real(arms):
    sys.path.insert(0, str(DA))
    import rq3_corrections as R3
    CM, _ = R3.imports()
    import anndata as ad, pandas as pd, shutil, subprocess
    from scipy import sparse
    z = np.load(OUT / 'H1_real_arms.npz', allow_pickle=True)
    tg = [str(t) for t in z['targets']]; n = int(z['n_cells'])
    ctrl, genes = CM.controls('H1'); genes = np.asarray(genes).astype(str)
    assert list(genes) == list(z['genes']); truth = CM.truth_path('H1')
    assert set(tg) == set(CM.labels('H1')) - {'non-targeting'}
    fm = np.asarray(ctrl.sum(0)).ravel() / ctrl.shape[0]; assert np.allclose(fm, z['fm'], rtol=1e-5)
    cache = R3.CACHE / 'cache_real_H1'; cache.mkdir(parents=True, exist_ok=True)
    for a in arms:
        out = OUT / 'real' / a
        if (out / 'result.json').exists():
            log(event='cached', arm=a); continue
        out.mkdir(parents=True, exist_ok=True); ta = time.time()
        blocks, Ts = R3.emit_share(ctrl, z[f'lfc_{a}'], z[f'Od_{a}'], n)
        X = sparse.vstack(blocks, format='csr'); del blocks
        pred = ad.AnnData(X=sparse.vstack([X, ctrl], format='csr'),
                          obs=pd.DataFrame({'target': pd.Categorical([*np.repeat(tg, n), *['non-targeting'] * ctrl.shape[0]])},
                                           index=[f'p{i}' for i in range(X.shape[0] + ctrl.shape[0])]),
                          var=pd.DataFrame(index=genes))
        pp = out / 'pred.h5ad'; pred.write_h5ad(pp); del pred, X
        log(event='emitted', arm=a, secs=round(time.time() - ta, 1))
        cmd = [str(R3.CE), 'run', '-ap', str(pp), '-ar', str(truth), '--preset', 'vcc2026', '--pert-col', 'target',
               '--control', 'non-targeting', '--input-type', 'counts', '-o', str(out / 'cell_eval2'), '--cache-real',
               str(cache), '--set', 'de.backend=pdex', '--set', 'num_threads=4', '--set', 'gather_threads=4']
        env = dict(os.environ, NUMBA_NUM_THREADS='4', OMP_NUM_THREADS='4')
        with (out / 'scorer.log').open('w') as g:
            subprocess.run(cmd, stdout=g, stderr=subprocess.STDOUT, env=env, check=True)
        agg = pd.read_csv(out / 'cell_eval2/agg_results.csv', index_col=0)
        raw = {k: float(agg.loc['mean', k]) for k in agg.columns
               if np.isfinite(pd.to_numeric(agg.loc['mean', k], errors='coerce'))}
        rl = pd.read_csv(out / 'cell_eval2/results.csv')
        rl[rl['metric'] == 'pds_cosine'].set_index('perturbation')['value'].to_csv(out / 'pds_per_target.csv')
        shutil.copy(out / 'cell_eval2/agg_results.csv', out / 'agg_results.csv')
        res = dict(arm=a, raw=raw, pds_real=raw.get('pds_cosine'), seconds=round(time.time() - ta, 1),
                   emission='Task 17 emit_share: n=300 random controls (rng 0), pinned sums round(fm exp(L) n); PDS only')
        (out / 'result.json').write_text(json.dumps(res, indent=1)); pp.unlink()
        shutil.rmtree(out / 'cell_eval2', ignore_errors=True)
        log(event='real', arm=a, pds=res['pds_real'], secs=res['seconds'])


def cmd_collect():
    import pandas as pd
    res = {}
    for a in ('pipeline', 'oracle'):
        f = OUT / 'real' / a / 'result.json'
        if f.exists():
            res[a] = json.loads(f.read_text())
    rs = dict(status='pending', arms={a: r['pds_real'] for a, r in res.items()},
              note='PDS with the Task 17 emission; DE metrics/composite not comparable; does not change the H1 verdict')
    t17 = RESULTS / 'rq3_corr/real/M3_dest/result.json'
    if t17.exists() and 'pipeline' in res:
        rs['pipeline_vs_T17_M3_dest'] = res['pipeline']['pds_real'] - json.loads(t17.read_text())['pds_real']
    if len(res) == 2:
        p = pd.read_csv(OUT / 'real/pipeline/pds_per_target.csv', index_col=0).iloc[:, 0]
        o = pd.read_csv(OUT / 'real/oracle/pds_per_target.csv', index_col=0).iloc[:, 0].reindex(p.index)
        d = (o - p).to_numpy(); rng = np.random.default_rng(REAL_BOOT_SEED)
        bs = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT)])
        rs.update(status='done', oracle_minus_pipeline=float(d.mean()),
                  ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], n_targets=int(len(d)))
    f = RESULTS / 'rq1_h1.json'
    j = json.loads(f.read_text()) if f.exists() else {}
    j['real_scorer'] = rs
    f.write_text(json.dumps(j, indent=1, default=float))
    log(event='collect', **rs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['analyse', 'build', 'real', 'collect'])
    ap.add_argument('arms', nargs='*')
    ap.add_argument('--panels', nargs='*', default=PANELS7)
    ap.add_argument('--boot', type=int, default=N_BOOT)
    a = ap.parse_args()
    if a.cmd == 'analyse':
        cmd_analyse(a.panels, a.boot)
    elif a.cmd == 'build':
        cmd_build()
    elif a.cmd == 'real':
        cmd_real(a.arms or ['oracle', 'pipeline'])
    else:
        cmd_collect()


if __name__ == '__main__':
    main()
