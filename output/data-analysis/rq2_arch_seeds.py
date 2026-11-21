"""Task 10 collector (plan input/ideas/zero-shot-limits-context-complexity-plan.md): seed replication of six screen
architectures on K562 + HepG2, official cell-eval2 vcc2026 composite, vs the refit_v2 transfer baseline.

ASSUMPTIONS (stated; checked where checkable):
A1 Composite = `avg_score` (from_replicate column) of `cell-eval2 score` against each line's frozen real bundle;
   six metrics = the six normalised components in the same file (raw six from receipt.json also kept). This is the
   quantity FINDINGS.md ranks on.
A2 Unit of replication = training seed s in {0,1,2}; seed s also sets the cell-reconstruction seed (as in seed 0).
   delta(arch, line, s) = composite(arch_correct, s) - composite(base_refit_v2, s), paired at the same s.
A3 Seeds are independent replicates (separate torch/numpy RNG streams; same data). Lines are NOT exchangeable
   replicates (K562 has no measurable context specificity, HepG2 has); the two-line mean is the pre-specified
   primary quantity (as in FINDINGS / Task 18), per-line values are secondary.
A4 95 % CI = Student t on the n = 3 seed values (df = 2, t = 4.303); normality of seed-to-seed variation is
   ASSUMED and cannot be tested with n = 3 (reported as such). A sign count (k/3 seeds > 0) is given as a
   distribution-free descriptor; with n = 3 an exact sign test cannot reach p < 0.05 (min one-sided p = 0.125).
A5 No multiplicity correction: descriptive consistency check for H2c (plan Task 10), not a hypothesis test.
A6 Reproduction gates: (a) base_refit_v2 seed-0 composite re-scored with the regenerated K562 real panel must equal
   the stored seed-0 value to 1e-6 (else the protocol is not reproduced); (b) additive/moe seed-0 retrained under
   the current environment are compared with the stored seed-0 harness PDS / composite (GPU nondeterminism and the
   current torch version may differ; reported, not gated).
Seed 0 values are read from the existing reports (never rewritten). Missing seeds are reported, not imputed.
"""
import json, math
from pathlib import Path
import numpy as np, pandas as pd

R = Path('/home/wonmin-kim/AIVC/vcc2026_data/reports')
O = Path('/home/wonmin-kim/AIVC/output/data-analysis/results')
A = O / 'rq2_arch'
ARCHS = json.loads((A / 'selection.json').read_text())['archs']
LINES = ['K562', 'HepG2']
SEED0 = {'K562': R / 'arch_screen_2026-09-22', 'HepG2': R / 'arch_screen_2026-09-22_HepG2'}
SIX = ['pds_cosine', 'expr_mse_unbiased_capped_norm', 'de_wilcoxon_lfc_nmae',
       'de_wilcoxon_direction_fidelity_yield_raw', 'de_wilcoxon_direction_reach_raw', 'de_wilcoxon_sig_jaccard']
T975_DF2 = 4.302652729911275


def src(line, seed):
    return SEED0[line] if seed == 0 else A / f's{seed}_{line}'


def read_run(base, line, tag, seed):
    d = base / 'preset' / 'runs' / line / tag / f'seed{seed}'
    f = d / 'scored.csv'
    if not f.exists():
        return None
    t = pd.read_csv(f); col = 'metric' if 'metric' in t.columns else t.columns[0]
    mm = {str(a): float(b) for a, b in zip(t[col], t['from_replicate']) if pd.notna(b)}
    rec = json.loads((d / 'receipt.json').read_text())
    row = dict(composite=mm['avg_score'], **{f'n_{k}': mm.get(k) for k in SIX},
               **{f'raw_{k}': rec['raw_six'].get(k) for k in SIX}, pds_parity_error=rec.get('pds_parity_error'),
               run_dir=str(d))
    return row


def harness(base, arch):
    p = base / 'results.json'
    if not p.exists():
        return {}
    r = json.loads(p.read_text()).get(arch, {})
    return dict(best_step=r.get('best_step'), harness_pds_correct=r.get('pds_correct'),
                harness_pds_correct_a025=r.get('pds_correct_a0.25'), resid_cos_correct=r.get('resid_cos_correct'),
                resid_cos_wrong=r.get('resid_cos_wrong'))


def tci(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float); n = len(v)
    if n == 0:
        return dict(n=0)
    m = float(v.mean()); sd = float(v.std(ddof=1)) if n > 1 else float('nan')
    out = dict(n=n, mean=m, sd=sd, values=[float(x) for x in v], n_pos=int((v > 0).sum()))
    if n == 3:
        h = T975_DF2 * sd / math.sqrt(n); out.update(ci_lo=m - h, ci_hi=m + h, ci_excludes_0=bool(m - h > 0 or m + h < 0),
                                                    beats_baseline_ci=bool(m - h > 0))
    return out


def main():
    rows = []
    for line in LINES:
        for seed in (0, 1, 2):
            b = src(line, seed); base = read_run(b, line, 'base_refit_v2', seed)
            if base is not None:
                rows.append(dict(arch='base_refit_v2', line=line, seed=seed, base_composite=base['composite'],
                                 delta_composite=0.0, **base))
            for arch in ARCHS:
                r = read_run(b, line, f'{arch}_correct', seed)
                if r is None:
                    continue
                rows.append(dict(arch=arch, line=line, seed=seed,
                                 base_composite=None if base is None else base['composite'],
                                 delta_composite=None if base is None else r['composite'] - base['composite'],
                                 **r, **harness(b, arch)))
    df = pd.DataFrame(rows)
    df.to_csv(O / 'rq2_arch_seeds.csv', index=False)

    summ = dict(plan='input/ideas/zero-shot-limits-context-complexity-plan.md', task=10, archs=ARCHS,
                ci='Student t, n = 3 seeds, df = 2 (normality assumed, untestable at n = 3)', per_arch={})
    have = df[df.arch != 'base_refit_v2']
    for arch in ARCHS:
        d = have[have.arch == arch]; e = {}
        for line in LINES:
            e[line] = tci(d[d.line == line].sort_values('seed').delta_composite.tolist())
        piv = d.pivot_table(index='seed', columns='line', values='delta_composite')
        two = piv.dropna().mean(axis=1) if set(LINES) <= set(piv.columns) else pd.Series(dtype=float)
        e['two_line_mean'] = tci(two.tolist()); e['two_line_seeds'] = [int(s) for s in two.index]
        s0 = d[d.seed == 0].set_index('line').delta_composite.to_dict()
        e['seed0_two_line_mean'] = float(np.mean([s0[l] for l in LINES])) if all(l in s0 for l in LINES) else None
        summ['per_arch'][arch] = e
    base = df[df.arch == 'base_refit_v2']
    summ['baseline_composite'] = {l: tci(base[base.line == l].sort_values('seed').composite.tolist()) for l in LINES}
    complete = all(summ['per_arch'][a]['two_line_mean'].get('n') == 3 for a in ARCHS)
    summ['complete'] = complete
    summ['n_runs_new'] = int(((df.seed > 0)).sum())
    beats = {a: [k for k in ('K562', 'HepG2', 'two_line_mean') if summ['per_arch'][a][k].get('beats_baseline_ci')] for a in ARCHS}
    summ['beats_baseline_ci_excl_0'] = {a: v for a, v in beats.items() if v}
    summ['any_learned_beats_baseline_two_line'] = any('two_line_mean' in v for v in beats.values())
    summ['any_learned_beats_baseline_any_line'] = any(beats.values())
    tl = {a: summ['per_arch'][a]['two_line_mean'].get('mean') for a in ARCHS if summ['per_arch'][a]['two_line_mean'].get('n') == 3}
    sl = {(a, l): summ['per_arch'][a][l].get('mean') for a in ARCHS for l in LINES if summ['per_arch'][a][l].get('n') == 3}
    if tl:
        ba = max(tl, key=tl.get); bs = max(sl, key=sl.get)
        summ['task18_best_arch_gain'] = dict(
            primary_two_line_seed_mean=dict(arch=ba, gain=tl[ba], **{k: summ['per_arch'][ba]['two_line_mean'].get(k) for k in ('sd', 'ci_lo', 'ci_hi')}),
            sensitivity_single_line_seed_mean=dict(arch=bs[0], line=bs[1], gain=sl[bs], **{k: summ['per_arch'][bs[0]][bs[1]].get(k) for k in ('sd', 'ci_lo', 'ci_hi')}),
            seed0_values_replaced=dict(primary='additive +0.0036 (two-line mean)', sensitivity='lora_hyper HepG2 +0.0169'),
            caveat='Seed means cover only the plan\'s six replicated architectures; LoRA hypernetwork (seed-0 single-line '
                    'best, HepG2 +0.0169) was not in the six and remains seed-0 only. Scorer setup: K562/HepG2 panels, '
                    'refit_v2 base, alpha 1, ext_preset_score emission (differs from the H1 panel / F2 emission of Task 18).')
    # reproduction gates
    rp = A / 's0repro_K562'; rep = {}
    b0 = read_run(rp, 'K562', 'base_refit_v2', 0); s0 = read_run(SEED0['K562'], 'K562', 'base_refit_v2', 0)
    if b0:
        rep['base_refit_v2_seed0_composite'] = dict(stored=s0['composite'], rescored=b0['composite'],
                                                    abs_diff=abs(b0['composite'] - s0['composite']),
                                                    pass_1e6=abs(b0['composite'] - s0['composite']) <= 1e-6)
    for arch in ('additive', 'moe'):
        r = read_run(rp, 'K562', f'{arch}_correct', 0); st = read_run(SEED0['K562'], 'K562', f'{arch}_correct', 0)
        h_new = harness(rp, arch); h_old = harness(SEED0['K562'], arch)
        if h_new:
            rep[arch] = dict(stored_harness_pds=h_old.get('harness_pds_correct'), new_harness_pds=h_new.get('harness_pds_correct'),
                             stored_best_step=h_old.get('best_step'), new_best_step=h_new.get('best_step'),
                             stored_composite=st['composite'] if st else None, new_composite=r['composite'] if r else None)
    summ['seed0_reproduction'] = rep
    (O / 'rq2_arch_seeds_summary.json').write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps(dict(event='collected', rows=len(df), complete=complete,
                          two_line={a: (round(v, 4) if v is not None else None) for a, v in tl.items()}), indent=None))
    plot(df)


def plot(df):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    d = df[df.arch != 'base_refit_v2'].dropna(subset=['delta_composite'])
    if d.empty:
        return
    col = {'K562': '#2a78d6', 'HepG2': '#eb6834'}
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=150)
    for j, line in enumerate(LINES):
        for i, arch in enumerate(ARCHS):
            v = d[(d.arch == arch) & (d.line == line)].delta_composite.to_numpy()
            x = i + (j - 0.5) * 0.3
            ax.scatter(np.full(len(v), x), v, s=18, color=col[line], alpha=0.55, zorder=3, linewidths=0)
            if len(v) == 3:
                m, h = v.mean(), T975_DF2 * v.std(ddof=1) / math.sqrt(3)
                ax.errorbar(x, m, yerr=h, fmt='o', ms=6, color=col[line], lw=2, capsize=0, zorder=4,
                            label=line if i == 0 else None)
    ax.axhline(0, color='#888888', lw=1, zorder=1)
    ax.set_xticks(range(len(ARCHS))); ax.set_xticklabels([a.replace('_', ' ') for a in ARCHS], rotation=20, ha='right')
    ax.set_ylabel('composite − refit_v2 (paired seed)')
    ax.set_title('Six architectures, seeds 0–2: mean ± t-95% CI (dots = seeds)', fontsize=10)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', color='#e5e5e5', lw=0.6); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(O / 'rq2_arch_seeds.png'); plt.close(fig)


if __name__ == '__main__':
    main()
