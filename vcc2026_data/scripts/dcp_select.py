"""Nested alpha / arm selection, then the final arm map for official preset scoring.

The alpha for line L is chosen by the mean harness PDS of the OTHER four lines, never L's own.
This is a cheap development metric used only to pick which file to score (Part 18 permits that);
the final ranking is the official VCC2026 preset, computed downstream.
"""
from pathlib import Path
import argparse, json, statistics
import numpy as np

ROOT = Path('/home/wonmin-kim/AIVC/vcc2026_data/reports/drug_context_pretrain_2026-09-20')
LINES = ('K562', 'RPE1', 'HepG2', 'Jurkat', 'KOLF')
ALPHAS = (0.25, 0.5, 1.0)
FAMILIES = ('scratch', 'tahoe_frozen', 'tgse_frozen', 'tgse_partial',
            'gse_partial', 'joint', 'joint_ft', 'ENS_plain', 'ENS_gated')
FIXED = ('refit_alpha0',)


def load():
    rows = json.loads((ROOT / 'scores.json').read_text())
    gp = ROOT / 'gate_scores.json'
    if gp.exists():
        rows += json.loads(gp.read_text())
    tab = {}
    for r in rows:
        tab.setdefault((r['line'], r['tag']), []).append(r['pds'])
    return {k: statistics.mean(v) for k, v in tab.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(ROOT / 'arm_map.json'))
    a = ap.parse_args()
    pds = load()
    lines = [l for l in LINES if (l, 'refit_alpha0') in pds]
    arm_map = {l: list(FIXED) for l in lines}
    chosen = {}
    for fam in FAMILIES:
        have = [al for al in ALPHAS if all((l, f'{fam}_a{al:g}') in pds for l in lines)]
        if not have:
            continue
        chosen[fam] = {}
        for l in lines:
            others = [x for x in lines if x != l]
            best = max(have, key=lambda al: statistics.mean(
                pds[(o, f'{fam}_a{al:g}')] - pds[(o, 'refit_alpha0')] for o in others))
            chosen[fam][l] = best
            arm_map[l].append(f'{fam}_a{best:g}')
    # context-specificity controls carry whatever tag the run produced
    for l in lines:
        for t in sorted({k[1] for k in pds if k[0] == l and k[1].startswith('BEST_')}):
            arm_map[l].append(t)
    json.dump(dict(arm_map=arm_map, alpha_choice=chosen,
                   rule='alpha for line L maximises mean (pds - refit) over the other four lines'),
              open(a.out, 'w'), indent=2)
    print(f"{'family':16s} " + ' '.join(f'{l:>8s}' for l in lines))
    for fam, d in chosen.items():
        print(f"{fam:16s} " + ' '.join(f'{d[l]:8g}' for l in lines))
    print('\nper-line arms:')
    for l in lines:
        print(f'  {l}: {len(arm_map[l])} -> {arm_map[l]}')
    print('\ntotal preset runs (x3 seeds):', sum(len(v) for v in arm_map.values()) * 3)


if __name__ == '__main__':
    main()
