"""Task 10 (plan zero-shot-limits-context-complexity-plan): seed-s rerun of the 2026-09-22 architecture screen.

Runs vcc2026_data/scripts/arch_screen.py UNCHANGED except for its module-level SEED (seed 0 in the original
screen) and its output folder. SEED controls (i) torch.manual_seed before model construction, (ii) Z.train's
torch/numpy minibatch RNG and (iii) the harness-PDS reconstruction seed stored in scores.json (which the official
preset step later re-uses as its reconstruction seed and PDS-parity key). Data preparation, bases, context space,
early-stopping rule and evaluation are fixed seeds (20260922) exactly as in seed 0.

Assumptions (stated, see knowledge/synthesis/rq2_arch_seeds-assumptions.md):
- "same protocol as seed 0" = same driver code, same data, same alpha-1 correct-context official preset; only the
  seed index changes, and it changes training AND reconstruction together (as in seed 0, where both were 0).
- The transfer baseline (refit_v2) is deterministic in training; its official composite is re-scored at the same
  reconstruction seed s so every architecture-minus-baseline difference is paired at seed s.
Usage: python rq2_arch_seed_screen.py --seed S --line K562|HepG2 --out DIR [--archs a,b,...] [--device cuda]
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, '/home/wonmin-kim/AIVC/vcc2026_data/scripts')
import arch_screen as A

SIX = ['additive', 'moe', 'bilinear', 'full_hyper', 'token_transformer', 'cross_attn']

ap = argparse.ArgumentParser()
ap.add_argument('--seed', type=int, required=True)
ap.add_argument('--line', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--archs', default=','.join(SIX))
ap.add_argument('--device', default='cuda')
a = ap.parse_args()
Path(a.out).mkdir(parents=True, exist_ok=True)
A.SEED = a.seed
print(f'{{"event": "seed_wrapper", "seed": {A.SEED}, "line": "{a.line}", "out": "{a.out}", "archs": "{a.archs}"}}', flush=True)
sys.argv = ['arch_screen.py', '--archs', a.archs, '--device', a.device, '--line', a.line, '--out', a.out]
A.main()
