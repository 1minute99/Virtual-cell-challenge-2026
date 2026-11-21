#!/usr/bin/env python3
"""Ingest the Feng TargetedScreen arm (19 iPSC lines x 444 targets) with the frozen
`ingest_feng.py` machinery and two arm-specific fixes, applied by monkeypatching rather than
by editing the historical script:

  * barcodes: the RNA header spells a cell `P4_I73_<BC>-1` while the metadata says
    `PC-P4-D3_I73_<BC>-1` (pool, day, inlet). Inlet numbers are unique within the arm, so both
    are normalised to `P{pool}_I{inlet}_<BC>-1`; the run still aborts below 95 % alignment.
  * controls: this arm's non-targeting guides are labelled `NonTarget_<protospacer>`, which
    the frozen NTC regex (written for `non-targeting`/`NTC`) does not match. They are
    confirmed non-targeting controls (8,241 cells, ~490 per line) and are labelled as such.
Outputs: pseudobulk stores `feng_targetedscreen_target` (pooled) and
`feng_targetedscreen_target_line` (target x line) plus their effects (matched control = the
same line's NonTarget cells; min 10 perturbed / 50 control cells).
"""
import re, sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import ingest_feng as IF

_META = re.compile(r'^(?:MP-)?PC-P(\d+)-D\d+_I(\d+)_(.+)$')
_HEAD = re.compile(r'^P(\d+)_I(\d+)_(.+)$')


def norm_bc(s):
    s = str(s).replace('.', '-')
    m = _META.match(s)
    if m:
        return f'P{m.group(1)}_I{m.group(2)}_{m.group(3)}'
    m = _HEAD.match(s)
    if m:
        return f'P{m.group(1)}_I{m.group(2)}_{m.group(3)}'
    return IF.norm_bc(s)


IF.norm_bc = norm_bc
IF.NTC_EXACT = re.compile(
    r"^(non[-_]?targeting|non[-_]?target|ntc|scrambled?|safe[-_]?harbou?r|no[-_]?target)"
    r"([-_]?\d+)?$", re.I)

if __name__ == '__main__':
    sys.argv = [sys.argv[0], '--experiments', 'TargetedScreen'] + sys.argv[1:]
    IF.main()
