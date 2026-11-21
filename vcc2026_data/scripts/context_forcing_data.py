"""Genetic-only, fold-local coordinates for context forcing and direct models."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import simple_residual_data as D

OUT = D.ROOT / 'reports/context_forcing_2026-09-16'
TABLE = Path(os.environ['TV2_BASE_TABLE']) if os.environ.get('TV2_BASE_TABLE') else D.OUT / 'data/table'   # 2026-10-27: y/controls deleted 10-11; S4 rebuild reads only features/genes (stub dir)
LINES = ('K562', 'RPE1', 'HepG2', 'Jurkat', 'KOLF')


def pair_keys(targets, sources):
    keys = []
    for target, source in zip(targets, sources):
        source = str(source)
        if source.startswith('jiang/'):
            family = '/'.join(source.split('/')[:2])
        elif source.startswith('rn/'):
            family = 'rn'
        elif source.startswith('xl/orion:'):
            family = 'orion'
        else:
            family = source  # No uncontrolled cross-study or stimulus contrasts.
        keys.append(family + '|' + str(target))
    return np.asarray(keys)


def permitted_ids(table, excluded=(), targets=None):
    ids = D.scope(table, excluded_bios=excluded, modality=0)
    if targets is not None:
        ids = ids[np.isin(table.targets[ids], list(targets))]
    return ids


def pilot_targets(table):
    import vcclib as V
    old = json.loads((D.OUT/'mock/selection.json').read_text())
    selected = {D.canonical_target(t) for t in old['targets'] + V.official_targets()}
    extras = sorted(set(table.targets[table.modality == 0]) - selected)
    selected.update(np.random.default_rng(20260917).choice(extras, min(768, len(extras)), replace=False))
    return sorted(selected)


def _center_labels(values, targets, bios, source_keys, compatible_keys):
    """Use one target panel for every study/replicate in a paired family."""
    centered = np.array(values, copy=True)
    targets = np.asarray(targets).astype(str)
    families = np.array([str(key).rsplit('|', 1)[0] for key in compatible_keys])
    support = {}
    for family in np.unique(families):
        family_rows = np.flatnonzero(families == family)
        groups = {key: family_rows[source_keys[family_rows] == key]
                  for key in np.unique(source_keys[family_rows])}
        biological_count = len(set(bios[family_rows]))
        target_sets = {key: set(targets[rows]) for key, rows in groups.items()}
        common = sorted(set.intersection(*target_sets.values())) if biological_count >= 2 else None
        if common is not None and not common:
            raise ValueError(f'Empty shared target panel for cross-biology family {family}')
        support[str(family)] = dict(
            mode='shared_equal_target_means' if common is not None else 'own_group_mean',
            biological_group_count=biological_count, source_group_count=len(groups),
            shared_target_count=len(common) if common is not None else None,
            shared_targets_sha256=hashlib.sha256('\n'.join(common).encode()).hexdigest()
                if common is not None else None,
            group_target_counts={str(key): len(panel) for key, panel in target_sets.items()})
        for rows in groups.values():
            if common is None:
                center = values[rows].mean(0)
            else:
                # Average duplicate rows within target, then targets equally.
                unique, inverse = np.unique(targets[rows], return_inverse=True)
                sums = np.zeros((len(unique), values.shape[1]), dtype=np.float64)
                np.add.at(sums, inverse, values[rows])
                means = sums / np.bincount(inverse)[:, None]
                positions = {target: i for i, target in enumerate(unique)}
                center = means[[positions[target] for target in common]].mean(0)
            centered[rows] -= center.astype(centered.dtype, copy=False)
    return centered, support


@dataclass
class Prepared:
    features: np.ndarray
    controls: np.ndarray
    baseline: np.ndarray
    y: np.ndarray
    weights: np.ndarray
    bios: np.ndarray
    pair_keys: np.ndarray
    ids: np.ndarray
    output_basis: np.ndarray
    scale: float
    control_mean: np.ndarray
    control_basis: np.ndarray
    control_scale: np.ndarray
    metadata: dict

    def query(self, controls, baseline=None):
        control = np.asarray(controls)
        h = np.clip(((control-self.control_mean) @ self.control_basis)/self.control_scale, -5, 5)
        if self.metadata['route'] == 'direct':
            a = np.zeros((len(control), self.output_basis.shape[1]), dtype=np.float32)
        else:
            if baseline is None:
                raise ValueError('Residual query requires its permitted-source anchor')
            a = np.asarray(baseline) @ self.output_basis / self.scale
        return h.astype(np.float32), a.astype(np.float32)

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.with_suffix('.tmp').open('wb') as f:
            np.savez(f, **{k:v for k,v in vars(self).items() if k != 'metadata'},
                     metadata=np.array(json.dumps(self.metadata)))
        path.with_suffix('.tmp').replace(path)

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as z:
            fields = {k:z[k] for k in z.files if k != 'metadata'}
            fields['metadata'] = json.loads(str(z['metadata']))
        fields['scale'] = float(fields['scale'])
        return cls(**fields)


def prepare(table, ids, route='residual', rank=32, anchors=None):
    ids = np.asarray(ids, dtype=np.int64)
    if route not in ('residual', 'direct'):
        raise ValueError('Unknown route')
    if len(ids) < 2 or len(set(ids)) != len(ids) or np.any(table.modality[ids] != 0):
        raise ValueError('Preparation requires unique permitted genetic rows')
    weights = D.balanced_weights(table, ids)
    gy = np.array(table.y[ids], copy=True)
    if not np.isfinite(gy).all():
        raise ValueError('Nonfinite response in permitted source rows')
    if route == 'direct':
        if anchors is not None:
            raise ValueError('Direct preparation forbids an empirical anchor')
        ga = np.zeros_like(gy)
    else:
        ga = D.training_anchors(table, ids, rank=32) if anchors is None else np.asarray(anchors)
        if ga.shape != gy.shape or not np.isfinite(ga).all():
            raise ValueError('Anchor shape or values invalid')
    response = gy-ga
    design = response-np.sum(response*weights[:,None], axis=0)
    design *= np.sqrt(weights*len(weights))[:,None]
    basis = D.basis(design, rank)
    projected = response @ basis
    scale = max(float(np.sqrt(np.sum(weights[:,None]*projected**2)/basis.shape[1])), 1e-4)
    bios = np.array([D.canonical_bio(b) for b in table.bios[ids]])
    keys = np.array([b+'|'+str(s) for b,s in zip(bios, table.sources[ids])])
    _, first = np.unique(keys, return_index=True)
    untreated = np.asarray(table.controls[ids[first]], dtype=np.float32)
    control_mean = untreated.mean(0)
    control_basis = D.basis(untreated-control_mean, min(32, max(len(untreated)-1, 1)))
    cp = (untreated-control_mean) @ control_basis
    control_scale = np.maximum(np.sqrt(np.mean(cp**2, axis=0)), 1e-3)
    controls = np.clip(((table.controls[ids]-control_mean) @ control_basis)/control_scale, -5, 5)
    compatible = pair_keys(table.targets[ids], table.sources[ids])
    # Both main and pair objectives use the same source-only centering. Within
    # a compatible cross-biology family, every study/replicate uses the same
    # shared target panel, preventing unequal assay panels from shifting pairs.
    y, centering_support = _center_labels(projected/scale, table.targets[ids],
                                          bios, keys, compatible)
    metadata = dict(route=route, rank=int(basis.shape[1]), requested_rank=int(rank),
        rows=len(ids), targets=len(set(table.targets[ids])), bios=sorted(set(bios)),
        anchor='none' if route == 'direct' else 'leave-destination-biology-out source-refitted rank32',
        basis='source genetic effect' if route == 'direct' else 'source genetic residual',
        control_profiles=len(untreated), control_rank=control_basis.shape[1],
        label_centering='shared equal-target panel within each compatible cross-biology family; otherwise own group mean; no response mean added at prediction',
        label_centering_support=centering_support,
        native_effect_energy=float(np.sum(weights[:,None]*(gy @ basis)**2)/max(np.sum(weights[:,None]*gy**2),1e-20)),
        objective_energy=float(np.sum(weights[:,None]*projected**2)/max(np.sum(weights[:,None]*response**2),1e-20)))
    return Prepared(np.array(table.features[ids], dtype=np.float32), controls.astype(np.float32),
        (ga @ basis / scale).astype(np.float32), y.astype(np.float32), weights,
        bios, compatible, ids.copy(), basis,
        scale, control_mean, control_basis, control_scale, metadata)


def full_delta(prepared, latent, columns, n_features=18533, empirical=None, alpha=1.):
    latent = np.asarray(latent, dtype=np.float64)
    latent = latent-latent.mean(0, keepdims=True)
    effect = (latent*prepared.scale) @ prepared.output_basis.astype(np.float64).T
    columns = np.asarray(columns, dtype=int)
    if prepared.metadata['route'] == 'direct':
        if empirical is not None:
            raise ValueError('Direct prediction forbids empirical mean transfer')
        delta = np.zeros((len(latent), n_features), dtype=np.float32)
    else:
        if empirical is None or empirical.shape != (len(latent), n_features):
            raise ValueError('Residual prediction requires full matching baseline')
        delta = np.asarray(empirical, dtype=np.float32).copy()
    delta[:, columns] += (float(alpha)*effect).astype(np.float32)
    if not np.isfinite(delta).all():
        raise FloatingPointError('Nonfinite prediction')
    return delta
