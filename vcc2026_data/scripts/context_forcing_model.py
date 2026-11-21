"""Deterministic latent context models with optional empirical-baseline input.

The module is deliberately data agnostic.  It consumes fold-local arrays and
never reads response tables, scorer artifacts, or validation data.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import math
import os
from pathlib import Path
import time

import numpy as np
import torch
from torch import nn


@dataclass(frozen=True)
class Config:
    route: str = "residual"
    family: str = "mlp"
    rank: int = 32
    context_dim: int = 16
    hidden: int = 64
    use_baseline: bool = True
    energy: float = 0.0
    weight_decay: float = 0.001
    pair_weight: float = 0.0
    steps: int = 3000
    batch_size: int = 256
    lr: float = 0.003
    mode: str = "correct"
    baseline_dropout: float = 0.0
    coverage_pass: bool = False

    def __post_init__(self):
        if self.route not in ("residual", "direct"):
            raise ValueError(f"Unknown route: {self.route}")
        if self.family not in ("mlp", "bilinear", "interaction"):
            raise ValueError(f"Unknown family: {self.family}")
        if self.mode not in ("correct", "none", "context_shuffle"):
            raise ValueError(f"Unknown context mode: {self.mode}")
        if min(self.rank, self.context_dim, self.hidden, self.steps,
               self.batch_size) <= 0:
            raise ValueError("rank, context_dim, hidden, steps, and batch_size must be positive")
        if self.lr <= 0 or self.weight_decay < 0 or self.energy < 0 or self.pair_weight < 0:
            raise ValueError("optimization settings must be nonnegative and lr must be positive")
        if not 0 <= self.baseline_dropout <= 1:
            raise ValueError("baseline_dropout must be between zero and one")


class Net(nn.Module):
    """A gene path plus a context path in a fold-local latent output space.

    ``mlp`` has a gene-only output plus a strict context interaction.  The
    ``bilinear`` and ``interaction`` families are strict: their complete output
    is exactly zero when centered context is zero.  Direct routes never create
    or inspect an empirical-baseline input path.
    """

    def __init__(self, feature_dim: int, output_rank: int, control_rank: int,
                 config: Config):
        super().__init__()
        if min(feature_dim, output_rank, control_rank) <= 0:
            raise ValueError("Network dimensions must be positive")
        if output_rank != config.rank:
            raise ValueError("output rank does not match Config.rank")
        self.feature_dim = int(feature_dim)
        self.output_rank = int(output_rank)
        self.control_rank = int(control_rank)
        self.config = config
        self.baseline_used = bool(config.route == "residual" and config.use_baseline)
        gene_inputs = self.feature_dim + (self.output_rank if self.baseline_used else 0)
        self.gene = nn.Sequential(nn.Linear(gene_inputs, config.hidden), nn.GELU())
        self.context = nn.Linear(self.control_rank, config.context_dim, bias=False)
        self.context_hidden = nn.Linear(config.context_dim, config.hidden, bias=False)

        if config.family == "mlp":
            self.global_head = nn.Sequential(
                nn.Linear(config.hidden, config.hidden), nn.GELU(),
                nn.Linear(config.hidden, self.output_rank, bias=False))
        else:
            self.global_head = None

        if config.family == "bilinear":
            self.interaction_body = nn.Identity()
            self.interaction_head = nn.Linear(config.hidden, self.output_rank, bias=False)
        else:
            self.interaction_body = nn.Sequential(
                nn.Linear(config.hidden * 2, config.hidden, bias=False), nn.GELU())
            self.interaction_head = nn.Linear(config.hidden, self.output_rank, bias=False)

        # A zero decoder strands both sides of a multiplicative interaction at
        # initialization.  Small nonzero heads keep outputs controlled while
        # allowing gradients through the gene and context encoders immediately.
        nn.init.normal_(self.interaction_head.weight, std=0.01)
        if self.global_head is not None:
            nn.init.normal_(self.global_head[-1].weight, std=0.01)

    def _gene_hidden(self, features, baseline):
        if features.ndim != 2 or features.shape[1] != self.feature_dim:
            raise ValueError("Feature dimension mismatch")
        if self.baseline_used:
            if baseline is None:
                raise ValueError("Residual route with use_baseline=True requires baseline")
            if baseline.ndim != 2 or baseline.shape != (len(features), self.output_rank):
                raise ValueError("Baseline dimension mismatch")
            inputs = torch.cat([features, baseline], dim=-1)
        else:
            # In particular, direct models do not even validate values from an
            # optional baseline tensor; callers may pass one for a shared API.
            inputs = features
        return self.gene(inputs)

    def encode_context(self, controls, zero_context=False):
        if controls.ndim != 2 or controls.shape[1] != self.control_rank:
            raise ValueError("Control dimension mismatch")
        raw = torch.zeros_like(controls) if zero_context else controls
        return torch.tanh(self.context_hidden(torch.tanh(self.context(raw))))

    def forward(self, features, controls, baseline=None, zero_context=False):
        if len(features) != len(controls):
            raise ValueError("Feature and control row counts differ")
        gene = self._gene_hidden(features, baseline)
        context = self.encode_context(
            controls, zero_context=bool(zero_context or self.config.mode == "none"))
        product = gene * context
        if self.config.family == "bilinear":
            interaction = self.interaction_head(product)
        else:
            # Every input to this branch is exactly zero for zero context;
            # bias-free layers and GELU(0)=0 preserve the exact-zero guarantee.
            interaction = self.interaction_head(
                self.interaction_body(torch.cat([product, context], dim=-1)))
        if self.global_head is None:
            return interaction
        return self.global_head(gene) + interaction


class _PairSampler:
    """Uniform compatible-key/group sampling without expanding row pairs.

    Indexing is one pass over rows followed by deterministic sorting of keys
    and group names.  Sampling is reproducible for a fixed NumPy generator;
    each key is equally likely, then an ordered pair of distinct biology groups
    is uniform, then one row is uniform within each chosen group.
    """

    def __init__(self, pair_keys, bios):
        keys = np.asarray(pair_keys).astype(str)
        groups = np.asarray(bios).astype(str)
        if keys.ndim != 1 or groups.ndim != 1 or len(keys) != len(groups):
            raise ValueError("pair_keys and bios must be aligned vectors")
        buckets = {}
        for row, (key, group) in enumerate(zip(keys.tolist(), groups.tolist())):
            by_group = buckets.get(key)
            if by_group is None:
                by_group = {}
                buckets[key] = by_group
            rows = by_group.get(group)
            if rows is None:
                rows = []
                by_group[group] = rows
            rows.append(row)

        eligible_keys = []
        key_group_starts = []
        key_group_counts = []
        group_rows = []
        stored_rows = 0
        for key in sorted(buckets):
            by_group = buckets[key]
            if len(by_group) < 2:
                continue
            row_arrays = tuple(np.asarray(by_group[group], dtype=np.int64)
                               for group in sorted(by_group))
            eligible_keys.append(key)
            key_group_starts.append(len(group_rows))
            key_group_counts.append(len(row_arrays))
            group_rows.extend(row_arrays)
            stored_rows += sum(len(rows) for rows in row_arrays)
        self._eligible_keys = tuple(eligible_keys)
        self._key_group_starts = np.asarray(key_group_starts, dtype=np.int64)
        self._key_group_counts = np.asarray(key_group_counts, dtype=np.int64)
        self._group_row_lengths = np.asarray([len(rows) for rows in group_rows],
                                             dtype=np.int64)
        self._group_row_starts = np.empty(len(group_rows), dtype=np.int64)
        if group_rows:
            self._group_row_starts[0] = 0
            if len(group_rows) > 1:
                np.cumsum(self._group_row_lengths[:-1],
                          out=self._group_row_starts[1:])
            self._pooled_rows = np.concatenate(group_rows)
        else:
            self._pooled_rows = np.empty(0, dtype=np.int64)
        self.stored_rows = int(stored_rows)

    @property
    def eligible_keys(self):
        return self._eligible_keys

    def sample(self, rng: np.random.Generator, size: int):
        if not self._eligible_keys:
            raise ValueError("Pair loss requires compatible cross-biology pair keys")
        if size <= 0:
            raise ValueError("Pair sample size must be positive")
        chosen = rng.integers(0, len(self._eligible_keys), size=size)
        group_counts = self._key_group_counts[chosen]
        first_local = rng.integers(0, group_counts)
        other_local = rng.integers(0, group_counts - 1)
        second_local = other_local + (other_local >= first_local)
        first_group = self._key_group_starts[chosen] + first_local
        second_group = self._key_group_starts[chosen] + second_local
        left_offsets = rng.integers(0, self._group_row_lengths[first_group])
        right_offsets = rng.integers(0, self._group_row_lengths[second_group])
        left = self._pooled_rows[self._group_row_starts[first_group] + left_offsets]
        right = self._pooled_rows[self._group_row_starts[second_group] + right_offsets]
        return left, right


def _pair_contrast_mse(prediction, target, left, right):
    pred_difference = prediction[left] - prediction[right]
    target_difference = target[left] - target[right]
    return torch.mean((pred_difference - target_difference) ** 2)


def _array_sha(value):
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256(str((array.shape, array.dtype.str)).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def _finite_array(name, value, ndim=2):
    array = np.asarray(value)
    if array.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}-dimensional array")
    if not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite numeric values")
    return array


def _validate_prepared(prepared, config):
    required = ("features", "controls", "baseline", "y", "weights", "bios",
                "pair_keys", "ids")
    missing = [name for name in required if not hasattr(prepared, name)]
    if missing:
        raise ValueError(f"Prepared data missing fields: {missing}")
    features = _finite_array("features", prepared.features)
    controls = _finite_array("controls", prepared.controls)
    baseline = _finite_array("baseline", prepared.baseline)
    target = _finite_array("y", prepared.y)
    n = len(features)
    if n == 0 or any(len(value) != n for value in (controls, baseline, target)):
        raise ValueError("Prepared row counts differ or are empty")
    if target.shape[1] != config.rank or baseline.shape[1] != config.rank:
        raise ValueError("Prepared target/baseline rank does not match Config.rank")
    weights = _finite_array("weights", prepared.weights, ndim=1).astype(np.float64)
    if len(weights) != n or np.any(weights <= 0):
        raise ValueError("weights must be positive and aligned with rows")
    bios = np.asarray(prepared.bios).astype(str)
    pair_keys = np.asarray(prepared.pair_keys).astype(str)
    ids = np.asarray(prepared.ids)
    if bios.ndim != 1 or pair_keys.ndim != 1 or ids.ndim != 1 or any(
            len(value) != n for value in (bios, pair_keys, ids)):
        raise ValueError("bios, pair_keys, and ids must be aligned vectors")
    if not np.issubdtype(ids.dtype, np.integer):
        raise ValueError("ids must be integers")
    if len(np.unique(ids)) != n:
        raise ValueError("ids must be unique")
    return features, controls, baseline, target, weights, bios, pair_keys, ids


def _derange_groups(controls, bios, rng):
    groups = sorted(set(bios.tolist()))
    if len(groups) < 2:
        raise ValueError("context_shuffle requires at least two biological groups")
    order = rng.permutation(len(groups))
    while np.any(order == np.arange(len(groups))):
        order = rng.permutation(len(groups))
    means = {group: controls[bios == group].mean(axis=0) for group in groups}
    shuffled = np.array(controls, copy=True)
    mapping = {}
    for destination_index, destination in enumerate(groups):
        source = groups[int(order[destination_index])]
        shuffled[bios == destination] = means[source]
        mapping[destination] = source
    return shuffled, mapping


def _to_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise FloatingPointError("Nonfinite training diagnostic")
    return result


@dataclass
class Fitted:
    module: Net
    config: Config
    seed: int
    receipt: dict

    def predict(self, features, controls, baseline=None, zero_context=False):
        device = next(self.module.parameters()).device
        f = _finite_array("features", features)
        c = _finite_array("controls", controls)
        if len(f) != len(c):
            raise ValueError("Prediction row counts differ")
        # Direct and baseline-free residual routes intentionally ignore the
        # supplied baseline, including arbitrary values from shared callers.
        if self.module.baseline_used:
            if baseline is None:
                raise ValueError("This fitted residual model requires baseline")
            a = _finite_array("baseline", baseline)
            if a.shape != (len(f), self.config.rank):
                raise ValueError("Prediction baseline dimension mismatch")
            base_tensor = torch.as_tensor(np.array(a, copy=True), dtype=torch.float32,
                                          device=device)
        else:
            base_tensor = None
        self.module.eval()
        with torch.no_grad():
            prediction = self.module(
                torch.as_tensor(np.array(f, copy=True), dtype=torch.float32, device=device),
                torch.as_tensor(np.array(c, copy=True), dtype=torch.float32, device=device),
                base_tensor,
                zero_context=bool(zero_context or self.config.mode == "none"))
        result = prediction.detach().cpu().numpy()
        if result.shape != (len(f), self.config.rank) or not np.isfinite(result).all():
            raise FloatingPointError("Nonfinite or malformed model prediction")
        return result

    def save(self, path):
        state = {key:value.detach().cpu() for key, value in self.module.state_dict().items()}
        if not all(torch.isfinite(value).all() for value in state.values()):
            raise FloatingPointError("Refusing to save nonfinite model tensors")
        _validate_receipt(self.receipt)
        if int(self.receipt.get("seed", self.seed)) != int(self.seed):
            raise ValueError("Fitted seed does not match receipt seed")
        payload = dict(state=state, config=asdict(self.config), seed=int(self.seed),
                       receipt=self.receipt, feature_dim=self.module.feature_dim,
                       output_rank=self.module.output_rank,
                       control_rank=self.module.control_rank)
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        torch.save(payload, destination)

    @classmethod
    def load(cls, path, device=None):
        saved = torch.load(path, map_location="cpu", weights_only=False)
        config = Config(**saved["config"])
        if int(saved["output_rank"]) != config.rank:
            raise ValueError("Saved output rank does not match saved config")
        module = Net(int(saved["feature_dim"]), int(saved["output_rank"]),
                     int(saved["control_rank"]), config)
        module.load_state_dict(saved["state"], strict=True)
        if not all(torch.isfinite(value).all() for value in module.state_dict().values()):
            raise FloatingPointError("Loaded model contains nonfinite tensors")
        receipt = saved["receipt"]
        _validate_receipt(receipt)
        if int(receipt.get("seed", saved["seed"])) != int(saved["seed"]):
            raise ValueError("Saved seed does not match receipt seed")
        target = device or "cpu"
        module.to(target)
        return cls(module=module, config=config, seed=int(saved["seed"]), receipt=receipt)


def _validate_receipt(value):
    if isinstance(value, dict):
        for item in value.values():
            _validate_receipt(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _validate_receipt(item)
    elif isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        raise FloatingPointError("Receipt contains a nonfinite number")


def fit(prepared, config: Config, seed=0, device=None):
    """Fit a deterministic latent predictor on permitted genetic rows only."""
    start = time.monotonic()
    arrays = _validate_prepared(prepared, config)
    features, controls, baseline, target, weights, bios, pair_keys, ids = arrays
    pair_sampler = _PairSampler(pair_keys, bios)
    if config.pair_weight > 0 and not pair_sampler.eligible_keys:
        raise ValueError("Pair loss requires compatible cross-biology pair keys")

    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))
    torch.set_num_threads(2)
    target_device = device or ("cuda" if torch.cuda.is_available() else "cpu")

    shuffle_mapping = {}
    training_controls = np.asarray(controls, dtype=np.float32)
    if config.mode == "context_shuffle":
        training_controls, shuffle_mapping = _derange_groups(
            training_controls, bios, np.random.default_rng(38111 + int(seed)))
    elif config.mode == "none":
        training_controls = np.zeros_like(training_controls)

    module = Net(features.shape[1], target.shape[1], controls.shape[1], config).to(target_device)
    f = torch.as_tensor(np.array(features, copy=True), dtype=torch.float32, device=target_device)
    c = torch.as_tensor(np.array(training_controls, copy=True), dtype=torch.float32,
                        device=target_device)
    a = torch.as_tensor(np.array(baseline, copy=True), dtype=torch.float32, device=target_device)
    y = torch.as_tensor(np.array(target, copy=True), dtype=torch.float32, device=target_device)
    optimizer = torch.optim.AdamW(module.parameters(), lr=config.lr,
                                  weight_decay=config.weight_decay)

    probabilities = weights / weights.sum()
    sample_rng = np.random.default_rng(19317 + int(seed))
    coverage_rng = np.random.default_rng(829317 + int(seed))
    dropout_rng = np.random.default_rng(319317 + int(seed))
    pair_rng = np.random.default_rng(519317 + int(seed))
    probe_rng = np.random.default_rng(619317 + int(seed))
    n = len(f)
    probe = np.sort(probe_rng.choice(n, size=min(n, 512), replace=False))
    if pair_sampler.eligible_keys:
        probe_left, probe_right = pair_sampler.sample(probe_rng, min(512, max(32, n)))
    else:
        probe_left = probe_right = np.empty(0, dtype=np.int64)

    covered = np.zeros(n, dtype=bool)
    counts = dict(coverage_batches=0, coverage_examples=0, sampled_batches=0,
                  sampled_examples=0, pair_batches=0, pair_examples=0,
                  optimizer_steps=0)
    coverage_sizes = []
    sampled_sizes = []
    trace = []

    def baseline_for(row_ids):
        chosen = a[row_ids]
        if not module.baseline_used or config.baseline_dropout == 0:
            return chosen if module.baseline_used else None
        if config.baseline_dropout == 1:
            return torch.zeros_like(chosen)
        keep = dropout_rng.random(len(row_ids)) >= config.baseline_dropout
        mask = torch.as_tensor(keep[:, None], dtype=torch.float32, device=target_device)
        return chosen * mask

    def pair_value(size):
        left, right = pair_sampler.sample(pair_rng, size)
        li = torch.as_tensor(left, dtype=torch.long, device=target_device)
        ri = torch.as_tensor(right, dtype=torch.long, device=target_device)
        pred_left = module(f[li], c[li], baseline_for(left),
                           zero_context=config.mode == "none")
        pred_right = module(f[ri], c[ri], baseline_for(right),
                            zero_context=config.mode == "none")
        return torch.mean(((pred_left - pred_right) - (y[li] - y[ri])) ** 2)

    def update(row_ids, phase):
        row_ids = np.asarray(row_ids, dtype=np.int64)
        ix = torch.as_tensor(row_ids, dtype=torch.long, device=target_device)
        optimizer.zero_grad(set_to_none=True)
        prediction = module(f[ix], c[ix], baseline_for(row_ids),
                            zero_context=config.mode == "none")
        mse = torch.mean((prediction - y[ix]) ** 2)
        total = mse + config.energy * torch.mean(prediction ** 2)
        if config.pair_weight > 0:
            contrast = pair_value(len(row_ids))
            total = total + config.pair_weight * contrast
            counts["pair_batches"] += 1
            counts["pair_examples"] += len(row_ids)
        if not torch.isfinite(total):
            raise FloatingPointError("Nonfinite context-forcing loss")
        total.backward()
        nn.utils.clip_grad_norm_(module.parameters(), 5.0)
        optimizer.step()
        counts["optimizer_steps"] += 1
        covered[row_ids] = True
        if phase == "coverage":
            counts["coverage_batches"] += 1
            counts["coverage_examples"] += len(row_ids)
            coverage_sizes.append(int(len(row_ids)))
        else:
            counts["sampled_batches"] += 1
            counts["sampled_examples"] += len(row_ids)
            sampled_sizes.append(int(len(row_ids)))

    def measure(phase, step):
        module.eval()
        with torch.no_grad():
            pi = torch.as_tensor(probe, dtype=torch.long, device=target_device)
            prediction = module(f[pi], c[pi], a[pi] if module.baseline_used else None,
                                zero_context=config.mode == "none")
            train_mse = torch.mean((prediction - y[pi]) ** 2)
            output_rms = torch.sqrt(torch.mean(prediction ** 2))
            context = module.encode_context(c[pi], zero_context=config.mode == "none")
            context_rms = torch.sqrt(torch.mean(context ** 2))
            if len(probe_left):
                li = torch.as_tensor(probe_left, dtype=torch.long, device=target_device)
                ri = torch.as_tensor(probe_right, dtype=torch.long, device=target_device)
                pl = module(f[li], c[li], a[li] if module.baseline_used else None,
                            zero_context=config.mode == "none")
                pr = module(f[ri], c[ri], a[ri] if module.baseline_used else None,
                            zero_context=config.mode == "none")
                contrast = torch.mean(((pl - pr) - (y[li] - y[ri])) ** 2)
                contrast_mse = _to_float(contrast)
            else:
                contrast_mse = None
        trace.append(dict(phase=phase, step=int(step),
                          optimizer_step=int(counts["optimizer_steps"]),
                          train_mse=_to_float(train_mse),
                          contrast_mse=contrast_mse,
                          output_rms=_to_float(output_rms),
                          context_rms=_to_float(context_rms)))
        module.train()

    measure("initial", 0)
    if config.coverage_pass:
        order = coverage_rng.permutation(n)
        for lo in range(0, n, config.batch_size):
            update(order[lo:lo + config.batch_size], "coverage")
        measure("coverage", counts["optimizer_steps"])

    log_every = max(1, config.steps // 8)
    for step in range(config.steps):
        batch = sample_rng.choice(n, size=config.batch_size, replace=True, p=probabilities)
        update(batch, "sampled")
        if (step + 1) % log_every == 0 or step + 1 == config.steps:
            measure("sampled", step + 1)

    module.eval()
    receipt = dict(
        **counts,
        examples_seen=int(counts["coverage_examples"] + counts["sampled_examples"]),
        unique_rows_seen=int(covered.sum()),
        coverage_batch_sizes=coverage_sizes,
        sampled_batch_sizes=sampled_sizes,
        training_rows=n,
        training_ids_sha256=_array_sha(ids.astype(np.int64)),
        feature_dim=int(features.shape[1]), control_rank=int(controls.shape[1]),
        output_rank=int(target.shape[1]), parameters=sum(p.numel() for p in module.parameters()),
        route=config.route, family=config.family, empirical_baseline_used=module.baseline_used,
        direct_empirical_baseline_used=False,
        mode=config.mode, context_shuffle=shuffle_mapping,
        pair_eligible_keys=len(pair_sampler.eligible_keys),
        pair_rows_stored=pair_sampler.stored_rows,
        pair_sampling="uniform compatible key, then uniform distinct biology groups; no expanded pairs",
        fixed_probe_rows=int(len(probe)), fixed_probe_pairs=int(len(probe_left)),
        trace=trace, seed=int(seed), device=str(target_device),
        coverage_pass=bool(config.coverage_pass), baseline_dropout=float(config.baseline_dropout),
        deterministic_algorithms=torch.are_deterministic_algorithms_enabled(),
        allow_tf32=False, cublas_workspace=os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
        torch_version=torch.__version__, cuda_version=torch.version.cuda,
        seconds=float(time.monotonic() - start),
    )
    _validate_receipt(receipt)
    return Fitted(module=module, config=config, seed=int(seed), receipt=receipt)
