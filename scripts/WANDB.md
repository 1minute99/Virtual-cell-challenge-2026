# AIVC training logs in Weights & Biases

Default destination is set in `wandb_config.json` (only entity/project; no credentials).
Credentials are already saved in the user's standard W&B login store. To change login,
run `/home/wonmin-kim/miniconda3/envs/ct_transfer/bin/wandb login --relogin --verify`
in your own terminal and enter the key at its hidden prompt. Never put the key in a
command argument, chat, report, or repository file.

For future runners using `context_forcing_model.fit` or `b19_optimization_model.fit`:

```bash
scripts/run_detached.sh experiment-01 /absolute/path/to/new/report/logs \
  /home/wonmin-kim/miniconda3/envs/ct_transfer/bin/python \
  /home/wonmin-kim/AIVC/vcc2026_data/scripts/run_with_wandb.py \
  /absolute/path/to/new_runner.py --your-arguments
```

The wrapper replaces only those two trainers in memory with W&B-enabled versions.
Frozen historical source files and results remain unchanged. Do not rerun completed
studies just to generate charts. Other training frameworks require explicit integration.
New Python runners can instead import `context_forcing_wandb_model` or `b19_wandb_model`
and call `fit(..., tracking={'name': 'descriptive-name', 'group': 'study-name'})`.
Config includes the seed and model/training settings. Use per-fold names/groups in
new runners to distinguish fits. One W&B run is created per actual fit.

Logged metrics:

- `train/batch_loss`: actual penalized minibatch loss, first update and every50 updates.
- `train/batch_mse`: standard trainer's unpenalized minibatch MSE.
- `probe/mse`: fixed training-subset MSE at the trainer's existing measurement points.
- `probe/penalized_loss`: optimization trainer's fixed-probe objective.
- `train/lr`, `train/phase`, and `train/optimizer_step`.
- Correction/context norms; optimization trainer also logs whole-source objectives at checkpoints.

Probe metrics are training diagnostics, not validation. Neither fit function performs
held-out scoring internally, so no validation curve is invented. For a new runner's
separate evaluation stage, use `aivc_wandb.TrainingLog` and log real metrics as
`validation/pds` or another accurately named field. Keep source-inner selection
scores separate from outer evaluation; evaluation logging is not automatically
intercepted by the wrapper.

Scalar JSONL backups and run IDs/URLs live in `vcc2026_data/wandb_runs/ID/`.
No model files, biological matrices, source code or terminal output are uploaded by
this integration. W&B also collects its standard system metadata/metrics.
Initialization failures stop the run before training rather than silently pretending
logging worked. Each fit gets a fresh ID; this integration does not resume model or
optimizer state. Inspect saved scientific checkpoints before restarting a job.

For offline operation, pass `tracking={'mode':'offline'}` to `fit`, or explicitly
launch `/usr/bin/env WANDB_MODE=offline ...` inside the detached job. Shell exports
are not automatically inherited by user-systemd services. Later sync the specific
W&B offline run directory with the training environment's `wandb sync` command.

Verification: `python -m pytest -q vcc2026_data/scripts/test_wandb_training.py`.
Both trainers are checked against the originals for bit-exact final weights.
`wandb_smoke_test.py` is a tiny synthetic online setup test, not a challenge result.
