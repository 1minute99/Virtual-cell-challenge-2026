# Long AIVC runs

Run from the intended working directory, using the absolute executable for the correct environment:

```bash
scripts/run_detached.sh study-stage-01 /home/wonmin-kim/AIVC/vcc2026_data/reports/MY_STUDY/logs /home/wonmin-kim/miniconda3/envs/ct_transfer/bin/python /absolute/path/to/runner.py --stage STAGE
systemctl --user status aivc-study-stage-01.service
systemctl --user show aivc-study-stage-01.service -p Result -p ExecMainStatus -p MemoryPeak
```

Output goes directly to `LOG_DIRECTORY/NAME.log`. Use `tail -n 50` to inspect it without dumping entire logs into VS Code. Launch success only means the process started; check `Result` and `ExecMainStatus` for completion. Completed units remain loaded, and failed units remain inspectable. Reusing a loaded name is blocked. Inspect saved stage receipts before using a new name or restarting anything. There is no automatic retry.

Jobs share `aivc-jobs.slice`: RAM throttles above 40GiB and is capped at48GiB; swap is capped at2GiB; total CPU is capped at8 cores' equivalent. A job exceeding the hard cap can fail: use its checkpoints to resume after investigating, rather than raising limits blindly. GPU memory is not capped; serialize GPU-heavy jobs. The launcher sets common BLAS/OpenMP thread counts to4. A user service does not inherit an activated shell environment; use the absolute interpreter, or explicitly pass required non-secret environment settings with `/usr/bin/env`.

The slice is installed at `~/.config/systemd/user/aivc-jobs.slice`; its source is `scripts/aivc-jobs.slice`. On another machine install it and run `systemctl --user daemon-reload` before using the launcher. Jobs survive editor/terminal closure, but not reboot and potentially not logout. This is not a scheduler or a checkpoint implementation.

Stop only a job you intend to cancel with `systemctl --user stop UNIT`. After reviewing its outcome, release its name with `systemctl --user stop UNIT` and, for failed jobs, `systemctl --user reset-failed UNIT`.

Integration verification: `python3 scripts/test_run_detached.py` (requires a running user systemd session).
