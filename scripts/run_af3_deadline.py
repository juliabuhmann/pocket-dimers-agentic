#!/usr/bin/env python3
"""Wrapper around scripts/run_af3_wrapper.py (imported, unmodified) that stops LAUNCHING new fold jobs after a
wall-clock deadline. A job already running is finished. Deadline: env AF3_LAUNCH_DEADLINE (ISO 8601 with offset,
e.g. 2026-10-30T04:30:00+01:00). Used by scripts/cofold_run_hetero_ext.sh. All CLI flags pass through to AF3."""
import datetime as dt
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_af3_wrapper as w  # noqa: E402  (patches bf16 config; exposes run_alphafold as w.ra)

DEADLINE = dt.datetime.fromisoformat(os.environ["AF3_LAUNCH_DEADLINE"])
_orig = w.ra.process_fold_input


def gated(**kw):
    now = dt.datetime.now(dt.timezone.utc)
    name = kw["fold_input"].name
    if now >= DEADLINE:
        print(f"[deadline] {now.isoformat()} >= {DEADLINE.isoformat()}: NOT launching {name}", flush=True)
        return None
    print(f"[deadline] launching {name} at {now.isoformat()}", flush=True)
    return _orig(**kw)


w.ra.process_fold_input = gated

if __name__ == "__main__":
    w.ra.flags.mark_flags_as_required(["output_dir"])
    w.app.run(w.ra.main)
