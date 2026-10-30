#!/usr/bin/env python3
"""Thin wrapper around AF3 run_alphafold.py that lets us override global_config.bfloat16
(env AF3_BFLOAT16 = all|intermediate|none) for non-bf16 GPUs. All CLI flags pass through."""
import os, sys
sys.path.insert(0, os.environ["AF3_SRC"])
import run_alphafold as ra
from absl import app

_orig = ra.make_model_config
def patched(**kw):
    cfg = _orig(**kw)
    mode = os.environ.get("AF3_BFLOAT16")
    if mode:
        cfg.global_config.bfloat16 = mode
        print(f"[wrapper] global_config.bfloat16 = {mode}", flush=True)
    return cfg
ra.make_model_config = patched

if __name__ == "__main__":
    ra.flags.mark_flags_as_required(['output_dir'])
    app.run(ra.main)
