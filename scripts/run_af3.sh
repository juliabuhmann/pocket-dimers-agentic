#!/usr/bin/env bash
# usage: run_af3.sh INPUT.json OUTDIR GPU_ID [extra run_alphafold.py flags]
# Turing (cc 7.5) settings: XLA attention, custom-kernel-fusion-rewriter disabled (AF3 docs, known_issues.md).
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
source "$WS/scripts/env_threads.sh"
source "$WS/conda_envs/af3/bin/activate"
IN=$1; OUT=$2; GPU=$3; shift 3
export CUDA_VISIBLE_DEVICES=$GPU
export XLA_FLAGS=${XLA_FLAGS_OVERRIDE:-"--xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_autotune_level=0"}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export AF3_SRC="$WS/tools/alphafold3"
python "$WS/scripts/run_af3_wrapper.py" --json_path="$IN" --output_dir="$OUT" \
  --model_dir="$WS/data/af3/models" --norun_data_pipeline \
  --flash_attention_implementation=xla --jax_compilation_cache_dir="$WS/data/af3/jax_cache" "$@"
