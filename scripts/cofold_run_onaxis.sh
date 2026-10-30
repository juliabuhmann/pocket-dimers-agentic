#!/usr/bin/env bash
# Run AF3 for the deep-dive #1 cases (results/cofold/onaxis/cases.tsv) on two GPUs, one process per GPU,
# each process looping over its own input directory (model loaded and XLA buckets compiled once per process).
# usage: cofold_run_onaxis.sh GPU_A GPU_B [--no-split]   (GPU 0 is faulty: never use it; see docs/COFOLDING.md)
#   --no-split: run whatever is already in the job dirs (e.g. after scripts/cofold_seed1_onaxis.py)
# Job split: scripts/cofold_split_jobs.py (case priority order, greedy balance on estimated runtime).
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
GA=${1:?gpu}; GB=${2:?gpu}
[[ "$GA" == 0 || "$GB" == 0 ]] && { echo "GPU 0 is faulty, refusing"; exit 1; }
source "$WS/scripts/env_threads.sh"
[[ "${3:-}" == "--no-split" ]] || "$WS/.venvs/part1/bin/python" "$WS/scripts/cofold_split_jobs.py" "$GA" "$GB"
source "$WS/conda_envs/af3/bin/activate"
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_autotune_level=0"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export AF3_SRC="$WS/tools/alphafold3"
for G in "$GA" "$GB"; do
  ( CUDA_VISIBLE_DEVICES=$G python "$WS/scripts/run_af3_wrapper.py" --input_dir="$WS/data/cofold/onaxis/jobs_gpu$G" \
      --output_dir="$WS/data/cofold/onaxis/af3" --model_dir="$WS/data/af3/models" --norun_data_pipeline \
      --flash_attention_implementation=xla --jax_compilation_cache_dir="$WS/data/af3/jax_cache" --force_output_dir \
      >> "$WS/logs/cofold_onaxis/af3_gpu$G.log" 2>&1; echo "EXIT $?" >> "$WS/logs/cofold_onaxis/af3_gpu$G.log" ) &
done
wait
