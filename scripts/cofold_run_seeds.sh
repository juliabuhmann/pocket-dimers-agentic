#!/usr/bin/env bash
# Run the replicate-seed AF3 jobs (scripts/cofold_seeds_jobs.py) on two GPUs, one process per GPU.
# Same AF3 settings as scripts/cofold_run_onaxis.sh / cofold_run_dimeronly.sh; only the input dir / output dir differ.
# Restart loop: a long-lived AF3 process hit RESOURCE_EXHAUSTED (OOM) on the ~1100-token dimers after many different
# buckets had been compiled in the same process (both GPUs, seed-1 SGPL1/PCCB). On exit the job dirs are rebuilt
# (finished jobs skipped) and fresh processes started, up to MAX_BATCHES, until nothing is left or DEADLINE_EPOCH passes.
# A job that is first in a fresh process and still unfinished afterwards is written to data/cofold/seeds/skip.txt.
# usage: cofold_run_seeds.sh GPU [GPU ...] [--seeds 1 2]   (GPU 0 is faulty: never use it)
set -uo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
GPUS=(); while [[ $# -gt 0 && "$1" != --* ]]; do GPUS+=("$1"); shift; done   # any number of GPU ids, then optional flags
[[ ${#GPUS[@]} -ge 1 ]] || { echo "need GPU ids"; exit 1; }
for g in "${GPUS[@]}"; do [[ "$g" == 0 ]] && { echo "GPU 0 is faulty, refusing"; exit 1; }; done
DEADLINE_EPOCH=${DEADLINE_EPOCH:-$(date -d "2026-10-30 03:13" +%s)}
MAX_BATCHES=${MAX_BATCHES:-8}
source "$WS/scripts/env_threads.sh"
source "$WS/scripts/af3_local_env.sh"   # node-local AF3 env/weights (stall mitigation)
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_autotune_level=0"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
for b in $(seq 1 $MAX_BATCHES); do
  [[ $(date +%s) -ge $DEADLINE_EPOCH ]] && { echo "deadline reached before batch $b"; break; }
  python "$WS/scripts/cofold_seeds_jobs.py" "${GPUS[@]}" "$@" | tee "$WS/logs/cofold_seeds/batch$b.jobs"
  grep -q "^jobs 0 " "$WS/logs/cofold_seeds/batch$b.jobs" && break
  for G in "${GPUS[@]}"; do
    ( CUDA_VISIBLE_DEVICES=$G python "$WS/scripts/run_af3_wrapper.py" --input_dir="$WS/data/cofold/seeds/jobs_gpu$G" \
        --output_dir="$WS/data/cofold/seeds/af3" --model_dir="$AF3_MODEL_DIR" --norun_data_pipeline \
        --flash_attention_implementation=xla --jax_compilation_cache_dir="$AF3_JAX_CACHE" --force_output_dir \
        >> "$WS/logs/cofold_seeds/af3_gpu$G.log" 2>&1; echo "EXIT $? batch $b" >> "$WS/logs/cofold_seeds/af3_gpu$G.log" ) &
  done
  wait
  for G in "${GPUS[@]}"; do   # first job of this batch's dir, if unfinished -> skip list
    f=$(ls "$WS/data/cofold/seeds/jobs_gpu$G" | head -1); [[ -z "$f" ]] && continue
    n=${f#*_}; n=${n%.json}
    [[ -f "$WS/data/cofold/seeds/af3/$n/${n}_ranking_scores.csv" ]] || echo "$n" >> "$WS/data/cofold/seeds/skip.txt"
  done
done
