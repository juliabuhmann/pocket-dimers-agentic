#!/usr/bin/env bash
# Run AF3 for the EXPLORATORY heterodimer MBS extension (results/cofold/hetero_ext/cases.tsv; inputs written by
# scripts/cofold_select_hetero_ext.py) in ONE process on ONE GPU, file-name (= priority) order, finished jobs skipped.
# Same AF3 settings as scripts/cofold_run_h4.sh; adds a launch deadline (scripts/run_af3_deadline.py).
# usage: cofold_run_hetero_ext.sh GPU DEADLINE_ISO     e.g. 1 2026-10-30T04:30:00+01:00   (GPU 0 is faulty)
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
G=${1:?gpu}; export AF3_LAUNCH_DEADLINE=${2:?deadline}
[[ "$G" == 0 ]] && { echo "GPU 0 is faulty, refusing"; exit 1; }
source "$WS/scripts/env_threads.sh"
IN="$WS/data/cofold/hetero_ext/inputs"; OUT="$WS/data/cofold/hetero_ext/af3"; TODO="$WS/data/cofold/hetero_ext/todo_gpu$G"
mkdir -p "$OUT" "$WS/logs/cofold_hetero_ext"; rm -rf "$TODO"; mkdir -p "$TODO"
for f in "$IN"/*.json; do
  n=$(basename "$f" .json); n=${n#*_}
  [[ -s "$OUT/$n/${n}_ranking_scores.csv" ]] || ln -s "$f" "$TODO/"
done
source "$WS/scripts/af3_local_env.sh"   # node-local AF3 env/weights (stall mitigation)
export CUDA_VISIBLE_DEVICES=$G
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_autotune_level=0"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
python "$WS/scripts/run_af3_deadline.py" --input_dir="$TODO" --output_dir="$OUT" --model_dir="$AF3_MODEL_DIR" \
  --norun_data_pipeline --flash_attention_implementation=xla --jax_compilation_cache_dir="$AF3_JAX_CACHE" \
  --force_output_dir
