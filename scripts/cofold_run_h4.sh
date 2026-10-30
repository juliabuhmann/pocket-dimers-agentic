#!/usr/bin/env bash
# Run AF3 for deep dive #2 (results/cofold/h4/cases.tsv; inputs data/cofold/h4/inputs, written by cofold_select_h4.py)
# in ONE process on ONE GPU, looping over the input dir in file-name (= priority) order; finished jobs are skipped.
# usage: cofold_run_h4.sh GPU      (this deep dive was assigned GPU 1; GPU 0 is faulty, GPUs 2-3 belong to deep dive #1)
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
G=${1:?gpu}
[[ "$G" == 0 ]] && { echo "GPU 0 is faulty, refusing"; exit 1; }
source "$WS/scripts/env_threads.sh"
IN="$WS/data/cofold/h4/inputs"; OUT="$WS/data/cofold/h4/af3"; TODO="$WS/data/cofold/h4/todo_gpu$G"
mkdir -p "$OUT" "$WS/logs/cofold_h4"; rm -rf "$TODO"; mkdir -p "$TODO"
for f in "$IN"/*.json; do   # skip jobs with an existing ranking_scores.csv
  n=$(basename "$f" .json); n=${n#*_}
  [[ -s "$OUT/$n/${n}_ranking_scores.csv" ]] || ln -s "$f" "$TODO/"
done
source "$WS/conda_envs/af3/bin/activate"
export CUDA_VISIBLE_DEVICES=$G
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_autotune_level=0"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export AF3_SRC="$WS/tools/alphafold3"
python "$WS/scripts/run_af3_wrapper.py" --input_dir="$TODO" --output_dir="$OUT" --model_dir="$WS/data/af3/models" \
  --norun_data_pipeline --flash_attention_implementation=xla --jax_compilation_cache_dir="$WS/data/af3/jax_cache" \
  --force_output_dir
