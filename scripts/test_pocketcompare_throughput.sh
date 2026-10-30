#!/usr/bin/env bash
# Throughput of `pocketcompare describe` (BURROW + fingerprints) on 100 predicted dimers (dimer_batch.txt).
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
source "$WS/scripts/env_threads.sh"; source "$WS/.venvs/pockets/bin/activate"
export NUMBA_CACHE_DIR="$WS/.cache/numba" TRITON_CACHE_DIR="$WS/.cache/triton"
OUT="$WS/results/tool_tests/pocketcompare"; mkdir -p "$OUT"
LIST="$WS/data/interim/tool_tests/dimer_batch.txt"
WORKERS="${WORKERS:-8}"; GPU="${GPU:-cuda:2}"   # cuda:0 fails scripts/check_gpus.py on this node
mapfile -t FILES < "$LIST"
run() { local name="$1"; shift; local t0=$(date +%s%N)
  pocketcompare describe "${FILES[@]}" --af-mode --workers "$WORKERS" -o "$OUT/$name.parquet" "$@" > "$OUT/$name.log" 2>&1
  printf "%-24s wall %d ms\n" "$name" $(( ($(date +%s%N) - t0) / 1000000 )); }
RUNS="${RUNS:-fast gpu}"
[[ " $RUNS " == *" fast "* ]] && run dimers100_fast_cpu_w8 --backend fast
[[ " $RUNS " == *" gpu "* ]] && run dimers100_gpu_w8      --device "$GPU"
true
