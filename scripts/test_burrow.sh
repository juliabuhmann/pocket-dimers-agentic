#!/usr/bin/env bash
# Smoke/timing test of burrow on test inputs (run scripts/fetch_tool_test_inputs.py first).
# Outputs: results/tool_tests/burrow/<case>/  (burrow CSVs + timing.csv), logs/test_burrow.log
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
source "$WS/scripts/env_threads.sh"
source "$WS/.venvs/pockets/bin/activate"
export BURROW_MODELS="$WS/tools/burrow/models" NUMBA_CACHE_DIR="$WS/.cache/numba" TRITON_CACHE_DIR="$WS/.cache/triton"
IN="$WS/data/interim/tool_tests"; OUT="$WS/results/tool_tests/burrow"
GPU="${GPU:-cuda:2}"          # GPU used for the GPU runs (cuda:0 on this node FAILS scripts/check_gpus.py)
SECTION="${SECTION:-all}"     # all | gpu_forced
WORKERS="${WORKERS:-8}"       # compute policy: <= 8 workers for tests
mkdir -p "$OUT"

run() {  # name, then burrow args
  local name="$1"; shift
  local t0=$(date +%s%N)
  burrow predict "$@" --out "$OUT/$name" > "$OUT/$name.log" 2>&1
  printf "%-28s wall %d ms\n" "$name" $(( ($(date +%s%N) - t0) / 1000000 ))
}

if [[ "$SECTION" == all ]]; then
# single structures, CPU, 1 worker (wall time includes Python start-up + Numba cache load)
run afdb_monomer_cpu  "$IN/afdb/AF-P31946-F1-model_v6.cif" --af-mode --workers 1
run homodimer_cpu     "$IN/dimers/AF-0000000065760335.pdb" --af-mode --workers 1
run heterodimer_cpu   "$IN/dimers/AF-0000000203452075.pdb" --af-mode --workers 1
run pdb_raw_cpu       "$IN/pdb/1A28.cif" --workers 1          # experimental: NO --af-mode
run pdb_stripped_cpu  "$IN/pdb/1A28_protein.cif" --workers 1
# same on GPU
run homodimer_gpu     "$IN/dimers/AF-0000000065760335.pdb" --af-mode --workers 1 --device "$GPU"
# throughput: 100 dimers (50 homo + 50 hetero)
run batch100_cpu_w8   "@$IN/dimer_batch.txt" --af-mode --workers "$WORKERS" --no-points
run batch100_gpu_w8   "@$IN/dimer_batch.txt" --af-mode --workers "$WORKERS" --no-points --device "$GPU"
# throughput on 1000 dimers (500 homo + 500 hetero); BURROW_FAST_PROFILE records whether the GPU path was used
export BURROW_FAST_PROFILE="$OUT/batch1000_gpu_w8.profile.jsonl"; rm -f "$BURROW_FAST_PROFILE"
run batch1000_gpu_w8  "@$IN/dimer_batch_large.txt" --af-mode --workers "$WORKERS" --no-points --device "$GPU"
unset BURROW_FAST_PROFILE
run batch1000_cpu_w8  "@$IN/dimer_batch_large.txt" --af-mode --workers "$WORKERS" --no-points
fi
# Forced GPU path: by default burrow runs a hybrid scheduler (CPU workers process structures while the GPU
# initialises, ~30-80 s); with few structures the GPU never gets work. BURROW_FAST_HYBRID=0 sends all to the GPU.
export BURROW_FAST_HYBRID=0
export BURROW_FAST_PROFILE="$OUT/batch1000_gpuforced_w8.profile.jsonl"; rm -f "$BURROW_FAST_PROFILE"
run batch1000_gpuforced_w8 "@$IN/dimer_batch_large.txt" --af-mode --workers "$WORKERS" --no-points --device "$GPU"
export BURROW_FAST_PROFILE="$OUT/homodimer_gpuforced.profile.jsonl"; rm -f "$BURROW_FAST_PROFILE"
run homodimer_gpuforced    "$IN/dimers/AF-0000000065760335.pdb" --af-mode --workers 1 --device "$GPU"
unset BURROW_FAST_PROFILE BURROW_FAST_HYBRID
