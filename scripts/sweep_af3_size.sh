#!/usr/bin/env bash
# Size sweep: homodimer + ATP, 1 seed, one GPU, sequential. usage: sweep_af3_size.sh GPU ACC [ACC...]
# Appends to logs/af3_size_sweep.tsv: acc  total_tokens_est  status  inference_s  peak_gpu_MiB
set -uo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
GPU=$1; shift
source "$WS/scripts/env_threads.sh"
for ACC in "$@"; do
  L=$(echo $ACC | tr A-Z a-z)
  python3 "$WS/scripts/make_af3_input.py" --acc $ACC $ACC --ligand-ccd ATP --out "$WS/data/af3/inputs/${L}_x2_ATP.json" || continue
  OUT="$WS/results/cofold/sweep_${L}"; LOG="$WS/logs/af3_sweep_${L}.log"; rm -rf "$OUT"
  "$WS/scripts/run_af3.sh" "$WS/data/af3/inputs/${L}_x2_ATP.json" "$OUT" $GPU > "$LOG" 2>&1 &
  PID=$!
  PEAK=0
  while kill -0 $PID 2>/dev/null; do
    M=$(nvidia-smi -i $GPU --query-gpu=memory.used --format=csv,noheader,nounits); [ "$M" -gt "$PEAK" ] && PEAK=$M; sleep 3
  done
  wait $PID; RC=$?
  T=$(grep -oP "Running model inference with seed \d+ took \K[0-9.]+" "$LOG" | tail -1)
  S=$([ $RC -eq 0 ] && echo OK || echo "FAIL(rc=$RC)")
  TOK=$(grep -oP "protein tokens=\K[0-9]+" "$LOG" || true)
  printf "%s\t%s\t%s\t%s\n" "$ACC" "$S" "${T:-NA}" "$PEAK" >> "$WS/logs/af3_size_sweep.tsv"
done
