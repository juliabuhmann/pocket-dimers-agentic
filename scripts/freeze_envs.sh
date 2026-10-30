#!/usr/bin/env bash
# Writes envs/<name>.txt for every env in .venvs/* and conda_envs/*; skips an env after 120 s.
cd "$(dirname "$0")/.."
mkdir -p envs
for d in .venvs/* conda_envs/*; do
  n=$(basename "$d"); [ "$(dirname "$d")" = conda_envs ] && n="conda_$n"
  if [ -d "$d/conda-meta" ]; then
    timeout 120 conda list -p "$d" --explicit > "envs/$n.txt" 2> "envs/$n.err" || echo "SKIPPED/FAILED $d" >> envs/FREEZE_NOTES.txt
  else
    timeout 120 uv pip freeze --python "$d/bin/python" > "envs/$n.txt" 2> "envs/$n.err" || timeout 120 "$d/bin/python" -m pip freeze > "envs/$n.txt" 2> "envs/$n.err" || echo "SKIPPED/FAILED $d" >> envs/FREEZE_NOTES.txt
  fi
  [ -s "envs/$n.err" ] || rm -f "envs/$n.err"
done
