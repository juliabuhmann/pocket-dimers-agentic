#!/usr/bin/env bash
# Copy the AF3 env, AF3 source, weights and JAX cache from the workspace (Lustre) to node-local /tmp. Resumable.
# (cp -a of the 6.8 GB env took ~15 min from Lustre.) A killed copy can leave one partial file: delete files newer than the start time first.
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}; L=/tmp/${USER}_af3
mkdir -p "$L/data/af3"
cp -a --update=none "$WS/conda_envs/af3" "$L/" 2>/dev/null || cp -a --update=none "$WS/conda_envs/af3/." "$L/af3/"
cp -a --update=none "$WS/tools/alphafold3" "$L/alphafold3_src"
cp -a --update=none "$WS/data/af3/models" "$L/data/af3/"
cp -a --update=none "$WS/data/af3/jax_cache" "$L/jax_cache"
