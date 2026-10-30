#!/usr/bin/env bash
# Install AlphaFold3 (google-deepmind, commit pinned below) and Boltz-2 into workspace uv venvs;
# decompress AF3 weights; fetch Boltz-2 assets. Single process, no parallel pools.
# Commands as actually run during setup on 2026-10-29; see docs/COFOLDING.md for verified status.
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
AF3_COMMIT=8f8ac383afba995f5953182d39a560ffb8bbbfd3
WEIGHTS_ZST=${AF3_WEIGHTS_ZST:-/cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst}
export UV_CACHE_DIR=$WS/.uv_cache
source "$WS/scripts/env_threads.sh"
cd "$WS"

# --- AF3 ---
[ -d tools/alphafold3 ] || git clone https://github.com/google-deepmind/alphafold3 tools/alphafold3
git -C tools/alphafold3 checkout "$AF3_COMMIT"
uv venv --python 3.12 conda_envs/af3
( source conda_envs/af3/bin/activate
  cd tools/alphafold3
  uv pip install .          # builds C++ ext; pulls jax 0.10.2 + CUDA 12.9 wheels (~10 min)
  build_data                # writes CCD msgpack into the installed package
)
mkdir -p data/af3/models data/af3/jax_cache
[ -f data/af3/models/af3.bin ] || conda_envs/af3/bin/python -c "
import zstandard
with open('$WEIGHTS_ZST','rb') as f, open('data/af3/models/af3.bin','wb') as o:
    zstandard.ZstdDecompressor().copy_stream(f,o)"

# --- Boltz-2 (fallback / comparison) ---
uv venv --python 3.12 conda_envs/boltz
uv pip install --python conda_envs/boltz/bin/python boltz
scripts/download_boltz_assets.sh     # weights/CCD from HuggingFace into data/boltz/cache (~4 GB; mols.tar extraction is slow)

# --- sanity: GPU numerics (GPU 0 on this node FAILED this test, see docs/COFOLDING.md) ---
for G in 0 1 2 3; do echo "GPU $G"; CUDA_VISIBLE_DEVICES=$G conda_envs/af3/bin/python scripts/test_gpu_basic.py | head -3; done
