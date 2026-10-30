# Source to use the in-house tools from the workspace (after scripts/setup_tools.sh).
#   source scripts/env_tools.sh pockets   -> burrow + pocketcompare venv
#   source scripts/env_tools.sh dock      -> adock conda env
# Keeps every cache inside the workspace (defaults would write to $HOME or /tmp).
_WS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$_WS/scripts/env_threads.sh"                       # compute policy: 1 BLAS/OMP thread per process
export NUMBA_CACHE_DIR="$_WS/.cache/numba"
export TRITON_CACHE_DIR="$_WS/.cache/triton"               # default ~/.triton
export BURROW_MODELS="$_WS/tools/burrow/models"
export ADOCK_WORK="$_WS/.cache/adock"                      # default /tmp/adock
export ADOCK_CPU_CACHE="$_WS/.cache/adock/cpukernel"
export ADOCK_PDB_MIRROR="${ADOCK_PDB_MIRROR:-$_WS/.cache/adock/pdb_mirror}"
export CUPY_CACHE_DIR="$_WS/.cache/cupy"                   # default ~/.cupy
# GPU 0 of this node returns wrong results (scripts/check_gpus.py, 2026-10-29): never expose it.
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1,2,3}"
case "${1:-pockets}" in
  pockets) source "$_WS/.venvs/pockets/bin/activate" ;;
  dock)    export PATH="$_WS/conda_envs/dock/bin:$PATH" ;;
esac
unset _WS
