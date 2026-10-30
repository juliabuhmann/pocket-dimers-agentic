#!/usr/bin/env bash
# Clone and install the three in-house tools (burrow, pocketcompare, adock) into the workspace.
# Re-runnable: existing clones are fetched and checked out at the pinned commit.
#
# Usage: bash scripts/setup_tools.sh [burrow_pc|adock|all]   (default: all)
#
# Layout (all under the workspace, nothing under $HOME):
#   tools/{burrow,pocketcompare,adock}   git clones at the commits pinned below
#   .venvs/pockets                       uv venv (Python 3.12): burrow + pocketcompare (+gpu, +search extras)
#   conda_envs/dock                      conda prefix env for adock (from tools/adock/envs/dock.yml)
#   .cache/{uv,conda_pkgs,numba,triton,cupy,adock}  caches (see scripts/env_tools.sh)
# Afterwards: source scripts/env_tools.sh pockets|dock ; python scripts/check_gpus.py (GPU 0 of this node is faulty)
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="$WS/tools"
WHAT="${1:-all}"

# Pinned commits (recorded 2026-10-29 from the default branch 'main' of each repository).
BURROW_URL="https://github.com/zoidpuff/burrow";             BURROW_COMMIT="b2da7abfbbc56d4be4850dfd44e02e55c7b5a648"
POCKETCOMPARE_URL="https://github.com/zoidpuff/pocketcompare"; POCKETCOMPARE_COMMIT="efcdf3ebb47bd8e9b8fe431aad51d41d466f4b76"
ADOCK_URL="https://github.com/juliabuhmann/adock";           ADOCK_COMMIT="a7f28535f4bbed06d1c5e7f09a5b942c84835d08"

PY_VERSION="3.12"   # burrow/pocketcompare were developed and pinned on Python 3.12.3
export UV_CACHE_DIR="$WS/.cache/uv"
export UV_PYTHON_INSTALL_DIR="$WS/.cache/uv-python"
export CONDA_PKGS_DIRS="$WS/.cache/conda_pkgs"
mkdir -p "$TOOLS" "$WS/.venvs" "$WS/conda_envs" "$UV_CACHE_DIR" "$CONDA_PKGS_DIRS"

clone_at() {  # url dir commit
  local url="$1" dir="$2" commit="$3"
  if [[ ! -d "$dir/.git" ]]; then git clone -q "$url" "$dir"; fi
  git -C "$dir" fetch -q origin
  git -C "$dir" checkout -q "$commit"
  echo "$(basename "$dir") @ $(git -C "$dir" rev-parse HEAD)"
}

if [[ "$WHAT" == all || "$WHAT" == burrow_pc ]]; then
  clone_at "$BURROW_URL" "$TOOLS/burrow" "$BURROW_COMMIT"
  clone_at "$POCKETCOMPARE_URL" "$TOOLS/pocketcompare" "$POCKETCOMPARE_COMMIT"
  VENV="$WS/.venvs/pockets"
  [[ -x "$VENV/bin/python" ]] || uv venv --python "$PY_VERSION" "$VENV"
  # burrow first (not on PyPI; pocketcompare depends on burrow-pockets==1.0.0)
  uv pip install --python "$VENV/bin/python" -e "$TOOLS/burrow[test,gpu]"
  uv pip install --python "$VENV/bin/python" -e "$TOOLS/pocketcompare[search,test]"
  # analysis helpers used by our wrapper scripts (not needed by the tools themselves)
  uv pip install --python "$VENV/bin/python" gemmi pyarrow
  mkdir -p "$WS/envs"; uv pip freeze --python "$VENV/bin/python" > "$WS/envs/pockets_freeze.txt"   # record resolved versions
  echo "pockets venv ready: source scripts/env_tools.sh pockets"
fi

if [[ "$WHAT" == all || "$WHAT" == adock ]]; then
  clone_at "$ADOCK_URL" "$TOOLS/adock" "$ADOCK_COMMIT"
  PREFIX="$WS/conda_envs/dock"
  if [[ ! -x "$PREFIX/bin/python" ]]; then
    # Same as tools/adock/scripts/install_env_dock.sh --from-yml, but with the prefix under the workspace.
    conda env create -y -p "$PREFIX" -f "$TOOLS/adock/envs/dock.yml"
  fi
  "$PREFIX/bin/python" -m pip install --no-deps -e "$TOOLS/adock"
  echo "adock env ready: source scripts/env_tools.sh dock"
fi
