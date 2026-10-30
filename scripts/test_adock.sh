#!/usr/bin/env bash
# Smoke/timing test of adock: env check, 1STP (streptavidin/biotin) fixture from the README, docking on CPU
# (v2_default, v2_fast; GPU v2_fast/v2_fitted only with RUN_GPU=1, see below), then top-1 RMSD vs crystal.
# Outputs: results/tool_tests/adock/, logs/test_adock.log
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
source "$WS/scripts/env_threads.sh"
export PATH="$WS/conda_envs/dock/bin:$PATH"
export ADOCK_WORK="$WS/.cache/adock"                       # default would be /tmp/adock
export ADOCK_PDB_MIRROR="$ADOCK_WORK/pdb_mirror"           # minimal mirror with the 3 wwPDB files the tests use
export ADOCK_CPU_CACHE="$ADOCK_WORK/cpukernel"
export CUPY_CACHE_DIR="$WS/.cache/cupy"                    # default ~/.cupy
GPU_ID="${GPU_ID:-2}"                                       # GPU 0 fails scripts/check_gpus.py on this node
export CUDA_VISIBLE_DEVICES="$GPU_ID"
AD="$WS/tools/adock"; OUT="$WS/results/tool_tests/adock"; mkdir -p "$OUT" "$ADOCK_WORK"
cd "$AD"

FX="$ADOCK_WORK/smoke_test/1stp"
if [[ ! -s "$FX/box.json" ]]; then
  python scripts/check_env_dock.py > "$OUT/check_env.log" 2>&1 || true   # numba CUDA FAIL is expected (README)
  python scripts/fetch_test_mirror.py > "$OUT/fetch_mirror.log" 2>&1
  # writes $FX/{receptor.pdb,ligand_input.sdf,ligand_xtal.sdf,box.json}
  python scripts/smoke_test_baselines.py --prepare-only > "$OUT/prepare.log" 2>&1
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS_DOCK:-8}"   # CPU engine is OpenMP; env_threads.sh would make it 1 thread

# All dockings in ONE Python process: importing the dock env (torch, cupy, meeko, openmm) from the workspace
# filesystem took 1-5 min per process on 2026-10-29, which would dominate per-call CLI timings.
python - "$FX" "$OUT" <<'EOF1'
import sys, time, contextlib, io, json
from adock import dock
fx, out = sys.argv[1], sys.argv[2]
import os
# CPU engine (OpenMP). v2_fitted (final holo preset) has no CPU path (README).
runs = [("1stp_v2_default_cpu", "v2_default", "cpu"), ("1stp_v2_fast_cpu", "v2_fast", "cpu")]
# GPU engine: on this node (Quadro RTX 6000, sm_75) the CUDA MC kernel did not finish within 15-30 min
# on 2026-10-29 (README: developed on Blackwell sm_120). Opt in with RUN_GPU=1.
if os.environ.get("RUN_GPU") == "1":
    runs += [("1stp_v2_fast_gpu", "v2_fast", "cuda"), ("1stp_v2_fitted_gpu", "v2_fitted", "cuda")]
for name, mode, dev in runs:
    buf = io.StringIO(); t = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        dock.main(["--receptor", f"{fx}/receptor.pdb", "--ligand", f"{fx}/ligand_input.sdf", "--box", f"{fx}/box.json",
                   "--out", f"{out}/{name}.sdf", "--mode", mode, "--device", dev, "--seed", "0"])
    open(f"{out}/{name}.log", "w").write(buf.getvalue())
    print(f"{name:22s} wall {time.perf_counter() - t:.1f} s (in-process, imports excluded)", flush=True)
EOF1

# top-1 / best-of-9 symmetry-corrected heavy-atom RMSD vs the crystal pose (adock.rmsd.sym_rmsd)
python - "$OUT" "$FX/ligand_xtal.sdf" <<'EOF2' | tee "$OUT/rmsd.tsv"
import sys, json
from pathlib import Path
from rdkit import Chem
from adock.rmsd import sym_rmsd
out, xtal = Path(sys.argv[1]), sys.argv[2]
ref = Chem.MolFromMolFile(xtal)
print("run\ttop1_rmsd\tbest_rmsd\tn_poses\ttop_affinity\ttimings")
for f in sorted(out.glob("1stp_*.sdf")):
    poses = [m for m in Chem.SDMolSupplier(str(f)) if m is not None]
    r = sym_rmsd(ref, poses)
    last = [l for l in (out / (f.stem + ".log")).read_text().splitlines() if l.startswith("{")]
    info = json.loads(last[-1]) if last else {}
    aff = info.pop("top_affinity", None)
    print(f"{f.stem}\t{r[0]:.2f}\t{min(r):.2f}\t{len(poses)}\t{aff}\t{json.dumps(info)}")
EOF2
