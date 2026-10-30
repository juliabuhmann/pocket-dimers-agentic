# Co-folding (protein dimer + ligand) on this node

Status written 2026-10-29. All numbers come from the logs in `logs/` (single runs, n=1, not replicated).
Structure quality was only checked via the model's own confidence scores, not against experimental structures.

## Summary
- AF3 (google-deepmind, commit `8f8ac38`) installs and runs on the Quadro RTX 6000 (Turing, cc 7.5) with the settings below. Predictions are numerically sane on GPUs 1-3.
- **GPU 0 is broken**: host-to-device round trip returns corrupted data (`scripts/test_gpu_basic.py`: `roundtrip False`, `x+1` max error 0.5). On GPU 0 AF3 produced all-NaN coordinates (logs `af3_tpi1_pga_try2..5*`) after long debugging. Use GPU 1-3 only. Please report GPU 0 to the admins.
- Boltz-2 also installs and runs (GPU 1), with precomputed MSAs given as a3m.
- Existing envs in `/path/to/conda/envs`: no standalone alphafold3 or boltz env (intellifold-host and openfold3 contain a vendored `alphafold3` 3.0.2 lib from IntelliFold, not used here).

## Install (exact commands: `scripts/setup_cofold.sh`)
```
bash scripts/setup_cofold.sh          # AF3 + Boltz envs under conda_envs/{af3,boltz}, weights, assets
```
Key points: `uv venv --python 3.12 conda_envs/af3; cd tools/alphafold3; uv pip install .; build_data`
(needs internet; ~10 min). Weights are decompressed from `/cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst`
to `data/af3/models/af3.bin` (1.1 GB; with the env's python, system python has no `zstandard`).
Boltz: `uv pip install boltz` (torch 2.15+cu132 wheel, works with driver 580). Boltz's own downloader hung
(model-gateway.boltz.bio), so `scripts/download_boltz_assets.sh` fetches weights/CCD from HuggingFace.

## Turing (no bf16) settings for AF3 (`scripts/run_af3.sh`)
- `--flash_attention_implementation=xla` (triton/cudnn kernels need Ampere+; tokamax logs `NotImplementedError: Not supported on Quadro RTX 6000` and falls back, harmless).
- `XLA_FLAGS=--xla_disable_hlo_passes=custom-kernel-fusion-rewriter` (AF3 docs: required for CUDA cc 7.x; `docs/known_issues.md`).
- `XLA_PYTHON_CLIENT_PREALLOCATE=false`, `--norun_data_pipeline`, JAX compilation cache in `data/af3/jax_cache`.
- bf16 weights run as-is (XLA emulates); forcing float32 (`AF3_BFLOAT16=none` in `scripts/run_af3_wrapper.py`) was not needed once a healthy GPU was used (the float32 test was on broken GPU 0, so it is uninformative).
- Do NOT set `--xla_gpu_autotune_level=0` on a healthy GPU; not needed. (Default autotuning failed with a cuBLASLt redzone error only on broken GPU 0.)

## Inputs from cached MSAs
`/cluster/project/beltrao/shared/25.05_alphafold3_msas_human/<acc>_data.json.gz` holds one protein per file with
`unpairedMsa`, `pairedMsa` (present for the files checked, checked 24 accessions, all non-empty) and `templates`.
```
python scripts/make_af3_input.py --acc P60174 P60174 --ligand-ccd PGA --out data/af3/inputs/tpi1_pga.json
scripts/run_af3.sh data/af3/inputs/tpi1_pga.json results/cofold/af3_tpi1_pga 1     # GPU id 1
python scripts/make_boltz_input.py --acc P60174 P60174 --ligand-ccd PGA --outdir data/boltz/inputs/tpi1_pga --name tpi1_pga
```
`make_af3_input.py` checks the cached sequence against the UniProt REST sequence (aborts on mismatch), merges identical
accessions into one entity (ids A,B), supports `--ligand-ccd` / `--ligand-smiles`, `--seeds`, `--no-templates`.
Caveats: the cached pairedMsa comes from a monomer search (UniProt hits), AF3 pairs by species at featurisation; for
a heterodimer no true cross-chain pairing beyond that is computed. Cached MSA queries may cover only a sub-sequence
(then the sequence check fails against UniProt; use `--skip-uniprot-check` knowingly). Boltz input uses only unpairedMsa.

## Test case
Human TPI1 (P60174, 249 aa) homodimer + PGA (2-phosphoglycolate), 1 seed, 5 diffusion samples, 10 recycles (AF3 defaults).
AF3: ranking_score 0.95, iptm 0.95, no clash (`results/cofold/af3_tpi1_pga/`). Ligand chain ptm 0.03 is expected for a tiny ligand.
Boltz-2: confidence 0.96, ligand_iptm 0.95 (`results/cofold/boltz_tpi1_pga/`), 1 diffusion sample, 3 recycles (boltz defaults).
Not compared against a reference structure.

## Runtime table (GPU 1, one process, homodimer + ATP unless noted, AF3 1 seed = 5 samples)
GPU memory = max of `nvidia-smi` memory.used sampled every 3-5 s (allocator reserve, upper bound of need).
"Inference" includes XLA compilation of that token bucket (not separable here; warm JAX cache saved ~40 s at 507 tokens).
| Engine | Case | Tokens (padded bucket) | Inference | Wall | Peak GPU | Result |
|---|---|---|---|---|---|---|
| AF3 | TPI1 x2 + PGA | 507 (512) | 174 s (135 s warm cache) | ~190 s | 6.6 GB | OK |
| AF3 | P35557 x2 + ATP | 961 (1024) | 710 s | ~750 s | 12.4 GB | OK |
| AF3 | P53350 x2 + ATP | 1237 (1280) | 1330 s | ~1380 s | 18.2 GB | OK |
| AF3 | P07900 x2 + ATP | 1495 (1536) | - | - | OOM (16 GiB single alloc) | FAIL |
| AF3 + unified memory | P07900 x2 + ATP | 1495 (1536) | 2093 s | 2752 s | ~18 GB + host spill | OK |
| AF3 | P10721 x2, P00533 x2, P35968 x2 + ATP | 1983 / 2451 / 2743 | - | - | OOM | FAIL (unified memory not tried) |
| Boltz-2 | TPI1 x2 + PGA | 507 | predict 142 s | 922-1016 s | 11.3 GB | OK |

Unified memory: `TF_FORCE_UNIFIED_MEMORY=true XLA_CLIENT_MEM_FRACTION=3.2` (AF3 docs, performance.md).
Boltz wall time is dominated by ~13 min startup (reading the 45k-file `mols/` dir on Lustre; I/O wait), paid once per `boltz predict` call, so batch many YAMLs in one call (pass a directory) or copy `mols/` to node-local disk.
Max AF3 size without unified memory: bucket 1280 (<=~1280 tokens incl. ligand atoms); with unified memory 1536 works but is slow. Max sizes for Boltz not tested.
Ranking scores of the larger runs (0.36-0.83) are for arbitrary protein pairs (not known dimers); not a quality benchmark.

## Threads / compute policy
Policy: `docs/COMPUTE_POLICY.md`; every run is in `logs/compute_ledger.txt`. One GPU process at a time, `scripts/env_threads.sh` sourced.
Measured `ps -L -p PID | wc -l` (AF3 python, during inference, after sourcing env_threads.sh): **253 threads**.
Boltz-2 (`--num_workers 0`): ~5 threads sampled during prediction (sampling only the main process; startup 1 thread).
Note `ps -L` on the `run_af3.sh` wrapper PID shows 1, use the `python` child.

## Limitations
- GPU 0 faulty; avoid. Always run `scripts/test_gpu_basic.py` first on a new node/GPU.
- XLA attention (no flash) makes AF3 slow/memory hungry: ~O(N^2); >1280 tokens needs unified memory and is ~35 min+ per seed.
- Per-bucket JIT compile is repeated per process (cache helps only modestly).
- Cached MSAs are human-only, monomer-derived; no cross-species or true paired MSAs for heterodimers; templates from the cached files (date cutoff of 25.05 build).
- Ligand tokens count per heavy atom (ATP = 31), so large ligands eat token budget.
- Single runs only; no timing variance, no accuracy validation against PDB.
- AF3 weights licence: AF3 terms of use apply (output TERMS_OF_USE.md is written with each run).

## Files
`scripts/{setup_cofold.sh,make_af3_input.py,make_boltz_input.py,run_af3.sh,run_af3_wrapper.py,sweep_af3_size.sh,download_boltz_assets.sh,test_gpu_basic.py,test_gpu_numerics.py,debug_af3_nan.py}`;
logs: `logs/af3_*.log`, `logs/boltz_*.log`, `logs/af3_size_sweep.tsv`, failed attempts `logs/af3_tpi1_pga_try*`.
