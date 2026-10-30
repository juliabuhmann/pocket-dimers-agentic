# In-house tools: burrow, pocketcompare, adock

Status: installed and smoke-tested on 2026-10-29 on this node (64 CPUs, 4x Quadro RTX 6000 / Turing sm_75,
shared with other agents; load average 20-70 during the tests). All numbers below were measured here unless they
are marked "(README)". Timings on a shared node vary by about ±30 % or more. None of the three tools has been
validated by a human; their READMEs say they were built by AI agents.

Reproduce:
```bash
bash scripts/setup_tools.sh                    # clone at pinned commits, create .venvs/pockets and conda_envs/dock
.venvs/pockets/bin/python scripts/check_gpus.py   # GPU sanity check (GPU 0 fails, see Gotchas)
.venvs/pockets/bin/python scripts/fetch_tool_test_inputs.py   # test structures -> data/interim/tool_tests/
bash scripts/test_burrow.sh                    # -> results/tool_tests/burrow/, logs/test_burrow.log
.venvs/pockets/bin/python scripts/test_pocketcompare.py        # -> results/tool_tests/pocketcompare/report.json
bash scripts/test_pocketcompare_throughput.sh  # 100 dimers, CPU fast backend vs GPU
bash scripts/test_adock.sh                     # 1STP fixture docking -> results/tool_tests/adock/
.venvs/pockets/bin/python scripts/summarize_pocket_chains.py <burrow_out_dir|describe.parquet> -o chains.tsv
source scripts/env_tools.sh pockets|dock       # use the tools interactively (sets caches, thread caps, GPUs)
```

## 1. Versions

| tool | repository | commit (pinned in `scripts/setup_tools.sh`) | package version | env |
|---|---|---|---|---|
| burrow | github.com/zoidpuff/burrow (main) | `b2da7abfbbc56d4be4850dfd44e02e55c7b5a648` (code = release commit d4ed36c, 2026-10-27; later commits are docs) | burrow-pockets 1.0.0 | `.venvs/pockets` (uv, Python 3.12.3; torch 2.6.0+cu124, triton 3.2.0) |
| pocketcompare | github.com/zoidpuff/pocketcompare (main) | `efcdf3ebb47bd8e9b8fe431aad51d41d466f4b76` (code = 0ca2855, 2026-10-28; later commits are docs) | 1.0.0, method `pocketcompare-final:A4b:n1_linear\|raw` | same venv (+faiss-cpu 1.14.3) |
| adock | github.com/juliabuhmann/adock (main) | `a7f28535f4bbed06d1c5e7f09a5b942c84835d08` | adock 0.1.0 | `conda_envs/dock` (from `tools/adock/envs/dock.yml`: Python 3.11, torch 2.8.0+cu128, cupy-cuda12x 14.2.0, CUDA 12.9) |

Resolved versions of the pockets venv: `envs/pockets_freeze.txt`.

### Which versions built the AFDB pocket atlas (`/path/to/pocketcompare_afdb_atlas/`)

- The atlas parquet files store a method record (schema metadata key `pocketcompare`):
  `pocketcompare_version 1.0.0`, `config_sha256 b3ea8fe9…a36a0c`, `af_mode true`, `unit merged`, `min_points 10`,
  `backend fast`, `device cuda:N`, `storage float16`. Our checkout's `configs/final_method.yaml` has the **same
  sha256**, and the BURROW model shipped inside pocketcompare (`models/burrow/burrow_all_default_s0.json`) is
  byte-identical to burrow's (sha256 `f33963e7…c4315`).
- Atlas filters (metadata key `pocketcompare_afdb`): merged pockets with `burrow_score >= 0.03` and lining
  `plddt_mean >= 70` were kept.
- No git commit is recorded in the atlas, and the pipeline code (`afdb/worker.py`, `afdb.sbatch`) and its env
  (`<atlas authors' private conda env>`) are not readable from here. The run
  started 2026-10-28 15:19, after the last code commit of both repositories.
- **Checked directly:** re-describing `AF-P31946-F1-model_v6` with our install gave the same 2 pocket ids, identical
  lining residues, and fingerprints within 4.3e-5 of the stored float16 values (float16 rounding).
  `scripts/test_pocketcompare.py`, `results/tool_tests/pocketcompare/report.json` (`atlas_check`).
- Atlas status (from `state/done/*.json`): tier 1 (46 reference proteomes, 555,846 models, 1,788,750 pockets) and
  tier 2 (Swiss-Prot rest, 436,886 models, 1,442,678 pockets) complete; tier 3 (rest of AFDB, 48,341 units) not
  started. No FAISS index is present (`index/` missing). Human (taxid 9606) has 90,560 kept pockets. The atlas is
  on scratch (15-day purge).

**To run identically on dimers:** `pocketcompare describe <files> --af-mode` (defaults: `--unit merged
--min-points 10`, no site, no score cut-off), then keep rows with `burrow_score >= 0.03` and `plddt_mean >= 70`.
Backend/device do not change the output (checked: CPU fast vs GPU, max |Δ| = 0 on 100 dimers).

## 2. burrow (pocket detection)

**What it does.** Fills the solvent shell with 1 Å grid points, describes them (ray-cast buriedness, chemistry),
scores points and pockets with two LightGBM models, outputs ranked pockets. No neural nets. CPU (Numba) and GPU
(Triton) backends with identical output (checked here: CPU vs GPU files identical on 1,000 dimers).

**CLI.**
```bash
source scripts/env_tools.sh pockets
burrow predict model.cif.gz --out out/ --af-mode                    # predicted model (pLDDT in B-factor column)
burrow predict 1abc.cif --out out/                                  # experimental structure: no --af-mode
burrow predict @list.txt --out out/ --af-mode --workers 8 --no-points
burrow predict @list.txt --out out/ --af-mode --workers 8 --device cuda:2   # one process per GPU
```
Options: `--preset default|fast|accurate`, `--ranker all|small|geometry` (`small` = non-ion ligands),
`--min-score none` (exhaustive, up to 50 pockets), `--output-rule v2|region`, `--af-ranker models/afr_final_afpdb_pl.json`.

**Inputs.** PDB or mmCIF (`.gz` ok), first model only, author chain ids and numbering. Any number of protein chains.
Only protein residues are used: ligands, waters, ions and nucleic acids are ignored. **Ligand stripping is not
needed**: raw 1A28 mmCIF (2 chains + progesterone + waters) and a gemmi-stripped copy gave byte-identical output
(rank 1 and 2 = the progesterone sites of chains B and A, scores 0.62 / 0.60).

**Outputs** (per input `<stem>`): `<stem>.csv` with `rank,score,center_x,center_y,center_z,residues`
(`residues` = space-separated `chain:resnum[icode]`); `<stem>.points.csv` (`rank,x,y,z` grid points, skip with
`--no-points`); `timing.csv` (`stem,status,n_atoms,n_pockets,total_ms,parse_ms,features_ms,point_score_ms,
pockets_ms,rank_ms,device,note`; the header line is truncated, so pass column names when reading it).
`score` = ranker probability-like score in 0-1, "comparable across proteins" (README) but not calibrated.
Default list: pockets with score < 0.03 dropped once 3 are listed. Python API with all intermediate data (points,
76 pocket descriptors incl. `n_chains`): `pocketdetect.burrow.export.pocket_representations(path, af_mode=True)`.

**Multi-chain behaviour.**
- Runs natively on whole complexes: grid, features and ranking are computed on the assembly, so interface
  pockets are found as single pockets. The predicted human dimers (PDB format, chains A and B, pLDDT in
  B-factors) run without any preprocessing.
- Chain attribution is only implicit: the `residues` column carries the chain of each lining residue. There is no
  chain column in the CSV. `scripts/summarize_pocket_chains.py` adds `chains`, `n_chains`, `res_per_chain`,
  `min_chain_frac` and `n_chains_ge3` (chains contributing ≥ 3 residues; the interface rule of burrow's dimer
  benchmark). The ranker uses the same `n_chains` descriptor internally.
- 1,000 human dimers (first 500 homo + 500 hetero by file name), default output: 13,321 pockets (13.3 per dimer);
  29.8 % have lining residues on both chains, 19.7 % have ≥ 3 residues on each chain (homodimers 20.5 %,
  heterodimers 18.7 %); among rank-1 pockets 24.4 % (homo) and 18.0 % (hetero) are interface pockets by the ≥ 3
  rule (`results/tool_tests/burrow/batch1000_chains.tsv`).
- Benchmarks (README / docs/REPORT.md, not re-run here): on 400 PDB assemblies BURROW is far ahead of P2Rank on
  single-chain sites in assemblies but **not better on interface sites** (−3.7 points, CI includes 0;
  hetero-interface −6.1, n = 66). Running chains in isolation still recovers many interface sites (interface DCC12
  69.6 vs 84.0 on full assemblies by a different protocol; see docs/DIMER_BENCHMARK.md).
- Homodimers: symmetric copies of a site appear as two separate pockets (one per chain or per interface copy);
  collapse them yourself (burrow's benchmark used residue-set Jaccard ≥ 0.5 after mapping chains to entities).

**Measured runtimes** (AF mode, default preset, `scripts/test_burrow.sh`):

| case | result |
|---|---|
| single structure, 1 CPU worker, wall incl. start-up | 16.5-18.5 s (start-up + per-worker warm-up ≈ 15 s; first-ever call compiles Numba ≈ 80 s) |
| warm, per structure, 1 worker | 150-165 ms for the single-chain examples (1a28_A, AF-P00918); 345-560 ms for dimers; the first structure of each worker takes 2.5-3.2 s (warm-up) |
| 1,000 human dimers (mean 5,159 heavy atoms), CPU, 8 workers | 82 s wall incl. start-up (12.3 structures/s); per structure median 449 ms, p90 961 ms, max 2.7 s (≈ 1.9 structures/s per worker) |
| same, GPU (`--device cuda:2`, 8 prep workers, `BURROW_FAST_HYBRID=0`) | 84 s wall, of which ≈ 48 s start-up (torch import, Triton warm-up); then 1,000 structures in 28.5 s (≈ 35/s, limited by CPU prep: GPU featurize 17 s, waiting for prep 10 s) |
| same, GPU with default hybrid scheduler | the GPU never received work: CPU workers finished all 1,000 during the 70-80 s GPU initialisation |

Extrapolation (not measured): the 10,614 human dimers (4,821 homo + 5,793 hetero) take ≈ 15 min with 8 CPU
workers or ≈ 6 min on one GPU with 8 prep workers.

**Threads per process** (`ps -L`): main process 4-5 threads and workers 1 thread each with
`scripts/env_threads.sh`; **without it the main process starts 64 threads** (BLAS/OpenMP pools). GPU runner: 5.

## 3. pocketcompare (pocket similarity)

**What it does.** Runs BURROW itself (it does not read burrow's CSV), merges touching region sub-pockets, and
describes each merged pocket's buried core (all buried grid points within 8 Å of the buriedness-weighted centre,
"sph8") by a 627-d descriptor (typed cavity histograms, pair-distance distributions, shape moments, lining
chemistry), mapped by a frozen linear head to a **384-d unit vector**. No sequence, fold or ligand input.

**CLI.**
```bash
pocketcompare describe files/ --af-mode --workers 8 --backend fast -o pockets.parquet   # all pockets
pocketcompare describe x.cif --residues A:43,B:45 -o site.parquet            # a given site (any chains)
pocketcompare describe x.cif --ligand lig.sdf --locate burrow -o site.parquet # BURROW pocket at a known ligand
pocketcompare compare a.parquet b.parquet --top 5 [--matrix S.npy]
pocketcompare index pockets.parquet -o idx/ [--kind ivfpq]                    # FAISS
pocketcompare search q.parquet idx/ -k 50 [--exclude-self] -o hits.tsv
```
`--chains A,B` restricts chains. Python: `pocketcompare.method.describe_all(path, af_mode=True)` returns
(vectors, table) but **without `pocket_id` and pLDDT columns** (the CLI adds them); read CLI output with
`pocketcompare.fpindex.read_table(path)` → (metadata, float32 (n, 384), method record).

**Output.** One row per merged pocket: `pocket_id` (`<stem>:m<k>`, k = merged index; ids are not contiguous after
the 10-point filter), `burrow_rank`, `burrow_score`, `members`, `residues` (`chain:resnum`), `center_*`,
`volume_A3`, `site_volume_A3`, `sph8_fallback`, `ok`, `plddt_mean/min/frac_lt70` (AF mode), `fingerprint`;
plus `<out>.structures.parquet` (status, timings). No score cut-off is applied.

**Score.** Inner product of two fingerprints = mean cosine of 3 heads, range −1 to 1 (self = 1). Higher = more
similar; trained to rank pockets binding similar ligands in unrelated proteins. Not calibrated, no validated
cut-off. Background (README, 116 k human AF pockets): random pairs median 0.11, 99th pct 0.62, 99.9th 0.77;
nearest neighbour median 0.82. Validated ROC-AUC ≈ 0.70-0.75 (TOUGH-M1 / temporal set); no gain on ProSPECCTs;
**PPI-interface pockets are recognised poorly across families** (README limitation). Use ranks, per-query
z-scores, or size-matched backgrounds rather than raw thresholds.

**Measured runtimes.**

| step | result |
|---|---|
| `describe` single structure via CLI (ref backend, cold) | 17-19 s wall; in-process 7.7 s for a dimer of which 4.3 s parse/warm-up, 3.3 s BURROW, 0.1 s fingerprinting |
| `describe` 100 dimers, `--backend fast`, 8 CPU workers | 55 s wall (37.5 s inside the run, 0.375 s/structure aggregate); per structure median 0.73 s (BURROW 0.54 s, fingerprints 0.09 s); 984 pockets |
| `describe` 100 dimers, `--device cuda:2`, 8 workers | 561 s wall: slower than CPU on this node (median 0.48 s/structure but long start-up/parse tails up to 43 s). Fingerprints identical to CPU. Not recommended at this batch size; README reports 380 structures/s with 4x RTX 4090 at AFDB scale |
| compare | a dot product; README: 0.1-0.7 µs/pair batched. Our small-matrix Python loop: 60 µs/pair (overhead) |
| all-vs-atlas search | human tier-1 atlas (90,560 pockets): load fingerprints 19 s, exact `IndexFlatIP` search 6.6 ms/query on 1 thread. Top hits of the 14-3-3 dimer pockets were the paralogous 14-3-3 monomer pockets (cosine 0.92-0.96) |

All-vs-all against the atlas is efficient with FAISS: tiers 1-2 (3.23 M pockets × 768 B fp16 ≈ 2.5 GB, 5 GB as
float32) fit in RAM for an exact flat index (≈ 0.25 s/query single-thread, extrapolated linearly from 6.6 ms per
90 k, or batched matrix products on GPU). For the full atlas use `--kind ivfpq` (README: 0.77 ms/query with
re-rank, recall@10 0.987). No index has been built for the atlas yet.

**Multi-chain.** Works on multi-chain files (all chains by default); interface pockets get fingerprints like any
other pocket, with residues from both chains in `residues`. The atlas contains **monomers only**, so dimer
interface pockets have no counterpart there. The same site can be merged/split differently in monomer and dimer
context (14-3-3 beta: dimer pockets m1/m2 vs monomer m0/m1 scored only 0.63-0.71), so match monomer and dimer
pockets by residue overlap, not by `pocket_id` or by fingerprint alone.

**Threads.** With `env_threads.sh`: main 5, workers 1 each.

## 4. adock (docking)

**What it does.** Vina-family docking engine (Vina/Vinardo scoring re-implemented, plus an 18-weight fitted
scoring function) with a CUDA search kernel (cupy NVRTC) and an OpenMP CPU kernel (g++ compiled at first use).
Needs a box; no pocket detection. README: on holo PoseBusters V2 the final preset (v2_fitted) reaches 58.6 %
top-1 ≤ 2 Å; on AlphaFold receptors success falls to ≈ 11 % (vs 50.5 % on the experimental receptor), and GNINA
falls similarly. Treat docking into predicted dimers as low-confidence pose generation.

**CLI.**
```bash
source scripts/env_tools.sh dock
python -m adock.dock --receptor rec.pdb --ligand lig.sdf --box box.json --out poses.sdf --mode v2_fast --device cpu --seed 0
python scripts/run_adock.py --set <records> ...   # batch runner over prepared record sets (tools/adock/README)
```
Presets (`config/adock.json`): `v2_fast`, `v2_default`, `v2_accurate`, `v2_fitted` (final holo method, **GPU only**).
`--af-mode` values other than `none` and the rerank/consensus presets need modules not included in the release.
Python API for batches in one process: `adock.dock.dock_batch(complexes, mode=..., device=..., timings={})`.

**Inputs.** Receptor: **PDB format** (own reader; mmCIF must be converted), heavy atoms typed by element/bonds,
hydrogens optional, waters dropped, atoms cropped to the box, so multi-chain receptors are fine. Ligand: SDF with
hydrogens (3D input conformer). Box: JSON `{"center": [x, y, z], "size": [sx, sy, sz]}` in Å; a burrow/
pocketcompare pocket gives the centre, and the extent of its `.points.csv` grid points + margin gives the size
(our choice; adock's own record boxes come from ligands). **Outputs:** SDF poses (best first; up to 9 after 2 Å clustering, only 1 for 1STP) and a JSON line
on stdout with `top_affinity` (Vina-like kcal/mol) and stage timings.

**Measured (1STP streptavidin/biotin fixture from the README, `scripts/test_adock.sh`, CPU engine, 8 OpenMP threads):**

| run | top-1 RMSD vs crystal | affinity | time |
|---|---|---|---|
| v2_default, CPU, first call in process | 0.58 Å | −7.51 | 183 s (prep 150 s: process-pool start-up and imports on the slow filesystem; grid 13 s; search 18 s) |
| v2_fast, CPU, second call (warm) | 0.73 Å | −7.48 | 2.1 s (prep 0.6 s, grid 0.4 s, search 1.0 s) |
| an earlier separate CLI call, v2_default CPU | same pose (−7.51) | | 635 s wall, of which ≈ 480 s imports/compilation |

**GPU engine does not work on this node.** On Quadro RTX 6000 (Turing, sm_75; GPUs 2 and 3, which pass
`scripts/check_gpus.py`) the CUDA Monte-Carlo kernel (`adock/search.py` → `run_mc`) ran at 100 % GPU utilisation for
15-30 min on the single 1STP complex without finishing (killed; stack dump in `dock_problem_kernel`). The CPU
engine finishes the same search in 1-18 s. The repository was developed and tested on a Blackwell RTX PRO 6000
(sm_120); sm_75 was never tested by its authors. The env check passes for torch, cupy and openmm on these GPUs.
Consequences: `v2_fitted` (GPU-only) is unavailable here; use `v2_default`/`v2_fast` on CPU, or fix/test the
kernel on sm_75, or run on another node. GPU docking throughput could therefore not be measured (README: ≈ 3 s
per complex warm on Blackwell, 46-54× faster than GNINA).

**Imports are very slow from the workspace (Lustre):** the env has 70,347 files; `import torch` took 192 s,
`import adock.dock` 312 s cold and 68 s warm on 2026-10-29. Run all dockings of a batch in one process. If this
stays a problem, a node-local copy of the env (e.g. `conda-pack` into `/tmp`, 1.7 TB local disk; what the adock
authors did) would help; this was not done because environments must live in the workspace.

**Threads.** The CPU kernel uses `OMP_NUM_THREADS` (set to 8 for the test; `env_threads.sh` would make it 1);
`prep_workers` = 8 processes in `config/adock.json`. Main process during GPU runs: 4 threads.

## 5. Druggability

None of the three tools computes a druggability score. burrow's `score` is a ligand-binding-site likelihood
(trained on BioLiP sites, **including ions** with the default `--ranker all`; `--ranker small` scores non-ion
small-molecule sites, P2Rank-style). pocketcompare has no druggability output. adock gives docking affinities
(Vina-like kcal/mol), not druggability. If needed, options are an external tool (e.g. fpocket's druggability
score, PockDrug, DoGSiteScorer) or a proxy from burrow's `--ranker small` score, volume and buriedness.

## 6. Gotchas

1. **GPU 0 of this node is faulty**: fp64 matmul off by ~2 % and corrupted int64 transfers
   (`scripts/check_gpus.py`); burrow's GPU kernels hit device-side asserts on it. GPUs 1-3 pass. Use
   `CUDA_VISIBLE_DEVICES=1,2,3` (set by `scripts/env_tools.sh`). burrow's GPU path swallows warm-up errors and
   falls back to the CPU silently; check `timing.csv` column `device` (or `BURROW_FAST_PROFILE`, `n_cpu_path`).
2. **burrow on GPU does not help small batches**: the hybrid scheduler gives work to CPU workers while the GPU
   initialises (≈ 50-80 s here). Use `BURROW_FAST_HYBRID=0` and large batches if you want the GPU.
3. **Start-up cost dominates single calls** (≈ 15-18 s for burrow/pocketcompare; 1-5 min for adock imports from the
   workspace filesystem). Always batch: give lists (`@list.txt`) or directories, or call the Python APIs in one
   long-lived process.
4. **Caches default to $HOME or /tmp**: Numba (`NUMBA_CACHE_DIR`), Triton (`TRITON_CACHE_DIR`, default
   `~/.triton`), CuPy (`CUPY_CACHE_DIR`, default `~/.cupy`), adock (`ADOCK_WORK`, `ADOCK_CPU_CACHE`, default
   `/tmp/adock`). `scripts/env_tools.sh` points all of them into `.cache/`. Note: `~/.triton` was already
   written by our first burrow GPU runs before this was set.
5. **`--af-mode` only for predicted models** (B-factor = pLDDT). Never on crystal structures.
6. The heterodimer PDB files' `TITLE` names only one protein; use the `DBREF` lines (chain → UniProt) for both
   partners. Model ids (`AF-0000000203452075`) are not UniProt accessions.
7. `burrow predict` output (v2 list with score floor) and `pocketcompare describe` pockets (merged region pockets,
   no floor) are different pocket sets with different numbering. For atlas-compatible pockets use
   `pocketcompare describe` and the atlas filters.
8. burrow's `timing.csv` header has fewer fields than its rows (`backend_or_timings...`); pass explicit column
   names when parsing.
9. The atlas covers AFDB tiers 1-2 only (no tier 3 yet) and lives on scratch (purged after 15 days): copy what we need.
10. Compute policy: `source scripts/env_threads.sh` before any run (else 64 threads per main process), ≤ 8
    workers for tests, register long jobs in `logs/compute_ledger.txt`.
