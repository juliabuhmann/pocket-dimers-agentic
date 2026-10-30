# Random seeds and RNG usage

Compiled by grep over `scripts/`, `src/`, `tests/` (patterns: np.random, default_rng, random_state, set.seed, modelSeeds, manual_seed, PRNGKey, --seed). Status: first-pass inventory, line-level semantics not independently verified. Project convention (DESIGN.md section 5): seeds 0, 1, 2; seed 0 primary.

## Statistical resampling (numpy `default_rng`)

| Script | Seed variable | Value | Purpose |
|---|---|---|---|
| scripts/pdbsite_common.py | `SEED`, `SEEDS`, `N_BOOT` | 0; (0,1,2); 1000 | shared constants for pdbsite_* scripts |
| scripts/pdbsite_sample.py | `--seed` (default `SEED`) | 0 | stratified sample of PDB sites (`rng.random`, `.sample(random_state=)`, shuffle) |
| scripts/pdbsite_compare.py | `seed` loop over `SEEDS` (all_weighted stratum), `SEEDS[:1]` otherwise | 0,1,2 / 0 | cluster bootstrap |
| scripts/pdbsite_crop_sensitivity.py, pdbsite_figures.py | reads `seed` column / `seed == 0` | 0 | filter only, no RNG |
| scripts/q22_common.py | `SEEDS`, `N_BOOT`, `cluster_boot(seed=0)` | (0,1,2); 1000; 0 | cluster bootstrap |
| scripts/q22_compare.py | `seed` arg (default 0); line 189 `default_rng(0)` | 0 (and SEEDS) | bootstrap of rates/ratios |
| scripts/q22_novelty.py | `seed` arg of `size_matched` | SEEDS | volume-matched bootstrap |
| scripts/q22_dimer_only_control.py | `excess(seed=0)` | 0 | bootstrap |
| scripts/q22_figures.py | filter `seed == 0` | 0 | no RNG |
| scripts/q24_conservation.py | `N_BOOT, SEED` | 1000, 0 | random pick of ok sites, bootstrap |
| scripts/q24_stats.py | imports `SEED` | 0 | bootstrap weights |
| scripts/h4_common.py | `SEEDS`, `N_BOOT`, `N_PERM` | [0,1,2]; 1000; 10000 | H4 constants |
| scripts/h4_stats.py | `seed` args (default 0), loops over SEEDS | 0,1,2 | bootstrap, within-cluster label permutation, sign flips |
| scripts/h4_recovery.py | `default_rng([sd,row_id,ai,cond_idx])` for sd in SEEDS; `default_rng([0,row_id,ai])` | 0,1,2; 0 | random size-matched null patches; random-rotation negative control |
| scripts/h5_analysis.py | `SEEDS`, `N_BOOT`; `boot_rates(seed=)` | [0,1,2]; 1000 | protein-level multinomial bootstrap; shuffle-seed negative controls passed to R |
| scripts/h5_glmm.R | arg 5 `shuffle_seed` (default -1 = no shuffle) | 0,1,2 from h5_analysis | `set.seed(seed)`, permute categories within protein (negative control) |
| scripts/part1_analysis.py | `SEEDS`, `N_BOOT` | [0,1,2]; 1000 | cluster bootstrap and differences |
| scripts/fpocket_compare.py | `SEED`, `N_BOOT` | 0; 1000 | bootstrap |
| scripts/fpocket_afdimer.py | `SEED` | 0 | random matched-candidate pick (`rng.permutation`) |
| scripts/fpocket_subsample.py | `SEED`, `N_PER_STRATUM` | 0; 500 | per-stratum subsample |
| scripts/fpocket_pdbsite.py | `sample(frac=1, random_state=0)` | 0 | processing order (inline literal) |
| scripts/fpocket_figures.py | `default_rng(0)` | 0 | plot jitter only (inline literal) |
| scripts/iface_lm_analyze.py | `SEED, K_FOLDS` | 0; 5 | CV fold assignment, bootstrap |
| scripts/ligysis_03_summarise.py | `SEED`, `B` | 0; 2000 | bootstrap |
| scripts/ligysis_04_burrow_exp.py | `SEED`, `SAMPLE_N` | 0; 200 | representative sample |
| scripts/run_part1_census.py / src/part1/config.py | `SAMPLE_SEED` | 0 | 1% id sample (`rng.choice`) |
| scripts/annotate_pocket_chains.py | `.sample(n, random_state=0)` | 0 | inline literal, subsample |
| scripts/summarize_pocket_run.py | `.sample(5, random_state=1)` | 1 | inline literal, display sample |

## Structure prediction / docking

| Script | Seed variable | Value | Purpose |
|---|---|---|---|
| scripts/cofold_select_onaxis.py, cofold_select_h4.py, cofold_select_dimeronly.py | `modelSeeds` | [0] | AF3 seed in generated input JSON (inline literal) |
| scripts/cofold_seeds_jobs.py, cofold_run_seeds.sh | `PROJECT_SEEDS`; `--seeds` default | [0,1,2]; 1 2 | replicate AF3 seeds |
| scripts/cofold_seed1_onaxis.py | `SEED` | 1 | extra AF3 seed, on-axis dimer runs |
| scripts/cofold_metrics_seeds.py | `SEEDS` | [0,1,2] | aggregation only |
| scripts/cofold_metrics_hetero_ext.py | `SEED` | 0 | check use in script (defined, usage not inspected) |
| scripts/make_af3_input.py | `--seeds` | default [1] | `modelSeeds` for ad hoc inputs (differs from project default 0) |
| scripts/test_adock.sh | `--seed` | 0 | adock docking test |
| scripts/debug_af3_nan.py | `jax.random.PRNGKey` | 1 | debug script |
| scripts/test_gpu_basic.py, test_gpu_numerics.py, check_gpus.py | `default_rng(0)`, `torch.manual_seed(0)` | 0 | GPU diagnostics |

## Flags

- No unseeded `np.random.*` / `random.*` / `set.seed` call was found in the scanned files. All `default_rng` calls carry an explicit seed.
- Seeds hard-coded as literals rather than named constants: fpocket_figures.py, fpocket_pdbsite.py, annotate_pocket_chains.py, summarize_pocket_run.py, h4_recovery.py (rotation control), test_gpu_*.py, debug_af3_nan.py.
- `make_af3_input.py` defaults to `modelSeeds=[1]` while project runs use 0: check which was used for any given input JSON.
- `cofold_metrics_hetero_ext.py` defines `SEED = 0`; not checked whether it is used.
- Residual non-determinism not controlled by these seeds: AF3 (JAX, GPU kernels; Turing cc 7.5 per docs/COFOLDING.md), fpocket, BURROW/torch inference, adock GPU docking (only `--seed 0` set in the test script), multi-threaded fetches. Pandas `.sample`/`shuffle` outputs depend on pandas version (see envs/).
- Tools in `tools/` (burrow, pocketcompare, adock, alphafold3) were not scanned.
