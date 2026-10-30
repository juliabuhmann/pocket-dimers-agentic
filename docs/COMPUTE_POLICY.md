# Compute policy (all agents must follow)

Hard limit: `ulimit -u` = 5120 threads+processes for our user (shared by every agent, shell and job).
A previous session crashed when four agents started pools at once (22 + 32 workers + AF3 queue).

Rules
1. `source scripts/env_threads.sh` before any multi-process job (sets OMP/BLAS/etc. threads = 1 per worker).
2. CPU worker budget across ALL concurrent jobs: ≤ 48 workers total. Before launching, check
   `ps -u $(id -u) -L --no-headers | wc -l` (current threads) and `cat logs/compute_ledger.txt`.
   Do not launch if current threads > 2500.
3. Register every long job (>5 min or >4 workers) in `logs/compute_ledger.txt`:
   `<start ISO time> | <agent/task> | <n workers> | <GPU ids> | <pid> | <log path>` and mark it done when finished.
4. Default per job: ≤ 16 workers unless the PI allots more. GPU jobs: one process per GPU, pin with CUDA_VISIBLE_DEVICES.
   JAX/AF3: also set XLA_FLAGS threads modestly; torch: torch.set_num_threads(small).
5. Use `concurrent.futures`/`multiprocessing` with explicit max_workers; never `Pool()` without a size; never nested pools.
6. Long jobs: run with nohup, log to logs/, make them resumable (skip outputs that exist).
