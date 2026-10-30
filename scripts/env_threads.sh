# Source before ANY parallel job: caps per-process thread pools.
# Reason: account limit ulimit -u = 5120 threads+processes; a prior run crashed when
# several agents launched pools concurrently (each numpy/torch worker spawns ~64 BLAS/OMP threads).
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-1}
export TF_NUM_INTRAOP_THREADS=${TF_NUM_INTRAOP_THREADS:-1}
export TF_NUM_INTEROP_THREADS=${TF_NUM_INTEROP_THREADS:-1}
export POLARS_MAX_THREADS=${POLARS_MAX_THREADS:-1}
