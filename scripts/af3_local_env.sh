# Source to run AF3 from NODE-LOCAL copies (/tmp) instead of Lustre (LESSONS rule 7: Lustre import stalls).
# Staging (rerunnable): scripts/stage_af3_local.sh. Status: first attempt, test pending.
L=/tmp/${USER}_af3
export PATH="$L/af3/bin:$PATH" VIRTUAL_ENV="$L/af3"
export AF3_SRC="$L/alphafold3_src"
export AF3_MODEL_DIR="$L/data/af3/models" AF3_JAX_CACHE="$L/jax_cache"
