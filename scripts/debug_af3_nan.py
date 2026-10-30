#!/usr/bin/env python3
"""Debug: run one AF3 inference and report which raw outputs contain NaN/inf. usage: debug_af3_nan.py INPUT.json"""
import os, sys, pathlib
sys.path.insert(0, os.environ["AF3_SRC"])
import numpy as np, jax
import run_alphafold as ra
import tokamax, jax.numpy as jnp
if os.environ.get("AF3_GLU_PURE"):
    def _glu(x, weights, activation=None, precision=None, **kw):
        h = jnp.einsum('...i,iof->...of', x, weights)
        return activation(h[..., 0, :]) * h[..., 1, :]
    tokamax.gated_linear_unit = _glu
from alphafold3.data import featurisation
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components

fi = list(folding_input.load_fold_inputs_from_path(pathlib.Path(sys.argv[1])))[0]
impl = os.environ.get("AF3_ATTN", "xla")
cfg = ra.make_model_config(flash_attention_implementation=impl,
                           num_diffusion_samples=int(os.environ.get("AF3_NSAMP", 1)),
                           num_recycles=int(os.environ.get("AF3_NREC", 1)))
if os.environ.get("AF3_BFLOAT16"):
    cfg.global_config.bfloat16 = os.environ["AF3_BFLOAT16"]
runner = ra.ModelRunner(config=cfg, device=jax.local_devices()[0], model_dir=pathlib.Path(os.environ["AF3_MODEL_DIR"]))
ex = featurisation.featurise_input(fold_input=fi, buckets=None, ccd=chemical_components.Ccd(user_ccd=fi.user_ccd), verbose=False)[0]
res = runner.run_inference(ex, jax.random.PRNGKey(1))
def walk(d, p=""):
    for k, v in d.items():
        if isinstance(v, dict): walk(v, p + k + "/")
        else:
            a = np.asarray(v)
            if a.dtype.kind == "f":
                print(f"{p}{k}: shape={a.shape} nan={np.isnan(a).mean():.3f} inf={np.isinf(a).mean():.3f}")
walk(res)
