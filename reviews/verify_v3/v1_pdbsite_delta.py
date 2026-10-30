"""Independent recomputation of weighted Cliff's delta (MBS vs SBS) at PDB ligand sites, uncropped.
Weights = units_available / sampled per stratum (from sample_strata tsv); cluster30 bootstrap (1000, seed 0).
Written by verifier; does not import producer code."""
from pathlib import Path
import numpy as np, pandas as pd, sys
ROOT = str(Path(__file__).resolve().parents[2])
SITES = f'{ROOT}/data/interim/pdbsite/sites.parquet'
STRATA = f'{ROOT}/results/part1/pocket_descriptors/sample_strata_n12000_seed0.tsv'
N_BOOT = 1000
SEED = 0
METRICS = ['burrow_score', 'volume_A3', 'merged_volume_A3', 'merged_n_lining', 'merged_frac_charged',
           'merged_frac_hydrophobic', 'frac_charged', 'frac_hydrophobic', 'merged_frac_aromatic', 'merged_buriedness',
           'merged_burrow_score']

def wdelta(x, wx, y, wy):
    o = np.argsort(y); ys = y[o]; cw = np.concatenate([[0], np.cumsum(wy[o])])
    lo = np.searchsorted(ys, x, 'left'); hi = np.searchsorted(ys, x, 'right')
    less = cw[lo]; greater = cw[-1] - cw[hi]
    return float(np.sum(wx * (less - greater)) / (wx.sum() * cw[-1]))

df = pd.read_parquet(SITES)
st = pd.read_csv(STRATA, sep='\t', comment='#')
wmap = dict(zip(st.stratum, st.units_available / st.sampled))
print('stratum values:', df.stratum.unique()[:8], file=sys.stderr)
df['wv'] = df.stratum.map(wmap)
assert df.wv.notna().all(), df.stratum[df.wv.isna()].unique()
print('my weights vs producer w (max abs diff, after normalising by mean):',
      np.nanmax(np.abs(df.wv / df.wv.mean() - df.w / df.w.mean())), file=sys.stderr)
df = df[df.status.eq('ok')] if 'status' in df and (df.status == 'ok').any() else df
print('n sites', len(df), 'MBS', int(df.y_mbs.sum()), 'clusters', df.cluster30.nunique(), file=sys.stderr)
cl_codes, cl_idx = np.unique(df.cluster30.values, return_inverse=True)
rng = np.random.default_rng(SEED)
mults = np.stack([np.bincount(rng.integers(0, len(cl_codes), len(cl_codes)), minlength=len(cl_codes))[cl_idx]
                  for _ in range(N_BOOT)])
rows = []
for m in METRICS:
    d = df[df[m].notna()]
    keep = df[m].notna().values
    g1 = d.y_mbs.values == 1
    v = d[m].values.astype(float); w = d.wv.values
    est = wdelta(v[g1], w[g1], v[~g1], w[~g1])
    est_u = wdelta(v[g1], np.ones(g1.sum()), v[~g1], np.ones((~g1).sum()))
    bs = []
    for b in range(N_BOOT):
        mb = mults[b][keep] * w
        bs.append(wdelta(v[g1], mb[g1], v[~g1], mb[~g1]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append(dict(metric=m, n_mbs=int(g1.sum()), n_sbs=int((~g1).sum()), delta_w=round(est, 4),
                     lo=round(lo, 3), hi=round(hi, 3), delta_unw=round(est_u, 4)))
    print(rows[-1], file=sys.stderr)
pd.DataFrame(rows).to_csv(f'{ROOT}/reviews/verify_v3/v1_pdbsite_delta.tsv', sep='\t', index=False)
