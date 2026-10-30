"""Verifier: fpocket druggability MBS vs SBS from per-site table; weights units_available / n_ok per stratum."""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = str(Path(__file__).resolve().parents[2])
THR = 0.5  # fpocket conventional druggability cutoff (Schmidtke & Barril 2010)
d = pd.read_parquet(f'{ROOT}/data/interim/fpocket/pdbsite_sites.parquet')
d = d[d.status == 'ok']
st = pd.read_csv(f'{ROOT}/results/part1/pocket_descriptors/sample_strata_n12000_seed0.tsv', sep='\t', comment='#')
avail = dict(zip(st.stratum, st.units_available))
nok = d.stratum.value_counts()
d['wv'] = d.stratum.map(lambda s: avail[s] / nok[s])
print('n ok', len(d), 'MBS', int(d.y_mbs.sum()), 'clusters', d.cluster30.nunique(), 'weight check vs producer w (normalised):',
      float(np.max(np.abs(d.wv / d.wv.mean() - d.w / d.w.mean()))))
def wdelta(x, wx, y, wy):
    o = np.argsort(y); ys = y[o]; cw = np.concatenate([[0], np.cumsum(wy[o])])
    lo = np.searchsorted(ys, x, 'left'); hi = np.searchsorted(ys, x, 'right')
    return float(np.sum(wx * (cw[lo] - (cw[-1] - cw[hi]))) / (wx.sum() * cw[-1]))
g1 = d.y_mbs.values == 1; v = d.fp_drug.values; w = d.wv.values; gt = (v > THR).astype(float)
def stats(mw):
    return (wdelta(v[g1], mw[g1], v[~g1], mw[~g1]), np.sum(mw[g1] * gt[g1]) / mw[g1].sum(), np.sum(mw[~g1] * gt[~g1]) / mw[~g1].sum())
est = stats(w)
cl, idx = np.unique(d.cluster30.values, return_inverse=True)
rng = np.random.default_rng(0)
bs = np.array([stats(np.bincount(rng.integers(0, len(cl), len(cl)), minlength=len(cl))[idx] * w) for _ in range(1000)])
print('weighted delta %.4f CI [%.3f, %.3f]' % (est[0], *np.percentile(bs[:, 0], [2.5, 97.5])))
print('frac>0.5 MBS %.4f SBS %.4f diff %.4f CI [%.3f, %.3f]' % (est[1], est[2], est[1] - est[2], *np.percentile(bs[:, 1] - bs[:, 2], [2.5, 97.5])))
print('unweighted delta %.4f' % wdelta(v[g1], np.ones(g1.sum()), v[~g1], np.ones((~g1).sum())))
print('medians', np.median(v[g1]), np.median(v[~g1]))
