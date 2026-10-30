"""Verifier extras: LIGYSIS full-atlas recall (A denominator, acc level) with protein bootstrap; fpocket drug-like stratum delta;
hetero-ext pre-cutoff complex count."""
from pathlib import Path
import pandas as pd, numpy as np
ROOT = str(Path(__file__).resolve().parents[2])
a = pd.read_csv(f'{ROOT}/results/benchmarks/ligysis/ligysis_site_hits.tsv', sep='\t')
a = a[a.level == 'acc']
hit = a.hit_dcc12_topN2.astype(str).eq('True').values  # NaN (no pocket / no centre) counted as miss (denominator A)
print('A: n sites', len(a), 'proteins', a.ACC.nunique(), 'recall %.4f' % hit.mean())
cl, idx = np.unique(a.ACC.values, return_inverse=True)
S = np.bincount(idx, hit); C = np.bincount(idx)
rng = np.random.default_rng(0); bs = []
for _ in range(1000):
    i = rng.integers(0, len(cl), len(cl)); bs.append(S[i].sum() / C[i].sum())
print('protein bootstrap CI', np.round(100 * np.percentile(bs, [2.5, 97.5]), 1), '(producer used protein-cluster bootstrap)')
c = a[a.n_pockets_atlas.fillna(0).gt(0) & a.has_centre.astype(str).eq('True')]
print('C: n', len(c), 'recall %.4f' % c.hit_dcc12_topN2.astype(str).eq('True').mean())
f = pd.read_parquet(f'{ROOT}/data/interim/fpocket/pdbsite_sites.parquet'); f = f[(f.status == 'ok') & (f.stratum_lig == 'drug_like')]
x = f[f.y_mbs == 1].fp_drug.values; y = f[f.y_mbs == 0].fp_drug.values
print('fpocket drug-like unweighted delta %.3f' % ((np.sign(x[:, None] - y[None]).sum()) / (len(x) * len(y))), 'frac>0.5', (x > 0.5).mean(), (y > 0.5).mean())
p = pd.read_csv(f'{ROOT}/results/cofold/hetero_ext/paired.tsv', sep='\t')
print('hetero paired rows', len(p), 'pair_complex_pre_cutoff True:', int(p.pair_complex_pre_cutoff.astype(str).eq('True').sum()))
