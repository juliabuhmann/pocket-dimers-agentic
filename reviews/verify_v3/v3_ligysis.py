"""Verifier: paired LIGYSIS recall, experimental chains vs AFDB atlas (seg-level N), from per-site hit tables."""
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import binomtest
ROOT = str(Path(__file__).resolve().parents[2] / 'results/benchmarks/ligysis')
DCC_THR = 12.0  # LIGYSIS recommended DCC threshold (Utges & Barton 2024)
e = pd.read_csv(f'{ROOT}/ligysis_exp_site_hits.tsv', sep='\t')
a = pd.read_csv(f'{ROOT}/ligysis_site_hits.tsv', sep='\t')
print('exp rows', len(e), 'chains', e.rep_chain.nunique(), 'afdb levels', a.level.value_counts().to_dict())
for lvl in ['seg', 'acc']:
    m = e.merge(a[a.level == lvl], on='SITE_NAME', suffixes=('_e', '_a'))
    for c in ['hit_dcc12_topN2', 'hit_dcc4_topN2', 'hit_ovl_topN2']:
        mm = m[m[c + '_e'].notna() & m[c + '_a'].notna()]
        x = mm[c + '_e'].astype(str).eq('True'); y = mm[c + '_a'].astype(str).eq('True')
        b, cc = int((x & ~y).sum()), int((~x & y).sum())
        print(lvl, c, 'n', len(mm), 'chains', mm.rep_chain.nunique(), 'exp %.4f afdb %.4f' % (x.mean(), y.mean()),
              'disc', b, cc, 'exact McNemar p %.2g' % binomtest(b, b + cc).pvalue)
# internal consistency: dcc12 hit flag == min_dcc_topN2 <= 12
for name, t in [('exp', e), ('afdb', a)]:
    t = t[t.min_dcc_topN2.notna()]
    print(name, 'hit_dcc12_topN2 vs min_dcc<=12 disagreements:', int((t.hit_dcc12_topN2.astype(str).eq('True') != (t.min_dcc_topN2 <= DCC_THR)).sum()))
