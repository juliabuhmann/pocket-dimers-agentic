"""Verifier: co-folding seed success counts and hetero-extension paired tests, from metrics tables."""
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import binomtest, spearmanr
ROOT = str(Path(__file__).resolve().parents[2] / 'results/cofold')
THR = 2.0  # success = ligand RMSD <= 2 A (docs/COFOLD_HETERO_EXT.md, docs/COFOLDING.md)
for exp in ['onaxis', 'dimeronly']:
    m = pd.read_csv(f'{ROOT}/{exp}/seeds_metrics.tsv', sep='\t')
    m['succ'] = m.lig_rmsd_sym <= THR
    print(exp, 'flag vs recomputed disagreements:', int((m.succ != m.success_sym.astype(bool)).sum()) if 'success_sym' in m else 'na',
          'NaN rmsd rows:', int(m.lig_rmsd_sym.isna().sum()))
    print(m.groupby(['group', 'cond', 'seed']).succ.agg(['sum', 'count']).to_string())
    # per-case agreement across available seeds
    agg = m.groupby(['group', 'cond', 'case_id']).succ.agg(['nunique', 'count'])
    print('cases with consistent outcome across seeds:', agg.groupby(['group', 'cond']).apply(lambda t: f"{int((t['nunique']==1).sum())}/{len(t)}").to_dict())
    s0 = m[m.seed == 0].pivot_table(index=['group', 'case_id'], columns='cond', values=['lig_rmsd_sym', 'succ'], aggfunc='first')
    for g, t in s0.groupby(level=0):
        d, mo = t[('succ', 'dimer')].astype(bool), t[('succ', 'monomer')].astype(bool)
        b, c = int((d & ~mo).sum()), int((~d & mo).sum())
        lower = int((t[('lig_rmsd_sym', 'dimer')] < t[('lig_rmsd_sym', 'monomer')]).sum()); n = len(t)
        print(exp, g, 'seed0 n', n, 'dimer', int(d.sum()), 'single', int(mo.sum()), 'discordant', b, c,
              'McNemar exact p %.4g' % (binomtest(b, b + c).pvalue if b + c else 1), 'dimer lower RMSD %d/%d sign p %.4g' % (lower, n, binomtest(lower, n).pvalue))
# hetero extension
h = pd.read_csv(f'{ROOT}/hetero_ext/metrics.tsv', sep='\t')
print('hetero rows', len(h), h.status.value_counts().to_dict())
h = h[h.status == 'ok']
h['succ'] = h.rmsd_ver <= THR
print('flag vs recomputed disagreements:', int((h.succ != h.success_ver.astype(bool)).sum()), 'prod vs ver success disagreements:', int((h.success_ver != h.success_prod).sum()))
p = h.pivot_table(index=['case_id', 'origin'], columns='cond', values=['rmsd_ver', 'succ', 'partner_share'], aggfunc='first').dropna()
d, mo = p[('succ', 'dimer')].astype(bool), p[('succ', 'monomer')].astype(bool)
b, c = int((d & ~mo).sum()), int((~d & mo).sum()); n = len(p)
lower = int((p[('rmsd_ver', 'dimer')] < p[('rmsd_ver', 'monomer')]).sum())
print('hetero paired n', n, 'by origin', p.reset_index().origin.value_counts().to_dict(), 'dimer', int(d.sum()), 'single', int(mo.sum()),
      'disc', b, c, 'McNemar exact p %.4f' % binomtest(b, b + c).pvalue, 'dimer lower %d/%d sign p %.5f' % (lower, n, binomtest(lower, n).pvalue))
gain = np.log10(p[('rmsd_ver', 'monomer')] / p[('rmsd_ver', 'dimer')])
print('spearman partner_share vs log gain: rho %.3f p %.3f' % spearmanr(p[('partner_share', 'dimer')], gain))
