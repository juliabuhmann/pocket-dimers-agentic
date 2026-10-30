"""EXPLORATORY / UNVERIFIED. Paired AFDB-atlas vs experimental-chain BURROW recall on the same sampled LIGYSIS sites (segment-level N)."""
import pandas as pd
O = 'results/benchmarks/ligysis/'
a = pd.read_csv(O + 'ligysis_site_hits.tsv', sep='\t'); a = a[(a.level == 'seg') & (a.has_centre == True) & (a.n_pockets_atlas > 0)]
e = pd.read_csv(O + 'ligysis_exp_site_hits.tsv', sep='\t').dropna(subset=['hit_dcc12_topN2'])
m = e.merge(a, on='SITE_NAME', suffixes=('_exp', '_afdb'))
rows = []
for c in ['hit_dcc12_topN2', 'hit_dcc4_topN2', 'hit_ovl_topN2']:
    rows.append(dict(criterion=c, n_sites=len(m), n_chains=m.rep_chain.nunique(), recall_exp=m[c + '_exp'].astype(bool).mean(),
                     recall_afdb=m[c + '_afdb'].astype(bool).mean(), exp_only=int((m[c + '_exp'].astype(bool) & ~m[c + '_afdb'].astype(bool)).sum()),
                     afdb_only=int((~m[c + '_exp'].astype(bool) & m[c + '_afdb'].astype(bool)).sum())))
r = pd.DataFrame(rows); r.to_csv(O + 'ligysis_paired_exp_vs_afdb.tsv', sep='\t', index=False); print(r.round(3).to_string())
print('exp sample sites with computable centre+pockets:', len(e), 'of 480')
