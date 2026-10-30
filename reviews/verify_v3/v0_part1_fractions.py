"""Verifier recomputation of Part 1 MBS fractions (census v2) from data/part1/ligand_sites.parquet.
Definitions from DESIGN.md §0, docs/PART1_PIPELINE.md (primary set) and docs/PART1_RESULTS.md §2 (units). No producer code."""
from pathlib import Path
import pandas as pd, numpy as np, pyarrow.parquet as pq, sys
ROOT = str(Path(__file__).resolve().parents[2])
cols = ['pdb_id', 'ccd_key', 'assembly_type', 'marsh_class', 'is_mbs', 'in_primary_set', 'is_photosynthetic',
        'contacts_antibody', 'cluster30', 'site_group', 'resolution', 'is_dup_copy', 'is_bound', 'is_glycan',
        'is_additive', 'biolip_covered', 'ligand_class']
d = pq.read_table(f'{ROOT}/data/part1/ligand_sites.parquet', columns=cols).to_pandas()
out = []
def rec(k, v, n=None):
    out.append((k, round(100 * v, 2), n)); print(k, round(100 * v, 2), n)
# independent check of in_primary_set definition
mine = (~d.is_dup_copy) & d.is_bound & (~d.is_glycan) & (~d.is_additive.fillna(True)) & d.biolip_covered
print('primary-set flag vs own rule disagreements:', int((mine != d.in_primary_set).sum()), 'n_primary', int(d.in_primary_set.sum()))
p = d[d.in_primary_set]
print(p.assembly_type.value_counts().to_dict())
for at in ['homomer', 'heteromer']:
    s = p[p.assembly_type == at]
    rec(f'instance {at} all', s.is_mbs.mean(), len(s))
    for mc in ['cofactor', 'metal', 'small_molecule']:
        t = s[s.marsh_class == mc]; rec(f'instance {at} {mc}', t.is_mbs.mean(), len(t))
    t = s[~s.is_photosynthetic.fillna(False)]; rec(f'instance {at} excl photosynthetic', t.is_mbs.mean(), len(t))
# sites / units
o = p[p.assembly_type.isin(['homomer', 'heteromer']) & ~p.contacts_antibody.fillna(False) & p.cluster30.notna()].copy()
print('oligomer rows', len(o), 'clusters', o.cluster30.nunique())
sites = o.groupby(['pdb_id', 'site_group']).agg(mbs=('is_mbs', 'max'), cluster30=('cluster30', 'first'),
        ccd_key=('ccd_key', 'first'), at=('assembly_type', 'first'), res=('resolution', 'first')).reset_index()
print('sites', len(sites), 'site groups spanning >1 ccd_key:', int(o.groupby(['pdb_id', 'site_group']).ccd_key.nunique().gt(1).sum()))
rec('sites all (raw)', sites.mbs.mean(), len(sites))
# units: (cluster30, ccd_key); representative entry = best resolution, ties earliest pdb id
e = o.groupby(['cluster30', 'ccd_key', 'pdb_id']).agg(res=('resolution', 'min'), mbs=('is_mbs', 'max'),
                                                       at=('assembly_type', 'first')).reset_index()
e['res_s'] = e.res.fillna(np.inf)
e = e.sort_values(['cluster30', 'ccd_key', 'res_s', 'pdb_id'])
u = e.groupby(['cluster30', 'ccd_key']).head(1)
print('units', len(u), u['at'].value_counts().to_dict())
rec('units all', u.mbs.mean(), len(u))
for at in ['homomer', 'heteromer']:
    rec(f'units {at}', u[u['at'] == at].mbs.mean(), int((u['at'] == at).sum()))
# compare to producer unit table
up = pd.read_parquet(f'{ROOT}/data/part1/analysis/units_primary_v2.parquet')
m = u.merge(up[['cluster30', 'ccd_key', 'pdb_id', 'y_mbs']], on=['cluster30', 'ccd_key'], how='outer', indicator=True)
print('merge with producer units:', m._merge.value_counts().to_dict(), 'rep pdb mismatches', int((m.pdb_id_x != m.pdb_id_y).sum()),
      'y_mbs mismatches', int((m.mbs.astype(float) != m.y_mbs.astype(float)).sum()))
# cluster bootstrap CIs for sites/units (1000, seed 0)
rng = np.random.default_rng(0)
for name, t, col in [('sites', sites, 'mbs'), ('units', u, 'mbs')]:
    g = t.groupby('cluster30')[col].agg(['sum', 'count'])
    S, C = g['sum'].values, g['count'].values
    bs = []
    for _ in range(1000):
        i = rng.integers(0, len(g), len(g)); bs.append(S[i].sum() / C[i].sum())
    print(name, 'CI', np.round(100 * np.percentile(bs, [2.5, 97.5]), 1))
pd.DataFrame(out, columns=['quantity', 'pct', 'n']).to_csv(f'{ROOT}/reviews/verify_v3/v0_part1_fractions.tsv', sep='\t', index=False)
lo = m[m._merge == 'left_only']
print(lo[['cluster30', 'ccd_key', 'pdb_id_x', 'mbs', 'at']].head(15).to_string())
print('left_only at:', lo['at'].value_counts().to_dict(), 'mbs rate', lo.mbs.mean())
x = o[o.set_index(['cluster30', 'ccd_key']).index.isin(lo.set_index(['cluster30', 'ccd_key']).index)]
print(x[['pdb_id', 'ccd_key', 'assembly_type', 'marsh_class', 'ligand_class', 'resolution']].drop_duplicates().head(15).to_string())
print('ligand_class of left_only:', x.ligand_class.value_counts().to_dict())
u2 = u[~u.set_index(['cluster30', 'ccd_key']).index.isin(lo.set_index(['cluster30', 'ccd_key']).index)]
print('after dropping left_only: all', round(100*u2.mbs.mean(), 2), {a: round(100*u2[u2['at']==a].mbs.mean(), 2) for a in ['homomer', 'heteromer']})
# Variant: unit cluster taken per collapsed site (one cluster30 per site), not per instance
ncl = o.groupby(['pdb_id', 'site_group']).cluster30.nunique()
print('sites whose instances span >1 cluster30:', int((ncl > 1).sum()))
for how in ['first', 'min']:
    sc = o.groupby(['pdb_id', 'site_group']).cluster30.agg(how).rename('c').reset_index()
    s2 = sites.drop(columns='cluster30').merge(sc, on=['pdb_id', 'site_group'])
    e2 = s2.groupby(['c', 'ccd_key', 'pdb_id']).agg(res=('res', 'min'), mbs=('mbs', 'max'), at=('at', 'first')).reset_index()
    e2['rs'] = e2.res.fillna(np.inf)
    u3 = e2.sort_values(['c', 'ccd_key', 'rs', 'pdb_id']).groupby(['c', 'ccd_key']).head(1)
    print('site-level cluster', how, 'units', len(u3), 'all', round(100*u3.mbs.mean(), 2),
          {a: (int((u3['at']==a).sum()), round(100*u3[u3['at']==a].mbs.mean(), 2)) for a in ['homomer', 'heteromer']})
