# Audit task 1: independent recomputation of MBS fractions from ligand_sites.parquet
import sys, pandas as pd, numpy as np
cols=['pdb_id','ligand_uid','assembly_type','ccd','is_dup_copy','is_bound','is_glycan','is_additive','biolip_status',
 'in_primary_set','is_mbs','is_mbs_min2','is_mbs_biolip2','is_mbs_4A','marsh_class','ligand_class','chain_counts','n_chains_contacted']
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=cols)
d.to_pickle(sys.argv[1]+'/ls.pkl')
print(len(d))
# own primary-set definition
mine=(~d.is_dup_copy)&d.is_bound&(~d.is_glycan)&(~d.is_additive.fillna(False).astype(bool))
print('primary own',mine.sum(),'flag',d.in_primary_set.sum(),'xor',(mine!=d.in_primary_set).sum())
# own MBS from chain_counts (chains with >=1 res); note: chain_counts may include peptides? check
def nch(cc,k=1): return sum(1 for x in cc if x['n_res']>=k) if cc is not None else 0
d['my_nch']=d.chain_counts.map(nch); d['my_nch2']=d.chain_counts.map(lambda c:nch(c,2))
d['my_mbs']=d.my_nch>=2
print('mbs mismatch vs is_mbs', (d.my_mbs!=d.is_mbs).sum(), 'nch mismatch', (d.my_nch!=d.n_chains_contacted).sum(), 'min2 mismatch',((d.my_nch2>=2)!=d.is_mbs_min2).sum())
p=d[d.in_primary_set]
t=p.groupby(['assembly_type','marsh_class']).agg(n=('my_mbs','size'),frac=('my_mbs','mean'),min2=('is_mbs_min2','mean'))
print(t.round(4)); print(p.groupby('assembly_type').my_mbs.agg(['size','mean']).round(4))
# before additive: which denominator?
b=d[(~d.is_dup_copy)&d.is_bound&(~d.is_glycan)]
print('before-add (dup/bound/glycan removed)',b.groupby('assembly_type').my_mbs.agg(['size','mean']).round(4))
