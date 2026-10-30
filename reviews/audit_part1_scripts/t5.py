import pandas as pd, numpy as np
cols=['pdb_id','ligand_uid','ccd','resolution','release_date','assembly_type','cluster30','in_primary_set','is_h1_rep_entry','is_h1_rep_instance','is_mbs','is_mbs_min2',
 'is_on_axis','on_axis_jaccard','c2_axis_dist','chainpair_rot_angle','sym_group','is_dup_copy','is_bound','biolip_status','chain_counts','n_protein_chains','mbs_contrib_symmetry_copy']
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=cols)
p=d[d.in_primary_set]
# independent H1 rep recomputation
k=p.dropna(subset=['cluster30'])
best=k.sort_values(['resolution','pdb_id']).groupby(['cluster30','ccd']).pdb_id.first()
flag=k[k.is_h1_rep_entry].groupby(['cluster30','ccd']).pdb_id.agg(lambda s:sorted(set(s)))
cmp=pd.DataFrame({'mine':best,'flag':flag})
cmp['ok']=[isinstance(f,list) and len(f)==1 and f[0]==m for m,f in zip(cmp.mine,cmp.flag)]
print('H1 groups', len(cmp), 'agree', round(cmp.ok.mean(),4), 'groups with >1 flagged entry', cmp.flag.map(lambda f:isinstance(f,list) and len(f)>1).sum(), 'no flag', cmp.flag.isna().sum())
print(cmp[~cmp.ok].head(5).to_string())
ri=k[k.is_h1_rep_instance].groupby(['cluster30','ccd']).size(); print('rep instances per group', ri.value_counts().head().to_dict(), 'total', ri.sum())
print('primary with null cluster30', p.cluster30.isna().sum())
# H1 rep ignores assembly type? groups spanning types
ty=k.groupby(['cluster30','ccd']).assembly_type.nunique(); print('cluster30xCCD groups spanning >1 assembly type', (ty>1).sum(), 'of', len(ty))
show=['ligand_uid','assembly_type','in_primary_set','is_dup_copy','biolip_status','is_mbs','is_mbs_min2','is_on_axis','on_axis_jaccard','c2_axis_dist','chainpair_rot_angle','sym_group','is_h1_rep_entry','chain_counts']
cases={'1hvr':['XK2'],'2rox':['T44'],'3tct':['3MI'],'1bm7':None,'3grs':['FAD'],'1hvy':['UMP','D16'],'2tsc':None,'1ql3':None,'4hhb':['HEM']}
pd.set_option('display.width',300); pd.set_option('display.max_colwidth',70)
for pdb,cc in cases.items():
    x=d[d.pdb_id==pdb]
    if cc: x=x[x.ccd.isin(cc)]
    else: x=x[x.is_bound & ~x.is_dup_copy]
    if len(x)==0: print(pdb,'absent'); continue
    x=x.assign(chain_counts=x.chain_counts.map(lambda c:';'.join(f"{a['chain']}:{a['n_res']}" for a in c)), sym_group=x.sym_group.str[-6:])
    print(x[show].round(2).to_string(index=False))
