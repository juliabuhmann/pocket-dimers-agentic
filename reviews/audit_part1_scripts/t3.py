# Audit task 3: additive / class label spot-check + aggregate checks
import sys, pandas as pd
SP=sys.argv[1]; SEED=0; N=20
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['ligand_uid','pdb_id','ccd','is_dup_copy','is_bound','is_glycan','is_additive',
  'biolip_status','additive_label_source','is_ligand_of_interest','in_primary_set','ligand_class','marsh_class','in_biolip_artifact_list','plinder_artifact','assembly_type','is_mbs'])
lc=pd.read_parquet('data/part1/ligand_classes.parquet',columns=['ccd','name'])
d=d.merge(lc,on='ccd',how='left')
pre=d[~d.is_dup_copy & d.is_bound & ~d.is_glycan]
show=['ligand_uid','ccd','name','ligand_class','biolip_status','additive_label_source','is_ligand_of_interest','in_biolip_artifact_list']
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',40)
print('PRIMARY sample'); print(pre[pre.in_primary_set].sample(N,random_state=SEED)[show].to_string(index=False))
print('ADDITIVE sample'); print(pre[pre.is_additive.fillna(False).astype(bool)].sample(N,random_state=SEED)[show].to_string(index=False))
p=pre[pre.in_primary_set]
print('primary by status/source'); print(p.groupby(['biolip_status','additive_label_source'],dropna=False).size())
print('primary on BioLiP artifact list frac', p.in_biolip_artifact_list.mean().round(4), 'n', p.in_biolip_artifact_list.sum())
top=p.ccd.value_counts().head(25); print('top primary CCDs'); print(pd.DataFrame({'n':top,'name':lc.set_index('ccd').name.reindex(top.index).str[:30],'on_list':p.groupby('ccd').in_biolip_artifact_list.first().reindex(top.index)}).to_string())
for c in ['GOL','SO4','EDO','PEG','ACT','CL','DMS','MPD','PO4','FMT','NA','MG','ZN','CA']:
    s=pre[pre.ccd==c]; print(c, len(s), 'kept', s.in_primary_set.mean().round(3))
add=pre[pre.is_additive.fillna(False).astype(bool)]
for c in ['FOL','FAD','NAD','NAP','HEM','ATP','ADP','ANP','SAM','SAH','PLP','FMN','COA','HEC','CLA','BCL','ZN','FE','MG','MN','CA','SF4','FES','BTN','TPP','GTP','GDP','STI','IMN']:
    s=pre[pre.ccd==c]
    if len(s): print(f'{c}: n={len(s)} excluded_as_additive={1-s.in_primary_set.mean():.3f} loi_among_excl={s[~s.in_primary_set].is_ligand_of_interest.mean():.2f}')
print('excluded but wwPDB ligand-of-interest:', add.is_ligand_of_interest.sum(), 'of', len(add))
print('drug_like excluded frac', 1-pre[pre.ligand_class=='drug_like'].in_primary_set.mean())
print('cofactor excluded frac', 1-pre[pre.ligand_class=='cofactor'].in_primary_set.mean())
