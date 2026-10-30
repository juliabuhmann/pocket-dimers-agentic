import sys, pandas as pd
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['pdb_id','ccd','is_dup_copy','is_bound','is_glycan','is_additive','biolip_status','in_primary_set','ligand_class','in_biolip_artifact_list','plinder_artifact','biolip_list_geom_artifact'])
lc=pd.read_parquet('data/part1/ligand_classes.parquet',columns=['ccd','name','chembl_max_phase'])
pre=d[~d.is_dup_copy & d.is_bound & ~d.is_glycan].copy()
ex=pre[~pre.in_primary_set]
t=ex[~ex.in_biolip_artifact_list].ccd.value_counts().head(20)
print('top excluded CCDs NOT on BioLiP list'); print(pd.DataFrame({'n':t}).join(lc.set_index('ccd')[['name']]).assign(name=lambda x:x.name.str[:35]).to_string())
t=ex[ex.ligand_class=='drug_like'].ccd.value_counts().head(15)
print('top excluded drug_like'); print(pd.DataFrame({'n':t}).join(lc.set_index('ccd')[['name','chembl_max_phase']]).assign(name=lambda x:x.name.str[:35]).to_string())
for c in ['GOL','SO4','CL','PO4','EDO']:
    s=pre[(pre.ccd==c)&pre.in_primary_set]; print(c,'kept by status', s.biolip_status.value_counts().to_dict())
nc=pre[pre.biolip_status=='not_covered']
print('not_covered: kept frac', nc.in_primary_set.mean().round(3), ' covered kept frac', pre[pre.biolip_status!='not_covered'].in_primary_set.mean().round(3))
for c in ['GOL','SO4','EDO','CL','NA','MG']:
    a=pre[(pre.ccd==c)]
    print(c,'kept frac covered', round(a[a.biolip_status!='not_covered'].in_primary_set.mean(),3),'not_covered', round(a[a.biolip_status=='not_covered'].in_primary_set.mean(),3))
dl=pre[pre.ligand_class=='drug_like'].ccd.value_counts().head(25)
print('top drug_like CCDs (all)'); print(pd.DataFrame({'n':dl}).join(lc.set_index('ccd')[['name','chembl_max_phase']]).assign(name=lambda x:x.name.str[:30]).to_string())
