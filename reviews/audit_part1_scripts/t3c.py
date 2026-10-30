import pandas as pd
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['ccd','in_primary_set','ligand_class','assembly_type','is_mbs','biolip_status','marsh_class','is_h1_rep_instance'])
lc=pd.read_parquet('data/part1/ligand_classes.parquet',columns=['ccd','druglike_source','is_druglike_chembl','chembl_max_phase','n_heavy'])
p=d[d.in_primary_set].merge(lc,on='ccd',how='left')
dl=p[p.ligand_class=='drug_like']
print('primary drug_like instances', len(dl), dl.druglike_source.value_counts().to_dict())
print('distinct CCDs by source', dl.groupby('druglike_source').ccd.nunique().to_dict())
h=p[p.is_h1_rep_instance & (p.ligand_class=='drug_like')]; print('H1 reps drug_like', len(h), h.druglike_source.value_counts().to_dict())
print('ccd_classes druglike_source overall', lc.druglike_source.value_counts(dropna=False).to_dict())
# MBS by class among oligomers, chembl vs ro5
o=p[p.assembly_type!='monomer']
print(o[o.ligand_class=='drug_like'].groupby('druglike_source').is_mbs.agg(['size','mean']).round(3))
# not_covered fraction by assembly type
print(p.groupby('assembly_type').apply(lambda x:(x.biolip_status=='not_covered').mean()).round(3))
