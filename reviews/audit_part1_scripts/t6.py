import pandas as pd, numpy as np
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['pdb_id','ccd','is_dup_copy','is_bound','is_glycan','biolip_status','release_date','method_class','n_components','is_covalent','marsh_class','ligand_class','is_ligand_of_interest','n_protein_chains'])
b=pd.read_csv('data/raw/part1/BioLiP.txt.gz',sep='\t',header=None,usecols=[0],dtype=str,quoting=3,na_filter=False)[0].str.lower()
bp=set(b)
pre=d[~d.is_dup_copy & d.is_bound & ~d.is_glycan].copy()
pre['entry_in_biolip']=pre.pdb_id.isin(bp)
cov=pre[pre.biolip_status!='not_covered']
x=cov[~cov.entry_in_biolip]
print('covered bound instances in entries absent from BioLiP:', len(x), 'entries', x.pdb_id.nunique(), 'by method', x.groupby('method_class').pdb_id.nunique().to_dict())
print(' classes there', x.ligand_class.value_counts().head(8).to_dict(), 'LOI', x.is_ligand_of_interest.sum())
print(' release year', x.drop_duplicates('pdb_id').release_date.str[:4].value_counts().sort_index().tail(6).to_dict())
print(' n_protein_chains median', x.drop_duplicates('pdb_id').n_protein_chains.median(), 'vs in-biolip', cov[cov.entry_in_biolip].drop_duplicates('pdb_id').n_protein_chains.median())
for k in ['n_components','is_covalent']:
    print(k, cov.groupby(cov[k] if k=='is_covalent' else (cov[k]>1)).apply(lambda g:(g.biolip_status=='relevant').mean()).round(3).to_dict())
# relevance for classes that should be nearly always relevant
for c in ['HEM','FAD','NAD','ATP','CLA','SF4','ZN']:
    g=cov[cov.ccd==c]; print(c,'relevant frac', round((g.biolip_status=='relevant').mean(),3), 'n', len(g), 'missing-entry frac among additive', round((~g[g.biolip_status=='additive'].entry_in_biolip).mean(),3))
