import pandas as pd, numpy as np
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['pdb_id','ccd','in_primary_set','assembly_type','n_protein_chains','sym_group','contacts','is_mbs','marsh_class','is_on_axis'],
   filters=[('in_primary_set','=',True)])
hd=d[(d.assembly_type=='homomer')&(d.n_protein_chains==2)]
g=hd.groupby(['pdb_id','ccd'])
pairs=[]
for (pdb,ccd),x in g:
    if len(x)!=2 or x.sym_group.nunique()==1: continue
    s=[set((c['unp_acc'],c['unp_pos']) for c in cs if c['kind']=='protein') for cs in x.contacts]
    j=len(s[0]&s[1])/max(1,len(s[0]|s[1])); pairs.append((pdb,ccd,j,x.is_mbs.sum()))
pp=pd.DataFrame(pairs,columns=['pdb','ccd','jac','n_mbs'])
print('non-collapsed 2-instance homodimer pairs', len(pp), 'Jaccard>=0.5:', round((pp.jac>=0.5).mean(),3), 'Jaccard>=0.8:', round((pp.jac>=0.8).mean(),3), 'jac==0', round((pp.jac==0).mean(),3))
print('of which both MBS', (pp.n_mbs==2).sum(), 'and jac>=0.5', ((pp.n_mbs==2)&(pp.jac>=0.5)).sum())
# missing UniProt positions in contacts (unp_pos None) -> hash on what?
allc=hd.contacts.explode().dropna()
print('homodimer contacts with unp_pos null frac', round(allc.map(lambda c:c['unp_pos'] is None).mean(),4))
ho=d[d.assembly_type=='homomer']; lc=pd.read_parquet('data/part1/ligand_classes.parquet',columns=['ccd','name']).set_index('ccd').name
chl=set(c for c in d.ccd.unique() if isinstance(lc.get(c),str) and ('CHLOROPHYLL' in lc[c].upper() or 'PHEOPHYTIN' in lc[c].upper()))
hc=ho[ho.marsh_class=='cofactor']; print('homomer cofactor chl frac', round(hc.ccd.isin(chl).mean(),3), 'MBS excl chl', round(hc[~hc.ccd.isin(chl)].is_mbs.mean(),3))
