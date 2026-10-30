# Audit task 4: pathologies
import pandas as pd, numpy as np
d=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['pdb_id','ccd','in_primary_set','ligand_class','marsh_class','assembly_type','is_mbs','n_protein_chains',
   'mbs_contrib_same_entity','mbs_contrib_same_protein','sym_group','is_on_axis','cluster30','is_h1_rep_instance','chain_counts','method_class'])
e=pd.read_parquet('data/part1/entries.parquet',columns=['pdb_id','assembly_type','skip_reason','n_distinct_proteins','n_chains_without_uniprot','any_chimeric','n_protein_chains_asm'])
lc=pd.read_parquet('data/part1/ligand_classes.parquet',columns=['ccd','name'])
p=d[d.in_primary_set]
t=p.pdb_id.value_counts().head(10)
top=p[p.pdb_id.isin(t.index)].groupby('pdb_id').agg(n=('ccd','size'),type=('assembly_type','first'),top_ccd=('ccd',lambda s:s.value_counts().index[0]),mbs=('is_mbs','mean'))
print('top10 entries', top.sort_values('n',ascending=False).round(2).to_string()); print('share top10', round(t.sum()/len(p),4), 'top100', round(p.pdb_id.value_counts().head(100).sum()/len(p),4))
for at in ['homomer','heteromer']:
    vc=p[p.assembly_type==at].pdb_id.value_counts(); print(at,'top1% entries share', round(vc.head(int(len(vc)*0.01)).sum()/vc.sum(),3), 'n entries', len(vc))
# heteromer heme / cofactor composition
hh=p[(p.assembly_type=='heteromer')&(p.ligand_class=='heme_porphyrin')]
print('heteromer heme_porphyrin n', len(hh)); print(hh.ccd.value_counts().head(8).to_dict())
CHL={'CLA','CHL','BCL','BCB','F6C','CL7','CL0','PHO','BPH','BPB','KC1','CHA','G9R'}  # chlorophyll/pheophytin codes; inspect names below
names=lc.set_index('ccd').name
chl=set(c for c in p.ccd.unique() if isinstance(names.get(c),str) and ('CHLOROPHYLL' in names[c].upper() or 'PHEOPHYTIN' in names[c].upper()))
print('chlorophyll-like codes found by name', sorted(chl)[:30])
hc=p[(p.assembly_type=='heteromer')&(p.marsh_class=='cofactor')]
print('heteromer cofactor n',len(hc),'chl frac',round(hc.ccd.isin(chl).mean(),3),'MBS all',round(hc.is_mbs.mean(),3),'MBS excl chl',round(hc[~hc.ccd.isin(chl)].is_mbs.mean(),3))
for at in ['homomer','heteromer']:
    x=p[(p.assembly_type==at)]
    r=x[x.is_h1_rep_instance]; print(at,'H1-rep MBS by marsh', r.groupby('marsh_class').is_mbs.agg(['size','mean']).round(3).to_dict())
# same-UniProt homomer with different entities contributing
hm=p[(p.assembly_type=='homomer')&p.is_mbs]
print('homomer MBS instances with contributors from different entities', (hm.mbs_contrib_same_entity==False).sum(), 'of', len(hm), 'entries', hm[hm.mbs_contrib_same_entity==False].pdb_id.nunique())
print('homomer MBS same_protein False', (hm.mbs_contrib_same_protein==False).sum())
ee=e[e.skip_reason.isna()]
print('entries types', ee.assembly_type.value_counts().to_dict())
print('heteromers with chains lacking UniProt', ((ee.assembly_type=='heteromer')&(ee.n_chains_without_uniprot>0)).sum())
print('any_chimeric by type', ee.groupby('assembly_type').any_chimeric.sum().to_dict())
# symmetric copy collapse in homodimers
hd=p[(p.assembly_type=='homomer')&(p.n_protein_chains==2)]
g=hd.groupby(['pdb_id','ccd']).agg(n=('sym_group','size'),nsg=('sym_group','nunique'))
g2=g[g.n==2]
print('homodimer (pdb,ccd) with exactly 2 instances:', len(g2), 'collapsed to 1 sym_group:', round((g2.nsg==1).mean(),3))
hdm=hd[hd.is_mbs & ~hd.is_on_axis.fillna(False)]
gm=hdm.groupby(['pdb_id','ccd']).agg(n=('sym_group','size'),nsg=('sym_group','nunique'))
gm2=gm[gm.n==2]; print('homodimer off-axis MBS pairs (n=2):',len(gm2),'collapsed:',round((gm2.nsg==1).mean(),3))
print('homomer primary: instances', (p.assembly_type=='homomer').sum(), 'distinct sym_groups', p[p.assembly_type=='homomer'].sym_group.nunique())
