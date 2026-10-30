import importlib.util, glob, numpy as np, pandas as pd, gemmi
sp=glob.glob('.venvs/q22/lib/python*/site-packages/pydssp/pydssp_numpy.py')[0]
spec=importlib.util.spec_from_file_location('pn',sp); pn=importlib.util.module_from_spec(spec); spec.loader.exec_module(pn)
c=pd.read_csv('results/part2/q23a_curated.tsv',sep='\t',comment='#')
for _,r in c.sample(5,random_state=0).iterrows():
    p=f"/cluster/scratch/gmagnusson/dimers/afdb_human_dimers/{'homodimers' if r.kind=='homo' else 'heterodimers'}/{r.modelEntityId}.pdb"
    st=gemmi.read_structure(p); hel={}
    for ch in st[0]:
        xyz=[];nums=[]
        for res in ch:
            try: xyz.append([list(res[a][0].pos.tolist()) for a in ('N','CA','C','O')]); nums.append(res.seqid.num)
            except Exception: pass
        oh=pn.assign(np.array(xyz)); 
        for i,n in enumerate(nums): hel[(ch.name,n)]=oh[i,1]
    v=[hel.get((x.split(':')[0],int(x.split(':')[1])),np.nan) for x in r.residues.split()]
    print(r.pocket_id,'helix frac mine %.3f stored %.3f'%(np.nanmean(v),r.helix_frac_lining),'<0.7:',np.nanmean(v)<0.7)
