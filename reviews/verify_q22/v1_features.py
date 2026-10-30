"""Independent verification (Q2.2): sample 1000 models (seed 0), per-residue rSASA (isolated chain), interface (<=5 A heavy atom),
and write partner-translated (+20 A) PDBs for a BURROW control. Status: first attempt."""
from pathlib import Path
import sys, os, numpy as np, pandas as pd, gemmi, freesasa
from scipy.spatial import cKDTree
from multiprocessing import Pool
WS=str(Path(__file__).resolve().parents[2])
OUT=f'{WS}/reviews/verify_q22/work'; os.makedirs(OUT+'/shift',exist_ok=True)
ROOT='/cluster/scratch/gmagnusson/dimers/afdb_human_dimers'
SEED=0; N=1000; IF_CUT=5.0; SHIFT=20.0; GAP=15.0  # shift along centroid axis from 20 A in 5 A steps until min inter-chain distance >=15 A
RAD={'C':1.70,'N':1.55,'O':1.52,'S':1.80,'SE':1.90}   # own element radii (Bondi-like)
MAXASA={'ALA':129,'ARG':274,'ASN':195,'ASP':193,'CYS':167,'GLN':225,'GLU':223,'GLY':104,'HIS':224,'ILE':197,'LEU':201,'LYS':236,'MET':224,'PHE':240,'PRO':159,'SER':155,'THR':172,'TRP':285,'TYR':263,'VAL':174}  # Tien et al. 2013 theoretical
def sasa_res(coords,radii,resid,nres):
    r=freesasa.calcCoord(coords.ravel().tolist(),radii.tolist())
    a=np.array([r.atomArea(i) for i in range(len(radii))])
    return np.bincount(resid,weights=a,minlength=nres)
def work(a):
    mid,kind=a
    p=f'{ROOT}/{kind}dimers/{mid}.pdb'
    st=gemmi.read_structure(p); st.remove_hydrogens(); m=st[0]
    chs=[c for c in m]; 
    if len(chs)!=2: return None
    dat=[]
    for c in chs:
        X=[];R=[];ri=[];res=[]
        for k,r in enumerate(c):
            if r.name not in MAXASA: continue
            res.append((c.name,r.seqid.num,r.name))
            for at in r:
                X.append(at.pos.tolist()); R.append(RAD.get(at.element.name.upper(),1.7)); ri.append(len(res)-1)
        dat.append((np.array(X),np.array(R),np.array(ri),res))
    rows=[]
    tr=[cKDTree(d[0]) for d in dat]
    for i in (0,1):
        X,R,ri,res=dat[i]; o=1-i
        s=sasa_res(X,R,ri,len(res))
        d,_=tr[o].query(X)
        isif=np.zeros(len(res),bool); isif[np.unique(ri[d<=IF_CUT])]=True
        for j,(cn,num,nm) in enumerate(res):
            rows.append((mid,cn,num,nm,s[j]/MAXASA[nm],isif[j]))
    # shifted copy
    cA=dat[0][0].mean(0); cB=dat[1][0].mean(0); v=(cB-cA); v/=np.linalg.norm(v)
    SH=SHIFT
    while tr[0].query(dat[1][0]+SH*v)[0].min()<GAP and SH<400: SH+=5.0
    dmin=tr[0].query(dat[1][0]+SH*v)[0].min()
    st2=gemmi.read_structure(p); ch=st2[0][1]
    for r in ch:
        for at in r: at.pos=gemmi.Position(at.pos.x+SH*v[0],at.pos.y+SH*v[1],at.pos.z+SH*v[2])
    st2.write_pdb(f'{OUT}/shift/{mid}.pdb')
    return rows,(mid,dmin,SH)
if __name__=='__main__':
    d=pd.read_parquet(f'{WS}/data/part2/human_dimer_pockets_annotated.parquet',columns=['modelEntityId','kind'])
    m=d.drop_duplicates('modelEntityId').sample(N,random_state=SEED).reset_index(drop=True)
    m.to_csv(OUT+'/sample.csv',index=False)
    with Pool(6) as pool: res=pool.map(work,[(a,b) for a,b in zip(m.modelEntityId,m.kind)],chunksize=10)
    res=[r for r in res if r]
    pd.DataFrame([x for r in res for x in r[0]],columns=['mid','chain','resnum','resname','rsasa_dimer_chain_free','is_if']).to_parquet(OUT+'/res.parquet')
    pd.DataFrame([r[1] for r in res],columns=['mid','min_dist_after_shift','shift']).to_csv(OUT+'/shift_dist.csv',index=False)
    print(len(res))
