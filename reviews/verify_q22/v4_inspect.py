"""Independent verification (4): inspect random Q2.3a entries. Status: first attempt."""
from pathlib import Path
import numpy as np, pandas as pd, gemmi, freesasa, math
from scipy.spatial import cKDTree
import importlib.util,sys
WS=str(Path(__file__).resolve().parents[2])
spec=importlib.util.spec_from_file_location('v1',f'{WS}/reviews/verify_q22/v1_features.py'); v1=importlib.util.module_from_spec(spec); spec.loader.exec_module(v1)
t=pd.read_csv(f'{WS}/results/part2/q23a_dimer_only_ligandable.tsv',sep='\t',comment='#')
print(len(t))
ro=t.sample(10,random_state=0)
named=t[t.gene.isin(['RGS4','MYCN','MAX'])]
sel=pd.concat([named,ro]).drop_duplicates('accession')
full=pd.read_parquet(f'{WS}/data/part2/human_dimer_pockets_annotated.parquet')
fp=full.set_index('pocket_id')
def helix(phi,psi): return (-100<phi<-30) and (-80<psi<-5)
out=[]
for _,r in sel.iterrows():
    mid=r.best_modelEntityId; kind=r.best_kind
    st=gemmi.read_structure(f'{v1.ROOT}/{kind}dimers/{mid}.pdb'); st.remove_hydrogens(); m=st[0]
    res=r.best_residues.split()
    info={}; bf={}; hel={}; 
    for c in m:
        pp=c.get_polymer(); 
        for i,rr in enumerate(c):
            bf[(c.name,rr.seqid.num)]=np.mean([a.b_iso for a in rr])
            prev=c[i-1] if i>0 else None; nxt=c[i+1] if i<len(c)-1 else None
            if prev is not None and nxt is not None:
                phi,psi=gemmi.calculate_phi_psi(prev,rr,nxt); hel[(c.name,rr.seqid.num)]=helix(math.degrees(phi),math.degrees(psi)) if not (math.isnan(phi) or math.isnan(psi)) else False
    L=[(x.split(':')[0],int(x.split(':')[1])) for x in res]
    pl=[bf[k] for k in L]; hf=np.mean([hel.get(k,False) for k in L])
    # whole-chain plddt & disorder
    chs={c.name:np.array([np.mean([a.b_iso for a in rr]) for rr in c]) for c in m}
    # dSASA
    X=[];R=[];ri=[];ch=[];nres={}
    dat=[]
    for c in m:
        Xc=[];Rc=[];rc=[];n=0
        for rr in c:
            if rr.name not in v1.MAXASA: continue
            for a in rr: Xc.append(a.pos.tolist()); Rc.append(v1.RAD.get(a.element.name.upper(),1.7)); rc.append(n)
            n+=1
        dat.append((np.array(Xc),np.array(Rc),np.array(rc),n))
    sa=[v1.sasa_res(*d_) for d_ in [(d[0],d[1],d[2],d[3]) for d in dat]]
    Xall=np.vstack([d[0] for d in dat]);Rall=np.concatenate([d[1] for d in dat]);rall=np.concatenate([dat[0][2],dat[1][2]+dat[0][3]])
    sc=v1.sasa_res(Xall,Rall,rall,dat[0][3]+dat[1][3])
    bsa=(sa[0].sum()+sa[1].sum()-sc.sum())
    # buried interface area in A^2 (all, using raw areas)
    # chain-level: lengths
    nA=sum(1 for k in L if k[0]=='A'); nB=sum(1 for k in L if k[0]=='B')
    lowfrac={c:float((chs[c]<70).mean()) for c in chs}
    pid=r.best_pocket_id; p=fp.loc[pid] if pid in fp.index else None
    info=dict(acc=r.accession,gene=r.gene,partner=r.best_partner,kind=kind,vol=r.best_volume_A3,score=round(r.best_burrow_score,2),ipSAE=round(r.best_max_ipSAE,2),
      nA=nA,nB=nB,lin_plddt_mean=round(np.mean(pl)),lin_plddt_min=round(min(pl)),helix_frac_lining=round(hf,2),BSA=round(bsa),chainlen=[len(chs[c]) for c in chs],chain_frac_plddt_lt70=[round(lowfrac[c],2) for c in chs],oc_iso=round(r.best_max_overlap_iso,2),res_span_A=(min([k[1] for k in L if k[0]=='A'],default=0),max([k[1] for k in L if k[0]=='A'],default=0)),res_span_B=(min([k[1] for k in L if k[0]=='B'],default=0),max([k[1] for k in L if k[0]=='B'],default=0)),
      pocket_pts=None if p is None else int(p.n_points), chembl=r.chembl_status)
    out.append(info)
pd.set_option('display.width',250,'display.max_columns',40)
df=pd.DataFrame(out); print(df.to_string()); df.to_csv(f'{WS}/reviews/verify_q22/v4_inspect.tsv',sep='\t',index=False)
