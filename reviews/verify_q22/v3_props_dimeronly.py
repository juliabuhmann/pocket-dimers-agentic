"""Independent verification: (2) property deltas, (3) dimer-only fraction & noise baseline. Status: first attempt."""
from pathlib import Path
import numpy as np, pandas as pd, glob, os
from scipy.stats import mannwhitneyu
WS=str(Path(__file__).resolve().parents[2]); W=f'{WS}/reviews/verify_q22/work'
full=pd.read_parquet(f'{WS}/data/part2/human_dimer_pockets_annotated.parquet')
# ---- (3a) full-data dimer-only fraction recomputed from max_overlap_iso column (their overlap values)
f=full[full.atlas_pass&full.is_sym_rep]
marg=(f.n_res_A+f.n_res_B>0)&((f.n_res_A.between(1,2)&(f.n_res_B==0))|(f.n_res_B.between(1,2)&(f.n_res_A==0)))
f=f[~marg]
itf=f[(f.n_res_A>=3)&(f.n_res_B>=3)]
print('FULL: interface n',len(itf),'dimer-only(OC<.5) frac',(itf.max_overlap_iso<0.5).mean(),'NaN OC',itf.max_overlap_iso.isna().mean(),
      ' without sym-rep filter:',(full[full.atlas_pass&full.is_interface].max_overlap_iso<0.5).mean())
# ---- sample
res=pd.read_parquet(W+'/res.parquet'); res['key']=res.mid+'|'+res.chain+':'+res.resnum.astype(str)
IF=set(res.key[res.is_if]); aa=dict(zip(res.key,res.resname))
samp=set(res.mid)
d=full[full.modelEntityId.isin(samp)&full.atlas_pass].copy()
d['rl']=d.residues.str.split()
iso=[]
for k in ('homo','hetero'):
    for fn in sorted(glob.glob(f'{WS}/data/part2/pockets/human_isolated/{k}/chunk_*[0-9].parquet')):
        c=pd.read_parquet(fn,columns=['structure_stem','residues','burrow_score','plddt_mean'])
        c=c[c.structure_stem.str[:19].isin(samp)]; iso.append(c)
iso=pd.concat(iso); iso['mid']=iso.structure_stem.str[:19]; iso['rl']=iso.residues.str.split(); iso['ok']=(iso.burrow_score>=0.03)&(iso.plddt_mean>=70)
print('iso pockets in sample',len(iso),iso.ok.sum())
isosets={m:[set(r) for r in g.rl] for m,g in iso[iso.ok].groupby('mid')}
isosets_all={m:[set(r) for r in g.rl] for m,g in iso.groupby('mid')}
dsets={m:[set(r) for r in g.rl] for m,g in d.groupby('modelEntityId')}
oc=lambda a,b:len(a&b)/min(len(a),len(b))
def maxoc(s,lst): return max([oc(s,t) for t in lst]+[0.0])
d['oc_iso']=[maxoc(set(r),isosets.get(m,[])) for r,m in zip(d.rl,d.modelEntityId)]
d['oc_iso_all']=[maxoc(set(r),isosets_all.get(m,[])) for r,m in zip(d.rl,d.modelEntityId)]
d['kind2']=np.where((d.n_res_A>=3)&(d.n_res_B>=3),'interface',np.where((d.n_res_A==0)|(d.n_res_B==0),'single','marginal'))
r=d[d.is_sym_rep&(d.kind2!='marginal')]
for kd in ('interface','single'):
    x=r[r.kind2==kd]; print(f'SAMPLE {kd}: n={len(x)} new(OC<0.5 vs atlas-pass iso)={np.mean(x.oc_iso<0.5):.3f}; vs all iso={np.mean(x.oc_iso_all<0.5):.3f}; their col={np.nanmean(x.max_overlap_iso<0.5):.3f}')
# reverse noise: isolated pockets with no interface residue lost in dimer
iso['n_if']=[sum(f'{m}|{x}' in IF for x in rl) for m,rl in zip(iso.mid,iso.rl)]
io=iso[iso.ok].copy()
io['oc_d']=[maxoc(set(rl),dsets.get(m,[])) for m,rl in zip(io.mid,io.rl)]
print('ISO no-interface pockets lost (OC<0.5 vs dimer atlas-pass):',np.mean(io[io.n_if==0].oc_d<0.5).round(3),'n',(io.n_if==0).sum(),'; with >=3 IF residues',np.mean(io[io.n_if>=3].oc_d<0.5).round(3),'n',(io.n_if>=3).sum())
# score-matched null: single-chain dimer-only rate in lowest vs highest score quintile
q=pd.qcut(r.burrow_score,5,labels=False); r=r.assign(q=q,new=r.oc_iso<0.5)
print(r.groupby(['kind2','q']).new.mean().unstack().round(3))
# ---- (2) property deltas
def cd(a,b):
    u=mannwhitneyu(a,b).statistic; return 2*u/(len(a)*len(b))-1
def comp(rl,m):
    n=[aa.get(f'{m}|{x}','XXX') for x in rl]; L=max(len(n),1)
    return (sum(x in('ALA','VAL','LEU','ILE','MET','PHE','TRP','CYS') for x in n)/L, sum(x in('PHE','TRP','TYR','HIS') for x in n)/L, sum(x in('ASP','GLU','LYS','ARG') for x in n)/L)
r=r.copy(); cc=np.array([comp(rl,m) for rl,m in zip(r.rl,r.modelEntityId)]); r['hyd'],r['aro'],r['chg']=cc.T; r['nres']=r.rl.map(len)
I=r[r.kind2=='interface']; S=r[r.kind2=='single']
print('PROPS real (sample, sym reps, atlas_pass; n_int',len(I),'n_single',len(S),')')
for c in ['volume_A3','nres','burrow_score','hyd','aro','chg','plddt_mean']:
    print(f'  {c}: median I {I[c].median():.3f} S {S[c].median():.3f}  Cliff delta {cd(I[c],S[c]):+.3f}')
# control: shifted run; "former-interface" pockets (>=50% IF res, >=3) vs other, same-chain geometry
cp=pd.read_parquet(W+'/shift_pockets.parquet'); cp['mid']=cp.structure.str.extract(r'(AF-\d+)')[0]
cp=cp[(cp.burrow_score>=0.03)&(cp.plddt_mean>=70)].copy(); cp['rl']=cp.residues.str.split()
cp['nif']=[sum(f'{m}|{x}' in IF for x in rl) for m,rl in zip(cp.mid,cp.rl)]; cp['nres']=cp.rl.map(len)
cp['fi']=cp.nif/cp.nres
cc=np.array([comp(rl,m) for rl,m in zip(cp.rl,cp.mid)]); cp['hyd'],cp['aro'],cp['chg']=cc.T
A=cp[(cp.nif>=3)&(cp.fi>=0.5)]; B=cp[cp.nif==0]
print('PROPS control (shifted; former-interface pockets n',len(A),'vs no-IF-residue pockets n',len(B),')')
for c in ['volume_A3','nres','burrow_score','hyd','aro','chg']:
    print(f'  {c}: median F {A[c].median():.3f} N {B[c].median():.3f}  Cliff delta {cd(A[c],B[c]):+.3f}')
