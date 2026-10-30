import numpy as np, pandas as pd
d=pd.read_parquet('data/interim/pdbsite/sites.parquet')
rng=np.random.default_rng(12345)
def wcd(x1,w1,x0,w0):
    # weighted Cliff's delta via sorting + cumulative weights
    o=np.argsort(x0); s=x0[o]; cw=np.concatenate([[0],np.cumsum(w0[o])]); W0=cw[-1]
    lo=np.searchsorted(s,x1,'left'); hi=np.searchsorted(s,x1,'right')
    less=cw[lo]; greater=W0-cw[hi]
    return float((w1*(less-greater)).sum()/(w1.sum()*W0))
cols={'volume (merged_volume_A3)':'merged_volume_A3','charged (merged_frac_charged)':'merged_frac_charged',
 'hydrophobic (merged_frac_hydrophobic)':'merged_frac_hydrophobic','burrow_score (located member)':'burrow_score'}
cl=d.cluster30.values; ucl,inv=np.unique(cl,return_inverse=True)
idx=[np.where(inv==k)[0] for k in range(len(ucl))]
B=1000; out=[]
for name,c in cols.items():
    m=d[c].notna().values
    x=d[c].values; w=d.w.values; y=d.y_mbs.values
    pt=wcd(x[m&(y==1)],w[m&(y==1)],x[m&(y==0)],w[m&(y==0)])
    bs=[]
    for b in range(B):
        pick=rng.integers(0,len(ucl),len(ucl))
        ii=np.concatenate([idx[k] for k in pick]); ii=ii[m[ii]]
        g1=ii[y[ii]==1]; g0=ii[y[ii]==0]
        bs.append(wcd(x[g1],w[g1],x[g0],w[g0]))
    lo,hi=np.percentile(bs,[2.5,97.5]); out.append((name,pt,lo,hi)); print(name,round(pt,4),round(lo,3),round(hi,3),flush=True)
pd.DataFrame(out,columns=['metric','delta','lo','hi']).to_csv('reviews/verify_misc/v1a_cliffs.tsv',sep='\t',index=False)
