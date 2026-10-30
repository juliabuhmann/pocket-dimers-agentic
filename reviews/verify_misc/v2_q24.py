import numpy as np, pandas as pd
d=pd.read_parquet('data/interim/q24/q24_pocket_scores.parquet')
print(d.cls.value_counts().to_dict(), d['sub'].value_counts().to_dict(), d.unit.nunique(), d.kind.value_counts().to_dict())
# own quintiles over unique pockets
u=d.drop_duplicates('pocket_id')
edges=np.quantile(u.burrow_score,[.2,.4,.6,.8])
d['q']=np.searchsorted(edges,d.burrow_score,side='right')
print('agreement with stored score_q:',(d.q.astype(str).radd('Q')==d.score_q.str.replace('Q','Q').map(lambda s:s)).mean() if False else pd.crosstab(d.q,d.score_q).values.tolist())
d['y']=(d.oc_atlas>=0.5).astype(float)
def mdiff(x):
    a=x[x.cls=='interface']; b=x[x.cls=='single']
    pa=a.groupby('q').y.agg(['mean','size']); pb=b.groupby('q').y.mean()
    w=pa['size']/pa['size'].sum(); ok=pb.index.intersection(pa.index)
    return float((w[ok]/w[ok].sum()*(pa['mean'][ok]-pb[ok])).sum()), float(a.y.mean()), float(b.y.mean())
rng=np.random.default_rng(7)
units=d.unit.values; uu,inv=np.unique(units,return_inverse=True)
order=np.argsort(inv); bounds=np.searchsorted(inv[order],np.arange(len(uu)+1))
def boot(x_all, sel, B=1000):
    pt=mdiff(x_all[sel])[0]
    # pre-split by unit for speed
    sub=x_all[sel].reset_index(drop=True)
    g=sub.groupby('unit').indices; keys=list(g.keys()); res=[]
    for b in range(B):
        pick=rng.integers(0,len(keys),len(keys))
        ii=np.concatenate([g[keys[k]] for k in pick])
        res.append(mdiff(sub.iloc[ii])[0])
    return pt,*np.percentile(res,[2.5,97.5])
for name,sel in [('ALL',d.cls.isin(['interface','single'])),('homo',(d.kind=='homo')),('hetero',(d.kind=='hetero'))]:
    s=sel&d.cls.isin(['interface','single'])
    r=boot(d,s,500); raw=mdiff(d[s])
    print(name,'matched %.3f [%.3f, %.3f]'%r,'raw iface %.3f single %.3f'%raw[1:],'n',s.sum(),flush=True)
