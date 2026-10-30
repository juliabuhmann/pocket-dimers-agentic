import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import rankdata
f=pd.read_parquet('data/interim/q22/dimer_pocket_features_units.parquet')
m=pd.read_parquet('data/interim/iface_ligand_matched/pocket_match.parquet')
d=m.merge(f,on='pocket_id')
print(len(d), 'atlas_pass&sym_rep share',(d.atlas_pass&d.is_sym_rep).mean())
d=d[d.atlas_pass&d.is_sym_rep]
d['iface']=(d.n_res_A>=3)&(d.n_res_B>=3)
print('own iface vs is_interface agreement',(d.iface==d.is_interface).mean())
I=d[d.iface]
im=I[I.match_mbs_prim]; iu=I[~I.match_any]
D=pd.concat([im.assign(y=1),iu.assign(y=0)])
print('IM',len(im),im.unit.nunique(),'IU',len(iu),iu.unit.nunique())
def auc(y,s):
    r=rankdata(s); n1=y.sum(); n0=len(y)-n1
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def fit(X,y,lam=1.0):
    def nll(w):
        z=X@w[1:]+w[0]; return np.sum(np.logaddexp(0,z)-y*z)+0.5*lam*np.sum(w[1:]**2)
    return minimize(nll,np.zeros(X.shape[1]+1),method='L-BFGS-B').x
rng=np.random.default_rng(0)
def cv(cols,seed=0):
    X=D[cols].values.astype(float); y=D.y.values
    cl=D.unit.values; u=np.unique(cl); rng=np.random.default_rng(seed); fo=dict(zip(u,rng.permutation(len(u))%5))
    fold=np.array([fo[c] for c in cl]); oof=np.zeros(len(y)); fa=[]
    for k in range(5):
        tr=fold!=k; mu=X[tr].mean(0); sd=X[tr].std(0)+1e-9
        w=fit((X[tr]-mu)/sd,y[tr]); oof[~tr]=((X[~tr]-mu)/sd)@w[1:]+w[0]; fa.append(auc(y[~tr],oof[~tr]))
    return oof,fa
g=D.groupby('unit').indices; keys=list(g.keys())
for cols in (['frac_charged'],['burrow_score'],['frac_charged','burrow_score'],['volume_A3']):
    res=[]
    for s in range(3):
        oof,fa=cv(cols,s); res.append(auc(D.y.values,oof))
    oof,fa=cv(cols,0); bs=[]
    for b in range(500):
        ii=np.concatenate([g[keys[k]] for k in rng.integers(0,len(keys),len(keys))]); 
        yy=D.y.values[ii]
        if 0<yy.sum()<len(yy): bs.append(auc(yy,oof[ii]))
    print(cols,'AUC oof (3 fold-seeds)',np.round(res,3),'CI',np.round(np.percentile(bs,[2.5,97.5]),3),'fold range',round(min(fa),3),round(max(fa),3))
