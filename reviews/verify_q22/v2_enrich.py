"""Independent verification: enrichment real vs partner-shifted control. Status: first attempt."""
from pathlib import Path
import numpy as np, pandas as pd, os, sys
WS=str(Path(__file__).resolve().parents[2]); W=f'{WS}/reviews/verify_q22/work'
SURF=float(sys.argv[1]) if len(sys.argv)>1 else 0.2
res=pd.read_parquet(W+'/res.parquet'); res['key']=res.mid+'|'+res.chain+':'+res.resnum.astype(str)
res['surf']=res.rsasa_dimer_chain_free>=SURF
IF=set(res.key[res.is_if]); SURFK=set(res.key[res.surf])
real=pd.read_parquet(f'{WS}/data/part2/human_dimer_pockets_annotated.parquet')
samp=set(res.mid); real=real[real.modelEntityId.isin(samp)&real.atlas_pass].copy()
ctrl=pd.read_parquet(W+'/shift_pockets.parquet') if os.path.exists(W+'/shift_pockets.parquet') else None
def prep(df,mcol):
    df=df[(df.burrow_score>=0.03)&(df.plddt_mean>=70)].copy()
    df['mid']=df[mcol]; df['res']=df.residues.str.split()
    df['rk']=[[f'{m}|{r}' for r in rs] for m,rs in zip(df.mid,df.res)]
    df['nA']=df.res.map(lambda r:sum(x.startswith('A:') for x in r)); df['nB']=df.res.map(lambda r:sum(x.startswith('B:') for x in r))
    df['nif']=df['rk'].map(lambda k:sum(x in IF for x in k)); df['fif']=df.nif/df.res.map(len)
    return df
def stats(df,name):
    exp=len(IF&SURFK)/len(SURFK)
    p_def1=((df.nA>=3)&(df.nB>=3)).mean(); p_def2=((df.nif>=3)&(df.fif>=0.5)).mean()
    lin=set(k for ks in df['rk'] for k in ks)&SURFK
    rfrac=len(lin&IF)/len(lin)
    surf_rate_if=len(lin&IF)/len(IF&SURFK); surf_rate_non=len(lin-IF)/len(SURFK-IF)
    print(f'{name}: n_pockets={len(df)} exp_if_share_surface={exp:.3f} | pocket-level: def1(>=3/chain) {p_def1:.3f} -> enr {p_def1/exp:.2f}; def2(>=50% IF res,>=3) {p_def2:.3f} -> enr {p_def2/exp:.2f} | residue-level: if share of lining surf res {rfrac:.3f} -> enr {rfrac/exp:.2f}; rate if {surf_rate_if:.3f} vs non-if {surf_rate_non:.3f}')
    return df
print('surface cut',SURF,'n models',len(samp),'surf res',len(SURFK),'interface res',len(IF),'IF&surf',len(IF&SURFK))
real=stats(prep(real,'modelEntityId'),'REAL (their BURROW pockets)')
if ctrl is not None:
    ctrl['mid']=ctrl.structure.str.extract(r'(AF-\d+)')[0]
    c=prep(ctrl,'mid'); stats(c,'CONTROL (partner shifted >=15A gap)')
    print('models w/ control pockets',c.mid.nunique(),'real',real.mid.nunique())
