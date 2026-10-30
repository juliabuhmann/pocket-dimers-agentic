import pandas as pd, numpy as np, requests, glob, os
cols=['q','t','fident','alnlen','evalue','bits','qs','qe','ts','te','qlen','tlen','qseq','tseq'][:12]
def rd(p): 
    return pd.read_csv(p,sep='\t',header=None,usecols=range(6),names=['q','t','fident','alnlen','evalue','bits'])
rows=[]
for sp in sorted(set(os.path.basename(f).split('_')[1] for f in glob.glob('data/interim/q24/hits_*_fwd.tsv'))):
    f=rd(f'data/interim/q24/hits_{sp}_fwd.tsv'); r=rd(f'data/interim/q24/hits_{sp}_rev.tsv')
    f=f[f.evalue<=1e-3]; r=r[r.evalue<=1e-3]
    fb=f.sort_values(['q','bits'],ascending=[True,False]).drop_duplicates('q').set_index('q')
    rb=r.sort_values(['q','bits'],ascending=[True,False]).drop_duplicates('q').set_index('q')
    for h,row in fb.iterrows():
        if row.t in rb.index and rb.loc[row.t,'t']==h: rows.append((sp,h,row.t,row.fident,row.alnlen,row.evalue))
R=pd.DataFrame(rows,columns=['species','human','ortholog','fident','alnlen','evalue'])
print(R.groupby('species').size().to_dict())
cache={}
def gene(acc):
    if acc in cache: return cache[acc]
    try:
        t=requests.get(f'https://rest.uniprot.org/uniprotkb/{acc.split('-')[0]}.tsv?fields=gene_names,protein_name,organism_name',timeout=20).text.strip().split('\n')[1].split('\t')
    except Exception as e: t=['?','?','?']
    t=(t+['?','?','?'])[:3]
    cache[acc]=t; return t
def norm(g): return set(x.upper() for x in g.replace(';',' ').split())
def check(df):
    out=[]
    for _,r in df.iterrows():
        gh,ph,_=gene(r.human); go,po,oo=gene(r.ortholog)
        same=bool(norm(gh)&norm(go)); 
        out.append(dict(species=r.species,human=r.human,gene_h=gh,ortholog=r.ortholog,gene_o=go,prot_o=po[:60],fident=r.fident,alnlen=r.alnlen,same_gene=same))
    return pd.DataFrame(out)
five=R.sample(5,random_state=0); o5=check(five); print(o5.to_string()); o5.to_csv('reviews/verify_misc/v2b_rbh_5.tsv',sep='\t',index=False)
ext=R[R.species.isin(['mouse','rat','zebrafish','fly','worm','yeast','arabidopsis'])].sample(70,random_state=1)
oe=check(ext); print('extra 70 euk pairs: same gene name',oe.same_gene.mean()); print(oe.groupby('species').same_gene.agg(['mean','size']))
oe.to_csv('reviews/verify_misc/v2b_rbh_70.tsv',sep='\t',index=False)
