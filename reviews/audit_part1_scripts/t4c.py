# Estimate same-UniProt "homomers" whose protein entities differ (p66/p51-like) and heteromers with identical entities.
import pandas as pd, gemmi, difflib
from multiprocessing import Pool
N=300; SEED=0; MINLEN=30
e=pd.read_parquet('data/part1/entries.parquet',columns=['pdb_id','assembly_type','skip_reason','uniprot_accs','n_chains_without_uniprot'])
e=e[e.skip_reason.isna()]
def ents(pdb):
    b=gemmi.cif.read(f'/path/to/pdb_mirror/{pdb[1:3]}/{pdb}.cif.gz').sole_block()
    t=b.find('_entity_poly.',['entity_id','type','pdbx_seq_one_letter_code_can'])
    seqs=[r.str(2).replace('\n','').replace(';','').strip() for r in t if 'polypeptide' in r.str(1)]
    return pdb,[s for s in seqs if len(s)>=MINLEN]
def classify(seqs):
    if len(seqs)<2: return 'single_entity'
    s0=seqs[0]; worst=1.0; lr=1.0
    for s in seqs[1:]:
        lr=min(lr,min(len(s),len(s0))/max(len(s),len(s0)))
        worst=min(worst,difflib.SequenceMatcher(None,s0,s,autojunk=False).ratio())
    if worst>0.95 and lr>0.95: return 'multi_entity_near_identical'
    if lr<0.9: return 'multi_entity_length_differs'
    return 'multi_entity_seq_differs'
if __name__=='__main__':
    for at in ['homomer','heteromer']:
        s=e[e.assembly_type==at].sample(N,random_state=SEED)
        with Pool(8) as pl: out=pl.map(ents,s.pdb_id)
        c=pd.Series({p:classify(q) for p,q in out}); print(at, c.value_counts().to_dict())
        if at=='homomer': print('length-differs examples', list(c[c=='multi_entity_length_differs'].index[:8]))
        if at=='heteromer':
            print('heteromer w/ all-identical-ish entities', list(c[c.isin(['single_entity','multi_entity_near_identical'])].index[:10]))
