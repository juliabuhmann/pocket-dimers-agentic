# Audit task 2: independent contact recomputation for 40 stratified random primary-set ligand instances.
# Own implementation: gemmi assembly generators + numpy distances; BioLiP rule d <= vdW_i + vdW_j + TOL.
import sys, gzip, numpy as np, pandas as pd, gemmi
from multiprocessing import Pool
PDB_ROOT = '/cluster/scratch/jjaenes/pdb'
TOL = 0.5          # BioLiP tolerance (DESIGN §0)
MIN_PROT_LEN = 30  # eligible chain SEQRES length (DESIGN §0)
SEED = 0; N_PER = 10
SP = sys.argv[1]
d = pd.read_pickle(f'{SP}/ls.pkl')
full = pd.read_parquet('data/part1/ligand_sites.parquet',
    columns=['ligand_uid','label_asym','auth_seq','op','centroid','components','n_components'])
d = d.merge(full, on='ligand_uid')
p = d[d.in_primary_set]
strata = {'homomer_MBS': p[(p.assembly_type=='homomer')&p.is_mbs], 'homomer_SBS': p[(p.assembly_type=='homomer')&~p.is_mbs],
          'heteromer_MBS': p[(p.assembly_type=='heteromer')&p.is_mbs], 'monomer': p[p.assembly_type=='monomer']}
samp = pd.concat([v.sample(N_PER, random_state=SEED).assign(stratum=k) for k,v in strata.items()])

def vdw(el):
    r = gemmi.Element(el).vdw_r
    return r if r > 0 else gemmi.Element('C').vdw_r

def run(row):
    pdb = row['pdb_id']
    st = gemmi.read_structure(f'{PDB_ROOT}/{pdb[1:3]}/{pdb}.cif.gz')
    st.setup_entities(); st.remove_hydrogens(); st.remove_alternative_conformations()
    model = st[0]
    ent_ok = {}
    for e in st.entities:
        ok = e.entity_type == gemmi.EntityType.Polymer and e.polymer_type in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD) \
             and len(e.full_sequence) >= MIN_PROT_LEN
        for s in e.subchains: ent_ok[s] = ok
    # residues by subchain
    by_sub = {}
    for ch in model:
        for r in ch:
            by_sub.setdefault(r.subchain, []).append((ch.name, r))
    if st.assemblies:
        asm = st.assemblies[0]
        gens = [(g.subchains, list(g.operators)) for g in asm.generators]
    else:
        gens = [(list(by_sub), [None])]
    prot_xyz, prot_r, prot_key = [], [], []
    lig_cands = []
    for subs, opers in gens:
        for op in opers:
            tr = op.transform if op is not None else gemmi.Transform()
            oname = op.name if op is not None else '1'
            for s in subs:
                for chname, r in by_sub.get(s, []):
                    xyz = np.array([tr.apply(a.pos).tolist() for a in r])
                    if ent_ok.get(s, False):
                        prot_xyz.append(xyz); prot_r += [vdw(a.element.name) for a in r]
                        prot_key += [((s, oname), r.seqid.num, r.seqid.icode)]*len(r)
                    elif s == row['label_asym'] and str(r.seqid.num)+(r.seqid.icode.strip()) == str(row['auth_seq']) and r.name == row['ccd']:
                        lig_cands.append((oname, xyz, [vdw(a.element.name) for a in r]))
                    elif r.name in (row['components'] if row['components'] is not None else []) and row['n_components'] > 1:
                        lig_cands.append(('__comp__'+oname, xyz, [vdw(a.element.name) for a in r]))
    P = np.vstack(prot_xyz); PR = np.array(prot_r)
    main = [c for c in lig_cands if not c[0].startswith('__comp__')]
    if not main: return dict(ligand_uid=row['ligand_uid'], status='ligand_not_found')
    cen = np.array(row['centroid'])
    best = min(main, key=lambda c: np.linalg.norm(c[1].mean(0)-cen))
    L, LR = best[1], np.array(best[2])
    if row['n_components'] > 1:   # merge covalently attached components (atoms within 2.0 A), iteratively
        comps = [c for c in lig_cands if c[0] == '__comp__'+best[0]]
        changed = True
        while changed:
            changed = False
            for c in list(comps):
                if np.min(np.linalg.norm(c[1][:,None]-L[None], axis=2)) < 2.0:
                    L = np.vstack([L, c[1]]); LR = np.concatenate([LR, c[2]]); comps.remove(c); changed = True
    cdiff = float(np.linalg.norm(L.mean(0)-cen))
    D = np.linalg.norm(L[:,None,:]-P[None,:,:], axis=2)
    hit = D <= LR[:,None] + PR[None,:] + TOL
    pj = np.where(hit.any(0))[0]
    res = {}
    for j in pj: res[prot_key[j]] = res.get(prot_key[j], 0) + int(hit[:,j].sum())
    chains = {}
    for (ch, num, ic), n in res.items(): chains.setdefault(ch, []).append(n)
    nch1 = sum(1 for v in chains.values() if len(v) >= 1)
    nch2 = sum(1 for v in chains.values() if len(v) >= 2)
    return dict(ligand_uid=row['ligand_uid'], status='ok', centroid_diff=round(cdiff,2), n_lig_atoms=len(L),
                my_n_res=len(res), my_nch=nch1, my_mbs=nch1 >= 2, my_mbs_min2=nch2 >= 2,
                my_chain_counts=';'.join(f'{k[0]}/{k[1]}:{len(v)}' for k,v in sorted(chains.items())))

if __name__ == '__main__':
    rows = samp.to_dict('records')
    with Pool(8) as pool: out = pool.map(run, rows)
    o = pd.DataFrame(out).merge(samp[['ligand_uid','stratum','pdb_id','ccd','n_chains_contacted','is_mbs','is_mbs_min2','chain_counts','n_components']], on='ligand_uid')
    o['pipe_n_res'] = o.chain_counts.map(lambda c: sum(x['n_res'] for x in c))
    o['chain_counts'] = o.chain_counts.map(lambda c: ';'.join(f"{x['chain']}:{x['n_res']}" for x in c))
    o['agree_nch'] = o.my_nch == o.n_chains_contacted
    o['agree_mbs'] = o.my_mbs == o.is_mbs
    o.to_csv(f'{SP}/t2_out.tsv', sep='\t', index=False)
    print(o.groupby('stratum')[['agree_nch','agree_mbs']].mean())
    print('overall nch', o.agree_nch.mean(), 'mbs', o.agree_mbs.mean(), 'n_res agree', (o.my_n_res==o.pipe_n_res).mean())
    print(o.loc[~o.agree_nch | (o.my_n_res!=o.pipe_n_res), ['stratum','ligand_uid','n_components','centroid_diff','my_n_res','pipe_n_res','my_chain_counts','chain_counts']].to_string())
