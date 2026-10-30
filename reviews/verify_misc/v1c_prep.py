import numpy as np, pandas as pd, gemmi, os
d=pd.read_parquet('data/interim/pdbsite/sites.parquet')
c=d[(d.n_chains_asm>d.n_chains_kept)&(d.n_chains_asm<=24)&d.d_match]
print('sites where crop dropped chains:',(d.n_chains_asm>d.n_chains_kept).sum(),'of',len(d))
s=c.sample(10,random_state=0); out='reviews/verify_misc/crop'; os.makedirs(out,exist_ok=True)
for _,r in s.iterrows():
    st=gemmi.read_structure(f'/path/to/pdb_mirror/{r.pdb_id[1:3]}/{r.pdb_id}.cif.gz'); st.setup_entities()
    a=[x for x in st.assemblies if x.name==str(r.assembly_id)][0]
    m=gemmi.make_assembly(a,st[0],gemmi.HowToNameCopiedChain.AddNumber); m.remove_hydrogens()
    best=None
    for ch in m:
        for res in ch:
            if res.name==r.ccd and res.het_flag=='H':
                dd=np.linalg.norm(np.mean([[x.pos.x,x.pos.y,x.pos.z] for x in res],0)-np.array(r.centroid))
                if best is None or dd<best[0]: best=(dd,ch.name,res.seqid)
    # ligand pdb
    lig=gemmi.Structure(); lm=gemmi.Model('1'); lc=gemmi.Chain('Z')
    for ch in m:
        for res in ch:
            if ch.name==best[1] and res.seqid==best[2] and res.name==r.ccd: lc.add_residue(res)
    lm.add_chain(lc); lig.add_model(lm); lig.setup_entities(); lig.write_pdb(f'{out}/{r.site_id}_lig.pdb')
    # protein-only assembly
    nm=gemmi.Model('1')
    for ch in m:
        nc=gemmi.Chain(ch.name)
        for res in ch:
            if res.het_flag=='A' and not res.is_water(): nc.add_residue(res)
        if len(nc): nm.add_chain(nc)
    ns=gemmi.Structure(); ns.cell=st.cell; ns.spacegroup_hm=st.spacegroup_hm; ns.add_model(nm); ns.setup_entities()
    ns.make_mmcif_document().write_file(f'{out}/{r.site_id}_full.cif')
    print(r.site_id,r.pdb_id,r.ccd,int(r.n_chains_asm),int(r.n_chains_kept),flush=True)
s[['site_id','pdb_id','ccd']].to_csv(f'{out}/sites.tsv',sep='\t',index=False)
