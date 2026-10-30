import numpy as np, pandas as pd, gemmi
d=pd.read_parquet('data/interim/pdbsite/sites.parquet'); s=d.sample(20,random_state=0)
cc=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['ligand_uid','contacts']).set_index('ligand_uid')
for _,r in s.iloc[[0,1,2,3,5,7,8,13,14,16]].iterrows():
    st=gemmi.read_structure(f'/path/to/pdb_mirror/{r.pdb_id[1:3]}/{r.pdb_id}.cif.gz'); st.setup_entities()
    a=[x for x in st.assemblies if x.name==str(r.assembly_id)][0]
    m=gemmi.make_assembly(a,st[0],gemmi.HowToNameCopiedChain.AddNumber)
    ns=gemmi.NeighborSearch(m,st.cell,6).populate()
    best=None
    for ch in m:
        for res in ch:
            if res.name==r.ccd and res.het_flag=='H':
                dd=np.linalg.norm(np.mean([[x.pos.x,x.pos.y,x.pos.z] for x in res],0)-np.array(r.centroid))
                if best is None or dd<best[0]: best=(dd,res)
    lig=best[1]; res_={k:{} for k in ['vdw05','vdw05_noalt','4A']}
    for la in lig:
        for mk in ns.find_atoms(la.pos,'\0',radius=5.5):
            cra=mk.to_cra(m)
            if cra.residue.is_water() or cra.residue.het_flag=='H' and cra.residue.name!='MSE': continue
            t=gemmi.find_tabulated_residue(cra.residue.name)
            if t is None or not t.is_amino_acid(): continue
            dist=cra.atom.pos.dist(la.pos); k=(cra.chain.name,str(cra.residue.seqid))
            if dist<=la.element.vdw_r+cra.atom.element.vdw_r+0.5: res_['vdw05'][k]=1
            if dist<=4.0: res_['4A'][k]=1
    cs=cc.loc[r.ligand_uid,'contacts']; pr=[c for c in cs if c['kind']=='protein']
    print(r.pdb_id,r.ccd,'stored',len(pr),'n_vdw>=1',sum(c['n_vdw']>=1 for c in pr),'n_vdw>=2',sum(c['n_vdw']>=2 for c in pr),'mine vdw',len(res_['vdw05']),'4A',len(res_['4A']),'n_simple>=1',sum(c['n_simple']>=1 for c in pr),'kinds',set(c['kind'] for c in cs), 'nchain', len(set(c['chain'] for c in pr)))
