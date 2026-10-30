import numpy as np, pandas as pd, gemmi, os
d=pd.read_parquet('data/interim/pdbsite/sites.parquet')
s=d.sample(20,random_state=0)
HY=set('AVLIMFWC'); CH=set('DEKR')
rows=[]
import sys
CUT=float(sys.argv[1])
for _,r in s.iterrows():
    p=f'/cluster/scratch/jjaenes/pdb/{r.pdb_id[1:3]}/{r.pdb_id}.cif.gz'
    st=gemmi.read_structure(p); st.setup_entities()
    try:
        a=[x for x in st.assemblies if x.name==str(r.assembly_id)][0]
        m=gemmi.make_assembly(a,st[0],gemmi.HowToNameCopiedChain.AddNumber)
    except Exception as e:
        m=st[0]
    m.remove_hydrogens(); 
    # first altloc only
    for ch in m:
        for res in ch:
            pass
    ns=gemmi.NeighborSearch(m,st.cell,6).populate()
    cen=gemmi.Position(*r.centroid)
    best=None
    for ch in m:
        for res in ch:
            if res.name!=r.ccd or res.het_flag!='H': continue
            c=res.calculate_center_of_mass() if False else gemmi.calculate_center_of_mass(res) if hasattr(gemmi,'calculate_center_of_mass') else None
            xs=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in res]); cc=xs.mean(0)
            dd=np.linalg.norm(cc-np.array(r.centroid))
            if best is None or dd<best[0]: best=(dd,ch,res)
    dd,ch,lig=best
    resid={}
    for la in lig:
        if la.altloc not in ('\0','A'): continue
        for mk in ns.find_atoms(la.pos,'\0',radius=5.5):
            cra=mk.to_cra(m)
            if cra.residue.is_water() or cra.residue.het_flag=='H': continue
            if cra.atom.altloc not in ('\0','A'): continue
            if gemmi.find_tabulated_residue(cra.residue.name) is None or not gemmi.find_tabulated_residue(cra.residue.name).is_amino_acid(): continue
            dist=cra.atom.pos.dist(la.pos)
            if dist<=CUT:
                resid[(cra.chain.name,str(cra.residue.seqid))]=cra.residue.name
    codes=[gemmi.find_tabulated_residue(n).one_letter_code.upper() for n in resid.values()]
    n=len(codes)
    rows.append(dict(site=r.site_id,pdb=r.pdb_id,ccd=r.ccd,dist=dd,n_mine=n,n_stored=None,
      hyd_mine=np.mean([c in HY for c in codes]) if n else np.nan,chg_mine=np.mean([c in CH for c in codes]) if n else np.nan,
      hyd_st=r.contact_frac_hydrophobic,chg_st=r.contact_frac_charged))
o=pd.DataFrame(rows)
cc=pd.read_parquet('data/part1/ligand_sites.parquet',columns=['ligand_uid','contacts']).set_index('ligand_uid')
o['n_stored']=[sum(c['kind']=='protein' for c in cc.loc[u,'contacts']) if u in cc.index else np.nan for u in s.ligand_uid]
o['dh']=o.hyd_mine-o.hyd_st; o['dc']=o.chg_mine-o.chg_st
o.to_csv(f'reviews/verify_misc/v1b_contacts_{CUT}.tsv',sep='\t',index=False)
print(o.round(3).to_string())
print('max|dh|',o.dh.abs().max(),'max|dc|',o.dc.abs().max(),'n exact (both)',((o.dh.abs()<1e-6)&(o.dc.abs()<1e-6)).sum(), 'corr', o[['hyd_mine','hyd_st']].corr().iloc[0,1], o[['chg_mine','chg_st']].corr().iloc[0,1])
