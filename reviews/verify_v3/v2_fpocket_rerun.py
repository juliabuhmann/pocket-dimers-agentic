"""Verifier: rerun fpocket on 10 random ok sites (seed 0) from raw mmCIF; compare selected-pocket druggability.
Own implementation (gemmi assembly 1, polymer only, fpocket -f default, pocket = most ligand heavy atoms within 3 A
of its alpha-sphere centres, ties -> lower pocket number)."""
from pathlib import Path
import os, re, subprocess, sys, tempfile, string, time
import numpy as np, pandas as pd, gemmi
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor
ROOT = str(Path(__file__).resolve().parents[2])
PDB_MIRROR = '/path/to/pdb_mirror'   # mirror named in docs/FPOCKET_DRUGGABILITY.md / src/part1/config.py
FPOCKET = f'{ROOT}/conda_envs/fpocket/bin/fpocket'
N_SITES, SEED, OVERLAP_A = 10, 0, 3.0
OUTDIR = f'{ROOT}/reviews/verify_v3/fpocket_runs'

def run(row):
    t0 = time.time()
    st = gemmi.read_structure(f'{PDB_MIRROR}/{row.pdb_id[1:3]}/{row.pdb_id}.cif.gz')
    st.setup_entities(); st.remove_hydrogens(); st.remove_alternative_conformations()
    asm = st.assemblies[0]
    model = gemmi.make_assembly(asm, st[0], gemmi.HowToNameCopiedChain.Short)
    comps = set(row.components); cen = np.asarray(row.centroid, float)
    cands = []
    for ch in model:
        for r in ch:
            if r.name in comps and r.het_flag == 'H':
                xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r])
                cands.append((np.linalg.norm(xyz.mean(0) - cen), xyz))
    cands.sort(key=lambda t: t[0])
    # nearest residue, plus further component residues of a multi-CCD ligand within 8 A of the census centroid
    lig = [cands[0][1]] + [x for dd, x in cands[1:] if len(comps) > 1 and dd < 8.0]
    lig = np.vstack(lig)
    ns = gemmi.Structure(); ns.cell = st.cell; m = gemmi.Model('1')
    names = iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    for ch in model:
        nc = gemmi.Chain(next(names)) if True else None
        for r in ch:
            if r.entity_type == gemmi.EntityType.Polymer:
                nc.add_residue(r)
        if len(nc): m.add_chain(nc)
    ns.add_model(m); ns.setup_entities()
    wd = f'{OUTDIR}/{row.site_id}'; os.makedirs(wd, exist_ok=True)
    pdb = f'{wd}/s.pdb'; ns.write_pdb(pdb)
    subprocess.run([FPOCKET, '-f', pdb], check=True, capture_output=True, timeout=1200)
    info = open(f'{wd}/s_out/s_info.txt').read()
    drug = {int(p): float(v) for p, v in re.findall(r'Pocket (\d+) :.*?Druggability Score : \t?\s*([-\d.]+)', info, re.S)}
    best = None
    for p in sorted(drug):
        cen_s = [list(map(float, (l[30:38], l[38:46], l[46:54]))) for l in open(f'{wd}/s_out/pockets/pocket{p}_vert.pqr') if l.startswith('ATOM')]
        cen_s = np.array(cen_s)
        dmin = np.sqrt(((lig[:, None, :] - cen_s[None]) ** 2).sum(-1)).min(1)
        n = int((dmin <= OVERLAP_A).sum())
        if best is None or n > best[1]: best = (p, n)
    return dict(site_id=row.site_id, my_pocket=best[0], my_overlap=best[1], my_n_lig=len(lig), my_drug=drug[best[0]],
                my_n_pockets=len(drug), secs=round(time.time() - t0, 1))

if __name__ == '__main__':
    fp = pd.read_parquet(f'{ROOT}/data/interim/fpocket/pdbsite_sites.parquet')
    fp = fp[fp.status == 'ok'].sort_values('site_id').reset_index(drop=True)
    pick = fp.sample(N_SITES, random_state=SEED)
    s = pd.read_parquet(f'{ROOT}/data/interim/pdbsite/sites.parquet', columns=['site_id', 'components', 'centroid'])
    rows = pick[['site_id', 'pdb_id']].merge(s, on='site_id')
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(run, [SimpleNamespace(**r) for r in rows.to_dict("records")]))
    out = pd.DataFrame(res).merge(pick[['site_id', 'stratum', 'fp_pocket', 'fp_overlap_n', 'fp_n_lig_atoms', 'fp_drug', 'fp_n_pockets', 'fp_n_atoms']], on='site_id')
    out['drug_match'] = np.isclose(out.my_drug, out.fp_drug, atol=1e-3)
    out.to_csv(f'{ROOT}/reviews/verify_v3/v2_fpocket_rerun.tsv', sep='\t', index=False)
    print(out.to_string())
