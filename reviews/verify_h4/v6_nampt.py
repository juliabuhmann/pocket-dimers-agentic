"""Verifier: NAMPT (P43490) diagnostic. Superpose PDB ligand atoms into the dimer model and measure (i) clashes of the
transformed ligand with model heavy atoms (dimer vs the binding chain alone), (ii) per-chain share of binding residues.
Sites: the 5 seed-0 NAMPT sites + 2 spot-check misses from dcc_recompute.tsv. UNVERIFIED first attempt."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, gemmi
from scipy.spatial import cKDTree
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v2_dcc import parse_sifts, build_chains, read_model_ca, kabsch, LOCAL_R
OUT = Path(__file__).resolve().parent
CLASH = 2.5   # Å, ligand heavy atom to model heavy atom
R = pd.read_csv(OUT / "dcc_recompute.tsv", sep="\t")
S = R[(R.acc1 == "P43490")]
for _, r in S.iterrows():
    sm = parse_sifts(r.pdb_id); st, chains = build_chains(r.pdb_id)
    f = r.ligand_uid.split(":"); lig = []
    for c in st[0]:
        for res in c:
            if res.subchain == f[1] and res.name == f[3] and str(res.seqid.num) == f[2]:
                lig += [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H" and a.altloc in ("\0", "A")]
    lig = np.array(lig); cen = lig.mean(0)
    mst = gemmi.read_structure(r.model_path if "model_path" in r else f"/path/to/afdb_dimers/afdb_human_dimers/homodimers/{r.model}.pdb")
    mca = read_model_ca(f"/path/to/afdb_dimers/afdb_human_dimers/homodimers/{r.model}.pdb")
    a1, a2 = r.v_local_dimer_assign[0], r.v_local_dimer_assign[1]
    Pp, Q = [], []
    for pc, mc in [(r.chain1, a1), (r.chain2, a2)]:
        for ls, rn, xyz in chains[pc]:
            u = sm.get((pc.split("-")[0], ls))
            if u and u[1] in mca[mc]: Pp.append(xyz); Q.append(mca[mc][u[1]][1])
    Pp, Q = np.array(Pp), np.array(Q); sel = np.linalg.norm(Pp - cen, axis=1) <= LOCAL_R
    Rm, t = kabsch(Pp[sel], Q[sel]); L = lig @ Rm.T + t
    xyz = {ch.name: np.array([[a.pos.x, a.pos.y, a.pos.z] for res in ch for a in res if a.element.name != "H"]) for ch in mst[0]}
    allx = np.vstack(list(xyz.values()))
    d_all = cKDTree(allx).query(L)[0]
    per = {c: int((cKDTree(x).query(L)[0] < CLASH).sum()) for c, x in xyz.items()}
    print(f"{r.pdb_id} {r.ccd} n_lig_atoms={len(L)} assign={a1}{a2} bind1={r.n_bind_pair} | ligand atoms clashing (<{CLASH} A) with dimer: {int((d_all < CLASH).sum())}, per chain {per}; "
          f"min dist {d_all.min():.2f}; median nearest-atom dist {np.median(d_all):.2f}")
