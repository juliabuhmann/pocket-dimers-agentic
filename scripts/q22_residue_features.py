#!/usr/bin/env python
"""Per-residue features of every human dimer model (homo + hetero) for Q2.2.

Per residue: chain, resnum (= UniProt position), aa (1-letter), pLDDT (B-factor of CA),
rel_sasa_iso (FreeSASA relative total SASA with the chain alone, same coordinates),
rel_sasa_cpx (in the dimer), is_surface (rel_sasa_iso >= REL_SASA_SURFACE),
is_interface (any heavy atom within CONTACT_A of a partner-chain heavy atom).
Output: data/interim/q22/residues.parquet (one row per residue) + failures ledger.
Status: first attempt, not yet verified.
"""
import argparse, sys
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import DIMER_ROOT, KINDS, INTERIM, CONTACT_A, REL_SASA_SURFACE, load_annot

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
             "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}


def one(args):
    kind, stem = args
    import freesasa, gemmi
    from scipy.spatial import cKDTree
    freesasa.setVerbosity(freesasa.silent)
    path = str(DIMER_ROOT / KINDS[kind] / f"{stem}.pdb")
    try:
        st = gemmi.read_structure(path)
        st.remove_hydrogens(); st.remove_ligands_and_waters()
        m = st[0]
        rows, xyz, owner = [], {}, {}
        for ch in m:
            pts, own = [], []
            for r in ch:
                ca = r.find_atom("CA", "*")
                rows.append((ch.name, r.seqid.num, THREE2ONE.get(r.name, "X"), ca.b_iso if ca else np.nan))
                for a in r:
                    pts.append(a.pos.tolist()); own.append(r.seqid.num)
            xyz[ch.name] = np.array(pts); owner[ch.name] = np.array(own)
        res = pd.DataFrame(rows, columns=["chain", "resnum", "aa", "plddt"])
        chains = list(xyz)
        iface = set()
        if len(chains) == 2:
            a, b = chains
            ta, tb = cKDTree(xyz[a]), cKDTree(xyz[b])
            pairs = ta.query_ball_tree(tb, CONTACT_A)
            for i, js in enumerate(pairs):
                if js:
                    iface.add((a, int(owner[a][i])))
                    for j in js:
                        iface.add((b, int(owner[b][j])))
        # SASA: complex and each chain alone (same coordinates)
        cpx_res = freesasa.calc(freesasa.Structure(path)).residueAreas()
        sep = freesasa.structureArray(path, {"separate-chains": True, "separate-models": False})
        iso_res = {}
        for s in sep:
            ra = freesasa.calc(s).residueAreas()
            for ch, d in ra.items():
                iso_res[ch] = d
        def get(dd, ch, n, attr):
            try:
                v = getattr(dd[ch][str(n)], attr)
                return v if np.isfinite(v) else np.nan
            except KeyError:
                return np.nan
        res["rel_sasa_iso"] = [get(iso_res, c, n, "relativeTotal") for c, n in zip(res.chain, res.resnum)]
        res["rel_sasa_cpx"] = [get(cpx_res, c, n, "relativeTotal") for c, n in zip(res.chain, res.resnum)]
        res["sasa_iso"] = [get(iso_res, c, n, "total") for c, n in zip(res.chain, res.resnum)]
        res["sasa_cpx"] = [get(cpx_res, c, n, "total") for c, n in zip(res.chain, res.resnum)]
        res["is_surface"] = res.rel_sasa_iso >= REL_SASA_SURFACE
        res["is_interface"] = [(c, n) in iface for c, n in zip(res.chain, res.resnum)]
        res["modelEntityId"] = stem; res["kind"] = kind
        return res, None
    except Exception as e:  # ledger, never silent
        return None, (kind, stem, repr(e)[:300])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    d = load_annot()
    models = d[["kind", "modelEntityId"]].drop_duplicates().sort_values(["kind", "modelEntityId"])
    jobs = list(models.itertuples(index=False, name=None))
    if a.limit:
        jobs = jobs[:a.limit]
    frames, fails = [], []
    with Pool(a.workers) as p:
        for k, (res, err) in enumerate(p.imap_unordered(one, jobs, chunksize=4)):
            if res is not None:
                frames.append(res)
            else:
                fails.append(err)
            if k % 500 == 0:
                print(k, len(jobs), flush=True)
    out = pd.concat(frames, ignore_index=True)
    sfx = f"_limit{a.limit}" if a.limit else ""
    out.to_parquet(INTERIM / f"residues{sfx}.parquet")
    pd.DataFrame(fails, columns=["kind", "modelEntityId", "error"]).to_csv(INTERIM / f"residues_failures{sfx}.tsv", sep="\t", index=False)
    print(f"models={len(jobs)} ok={len(frames)} failed={len(fails)} residues={len(out)}")
    print(out.groupby("kind")[["is_surface", "is_interface"]].mean())


if __name__ == "__main__":
    main()
