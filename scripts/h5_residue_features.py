#!/usr/bin/env python
"""H5 step 1a: per-residue structural features for every human AFDB dimer model.

For each model (both chains): UniProt accession, position (AFDB full-length numbering = UniProt position),
residue name, pLDDT (CA B-factor), relative SASA in the dimer context (rsa_dimer) and in the isolated-chain
context (rsa_iso, same coordinates, partner removed), and an interface flag.

Definitions (H5 analysis plan, written before outcomes were looked at):
  INTERFACE_CUTOFF = 5.0 A: a residue is 'interface' if any of its heavy atoms is within 5.0 A of any heavy atom
    of the partner chain (task specification for H5; DESIGN.md v3 does not define a residue interface cutoff).
  SASA: freesasa (Lee-Richards, default resolution, ProtOr radii, probe 1.4 A); relative SASA = freesasa
    relativeTotal (reference = freesasa's Ala-X-Ala tripeptide max areas). Isolated-chain context via
    freesasa 'separate-chains'.
Status: first attempt, not independently verified.

Usage: .venvs/h5/bin/python scripts/h5_residue_features.py --workers 8
"""
import os
import argparse, sys, time, traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/cluster/scratch/gmagnusson/dimers")) / "afdb_human_dimers"
KINDS = {"homo": "homodimers", "hetero": "heterodimers"}
INTERFACE_CUTOFF = 5.0
OUT = WS / "data/part2/h5/residue_features"


def models():
    rows = []
    for kind, sub in KINDS.items():
        c = pd.read_csv(DIMER_ROOT / sub / "candidates.csv").drop_duplicates("modelEntityId")
        c["kind"] = kind
        c["path"] = [str(DIMER_ROOT / sub / f"{m}.pdb") for m in c.modelEntityId]
        rows.append(c)
    return pd.concat(rows, ignore_index=True)


def one(args):
    import freesasa, gemmi
    from scipy.spatial import cKDTree
    path, mid, acc = args
    freesasa.setVerbosity(freesasa.silent)
    st = gemmi.read_structure(path); st.remove_hydrogens()
    model = st[0]
    chains = [ch.name for ch in model]
    if len(chains) != 2:
        raise ValueError(f"expected 2 chains, got {chains}")
    # heavy-atom coordinates per chain
    coords, owner = {}, {}
    for ch in model:
        xyz, own = [], []
        for i, r in enumerate(ch):
            for a in r:
                xyz.append(a.pos.tolist()); own.append(i)
        coords[ch.name] = np.asarray(xyz); owner[ch.name] = np.asarray(own)
    # SASA
    s_d = freesasa.Structure(path)
    r_d = freesasa.calc(s_d).residueAreas()
    r_i = {}
    for s in freesasa.structureArray(path, {"separate-chains": True}):
        ra = freesasa.calc(s).residueAreas()
        r_i.update(ra)
    out = []
    for ci, ch in enumerate(model):
        partner = chains[1 - ci]
        tree = cKDTree(coords[partner])
        d, _ = tree.query(coords[ch.name], k=1, distance_upper_bound=INTERFACE_CUTOFF + 1e-6)
        n_res = len(ch)
        mind = np.full(n_res, np.inf)
        np.minimum.at(mind, owner[ch.name], d)
        for i, r in enumerate(ch):
            key = str(r.seqid.num) + r.seqid.icode.strip()
            ca = r.find_atom("CA", "*")
            ad, ai = r_d[ch.name].get(key), r_i[ch.name].get(key)
            out.append((mid, ch.name, acc[ci], r.seqid.num, r.name, ca.b_iso if ca else np.nan,
                        ad.relativeTotal if ad else np.nan, ai.relativeTotal if ai else np.nan,
                        ad.total if ad else np.nan, ai.total if ai else np.nan,
                        mind[i] if np.isfinite(mind[i]) else np.nan, bool(mind[i] <= INTERFACE_CUTOFF)))
    return out


def run_chunk(job):
    k, rows = job
    f = OUT / f"chunk_{k:04d}.parquet"
    if f.exists():
        return k, 0, []
    recs, fails = [], []
    for path, mid, a1, a2 in rows:
        try:
            recs.extend(one((path, mid, (a1, a2))))
        except Exception as e:
            fails.append((mid, repr(e)[:200]))
    df = pd.DataFrame(recs, columns=["modelEntityId", "chain", "accession", "position", "resname", "plddt",
                                     "rsa_dimer", "rsa_iso", "sasa_dimer", "sasa_iso", "min_dist_partner", "interface"])
    df.to_parquet(f, index=False)
    return k, len(rows), fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    m = models()
    m.drop_duplicates("modelEntityId").to_csv(OUT / "models.tsv", sep="\t", index=False)
    if a.limit:
        m = m.head(a.limit)
    rows = list(zip(m.path, m.modelEntityId, m.uniprotAccession_1, m.uniprotAccession_2))
    jobs = [(i // a.chunk, rows[i:i + a.chunk]) for i in range(0, len(rows), a.chunk)]
    t0 = time.time(); allf = []
    with ProcessPoolExecutor(a.workers) as ex:
        for k, n, fails in ex.map(run_chunk, jobs):
            allf += fails
            print(f"chunk {k} n={n} fails={len(fails)} t={time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(allf, columns=["modelEntityId", "error"]).to_csv(OUT / "failures.tsv", sep="\t", index=False)
    print("done", len(rows), "models; failures", len(allf))


if __name__ == "__main__":
    main()
