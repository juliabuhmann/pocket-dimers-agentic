#!/usr/bin/env python
"""fpocket at the sampled PDB ligand sites (EXPLORATORY, not pre-registered; docs/FPOCKET_DRUGGABILITY.md).

Sites: the 11,983 described sites of the PDB-site descriptor run (docs/PDB_SITE_DESCRIPTORS.md), read from
data/interim/pdbsite/sites_crop30.parquet (read-only; same site_id / ligand instance / weights / labels).
Per site:
 1. read <PDB mirror>/<id[1:3]>/<id>.cif.gz with gemmi, first altloc only, no H; build assembly 1 (first listed,
    as the census and scripts/pdbsite_describe.py) - WHOLE assembly, no crop;
 2. re-find the census ligand instance (CCD component het residues nearest the census centroid; same rule as
    scripts/pdbsite_describe.py build());
 3. keep only polymer residues (gemmi entity_type Polymer, so modified residues such as MSE stay); ligands, waters
    and additives are removed. Chains are renamed to one character for PDB format;
 4. fpocket -f <pdb> (default parameters);
 5. pick the pocket with the most ligand heavy atoms within OVERLAP_A of any of its alpha-sphere centres
    (ties: better fpocket rank = lower pocket number). Record its fpocket descriptors and the overlap.
Outputs (resumable): data/interim/fpocket/pdbsite/part_*.parquet
Run: source scripts/env_threads.sh; .venvs/part1/bin/python scripts/fpocket_pdbsite.py --workers 10 [--limit 50]
Status: first attempt, not yet verified.
"""
import argparse, glob, os, sys, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fpocket_common import WS, INTERIM, OVERLAP_A, run_fpocket, fpocket_version

PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/path/to/pdb_mirror"))     # same mirror as src/part1/config.py and scripts/pdbsite_common.py
SITES = WS / "data/interim/pdbsite/sites_crop30.parquet"
OUTDIR = INTERIM / "pdbsite"
PDB_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
MAX_PDB_ATOMS = 99999   # PDB-format atom serial limit


def build(row):
    import gemmi
    pdb = row["pdb_id"].lower()
    st = gemmi.read_structure(str(PDB_MIRROR / pdb[1:3] / f"{pdb}.cif.gz"))
    st.remove_alternative_conformations(); st.remove_hydrogens(); st.setup_entities()
    model = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber) if len(st.assemblies) else st[0]
    comps = set(row["components"]); cen = np.asarray(row["centroid"], float)
    cand = []
    for ch in model:
        for r in ch:
            if r.name in comps and r.het_flag == "H":
                xyz = np.array([a.pos.tolist() for a in r if a.element.name not in ("H", "D")])
                if len(xyz):
                    cand.append((float(np.linalg.norm(xyz.mean(0) - cen)), xyz))
    if not cand:
        raise ValueError("ligand not found in assembly")
    cand.sort(key=lambda t: t[0])
    if len(row["components"]) == 1:
        lig = cand[0][1]
    else:
        parts, n = [], 0
        for d, xyz in cand:
            parts.append(xyz); n += len(xyz)
            if n >= int(row["n_heavy_atoms_instance"]):
                break
        lig = np.vstack(parts)
    st2 = gemmi.Structure(); st2.cell = st.cell; st2.spacegroup_hm = st.spacegroup_hm
    m2 = gemmi.Model("1")
    n_ch, n_atoms = 0, 0
    for ch in model:
        pol = [r for r in ch if r.entity_type == gemmi.EntityType.Polymer and r.name != "HOH"]
        if not pol:
            continue
        if n_ch >= len(PDB_CHAIN_IDS):
            raise ValueError(f"more than {len(PDB_CHAIN_IDS)} polymer chains (PDB format limit)")
        nc = gemmi.Chain(PDB_CHAIN_IDS[n_ch]); n_ch += 1
        for r in pol:
            nc.add_residue(r); n_atoms += len(r)
        m2.add_chain(nc)
    if n_atoms > MAX_PDB_ATOMS:
        raise ValueError(f"{n_atoms} atoms > PDB-format limit")
    st2.add_model(m2); st2.setup_entities()
    return st2, lig, float(np.linalg.norm(lig.mean(0) - cen)), n_ch, n_atoms


def one(row):
    from scipy.spatial import cKDTree
    t0 = time.perf_counter()
    rec = {"site_id": row["site_id"], "status": "ok", "error": ""}
    try:
        st2, lig, md, n_ch, n_atoms = build(row)
        rec.update(fp_lig_match_dist=md, fp_n_lig_atoms=len(lig), fp_n_chains=n_ch, fp_n_atoms=n_atoms)
        with tempfile.TemporaryDirectory(prefix="fpb_") as td:
            p = Path(td) / "x.pdb"
            st2.write_minimal_pdb(str(p))
            t1 = time.perf_counter()
            pk, cen, _ = run_fpocket(p)
            rec["fpocket_s"] = time.perf_counter() - t1
        rec["fp_n_pockets"] = len(pk)
        best, best_n = None, 0
        for k in sorted(cen):
            d, _ = cKDTree(cen[k]).query(lig, k=1)
            n = int((d <= OVERLAP_A).sum())
            if n > best_n:
                best, best_n = k, n
        rec["fp_overlap_n"] = best_n
        rec["fp_overlap_frac"] = best_n / len(lig)
        if best is None:
            rec["status"] = "no_overlap"
        else:
            rec.update(pk[best])
            rec["fp_rank"] = best   # fpocket numbers pockets by its ranking score
    except Exception as e:  # noqa: BLE001 (ledger, never silent)
        rec.update(status="error", error=f"{type(e).__name__}: {e}"[:300])
    rec["total_s"] = time.perf_counter() - t0
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0, help="first N sites of a seed-0 shuffle (pilot)")
    ap.add_argument("--subsample", default="", help="TSV of site_id to run (stratified subsample), else all")
    ap.add_argument("--flush", type=int, default=200)
    ap.add_argument("--out", default=str(OUTDIR))
    a = ap.parse_args()
    assert a.workers <= 10, "task budget: <=10 CPU workers"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    print(fpocket_version(), flush=True)
    s = pd.read_parquet(SITES, columns=["site_id", "pdb_id", "components", "centroid", "n_heavy_atoms_instance"])
    s = s.sample(frac=1, random_state=0)          # seed-0 order so a pilot / partial run is a random subset
    if a.subsample:
        keep = set(pd.read_csv(a.subsample, sep="\t", comment="#").site_id)
        s = s[s.site_id.isin(keep)]
    if a.limit:
        s = s.head(a.limit)
    done = set()
    for f in glob.glob(str(out / "part_*.parquet")):
        done |= set(pd.read_parquet(f, columns=["site_id"]).site_id)
    rows = s[~s.site_id.isin(done)].to_dict("records")
    print(f"sites {len(s)} done {len(done)} todo {len(rows)}", flush=True)
    k0 = len(glob.glob(str(out / "part_*.parquet")))
    buf, t0, n = [], time.time(), 0
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for rec in ex.map(one, rows, chunksize=2):
            buf.append(rec); n += 1
            if len(buf) >= a.flush:
                pd.DataFrame(buf).to_parquet(out / f"part_{k0:04d}.parquet", index=False); k0 += 1; buf = []
                el = time.time() - t0
                print(f"{n}/{len(rows)} {el:.0f}s {el / n:.3f} s/site", flush=True)
    if buf:
        pd.DataFrame(buf).to_parquet(out / f"part_{k0:04d}.parquet", index=False)
    el = time.time() - t0
    print(f"finished {n} sites in {el:.0f}s ({el / max(n, 1):.3f} s/site aggregate)", flush=True)


if __name__ == "__main__":
    main()
