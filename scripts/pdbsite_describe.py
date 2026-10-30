#!/usr/bin/env python
"""BURROW pocket at each sampled PDB ligand site (pocketcompare describe, locate=burrow, CPU fast backend).

Per site (data/interim/pdbsite/sample_*.parquet from scripts/pdbsite_sample.py):
 1. read <PDB mirror>/<id[1:3]>/<id>.cif.gz with gemmi, first altloc only, no hydrogens; build the first listed
    assembly (as the census, src/part1/config.py ASSEMBLY_ID=None) with gemmi.make_assembly;
 2. find the census ligand instance in the assembly: het residues named as its CCD components, nearest to the census
    centroid (multi-component instances: nearest residues until the census heavy-atom count is reached); flag if the
    rebuilt centroid is > LIG_MATCH_A from the census centroid;
 3. keep all polymer chains of the assembly (primary, DESIGN: whole assembly 1); with --crop X (sensitivity only)
    keep only polymer chains with a heavy atom <= X Å of the ligand. Remove all ligands and waters, write a temp mmCIF;
 4. pocketcompare describe with the ligand heavy atoms as the site and --locate burrow (same code path as
    `pocketcompare describe x.cif --ligand lig.sdf --locate burrow`, ligand passed as coordinates in-process);
 5. lining chemistry with the residue classes of scripts/q22_common.py; buriedness = mean(1 - min(rel SASA, 1)) of
    lining residues in the described complex (whole assembly, or the cropped one) (FreeSASA, as scripts/q22_residue_features.py); lining residues per chain.
Outputs (resumable): data/interim/pdbsite/describe/part_*.parquet (one row per site, fingerprint as list<float32>)
Run with the pockets venv (source scripts/env_tools.sh pockets). FreeSASA is imported from the q22 venv (same Python 3.12).
Status: first attempt, not yet verified.
"""
import argparse, glob, os, sys, time, traceback
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdbsite_common import (WS, PDB_MIRROR, INTERIM, LIG_MATCH_A, INTERFACE_MIN_RES, HYDROPHOBIC, AROMATIC,
                            CHARGED)

Q22_SITE = sorted(glob.glob(str(WS / ".venvs/q22/lib/python3*/site-packages")))
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
             "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V", "MSE": "M"}   # MSE -> M (selenomethionine; our choice)
SASA_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
TMP = INTERIM / "tmp"
OUTDIR = INTERIM / "describe_uncropped"   # primary; the 30 Å crop sensitivity lives in describe_crop30


def init():
    from pocketcompare import cli
    cli._worker_init(single_thread=True)
    if Q22_SITE:
        sys.path.append(Q22_SITE[-1])


def resname_key(r):
    return f"{r.seqid.num}{r.seqid.icode.strip()}"


def build(row):
    import gemmi
    pdb = row["pdb_id"].lower()
    st = gemmi.read_structure(str(PDB_MIRROR / pdb[1:3] / f"{pdb}.cif.gz"))
    st.remove_alternative_conformations(); st.remove_hydrogens()
    st.setup_entities()
    if len(st.assemblies):
        model = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    else:
        model = st[0]
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
    match_d = float(np.linalg.norm(lig.mean(0) - cen))
    # polymer chains to keep (all, or within --crop Å of the ligand for the sensitivity run)
    from scipy.spatial import cKDTree
    lt = cKDTree(lig)
    st2 = gemmi.Structure(); st2.cell = st.cell; st2.spacegroup_hm = st.spacegroup_hm
    m2 = gemmi.Model(1)
    n_chains_asm, kept = 0, 0
    for ch in model:
        pol = [r for r in ch if r.het_flag != "H" and r.name != "HOH"]
        if not pol:
            continue
        n_chains_asm += 1
        xyz = np.array([a.pos.tolist() for r in pol for a in r])
        crop = row.get("crop") or 0
        if crop <= 0 or lt.query(xyz, k=1, distance_upper_bound=crop)[0].min() <= crop:
            nc = gemmi.Chain(ch.name)
            for r in pol:
                nc.add_residue(r)
            m2.add_chain(nc); kept += 1
    st2.add_model(m2)
    st2.setup_entities()
    return st2, lig, match_d, n_chains_asm, kept


def sasa_bur(st2, lining):
    """{(chain, resnum_str): 1 - min(relSASA,1)} for lining residues (FreeSASA on a temp PDB, chains renamed)."""
    import freesasa
    freesasa.setVerbosity(freesasa.silent)
    names = [c.name for c in st2[0]]
    if len(names) > len(SASA_CHAIN_IDS):
        return {}
    ren = dict(zip(names, SASA_CHAIN_IDS))
    s3 = st2.clone()
    for c in s3[0]:
        c.name = ren[c.name]
    s3.setup_entities()
    p = TMP / f"sasa_{os.getpid()}.pdb"
    s3.write_pdb(str(p))
    ra = freesasa.calc(freesasa.Structure(str(p))).residueAreas()
    out = {}
    for tok in lining:
        ch, _, num = tok.rpartition(":")
        try:
            v = ra[ren[ch]][num].relativeTotal
            out[tok] = 1 - min(v, 1.0) if np.isfinite(v) else np.nan
        except KeyError:
            out[tok] = np.nan
    p.unlink(missing_ok=True)
    return out


def lining_feats(lining, aa, pre):
    codes = [aa.get(t, "X") for t in lining]
    n = len(codes)
    f = {"n_lining": n, "n_lining_std": sum(c != "X" for c in codes),
         "frac_hydrophobic": float(np.mean([c in HYDROPHOBIC for c in codes])) if n else np.nan,
         "frac_aromatic": float(np.mean([c in AROMATIC for c in codes])) if n else np.nan,
         "frac_charged": float(np.mean([c in CHARGED for c in codes])) if n else np.nan}
    per = pd.Series([t.rpartition(":")[0] for t in lining], dtype=object).value_counts()
    f["lining_chain_counts"] = ";".join(f"{k}={v}" for k, v in per.items())
    f["n_lining_chains"] = int(len(per))
    f["frac_top_chain"] = float(per.iloc[0] / n) if n else np.nan
    f["frac_second_chain"] = float(per.iloc[1] / n) if len(per) > 1 else 0.0
    f["pocket_interface"] = bool((per >= INTERFACE_MIN_RES).sum() >= 2)
    return {pre + k: v for k, v in f.items()}


def one(row):
    from pocketcompare import cli
    t0 = time.perf_counter()
    rec = {"site_id": row["site_id"], "status": "ok", "error": ""}
    fp = None
    path = TMP / f"{row['site_id']}.cif"
    try:
        st2, lig, match_d, n_asm, kept = build(row)
        rec.update(lig_match_dist=match_d, lig_match_ok=match_d <= LIG_MATCH_A, n_lig_atoms_rebuilt=len(lig),
                   n_chains_asm=n_asm, n_chains_kept=kept)
        st2.make_mmcif_document().write_file(str(path))
        t1 = time.perf_counter()
        opt = {"af_mode": False, "chains": None, "backend": "fast", "site": lig, "locate": "burrow",
               "unit": "merged", "min_points": 10, "n1": False}
        fs, _, tim = cli.featurize(Path(path), opt)
        info, V, _ = cli.describe_featurized(fs, path, None, opt, tim)
        rec["describe_s"] = time.perf_counter() - t1
        r = info.iloc[0]
        for k in ("burrow_score", "burrow_rank", "volume_A3", "site_volume_A3", "n_points", "site_n_points",
                  "sph8_fallback", "burrow_index", "coverage", "merged_coverage", "merged_precision", "d_match", "n_site_points"):
            if k in info:
                v = r[k]
                rec[k] = v.item() if hasattr(v, "item") else v
        aa = {}
        for c in st2[0]:
            for res in c:
                aa[f"{c.name}:{resname_key(res)}"] = THREE2ONE.get(res.name, "X")
        lining = r["residues"].split()
        rec["residues"] = " ".join(lining)
        rec.update(lining_feats(lining, aa, ""))
        # merged pocket containing the located member (= the unit of `describe` without a site and of the AFDB
        # atlas / Q2.2 dimer pockets): score = max member score, volume and lining of the merged point set
        from pocketcompare import pockets as pk
        bi = int(r["burrow_index"])
        mp = next((m for m in pk.merge_pockets(fs, "region") if bi in [int(x) for x in m.members]), None)
        mlin = []
        if mp is not None:
            rec["merged_burrow_score"] = float(mp.score)
            rec["merged_volume_A3"] = float(fs.volume(np.asarray(mp.idx, np.int64)))
            rec["merged_n_points"] = int(len(mp.idx))
            rec["merged_n_members"] = int(len(mp.members))
            mlin = pk.sort_residues(list(mp.residues))
            rec["merged_residues"] = " ".join(mlin)
            rec.update(lining_feats(mlin, aa, "merged_"))
        try:
            b = sasa_bur(st2, sorted(set(lining) | set(mlin)))
            for pre, lst in (("", lining), ("merged_", mlin)):
                vals = [b[t] for t in lst if t in b and np.isfinite(b[t])]
                rec[pre + "buriedness"] = float(np.mean(vals)) if vals else np.nan
        except Exception as e:  # noqa: BLE001
            rec["buriedness"] = np.nan; rec["sasa_error"] = repr(e)[:200]
        fp = np.asarray(V[0], np.float32)
    except Exception as e:  # noqa: BLE001  (ledger, never silent)
        rec.update(status="error", error=f"{type(e).__name__}: {e}"[:300])
    finally:
        path.unlink(missing_ok=True)
    rec["total_s"] = time.perf_counter() - t0
    rec["fingerprint"] = fp.tolist() if fp is not None else None
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--flush", type=int, default=250)
    ap.add_argument("--out", default=str(OUTDIR))
    ap.add_argument("--crop", type=float, default=0.0, help="sensitivity: crop radius in Å (0 = whole assembly, primary)")
    a = ap.parse_args()
    assert a.workers <= 12, "compute policy for this task: <=12 CPU workers"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); TMP.mkdir(parents=True, exist_ok=True)
    smp = pd.read_parquet(a.sample)
    if a.limit:
        smp = smp.head(a.limit)
    done = set()
    for f in glob.glob(str(out / "part_*.parquet")):
        done |= set(pd.read_parquet(f, columns=["site_id"]).site_id)
    todo = smp[~smp.site_id.isin(done)]
    print(f"sample {len(smp)} done {len(done)} todo {len(todo)}", flush=True)
    cols = ["site_id", "pdb_id", "components", "centroid", "n_heavy_atoms_instance"]
    rows = todo[cols].assign(crop=a.crop).to_dict("records")
    k0 = len(glob.glob(str(out / "part_*.parquet")))
    buf, t0, n = [], time.time(), 0
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(a.workers, initializer=init) as pool:
        for rec in pool.imap_unordered(one, rows, chunksize=2):
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
