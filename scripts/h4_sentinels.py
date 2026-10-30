#!/usr/bin/env python
"""H4 known-answer sentinels (step 6). First attempt, NOT verified beyond what is printed.

1. Hand check, independent re-derivation: SENTINEL_PDB / SENTINEL_CCD (glutathione in human GSTP1, a C2 homodimer
   whose G-site is lined by both subunits; chosen from the benchmark). Recomputed WITHOUT the pipeline code path:
   gemmi reads the asymmetric unit (both chains present in the AU), the ligand centroid is taken from the mmCIF
   atoms, chains are superposed GLOBALLY (all mapped CA, both chains; the pipeline uses a local 15 Å frame) onto
   the dimer model through SIFTS XML UniProt numbering, and the DCC to the dimer pockets is recomputed.
   Pass: |DCC_independent - DCC_pipeline| <= HAND_TOL for the best chain assignment.
2. Negative control: ligand centroid randomly rotated about the model centre -> hit rate must not exceed the
   random-patch null by more than NEG_MARGIN (absolute), dimer and isolated, on the primary benchmark.
Writes results/part2/h4/<census>/sentinels.tsv; exits non-zero on failure.
"""
import argparse
import sys

import numpy as np
import pandas as pd

from h4_common import OUT_DATA, OUT_RES, PDB_MIRROR, DIMER_ROOT, POCKETS_DIMER, kabsch, parse_res, write_tsv
from fetch_sifts_xml import fetch_sifts_xml, parse_sifts_xml, CACHE

SENTINEL_PDB, SENTINEL_CCD, SENTINEL_ACC = "11gs", "GSH", "P09211"   # GSTP1-1 with glutathione (Oakley et al. 1997)
HAND_TOL = 1.5     # Å; local vs global superposition may differ by up to ~1 Å for a rigid homodimer
NEG_MARGIN = 0.05  # absolute hit-rate margin over the random-patch null


def hand_check(census):
    import gemmi
    rec = pd.read_parquet(OUT_DATA / census / "recovery.parquet")
    r = rec[(rec.pdb_id == SENTINEL_PDB) & (rec.ccd == SENTINEL_CCD) & (rec.acc1 == SENTINEL_ACC)]
    if r.empty:
        return dict(check="hand_check", status="NOT_IN_BENCHMARK")
    r = r.sort_values("ligand_uid").iloc[0]
    fetch_sifts_xml([SENTINEL_PDB])
    sx = parse_sifts_xml(CACHE / f"{SENTINEL_PDB}.xml.gz")
    sx = sx[sx.uniprot_acc == SENTINEL_ACC]
    amap = {(c, int(n)): int(p) for c, n, p in zip(sx.chain, sx.auth_resnum, sx.uniprot_pos) if pd.notna(n)}
    st = gemmi.read_structure(str(PDB_MIRROR / SENTINEL_PDB[1:3] / f"{SENTINEL_PDB}.cif.gz"))
    st.setup_entities()
    m = st[0]
    ca, lig = {}, []
    label_asym, auth_seq = r.ligand_uid.split(":")[1], r.ligand_uid.split(":")[2]  # uid = pdb:label_asym:auth_seq:ccd:op
    for ch in m:
        for res in ch:
            if res.name == SENTINEL_CCD and res.subchain == label_asym and str(res.seqid.num) == auth_seq:
                lig += [a.pos.tolist() for a in res if a.element.name != "H" and a.altloc in ("\0", "A")]
            a = res.find_atom("CA", "*")
            if a is not None and (ch.name, res.seqid.num) in amap and res.het_flag == "A":
                ca[(ch.name, amap[(ch.name, res.seqid.num)])] = np.array(a.pos.tolist())
    if not lig:  # ligand may sit on another auth chain; fall back to any copy with this CCD nearest the census centroid
        return dict(check="hand_check", status="LIGAND_NOT_FOUND")
    cen = np.mean(lig, 0)
    mdl = gemmi.read_structure(str(DIMER_ROOT / "homodimers" / f"{r.model}.pdb"))[0]
    mca = {(ch.name, res.seqid.num): np.array(res.find_atom("CA", "*").pos.tolist()) for ch in mdl for res in ch
           if res.find_atom("CA", "*") is not None}
    pk = pd.read_parquet(POCKETS_DIMER)
    pk = pk[pk.structure_stem == r.model].sort_values("burrow_score", ascending=False)
    ctr = pk[["center_x", "center_y", "center_z"]].to_numpy()
    c1, c2 = r.chain1_auth, r.chain2_auth
    best = None
    for m1, m2 in (("A", "B"), ("B", "A")):
        P, Q = [], []
        for pc, mc in ((c1, m1), (c2, m2)):
            for (c, p), x in ca.items():
                if c == pc and (mc, p) in mca:
                    P.append(x); Q.append(mca[(mc, p)])
        R, t, rmsd = kabsch(np.array(P), np.array(Q))
        x = cen @ R.T + t
        d = np.linalg.norm(ctr - x, axis=1)
        if best is None or d.min() < best[0]:
            best = (d.min(), int(np.argmin(d)) + 1, rmsd, len(P))
    k = int(r.n_sites_pair) + 2
    ok = abs(best[0] - r.dimer_nocollapse_dcc) <= HAND_TOL or abs(best[0] - r.dimer_dcc) <= HAND_TOL
    return dict(check="hand_check", status="PASS" if ok else "FAIL", pdb=SENTINEL_PDB, ccd=SENTINEL_CCD, model=r.model,
                is_mbs=bool(r.is_mbs), is_on_axis=r.is_on_axis, iface_jaccard=r.iface_jaccard,
                pipeline_dimer_dcc=r.dimer_dcc, pipeline_dimer_nocollapse_dcc=r.dimer_nocollapse_dcc, pipeline_iso_dcc=r.iso_dcc,
                pipeline_dimer_hit=bool(r.dimer_hit), pipeline_iso_hit=bool(r.iso_hit), pipeline_rmsd_local=r.rmsd_sup,
                independent_dcc_any_pocket=best[0], independent_nearest_pocket_rank=best[1], independent_global_rmsd=best[2],
                independent_n_ca=best[3], topk=k)


def negative_control(census):
    b = pd.read_parquet(OUT_DATA / census / "benchmark_filtered.parquet")
    b = b[b.lig_set == "organic"]
    out = []
    for c in ("dimer", "iso"):
        sh = b[f"{c}_shuf_hit"].mean(); nl = b[f"{c}_null0_hit"].mean(); real = b[f"{c}_hit"].mean()
        out.append(dict(check=f"negative_control_{c}", status="PASS" if sh <= nl + NEG_MARGIN else "FAIL", n=len(b),
                        shuffled_hit_rate=sh, null_seed0_hit_rate=nl, real_hit_rate=real))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    a = ap.parse_args()
    rows = [hand_check(a.census)] + negative_control(a.census)
    df = pd.DataFrame(rows)
    write_tsv(df, OUT_RES / a.census / "sentinels.tsv", {"census": a.census})
    pd.set_option("display.width", 250)
    print(df.T.to_string())
    sys.exit(1 if (df.status == "FAIL").any() else 0)
