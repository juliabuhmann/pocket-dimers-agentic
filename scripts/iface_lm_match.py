#!/usr/bin/env python
"""EXPLORATORY: label each dimer-model pocket by whether it matches a PDB ligand site (H4 benchmark rows).

Match (per pocket x site x chain assignment; homodimers try both A/B assignments, as h4_recovery.py):
  DCC: distance(pocket centre, ligand centroid superposed with the h4 local-frame rule) <= DCC_HIT (h4_common)
  OC : |L & P| / min(|L|, |P|) >= OC_MIN, L = PDB binding residues mapped to UniProt (= model numbering), P = pocket lining.
Uses ALL pockets of the model (not top-k, no collapse); symmetric-representative selection is done downstream.
Flags per pocket: matched to a primary-set (organic, interface-filtered) MBS / SBS site, and to ANY recovery row
(any ligand class, no interface filter, MBS or SBS), the latter used to define "unmatched" conservatively.
Output: data/interim/iface_ligand_matched/pocket_match.parquet. Status: first attempt, not verified.
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from h4_common import OUT_DATA, DCC_HIT, SUPERPOSE_RADIUS, SUPERPOSE_MIN_CA, DIMER_ROOT, KINDS, kabsch, parse_res
from h4_recovery import read_model
from q22_common import INTERIM

OC_MIN = 0.5  # task spec
OUTD = INTERIM.parent / "iface_ligand_matched"
FLAGS = ["mbs_prim", "sbs_prim", "any", "mbs_prim_dcc", "mbs_prim_oc", "sbs_prim_dcc", "sbs_prim_oc"]


def work(args):
    model, kind, P, sites, cas = args
    pr = [parse_res(s) for s in P.residues]
    ctr = P[["center_x", "center_y", "center_z"]].to_numpy(float)
    fl = {k: np.zeros(len(P), bool) for k in FLAGS}
    ca = read_model(DIMER_ROOT / KINDS[kind] / f"{model}.pdb")[0]
    for s in sites:
        assigns = [("A", "B"), ("B", "A")] if kind == "homo" else ([("A", "B")] if s["acc1"] == s["model_acc_A"] else [("B", "A")])
        cen = np.asarray(s["centroid"], float)
        for m1, m2 in assigns:
            A, B, near = [], [], []
            for role, mch in (("1", m1), ("2", m2)):
                key = (s["pdb_id"], s[f"chain{role}"], s[f"acc{role}"])
                if key not in cas:
                    continue
                pos, xyz = cas[key]
                for p, x in zip(pos, xyz):
                    y = ca.get((mch, p))
                    if y is not None:
                        A.append(x); B.append(y); near.append(np.linalg.norm(x - cen) <= SUPERPOSE_RADIUS)
            if len(A) < 3:
                continue
            A, B, near = np.array(A), np.array(B), np.array(near, bool)
            sel = near if near.sum() >= SUPERPOSE_MIN_CA else np.ones(len(A), bool)
            R, t, _ = kabsch(A[sel], B[sel])
            x = cen @ R.T + t
            dcc = np.linalg.norm(ctr - x, axis=1) <= DCC_HIT
            L = {(m1, p) for p in s["bind1"]} | {(m2, p) for p in s["bind2"]}
            oc = np.array([len(L & p) / min(len(L), len(p)) >= OC_MIN if L and p else False for p in pr])
            hit = dcc | oc
            fl["any"] |= hit
            if s["primary"]:
                k = "mbs_prim" if s["is_mbs"] else "sbs_prim"
                fl[k] |= hit; fl[k + "_dcc"] |= dcc; fl[k + "_oc"] |= oc
    out = P[["pocket_id"]].copy()
    for k, v in fl.items():
        out["match_" + k] = v
    return out


def main():
    rec = pd.read_parquet(OUT_DATA / "v2/recovery.parquet")
    bf = pd.read_parquet(OUT_DATA / "v2/benchmark_filtered.parquet", columns=["row_id", "lig_set"])
    prim = set(bf.row_id[bf.lig_set == "organic"])
    rec = rec[rec.superpose.notna()].copy()
    rec["primary"] = rec.row_id.isin(prim)
    print("site rows", len(rec), "primary", int(rec.primary.sum()), "models", rec.model.nunique(), flush=True)
    ca = pd.read_parquet(OUT_DATA / "v2/chain_ca.parquet")
    cas = {(p, c, a): (np.asarray(pos, int), np.asarray(x, float).reshape(-1, 3)) for p, c, a, pos, x in
           zip(ca.pdb_id, ca.chain, ca.acc, ca.pos, ca.xyz)}
    feats = pd.read_parquet(INTERIM / "dimer_pocket_features_units.parquet")
    feats = feats[feats.modelEntityId.isin(set(rec.model))]
    cols = ["row_id", "pdb_id", "chain1", "chain2", "acc1", "acc2", "model_acc_A", "bind1", "bind2", "bind1", "centroid", "is_mbs", "primary"]
    cols = list(dict.fromkeys(cols))
    pg = dict(list(feats.groupby("modelEntityId")))
    tasks = []
    for (m, kind), g in rec.groupby(["model", "kind"]):
        if m not in pg:
            continue
        sites = g[cols].to_dict("records")
        pdbs = set(g.pdb_id)
        tasks.append((m, kind, pg[m], sites, {k: v for k, v in cas.items() if k[0] in pdbs}))
    print("models with pockets", len(tasks), flush=True)
    res = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for i, r in enumerate(ex.map(work, tasks, chunksize=4)):
            res.append(r)
            if (i + 1) % 200 == 0:
                print(i + 1, flush=True)
    m = pd.concat(res)
    OUTD.mkdir(parents=True, exist_ok=True)
    m.to_parquet(OUTD / "pocket_match.parquet")
    print(m.drop(columns="pocket_id").sum())


if __name__ == "__main__":
    main()
