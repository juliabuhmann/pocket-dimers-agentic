#!/usr/bin/env python3
"""Metrics for deep dive #1 AF3 runs (results/cofold/onaxis/cases.tsv) against the PDB assembly-1 pose.

Status: first attempt, UNVERIFIED.

Rerun:  .venvs/part1/bin/python scripts/cofold_metrics_onaxis.py   -> results/cofold/onaxis/{metrics,metrics_samples}.tsv

Per prediction sample (5 per run) and for the top-ranked sample:
- lig_rmsd: heavy-atom RMSD of the predicted ligand vs the census reference instance (ligand_uid; picked in the
  assembly by centroid), after superposing protein CA (dimer: both chain mappings tried; monomer: onto either PDB
  chain), atoms matched by CCD atom name (AF3 CCD input keeps CCD names). Minimum over chain mappings and predicted
  ligand copies. lig_rmsd_hung = same but atoms matched by a Hungarian assignment within element (lower bound on a
  symmetry-corrected RMSD). Success = lig_rmsd <= SUCCESS_RMSD.
- n_res_chainA/B: residues per protein chain with a heavy atom within CONTACT_CUTOFF of the predicted ligand copy
  that matched the reference; contacts_both = both >= 1.
- c2_axis_dist_pred: ligand centroid distance to the predicted A->B superposition rotation axis (+ rotation angle).
- ca_rmsd, irmsd: global CA RMSD and interface backbone RMSD (residues with any heavy atom within IFACE_CUTOFF of
  the other chain in the PDB; DockQ-like iRMSD, DockQ itself not installed).
- site_dist: distance between predicted ligand centroid and the centroid (predicted coords) of CA atoms of the PDB
  site residues (residues within CONTACT_CUTOFF of the reference ligand; monomer: residues of the matched chain);
  site_frac: fraction of those site residues within CONTACT_CUTOFF of the predicted ligand.
- AF3 confidences: ranking_score, iptm, ptm, ligand chain_iptm, max ligand-protein chain_pair_iptm, has_clash.
- lig_rmsd_pocket: same, but superposed on CA atoms of reference residues within POCKET_CUTOFF of the reference
  ligand (all mapped chains), i.e. pose accuracy independent of the global inter-chain arrangement.
- finite: all predicted coordinates finite.
"""
import os
import json
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data/cofold/onaxis"
OUT = ROOT / "results/cofold/onaxis"
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/cluster/scratch/jjaenes/pdb"))  # = src/part1/config.py PDB_MIRROR
SUCCESS_RMSD = 2.0     # task definition (A)
CONTACT_CUTOFF = 4.0   # A, heavy-atom contact for "ligand contacts chain" / site residues
POCKET_CUTOFF = 10.0   # A, residues (any heavy atom) around the reference ligand used for the pocket superposition
IFACE_CUTOFF = 10.0    # A, CAPRI/DockQ interface definition for iRMSD
BB = ("N", "CA", "C", "O")


def kabsch(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    H = (X - cx).T @ (Y - cy)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cy - R @ cx


def rmsd(X, Y):
    return float(np.sqrt(((X - Y) ** 2).sum(1).mean()))


def heavy(res, alt=("\0", "A")):
    return {a.name: (a.element.name, np.array(a.pos.tolist())) for a in res if a.element.name != "H" and a.altloc in alt}


_REF = {}


def ref_structure(case):
    if case.case_id not in _REF:
        _REF[case.case_id] = _ref_structure(case)
    return _REF[case.case_id]


def _ref_structure(case):
    st = gemmi.read_structure(str(PDB_MIRROR / case.pdb_id[1:3] / f"{case.pdb_id}.cif.gz"))
    st.setup_entities()
    a = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    a.remove_hydrogens()
    chains = case.pchains.split(",")
    prot = {}
    for cn in chains:
        ch = a[cn]
        prot[cn] = {r.label_seq - 1 - case.seqres_lo: heavy(r) for r in ch.get_polymer() if r.label_seq is not None}
    ligs = []
    for ch in a:
        for r in ch:
            if r.name == case.ccd and r.het_flag == "H":
                h = heavy(r)
                if h:
                    ligs.append(h)
    cen = np.array(json.loads(case.centroid))
    ligs.sort(key=lambda h: np.linalg.norm(np.mean([v[1] for v in h.values()], 0) - cen))
    # dedupe altloc/duplicate copies by centroid
    uniq = []
    for h in ligs:
        c = np.mean([v[1] for v in h.values()], 0)
        if all(np.linalg.norm(c - np.mean([v[1] for v in u.values()], 0)) > 2.0 for u in uniq):
            uniq.append(h)
    return chains, prot, uniq


def pred_structure(path, ccd):
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    prot, ligs, allxyz = {}, [], []
    for ch in st[0]:
        for r in ch:
            allxyz += [a.pos.tolist() for a in r]
        if ch[0].name == ccd and ch[0].het_flag == "H":
            ligs.append(heavy(ch[0], alt=("\0",)))
        else:
            prot[ch.name] = {r.seqid.num - 1: heavy(r, alt=("\0",)) for r in ch}
    return prot, ligs, bool(np.isfinite(np.array(allxyz)).all())


def ca_pairs(P, R, atoms=("CA",), keys=None):
    X, Y = [], []
    for k in (keys if keys is not None else sorted(set(P) & set(R))):
        if k not in P or k not in R:
            continue
        for an in atoms:
            if an in P[k] and an in R[k]:
                X.append(P[k][an][1]); Y.append(R[k][an][1])
    return np.array(X), np.array(Y)


def lig_rmsds(pl, rl):
    """name-matched and Hungarian-by-element RMSD between two ligand atom dicts (already in the same frame)."""
    common = [n for n in rl if n in pl and pl[n][0] == rl[n][0]]
    rn = rmsd(np.array([pl[n][1] for n in common]), np.array([rl[n][1] for n in common])) if common else np.nan
    sq = []
    for el in {v[0] for v in rl.values()}:
        A = np.array([v[1] for v in pl.values() if v[0] == el])
        B = np.array([v[1] for v in rl.values() if v[0] == el])
        if len(A) == 0 or len(B) == 0:
            continue
        D = ((A[:, None] - B[None]) ** 2).sum(-1)
        i, j = linear_sum_assignment(D)
        sq += list(D[i, j])
    rh = float(np.sqrt(np.mean(sq))) if sq else np.nan
    return rn, rh, len(common), len(rl)


def transform(d, R, t):
    return {k: {n: (e, R @ x + t) for n, (e, x) in v.items()} for k, v in d.items()}


def residues_near(prot_chain, xyz, cut):
    keys = [k for k, atoms in prot_chain.items() for _ in atoms]
    if not keys:
        return set()
    P = np.array([x for atoms in prot_chain.values() for _, x in atoms.values()])
    from scipy.spatial import cKDTree
    hit = cKDTree(np.asarray(xyz)).query(P, distance_upper_bound=cut)[0] <= cut
    return {k for k, h in zip(keys, hit) if h}


def axis_dist(PA, PB, cen):
    X, Y = ca_pairs(PA, PB)
    R, t = kabsch(X, Y)
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    w, v = np.linalg.eig(R)
    u = np.real(v[:, np.argmin(np.abs(w - 1))]); u /= np.linalg.norm(u)
    s = float(u @ t)
    p = np.linalg.lstsq(np.eye(3) - R, t - s * u, rcond=None)[0]
    dv = cen - p
    return ang, float(np.linalg.norm(dv - (dv @ u) * u))


def evaluate(case, cond, model_path, conf):
    rchains, rprot, rligs = ref_structure(case)
    pprot, pligs, finite = pred_structure(model_path, case.ccd)
    out = dict(finite=finite, n_pred_lig=len(pligs), n_ref_lig=len(rligs))
    if not finite or not pligs or not rligs:
        return out
    rlig = rligs[0]
    rlig_xyz = np.array([v[1] for v in rlig.values()])
    pch = sorted(pprot)
    if cond.startswith("dimer"):
        maps = [dict(zip(pch, rchains)), dict(zip(pch, rchains[::-1]))]
    else:
        maps = [{pch[0]: rc} for rc in rchains]
    best = None
    for m in maps:
        X = np.vstack([ca_pairs(pprot[p], rprot[r])[0] for p, r in m.items()])
        Y = np.vstack([ca_pairs(pprot[p], rprot[r])[1] for p, r in m.items()])
        R, t = kabsch(X, Y)
        ca = rmsd(X @ R.T + t, Y)
        for li, pl in enumerate(pligs):
            plt_ = {n: (e, R @ x + t) for n, (e, x) in pl.items()}
            rn, rh, nm, nr = lig_rmsds(plt_, rlig)
            if best is None or rn < best["lig_rmsd"]:
                best = dict(lig_rmsd=rn, lig_rmsd_hung=rh, n_atoms_matched=nm, n_atoms_ref=nr, ca_rmsd=ca,
                            _m=m, _li=li, _R=R, _t=t)
    m, li = best.pop("_m"), best.pop("_li")
    best.pop("_R"); best.pop("_t")
    out.update(best)
    # pocket-aligned ligand RMSD: superpose on CA of reference residues within POCKET_CUTOFF of the reference ligand
    # (separates pocket/pose accuracy from errors in the inter-chain arrangement); min over mappings and copies
    pk = None
    for mm in maps:
        X, Y = [], []
        for p, r in mm.items():
            keys = sorted(residues_near(rprot[r], rlig_xyz, POCKET_CUTOFF))
            x, y = ca_pairs(pprot[p], rprot[r], ("CA",), keys)
            if len(x):
                X.append(x); Y.append(y)
        if not X or sum(len(x) for x in X) < 3:
            continue
        X, Y = np.vstack(X), np.vstack(Y)
        R, t = kabsch(X, Y)
        for pl_ in pligs:
            rn = lig_rmsds({n: (e, R @ x + t) for n, (e, x) in pl_.items()}, rlig)[0]
            pk = rn if pk is None or rn < pk else pk
    out["lig_rmsd_pocket"] = pk
    pl = pligs[li]
    pxyz = np.array([v[1] for v in pl.values()])
    pcen = pxyz.mean(0)
    for i, p in enumerate(pch):
        out[f"n_res_chain{'AB'[i]}"] = len(residues_near(pprot[p], pxyz, CONTACT_CUTOFF))
    out["contacts_both"] = bool(cond.startswith("dimer") and out["n_res_chainA"] >= 1 and out.get("n_res_chainB", 0) >= 1)
    # true site residues (reference) on mapped chains -> predicted coordinates
    site_ca, site_hit, site_n = [], 0, 0
    for p, r in m.items():
        site = residues_near(rprot[r], rlig_xyz, CONTACT_CUTOFF)
        near = residues_near(pprot[p], pxyz, CONTACT_CUTOFF)
        for k in site:
            if k in pprot[p] and "CA" in pprot[p][k]:
                site_ca.append(pprot[p][k]["CA"][1]); site_n += 1; site_hit += k in near
    if site_ca:
        out["site_dist"] = float(np.linalg.norm(pcen - np.mean(site_ca, 0)))
        out["site_frac"] = site_hit / site_n
    out["n_site_res"] = site_n
    if cond.startswith("dimer"):
        out["pred_rot_angle"], out["c2_axis_dist_pred"] = axis_dist(pprot[pch[0]], pprot[pch[1]], pcen)
        # interface backbone RMSD (reference interface)
        ra, rb = m[pch[0]], m[pch[1]]
        rxyz = {c: np.array([x for v in rprot[c].values() for _, x in v.values()]) for c in (ra, rb)}
        ia = residues_near(rprot[ra], rxyz[rb], IFACE_CUTOFF)
        ib = residues_near(rprot[rb], rxyz[ra], IFACE_CUTOFF)
        X1, Y1 = ca_pairs(pprot[pch[0]], rprot[ra], BB, sorted(ia))
        X2, Y2 = ca_pairs(pprot[pch[1]], rprot[rb], BB, sorted(ib))
        X, Y = np.vstack([X1, X2]), np.vstack([Y1, Y2])
        R, t = kabsch(X, Y)
        out["irmsd"] = rmsd(X @ R.T + t, Y)
        out["n_iface_res"] = len(ia) + len(ib)
    # confidences
    ids = conf["chain_ids"]
    lig_idx = [i for i, c in enumerate(ids) if c not in pprot]
    prot_idx = [i for i, c in enumerate(ids) if c in pprot]
    out.update(ranking_score=conf["ranking_score"], iptm=conf.get("iptm"), ptm=conf["ptm"],
               has_clash=conf["has_clash"],
               lig_chain_iptm=float(np.max([conf["chain_iptm"][i] for i in lig_idx])),
               lig_prot_pair_iptm=float(np.max([conf["chain_pair_iptm"][i][j] for i in lig_idx for j in prot_idx])))
    out["success"] = bool(out["lig_rmsd"] <= SUCCESS_RMSD)
    return out


def main():
    cases = pd.read_csv(OUT / "cases.tsv", sep="\t")
    rows = []
    for _, c in cases.iterrows():
        for cond in ("dimer", "monomer", "dimer_seed1"):   # dimer_seed1: optional extra seed (on-axis only)
            name = f"{c.case_id}_{cond}"  # AF3 keeps case in output names (checked)
            d = W / "af3" / name
            rs = d / f"{name}_ranking_scores.csv"
            if not rs.exists():
                if cond == "dimer_seed1":
                    continue
                rows.append(dict(case_id=c.case_id, cond=cond, status="missing")); continue
            r = pd.read_csv(rs).sort_values("ranking_score", ascending=False).reset_index(drop=True)
            for rank, x in r.iterrows():
                sd = d / f"seed-{int(x.seed)}_sample-{int(x['sample'])}"
                pre = f"{name}_seed-{int(x.seed)}_sample-{int(x['sample'])}"
                conf = json.loads((sd / f"{pre}_summary_confidences.json").read_text())
                try:
                    m = evaluate(c, cond, sd / f"{pre}_model.cif", conf)
                    st = "ok"
                except Exception as e:  # noqa: BLE001
                    m, st = {}, f"error:{e}"
                rows.append(dict(case_id=c.case_id, stratum=c.stratum, post_cutoff=c.post_cutoff, cond=cond,
                                 seed=int(x.seed), sample=int(x["sample"]), rank=rank, status=st, **m))
            print(c.case_id, cond, flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "metrics_samples.tsv", sep="\t", index=False, float_format="%.3f")
    ok = s[s.status == "ok"]
    top = ok[ok["rank"] == 0].set_index(["case_id", "cond"])
    spread = ok.groupby(["case_id", "cond"]).agg(lig_rmsd_min5=("lig_rmsd", "min"), lig_rmsd_med5=("lig_rmsd", "median"),
                                                 lig_rmsd_max5=("lig_rmsd", "max"), n_success5=("success", "sum"),
                                                 n_samples=("lig_rmsd", "size"))
    met = top.join(spread).reset_index().drop(columns=["rank", "status"])
    miss = s[s.status != "ok"].groupby(["case_id", "cond"]).status.first().reset_index()
    met = pd.concat([met, miss[~miss.set_index(["case_id", "cond"]).index.isin(top.index)]], ignore_index=True)
    met = cases[["case_id", "stratum", "pdb_id", "ccd", "uniprot", "post_cutoff", "release_date"]].merge(
        met.drop(columns=["stratum", "post_cutoff"], errors="ignore"), on="case_id", how="right")
    met.to_csv(OUT / "metrics.tsv", sep="\t", index=False, float_format="%.3f")
    # summary: success (top-ranked) by stratum x condition, overall and split by AF3 training cutoff
    ok_m = met[met.lig_rmsd.notna()].copy()
    ok_m["success"] = ok_m.success.astype(bool)
    ok_m["any5"] = ok_m.n_success5 > 0
    rows = []
    for split, g0 in [("all", ok_m)] + [(f"post_cutoff={v}", g) for v, g in ok_m.groupby("post_cutoff")]:
        for (st, cond), g in g0.groupby(["stratum", "cond"]):
            rows.append(dict(split=split, stratum=st, cond=cond, n=len(g), n_success_top=int(g.success.sum()),
                             n_success_any_of_5=int(g.any5.sum()),
                             n_success_pocket=int((g.lig_rmsd_pocket <= SUCCESS_RMSD).sum()), median_lig_rmsd=g.lig_rmsd.median(),
                             n_contacts_both=int(g.contacts_both.astype(bool).sum()) if cond.startswith("dimer") else np.nan,
                             median_site_dist=g.site_dist.median(), median_c2_axis_dist_pred=g.get("c2_axis_dist_pred", pd.Series(dtype=float)).median()))
    pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False, float_format="%.2f")
    print(met[["case_id", "cond", "lig_rmsd", "success", "contacts_both", "site_dist", "ranking_score"]].to_string())


if __name__ == "__main__":
    main()
