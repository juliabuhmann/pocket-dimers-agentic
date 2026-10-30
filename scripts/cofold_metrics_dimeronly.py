#!/usr/bin/env python3
"""Metrics for deep dive #3 AF3 runs (dimer-only pockets; results/cofold/dimeronly/cases.tsv) against the PDB
assembly-1 pose. Adapted from scripts/cofold_metrics_onaxis.py (geometry helpers imported from there).

Status: first attempt, UNVERIFIED.

Rerun:  .venvs/part1/bin/python scripts/cofold_metrics_dimeronly.py
        -> results/cofold/dimeronly/{metrics,metrics_samples,summary}.tsv

Per prediction sample (5 per run) and for the top-ranked sample:
- Reference = the census ligand instance (ligand_uid, matched by CCD + centroid) in assembly 1 of the PDB entry, and
  the two reference chains of the pair (cases.tsv pchains; per-chain SEQRES offsets in case_chains.json).
- Superposition: protein CA of all mapped chains (dimer: both chain mappings tried for homodimers; monomer: onto
  the reference chain(s) of the same accession). Minimum over mappings and predicted ligand copies.
- lig_rmsd_sym: symmetry-corrected heavy-atom RMSD (spyrmsd 0.9 symmrmsd, graph isomorphisms of the CCD bond
  graph, no extra minimisation), over CCD atom names present in both the reference and the prediction; bonds from
  the RCSB CCD definition (https://files.rcsb.org/ligands/view/<CCD>.cif, cached). lig_rmsd_name = plain name-matched
  RMSD (as in the on-axis analysis). Success = lig_rmsd_sym <= SUCCESS_RMSD.
- lig_rmsd_sym_pocket: same, superposed on CA atoms of reference residues within POCKET_CUTOFF of the reference
  ligand (independent of the global inter-chain arrangement).
- dcc_pocket: distance between predicted and reference ligand centroids after the pocket superposition.
- dcc_curated: distance between the predicted ligand centroid and the centroid of the CA atoms (predicted coords) of
  the curated dimer-only pocket residues (results/part2/q23a_curated.tsv, mapped to the AF3 construct); dimer: both
  role assignments for homodimers, min over ligand copies; monomer: residues of the pocket on that one chain only.
- n_res_chainA/B, contacts_both: residues per protein chain within CONTACT_CUTOFF of the matched predicted copy.
- n_ref_copies / n_ref_copies_recovered: bound reference copies of the CCD on the chain pair and how many of them
  have a predicted copy with sym RMSD <= SUCCESS_RMSD under the global superposition (dimer runs).
- site_frac: fraction of reference site residues (within CONTACT_CUTOFF of the reference ligand, mapped chains) that
  are within CONTACT_CUTOFF of the predicted ligand.
- ref_site_major/minor, ref_site_minor_frac: reference site residues (within CONTACT_CUTOFF of the reference
  ligand) on the more / less contacted chain of the pair; minor_frac = share of the site on the partner chain.
- ca_rmsd, irmsd (dimer): global CA RMSD and interface backbone RMSD (10 A CAPRI interface; DockQ not installed).
- AF3 confidences: ranking_score, iptm, ptm, has_clash, ligand chain_iptm, max ligand-protein chain_pair_iptm.
"""
import json
import sys
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from spyrmsd import rmsd as spy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cofold_metrics_onaxis as OM  # noqa: E402  reuse kabsch, rmsd, heavy, ca_pairs, residues_near, axis_dist, pred_structure
from cofold_metrics_onaxis import BB, PDB_MIRROR, ca_pairs, heavy, kabsch, pred_structure, residues_near, rmsd  # noqa: E402

W = ROOT / "data/cofold/dimeronly"
OUT = ROOT / "results/cofold/dimeronly"
SUCCESS_RMSD = 2.0     # A, task definition (as on-axis)
CONTACT_CUTOFF = 4.0   # A
POCKET_CUTOFF = 10.0   # A
IFACE_CUTOFF = 10.0    # A
DEDUP = 2.0            # A
CCD_URL = "https://files.rcsb.org/ligands/view/{}.cif"
ANUM = {e: gemmi.Element(e).atomic_number for e in ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Se", "B", "Fe", "Mg")}


# ------------------------------------------------------------------ CCD bond graph + symmetry RMSD
_CCD = {}


def ccd_bonds(ccd):
    if ccd in _CCD:
        return _CCD[ccd]
    p = W / "ccd" / f"{ccd}.cif"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(CCD_URL.format(ccd), timeout=60) as r:
            p.write_bytes(r.read())
    b = gemmi.cif.read(str(p)).sole_block()
    t = b.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"])
    _CCD[ccd] = {(r[0].strip('"'), r[1].strip('"')) for r in t}
    return _CCD[ccd]


def sym_rmsd(pl, rl, ccd):
    """symmetry-corrected RMSD (spyrmsd) on common CCD atom names; pl, rl = {name: (element, xyz)} same frame."""
    names = [n for n in rl if n in pl and pl[n][0] == rl[n][0]]
    if len(names) < 3:
        return np.nan
    idx = {n: i for i, n in enumerate(names)}
    A = np.zeros((len(names), len(names)), dtype=int)
    for a, b in ccd_bonds(ccd):
        if a in idx and b in idx:
            A[idx[a], idx[b]] = A[idx[b], idx[a]] = 1
    z = np.array([gemmi.Element(rl[n][0]).atomic_number for n in names])
    X = np.array([rl[n][1] for n in names]); Y = np.array([pl[n][1] for n in names])
    try:
        return float(spy.symmrmsd(X, Y, z, z, A, A, center=False, minimize=False))
    except Exception:  # noqa: BLE001  (disconnected subgraph etc.: fall back to name-matched RMSD)
        return rmsd(X, Y)


def name_rmsd(pl, rl):
    c = [n for n in rl if n in pl and pl[n][0] == rl[n][0]]
    return rmsd(np.array([pl[n][1] for n in c]), np.array([rl[n][1] for n in c])) if c else np.nan


def tr(d, R, t):
    return {n: (e, R @ x + t) for n, (e, x) in d.items()}


def cen(d):
    return np.mean([x for _, x in d.values()], 0)


# ------------------------------------------------------------------ reference
_REF = {}


def ref_structure(case, chains_info):
    if case.case_id in _REF:
        return _REF[case.case_id]
    st = gemmi.read_structure(str(PDB_MIRROR / case.pdb_id[1:3] / f"{case.pdb_id}.cif.gz"))
    st.setup_entities()
    a = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    a.remove_hydrogens()
    prot = {}
    for ci in chains_info:
        lo = ci["seqres_lo"]
        prot[ci["chain"]] = {r.label_seq - 1 - lo: heavy(r) for r in a[ci["chain"]].get_polymer() if r.label_seq is not None}
    ligs = []
    for ch in a:
        for r in ch:
            if r.name == case.ccd and r.het_flag == "H":
                h = heavy(r)
                if h:
                    ligs.append(h)
    c0 = np.array(json.loads(case.centroid))
    ligs.sort(key=lambda h: np.linalg.norm(cen(h) - c0))
    uniq = []
    for h in ligs:
        if all(np.linalg.norm(cen(h) - cen(u)) > DEDUP for u in uniq):
            uniq.append(h)
    # bound to the pair chains, same rule as the selection: >= 1 contact residue and a majority of the contacting
    # protein residues (any chain of the assembly) on the pair chains
    ns = gemmi.NeighborSearch(a, st.cell, 6).populate()
    pair = set(prot)
    bound = []
    for h in uniq:
        res = set()
        for _, x in h.values():
            for mk in ns.find_atoms(gemmi.Position(*x), "\0", radius=CONTACT_CUTOFF):
                cra = mk.to_cra(a)
                if cra.residue.het_flag != "H" and cra.residue.entity_type == gemmi.EntityType.Polymer:
                    res.add((cra.chain.name, cra.residue.seqid.num))
        on = sum(c in pair for c, _ in res)
        if on and 2 * on > len(res):
            bound.append(h)
    _REF[case.case_id] = (prot, bound)
    return prot, bound


def evaluate(case, cond, model_path, conf, chains_info, pocket):
    rprot, rligs = ref_structure(case, chains_info)
    pprot, pligs, finite = pred_structure(model_path, case.ccd)
    out = dict(finite=finite, n_pred_lig=len(pligs), n_ref_copies=len(rligs))
    if not finite or not pligs or not rligs:
        return out
    rlig = rligs[0]
    rlig_xyz = np.array([v[1] for v in rlig.values()])
    rch = [c["chain"] for c in chains_info]
    racc = {c["chain"]: c["acc"] for c in chains_info}
    pch = sorted(pprot)
    homo = case.kind == "homo"
    if cond == "dimer":
        maps = [dict(zip(pch, rch))] + ([dict(zip(pch, rch[::-1]))] if homo else [])
    else:
        maps = [{pch[0]: rc} for rc in rch if racc[rc] == racc[case.mono_chain]]
    best = None
    for m in maps:
        X = np.vstack([ca_pairs(pprot[p], rprot[r])[0] for p, r in m.items()])
        Y = np.vstack([ca_pairs(pprot[p], rprot[r])[1] for p, r in m.items()])
        R, t = kabsch(X, Y)
        ca = rmsd(X @ R.T + t, Y)
        for li, pl in enumerate(pligs):
            plt_ = tr(pl, R, t)
            rs = sym_rmsd(plt_, rlig, case.ccd)
            if best is None or rs < best["lig_rmsd_sym"]:
                best = dict(lig_rmsd_sym=rs, lig_rmsd_name=name_rmsd(plt_, rlig), ca_rmsd=ca, _m=m, _li=li, _R=R, _t=t)
    m, li, R, t = best.pop("_m"), best.pop("_li"), best.pop("_R"), best.pop("_t")
    out.update(best)
    out["n_atoms_ref"] = len(rlig)
    if cond == "dimer":
        rec = 0
        for rl in rligs:
            rec += min(sym_rmsd(tr(pl, R, t), rl, case.ccd) for pl in pligs) <= SUCCESS_RMSD
        out["n_ref_copies_recovered"] = int(rec)
    # pocket superposition
    pk, dcc = None, None
    for mm in maps:
        X, Y = [], []
        for p, r in mm.items():
            keys = sorted(residues_near(rprot[r], rlig_xyz, POCKET_CUTOFF))
            x, y = ca_pairs(pprot[p], rprot[r], ("CA",), keys)
            if len(x):
                X.append(x); Y.append(y)
        if not X or sum(len(x) for x in X) < 3:
            continue
        R2, t2 = kabsch(np.vstack(X), np.vstack(Y))
        for pl_ in pligs:
            q = tr(pl_, R2, t2)
            rs = sym_rmsd(q, rlig, case.ccd)
            if pk is None or rs < pk:
                pk, dcc = rs, float(np.linalg.norm(cen(q) - cen(rlig)))
    out["lig_rmsd_sym_pocket"], out["dcc_pocket"] = pk, dcc
    pl = pligs[li]
    pxyz = np.array([v[1] for v in pl.values()])
    pcen = pxyz.mean(0)
    for i, p in enumerate(pch):
        out[f"n_res_chain{'AB'[i]}"] = len(residues_near(pprot[p], pxyz, CONTACT_CUTOFF))
    out["contacts_both"] = bool(cond == "dimer" and out["n_res_chainA"] >= 1 and out.get("n_res_chainB", 0) >= 1)
    site_hit = site_n = 0
    for p, r in m.items():
        site = residues_near(rprot[r], rlig_xyz, CONTACT_CUTOFF)
        near = residues_near(pprot[p], pxyz, CONTACT_CUTOFF)
        site_n += sum(k in pprot[p] for k in site); site_hit += sum(k in near for k in site)
    out["site_frac"] = site_hit / site_n if site_n else np.nan
    # reference site size on each reference chain of the pair (minor = share of the site on the less-contacted chain;
    # the part of the site a single-chain run cannot reproduce)
    rs_n = sorted((len(residues_near(rprot[r], rlig_xyz, CONTACT_CUTOFF)) for r in rch), reverse=True)
    out["ref_site_major"], out["ref_site_minor"] = rs_n[0], rs_n[1]
    out["ref_site_minor_frac"] = rs_n[1] / sum(rs_n) if sum(rs_n) else np.nan
    # curated dimer-only pocket centre (construct indices via per-chain UniProt mapping)
    u2i = {c["chain"]: {u: i for i, u in enumerate(c["mapping"]) if u is not None} for c in chains_info}
    pred_of_ref = {r: p for p, r in m.items()}
    roles = [("A", "B"), ("B", "A")] if homo else [("A", "B")]
    best_d = None
    for ra, rb in roles:
        cas = []
        for role, rc in zip((ra, rb), rch):
            if rc not in pred_of_ref:
                continue
            pc_ = pprot[pred_of_ref[rc]]
            for pos in pocket[role][1]:
                i = u2i[rc].get(pos - 1)
                if i is not None and i in pc_ and "CA" in pc_[i]:
                    cas.append(pc_[i]["CA"][1])
        if not cas:
            continue
        pc0 = np.mean(cas, 0)
        for pl_ in pligs:
            d = float(np.linalg.norm(cen(pl_) - pc0))
            best_d = d if best_d is None or d < best_d else best_d
    out["dcc_curated"] = best_d
    if cond == "dimer":
        out["pred_rot_angle"], out["c2_axis_dist_pred"] = OM.axis_dist(pprot[pch[0]], pprot[pch[1]], pcen)
        ra, rb = m[pch[0]], m[pch[1]]
        rxyz = {c: np.array([x for v in rprot[c].values() for _, x in v.values()]) for c in (ra, rb)}
        ia = residues_near(rprot[ra], rxyz[rb], IFACE_CUTOFF)
        ib = residues_near(rprot[rb], rxyz[ra], IFACE_CUTOFF)
        X1, Y1 = ca_pairs(pprot[pch[0]], rprot[ra], BB, sorted(ia))
        X2, Y2 = ca_pairs(pprot[pch[1]], rprot[rb], BB, sorted(ib))
        X, Y = np.vstack([X1, X2]), np.vstack([Y1, Y2])
        R3, t3 = kabsch(X, Y)
        out["irmsd"] = rmsd(X @ R3.T + t3, Y)
    ids = conf["chain_ids"]
    lig_idx = [i for i, c in enumerate(ids) if c not in pprot]
    prot_idx = [i for i, c in enumerate(ids) if c in pprot]
    out.update(ranking_score=conf["ranking_score"], iptm=conf.get("iptm"), ptm=conf["ptm"], has_clash=conf["has_clash"],
               lig_chain_iptm=float(np.max([conf["chain_iptm"][i] for i in lig_idx])),
               lig_prot_pair_iptm=float(np.max([conf["chain_pair_iptm"][i][j] for i in lig_idx for j in prot_idx])))
    out["success"] = bool(out["lig_rmsd_sym"] <= SUCCESS_RMSD)
    return out


def main():
    cases = pd.read_csv(OUT / "cases.tsv", sep="\t")
    chains = json.loads((W / "case_chains.json").read_text())
    pockets = json.loads((W / "case_pockets.json").read_text())
    rows = []
    for _, c in cases.iterrows():
        for cond in ("dimer", "monomer"):
            name = f"{c.case_id}_{cond}"
            d = W / "af3" / name
            rs = d / f"{name}_ranking_scores.csv"
            if not rs.exists():
                rows.append(dict(case_id=c.case_id, cond=cond, status="missing")); continue
            r = pd.read_csv(rs).sort_values("ranking_score", ascending=False).reset_index(drop=True)
            for rank, x in r.iterrows():
                sd = d / f"seed-{int(x.seed)}_sample-{int(x['sample'])}"
                pre = f"{name}_seed-{int(x.seed)}_sample-{int(x['sample'])}"
                conf = json.loads((sd / f"{pre}_summary_confidences.json").read_text())
                try:
                    mt = evaluate(c, cond, sd / f"{pre}_model.cif", conf, chains[c.case_id], pockets[c.case_id])
                    st = "ok"
                except Exception as e:  # noqa: BLE001
                    mt, st = {}, f"error:{e}"
                rows.append(dict(case_id=c.case_id, cond=cond, seed=int(x.seed), sample=int(x["sample"]), rank=rank,
                                 status=st, **mt))
            print(c.case_id, cond, flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "metrics_samples.tsv", sep="\t", index=False, float_format="%.3f")
    ok = s[s.status == "ok"]
    top = ok[ok["rank"] == 0].set_index(["case_id", "cond"])
    spread = ok.groupby(["case_id", "cond"]).agg(lig_rmsd_sym_min5=("lig_rmsd_sym", "min"),
                                                 lig_rmsd_sym_med5=("lig_rmsd_sym", "median"),
                                                 n_success5=("success", "sum"), n_samples=("lig_rmsd_sym", "size"))
    met = top.join(spread).reset_index().drop(columns=["rank", "status"])
    miss = s[s.status != "ok"].groupby(["case_id", "cond"]).status.first().reset_index()
    met = pd.concat([met, miss[~miss.set_index(["case_id", "cond"]).index.isin(top.index)]], ignore_index=True)
    met = cases[["case_id", "kind", "pdb_id", "ccd", "ligand_class", "post_cutoff", "release_date", "n_lig_dimer_run"]].merge(
        met, on="case_id", how="right")
    met.to_csv(OUT / "metrics.tsv", sep="\t", index=False, float_format="%.3f")
    g = met[met.lig_rmsd_sym.notna()]
    summ = []
    for split, g0 in [("all", g)] + [(f"post_cutoff={v}", x) for v, x in g.groupby("post_cutoff")]:
        for cond, x in g0.groupby("cond"):
            summ.append(dict(split=split, cond=cond, n=len(x), n_success_top=int(x.success.astype(bool).sum()),
                             n_success_any_of_5=int((x.n_success5 > 0).sum()),
                             n_success_pocket=int((x.lig_rmsd_sym_pocket <= SUCCESS_RMSD).sum()),
                             median_lig_rmsd_sym=x.lig_rmsd_sym.median(), median_dcc_pocket=x.dcc_pocket.median(),
                             median_dcc_curated=x.dcc_curated.median(),
                             n_contacts_both=int(x.contacts_both.astype(bool).sum()) if cond == "dimer" else np.nan))
    pd.DataFrame(summ).to_csv(OUT / "summary.tsv", sep="\t", index=False, float_format="%.2f")
    print(met[["case_id", "cond", "lig_rmsd_sym", "lig_rmsd_name", "lig_rmsd_sym_pocket", "dcc_pocket", "dcc_curated",
               "contacts_both", "ca_rmsd", "ranking_score"]].to_string())


if __name__ == "__main__":
    main()
