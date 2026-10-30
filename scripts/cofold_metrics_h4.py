#!/usr/bin/env python3
"""Metrics for deep dive #2 AF3 runs (results/cofold/h4/cases.tsv) against the PDB assembly-1 pose.

Status: first attempt, UNVERIFIED.

Rerun:  .venvs/part1/bin/python scripts/cofold_metrics_h4.py  -> results/cofold/h4/{metrics,metrics_samples,summary}.tsv

Reuses (imports, unmodified) helpers from deep dive #1 (scripts/cofold_metrics_onaxis.py): kabsch, rmsd, heavy,
lig_rmsds (CCD-name-matched RMSD and Hungarian-within-element RMSD), residues_near. Extended here to heterodimers
(two entities, per-chain SEQRES offsets) and to the case-specific PDB chain pair chosen by cofold_select_h4.py.

Per sample (5 per run) and for the top-ranked sample:
- lig_rmsd: SYMMETRY-CORRECTED heavy-atom RMSD (min over ligand graph automorphisms, sym_rmsd) of a predicted
  ligand copy vs the census reference instance, after superposing protein
  CA atoms of the pair (dimer: for homodimers both chain mappings are tried; monomer: onto the PDB chain with the
  most contacts, for homodimers also onto the other copy). Minimum over chain mappings and predicted copies.
  lig_rmsd_name: atoms matched by CCD atom name only (no symmetry correction). lig_rmsd_hung: Hungarian-within-element match (lower bound for symmetry-corrected RMSD).
  success = lig_rmsd <= SUCCESS_RMSD.
- dcc: distance between the predicted (superposed) ligand centroid and the PDB ligand centroid (pocket-level).
  dcc_hit = dcc <= DCC_HIT (same 4 A threshold as H4).
- contacts_both: dimer run, the matched predicted copy has >= 1 heavy atom within CONTACT_CUTOFF of both chains.
- site_frac: fraction of PDB site residues (within CONTACT_CUTOFF of the reference ligand) that are within
  CONTACT_CUTOFF of the predicted ligand (residues of the chains present in the run); site_frac_full: same
  numerator over the site residues of BOTH PDB chains (monomer run cannot exceed its own chain's share).
- ca_rmsd (superposition fit), ranking_score, iptm, ptm, has_clash, ligand chain_iptm and max ligand-protein
  chain_pair_iptm from AF3 summary confidences.
"""
import os
import json
import sys
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from cofold_metrics_onaxis import heavy, kabsch, lig_rmsds, residues_near, rmsd  # noqa: E402

W = ROOT / "data/cofold/h4"
OUT = ROOT / "results/cofold/h4"
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/cluster/scratch/jjaenes/pdb"))  # = src/part1/config.py PDB_MIRROR
SUCCESS_RMSD = 2.0   # task definition
DCC_HIT = 4.0        # H4 pocket-level threshold (scripts/h4_common.py DCC_HIT)
CONTACT_CUTOFF = 4.0


MAX_AUTOMORPHISMS = 20000   # cap on enumerated ligand graph automorphisms


def sym_rmsd(pl, rl):
    """Symmetry-corrected RMSD: min over graph automorphisms of the reference ligand (bonds inferred from the
    reference geometry: distance <= sum of covalent radii + 0.45 A; nodes labelled by element)."""
    names = [n for n in rl if n in pl and pl[n][0] == rl[n][0]]
    if not names:
        return np.nan
    G = nx.Graph()
    for n in names:
        G.add_node(n, el=rl[n][0])
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            cut = gemmi.Element(rl[a][0]).covalent_r + gemmi.Element(rl[b][0]).covalent_r + 0.45
            if np.linalg.norm(rl[a][1] - rl[b][1]) <= cut:
                G.add_edge(a, b)
    R = np.array([rl[n][1] for n in names])
    best = np.inf
    gm = nx.algorithms.isomorphism.GraphMatcher(G, G, node_match=lambda x, y: x["el"] == y["el"])
    for k, iso in enumerate(gm.isomorphisms_iter()):
        if k >= MAX_AUTOMORPHISMS:
            break
        P = np.array([pl[iso[n]][1] for n in names])
        best = min(best, rmsd(P, R))
    return float(best)


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
    prot = {}
    for role, cn, lo in (("A", case.chainA_pdb, case.loA), ("B", case.chainB_pdb, case.loB)):
        prot[role] = {r.label_seq - 1 - lo: heavy(r) for r in a[cn].get_polymer() if r.label_seq is not None}
    cen = np.array(json.loads(case.centroid))
    best, bd = None, 1e9
    for ch in a:
        for r in ch:
            if r.name == case.ccd and r.het_flag == "H":
                h = heavy(r)
                if h:
                    d = np.linalg.norm(np.mean([v[1] for v in h.values()], 0) - cen)
                    if d < bd:
                        best, bd = h, d
    return prot, best


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


def pairs(P, R):
    k = sorted(set(P) & set(R))
    k = [x for x in k if "CA" in P[x] and "CA" in R[x]]
    return np.array([P[x]["CA"][1] for x in k]), np.array([R[x]["CA"][1] for x in k])


def evaluate(case, cond, model_path, conf):
    rprot, rlig = ref_structure(case)
    pprot, pligs, finite = pred_structure(model_path, case.ccd)
    out = dict(finite=finite, n_pred_lig=len(pligs))
    if not finite or not pligs or rlig is None:
        return out
    rxyz = np.array([v[1] for v in rlig.values()])
    rcen = rxyz.mean(0)
    pch = sorted(pprot)       # AF3 chain ids: A (and B); A = PDB chain with most contacts (chainA_pdb)
    homo = case.accA == case.accB
    if cond == "dimer":
        maps = [{pch[0]: "A", pch[1]: "B"}] + ([{pch[0]: "B", pch[1]: "A"}] if homo else [])
    else:
        maps = [{pch[0]: "A"}] + ([{pch[0]: "B"}] if homo else [])
    best = None
    for m in maps:
        XY = [pairs(pprot[p], rprot[r]) for p, r in m.items()]
        X, Y = np.vstack([x for x, _ in XY]), np.vstack([y for _, y in XY])
        R, t = kabsch(X, Y)
        ca = rmsd(X @ R.T + t, Y)
        for li, pl in enumerate(pligs):
            plt_ = {n: (e, R @ x + t) for n, (e, x) in pl.items()}
            rn, rh, nm, nr = lig_rmsds(plt_, rlig)
            rs = sym_rmsd(plt_, rlig)
            key = rs if np.isfinite(rs) else rh
            if best is None or key < best[0]:
                best = (key, dict(lig_rmsd=rs, lig_rmsd_name=rn, lig_rmsd_hung=rh, n_atoms_matched=nm, n_atoms_ref=nr, ca_rmsd=ca,
                                  n_ca_sup=len(X)), m, plt_)
    _, res, m, plt_ = best
    out.update(res)
    pxyz = np.array([v[1] for v in plt_.values()])
    out["dcc"] = float(np.linalg.norm(pxyz.mean(0) - rcen))
    out["dcc_hit"] = bool(out["dcc"] <= DCC_HIT)
    # map the predicted protein into the PDB frame with the chosen mapping (contacts are frame-invariant)
    XY = [pairs(pprot[p], rprot[r]) for p, r in m.items()]
    R, t = kabsch(np.vstack([x for x, _ in XY]), np.vstack([y for _, y in XY]))
    pprot_t = {p: {k: {n: (e, R @ x + t) for n, (e, x) in v.items()} for k, v in pprot[p].items()} for p in pch}
    nres = {m[p]: len(residues_near(pprot_t[p], pxyz, CONTACT_CUTOFF)) for p in pch}
    out["n_res_A"] = nres.get("A", 0)
    out["n_res_B"] = nres.get("B", 0)
    out["contacts_both"] = bool(cond == "dimer" and out["n_res_A"] >= 1 and out["n_res_B"] >= 1)
    site_n = site_hit = 0
    for p in pch:
        site = residues_near(rprot[m[p]], rxyz, CONTACT_CUTOFF)
        near = residues_near(pprot_t[p], pxyz, CONTACT_CUTOFF)
        site_n += len(site); site_hit += len(site & near)
    out["site_frac"] = site_hit / site_n if site_n else np.nan
    site_all = sum(len(residues_near(rprot[r], rxyz, CONTACT_CUTOFF)) for r in ("A", "B"))
    out["site_frac_full"] = site_hit / site_all if site_all else np.nan
    ids = conf["chain_ids"]
    lig_idx = [i for i, c in enumerate(ids) if c not in pprot]
    prot_idx = [i for i, c in enumerate(ids) if c in pprot]
    out.update(ranking_score=conf["ranking_score"], iptm=conf.get("iptm"), ptm=conf["ptm"], has_clash=conf["has_clash"],
               lig_chain_iptm=float(np.max([conf["chain_iptm"][i] for i in lig_idx])),
               lig_prot_pair_iptm=float(np.max([conf["chain_pair_iptm"][i][j] for i in lig_idx for j in prot_idx])))
    out["success"] = bool(out["lig_rmsd"] <= SUCCESS_RMSD) if np.isfinite(out["lig_rmsd"]) else False
    return out


def main():
    cases = pd.read_csv(OUT / "cases.tsv", sep="\t")
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
                pre = f"{name}_seed-{int(x.seed)}_sample-{int(x['sample'])}"
                sd = d / f"seed-{int(x.seed)}_sample-{int(x['sample'])}"
                conf = json.loads((sd / f"{pre}_summary_confidences.json").read_text())
                try:
                    m, st = evaluate(c, cond, sd / f"{pre}_model.cif", conf), "ok"
                except Exception as e:  # noqa: BLE001
                    m, st = {}, f"error:{e}"
                rows.append(dict(case_id=c.case_id, cond=cond, seed=int(x.seed), sample=int(x["sample"]), rank=rank,
                                 status=st, **m))
            print(c.case_id, cond, flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "metrics_samples.tsv", sep="\t", index=False, float_format="%.3f")
    ok = s[s.status == "ok"]
    top = ok[ok["rank"] == 0].set_index(["case_id", "cond"])
    spread = ok.groupby(["case_id", "cond"]).agg(lig_rmsd_min5=("lig_rmsd", "min"), lig_rmsd_med5=("lig_rmsd", "median"),
                                                 n_success5=("success", "sum"), dcc_min5=("dcc", "min"),
                                                 n_samples=("lig_rmsd", "size"))
    met = top.join(spread).reset_index().drop(columns=["rank"])
    miss = s[(s.status != "ok") & (s.get("rank", 0).fillna(0) == 0)][["case_id", "cond", "status"]]
    met = pd.concat([met, miss[~miss.set_index(["case_id", "cond"]).index.isin(top.index)]], ignore_index=True)
    keep = ["case_id", "stratum", "forced", "pdb_id", "ccd", "ligand_class", "kind", "accA", "accB", "release_date",
            "post_cutoff", "post_cutoff_pair", "apo_dimer_hit", "apo_iso_hit", "apo_dimer_dcc", "fm_dimer",
            "n_lig_dimer_run", "tokens_dimer"]
    met = cases[keep].merge(met, on="case_id", how="right")
    met.to_csv(OUT / "metrics.tsv", sep="\t", index=False, float_format="%.3f")

    # summary by stratum x cond, overall and split by release relative to the AF3 cutoff
    v = met[met.status == "ok"].copy()
    for col in ("success", "dcc_hit", "contacts_both"):
        v[col] = v[col].astype(bool)
    rows = []
    for split, g0 in [("all", v), ("post_cutoff", v[v.post_cutoff]), ("pre_cutoff", v[~v.post_cutoff])]:
        for st, g1 in [("all", g0)] + list(g0.groupby("stratum")):
            for cond, g in g1.groupby("cond"):
                rows.append(dict(split=split, stratum=st, cond=cond, n=len(g), n_success=int(g.success.sum()),
                                 n_success_any5=int((g.n_success5 > 0).sum()), n_dcc_hit=int(g.dcc_hit.sum()),
                                 n_contacts_both=int(g.contacts_both.sum()) if cond == "dimer" else np.nan,
                                 median_lig_rmsd=g.lig_rmsd.median(), median_dcc=g.dcc.median(),
                                 median_ranking=g.ranking_score.median(),
                                 median_lig_iptm=g.lig_prot_pair_iptm.median()))
    pd.DataFrame(rows).to_csv(OUT / "summary.tsv", sep="\t", index=False, float_format="%.3f")
    # paired A vs B per case
    pv = v.pivot(index="case_id", columns="cond", values="lig_rmsd").dropna()
    print("paired cases", len(pv), "dimer better", int((pv.dimer < pv.monomer).sum()))
    print(met[["case_id", "cond", "lig_rmsd", "success", "dcc", "contacts_both", "ranking_score",
               "lig_prot_pair_iptm"]].to_string())


if __name__ == "__main__":
    main()
