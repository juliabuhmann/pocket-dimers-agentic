"""Independent re-scoring of AF3 co-folding top-ranked predictions (verifier code; does not import cofold_* scripts).

First attempt -- outputs are NOT verified beyond the checks reported in reviews/verify_cofold.md.

For every (case, condition, seed) AF3 run in the on-axis, dimer-only and H4 batches:
  * pick the top-ranked sample from <job>_ranking_scores.csv (max ranking_score for the requested seed);
  * reference = PDB assembly 1 (gemmi make_assembly, AddNumber naming, H removed, first altloc);
    the two reference protein chains are the case table's chain columns (input definition, not a result);
  * residue mapping: global sequence alignment (gemmi) of each predicted chain to each candidate reference chain's
    *observed* residues; only identical aligned residues are paired;
  * superposition: Kabsch on all paired CA atoms of all predicted protein chains (global CA);
  * protein chain assignment: every assignment of predicted chains to reference chains with sequence identity
    >= MIN_ID_ASSIGN is tried; reported RMSD = min over assignments (protein-symmetry equivalent), and the RMSD
    under the CA-RMSD-best assignment is also kept;
  * ligand: CCD heavy-atom bond graph from the RCSB CCD file (fetched, cached); all element-preserving graph
    automorphisms (networkx VF2) are applied; RMSD over atoms present in both copies; min over
    (predicted copy, reference copy) pairs, reference copies = residues with that CCD id within REF_COPY_CUT of
    either reference protein chain.
Sanity: finite coordinates, min ligand-protein heavy-atom distance (clash < CLASH_CUT), duplicate predicted copies
(centroid distance < DUP_CUT), and partner share of the PDB site (BioLiP-like vdW+0.5 A, >=2 atom pairs).
"""
import argparse
import itertools
import json
import sys
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PDB_MIRROR = Path("/path/to/pdb_mirror")
CCD_URL = "https://files.rcsb.org/ligands/download/{}.cif"
CCD_CACHE = ROOT / "reviews/verify_cofold/ccd_cache"   # small CIFs, re-fetchable
OUT = ROOT / "reviews/verify_cofold"

SUCCESS_CUT = 2.0        # A, task definition
CLASH_CUT = 1.5          # A, task definition
DUP_CUT = 2.0            # A, centroid distance between two predicted copies treated as duplicate
REF_COPY_CUT = 5.0       # A, reference ligand copy must have a heavy atom this close to a reference protein chain
MIN_ID_ASSIGN = 0.9      # identity over aligned pairs for a predicted->reference chain assignment
BIOLIP_TOL = 0.5         # A, BioLiP: d <= vdW_i + vdW_j + 0.5
BIOLIP_MIN_PAIRS = 2     # BioLiP: >= 2 atom pairs per residue
MAX_AUTOMORPHISMS = 200000

BATCHES = {
    "onaxis": dict(cases="results/cofold/onaxis/cases.tsv", af3="data/cofold/onaxis/af3",
                   prod="results/cofold/onaxis/metrics.tsv", prod_rmsd="lig_rmsd"),
    "dimeronly": dict(cases="results/cofold/dimeronly/cases.tsv", af3="data/cofold/dimeronly/af3",
                      prod="results/cofold/dimeronly/metrics.tsv", prod_rmsd="lig_rmsd_sym"),
    "h4": dict(cases="results/cofold/h4/cases.tsv", af3="data/cofold/h4/af3",
               prod="results/cofold/h4/metrics.tsv", prod_rmsd="lig_rmsd"),
}
SEED_DIRS = ["data/cofold/{batch}/af3", "data/cofold/seeds/af3"]


def ref_chain_names(batch, row):
    if batch == "h4":
        return [row.chainA_pdb, row.chainB_pdb]
    return row.pchains.split(",")


def centroid_of(row):
    return np.array(json.loads(row.centroid), float)


# ------------------------------------------------------------------ CCD graph
def ccd_graph(ccd):
    CCD_CACHE.mkdir(parents=True, exist_ok=True)
    f = CCD_CACHE / f"{ccd}.cif"
    if not f.exists():
        with urllib.request.urlopen(CCD_URL.format(ccd), timeout=60) as r:
            f.write_bytes(r.read())
    doc = gemmi.cif.read(str(f))
    cc = gemmi.make_chemcomp_from_block(doc[0])
    G = nx.Graph()
    for a in cc.atoms:
        if not a.is_hydrogen():
            G.add_node(a.id, el=a.el.name)
    for b in cc.rt.bonds:
        if b.id1.atom in G and b.id2.atom in G:
            G.add_edge(b.id1.atom, b.id2.atom)
    return G


def automorphisms(G):
    gm = nx.algorithms.isomorphism.GraphMatcher(G, G, node_match=lambda a, b: a["el"] == b["el"])
    out = []
    for m in gm.isomorphisms_iter():
        out.append(m)
        if len(out) >= MAX_AUTOMORPHISMS:
            return out, True
    return out, False


# ------------------------------------------------------------------ structures
def one_letter(resname):
    t = gemmi.find_tabulated_residue(resname)
    return t.one_letter_code.upper() if t and t.one_letter_code != " " else "X"


def protein_chain(ch):
    """List of (one_letter, {atom_name: xyz}) for polymer residues with a CA."""
    out = []
    for r in ch:
        if r.het_flag == "H" and not gemmi.find_tabulated_residue(r.name):
            continue
        t = gemmi.find_tabulated_residue(r.name)
        if t is None or not t.is_amino_acid():
            continue
        atoms = {a.name: np.array(a.pos.tolist()) for a in r if a.element.name not in ("H", "D")}
        if "CA" in atoms:
            out.append((one_letter(r.name), atoms, r.seqid.num))
    return out


def ligand_copies(model, ccd, chains=None):
    out = []
    for ch in model:
        if chains is not None and ch.name not in chains:
            continue
        for r in ch:
            if r.name == ccd:
                atoms = {a.name: np.array(a.pos.tolist()) for a in r if a.element.name not in ("H", "D")}
                els = {a.name: a.element for a in r if a.element.name not in ("H", "D")}
                out.append(dict(chain=ch.name, seq=r.seqid.num, atoms=atoms, els=els))
    return out


_REF = {}


def load_ref(pdb_id):
    if pdb_id in _REF:
        return _REF[pdb_id]
    st = gemmi.read_structure(str(PDB_MIRROR / pdb_id[1:3] / f"{pdb_id}.cif.gz"))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    a = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    _REF[pdb_id] = a
    return a


def load_pred(path):
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    return st[0]


def align_pairs(pred_res, ref_res):
    """Indices (i_pred, i_ref) of identical aligned residues; identity over aligned columns."""
    s1 = "".join(x[0] for x in pred_res)
    s2 = "".join(x[0] for x in ref_res)
    al = gemmi.align_string_sequences(list(s1), list(s2), [])
    g1, g2 = al.add_gaps(s1, 1), al.add_gaps(s2, 2)
    i = j = 0
    pairs, nal, nid = [], 0, 0
    for c1, c2 in zip(g1, g2):
        if c1 != "-" and c2 != "-":
            nal += 1
            if c1 == c2:
                nid += 1
                pairs.append((i, j))
        if c1 != "-":
            i += 1
        if c2 != "-":
            j += 1
    return pairs, (nid / nal if nal else 0.0)


def kabsch(P, Q):
    """R, t minimising |P R + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, qc - pc @ R


def lig_rmsd(pred, ref, autos, R, t):
    """Min RMSD over CCD automorphisms (atoms present in both). Returns (rmsd, n_atoms, rmsd_name_matched)."""
    P = {k: v @ R + t for k, v in pred["atoms"].items()}
    best, nbest, name_rmsd = np.inf, 0, np.nan
    for m in autos:
        common = [k for k in P if k in m and m[k] in ref["atoms"]]
        if not common:
            continue
        d = np.array([P[k] - ref["atoms"][m[k]] for k in common])
        r = float(np.sqrt((d ** 2).sum(1).mean()))
        if all(m[k] == k for k in m):
            name_rmsd = r
        if r < best:
            best, nbest = r, len(common)
    return best, nbest, name_rmsd


def biolip_site(lig, chains_res):
    """Binding residues per chain name: >= BIOLIP_MIN_PAIRS atom pairs with d <= vdW_i+vdW_j+tol."""
    L = np.array(list(lig["atoms"].values()))
    Lr = np.array([lig["els"][k].vdw_r for k in lig["atoms"]])
    out = {}
    for cname, ch in chains_res.items():
        n = 0
        for r in ch:
            hits = 0
            for a in r:
                if a.element.name in ("H", "D"):
                    continue
                d = np.linalg.norm(L - np.array(a.pos.tolist()), axis=1)
                hits += int((d <= Lr + a.element.vdw_r + BIOLIP_TOL).sum())
                if hits >= BIOLIP_MIN_PAIRS:
                    break
            if hits >= BIOLIP_MIN_PAIRS:
                t = gemmi.find_tabulated_residue(r.name)
                if t is not None and t.is_amino_acid():
                    n += 1
        out[cname] = n
    return out


# ------------------------------------------------------------------ one job
def top_sample(jobdir, job, seed):
    rs = pd.read_csv(jobdir / f"{job}_ranking_scores.csv")
    rs = rs[rs["seed"] == seed]
    if rs.empty:
        return None, None
    best = rs.loc[rs.ranking_score.idxmax()]
    s, k = int(best.seed), int(best["sample"])
    return jobdir / f"seed-{s}_sample-{k}" / f"{job}_seed-{s}_sample-{k}_model.cif", float(best.ranking_score)


def score_job(args):
    batch, row, cond, seed, jobdir, job = args
    rec = dict(batch=batch, case_id=row.case_id, cond=cond, seed=seed, job=job)
    try:
        path, rscore = top_sample(jobdir, job, seed)
        if path is None or not path.exists():
            rec["status"] = "missing"
            return rec
        rec.update(model=str(path.relative_to(ROOT)), ranking_score=rscore)
        pm = load_pred(path)
        ref = load_ref(row.pdb_id)
        rnames = ref_chain_names(batch, row)
        rch = {n: ref.find_chain(n) for n in rnames}
        if any(c is None for c in rch.values()):
            rec["status"] = f"ref_chain_missing:{rnames}"
            return rec
        rres = {n: protein_chain(c) for n, c in rch.items()}
        pres = {ch.name: protein_chain(ch) for ch in pm if len(protein_chain(ch)) > 0}
        # finite
        allxyz = np.array([a.pos.tolist() for ch in pm for r in ch for a in r])
        rec["finite"] = bool(np.isfinite(allxyz).all())
        # candidate assignments
        cand = {}
        for pn, pr in pres.items():
            for rn, rr in rres.items():
                pairs, ident = align_pairs(pr, rr)
                if ident >= MIN_ID_ASSIGN and len(pairs) >= 20:
                    cand[(pn, rn)] = pairs
        pnames = sorted(pres)
        assigns = []
        for perm in itertools.permutations(rnames, len(pnames)):
            if all((p, r) in cand for p, r in zip(pnames, perm)):
                assigns.append(dict(zip(pnames, perm)))
        if not assigns:
            rec["status"] = "no_chain_assignment"
            return rec
        G = ccd_graph(row.ccd)
        autos, trunc = automorphisms(G)
        rec["n_automorphisms"], rec["auto_truncated"] = len(autos), trunc
        plig = ligand_copies(pm, row.ccd)
        rlig_all = ligand_copies(ref, row.ccd)
        # reference copies near either reference protein chain
        prot_xyz = np.array([a for ch in rres.values() for _, atoms, _ in ch for a in atoms.values()])
        rlig = [c for c in rlig_all
                if np.min(np.linalg.norm(prot_xyz[None] - np.array(list(c["atoms"].values()))[:, None], axis=2))
                <= REF_COPY_CUT]
        rec["n_pred_lig"], rec["n_ref_lig"] = len(plig), len(rlig)
        if not plig or not rlig:
            rec["status"] = "no_ligand"
            return rec
        results = []
        for asg in assigns:
            P, Q = [], []
            for pn, rn in asg.items():
                for i, j in cand[(pn, rn)]:
                    P.append(pres[pn][i][1]["CA"])
                    Q.append(rres[rn][j][1]["CA"])
            P, Q = np.array(P), np.array(Q)
            R, t = kabsch(P, Q)
            ca = float(np.sqrt(((P @ R + t - Q) ** 2).sum(1).mean()))
            best = (np.inf, None, None, 0, np.nan)
            for pi, pl in enumerate(plig):
                for ri, rl in enumerate(rlig):
                    r, n, nm = lig_rmsd(pl, rl, autos, R, t)
                    if r < best[0]:
                        best = (r, pi, ri, n, nm)
            results.append(dict(asg=asg, ca=ca, n_ca=len(P), lig=best[0], pi=best[1], ri=best[2],
                                n_atoms=best[3], name_rmsd=best[4]))
        bmin = min(results, key=lambda x: x["lig"])
        bca = min(results, key=lambda x: x["ca"])
        rec.update(status="ok", n_assign=len(assigns),
                   lig_rmsd=bmin["lig"], lig_rmsd_name=bmin["name_rmsd"], n_atoms=bmin["n_atoms"],
                   n_atoms_ccd=G.number_of_nodes(), ca_rmsd=bmin["ca"], n_ca=bmin["n_ca"],
                   assign=json.dumps(bmin["asg"]), pred_copy=bmin["pi"],
                   ref_copy=f'{rlig[bmin["ri"]]["chain"]}:{rlig[bmin["ri"]]["seq"]}',
                   lig_rmsd_cabest=bca["lig"], ca_rmsd_cabest=bca["ca"], assign_cabest=json.dumps(bca["asg"]))
        rec["success"] = bool(rec["lig_rmsd"] <= SUCCESS_CUT)
        # clash + duplicates on predicted model
        pprot = np.array([a for ch in pres.values() for _, atoms, _ in ch for a in atoms.values()])
        dmins = [float(np.min(np.linalg.norm(pprot[None] - np.array(list(c["atoms"].values()))[:, None], axis=2)))
                 for c in plig]
        rec["min_lig_prot_dist"] = min(dmins)
        rec["clash"] = bool(min(dmins) < CLASH_CUT)
        cents = [np.array(list(c["atoms"].values())).mean(0) for c in plig]
        dup = [float(np.linalg.norm(a - b)) for a, b in itertools.combinations(cents, 2)]
        rec["min_copy_centroid_dist"] = min(dup) if dup else np.nan
        rec["dup_ligand"] = bool(dup and min(dup) < DUP_CUT)
        rec["n_atoms_pred_vs_ccd"] = f'{min(len(c["atoms"]) for c in plig)}/{G.number_of_nodes()}'
    except Exception as e:  # report, do not hide
        rec["status"] = f"error:{type(e).__name__}:{e}"
    return rec


def site_share(args):
    batch, row = args
    rec = dict(batch=batch, case_id=row.case_id, pdb_id=row.pdb_id, ccd=row.ccd)
    ref = load_ref(row.pdb_id)
    rnames = ref_chain_names(batch, row)
    cen = centroid_of(row)
    copies = ligand_copies(ref, row.ccd)
    d = [np.linalg.norm(np.array(list(c["atoms"].values())).mean(0) - cen) for c in copies]
    k = int(np.argmin(d))
    rec["ref_copy"] = f'{copies[k]["chain"]}:{copies[k]["seq"]}'
    rec["ref_copy_centroid_offset"] = float(d[k])
    site = biolip_site(copies[k], {n: ref.find_chain(n) for n in rnames})
    v = sorted(site.values())
    rec["site_chain_counts"] = json.dumps(site)
    rec["n_site"] = sum(v)
    rec["n_minor"] = v[0]
    rec["partner_share"] = v[0] / sum(v) if sum(v) else np.nan
    return rec


def jobs_for(batch, seeds):
    cfg = BATCHES[batch]
    cases = pd.read_csv(ROOT / cfg["cases"], sep="\t")
    out = []
    for _, row in cases.iterrows():
        for cond in ("dimer", "monomer"):
            for seed in seeds:
                base = f"{row.case_id}_{cond}"
                names = [base] if seed == 0 else [f"{base}_seed{seed}"]
                for d in SEED_DIRS:
                    for nm in names:
                        jd = ROOT / d.format(batch=batch) / nm
                        if (jd / f"{nm}_ranking_scores.csv").exists():
                            out.append((batch, row, cond, seed, jd, nm))
                            break
                    else:
                        continue
                    break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=str(OUT / "rescored_seed0.tsv"))
    ap.add_argument("--sites", action="store_true")
    a = ap.parse_args()
    seeds = [int(s) for s in a.seeds.split(",")]
    # warm CCD cache serially (network)
    for b in BATCHES:
        for c in pd.read_csv(ROOT / BATCHES[b]["cases"], sep="\t").ccd.unique():
            ccd_graph(c)
    jobs = [j for b in BATCHES for j in jobs_for(b, seeds)]
    print(f"{len(jobs)} jobs", file=sys.stderr)
    with ProcessPoolExecutor(a.workers) as ex:
        recs = list(ex.map(score_job, jobs, chunksize=2))
    pd.DataFrame(recs).to_csv(a.out, sep="\t", index=False, float_format="%.3f")
    if a.sites:
        rows = [(b, r) for b in BATCHES for _, r in pd.read_csv(ROOT / BATCHES[b]["cases"], sep="\t").iterrows()]
        with ProcessPoolExecutor(a.workers) as ex:
            srec = list(ex.map(site_share, rows, chunksize=2))
        pd.DataFrame(srec).to_csv(OUT / "site_share.tsv", sep="\t", index=False, float_format="%.3f")


if __name__ == "__main__":
    main()
