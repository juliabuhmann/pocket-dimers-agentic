#!/usr/bin/env python
"""H4 / Q2.1 step 2-3, 5: interface correspondence + pocket recovery per benchmark site (first attempt, NOT verified).

Per site x model (status ok in sites.parquet) and per chain assignment (homodimer: both, hetero: fixed by accession):
  * local superposition: PDB CA atoms of the chain pair within SUPERPOSE_RADIUS of the ligand centroid, matched to
    the model by UniProt position (fallback: whole mapped chain pair, flagged); ligand centroid transformed;
  * interface Jaccard (UniProt numbering, chain-role aware) between the PDB interface and the model interface (5 Å);
  * DCC from the transformed centroid to the centres of the top-(n+2) pockets by burrow_score, for
      dimer      = dimer-model pockets (homodimers: symmetric copies collapsed, is_sym_rep),
      iso        = pooled pockets of the two isolated chains (same coordinates; homodimers: same collapse rule),
      *_nocollapse, *_atlas (atlas filter burrow_score>=0.03 & lining pLDDT>=70) = sensitivities;
  * residue metrics (best-F1 pocket of the top-(n+2)): recall, overlap coefficient, F1, pocket size;
  * null: size-matched random surface patches in the same model/context (seeds 0-2);
  * negative control: ligand centroid randomly rotated about the model centre of mass (seed 0).
Max over chain assignments is taken per metric. Secondary (step 5): residue metrics of the AFDB monomer atlas pockets
of the same proteins (pooled, top-(n+2)).
Usage: python scripts/h4_recovery.py --census v1 [--workers 6]
"""
import argparse
import glob
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from h4_common import (WS, DIMER_ROOT, KINDS, POCKETS_DIMER, POCKETS_ISO, ATLAS, OUT_DATA, DCC_HIT, TOPK_EXTRA,
                       IFACE_DIST, SUPERPOSE_RADIUS, SUPERPOSE_MIN_CA, SURFACE_NB_RADIUS, SURFACE_NB_MAX, SEEDS,
                       kabsch, parse_res, load_candidates)
from annotate_pocket_chains import sym_groups, ATLAS_SCORE, ATLAS_PLDDT

CONDS = ["dimer", "iso", "dimer_nocollapse", "iso_nocollapse", "dimer_occollapse", "iso_occollapse", "dimer_atlas", "iso_atlas"]
SYM_JACCARD = 0.5  # TOOLS.md §2: burrow dimer benchmark collapsed symmetric copies at residue-set Jaccard >= 0.5


def jaccard_sym_rep(df):
    """Homodimer symmetric-copy collapse: pockets i, j are copies if Jaccard(swap_AB(res_i), res_j) >= SYM_JACCARD
    (union-find); the highest-scoring pocket of each group is kept. Jaccard (not the overlap coefficient of
    scripts/annotate_pocket_chains.py) because the overlap coefficient lets a large interface pocket absorb distinct
    single-chain pockets (seen for GSTP1, AF-0000000065928178 m0 absorbing m1/m2)."""
    sets = list(df.rs)
    sw = [{({"A": "B", "B": "A"}.get(c, c), p) for c, p in x} for x in sets]
    n = len(sets); parent = list(range(n))
    def f(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            u = sw[i] | sets[j]
            if u and len(sw[i] & sets[j]) / len(u) >= SYM_JACCARD:
                parent[f(i)] = f(j)
    sc = df.burrow_score.values
    best = {}
    for i in range(n):
        r = f(i)
        if r not in best or sc[i] > sc[best[r]]:
            best[r] = i
    rep = np.zeros(n, bool); rep[list(best.values())] = True
    return rep


def pocket_sets(dim, iso, kind):
    """Return dict cond -> pockets (centres, residue sets, score) sorted by burrow_score desc."""
    out = {}
    dim = dim.copy().reset_index(drop=True); iso = iso.copy().reset_index(drop=True)
    dim["rs"] = dim.residues.map(parse_res); iso["rs"] = iso.residues.map(parse_res)
    if kind == "homo":
        dim["rep_j"] = jaccard_sym_rep(dim) if len(dim) else []
        iso["rep_j"] = jaccard_sym_rep(iso) if len(iso) else []
        iso["is_sym_rep"] = sym_groups(iso)[1] if len(iso) else []
    else:
        dim["rep_j"] = True; iso["rep_j"] = True; iso["is_sym_rep"] = True
    iso["atlas_pass"] = (iso.burrow_score >= ATLAS_SCORE) & (iso.plddt_mean >= ATLAS_PLDDT)
    out["dimer_nocollapse"] = dim
    out["iso_nocollapse"] = iso
    out["dimer"] = dim[dim.rep_j.astype(bool)]
    out["iso"] = iso[iso.rep_j.astype(bool)]
    out["dimer_occollapse"] = dim[dim.is_sym_rep.astype(bool)]
    out["iso_occollapse"] = iso[iso.is_sym_rep.astype(bool)]
    out["dimer_atlas"] = dim[dim.rep_j.astype(bool) & dim.atlas_pass.astype(bool)]
    out["iso_atlas"] = iso[iso.rep_j.astype(bool) & iso.atlas_pass]
    for k, v in out.items():
        v = v.sort_values("burrow_score", ascending=False, kind="stable")
        out[k] = dict(ctr=v[["center_x", "center_y", "center_z"]].to_numpy(float), rs=list(v.rs), score=v.burrow_score.to_numpy(),
                      vol=v.volume_A3.to_numpy())
    return out


def res_metrics(L, rs_list):
    best = (0.0, 0.0, 0.0, 0)
    for P in rs_list:
        i = len(L & P)
        if not i:
            continue
        rec, prec = i / len(L), i / len(P)
        f1 = 2 * rec * prec / (rec + prec)
        oc = i / min(len(L), len(P))
        if f1 > best[2]:
            best = (rec, oc, f1, len(P))
    return best


def read_model(path):
    import gemmi
    st = gemmi.read_structure(str(path))
    ca, cb, heavy, hres = {}, {}, [], []
    for ch in st[0]:
        for r in ch:
            p = r.seqid.num
            a = r.find_atom("CA", "*")
            if a is None:
                continue
            ca[(ch.name, p)] = np.array(a.pos.tolist())
            b = r.find_atom("CB", "*") or a
            cb[(ch.name, p)] = np.array(b.pos.tolist())
            for at in r:
                if at.element.name != "H":
                    heavy.append(at.pos.tolist()); hres.append((ch.name, p))
    return ca, cb, np.array(heavy), hres


def model_iface(heavy, hres):
    from scipy.spatial import cKDTree
    ch = np.array([c for c, _ in hres])
    out = set()
    A, B = np.where(ch == "A")[0], np.where(ch == "B")[0]
    if not len(A) or not len(B):
        return out
    tb = cKDTree(heavy[B])
    hits = tb.query_ball_point(heavy[A], IFACE_DIST)
    for i, h in zip(A, hits):
        if h:
            out.add(hres[i])
            for j in h:
                out.add(hres[B[j]])
    return out


def surface_context(cb):
    """Surface residues and CB arrays for dimer context (both chains) and isolated context (per chain)."""
    from scipy.spatial import cKDTree
    keys = list(cb)
    X = np.array([cb[k] for k in keys])
    chn = np.array([k[0] for k in keys])
    t = cKDTree(X)
    nb_dim = np.array([len(x) - 1 for x in t.query_ball_point(X, SURFACE_NB_RADIUS)])
    nb_iso = np.zeros(len(keys), int)
    for c in np.unique(chn):
        m = np.where(chn == c)[0]
        tc = cKDTree(X[m])
        nb_iso[m] = [len(x) - 1 for x in tc.query_ball_point(X[m], SURFACE_NB_RADIUS)]
    return dict(keys=keys, X=X, chn=chn, surf_dim=np.where(nb_dim < SURFACE_NB_MAX)[0], surf_iso=np.where(nb_iso < SURFACE_NB_MAX)[0])


def random_patches(ctx, sizes, rng, iso):
    """One random patch per size: seed = random surface residue, patch = seed + nearest residues (same chain if iso)."""
    surf = ctx["surf_iso"] if iso else ctx["surf_dim"]
    X, chn, keys = ctx["X"], ctx["chn"], ctx["keys"]
    ctrs, sets = [], []
    for s in sizes:
        i = rng.choice(surf)
        m = np.where(chn == chn[i])[0] if iso else np.arange(len(X))
        d = np.linalg.norm(X[m] - X[i], axis=1)
        sel = m[np.argsort(d)[: max(int(s), 1)]]
        ctrs.append(X[sel].mean(0)); sets.append({keys[j] for j in sel})
    return np.array(ctrs), sets


def random_rotation(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def process_model(args):
    model, kind, sites, cas, dim, iso = args
    try:
        return _process_model(model, kind, sites, cas, dim, iso)
    except Exception as ex:
        return [dict(model=model, row_id=s["row_id"], error=f"{type(ex).__name__}: {ex}") for s in sites]


def _process_model(model, kind, sites, cas, dim, iso):
    ca, cb, heavy, hres = read_model(DIMER_ROOT / KINDS[kind] / f"{model}.pdb")
    miface = model_iface(heavy, hres)
    ctx = surface_context(cb)
    com = heavy.mean(0)
    pk = pocket_sets(dim, iso, kind)
    out = []
    for s in sites:
        k = int(s["n_sites_pair"]) + TOPK_EXTRA
        cen = np.asarray(s["centroid"], float)
        assigns = [("A", "B"), ("B", "A")] if kind == "homo" else ([("A", "B")] if s["acc1"] == s["model_acc_A"] else [("B", "A")])
        best = defaultdict(lambda: None)
        recs = []
        for ai, (m1, m2) in enumerate(assigns):
            P, Q, near = [], [], []
            for role, mch in (("1", m1), ("2", m2)):
                pos, xyz = cas[(s["pdb_id"], s[f"chain{role}"], s[f"acc{role}"])]
                for p, x in zip(pos, xyz):
                    y = ca.get((mch, p))
                    if y is not None:
                        P.append(x); Q.append(y); near.append(np.linalg.norm(x - cen) <= SUPERPOSE_RADIUS)
            P, Q, near = np.array(P), np.array(Q), np.array(near, bool)
            if len(P) < 3:
                recs.append(dict(assign=ai, superpose="failed")); continue
            mode = "local" if near.sum() >= SUPERPOSE_MIN_CA else "pair_fallback"
            sel = near if mode == "local" else np.ones(len(P), bool)
            R, t, rmsd = kabsch(P[sel], Q[sel])
            x = cen @ R.T + t
            pdb_if = {(m1, p) for p in s["pdb_iface1"]} | {(m2, p) for p in s["pdb_iface2"]}
            jac = len(pdb_if & miface) / len(pdb_if | miface) if (pdb_if | miface) else 0.0
            L = {(m1, p) for p in s["bind1"]} | {(m2, p) for p in s["bind2"]}
            rec = dict(assign=ai, superpose=mode, n_ca_sup=int(sel.sum()), rmsd_sup=rmsd, iface_jaccard=jac, n_lig_res=len(L))
            rng_rot = np.random.default_rng([0, int(s["row_id"]), ai])
            xs = (x - com) @ random_rotation(rng_rot).T + com
            for c in CONDS:
                ctr = pk[c]["ctr"][:k]
                rs = pk[c]["rs"][:k]
                d = np.linalg.norm(ctr - x, axis=1) if len(ctr) else np.array([np.inf])
                rec[f"{c}_dcc"] = float(d.min())
                rec[f"{c}_hit"] = bool(d.min() <= DCC_HIT)
                rec[f"{c}_hit_rank"] = int(np.argmin(d)) + 1 if d.min() <= DCC_HIT else -1
                rec[f"{c}_npk"] = len(ctr)
                r_, o_, f_, n_ = res_metrics(L, rs) if L else (np.nan,) * 4
                rec.update({f"{c}_recall": r_, f"{c}_oc": o_, f"{c}_f1": f_, f"{c}_pksize": n_})
                if c in ("dimer", "iso"):
                    ds = np.linalg.norm(ctr - xs, axis=1) if len(ctr) else np.array([np.inf])
                    rec[f"{c}_shuf_hit"] = bool(ds.min() <= DCC_HIT)
                    sizes = [len(r) for r in rs]
                    for sd in SEEDS:
                        rng = np.random.default_rng([sd, int(s["row_id"]), ai, CONDS.index(c)])
                        if not sizes:
                            rec[f"{c}_null{sd}_hit"] = False; rec[f"{c}_null{sd}_f1"] = 0.0; continue
                        nc, ns = random_patches(ctx, sizes, rng, iso=(c == "iso"))
                        rec[f"{c}_null{sd}_hit"] = bool(np.linalg.norm(nc - x, axis=1).min() <= DCC_HIT)
                        rec[f"{c}_null{sd}_f1"] = res_metrics(L, ns)[2] if L else np.nan
            recs.append(rec)
        ok = [r for r in recs if r.get("superpose") != "failed"]
        row = dict(row_id=s["row_id"], model=model, n_assign_ok=len(ok))
        if ok:
            row["iface_jaccard"] = max(r["iface_jaccard"] for r in ok)
            row["superpose"] = ",".join(sorted({r["superpose"] for r in ok}))
            row["rmsd_sup"] = min(r["rmsd_sup"] for r in ok)
            row["n_ca_sup"] = max(r["n_ca_sup"] for r in ok)
            row["n_lig_res"] = ok[0]["n_lig_res"]
            for key in ok[0]:
                if key in ("assign", "superpose", "n_ca_sup", "rmsd_sup", "iface_jaccard", "n_lig_res"):
                    continue
                vals = [r[key] for r in ok]
                if key.endswith("_dcc"):
                    row[key] = min(vals)
                elif key.endswith("_hit_rank"):
                    pos = [v for v in vals if v > 0]
                    row[key] = min(pos) if pos else -1
                elif key.endswith("_hit"):
                    row[key] = any(vals)
                elif key.endswith("_f1"):
                    row[key] = np.nanmax(vals) if not all(pd.isna(vals)) else np.nan
                elif key.endswith(("_recall", "_oc", "_pksize", "_npk")):
                    # take the values from the assignment with the best F1 of the same condition
                    c = key.rsplit("_", 1)[0]
                    f1s = [r.get(f"{c}_f1", -1) for r in ok]
                    j = int(np.nanargmax([(-1 if pd.isna(v) else v) for v in f1s])) if key[-4:] != "_npk" else 0
                    row[key] = ok[j][key]
                else:
                    row[key] = vals[0]
        out.append(row)
    return out


def load_iso():
    fs = sorted(glob.glob(str(POCKETS_ISO / "*" / "chunk_[0-9][0-9][0-9][0-9].parquet")))
    cols = ["structure_stem", "burrow_score", "center_x", "center_y", "center_z", "residues", "volume_A3", "plddt_mean"]
    d = pd.concat([pd.read_parquet(f, columns=cols) for f in fs], ignore_index=True)
    d["model"] = d.structure_stem.str.rsplit("_", n=1).str[0]
    return d


def atlas_metrics(sites, census):
    """Step 5 (secondary): residue metrics of AFDB monomer atlas pockets (tier 1, fragment 1) of the same proteins."""
    accs = sorted(set(sites.acc1) | set(sites.acc2))
    fs = sorted(glob.glob(str(ATLAS / "pockets/tier1/*.parquet")))
    cols = ["accession", "fragment", "burrow_score", "residues"]
    parts = []
    for f in fs:
        t = pq.read_table(f, columns=cols, filters=[("accession", "in", accs), ("fragment", "==", 1)])
        if t.num_rows:
            parts.append(t.to_pandas())
    at = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
    at["pos"] = at.residues.map(lambda s: {p for _, p in parse_res(s)})
    by = {a: g.sort_values("burrow_score", ascending=False) for a, g in at.groupby("accession")}
    rows = []
    for s in sites.itertuples(index=False):
        k = int(s.n_sites_pair) + TOPK_EXTRA
        prots = [s.acc1] if s.acc1 == s.acc2 else [s.acc1, s.acc2]
        pool = pd.concat([by[a].assign(acc=a) for a in prots if a in by], ignore_index=True) if any(a in by for a in prots) else None
        r = dict(row_id=s.row_id, atlas_has_model=pool is not None)
        if pool is None:
            rows.append(r); continue
        pool = pool.sort_values("burrow_score", ascending=False).head(k)
        L1, L2 = set(s.bind1), set(s.bind2)
        # chain-role sets: a monomer pocket of acc covers the role(s) with that accession
        Ls = {("1", p) for p in L1} | {("2", p) for p in L2}
        cand = []
        for a, ps in zip(pool.acc, pool.pos):
            for role, acc in (("1", s.acc1), ("2", s.acc2)):
                if acc == a:
                    cand.append({(role, p) for p in ps})
        r["atlas_recall"], r["atlas_oc"], r["atlas_f1"], r["atlas_pksize"] = res_metrics(Ls, cand) if Ls else (np.nan,) * 4
        r["atlas_npk"] = len(pool)
        rows.append(r)
    return pd.DataFrame(rows), len(at), at.accession.nunique()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, choices=["v1", "v2"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--models", nargs="*")
    a = ap.parse_args()
    od = OUT_DATA / a.census
    sites = pd.read_parquet(od / "sites.parquet")
    sites = sites[sites.status == "ok"].reset_index(drop=True)
    sites["row_id"] = np.arange(len(sites))
    # n = number of ligand sites (collapsed site groups, same ligand set) of this chain pair in this entry
    sites["pair_key"] = sites.pdb_id + "|" + sites.model + "|" + sites.apply(lambda r: "+".join(sorted([r.chain1, r.chain2])), axis=1)
    sites["n_sites_pair"] = sites.groupby(["pair_key", "lig_set"]).site_group.transform("nunique")
    cand = load_candidates().set_index("modelEntityId")
    sites["model_acc_A"] = sites.model.map(cand.uniprotAccession_1)
    sites["model_acc_B"] = sites.model.map(cand.uniprotAccession_2)
    for k in ("ipTM", "max_ipSAE", "max_pDockQ2"):
        sites[k] = sites.model.map(cand[k])
    if a.models:
        sites = sites[sites.model.isin(a.models)]
    ca = pd.read_parquet(od / "chain_ca.parquet")
    cas = {(p, c, acc): (np.asarray(pos, int), np.asarray(x, float).reshape(-1, 3)) for p, c, acc, pos, x in
           zip(ca.pdb_id, ca.chain, ca.acc, ca.pos, ca.xyz)}
    dim = pd.read_parquet(POCKETS_DIMER, columns=["structure_stem", "burrow_score", "center_x", "center_y", "center_z", "residues",
                                                  "volume_A3", "is_sym_rep", "atlas_pass"])
    dim = dim[dim.structure_stem.isin(set(sites.model))]
    iso = load_iso()
    iso = iso[iso.model.isin(set(sites.model))]
    dg = dict(list(dim.groupby("structure_stem"))); ig = dict(list(iso.groupby("model")))
    empty_d, empty_i = dim.iloc[:0], iso.iloc[:0]
    tasks = []
    missing = 0
    for (m, kind), g in sites.groupby(["model", "kind"]):
        if m not in dg:
            missing += len(g)
        recs = g.to_dict("records")
        sub = {k: v for k, v in cas.items() if k[0] in set(g.pdb_id)}
        tasks.append((m, kind, recs, sub, dg.get(m, empty_d), ig.get(m, empty_i)))
    print(f"{len(sites)} site rows, {len(tasks)} models, rows on models without dimer pockets: {missing}", flush=True)
    res = []
    t = time.time()
    with ProcessPoolExecutor(max_workers=min(a.workers, 6)) as ex:
        for i, r in enumerate(ex.map(process_model, tasks, chunksize=2)):
            res += r
            if (i + 1) % 200 == 0:
                print(f"{i+1}/{len(tasks)} models {time.time()-t:.0f}s", flush=True)
    rec = pd.DataFrame(res)
    out = sites.merge(rec, on=["row_id", "model"], how="left")
    am, n_at, n_acc = atlas_metrics(out, a.census)
    out = out.merge(am, on="row_id", how="left")
    print(f"atlas pockets loaded {n_at} for {n_acc} accessions", flush=True)
    out["model_has_dimer_pockets"] = out.model.isin(set(dg))
    out.to_parquet(od / ("recovery.parquet" if not a.models else "recovery_test.parquet"))
    print("errors:", out["error"].notna().sum() if "error" in out else 0)


if __name__ == "__main__":
    main()
