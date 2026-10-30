"""Independent verifier: recompute DCC / hit for sampled H4 sites from raw coordinates.
Own SIFTS-XML parsing, own assembly-1 chain construction (gemmi), own Kabsch superposition, own pocket top-k + symmetric
collapse. Does not import scripts/h4_*.py or src/part1. Status: verifier's first attempt, UNVERIFIED.
Inputs: site table (producer), raw mmCIF mirror, SIFTS XML cache, AFDB dimer models, burrow pocket tables."""
import gzip, glob, sys, json, re
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd, gemmi

WS = Path(__file__).resolve().parents[2]
OUT = WS / "reviews/verify_h4"
SITE_TABLE = WS / "data/part2/h4/v2/benchmark_filtered.parquet"
PDB_MIRROR = Path("/cluster/scratch/jjaenes/pdb")
SIFTS_XML = WS / "data/raw/sifts_xml"
POCKETS = {"dimer": WS / "data/part2/pockets/human", "iso": WS / "data/part2/pockets/human_isolated"}
DCC_CUT = 4.0          # DESIGN Q2.1 primary metric
TOPK_EXTRA = 2         # DESIGN: top-(n+2)
LOCAL_R = 15.0         # producer's documented local-frame radius (used only for the 'local' variant)
SYM_JACC = 0.5         # producer's documented collapse rule (swapped-residue Jaccard)
SEED = 0
N_WORKERS = 6

def parse_sifts(pdb):
    """(auth_chain, label_seq) -> (uniprot acc, uniprot pos, one-letter)"""
    m = {}
    with gzip.open(SIFTS_XML / f"{pdb}.xml.gz") as fh:
        for _, el in ET.iterparse(fh):
            if el.tag.endswith("residue") and el.get("dbSource") == "PDBe":
                ls = int(el.get("dbResNum")); ch = None; up = None
                for x in el:
                    if not x.tag.endswith("crossRefDb"): continue
                    if x.get("dbSource") == "PDB": ch = x.get("dbChainId")
                    elif x.get("dbSource") == "UniProt" and up is None and x.get("dbResNum") not in (None, "null"):
                        up = (x.get("dbAccessionId"), int(x.get("dbResNum")), x.get("dbResName"))
                if ch is not None and up is not None: m[(ch, ls)] = up
                el.clear()
    return m

def build_chains(pdb):
    """assembly-1 polymer chains named auth or auth-<oper> (non-identity); returns {cid: [(label_seq, resname, CA xyz)]}, ligand atoms"""
    st = gemmi.read_structure(str(PDB_MIRROR / pdb[1:3] / f"{pdb}.cif.gz"))
    st.setup_entities()
    model = st[0]
    chains = {}
    gens = []
    if len(st.assemblies):
        for g in st.assemblies[0].generators:
            gens.append((set(g.subchains), set(g.chains), [(o.name, o.transform) for o in g.operators]))
    else:
        gens = [(None, None, [("1", gemmi.Transform())])]
    seen = set()
    for subs, auths, ops in gens:
        for oname, tr in ops:
            ident = np.allclose(np.array(tr.mat.tolist()), np.eye(3)) and np.allclose(tr.vec.tolist(), 0)
            for ch in model:
                # iterate polymer subchain of this chain
                poly = ch.get_polymer()
                if len(poly) == 0: continue
                sub = poly.subchain_id()
                if subs is not None and sub not in subs and not (auths and ch.name in auths and not subs): continue
                cid = ch.name if ident else f"{ch.name}-{oname}"
                if cid in seen: continue
                seen.add(cid)
                res = []
                for r in poly:
                    ca = r.find_atom("CA", "*")
                    if ca is None or r.label_seq is None: continue
                    p = tr.apply(ca.pos) if not ident else ca.pos
                    res.append((r.label_seq, r.name, np.array([p.x, p.y, p.z])))
                chains[cid] = res
    return st, chains

def ligand_centroid(st, uid):
    f = uid.split(":"); ch, num, ccd = f[1], f[2], f[3]
    alt = f[5][3:] if len(f) > 5 and f[5].startswith("alt") else None
    xyz = []
    for c in st[0]:   # uid field 2 is label_asym_id (subchain), field 3 auth seq num
        for r in c:
            if r.subchain == ch and r.name == ccd and str(r.seqid.num) + (r.seqid.icode.strip() or "") == num:
                xyz += [[a.pos.x, a.pos.y, a.pos.z] for a in r if a.element.name != "H"
                        and (a.altloc == "\0" or alt is None and a.altloc == "A" or a.altloc == alt)]
    return np.array(xyz).mean(0) if xyz else None

def read_model_ca(path):
    st = gemmi.read_structure(str(path)); out = {}
    for ch in st[0]:
        seq = {}
        for r in ch:
            ca = r.find_atom("CA", "*")
            if ca is not None: seq[r.seqid.num] = (gemmi.find_tabulated_residue(r.name).one_letter_code.upper(), np.array([ca.pos.x, ca.pos.y, ca.pos.z]))
        out[ch.name] = seq
    return out

def kabsch(P, Q):
    """returns R, t minimising |R P + t - Q|"""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc

def parse_res(s):
    return frozenset(tuple(x.split(":")) for x in s.split()) if isinstance(s, str) and s else frozenset()

def collapse(pk, homo):
    """symmetric-copy collapse for homodimers: swapped-residue Jaccard >= SYM_JACC, union-find, keep highest burrow_score"""
    pk = pk.sort_values("burrow_score", ascending=False).reset_index(drop=True)
    if not homo or len(pk) < 2: return pk
    sets = [parse_res(s) for s in pk.residues]
    swap = [frozenset((("B" if c == "A" else "A" if c == "B" else c), n) for c, n in s) for s in sets]
    par = list(range(len(pk)))
    def f(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for i in range(len(pk)):
        for j in range(i + 1, len(pk)):
            u = len(swap[i] | sets[j])
            if u and len(swap[i] & sets[j]) / u >= SYM_JACC: par[f(j)] = f(i)
    keep = [i for i in range(len(pk)) if f(i) == i]   # root is the lower index = higher score (union always to i<j root)
    return pk.iloc[sorted(set(f(i) for i in range(len(pk))))].reset_index(drop=True)

def score_site(row, pk_dimer, pk_iso):
    r = dict(row)
    pdb = row["pdb_id"]
    sm = parse_sifts(pdb)
    st, chains = build_chains(pdb)
    cen = ligand_centroid(st, row["ligand_uid"])
    r["v_centroid_shift"] = float(np.linalg.norm(cen - np.asarray(row["centroid"]))) if cen is not None else np.nan
    if cen is None: cen = np.asarray(row["centroid"], float)
    mdl = read_model_ca(row["model_path"])
    homo = row["kind"] == "homo"
    accs = {"A": row["model_acc_A"], "B": row["model_acc_B"]}
    assigns = [("A", "B"), ("B", "A")] if homo else [("A", "B") if accs["A"] == row["acc1"] else ("B", "A")]
    best = {}
    for mc1, mc2 in assigns:
        P, Q, nid, ntot = [], [], 0, 0
        for pc, mc, acc in [(row["chain1"], mc1, row["acc1"]), (row["chain2"], mc2, row["acc2"])]:
            auth = pc.split("-")[0].rstrip("'")
            for ls, rn, xyz in chains.get(pc, []):
                u = sm.get((auth, ls))
                if u is None or u[0] != acc or u[1] not in mdl[mc]: continue
                ntot += 1; nid += (mdl[mc][u[1]][0] == u[2])
                P.append(xyz); Q.append(mdl[mc][u[1]][1])
        P, Q = np.array(P), np.array(Q)
        if len(P) < 10: continue
        out = {}
        for frame in ("global", "local"):
            if frame == "local":
                sel = np.linalg.norm(P - cen, axis=1) <= LOCAL_R
                if sel.sum() < 10: sel = np.ones(len(P), bool)
            else:
                sel = np.ones(len(P), bool)
            R, t = kabsch(P[sel], Q[sel])
            rmsd = float(np.sqrt(((P[sel] @ R.T + t - Q[sel]) ** 2).sum(1).mean()))
            c = R @ cen + t
            o = dict(rmsd=rmsd, nca=int(sel.sum()), seqid=nid / max(ntot, 1))
            for lab, pk in (("dimer", pk_dimer), ("iso", pk_iso)):
                k = int(row["n_sites_pair"]) + TOPK_EXTRA
                cent = pk[["center_x", "center_y", "center_z"]].values
                dcc = np.linalg.norm(cent - c, axis=1) if len(pk) else np.array([])
                topd = dcc[:k]
                o[f"{lab}_dcc_topk"] = float(topd.min()) if len(topd) else np.nan
                o[f"{lab}_rank_topk"] = int(np.argmin(topd)) + 1 if len(topd) else -1
                o[f"{lab}_dcc_any"] = float(dcc.min()) if len(dcc) else np.nan
                o[f"{lab}_rank_any"] = int(np.argmin(dcc)) + 1 if len(dcc) else -1
                o[f"{lab}_npk_all"] = len(pk); o[f"{lab}_k"] = min(k, len(pk))
                o[f"{lab}_hit"] = bool(len(topd) and topd.min() <= DCC_CUT)
                if len(dcc):
                    j = int(np.argmin(dcc)); res = parse_res(pk.residues.iloc[j])
                    o[f"{lab}_nearest_size"] = len(res); o[f"{lab}_nearest_chains"] = "".join(sorted({x[0] for x in res}))
                    o[f"{lab}_nearest_vol"] = float(pk.volume_A3.iloc[j]); o[f"{lab}_nearest_score_rank"] = j + 1
                    # PDB binding residues (UniProt numbering) in assignment's model chains; compare with nearest pocket
                    bset = {(mc1, int(x)) for x in row["bind1"]} | {(mc2, int(x)) for x in row["bind2"]}
                    o[f"{lab}_nearest_recall"] = len(bset & {(a, int(b)) for a, b in res}) / max(len(bset), 1)
                    # per-pocket diagnostics (score rank, size, chains, recall of PDB binding residues, DCC)
                    pl = []
                    for jj in range(len(pk)):
                        rs = {(a, int(b)) for a, b in parse_res(pk.residues.iloc[jj])}
                        pl.append([jj + 1, round(float(pk.burrow_score.iloc[jj]), 3), len(rs), "".join(sorted({a for a, _ in rs})),
                                   round(len(bset & rs) / max(len(bset), 1), 2), round(float(dcc[jj]), 1)])
                    o[f"{lab}_pockets"] = json.dumps(pl); o[f"{lab}_nbind"] = len(bset)
            out[frame] = o
        best[(mc1, mc2)] = out
    if not best:
        r["v_status"] = "no_superposition"; return r
    r["v_status"] = "ok"
    for frame in ("global", "local"):
        # homodimer: max over assignments (hit if any assignment hits; report the assignment with min DCC)
        for lab in ("dimer", "iso"):
            ka = min(best, key=lambda a: best[a][frame][f"{lab}_dcc_topk"] if not np.isnan(best[a][frame][f"{lab}_dcc_topk"]) else 1e9)
            o = best[ka][frame]
            r[f"v_{frame}_{lab}_hit"] = any(best[a][frame][f"{lab}_hit"] for a in best)
            for key in ("dcc_topk", "rank_topk", "dcc_any", "rank_any", "npk_all", "k", "nearest_size", "nearest_chains",
                        "nearest_vol", "nearest_score_rank", "nearest_recall", "pockets", "nbind"):
                r[f"v_{frame}_{lab}_{key}"] = o.get(f"{lab}_{key}")
            r[f"v_{frame}_{lab}_assign"] = "".join(ka)
        r[f"v_{frame}_rmsd"] = best[ka][frame]["rmsd"]; r[f"v_{frame}_nca"] = best[ka][frame]["nca"]
        r[f"v_{frame}_seqid"] = best[ka][frame]["seqid"]
    return r

def load_pockets(stems, which):
    out = []
    for f in sorted(glob.glob(str(POCKETS[which] / "*" / "chunk_*[0-9].parquet"))):
        d = pd.read_parquet(f, columns=["pocket_id", "structure_stem", "burrow_score", "center_x", "center_y", "center_z", "residues", "volume_A3"])
        out.append(d[d.structure_stem.isin(stems)])
    return pd.concat(out)

def main():
    d = pd.read_parquet(SITE_TABLE)
    P = d[d.lig_set == "organic"].reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    picks = []
    # task 2: 10 MBS (dimer), 10 MBS (isolated), 10 SBS, disjoint, seed 0
    mbs = P.index[P.is_mbs].values; sbs = P.index[~P.is_mbs].values
    m20 = rng.choice(mbs, 20, replace=False); s10 = rng.choice(sbs, 10, replace=False)
    picks += [(i, "T2_MBS_dimer") for i in m20[:10]] + [(i, "T2_MBS_iso") for i in m20[10:]] + [(i, "T2_SBS") for i in s10]
    # task 3: 10 MBS dimer hits, 10 MBS dimer misses (seed 0, excluding already picked)
    used = set(i for i, _ in picks)
    h = np.setdiff1d(P.index[P.is_mbs & P.dimer_hit].values, list(used)); mi = np.setdiff1d(P.index[P.is_mbs & ~P.dimer_hit].values, list(used))
    picks += [(i, "T3_MBS_hit") for i in rng.choice(h, 10, replace=False)] + [(i, "T3_MBS_miss") for i in rng.choice(mi, 10, replace=False)]
    # NAMPT diagnostic: 5 NAMPT MBS sites (seed 0)
    nam = np.setdiff1d(P.index[P.is_mbs & (P.model_acc_A == "P43490")].values, [i for i, _ in picks])
    picks += [(i, "T4_NAMPT") for i in rng.choice(nam, 5, replace=False)]
    S = P.loc[[i for i, _ in picks]].copy(); S["sample"] = [s for _, s in picks]
    out_name = "dcc_recompute.tsv"
    if "--all-mbs" in sys.argv:   # full recomputation of every MBS site (for the failure-mode census)
        S = P.copy() if "--all" in sys.argv else P[P.is_mbs].copy()
        S["sample"] = "ALL"; out_name = "dcc_recompute_all.tsv"
    sub = {"homo": "homodimers", "hetero": "heterodimers"}
    S["model_path"] = [f"/cluster/scratch/gmagnusson/dimers/afdb_human_dimers/{sub[k]}/{m}.pdb" for k, m in zip(S.kind, S.model)]
    pd_d = load_pockets(set(S.model), "dimer")
    pd_i = load_pockets({f"{m}_{c}" for m in S.model for c in "AB"}, "iso")
    jobs = []
    for _, row in S.iterrows():
        homo = row.kind == "homo"
        pk_d = collapse(pd_d[pd_d.structure_stem == row.model], homo)
        pk_i = collapse(pd_i[pd_i.structure_stem.isin([f"{row.model}_A", f"{row.model}_B"])], homo)
        jobs.append((row.to_dict(), pk_d, pk_i))
    with ProcessPoolExecutor(N_WORKERS) as ex:
        res = list(ex.map(_wrap, jobs))
    R = pd.DataFrame(res)
    R.drop(columns=[c for c in R.columns if R[c].map(lambda x: isinstance(x, np.ndarray)).any()]).to_csv(OUT / out_name, sep="\t", index=False)
    print(R[["sample", "pdb_id", "ccd", "model", "v_status", "dimer_dcc", "v_local_dimer_dcc_topk", "v_global_dimer_dcc_topk", "dimer_hit", "v_local_dimer_hit", "v_global_dimer_hit",
             "iso_dcc", "v_local_iso_dcc_topk", "iso_hit", "v_local_iso_hit", "v_global_iso_hit", "v_centroid_shift", "v_local_rmsd", "v_local_seqid"]].round(2).to_string())

def _wrap(j):
    try:
        return score_site(*j)
    except Exception as e:
        r = dict(j[0]); r["v_status"] = f"error: {type(e).__name__}: {e}"; return r

if __name__ == "__main__":
    main()
