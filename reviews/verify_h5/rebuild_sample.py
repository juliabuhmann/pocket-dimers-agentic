"""Independent rebuild of H5 residue categories for a random protein sample (verifier code; not verified itself).
Run with .venvs/h5/bin/python (needs gemmi + freesasa). Does not import scripts/h5_*.py.
Categories (priority): 1 dimer-only pocket > 2 single-chain pocket > 3 interface (heavy atom <= INTERFACE_A to partner)
> 4 surface (rSASA_dimer >= RSA_SURF) > 5 buried. Pockets: atlas filter (burrow >= BURROW_MIN, pocket plddt_mean >= PLDDT_MIN),
interface pocket = >= IF_MIN_PER_CHAIN lining residues on each chain, dimer-only = max overlap coeff with any isolated-chain
pocket < OVERLAP_MAX. Lining residue kept only if its own pLDDT >= PLDDT_MIN.
"""
from pathlib import Path
import sys, os, glob, argparse
import numpy as np, pandas as pd, gemmi, freesasa
from concurrent.futures import ProcessPoolExecutor

ROOT = str(Path(__file__).resolve().parents[2])
DIMERS = "/path/to/afdb_dimers/afdb_human_dimers"
OUT = f"{ROOT}/reviews/verify_h5/out"
# DESIGN.md §0 / §1 H5 and docs/H5_RESULTS.md definitions
BURROW_MIN, PLDDT_MIN, IF_MIN_PER_CHAIN, OVERLAP_MAX = 0.03, 70.0, 3, 0.5
INTERFACE_A, RSA_SURF = 5.0, 0.2
N_SAMPLE, SEED = 300, 0
SYM_JACCARD = 0.5  # verifier choice for homodimer chain-swap collapse

def load_candidates():
    rows = []
    for kind, d in [("homo", "homodimers"), ("hetero", "heterodimers")]:
        c = pd.read_csv(f"{DIMERS}/{d}/candidates.csv")
        c["kind"] = kind; c["path"] = [f"{DIMERS}/{d}/{m}.pdb" for m in c.modelEntityId]
        rows.append(c)
    c = pd.concat(rows).drop_duplicates("modelEntityId")
    return c

def model_features(args):
    mid, path, accA, accB, keep = args
    st = gemmi.read_structure(path); st.remove_hydrogens(); st.setup_entities()
    m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    # SASA: dimer context and separate chains
    fs_d = freesasa.Structure(path)
    rd = freesasa.calc(fs_d).residueAreas()
    rel_iso = {}
    for s in freesasa.structureArray(path, {"separate-chains": True}):
        ra = freesasa.calc(s).residueAreas()
        for ch, dd in ra.items():
            for k, v in dd.items(): rel_iso[(ch, k)] = v.relativeTotal
    acc = {"A": accA, "B": accB}
    out = []
    for ch in m:
        if acc.get(ch.name) not in keep: continue
        for r in ch:
            ca = r.find_atom("CA", "*")
            pl = ca.b_iso if ca else np.mean([a.b_iso for a in r])
            mind = 99.0
            for a in r:
                for mk in ns.find_atoms(a.pos, "\0", radius=INTERFACE_A):
                    cra = mk.to_cra(m)
                    if cra.chain.name != ch.name:
                        d = cra.atom.pos.dist(a.pos)
                        if d < mind: mind = d
            k = str(r.seqid.num)
            ad = rd.get(ch.name, {}).get(k)
            out.append((mid, ch.name, acc[ch.name], r.seqid.num, r.name, pl,
                        ad.relativeTotal if ad else np.nan, rel_iso.get((ch.name, k), np.nan), mind))
    return out

def parse_res(s):
    return set(tuple((t.split(":")[0], int(t.split(":")[1]))) for t in s.split())

def overlap(a, b):
    return len(a & b) / min(len(a), len(b)) if a and b else 0.0

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    cand = load_candidates()
    accs = np.array(sorted(set(cand.uniprotAccession_1) | set(cand.uniprotAccession_2)))
    rng = np.random.default_rng(SEED)
    sample = set(rng.choice(accs, size=N_SAMPLE, replace=False))
    pd.Series(sorted(sample), name="accession").to_csv(f"{OUT}/sample_accessions.tsv", sep="\t", index=False)
    sub = cand[cand.uniprotAccession_1.isin(sample) | cand.uniprotAccession_2.isin(sample)]
    print(f"{len(accs)} proteins; sample {len(sample)}; models {len(sub)}", flush=True)
    jobs = [(r.modelEntityId, r.path, r.uniprotAccession_1, r.uniprotAccession_2, sample) for r in sub.itertuples()]
    feats = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, rows in enumerate(ex.map(model_features, jobs, chunksize=2)):
            feats.extend(rows)
            if i % 50 == 0: print("features", i, flush=True)
    F = pd.DataFrame(feats, columns=["model", "chain", "accession", "position", "resname", "plddt", "rsa_dimer", "rsa_iso", "min_dist"])
    F["interface"] = F.min_dist <= INTERFACE_A
    # pockets from raw describe output
    mids = set(sub.modelEntityId)
    dp = pd.concat([pd.read_parquet(f, columns=["pocket_id", "structure_stem", "residues", "burrow_score", "plddt_mean", "ok"])
                    for f in glob.glob(f"{ROOT}/data/part2/pockets/human/*/chunk_*[0-9].parquet")])
    dp = dp[dp.structure_stem.isin(mids)]
    ip = pd.concat([pd.read_parquet(f, columns=["pocket_id", "structure_stem", "residues"])
                    for f in glob.glob(f"{ROOT}/data/part2/pockets/human_isolated/*/chunk_*[0-9].parquet")])
    ip["model"] = ip.structure_stem.str.rsplit("_", n=1).str[0]
    ip = ip[ip.model.isin(mids)]
    iso = {m: [parse_res(s) for s in g.residues] for m, g in ip.groupby("model")}
    kind = dict(zip(sub.modelEntityId, sub.kind))
    plmap = {(r.model, r.chain, r.position): r.plddt for r in F.itertuples()}
    pk = []
    for r in dp.itertuples():
        res = parse_res(r.residues)
        nA = sum(1 for c, _ in res if c == "A"); nB = sum(1 for c, _ in res if c == "B")
        is_if = nA >= IF_MIN_PER_CHAIN and nB >= IF_MIN_PER_CHAIN
        ov = max([overlap(res, s) for s in iso.get(r.structure_stem, [])], default=0.0)
        pk.append(dict(pocket_id=r.pocket_id, model=r.structure_stem, res=res, burrow=r.burrow_score,
                       atlas=(r.burrow_score >= BURROW_MIN) and (r.plddt_mean >= PLDDT_MIN), is_if=is_if, ov=ov,
                       dimer_only=is_if and ov < OVERLAP_MAX, kind=kind[r.structure_stem]))
    PK = pd.DataFrame(pk)
    # homodimer symmetry collapse (verifier version): chain-swap Jaccard >= SYM_JACCARD -> group, rep = max burrow
    PK["sym_rep"] = True
    for m, g in PK[PK.kind == "homo"].groupby("model"):
        idx = list(g.index); parent = {i: i for i in idx}
        def f(i):
            while parent[i] != i: i = parent[i]
            return i
        sw = {i: {("B" if c == "A" else "A", p) for c, p in PK.at[i, "res"]} for i in idx}
        for x in idx:
            for y in idx:
                if x < y:
                    j = len(sw[x] & PK.at[y, "res"]) / max(1, len(sw[x] | PK.at[y, "res"]))
                    if j >= SYM_JACCARD: parent[f(x)] = f(y)
        comp = {}
        for i in idx: comp.setdefault(f(i), []).append(i)
        for members in comp.values():
            best = max(members, key=lambda i: PK.at[i, "burrow"])
            for i in members: PK.at[i, "sym_rep"] = (i == best)
    PK.drop(columns="res").to_parquet(f"{OUT}/sample_pockets.parquet")
    for variant, mask in [("nocollapse", PK.atlas), ("collapse", PK.atlas & PK.sym_rep)]:
        lab = {}
        for r in PK[mask].itertuples():
            code = 1 if r.dimer_only else 2
            for c, p in r.res:
                if plmap.get((r.model, c, p), 0) < PLDDT_MIN: continue
                k = (r.model, c, p); lab[k] = min(lab.get(k, 9), code)
        F[f"pk_{variant}"] = [lab.get((m, c, p), 9) for m, c, p in zip(F.model, F.chain, F.position)]
        cat = np.where(F[f"pk_{variant}"] < 9, F[f"pk_{variant}"], np.where(F.interface, 3, np.where(F.rsa_dimer >= RSA_SURF, 4, 5)))
        F[f"cat_{variant}"] = cat
    F.to_parquet(f"{OUT}/sample_model_residues.parquet")
    print("done", len(F), flush=True)

if __name__ == "__main__":
    main()
