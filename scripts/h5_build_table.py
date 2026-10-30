#!/usr/bin/env python
"""H5 step 1b/2: residue category table (per model-chain and aggregated per protein) + outcomes.

Inputs: data/part2/h5/residue_features/chunk_*.parquet (h5_residue_features.py),
        data/part2/human_dimer_pockets_annotated.parquet (annotate_pocket_chains.py),
        data/processed/variants_humsavar.parquet, data/processed/alphamissense_residue_mean.parquet,
        data/processed/sifts_chain_uniprot.parquet (PDB-complex flag, exploratory).
Categories (mutually exclusive, priority order; DESIGN.md v3 H5):
  1 dimer_only : lining residue (pLDDT >= LINING_PLDDT) of an atlas-pass pocket with >=3 lining residues on each
                 chain (is_interface) and max_overlap_iso < DIMER_ONLY_OVERLAP
  2 sc_pocket  : lining residue (pLDDT >= LINING_PLDDT) of any other atlas-pass pocket ("single-chain pocket";
                 includes interface pockets that also exist in the isolated chain, so that category 3 is
                 'interface residues outside any pocket' as the hypothesis states)
  3 interface  : heavy atom within 5 A of the partner chain (h5_residue_features.INTERFACE_CUTOFF)
  4 surface    : rsa_dimer >= SURFACE_RSA
  5 buried     : rest (excluded from models)
Homodimer pockets: only symmetry representatives (is_sym_rep) are used (DESIGN.md §0 collapse rule).
Per-protein aggregation (primary): residue takes its highest-priority category over all models and both chains
mapping to that protein; covariates (rsa_iso, plddt) = mean over the model-chain instances that carry that category.
Sensitivity: one model per protein = the model with the highest max_ipSAE (ties: modelEntityId order).
Status: first attempt, not independently verified.
"""
import glob, hashlib, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
FEAT = WS / "data/part2/h5/residue_features"
TAB = WS / "data/part2/h5/tables"
RES = WS / "results/part2/h5"
ATLAS_SCORE, ATLAS_PLDDT = 0.03, 70.0   # DESIGN.md §0 atlas filter (pocket level, as in annotate_pocket_chains.py)
LINING_PLDDT = 70.0                      # task spec: lining residues with pLDDT >= 70
DIMER_ONLY_OVERLAP = 0.5                 # DESIGN.md §0 (sensitivities 0.2, 0.8 -> exploratory columns)
SENS_OVERLAPS = (0.2, 0.8)
SURFACE_RSA = 0.2                        # task spec: other surface = relative SASA >= 0.2 (dimer context)
CATS = ["dimer_only", "sc_pocket", "interface", "surface", "buried"]
# PI-requested labelled sensitivities (added 2026-10-29 after Q2.2 found burrow differences away from the interface;
# requested after the primary H5 fit had been seen -> reported as SENSITIVITY/COMPARATOR, never replacing the primary):
BURROW_TOP_QUANTILE = 0.6   # (b) dimer-only restricted to pockets in the top two burrow_score quintiles (>= 60th pct of
                            #     atlas-pass symmetry-representative pockets); residues only in lower-quintile dimer-only
                            #     pockets are 'excluded' from that model
CATS_Q40 = ["dimer_only", "excluded", "sc_pocket", "interface", "surface", "buried"]
# comparator: interface pockets that DO exist in the isolated chain (is_interface & max_overlap_iso >= 0.5)
CATS_SPLIT = ["dimer_only", "ifpocket_preexisting", "sc_pocket", "interface", "surface", "buried"]
AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L",
           LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")


def provenance():
    sha = hashlib.sha256((WS / "DESIGN.md").read_bytes()).hexdigest()
    commit = subprocess.run(["git", "-C", str(WS), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {"design_sha256": sha, "git_commit": commit}


def load_features():
    fs = sorted(glob.glob(str(FEAT / "chunk_*.parquet")))
    f = pd.concat([pd.read_parquet(x) for x in fs], ignore_index=True)
    m = pd.read_csv(FEAT / "models.tsv", sep="\t")[["modelEntityId", "kind", "max_ipSAE", "ipTM", "max_pDockQ2"]]
    return f.merge(m, on="modelEntityId", how="left"), len(fs)


def pocket_residues(thr=DIMER_ONLY_OVERLAP):
    d = pd.read_parquet(WS / "data/part2/human_dimer_pockets_annotated.parquet")
    n0 = len(d)
    d = d[d.atlas_pass & d.is_sym_rep]
    qcut = float(d.burrow_score.quantile(BURROW_TOP_QUANTILE))
    rows = []
    for r in d.itertuples():
        do = bool(r.is_interface) and pd.notna(r.max_overlap_iso)
        flags = {t: do and r.max_overlap_iso < t for t in (thr,) + SENS_OVERLAPS}
        for tok in r.residues.split():
            ch, num = tok.split(":")
            rows.append((r.structure_stem, ch, int(num), flags[thr], flags[0.2], flags[0.8], bool(r.is_interface),
                         flags[thr] and r.burrow_score >= qcut, do and r.max_overlap_iso >= thr))
    cols = ["do", "do_s02", "do_s08", "pk_interface", "do_q40", "pre"]
    p = pd.DataFrame(rows, columns=["modelEntityId", "chain", "position"] + cols)
    p = p.groupby(["modelEntityId", "chain", "position"], as_index=False)[cols].max()
    return p, n0, len(d), qcut


def categorize(f, p):
    f = f.merge(p, on=["modelEntityId", "chain", "position"], how="left")
    inpk = f["do"].notna() & (f.plddt >= LINING_PLDDT)
    for col, docol in (("cat", "do"), ("cat_s02", "do_s02"), ("cat_s08", "do_s08")):
        c = np.full(len(f), 5, np.int8)
        c[(f.rsa_dimer >= SURFACE_RSA).values] = 4
        c[f.interface.values] = 3
        c[inpk.values] = 2
        c[(inpk & f[docol].fillna(False).astype(bool)).values] = 1
        f[col] = c
    b = lambda c: f[c].fillna(False).astype(bool).values
    base = np.full(len(f), 6, np.int8)
    base[(f.rsa_dimer >= SURFACE_RSA).values] = 5
    base[f.interface.values] = 4
    base[inpk.values] = 3
    q = base.copy(); q[inpk.values & b("do")] = 2; q[inpk.values & b("do_q40")] = 1; f["cat_q40"] = q
    sp = base.copy(); sp[inpk.values & b("pre")] = 2; sp[inpk.values & b("do")] = 1; f["cat_split"] = sp
    f["lining_low_plddt"] = f["do"].notna() & (f.plddt < LINING_PLDDT)
    f["sc_from_ifpocket"] = inpk & ~f["do"].fillna(False).astype(bool) & f.pk_interface.fillna(False).astype(bool)
    return f.drop(columns=["do", "do_s02", "do_s08", "pk_interface", "do_q40", "pre"])


def aggregate(f, catcol="cat"):
    """Per (accession, position): highest-priority (lowest code) category; covariates from instances with it."""
    best = f.groupby(["accession", "position"])[catcol].transform("min")
    g = f[f[catcol] == best].groupby(["accession", "position"], as_index=False).agg(
        cat=(catcol, "first"), rsa_iso=("rsa_iso", "mean"), rsa_dimer=("rsa_dimer", "mean"), plddt=("plddt", "mean"),
        resname=("resname", "first"), n_inst=("modelEntityId", "size"))
    return g


def outcomes(g):
    v = pd.read_parquet(WS / "data/processed/variants_humsavar.parquet")
    v = v[v.category == "P/LP"][["accession", "position", "wt"]]
    x = v.merge(g[["accession", "position", "resname"]], on=["accession", "position"], how="inner")
    x["match"] = x.wt == x.resname.map(AA3)
    stats = {"plp_rows_in_dimer_proteins": len(x), "plp_wt_mismatch": int((~x.match).sum())}
    ok = x[x.match].groupby(["accession", "position"]).size().rename("n_plp").reset_index()
    g = g.merge(ok, on=["accession", "position"], how="left")
    g["n_plp"] = g.n_plp.fillna(0).astype(int); g["y"] = (g.n_plp > 0).astype(int)
    am = pd.read_parquet(WS / "data/processed/alphamissense_residue_mean.parquet", columns=["accession", "position", "am_mean"])
    g = g.merge(am, on=["accession", "position"], how="left")
    return g, stats


def pdb_complex_models(models):
    s = pd.read_parquet(WS / "data/processed/sifts_chain_uniprot.parquet", columns=["pdb", "chain", "sp_primary"]).drop_duplicates()
    nch = s.groupby(["pdb", "sp_primary"]).chain.nunique()
    homo = set(nch[nch >= 2].reset_index().sp_primary)     # >=2 chains of the protein in one entry (ASU/SIFTS; approx.)
    pairs = set()
    for pdb, accs in s.groupby("pdb").sp_primary:
        u = sorted(set(accs))
        if 1 < len(u) <= 30:
            pairs.update((a, b) for i, a in enumerate(u) for b in u[i + 1:])
    out = {}
    for r in models.itertuples():
        a, b = r.uniprotAccession_1, r.uniprotAccession_2
        out[r.modelEntityId] = (a in homo) if a == b else (tuple(sorted((a, b))) in pairs)
    return out


def main():
    TAB.mkdir(parents=True, exist_ok=True); RES.mkdir(parents=True, exist_ok=True)
    prov = provenance()
    f, nchunks = load_features()
    p, npk0, npk, qcut = pocket_residues()
    f = categorize(f, p)
    models = pd.read_csv(FEAT / "models.tsv", sep="\t")
    pc = pdb_complex_models(models)
    f["pdb_complex"] = f.modelEntityId.map(pc).fillna(False).astype(bool)
    f.to_parquet(TAB / "model_residues.parquet", index=False)
    ledger = [("feature_chunks", nchunks), ("models_with_features", f.modelEntityId.nunique()),
              ("model_chain_residues", len(f)), ("pockets_annotated", npk0), ("pockets_atlas_pass_symrep", npk),
              ("pocket_lining_rows_model_chain", int((f.cat <= 2).sum() + f.lining_low_plddt.sum())),
              ("pocket_lining_rows_low_plddt_dropped", int(f.lining_low_plddt.sum())),
              ("models_with_pdb_complex", int(sum(pc.values()))), ("burrow_score_q60_cut_sensitivity", round(qcut, 5))]
    # best model per protein
    inst = f[["accession", "modelEntityId", "max_ipSAE"]].drop_duplicates().sort_values(["accession", "max_ipSAE", "modelEntityId"], ascending=[True, False, True])
    bestm = inst.drop_duplicates("accession")
    sets = {"primary": f, "sens_bestmodel": f.merge(bestm[["accession", "modelEntityId"]], on=["accession", "modelEntityId"]),
            "homo": f[f.kind == "homo"], "hetero": f[f.kind == "hetero"], "no_pdb_complex": f[~f.pdb_complex]}
    t = models.drop_duplicates("modelEntityId").max_ipSAE.quantile([1 / 3, 2 / 3]).values
    ledger += [("ipsae_tertile_cut1", round(t[0], 4)), ("ipsae_tertile_cut2", round(t[1], 4))]
    tert = np.digitize(f.max_ipSAE, t)
    for k in range(3):
        sets[f"ipsae_T{k + 1}"] = f[tert == k]
    for name, sub in sets.items():
        g = aggregate(sub)
        if name == "primary":
            for sc in ("cat_s02", "cat_s08", "cat_q40", "cat_split"):
                g = g.merge(aggregate(sub, sc)[["accession", "position", "cat"]].rename(columns={"cat": sc}), on=["accession", "position"])
            sfi = sub.groupby(["accession", "position"]).sc_from_ifpocket.max().reset_index()
            g = g.merge(sfi, on=["accession", "position"])
        g, st = outcomes(g)
        for cc in [c for c in ("cat", "cat_s02", "cat_s08", "cat_q40", "cat_split") if c in g]:
            g[cc] = pd.Categorical.from_codes(g[cc] - 1, {"cat_q40": CATS_Q40, "cat_split": CATS_SPLIT}.get(cc, CATS))
        g.to_parquet(TAB / f"protein_residues_{name}.parquet", index=False)
        ledger += [(f"{name}:proteins", g.accession.nunique()), (f"{name}:residues", len(g))]
        ledger += [(f"{name}:{k}", v) for k, v in st.items()]
        for c in CATS:
            gc = g[g.cat == c]
            ledger += [(f"{name}:{c}:residues", len(gc)), (f"{name}:{c}:plp_residues", int(gc.y.sum())),
                       (f"{name}:{c}:proteins_with_plp", gc[gc.y == 1].accession.nunique())]
        print(name, "done", flush=True)
    L = pd.DataFrame(ledger, columns=["item", "value"])
    for k, v in prov.items():
        L.loc[len(L)] = (k, v)
    L.to_csv(RES / "h5_ledger.tsv", sep="\t", index=False)
    print(L.to_string())


if __name__ == "__main__":
    main()
