#!/usr/bin/env python
"""Per-pocket lining features for dimer pockets and isolated-chain pockets (Q2.2, EXPLORATORY).

Joins pocket lining residues (`residues` = 'chain:resnum') with data/interim/q22/residues.parquet
(scripts/q22_residue_features.py) on (modelEntityId, chain, resnum).
Features: n_lining, frac_hydrophobic / frac_aromatic / frac_charged (residue classes in q22_common),
buriedness = mean(1 - min(rel_sasa,1)) of lining residues (rel SASA in the dimer for dimer pockets, in the
isolated chain for isolated pockets; also the dimer-context value for isolated pockets), n_lining_interface.
Outputs: data/interim/q22/dimer_pocket_features.parquet, data/interim/q22/iso_pocket_features.parquet
Status: first attempt, not yet verified.
"""
import glob, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import INTERIM, POCKETS_ISO, HYDROPHOBIC, AROMATIC, CHARGED, ATLAS_SCORE, ATLAS_PLDDT, load_annot

res = pd.read_parquet(INTERIM / "residues.parquet")
res["hyd"] = res.aa.isin(HYDROPHOBIC); res["aro"] = res.aa.isin(AROMATIC); res["chg"] = res.aa.isin(CHARGED)
res["bur_cpx"] = 1 - res.rel_sasa_cpx.clip(0, 1); res["bur_iso"] = 1 - res.rel_sasa_iso.clip(0, 1)
res = res.set_index(["modelEntityId", "chain", "resnum"])[["hyd", "aro", "chg", "bur_cpx", "bur_iso", "is_interface", "is_surface"]]


def lining_features(df, model_col, bur_col, prefix=""):
    ex = df[["pid", model_col, "residues"]].copy()
    ex["tok"] = ex.residues.str.split()
    ex = ex.explode("tok").dropna(subset=["tok"])
    ex["chain"] = ex.tok.str[0]; ex["resnum"] = pd.to_numeric(ex.tok.str[2:], errors="coerce")
    ex = ex.dropna(subset=["resnum"]); ex["resnum"] = ex.resnum.astype(int)
    j = ex.join(res, on=[model_col, "chain", "resnum"], how="left")
    for c in ["hyd", "aro", "chg", "is_interface", "is_surface"]:
        j[c] = j[c].astype(float)  # NaN when a lining residue is not in the residue table
    g = j.groupby("pid")
    out = pd.DataFrame({"n_lining": g.size(), "n_lining_matched": g.hyd.count(),
                        "frac_hydrophobic": g.hyd.mean(), "frac_aromatic": g.aro.mean(), "frac_charged": g.chg.mean(),
                        "buriedness": g[bur_col].mean(), "buriedness_dimer_ctx": g.bur_cpx.mean(),
                        "n_lining_interface": g.is_interface.sum()})
    return out


# dimer pockets
d = load_annot()
d["pid"] = np.arange(len(d))
f = lining_features(d, "modelEntityId", "bur_cpx")
d = d.join(f, on="pid")
d.drop(columns=["pid"]).to_parquet(INTERIM / "dimer_pocket_features.parquet")
print("dimer pockets", len(d), "unmatched lining frac", 1 - d.n_lining_matched.sum() / d.n_lining.sum())

# isolated-chain pockets
fs = []
for kind in ("homo", "hetero"):
    for fn in sorted(glob.glob(str(POCKETS_ISO / kind / "chunk_[0-9][0-9][0-9][0-9].parquet"))):
        x = pd.read_parquet(fn, columns=["pocket_id", "structure_stem", "burrow_score", "volume_A3", "site_volume_A3",
                                         "residues", "plddt_mean", "n_points"])
        x["kind"] = kind; fs.append(x)
iso = pd.concat(fs, ignore_index=True)
iso["modelEntityId"] = iso.structure_stem.str.rsplit("_", n=1).str[0]
iso["chain"] = iso.structure_stem.str.rsplit("_", n=1).str[1]
iso["atlas_pass"] = (iso.burrow_score >= ATLAS_SCORE) & (iso.plddt_mean >= ATLAS_PLDDT)
iso["pid"] = np.arange(len(iso))
iso = iso.join(lining_features(iso, "modelEntityId", "bur_iso"), on="pid").drop(columns=["pid"])
iso.to_parquet(INTERIM / "iso_pocket_features.parquet")
print("iso pockets", len(iso), "models", iso.modelEntityId.nunique(), "unmatched lining frac",
      1 - iso.n_lining_matched.sum() / iso.n_lining.sum())
