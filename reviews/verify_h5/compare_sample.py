"""Verifier: aggregate rebuilt sample to protein x position and compare with producer table (confusion matrix)."""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = str(Path(__file__).resolve().parents[2]); OUT = f"{ROOT}/reviews/verify_h5/out"
NAMES = {1: "dimer_only", 2: "sc_pocket", 3: "interface", 4: "surface", 5: "buried"}
F = pd.read_parquet(f"{OUT}/sample_model_residues.parquet")
samp = set(pd.read_csv(f"{OUT}/sample_accessions.tsv", sep="\t").accession)
prod = pd.read_parquet(f"{ROOT}/data/part2/h5/tables/protein_residues_primary.parquet")
prod = prod[prod.accession.isin(samp)].copy(); prod["cat"] = prod.cat.astype(str)
res = {}
for v in ["nocollapse", "collapse"]:
    c = f"cat_{v}"
    agg = F.groupby(["accession", "position"])[c].min().rename("code").reset_index()
    sel = F.merge(agg, on=["accession", "position"]); sel = sel[sel[c] == sel.code]
    cov = sel.groupby(["accession", "position"]).agg(rsa_iso=("rsa_iso", "mean"), plddt=("plddt", "mean"), resname=("resname", "first")).reset_index()
    agg = agg.merge(cov, on=["accession", "position"]); agg["cat"] = agg.code.map(NAMES)
    m = agg.merge(prod[["accession", "position", "cat", "rsa_iso", "plddt"]], on=["accession", "position"], how="outer", suffixes=("_ver", "_prod"), indicator=True)
    print(v, m._merge.value_counts().to_dict())
    m = m[m._merge == "both"]
    cm = pd.crosstab(m.cat_ver, m.cat_prod, margins=True)
    print(cm.to_string()); acc = (m.cat_ver == m.cat_prod).mean(); print("agreement", round(acc, 5))
    print("rsa_iso |diff|>0.01:", (abs(m.rsa_iso_ver - m.rsa_iso_prod) > 0.01).mean(), " plddt |diff|>0.5:", (abs(m.plddt_ver - m.plddt_prod) > 0.5).mean())
    cm.to_csv(f"{OUT}/confusion_{v}.tsv", sep="\t")
    agg.to_parquet(f"{OUT}/sample_protein_residues_{v}.parquet")
    res[v] = acc
