"""Verifier: disentangle amino-acid vs pLDDT attenuation of the H5 primary contrast (not verified itself)."""
from pathlib import Path
import sys, numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fits import fit, ROOT, OUT
P = pd.read_parquet(f"{ROOT}/data/part2/h5/tables/protein_residues_primary.parquet"); P["cat"] = P.cat.astype(str)
D = P[P.cat != "buried"].copy()
rows = []
rows += fit(D, "H decile-binned rsa_iso & pLDDT, no aa", bins=True)
rows += fit(D[D.plddt >= 70], "I pLDDT >= 70 only")
rows += fit(D[(D.plddt >= 70) & (D.plddt < 90)], "I2 70 <= pLDDT < 90")
rows += fit(D[D.plddt >= 90], "I3 pLDDT >= 90 + aa", aa=True)
D["plddt_lt90"] = (D.plddt < 90).astype(float)
for cls, aas in [("noCys", {"CYS"}), ("noCysGlyPro", {"CYS", "GLY", "PRO"}), ("noArom", {"PHE", "TYR", "TRP", "HIS"})]:
    rows += fit(D[~D.resname.isin(aas)], f"J drop {cls}")
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/fits_extra.tsv", sep="\t", index=False)
print(R[R.term == "dimer_only_vs_interface"].to_string())
c = D.groupby("cat").plddt.describe(); print(c)
