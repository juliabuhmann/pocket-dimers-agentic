"""Random spot-check of drug-like vs non-drug-like organic CCDs among the verifier's units (seed 0, 10 each)."""
from pathlib import Path
import pandas as pd
OUT = Path(__file__).resolve().parent / "out"
SEED, N = 0, 10  # DESIGN §5 primary seed; task asks for 10 per label
u = pd.read_parquet(OUT / "units_verifier.parquet")
dl = pd.read_parquet(Path(__file__).resolve().parents[2] / "data/processed/ccd_druglike.parquet",
                     columns=["ccd", "max_phase_sqlite", "max_pchembl", "druglike_reason"])
org = u[~u.fine.isin(["metal_ion", "metal_cluster", "inorganic"])].drop_duplicates("ccd_key")
a = org[org.fine == "drug_like"].sample(N, random_state=SEED)
b = org[org.fine.isin(["other_organic", "amino_acid"])].sample(N, random_state=SEED)  # non-drug-like organics that are not cofactor/heme/etc.
c = org[org.fine == "cofactor"].sample(N, random_state=SEED)
r = pd.concat([a.assign(set="drug_like"), b.assign(set="non_druglike_other"), c.assign(set="cofactor")])
r = r.merge(dl, left_on="ccd_key", right_on="ccd", how="left")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 55)
cols = ["set", "ccd_key", "name", "fine", "max_phase_sqlite", "max_pchembl", "druglike_reason"]
print(r[cols].to_string()); r[cols].to_csv(OUT / "label_spot_check.tsv", sep="\t", index=False)
