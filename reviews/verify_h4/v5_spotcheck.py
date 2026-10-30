"""Verifier: PI spot-check table for the seed-0 sample (10 MBS dimer hits, 10 MBS dimer misses, 5 NAMPT) from dcc_recompute.tsv.
UNVERIFIED first attempt; failure reasons are heuristic (see v4_failure.classify)."""
import json
from pathlib import Path
import pandas as pd
from v4_failure import classify
OUT = Path(__file__).resolve().parent
R = pd.read_csv(OUT / "dcc_recompute.tsv", sep="\t")
S = R[R["sample"].isin(["T3_MBS_hit", "T3_MBS_miss", "T4_NAMPT"])].copy()
S["reason_dimer"] = [classify(r, "dimer") for _, r in S.iterrows()]
S["reason_iso"] = [classify(r, "iso") for _, r in S.iterrows()]
def best(r, lab):
    pk = json.loads(r[f"v_local_{lab}_pockets"]); b = max(pk, key=lambda p: (p[4], -p[5]))
    return f"r{b[0]}/{b[2]}res/{b[3]}/rec{b[4]:.2f}/{b[5]}A"
S["dimer_best_pocket"] = [best(r, "dimer") for _, r in S.iterrows()]
S["iso_best_pocket"] = [best(r, "iso") for _, r in S.iterrows()]
cols = ["sample", "pdb_id", "ligand_uid", "ccd", "kind", "acc1", "acc2", "chain1", "chain2", "iface_jaccard", "n_sites_pair",
        "v_local_rmsd", "dimer_dcc", "v_local_dimer_dcc_topk", "v_local_dimer_hit", "iso_dcc", "v_local_iso_dcc_topk", "v_local_iso_hit",
        "dimer_best_pocket", "iso_best_pocket", "reason_dimer", "reason_iso"]
S[cols].to_csv(OUT / "spotcheck.tsv", sep="\t", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 60)
print(S[cols].round(2).to_string(index=False))
