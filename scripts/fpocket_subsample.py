#!/usr/bin/env python
"""Stratified subsample of the PDB-site sample for fpocket (EXPLORATORY; docs/FPOCKET_DRUGGABILITY.md).

Why: the 50-site pilot (logs/fpocket/pilot50.log) took 13.7 CPU-s/site on average (whole assemblies, heavy tail),
so all 11,983 sites (~4.6 h on 10 workers) did not fit the ~1.5 h compute budget of this task.
Rule: N_PER_STRATUM sites per (ligand class x MBS/SBS) stratum, drawn at random with SEED from
data/interim/pdbsite/sites_crop30.parquet. Pooled contrasts are re-weighted by units_available / n_subsampled.
Status: first attempt, not yet verified.
"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fpocket_common import WS, OUT

SEED = 0               # DESIGN.md §5: seed 0 primary
N_PER_STRATUM = 500    # sized from the pilot: 6 x 500 x 13.7 CPU-s / 9 workers ~ 76 min
s = pd.read_parquet(WS / "data/interim/pdbsite/sites_crop30.parquet", columns=["site_id", "stratum_lig", "y_mbs", "cluster30"])
sub = pd.concat([g.sample(n=min(N_PER_STRATUM, len(g)), random_state=SEED)
                 for _, g in s.sort_values("site_id").groupby(["stratum_lig", "y_mbs"])])
with open(OUT / "pdbsite_subsample_seed0.tsv", "w") as fh:
    fh.write(f"# status=UNVERIFIED EXPLORATORY fpocket subsample: {N_PER_STRATUM}/stratum seed={SEED} from sites_crop30.parquet\n")
    sub.to_csv(fh, sep="\t", index=False)
print(sub.groupby(["stratum_lig", "y_mbs"]).agg(n=("site_id", "size"), clusters=("cluster30", "nunique")))
