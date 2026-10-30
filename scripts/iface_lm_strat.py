#!/usr/bin/env python
"""EXPLORATORY: IM vs IU lining-chemistry difference after stratifying on burrow_score (IM quintile edges, IM-count weights).
Statistic: weighted mean over strata of (mean_IM - mean_IU); cluster bootstrap (paralog unit, seed 0). Not verified."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from q22_common import INTERIM, WS, cluster_boot, write_tsv

G = pd.read_parquet(INTERIM.parent / "iface_ligand_matched/groups.parquet")
D = G[G.grp.isin(["IM", "IU"])].reset_index(drop=True)
edges = np.quantile(D.loc[D.grp == "IM", "burrow_score"], [.2, .4, .6, .8])
D["stratum"] = np.digitize(D.burrow_score, edges)


def stat_fn(met):
    def f(x):
        num = den = 0.0
        for _, g in x.groupby("stratum"):
            a, b = g.loc[g.grp == "IM", met], g.loc[g.grp == "IU", met]
            if len(a) and len(b):
                num += len(a) * (a.mean() - b.mean()); den += len(a)
        return num / den
    return f


rows = []
for met in ["frac_charged", "frac_hydrophobic", "frac_aromatic", "volume_A3", "buriedness"]:
    e, lo, hi = cluster_boot(D, "unit", stat_fn(met), 1000, 0)
    rows.append(dict(metric=met, stratified_mean_diff_IM_minus_IU=e, lo=lo, hi=hi, n_IM=int((D.grp == "IM").sum()), n_IU=int((D.grp == "IU").sum()),
                     stratum_edges_burrow_score=",".join(f"{v:.3f}" for v in edges)))
r = pd.DataFrame(rows); write_tsv(r, WS / "results/part2/iface_ligand_matched/burrow_stratified_diff.tsv", "IM-IU mean diff within burrow_score strata"); print(r.round(3).to_string())
