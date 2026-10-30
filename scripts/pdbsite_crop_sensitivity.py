#!/usr/bin/env python
"""Cropped (30 Å, sensitivity) vs uncropped (whole assembly 1, primary) PDB-site pocket descriptors.

Inputs: data/interim/pdbsite/sites.parquet (uncropped) and sites_crop30.parquet (both from scripts/pdbsite_compare.py),
        results/part1/pocket_descriptors/{,crop30/}mbs_vs_sbs_cliffs_delta.tsv.
Outputs: results/part1/pocket_descriptors/crop_vs_uncropped_sites.tsv (per-site agreement, by whether the crop removed
         chains and by MBS/SBS), crop_vs_uncropped_deltas.tsv (pooled weighted MBS-vs-SBS deltas side by side).
Run with .venvs/part1/bin/python.  Status: first attempt, not yet verified.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdbsite_common import INTERIM, OUT, write_tsv

ABS_TOL = 1e-6           # equality tolerance for "unchanged"
SCORE_CHANGE = 0.05      # our choice: |Δ burrow_score| above which a site counts as changed


def jac(a, b):
    a, b = set(str(a).split()), set(str(b).split())
    return len(a & b) / len(a | b) if a | b else np.nan


u = pd.read_parquet(INTERIM / "sites.parquet")
c = pd.read_parquet(INTERIM / "sites_crop30.parquet")
cols = ["site_id", "burrow_score", "volume_A3", "merged_volume_A3", "merged_burrow_score", "d_match", "residues",
        "merged_frac_charged", "merged_frac_hydrophobic", "n_chains_kept", "n_chains_asm"]
m = u[cols + ["y_mbs"]].merge(c[cols], on="site_id", suffixes=("_u", "_c"))
m["crop_removed_chains"] = m.n_chains_kept_c < m.n_chains_asm_c
rows = []
for (rem, mbs), g in [((None, None), m)] + list(m.groupby(["crop_removed_chains", "y_mbs"])):
    dbs = (g.burrow_score_u - g.burrow_score_c).abs()
    rows.append(dict(crop_removed_chains="all" if rem is None else bool(rem), y_mbs="all" if mbs is None else int(mbs),
                     n=len(g), spearman_burrow_score=spearmanr(g.burrow_score_u, g.burrow_score_c)[0],
                     frac_burrow_score_changed_gt005=(dbs > SCORE_CHANGE).mean(), median_abs_d_burrow_score=dbs.median(),
                     frac_volume_unchanged=((g.volume_A3_u - g.volume_A3_c).abs() < ABS_TOL).mean(),
                     frac_dmatch_lost=(g.d_match_c.astype(bool) & ~g.d_match_u.astype(bool)).mean(),
                     frac_dmatch_gained=(~g.d_match_c.astype(bool) & g.d_match_u.astype(bool)).mean(),
                     median_jaccard_lining=np.nanmedian([jac(a, b) for a, b in zip(g.residues_u, g.residues_c)]),
                     spearman_merged_frac_charged=spearmanr(g.merged_frac_charged_u, g.merged_frac_charged_c, nan_policy="omit")[0]))
t = pd.DataFrame(rows)
write_tsv(t, OUT / "crop_vs_uncropped_sites.tsv", note="u=uncropped (primary), c=crop30 (sensitivity)")
print(t.to_string())

du = pd.read_csv(OUT / "mbs_vs_sbs_cliffs_delta.tsv", sep="\t", comment="#")
dc = pd.read_csv(OUT / "crop30/mbs_vs_sbs_cliffs_delta.tsv", sep="\t", comment="#")
k = ["metric", "stratum", "seed"]
x = du[k + ["n_g1", "n_g0", "cliffs_delta", "cd_lo", "cd_hi"]].merge(
    dc[k + ["n_g1", "n_g0", "cliffs_delta", "cd_lo", "cd_hi"]], on=k, suffixes=("_uncropped", "_crop30"))
x = x[x.seed == 0].drop(columns="seed")
x["diff"] = x.cliffs_delta_uncropped - x.cliffs_delta_crop30
write_tsv(x, OUT / "crop_vs_uncropped_deltas.tsv", note="MBS vs SBS Cliff's delta, uncropped (primary) vs crop30")
print(x[x.stratum == "all_weighted"].to_string())
