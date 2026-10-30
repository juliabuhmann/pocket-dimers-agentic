#!/usr/bin/env python
"""Figures for the PDB-site pocket descriptors (H2 secondary + EXPLORATORY). PNG, dpi 150.

Palette: validated reference categorical palette of the dataviz skill (slot 1 blue = PDB sites, slot 2 orange = AF
dimer models / second group), fixed order, as scripts/q22_figures.py.
Inputs: results/part1/pocket_descriptors/*.tsv, data/interim/pdbsite/sites.parquet (scripts/pdbsite_compare.py).
Run with .venvs/part1/bin/python.  Status: first attempt, not yet verified.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdbsite_common import OUT, FIG, INTERIM

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "legend.frameon": False})
LAB = {"volume_A3": "volume (merged pocket)", "burrow_score": "burrow_score (merged pocket)",
       "buriedness": "buriedness (1 - rel. SASA)", "n_lining": "# lining residues", "frac_hydrophobic": "hydrophobic frac.",
       "frac_aromatic": "aromatic frac.", "frac_charged": "charged frac. (DEKR)"}


def rd(name):
    return pd.read_csv(OUT / name, sep="\t", comment="#")


def forest_vs_af():
    c = rd("pdb_vs_afdimer_pattern.tsv")
    mets = list(dict.fromkeys(c.metric))[::-1]
    series = [(None, "af", ORANGE, "AF dimer models: interface vs single-chain pocket (Q2.2)"),
              ("PDB MBS vs SBS ligand sites", "pdb", BLUE, "PDB ligand sites: MBS vs SBS (weighted)"),
              ("PDB interface vs single-chain pocket at ligand sites", "pdb", AQUA,
               "PDB ligand sites: interface vs single-chain pocket (weighted)")]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axvline(0, color=INK2, lw=0.8)
    for k, (contrast, pre, col, lab) in enumerate(series):
        x = c[c.pdb_contrast == (contrast or c.pdb_contrast.iloc[0])].set_index("metric").loc[mets]
        est = x["af_dimer_delta"] if pre == "af" else x["pdb_delta"]
        y = np.arange(len(mets)) + 0.25 - 0.25 * k
        ax.errorbar(est, y, xerr=[est - x[f"{pre}_lo"], x[f"{pre}_hi"] - est], fmt="o", ms=5, color=col, ecolor=col,
                    elinewidth=2, capsize=0, label=lab)
    ax.set_yticks(np.arange(len(mets)), [LAB.get(m, m) for m in mets])
    ax.set_xlabel("Cliff's δ (95% cluster-bootstrap CI)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.legend(loc="upper center", bbox_to_anchor=(0.35, -0.16), fontsize=8)
    ax.set_title("Predicted-dimer pattern vs real PDB sites", fontsize=10, color=INK, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "pdbsite_vs_afdimer_forest.png", dpi=150); plt.close(fig)


def by_stratum():
    d = rd("mbs_vs_sbs_cliffs_delta.tsv")
    d = d[d.seed == 0]
    mets = ["burrow_score", "merged_volume_A3", "buriedness", "frac_charged", "frac_hydrophobic", "frac_aromatic"]
    strata = [("all_weighted", "all (weighted)"), ("drug_like", "drug-like"), ("cofactor", "cofactor"),
              ("other", "other organic"), ("homomer", "homomer"), ("heteromer", "heteromer")]
    fig, axes = plt.subplots(1, len(mets), figsize=(12, 2.9), sharey=True)
    for ax, m in zip(axes, mets):
        x = d[d.metric == m].set_index("stratum")
        y = np.arange(len(strata))[::-1]
        for yy, (s, _) in zip(y, strata):
            if s in x.index:
                r = x.loc[s]
                col = BLUE if s == "all_weighted" else INK2
                ax.plot([r.cd_lo, r.cd_hi], [yy, yy], color=col, lw=2, solid_capstyle="round")
                ax.plot(r.cliffs_delta, yy, "o", color=col, ms=5)
        ax.axvline(0, color=INK2, lw=0.8); ax.grid(axis="x", color=GRID, lw=0.6)
        t = "burrow_score (H2 secondary)" if m == "burrow_score" else m.replace("merged_", "").replace("_A3", "")
        ax.set_title(t, fontsize=9, color=INK)
        ax.set_yticks(y, [l for _, l in strata])
    axes[0].set_xlabel("Cliff's δ MBS vs SBS")
    fig.tight_layout(); fig.savefig(FIG / "pdbsite_mbs_vs_sbs_by_stratum.png", dpi=150); plt.close(fig)


def ecdf_burrow():
    s = pd.read_parquet(INTERIM / "sites.parquet")
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 2.9), sharey=True)
    for ax, st, lab in zip(axes, ["drug_like", "cofactor", "other"], ["drug-like", "cofactor", "other organic"]):
        for y, col, nm in ((1, BLUE, "MBS"), (0, ORANGE, "SBS")):
            v = np.sort(s[(s.stratum_lig == st) & (s.y_mbs == y)].burrow_score.dropna().values)
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=col, lw=2, label=f"{nm} (n={len(v)})")
        ax.set_title(lab, fontsize=9, color=INK); ax.set_xlabel("burrow_score of the pocket at the ligand")
        ax.grid(color=GRID, lw=0.6); ax.legend(fontsize=8, loc="lower right")
    axes[0].set_ylabel("cumulative fraction")
    fig.tight_layout(); fig.savefig(FIG / "pdbsite_burrow_score_ecdf.png", dpi=150); plt.close(fig)


def onaxis():
    d = rd("onaxis_vs_offaxis_mbs.tsv")
    d = d[(d.stratum == "all") & d.metric.isin(["burrow_score", "merged_volume_A3", "buriedness", "frac_charged",
                                                "frac_hydrophobic", "frac_aromatic", "n_lining"])].iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    y = np.arange(len(d))
    ax.axvline(0, color=INK2, lw=0.8)
    for yy, (_, r) in zip(y, d.iterrows()):
        ax.plot([r.cd_lo, r.cd_hi], [yy, yy], color=BLUE, lw=2, solid_capstyle="round")
        ax.plot(r.cliffs_delta, yy, "o", color=BLUE, ms=5)
    ax.set_yticks(y, [m.replace("merged_", "") for m in d.metric]); ax.grid(axis="x", color=GRID, lw=0.6)
    n = d.iloc[0]
    ax.set_xlabel(f"Cliff's δ on-axis (n={n.n_on}) vs off-axis (n={n.n_off}) homomer MBS")
    fig.tight_layout(); fig.savefig(FIG / "pdbsite_onaxis_vs_offaxis.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    forest_vs_af(); by_stratum(); ecdf_burrow(); onaxis()
    print("figures ->", FIG)
