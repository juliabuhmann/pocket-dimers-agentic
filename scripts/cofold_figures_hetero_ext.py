#!/usr/bin/env python3
"""Figure for the EXPLORATORY heterodimer MBS extension (reads results/cofold/hetero_ext/{paired,summary}.tsv).

Status: first attempt, UNVERIFIED (not visually signed off by the user).

Rerun:  .venvs/part1/bin/python scripts/cofold_figures_hetero_ext.py -> results/cofold/hetero_ext/figures/hetero_ext.png
Left: per-case ligand RMSD, single chain (x) vs heterodimer (y), verifier metric, top-ranked seed-0 sample.
Right: partner share of the PDB site vs partner gain = log10(RMSD single / RMSD dimer).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results/cofold/hetero_ext"
SUCCESS = 2.0
COL = {"original": "#2a78d6", "extension": "#eb6834"}   # dataviz reference palette slots 1-2 (validated, light)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DPI = 150


def main():
    w = pd.read_csv(D / "paired.tsv", sep="\t").dropna(subset=["rmsd_ver_dimer", "rmsd_ver_monomer"])  # complete pairs
    s = pd.read_csv(D / "summary.tsv", sep="\t").set_index("subset")
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    a = ax[0]
    lo, hi = 0.2, max(60, w[["rmsd_ver_dimer", "rmsd_ver_monomer"]].max().max() * 1.3)
    a.plot([lo, hi], [lo, hi], color=INK2, lw=1, ls="--", zorder=1)
    a.axhline(SUCCESS, color=GRID, lw=1.5, zorder=0); a.axvline(SUCCESS, color=GRID, lw=1.5, zorder=0)
    for o, g in w.groupby("origin"):
        for tier, mk in ((1, "o"), (2, "s")):
            h = g[g.tier == tier]
            if h.empty:
                continue
            a.scatter(h.rmsd_ver_monomer, h.rmsd_ver_dimer, s=42, marker=mk, color=COL[o], edgecolor="white",
                      linewidth=1, zorder=3, label=f"{o}, tier {tier} ligand (n={len(h)})")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lo, hi); a.set_ylim(lo, hi)
    a.set_xlabel("ligand RMSD, single chain + ligand (Å)"); a.set_ylabel("ligand RMSD, heterodimer + ligand (Å)")
    al = s.loc["all"]
    a.set_title(f"Paired RMSD (n={int(al.n)}): dimer only {int(al.dimer_only)}, single only {int(al.single_only)},"
                f"\nMcNemar exact p={al.mcnemar_exact_p:.2g}; dimer lower RMSD {int(al.rmsd_dimer_lower)}/"
                f"{int(al.rmsd_dimer_lower + al.rmsd_dimer_higher)}, sign p={al.sign_test_p:.2g}", fontsize=9, color=INK)
    a.text(hi * 0.8, 0.25, "below diagonal: dimer better", ha="right", color=INK2, fontsize=8)
    a.legend(frameon=False, fontsize=7.5, loc="upper left")
    b = ax[1]
    b.axhline(0, color=INK2, lw=1, ls="--", zorder=1)
    for o, g in w.groupby("origin"):
        for tier, mk in ((1, "o"), (2, "s")):
            h = g[g.tier == tier]
            b.scatter(h.partner_share, h.partner_gain_log10, s=42, marker=mk, color=COL[o], edgecolor="white",
                      linewidth=1, zorder=3)
    b.text(0.99, 0.02, "above 0: dimer lower RMSD", transform=b.transAxes, ha="right", color=INK2, fontsize=8)
    b.set_xlabel("partner share of PDB site (minor-chain binding residues / all)")
    b.set_ylabel("partner gain, log10(RMSD single / RMSD dimer)")
    b.set_title(f"Spearman ρ={al.spearman_share_vs_gain:.2f} (p={al.spearman_p:.2g}, n={int(al.n)})", fontsize=9,
                color=INK)
    b.grid(axis="y", color=GRID, lw=0.8); b.set_axisbelow(True)
    fig.suptitle("EXPLORATORY post-hoc extension: heterodimer MBS sites, AF3 seed 0, top-ranked sample (UNVERIFIED)",
                 fontsize=9.5, color=INK2)
    fig.tight_layout()
    (D / "figures").mkdir(exist_ok=True)
    fig.savefig(D / "figures/hetero_ext.png", dpi=DPI, facecolor="white")


if __name__ == "__main__":
    main()
