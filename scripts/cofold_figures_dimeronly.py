#!/usr/bin/env python3
"""Figures for deep dive #3 co-folding (dimer-only pockets; results/cofold/dimeronly/metrics*.tsv).
Status: first attempt, UNVERIFIED. Adapted from scripts/cofold_figures_onaxis.py.
-> results/cofold/dimeronly/figures/{lig_rmsd_by_case.png, pocket_placement.png}  (PNG, dpi 150)
Colors: dataviz-skill reference palette, categorical slots 1-2 (blue = dimer, orange = monomer)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/cofold/dimeronly"
FIG = R / "figures"
SUCCESS_RMSD = 2.0   # = scripts/cofold_metrics_dimeronly.py
COL = {"dimer": "#2a78d6", "monomer": "#eb6834"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
LOG_TICKS = [0.25, 0.5, 1, 2, 5, 10, 20, 40]


def style(ax, grid="x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def logx(ax):
    ax.set_xscale("log")
    ax.set_xticks(LOG_TICKS)
    ax.set_xticklabels([f"{t:g}" for t in LOG_TICKS])
    ax.minorticks_off()


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(R / "metrics.tsv", sep="\t")
    s = pd.read_csv(R / "metrics_samples.tsv", sep="\t")
    s = s[s.status == "ok"]
    m = m[m.lig_rmsd_sym.notna()].copy()
    cases = pd.read_csv(R / "cases.tsv", sep="\t").set_index("case_id")
    order = [c for c in cases.index if c in set(m.case_id)]
    y = {c: i for i, c in enumerate(order[::-1])}
    lab = {c: f"{c}{' *' if cases.post_cutoff[c] else ''}" for c in order}

    # 1) ligand RMSD per case, global vs pocket superposition
    fig, axes = plt.subplots(1, 2, figsize=(9, 0.42 * len(y) + 1.6), sharey=True)
    for ax, col, title in ((axes[0], "lig_rmsd_sym", "Global CA superposition"),
                           (axes[1], "lig_rmsd_sym_pocket", "Pocket CA superposition")):
        for cond, dy in (("dimer", 0.15), ("monomer", -0.15)):
            for cid, g in s[s.cond == cond].groupby("case_id"):
                if cid in y:
                    ax.plot([g[col].min(), g[col].max()], [y[cid] + dy] * 2, color=COL[cond], lw=1.2, alpha=0.5)
            mm = m[m.cond == cond]
            ax.scatter(mm[col], [y[c] + dy for c in mm.case_id], s=36, color=COL[cond], edgecolor="white", lw=1,
                       zorder=3, label=f"{cond} (top-ranked; line = range of 5 samples)")
        ax.axvline(SUCCESS_RMSD, color=INK2, lw=1, ls="--")
        logx(ax)
        style(ax)
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_xlabel("Symmetry-corrected ligand RMSD (Å, log)", fontsize=8, color=INK2)
    axes[0].set_yticks(list(y.values()))
    axes[0].set_yticklabels([lab[c] for c in y], fontsize=7.5)
    axes[0].legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(1.0, -0.13 - 0.6 / len(y)), ncol=2)
    fig.suptitle("AF3 seed 0, dimer-only pockets: ligand pose vs PDB (dashed = 2 Å; * = released after 2021-09-30)",
                 fontsize=9, color=INK, x=0.01, ha="left")
    fig.subplots_adjust(left=0.22, right=0.98, wspace=0.08, top=0.88)
    fig.savefig(FIG / "lig_rmsd_by_case.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 2) where does the ligand go: distance to curated pocket centre; contacts on chain A/B (top-ranked)
    fig, axes = plt.subplots(1, 2, figsize=(9, 0.42 * len(y) + 1.6), sharey=True)
    ax = axes[0]
    for cond, dy in (("dimer", 0.15), ("monomer", -0.15)):
        mm = m[m.cond == cond]
        ax.scatter(mm.dcc_curated, [y[c] + dy for c in mm.case_id], s=36, color=COL[cond], edgecolor="white", lw=1,
                   zorder=3, label=cond)
    style(ax)
    ax.set_xlabel("Ligand centroid to curated-pocket centre (Å)", fontsize=8, color=INK2)
    ax.set_title("Placement in the curated pocket", fontsize=9, color=INK, loc="left")
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels([lab[c] for c in y], fontsize=7.5)
    ax.legend(fontsize=7, frameon=False)
    ax = axes[1]
    d = m[m.cond == "dimer"].set_index("case_id")
    mo = m[m.cond == "monomer"].set_index("case_id")
    yy = np.array([y[c] for c in d.index])
    ax.barh(yy + 0.15, d.n_res_chainA, height=0.28, color=COL["dimer"], label="dimer: chain A")
    ax.barh(yy + 0.15, d.n_res_chainB, left=d.n_res_chainA + 0.3, height=0.28, color="#86b6ef", label="dimer: chain B")
    ax.barh(np.array([y[c] for c in mo.index]) - 0.15, mo.n_res_chainA, height=0.28, color=COL["monomer"],
            label="monomer: single chain")
    style(ax)
    ax.set_xlabel("Residues within 4 Å of the matched ligand copy", fontsize=8, color=INK2)
    ax.set_title("Contacts per chain (top-ranked)", fontsize=9, color=INK, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "pocket_placement.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", FIG)


if __name__ == "__main__":
    main()
