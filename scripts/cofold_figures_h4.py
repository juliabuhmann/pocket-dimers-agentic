#!/usr/bin/env python3
"""Figures for deep dive #2 co-folding (results/cofold/h4/metrics*.tsv). Status: first attempt, UNVERIFIED.
-> results/cofold/h4/figures/{lig_rmsd_by_case.png, rates_by_stratum.png}  (PNG, dpi 150)
Colors: dataviz-skill reference palette, categorical slots 1-2 (blue = dimer run A, orange = single chain run B),
same assignment as deep dive #1 (scripts/cofold_figures_onaxis.py)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/cofold/h4"
FIG = R / "figures"
SUCCESS_RMSD = 2.0
DCC_HIT = 4.0
COL = {"dimer": "#2a78d6", "monomer": "#eb6834"}
CLAB = {"dimer": "A: dimer + ligand", "monomer": "B: single chain + ligand"}
INK2, GRID = "#52514e", "#e4e3df"
ORDER = ["homo_miss", "homo_hit", "hetero"]
LABEL = {"homo_miss": "Homodimer, apo-dimer MISS", "homo_hit": "Homodimer, apo-dimer HIT", "hetero": "Heterodimer"}


def style(ax, grid="x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(R / "metrics.tsv", sep="\t")
    s = pd.read_csv(R / "metrics_samples.tsv", sep="\t")
    s = s[s.status == "ok"]
    m = m[m.lig_rmsd.notna()].copy()
    m["o"] = m.stratum.map({k: i for i, k in enumerate(ORDER)})
    cases = m.drop_duplicates("case_id").sort_values(["o", "post_cutoff", "case_id"]).set_index("case_id")
    y = {c: i for i, c in enumerate(cases.index[::-1])}

    # 1) per-case ligand RMSD, top-ranked sample + range of 5 samples
    fig, ax = plt.subplots(figsize=(7.5, 0.3 * len(y) + 1.6))
    for cond, dy in (("dimer", 0.15), ("monomer", -0.15)):
        for cid, g in s[s.cond == cond].groupby("case_id"):
            if cid in y:
                ax.plot([g.lig_rmsd.min(), g.lig_rmsd.max()], [y[cid] + dy] * 2, color=COL[cond], lw=1.2, alpha=0.5)
        mm = m[m.cond == cond]
        ax.scatter(mm.lig_rmsd, [y[c] + dy for c in mm.case_id], s=36, color=COL[cond], edgecolor="white", lw=1,
                   zorder=3, label=f"{CLAB[cond]} (dot = top-ranked, line = 5 samples)")
    ax.axvline(SUCCESS_RMSD, color=INK2, lw=1, ls="--")
    ax.text(SUCCESS_RMSD, 0.0, " 2 Å", transform=ax.get_xaxis_transform(), color=INK2, fontsize=7.5, va="bottom")
    ax.set_xscale("log")
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels([f"{c[len(cases.stratum[c]) + 1:]}{' *' if cases.post_cutoff[c] else ''}" for c in y], fontsize=7.5)
    prev = None
    for c in cases.index:
        st = cases.stratum[c]
        if st != prev:
            ax.axhline(y[c] + 0.5, color=GRID, lw=1)
            ax.text(0.995, y[c] + 0.42, LABEL[st], transform=ax.get_yaxis_transform(), ha="right", va="top",
                    fontsize=8, color=INK2, fontweight="bold")
            prev = st
    ax.set_xlabel("Ligand heavy-atom RMSD to PDB after protein superposition (Å, log)", fontsize=8.5)
    ax.set_title("AF3 co-folding at H4 multi-chain sites (seed 0; * = PDB entry released after 2021-09-30)",
                 fontsize=9, loc="left")
    style(ax)
    ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=1)
    fig.tight_layout()
    fig.savefig(FIG / "lig_rmsd_by_case.png", dpi=150)
    plt.close(fig)

    # 2) success (RMSD <= 2 A) and pocket-level DCC <= 4 A by stratum, A vs B
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
    for ax, (col, ttl) in zip(axs, (("success", f"Ligand RMSD ≤ {SUCCESS_RMSD:g} Å"),
                                    ("dcc_hit", f"Ligand centroid DCC ≤ {DCC_HIT:g} Å"))):
        x = np.arange(len(ORDER))
        for k, (cond, dx) in enumerate((("dimer", -0.18), ("monomer", 0.18))):
            fr, ns, ks = [], [], []
            for st in ORDER:
                g = m[(m.stratum == st) & (m.cond == cond)]
                ks.append(int(g[col].astype(bool).sum())); ns.append(len(g))
                fr.append(ks[-1] / len(g) if len(g) else 0)
            ax.bar(x + dx, fr, width=0.34, color=COL[cond], label=CLAB[cond], edgecolor="white", lw=2)
            for xi, f, kk, nn in zip(x + dx, fr, ks, ns):
                ax.text(xi, f + 0.02, f"{kk}/{nn}", ha="center", fontsize=7.5, color=INK2)
        ax.set_xticks(x)
        ax.set_xticklabels([LABEL[s].replace(", ", ",\n") for s in ORDER], fontsize=7.5)
        ax.set_title(ttl + " (top-ranked sample)", fontsize=9, loc="left")
        ax.set_ylim(0, 1.1)
        style(ax, grid="y")
    axs[0].set_ylabel("Fraction of cases", fontsize=8.5)
    axs[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "rates_by_stratum.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
