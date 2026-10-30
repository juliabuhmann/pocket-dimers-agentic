#!/usr/bin/env python
"""Q2.4 figures (PNG, dpi 150) from results/part2/q24/*.tsv. Status: first attempt, not yet verified.
Palette: dataviz skill reference categorical slots 1-4 (validated light mode; contrast WARN -> legend + direct labels)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WS = Path(__file__).resolve().parents[1]
OUT = WS / "results/part2/q24"
FIG = WS / "results/part2/figures"
DPI = 150
COL = {"interface": "#2a78d6", "single": "#eb6834", "dimer_only": "#1baf7a", "preexisting": "#eda100"}
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6})


def rd(name):
    return pd.read_csv(OUT / name, sep="\t", comment="#")


def species_order(T):
    s = T[T.divergence_mya.notna()].drop_duplicates("scope").sort_values(["divergence_mya", "scope"])
    return list(s.scope), list(s.divergence_mya)


def panel(ax, T, classes, order, div, null=True):
    x = np.arange(len(order))
    for k, cl in enumerate(classes):
        d = T[T.contrast == cl].set_index("scope").reindex(order)
        off = (k - (len(classes) - 1) / 2) * 0.18
        ax.errorbar(x + off, d.estimate, yerr=[d.estimate - d.lo, d.hi - d.estimate], fmt="o", ms=5, lw=1.5, color=COL[cl],
                    mec="white", mew=1, label=f"{cl.replace('_', '-')} (observed)")
        if null:
            ax.plot(x + off, d.null, "_", ms=9, mew=2, color=COL[cl], alpha=0.6, label=f"{cl.replace('_', '-')} (non-ortholog null)")
    ax.set_xticks(x); ax.set_xticklabels([f"{s}\n{d:.0f}" for s, d in zip(order, div)], fontsize=7.5)
    ax.set_ylim(0, 1); ax.set_ylabel("fraction of human pockets conserved\n(overlap ≥0.5 with an ortholog-dimer pocket)")
    ax.set_xlabel("species (divergence from human, Mya; bacteria assumed)")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    T = rd("q24_conservation.tsv")
    T = T[T.target.str.startswith("primary")]
    order, div = species_order(T)
    fig, axs = plt.subplots(2, 1, figsize=(9, 7.5), constrained_layout=True)
    panel(axs[0], T, ["interface", "single"], order, div)
    axs[0].set_title("Interface vs single-chain pockets", loc="left", fontsize=10, color=INK)
    panel(axs[1], T, ["dimer_only", "preexisting"], order, div)
    axs[1].set_title("Interface pockets: dimer-only vs pre-existing in the isolated chain", loc="left", fontsize=10, color=INK)
    for ax in axs:
        ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right")
    fig.suptitle("Q2.4 (exploratory, unverified): human dimer pocket conservation in ortholog dimers; 95% protein-cluster bootstrap CI",
                 fontsize=9, color=INK2, x=0.01, ha="left")
    fig.savefig(FIG / "q24_conservation_by_species.png", dpi=DPI); plt.close(fig)

    S = rd("q24_conservation_strata.tsv")
    for strat, xl in (("score_q", "burrow_score quintile (human pockets)"), ("id_bin", "ortholog pair identity, % (min of chains)")):
        scopes = list(dict.fromkeys(S.scope))
        fig, axs = plt.subplots(1, len(scopes), figsize=(2.3 * len(scopes), 3.2), sharey=True, constrained_layout=True)
        for ax, sc in zip(np.atleast_1d(axs), scopes):
            d0 = S[(S.scope == sc) & (S.stratum == strat)]
            levels = list(dict.fromkeys(d0.level))
            x = np.arange(len(levels))
            for k, cl in enumerate(["interface", "single"]):
                d = d0[d0.cls == cl].set_index("level").reindex(levels)
                ax.errorbar(x + (k - 0.5) * 0.2, d.estimate, yerr=[d.estimate - d.lo, d.hi - d.estimate], fmt="o-", ms=4, lw=1.2,
                            color=COL[cl], mec="white", label=cl)
                ax.plot(x + (k - 0.5) * 0.2, d.null, "_", ms=8, mew=2, color=COL[cl], alpha=0.6)
            ax.set_xticks(x); ax.set_xticklabels(levels, fontsize=7, rotation=45 if strat == "id_bin" else 0)
            ax.set_title(sc, fontsize=9, color=INK); ax.set_ylim(0, 1)
        np.atleast_1d(axs)[0].set_ylabel("fraction conserved (dash = null)")
        np.atleast_1d(axs)[0].legend(frameon=False, fontsize=7)
        fig.supxlabel(xl, fontsize=8.5, color=INK)
        fig.savefig(FIG / f"q24_conservation_by_{strat}.png", dpi=DPI); plt.close(fig)

    L = rd("q24_lining_identity.tsv")
    L = L[(L.metric == "lining_minus_chain") & L.scope.isin(order)]
    fig, ax = plt.subplots(figsize=(9, 3.6), constrained_layout=True)
    x = np.arange(len(order))
    for k, cl in enumerate(["interface", "single"]):
        d = L[L.cls == cl].set_index("scope").reindex(order)
        ax.errorbar(x + (k - 0.5) * 0.2, d["mean"], yerr=[d["mean"] - d.lo, d.hi - d["mean"]], fmt="o", ms=5, color=COL[cl], mec="white", label=cl)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"{s}\n{d:.0f}" for s, d in zip(order, div)], fontsize=7.5)
    ax.set_ylabel("lining-residue identity −\nchain alignment identity")
    ax.set_title("Lining-residue sequence conservation relative to the whole chain (exploratory, unverified)", loc="left", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "q24_lining_identity.png", dpi=DPI); plt.close(fig)


if __name__ == "__main__":
    main()
