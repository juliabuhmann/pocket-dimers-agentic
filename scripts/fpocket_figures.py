#!/usr/bin/env python
"""Figures for the fpocket druggability analysis (EXPLORATORY; docs/FPOCKET_DRUGGABILITY.md). PNG, dpi 150.

Inputs: data/interim/fpocket/pdbsite_sites.parquet (scripts/fpocket_compare.py), results/part1/fpocket/*.tsv,
data/interim/fpocket/afdimer_pockets.parquet (scripts/fpocket_afdimer.py).
Status: first attempt, not yet verified (figures not yet reviewed by a person).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fpocket_common import OUT, FIG, INTERIM, DRUG_THR

# dataviz skill reference palette (categorical slots 1-3, light mode) and text inks
C_MBS, C_SBS, C_3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "figure.facecolor": "white"})


def read(name):
    return pd.read_csv(OUT / name, sep="\t", comment="#")


def wecdf(v, w):
    o = np.argsort(v); v, w = np.asarray(v)[o], np.asarray(w)[o]
    return v, np.cumsum(w) / w.sum()


def fig_ecdf(d):
    panels = [("All sites (population-weighted)", d, True), ("Drug-like", d[d.stratum_lig == "drug_like"], False),
              ("Cofactor", d[d.stratum_lig == "cofactor"], False), ("Other organic", d[d.stratum_lig == "other"], False)]
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.3), sharey=True)
    for ax, (title, x, wt) in zip(axs, panels):
        for g, c, lab in ((1, C_MBS, "MBS"), (0, C_SBS, "SBS")):
            h = x[x.y_mbs == g]
            v, f = wecdf(h.fp_drug, h.w if wt else np.ones(len(h)))
            ax.step(v, f, where="post", color=c, lw=2, label=f"{lab} (n={len(h)})")
        ax.axvline(DRUG_THR, color=INK2, lw=1, ls="--")
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_xlabel("fpocket druggability score"); ax.set_xlim(0, 1)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    axs[0].set_ylabel("Cumulative fraction of sites")
    axs[0].text(DRUG_THR + 0.02, 0.97, f"{DRUG_THR} threshold", color=INK2, fontsize=7)
    fig.suptitle("fpocket druggability at PDB ligand sites, MBS vs SBS (EXPLORATORY, UNVERIFIED)", x=0.01, ha="left",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "fpocket_pdbsite_druggability_ecdf.png", dpi=150); plt.close(fig)


def fig_forest():
    cd = read("pdbsite_mbs_vs_sbs_cliffs_delta.tsv")
    oa = read("pdbsite_onaxis_vs_offaxis_mbs.tsv")
    names = {"fp_drug": "Druggability (primary)", "fp_score": "fpocket score", "fp_volume": "Volume",
             "fp_hydrophobicity": "Hydrophobicity score", "fp_polarity": "Polarity score"}
    strata = {"all_weighted": "All (weighted)", "all_unweighted": "All (unweighted)",
              "all_weighted_cover50": "All, pocket covers ≥50% of ligand", "drug_like": "Drug-like", "cofactor": "Cofactor",
              "other": "Other organic", "homomer": "Homomers", "heteromer": "Heteromers"}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    # left: druggability across strata (+ on- vs off-axis)
    x = cd[cd.metric == "fp_drug"].set_index("stratum").loc[list(strata)]
    rows = [(strata[s], r.cliffs_delta, r.cd_lo, r.cd_hi, C_MBS) for s, r in x.iterrows()]
    r = oa[oa.metric == "fp_drug"].iloc[0]
    rows.append(("Homomer MBS: on- vs off-axis", r.cliffs_delta, r.cd_lo, r.cd_hi, C_3))
    for ax, rr, title in ((axs[0], rows, "Druggability: MBS vs SBS by stratum"),
                          (axs[1], [(names[m], *cd[(cd.metric == m) & (cd.stratum == "all_weighted")][["cliffs_delta", "cd_lo", "cd_hi"]].iloc[0], C_MBS)
                                    for m in names], "All descriptors, all sites (weighted)")):
        y = np.arange(len(rr))[::-1]
        for yi, (lab, e, lo, hi, c) in zip(y, rr):
            ax.plot([lo, hi], [yi, yi], color=c, lw=2, solid_capstyle="round")
            ax.plot(e, yi, "o", color=c, ms=7, mec="white", mew=1.5)
        ax.set_yticks(y); ax.set_yticklabels([t[0] for t in rr])
        ax.axvline(0, color=INK2, lw=1)
        ax.set_xlabel("Cliff's δ, 95% cluster-bootstrap CI (MBS higher →)")
        ax.set_title(title, fontsize=9, loc="left", color=INK)
        ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "fpocket_pdbsite_forest.png", dpi=150); plt.close(fig)


def fig_afdimer():
    o = pd.read_parquet(INTERIM / "afdimer_pockets.parquet")
    groups = [("dimer_only_interface", "Dimer-only interface\n(q23a curated)", C_MBS),
              ("single_chain_matched", "Single-chain, same models\n(atlas-pass, random)", C_SBS),
              ("single_chain_strict", "Single-chain, same models\n(curated-list filters)", C_3)]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [1.3, 1]})
    rng = np.random.default_rng(0)
    for i, (g, lab, c) in enumerate(groups):
        v = o[(o.group == g) & o.fp_drug.notna()].fp_drug.values
        axs[0].scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=14, color=c, alpha=0.7, lw=0)
        axs[0].plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=2)
        axs[0].text(i, 1.03, f"n={len(v)}\n>{DRUG_THR}: {(v > DRUG_THR).mean():.0%}", ha="center", va="bottom",
                    fontsize=7.5, color=INK2)
        f = np.sort(v)
        axs[1].step(f, np.arange(1, len(f) + 1) / len(f), where="post", color=c, lw=2, label=lab.replace("\n", " "))
    axs[0].set_xticks(range(3)); axs[0].set_xticklabels([g[1] for g in groups], fontsize=8)
    axs[0].axhline(DRUG_THR, color=INK2, lw=1, ls="--"); axs[0].set_ylim(-0.03, 1.18)
    axs[0].set_ylabel("fpocket druggability\n(max-Jaccard mapped pocket)")
    axs[0].grid(axis="x", visible=False)
    axs[1].axvline(DRUG_THR, color=INK2, lw=1, ls="--"); axs[1].set_xlim(0, 1)
    axs[1].set_xlabel("fpocket druggability"); axs[1].set_ylabel("Cumulative fraction")
    axs[1].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.suptitle("fpocket druggability of BURROW pockets in AF dimer models (EXPLORATORY, UNVERIFIED); median = black bar",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "fpocket_afdimer_druggability.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    d = pd.read_parquet(INTERIM / "pdbsite_sites.parquet")
    fig_ecdf(d); fig_forest(); fig_afdimer()
    print("figures ->", FIG)
