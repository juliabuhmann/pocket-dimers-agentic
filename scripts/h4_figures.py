#!/usr/bin/env python
"""H4 figures (PNG, dpi 150) from results/part2/h4/<census>/*.tsv and data/part2/h4/<census>/benchmark_filtered.parquet.
First attempt, NOT verified. Palette: dataviz-skill reference categorical slots (validated: blue/orange/aqua, light mode).
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from h4_common import OUT_RES, OUT_FIG, OUT_DATA, DCC_HIT, provenance

COL = {"dimer": "#2a78d6", "iso": "#eb6834", "atlas": "#1baf7a", "null": "#8a8985"}
LAB = {"dimer": "Dimer model", "iso": "Isolated chains", "atlas": "AFDB monomer", "null": "Random patches (seed 0)"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DPI = 150


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)


def read(p):
    return pd.read_csv(p, sep="\t", comment="#")


def stamp(fig, census):
    pv = provenance()
    fig.text(0.01, 0.005, f"census {census} | DESIGN sha256 {pv['design_sha256'][:12]} | git {pv['git_commit']} | UNVERIFIED",
             fontsize=6.5, color=INK2)


def fig_rates(census):
    r = read(OUT_RES / census / "hit_rates.tsv")
    r = r[(r["filter"] == "filtered") & (r.lig_set == "organic")]
    kinds = [k for k in ("homo", "hetero") if k in set(r.kind)]
    fig, axes = plt.subplots(1, len(kinds), figsize=(5.2 * len(kinds), 4.2), sharey=True, squeeze=False)
    conds = [("dimer", "dimer"), ("iso", "iso"), ("dimer_null0", "null")]
    for ax, kind in zip(axes[0], kinds):
        for gi, st in enumerate(("MBS", "SBS")):
            for ci, (c, key) in enumerate(conds):
                row = r[(r.kind == kind) & (r.site_type == st) & (r.condition == c)]
                if row.empty:
                    continue
                row = row.iloc[0]
                x = gi + (ci - 1) * 0.26
                ax.bar(x, row.hit_rate, width=0.24, color=COL[key], label=LAB[key] if gi == 0 else None)
                ax.errorbar(x, row.hit_rate, yerr=[[row.hit_rate - row.ci_lo], [row.ci_hi - row.hit_rate]], color=INK, lw=1, capsize=2)
                ax.text(x, row.ci_hi + 0.01, f"{row.hit_rate:.2f}", ha="center", va="bottom", fontsize=7.5, color=INK2)
        ns = [r[(r.kind == kind) & (r.site_type == st)].n.max() for st in ("MBS", "SBS")]
        ax.set_xticks([0, 1]); ax.set_xticklabels([f"MBS\nn={ns[0]}", f"SBS\nn={ns[1]}"])
        ax.set_ylim(0, 1.12)
        ax.set_title(f"{kind}dimers", fontsize=11, color=INK, loc="left")
        style(ax)
    axes[0][0].set_ylabel(f"Hit rate (DCC ≤ {DCC_HIT:g} Å, top-(n+2))", color=INK2)
    fig.legend(*axes[0][0].get_legend_handles_labels(), frameon=False, fontsize=8, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.93))
    fig.suptitle("Recovery of PDB organic-ligand sites (interface-filtered, non-redundant; 95% cluster-bootstrap CI)", fontsize=10, color=INK)
    stamp(fig, census)
    fig.tight_layout(rect=(0, 0.02, 1, 0.88))
    fig.savefig(OUT_FIG / f"h4_hit_rates_{census}.png", dpi=DPI); plt.close(fig)


def fig_forest(census):
    t = read(OUT_RES / census / "h4_tests.tsv")
    t = t[~t.analysis.str.contains("seed [12]")].iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, 0.38 * len(t) + 1.2))
    y = np.arange(len(t))
    prim = t.role.eq("primary").values
    ax.errorbar(t.did, y, xerr=[t.did - t.did_lo, t.did_hi - t.did], fmt="none", ecolor=INK2, lw=1.2, capsize=2)
    ax.scatter(t.did, y, s=np.where(prim, 60, 30), color=np.where(prim, COL["dimer"], INK2), zorder=3)
    ax.axvline(0, color=INK2, lw=0.8, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([f"{a}  (n={n}, MBS {m})" for a, n, m in zip(t.analysis, t.n, t.n_mbs)], fontsize=7.5)
    ax.set_xlabel("DiD = dimer gain at MBS − dimer gain at SBS", color=INK2)
    style(ax); ax.yaxis.grid(False); ax.xaxis.grid(True, color=GRID)
    ax.set_title("H4 contrast, 95% cluster-bootstrap CI (primary in blue)", fontsize=10, color=INK, loc="right")
    stamp(fig, census)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT_FIG / f"h4_did_forest_{census}.png", dpi=DPI); plt.close(fig)


def fig_dcc(census):
    b = pd.read_parquet(OUT_DATA / census / "benchmark_filtered.parquet")
    b = b[b.lig_set == "organic"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
    for ax, (st, m) in zip(axes, (("MBS", True), ("SBS", False))):
        g = b[b.is_mbs == m]
        for c in ("dimer", "iso"):
            x = np.sort(np.clip(g[f"{c}_dcc"].values, 0, 30))
            ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=COL[c], lw=2, label=LAB[c])
        ax.axvline(DCC_HIT, color=INK2, lw=0.8, ls="--")
        ax.set_title(f"{st} (n={len(g)})", fontsize=10, color=INK, loc="left")
        ax.set_xlabel("Min DCC to top-(n+2) pockets (Å, clipped at 30)", color=INK2)
        style(ax)
    axes[0].set_ylabel("Cumulative fraction of sites", color=INK2)
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    stamp(fig, census)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT_FIG / f"h4_dcc_ecdf_{census}.png", dpi=DPI); plt.close(fig)


def fig_f1(census):
    r = read(OUT_RES / census / "residue_metrics.tsv")
    r = r[(r["filter"] == "filtered") & (r.lig_set == "organic") & (r.metric == "f1")]
    kinds = [k for k in ("homo", "hetero") if k in set(r.kind)]
    fig, axes = plt.subplots(1, len(kinds), figsize=(5 * len(kinds), 3.8), sharey=True, squeeze=False)
    conds = [("dimer[atlas subset]", "dimer"), ("iso[atlas subset]", "iso"), ("atlas", "atlas")]
    for ax, kind in zip(axes[0], kinds):
        for gi, st in enumerate(("MBS", "SBS")):
            for ci, (c, key) in enumerate(conds):
                row = r[(r.kind == kind) & (r.site_type == st) & (r.condition == c)]
                if row.empty:
                    continue
                row = row.iloc[0]
                x = gi + (ci - 1) * 0.26
                ax.bar(x, row["mean"], width=0.24, color=COL[key], label=LAB[key] if gi == 0 else None)
                ax.text(x, row["mean"] + 0.01, f"{row['mean']:.2f}", ha="center", va="bottom", fontsize=7.5, color=INK2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["MBS", "SBS"])
        ax.set_title(f"{kind}dimers", fontsize=11, color=INK, loc="left")
        style(ax)
    axes[0][0].set_ylabel("Mean residue F1 (best of top-(n+2))", color=INK2)
    axes[0][0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Residue overlap: dimer vs isolated chains vs AFDB monomer atlas (sites with an atlas model; secondary)",
                 fontsize=9.5, color=INK)
    stamp(fig, census)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    fig.savefig(OUT_FIG / f"h4_residue_f1_{census}.png", dpi=DPI); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    a = ap.parse_args()
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_rates, fig_forest, fig_dcc, fig_f1):
        f(a.census)
    print("figures written to", OUT_FIG)
