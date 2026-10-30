#!/usr/bin/env python
"""H5 figures (PNG, dpi 150): P/LP rate per category with protein-bootstrap CIs; forest plot of GLMM ORs.
Reads results/part2/h5/*.tsv. First attempt, not independently verified.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.ticker

WS = Path(__file__).resolve().parents[1]
RES = WS / "results/part2/h5"; FIG = WS / "results/part2/figures"
DPI = 150
# dataviz skill reference palette (light): series-1 blue, text, muted
BLUE, INK, MUTED, GRID, SURF = "#2a78d6", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
LABELS = {"dimer_only": "Dimer-only pocket", "sc_pocket": "Single-chain pocket", "interface": "Interface, no pocket (ref)",
          "surface": "Other surface"}


def read(name):
    return pd.read_csv(RES / name, sep="\t", comment="#")


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)


def rates_fig(prov):
    r = read("h5_category_rates.tsv").iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 3.4), facecolor=SURF); style(ax)
    y = np.arange(len(r))
    ax.hlines(y, r.ci_lo * 100, r.ci_hi * 100, color=BLUE, lw=2)
    ax.plot(r.plp_rate * 100, y, "o", ms=8, color=BLUE, mec=SURF, mew=2)
    for yi, row in zip(y, r.itertuples()):
        ax.text(row.ci_hi * 100, yi, f"  {row.plp_rate*100:.2f}%  ({row.n_plp:,}/{row.n_residues:,})", color=MUTED, fontsize=8, va="center")
    ax.set_yticks(y, [LABELS[c] for c in r.category])
    ax.set_xlabel("Residues with ≥1 P/LP missense variant (humsavar), %", color=INK)
    ax.set_xlim(0, r.ci_hi.max() * 100 * 1.6)
    ax.set_title("P/LP rate per residue category\n95% CI: protein bootstrap (1,000 resamples, seed 0)", color=INK, fontsize=10, loc="left")
    fig.text(0.01, 0.01, f"Unadjusted; buried residues excluded. {prov}", color=MUTED, fontsize=6)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIG / "h5_plp_rate_by_category.png", dpi=DPI); plt.close(fig)


def forest_fig(prov):
    m = read("h5_model_results.tsv")
    m = m[m.outcome == "plp"]
    a = m[m.term == "catdimer_only"].copy()
    a["row"] = a.analysis + " [" + a.label + "]"
    b = m[(m.analysis == "primary_plp") & m.term.str.startswith("cat")].copy()
    b["row"] = b.term.str[3:].map(LABELS) + " vs interface"
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.32 * len(a) + 1.6), facecolor=SURF,
                             gridspec_kw={"width_ratios": [1.4, 1]})
    for ax, d, title in ((axes[0], a, "Dimer-only pocket vs interface (all analyses)"),
                         (axes[1], b, "Primary model: all contrasts")):
        style(ax); d = d.iloc[::-1]; y = np.arange(len(d))
        col = [MUTED if (l == "NEG_CONTROL" or g is False or str(g) == "False") else BLUE for l, g in zip(d.label, d.power_gate_pass)]
        ax.hlines(y, d.OR_lo, d.OR_hi, color=col, lw=2)
        ax.scatter(d.OR, y, s=50, c=col, edgecolors=SURF, linewidths=1.5, zorder=3)
        ax.axvline(1, color=MUTED, lw=1, ls="--")
        ax.set_xscale("log"); ax.set_yticks(y, d.row, fontsize=8)
        ticks = [0.5, 0.7, 0.8, 1, 1.2, 1.5]
        lo_, hi_ = ax.get_xlim()
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([t for t in ticks if lo_ <= t <= hi_]))
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlabel("Odds ratio, GLMM (95% Wald CI, log)", color=INK)
        ax.set_title(title, color=INK, fontsize=10, loc="left")
    fig.text(0.01, 0.01, f"Grey = negative control or power gate failed (<30 proteins with P/LP per category). {prov}", color=MUTED, fontsize=6)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIG / "h5_forest.png", dpi=DPI); plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    prov = open(RES / "h5_category_rates.tsv").readline().strip("# \n")
    rates_fig(prov); forest_fig(prov)


if __name__ == "__main__":
    main()
