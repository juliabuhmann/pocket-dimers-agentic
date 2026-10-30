#!/usr/bin/env python
"""Part 1 figures from results/part1/analysis/*.tsv (PNG, dpi 150). Status: UNVERIFIED until independent re-derivation.

    .venvs/part1/bin/python scripts/part1_figures.py [--in results/part1/analysis] [--out results/part1/figures]
Palette: dataviz skill reference categorical slots 1-2 (blue/orange) + text greys.
"""
import argparse
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DPI = 150
COL = {"homomer": "#2a78d6", "heteromer": "#eb6834"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
REF = {"Q-BioLiP 28.8%": 0.288, "PLINDER 25%": 0.25}   # LITERATURE.md; instance-level, all ligand classes
FINE_ORDER = ["cofactor", "heme_porphyrin", "nucleotide_like", "drug_like", "other_organic", "lipid_fatty_acid",
              "sugar_glycan", "amino_acid", "metal_cluster", "metal_ion", "inorganic"]
MARSH_ORDER = ["cofactor", "small_molecule", "metal"]


def read(p):
    return pd.read_csv(p, sep="\t", comment="#")


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def fig_fractions(d, out, level="unit_nonredundant", filt="after_additive_filter"):
    d = d[(d.level == level) & (d["filter"] == filt)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), gridspec_kw={"width_ratios": [1, 2.4]})
    for ax, grp, col, order, title in [(axes[0], "assembly_type+marsh_class", "marsh_class", MARSH_ORDER, "Marsh 3-class"),
                                       (axes[1], "assembly_type+cls", "cls", FINE_ORDER, "Fine class")]:
        g = d[d.grouping == grp]
        order = [o for o in order if o in set(g[col])]
        y = np.arange(len(order))[::-1]
        for k, (at, off) in enumerate([("homomer", 0.17), ("heteromer", -0.17)]):
            gg = g[g.assembly_type == at].set_index(col).reindex(order)
            ax.errorbar(gg.mbs_frac, y + off, xerr=[gg.mbs_frac - gg.ci_lo, gg.ci_hi - gg.mbs_frac], fmt="o", ms=6,
                        color=COL[at], ecolor=COL[at], elinewidth=1.5, capsize=0, label=at, mec="white", mew=1)
            for yi, (_, r) in zip(y, gg.iterrows()):
                if pd.notna(r.n):
                    ax.text(1.0, yi + off, f"n={int(r.n):,}", fontsize=6.5, color=INK2, va="center", ha="right")
        for (lab, v), ls in zip(REF.items(), ("--", ":")):
            ax.axvline(v, color=INK2, lw=1, ls=ls, label=lab if ax is axes[1] else None)
        ax.set_yticks(y)
        ax.set_yticklabels([o.replace("_", " ") for o in order], color=INK, fontsize=8.5)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("MBS fraction (95% cluster-bootstrap CI)", color=INK2, fontsize=8.5)
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        style(ax)
    axes[1].legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(0.92, 0))
    fig.suptitle(f"MBS fraction by ligand class and assembly type ({level.replace('_', ' ')}, {filt.replace('_', ' ')})"
                 " - UNVERIFIED", color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("wrote", out)


def fig_forest(t, out):
    tests = [("H1", "H1: adjusted OR of MBS,\ncofactor vs drug-like"), ("H2", "H2: conditional OR of drug-like\noccupancy, MBS vs SBS"),
             ("H3", "H3: OR drug-like vs cofactor,\non- vs off-axis homomer MBS")]
    extra = {"H3_geom_c2": "H3", "H3_union_jaccard_or_geom": "H3"}
    t = t.copy()
    t["panel"] = t.test.map(lambda x: extra.get(x, x))
    # H3 geometric/union variants are shown for the primary analysis only (all rows are in hypothesis_tests.tsv)
    t = t[~t.test.isin(list(extra)) | (t.variant == "primary")]
    t["label"] = np.where(t.test.isin(list(extra)), t.variant + " [" + t.test.str.replace("H3_", "").str.replace(
        "_jaccard_or_geom", " J|geom") + "]", t.variant)
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.5), sharey=False)
    for ax, (key, title) in zip(axes, tests):
        g = t[t.panel == key].reset_index(drop=True)
        y = np.arange(len(g))[::-1]
        for yi, (_, r) in zip(y, g.iterrows()):
            prim = r.variant == "primary" and r.test == key
            c = COL["homomer"] if prim else INK2
            filled = bool(r.power_gate_pass)
            if pd.notna(r.estimate):
                ax.plot([r.ci_lo, r.ci_hi], [yi, yi], color=c, lw=2 if prim else 1.2)
                ax.plot(r.estimate, yi, "o" if not prim else "D", ms=7 if prim else 5.5, color=c,
                        mfc=c if filled else "white", mec=c)
        ax.axvline(1, color=INK, lw=0.8)
        ax.set_xscale("log")
        lo = np.nanmin(g.ci_lo.values); hi = np.nanmax(g.ci_hi.values)
        ticks = [x for x in (0.125, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 8, 16, 32, 64, 128, 256) if lo / 1.5 <= x <= hi * 1.5]
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_yticks(y)
        ax.set_yticklabels(g.label.str.replace("_", " "), fontsize=7.5, color=INK)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("odds ratio (95% CI, log scale)", fontsize=8.5, color=INK2)
        style(ax)
    fig.suptitle("Family A hypotheses: primary (blue diamond) and sensitivities; hollow = power gate failed (descriptive)"
                 " - UNVERIFIED", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(ROOT / "results/part1/analysis"))
    ap.add_argument("--out", default=str(ROOT / "results/part1/figures"))
    a = ap.parse_args()
    inp, out = Path(a.inp), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    d = read(inp / "mbs_fraction_descriptives_seed0.tsv")
    fig_fractions(d, out / "mbs_fraction_by_class_units.png")
    fig_fractions(d, out / "mbs_fraction_by_class_sites.png", level="site_symcollapsed")
    fig_fractions(d, out / "mbs_fraction_by_class_units_before_additive.png", filt="before_additive_filter")
    fig_forest(read(inp / "hypothesis_tests.tsv"), out / "forest_h1_h3.png")


if __name__ == "__main__":
    main()
