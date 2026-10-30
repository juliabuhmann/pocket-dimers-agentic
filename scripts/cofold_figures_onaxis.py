#!/usr/bin/env python3
"""Figures for deep dive #1 co-folding (results/cofold/onaxis/metrics*.tsv). Status: first attempt, UNVERIFIED.
-> results/cofold/onaxis/figures/{lig_rmsd_by_case.png, success_rates.png}  (PNG, dpi 150)
Colors: dataviz-skill reference palette, categorical slots 1-2 (blue = dimer, orange = monomer)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/cofold/onaxis"
FIG = R / "figures"
SUCCESS_RMSD = 2.0
COL = {"dimer": "#2a78d6", "monomer": "#eb6834"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ORDER = ["on_axis_druglike", "offaxis_cofactor", "sbs_druglike"]
LABEL = {"on_axis_druglike": "On-axis drug-like", "offaxis_cofactor": "Off-axis MBS cofactor", "sbs_druglike": "SBS drug-like"}


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(R / "metrics.tsv", sep="\t")
    s = pd.read_csv(R / "metrics_samples.tsv", sep="\t")
    s = s[(s.status == "ok") & s.cond.isin(list(COL))]
    m = m[m.lig_rmsd.notna() & m.cond.isin(list(COL))].copy()   # seed-1 extra runs are reported in the tables only
    m["o"] = m.stratum.map({k: i for i, k in enumerate(ORDER)})
    cases = m.drop_duplicates("case_id").sort_values(["o", "post_cutoff", "case_id"])
    y = {c: i for i, c in enumerate(cases.case_id[::-1])}
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(y) + 1.4))
    for cond, dy in (("dimer", 0.15), ("monomer", -0.15)):
        ss = s[s.cond == cond]
        for cid, g in ss.groupby("case_id"):
            if cid in y:
                ax.plot([g.lig_rmsd.min(), g.lig_rmsd.max()], [y[cid] + dy] * 2, color=COL[cond], lw=1.2, alpha=0.5)
        mm = m[m.cond == cond]
        ax.scatter(mm.lig_rmsd, [y[c] + dy for c in mm.case_id], s=36, color=COL[cond], edgecolor="white", lw=1,
                   zorder=3, label=f"{cond} (top-ranked; line = range of 5 samples)")
    ax.axvline(SUCCESS_RMSD, color=INK2, lw=1, ls="--")
    ax.set_xscale("log")
    ticks = [0.25, 0.5, 1, 2, 5, 10, 20, 40]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels([f"{c}{' *' if cases.set_index('case_id').post_cutoff[c] else ''}" for c in y], fontsize=7.5)
    # stratum separators
    prev = None
    for c in cases.case_id:
        st = cases.set_index("case_id").stratum[c]
        if st != prev:
            ax.text(1.0, y[c] + 0.45, LABEL[st], transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                    fontsize=8, color=INK, weight="bold")
            if prev is not None:
                ax.axhline(y[c] + 0.5, color=GRID, lw=1)
            prev = st
    ax.set_xlabel("Ligand RMSD vs PDB pose (Å, log scale; dashed = 2 Å)", color=INK, fontsize=9)
    ax.set_title("AF3 ligand placement: homodimer vs monomer input (* = released after 2021-09-30)",
                 fontsize=9.5, color=INK, loc="left")
    style(ax)
    ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "lig_rmsd_by_case.png", dpi=150)
    plt.close(fig)

    # success rates by stratum x condition
    t = m.groupby(["stratum", "cond"]).agg(k=("success", "sum"), n=("success", "size")).reset_index()
    fig, ax = plt.subplots(figsize=(6, 2.8))
    for i, st in enumerate(ORDER):
        for cond, dy in (("dimer", 0.18), ("monomer", -0.18)):
            r = t[(t.stratum == st) & (t.cond == cond)]
            if r.empty:
                continue
            k, n = int(r.k.iloc[0]), int(r.n.iloc[0])
            ax.barh(len(ORDER) - 1 - i + dy, k / n, height=0.32, color=COL[cond], label=cond if i == 0 else None)
            ax.text(k / n + 0.01, len(ORDER) - 1 - i + dy, f"{k}/{n}", va="center", fontsize=8, color=INK2)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([LABEL[x] for x in ORDER[::-1]], fontsize=8.5)
    ax.set_xlim(0, 1.1)
    ax.set_xlabel("Fraction of cases with top-ranked ligand RMSD ≤ 2 Å", fontsize=9, color=INK)
    style(ax)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "success_rates.png", dpi=150)
    plt.close(fig)
    print("wrote", FIG)


if __name__ == "__main__":
    main()
