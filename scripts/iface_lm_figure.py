#!/usr/bin/env python
"""Figure for the EXPLORATORY ligand-matched interface pocket analysis (reads results TSVs + oof.parquet)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from q22_common import WS, INTERIM

R = WS / "results/part2/iface_ligand_matched"
rd = lambda f: pd.read_csv(R / f, sep="\t", comment="#")
med, dl, auc = rd("group_medians.tsv"), rd("cliffs_delta.tsv"), rd("classifier_auc.tsv")
oof = pd.read_parquet(INTERIM.parent / "iface_ligand_matched/oof.parquet")
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
COL = {"IM": "#2a78d6", "IU": "#eb6834", "SM": "#1baf7a"}  # categorical slots 1-3
LAB = {"IM": "interface, ligand-matched", "IU": "interface, unmatched", "SM": "single-chain, ligand-matched"}
M = [("frac_charged", "charged fraction (DEKR)"), ("frac_hydrophobic", "hydrophobic fraction"), ("frac_aromatic", "aromatic fraction"),
     ("volume_A3", "volume (A^3)"), ("burrow_score", "burrow_score"), ("buriedness", "buriedness")]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
fig, axs = plt.subplots(2, 4, figsize=(16, 7.5), facecolor="#fcfcfb", constrained_layout=True)
def clean(ax):
    ax.set_facecolor("#fcfcfb"); ax.grid(axis="x", color=GRID, lw=.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
for ax, (m, t) in zip(axs.flat[:6], M):
    clean(ax)
    for i, g in enumerate(["IM", "IU", "SM"]):
        r = med[(med.group == g) & (med.metric == m)].iloc[0]
        ax.errorbar(r["median"], i, xerr=[[r["median"] - r.lo], [r.hi - r["median"]]], fmt="o", color=COL[g], ms=7, mec="#fcfcfb", mew=2, capsize=0, lw=2)
    ax.set_yticks(range(3)); ax.set_yticklabels([LAB[g] for g in ["IM", "IU", "SM"]] if m in ("frac_charged", "volume_A3") else [""] * 3)
    ax.set_ylim(2.6, -.6); ax.set_title(t, loc="left", fontsize=10); ax.set_xlabel("median, 95% CI")
ax = axs.flat[6]; clean(ax)
for j, (a, b, c) in enumerate([("IM", "IU", "#2a78d6"), ("IM", "SM", "#1baf7a"), ("IU", "SM", "#eb6834")]):
    for i, (m, t) in enumerate(M):
        r = dl[(dl.g1 == a) & (dl.g0 == b) & (dl.metric == m)].iloc[0]; y = i + (j - 1) * .22
        ax.errorbar(r.delta, y, xerr=[[r.delta - r.lo], [r.hi - r.delta]], fmt="o", color=c, ms=5, mec="#fcfcfb", mew=1.5, lw=1.5, label=f"{a} vs {b}" if i == 0 else None)
ax.axvline(0, color=MUTED, lw=1); ax.set_yticks(range(6)); ax.set_yticklabels([t.split(" (")[0] for _, t in M]); ax.set_ylim(5.6, -.6)
ax.set_xlabel("Cliff delta, 95% CI"); ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(.4, -.2), ncol=3); ax.set_title("effect sizes", loc="left", fontsize=10)
ax = axs.flat[7]; ax.set_facecolor("#fcfcfb")
from scipy.stats import rankdata
for key, c, ls in [("frac_charged+burrow_score", "#2a78d6", "-"), ("frac_charged", "#eb6834", "-"), ("burrow_score", "#1baf7a", "-"), ("volume_A3", "#52514e", "--")]:
    s, y = oof[key + "_oof"].to_numpy(), oof.y.to_numpy(); o = np.argsort(-s); yy = y[o]
    ax.plot(np.r_[0, np.cumsum(1 - yy) / (1 - y).sum()], np.r_[0, np.cumsum(yy) / y.sum()], color=c, ls=ls, lw=1.8,
            label=f"{key.replace('_', ' ')}  AUC {auc[auc.features == key].auc_oof.iloc[0]:.2f}")
ax.plot([0, 1], [0, 1], color=GRID, lw=1); ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(.5, -.2))
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate"); ax.set_title("IM vs IU classifier (OOF)", loc="left", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.suptitle("EXPLORATORY, unverified: predicted interface pockets in H4-benchmark dimer models, ligand-matched vs unmatched", x=.01, ha="left", fontsize=11)
fig.savefig(R / "iface_ligand_matched.png", dpi=150)
