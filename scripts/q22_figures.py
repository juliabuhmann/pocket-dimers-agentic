#!/usr/bin/env python
"""Figures for Q2.2 / Q2.3a (EXPLORATORY; PNG, dpi 150). Reads results/part2/q22/*.tsv and data/interim/q22/novelty.parquet.

Palette: validated reference categorical palette of the dataviz skill (slots 1-3: blue, orange, aqua), fixed order.
Status: first attempt, not yet verified (figures not yet visually inspected by the user).
"""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import OUT, FIG, INTERIM, WS, provenance

C = {"all": "#2a78d6", "homo": "#eb6834", "hetero": "#1baf7a"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "legend.frameon": False})
PV = provenance()
FOOT = f"EXPLORATORY · UNVERIFIED · design {PV['design_sha256'][:10]} · commit {PV['git_commit'][:8]}"


def rd(name):
    return pd.read_csv(OUT / name, sep="\t", comment="#")


def save(fig, name):
    fig.text(0.99, 1.04, FOOT, fontsize=6, color=INK2, ha="right", va="bottom")
    fig.savefig(FIG / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


LABEL = {"volume_A3": "volume (Å³)", "site_volume_A3": "site volume (Å³)", "burrow_score": "burrow_score",
         "plddt_mean": "lining pLDDT", "buriedness": "buriedness (1-rSASA)", "frac_hydrophobic": "hydrophobic frac.",
         "frac_aromatic": "aromatic frac.", "frac_charged": "charged frac.", "n_lining": "# lining residues",
         "max_excl": "max cosine to atlas (excl. self+paralogs)", "z_excl": "per-query z-score", "max_self": "max cosine to own atlas pockets"}


def forest(ax, df, title):
    df = df[(df.seed == 0) & df.stratum.isin(["all", "homo", "hetero"])]
    metrics = list(dict.fromkeys(df.metric))
    for k, s in enumerate(["all", "homo", "hetero"]):
        sub = df[df.stratum == s].set_index("metric").reindex(metrics)
        y = np.arange(len(metrics)) + (k - 1) * 0.22
        ax.errorbar(sub.cliffs_delta, y, xerr=[sub.cliffs_delta - sub.cd_lo, sub.cd_hi - sub.cliffs_delta],
                    fmt="o", ms=4, color=C[s], elinewidth=1.5, capsize=0, label=s)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(np.arange(len(metrics))); ax.set_yticklabels([LABEL.get(m, m) for m in metrics])
    ax.invert_yaxis(); ax.set_xlabel("Cliff's delta (95% cluster-bootstrap CI)"); ax.set_title(title, loc="left", fontsize=9)


# 1. property contrasts
p1, p2 = rd("q22_props.tsv"), rd("q22_props_dimer_only.tsv")
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
forest(axs[0], p1, "Interface vs single-chain pockets\n(>0: interface larger)")
forest(axs[1], p2, "Dimer-only vs pre-existing interface pockets\n(>0: dimer-only larger)")
axs[0].legend(loc="lower right", title=None)
save(fig, "q22_pocket_properties.png")

# 1b. ipSAE tertile forest (interface vs single-chain)
t = p1[p1.stratum.str.startswith("ipSAE")]
fig, ax = plt.subplots(figsize=(5.5, 4.2))
metrics = list(dict.fromkeys(t.metric))
for k, (s, col) in enumerate(zip(["ipSAE_T1_low", "ipSAE_T2_mid", "ipSAE_T3_high"], ["#86b6ef", "#3987e5", "#104281"])):
    sub = t[t.stratum == s].set_index("metric").reindex(metrics)
    y = np.arange(len(metrics)) + (k - 1) * 0.22
    ax.errorbar(sub.cliffs_delta, y, xerr=[sub.cliffs_delta - sub.cd_lo, sub.cd_hi - sub.cliffs_delta], fmt="o", ms=4,
                color=col, elinewidth=1.5, label=s.replace("ipSAE_", "ipSAE "))
ax.axvline(0, color=INK2, lw=0.8); ax.set_yticks(range(len(metrics))); ax.set_yticklabels([LABEL.get(m, m) for m in metrics])
ax.invert_yaxis(); ax.set_xlabel("Cliff's delta, interface vs single-chain"); ax.legend(loc="lower left")
ax.set_title("Interface vs single-chain by ipSAE tertile", loc="left", fontsize=9)
save(fig, "q22_props_by_ipsae.png")

# 2. expected vs observed
e = rd("q22_expected_rate.tsv")
qs = ["expected_iface_frac_of_surface_residues", "observed_iface_frac_of_pocket_lining_surface_residues",
      "observed_iface_frac_of_pockets_all_copies"]
ql = ["expected:\ninterface share of\nsurface residues", "observed:\ninterface share of\npocket-lining surface res.",
      "observed:\ninterface share\nof pockets"]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
for k, s in enumerate(["all", "homo", "hetero"]):
    sub = e[(e.stratum == s)].set_index("quantity").reindex(qs)
    x = np.arange(len(qs)) + (k - 1) * 0.26
    ax.bar(x, sub.estimate, width=0.24, color=C[s], label=s, edgecolor="#fcfcfb", linewidth=1)
    ax.errorbar(x, sub.estimate, yerr=[sub.estimate - sub.lo, sub.hi - sub.estimate], fmt="none", ecolor=INK, elinewidth=1)
ax.set_xticks(range(len(qs))); ax.set_xticklabels(ql, fontsize=8); ax.set_ylabel("fraction"); ax.legend()
ax.set_title("Interface pockets relative to interface surface area (atlas-pass pockets)", loc="left", fontsize=9)
save(fig, "q22_expected_vs_observed.png")

# 3. dimer-only fraction + destruction
f = rd("q22_dimer_only_frac.tsv"); de = rd("q22_destruction.tsv")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for k, s in enumerate(["all", "homo", "hetero"]):
    sub = f[f.stratum == s].sort_values("threshold_oc")
    x = np.arange(len(sub)) + (k - 1) * 0.26
    axs[0].bar(x, sub.frac, width=0.24, color=C[s], label=s, edgecolor="#fcfcfb")
    axs[0].errorbar(x, sub.frac, yerr=[sub.frac - sub.lo, sub.hi - sub.frac], fmt="none", ecolor=INK, elinewidth=1)
axs[0].set_xticks(range(3)); axs[0].set_xticklabels([f"OC < {t}" + (" (primary)" if t == 0.5 else "") for t in sorted(f.threshold_oc.unique())])
axs[0].set_ylabel("fraction of interface pockets that are dimer-only"); axs[0].legend()
axs[0].set_title("Dimer-only interface pockets (overlap-coefficient sensitivity)", loc="left", fontsize=9)
dd = de[de.stratum == "all"]
cts = ["partner_contact_ge3", "contact_1to2", "no_contact"]; ctl = ["≥3 lining res.\nat interface", "1-2 at\ninterface", "none at interface\n(control)"]
bottom = np.zeros(3)
for fate, col in [("lost", "#e34948"), ("shrunk", "#eda100"), ("retained", "#b7d3f6")]:
    v = dd[dd.fate == fate].set_index("contact").reindex(cts).frac.to_numpy()
    axs[1].bar(range(3), v, bottom=bottom, color=col, label=fate, edgecolor="#fcfcfb", linewidth=2, width=0.6)
    bottom += v
axs[1].set_xticks(range(3)); axs[1].set_xticklabels(ctl, fontsize=8); axs[1].set_ylabel("fraction of isolated-chain pockets")
axs[1].legend(loc="upper right", bbox_to_anchor=(1.28, 1)); axs[1].set_title("Fate of isolated-chain pockets in the dimer\n(no homodimer chain-swap; see q22_controls.png)", loc="left", fontsize=9)
save(fig, "q22_dimer_only_and_destruction.png")

# 4. novelty
nv = pd.read_parquet(INTERIM / "novelty.parquet")
S = nv[nv.atlas_pass & nv.is_sym_rep & ~(nv.is_interface_ge1 & ~nv.is_interface)].copy()
S["dec"] = pd.qcut(S.volume_A3.rank(method="first"), 10, labels=False)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for grp, col, lab in [(False, "#eb6834", "single-chain"), (True, "#2a78d6", "interface")]:
    g = S[S.is_interface == grp]
    xs = np.sort(g.max_excl.dropna()); axs[0].plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=col, lw=2, label=f"{lab} (n={len(xs)})")
    m = g.groupby("dec").max_excl.median(); axs[1].plot(m.index + 1, m.values, "-o", color=col, lw=2, ms=5, label=lab)
axs[0].set_xlabel("max cosine to human atlas (excl. self + ≥30% id paralogs)"); axs[0].set_ylabel("cumulative fraction"); axs[0].legend(loc="upper left")
axs[0].set_title("Similarity to monomer atlas pockets", loc="left", fontsize=9)
axs[1].set_xlabel("volume decile (1 = smallest)"); axs[1].set_ylabel("median max cosine"); axs[1].legend()
axs[1].set_title("Size-matched (volume deciles)\ncaveat: pocketcompare weak on PPI pockets", loc="left", fontsize=9)
save(fig, "q22_novelty.png")

# 5. Q2.3a ChEMBL status
q = pd.read_csv(WS / "results/part2/q23a_dimer_only_ligandable.tsv", sep="\t", comment="#")
q = q[q.in_atlas]
order = ["drug_target", "potent_ligands", "chembl_target", "not_in_chembl"]
v = q.chembl_status.value_counts().reindex(order).fillna(0)
fig, ax = plt.subplots(figsize=(5.5, 2.6))
ax.barh(range(4), v.values, color="#2a78d6", height=0.6); ax.set_yticks(range(4)); ax.set_yticklabels(order); ax.invert_yaxis()
for i, x in enumerate(v.values): ax.text(x, i, f" {int(x)}", va="center", color=INK, fontsize=8)
ax.set_xlabel("proteins"); ax.set_title(f"Q2.3a proteins with dimer-only drug-size pocket only (n={len(q)}): ChEMBL status", loc="left", fontsize=9)
save(fig, "q22_q23a_chembl_status.png")
print("figures written")

# 6. dimer-only null control + destruction control
ca, cb = rd("q22_dimer_only_control.tsv"), rd("q22_destruction_control.tsv")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
a = ca[ca.stratum == "all"].sort_values("threshold_oc")
x = np.arange(len(a))
axs[0].bar(x - 0.15, a.rate_g1, width=0.28, color="#2a78d6", label="interface pockets")
axs[0].bar(x + 0.15, a.rate_g0, width=0.28, color="#eb6834", label="single-chain pockets (null)")
for xi, (r1, r0, e, lo_, hi_) in enumerate(zip(a.rate_g1, a.rate_g0, a.score_matched_diff, a.sm_lo, a.sm_hi)):
    axs[0].text(xi, max(r1, r0) + 0.03, f"Δ {e:+.2f} [{lo_:.2f}, {hi_:.2f}]", ha="center", fontsize=7, color=INK)
axs[0].set_xticks(x); axs[0].set_xticklabels([f"OC < {t}" for t in a.threshold_oc]); axs[0].set_ylim(0, 1)
axs[0].set_ylabel("fraction without isolated-chain counterpart"); axs[0].legend(loc="upper left")
axs[0].set_title("'New in dimer' rate vs single-chain null\n(Δ = burrow_score-quintile-matched excess, 95% CI)", loc="left", fontsize=9)
b = cb.set_index("stratum").reindex(["all", "homo", "hetero"])
x = np.arange(3)
axs[1].bar(x - 0.15, b.rate_g1, width=0.28, color="#2a78d6", label="≥3 lining residues at interface")
axs[1].bar(x + 0.15, b.rate_g0, width=0.28, color="#eb6834", label="no interface contact (null)")
for xi, (r1, r0, e, lo_, hi_) in enumerate(zip(b.rate_g1, b.rate_g0, b.score_matched_diff, b.sm_lo, b.sm_hi)):
    axs[1].text(xi, max(r1, r0) + 0.03, f"Δ {e:+.2f} [{lo_:.2f}, {hi_:.2f}]", ha="center", fontsize=7, color=INK)
axs[1].set_xticks(x); axs[1].set_xticklabels(["all", "homo", "hetero"]); axs[1].set_ylim(0, 0.8)
axs[1].set_ylabel("fraction of isolated-chain pockets lost in dimer"); axs[1].legend(loc="upper left")
axs[1].set_title("Isolated-chain pocket loss (OC < 0.5, homodimer chain-swap matched)", loc="left", fontsize=9)
save(fig, "q22_controls.png")
print("control figure written")
