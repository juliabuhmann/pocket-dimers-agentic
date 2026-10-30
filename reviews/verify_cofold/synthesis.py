"""Exploratory synthesis: partner gain vs partner share of the PDB site (plan pre-declared in docs/COFOLD_SYNTHESIS.md).

First attempt; not verified by a second party. Inputs are the verifier's own re-scored tables:
  reviews/verify_cofold/rescored_seed0.tsv, reviews/verify_cofold/rescored_seeds12.tsv (if present),
  reviews/verify_cofold/site_share.tsv
Outputs: results/cofold/synthesis/{cases,stats}.tsv, results/cofold/synthesis/partner_gain_vs_share.png
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / "reviews/verify_cofold"
OUT = ROOT / "results/cofold/synthesis"
N_BOOT = 2000            # pre-declared
RNG_SEED = 0
SUCCESS_CUT = 2.0
# Sensitivity exclusions (pre-declared): ligand-induced dimers and probably non-biological 2-chain assemblies
EXCLUDE = ["on_8cza_ZN1", "on_7mr7_ZN1", "on_4mtz_0BU", "sbs_8ivf_HC3", "sbs_7era_JBC", "sbs_8ii3_PKK"]
# dataviz reference palette, categorical slots 1-5 in fixed order (validated with validate_palette.js, light mode)
CLASS_STYLE = {
    "on-axis drug (onaxis)": ("#2a78d6", "o"),
    "off-axis cofactor (onaxis)": ("#eb6834", "s"),
    "SBS drug (onaxis)": ("#1baf7a", "^"),
    "dimer-only pocket": ("#eda100", "D"),
    "H4 MBS site": ("#e87ba4", "v"),
}
ONAXIS_CLASS = {"on_axis_druglike": "on-axis drug (onaxis)", "offaxis_cofactor": "off-axis cofactor (onaxis)",
                "sbs_druglike": "SBS drug (onaxis)"}


def case_table(resc):
    d = resc[resc.cond == "dimer"].set_index(["batch", "case_id"])
    m = resc[resc.cond == "monomer"].set_index(["batch", "case_id"])
    agg = lambda x: x.groupby(level=[0, 1]).agg(rmsd=("lig_rmsd", "mean"), logr=("logr", "mean"),
                                                 succ=("success", "mean"), n_seeds=("seed", "nunique"))
    d, m = agg(d), agg(m)
    t = d.join(m, lsuffix="_dimer", rsuffix="_single", how="inner").reset_index()
    t["gain_rmsd"] = t.rmsd_single - t.rmsd_dimer
    t["gain_log10"] = t.logr_single - t.logr_dimer
    t["success_diff"] = t.succ_dimer - t.succ_single
    t["dimer_only"] = ((t.succ_dimer > 0.5) & (t.succ_single <= 0.5)).astype(int)
    return t


def stats(t, label):
    rng = np.random.default_rng(RNG_SEED)
    x = t.partner_share.to_numpy()
    out = []
    for y in ["gain_rmsd", "gain_log10", "success_diff"]:
        rho = spearmanr(x, t[y]).statistic
        bs = []
        for _ in range(N_BOOT):
            i = rng.integers(0, len(t), len(t))
            if np.ptp(x[i]) == 0 or np.ptp(t[y].to_numpy()[i]) == 0:
                continue
            bs.append(spearmanr(x[i], t[y].to_numpy()[i]).statistic)
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        out.append(dict(subset=label, n=len(t), test=f"spearman {y}~partner_share", estimate=rho, ci_lo=lo, ci_hi=hi,
                        p_one_sided_boot=float(np.mean(np.array(bs) <= 0)), n_boot_ok=len(bs)))
    y = t.dimer_only.to_numpy()

    def fit(xx, yy):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            r = sm.Logit(yy, sm.add_constant(xx)).fit(disp=0, maxiter=200)
        return r.params[1]
    try:
        b = fit(x, y)
    except Exception:
        b = np.nan
    bs, bad = [], 0
    for _ in range(N_BOOT):
        i = rng.integers(0, len(t), len(t))
        try:
            bs.append(fit(x[i], y[i]))
        except Exception:
            bad += 1
    lo, hi = (np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan))
    out.append(dict(subset=label, n=len(t), test="logit dimer_only~partner_share (slope per unit share)",
                    estimate=b, ci_lo=lo, ci_hi=hi,
                    p_one_sided_boot=float(np.mean(np.array(bs) <= 0)) if bs else np.nan,
                    n_boot_ok=len(bs), n_boot_failed=bad, n_events=int(y.sum())))
    return out


def plot(t, path, title):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), dpi=150)
    for ax, y, yl in [(axes[0], "gain_rmsd", "Partner gain: single − dimer ligand RMSD (Å)"),
                      (axes[1], "gain_log10", "Partner gain: log10(single / dimer RMSD)")]:
        ax.axhline(0, color="#52514e", lw=0.8, zorder=0)
        for cls, (col, mk) in CLASS_STYLE.items():
            s = t[t.cls == cls]
            ax.scatter(s.partner_share, s[y], c=col, marker=mk, s=46, edgecolors="#fcfcfb", linewidths=1.2,
                       label=f"{cls} (n={len(s)})", zorder=3)
        for _, r in t.iterrows():
            if (y == "gain_rmsd" and abs(r[y]) > 10) or (y == "gain_log10" and (abs(r[y]) > 0.6 or r.dimer_only)):
                ax.annotate(r.label, (r.partner_share, r[y]), fontsize=6.5, color="#52514e",
                            xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("Partner share of PDB site (minor-chain residues / all)")
        ax.set_ylabel(yl)
        ax.set_xlim(-0.03, 0.53)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].set_yscale("symlog", linthresh=2)
    axes[1].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    r0 = pd.read_csv(V / "rescored_seed0.tsv", sep="\t")
    parts = [r0]
    f12 = V / "rescored_seeds12.tsv"
    if f12.exists():
        parts.append(pd.read_csv(f12, sep="\t"))
    resc = pd.concat(parts)
    resc = resc[resc.status == "ok"].copy()
    resc["logr"] = np.log10(resc.lig_rmsd.clip(lower=0.1))
    share = pd.read_csv(V / "site_share.tsv", sep="\t")
    strata = {}
    for b in ["onaxis", "dimeronly", "h4"]:
        c = pd.read_csv(ROOT / f"results/cofold/{b}/cases.tsv", sep="\t")
        for _, r in c.iterrows():
            strata[(b, r.case_id)] = (ONAXIS_CLASS[r.stratum] if b == "onaxis" else
                                      "dimer-only pocket" if b == "dimeronly" else "H4 MBS site")
    stats_rows = []
    for tag, sub in [("seed0", resc[resc.seed == 0]), ("seeds012_pooled", resc)]:
        t = case_table(sub).merge(share, on=["batch", "case_id"])
        if tag == "seeds012_pooled":
            t = t[(t.n_seeds_dimer > 1) & (t.n_seeds_single > 1)]
            if len(t) < 5:
                continue
        t["cls"] = [strata[(b, c)] for b, c in zip(t.batch, t.case_id)]
        t["label"] = t.case_id.str.replace("^(on_|offaxis_|sbs_|homo_miss_|homo_hit_|hetero_)", "", regex=True)
        t.to_csv(OUT / f"cases_{tag}.tsv", sep="\t", index=False, float_format="%.3f")
        stats_rows += stats(t, f"{tag}:all")
        if tag == "seed0":
            for b in ["onaxis", "dimeronly", "h4"]:
                stats_rows += stats(t[t.batch == b], f"{tag}:{b}")
            stats_rows += stats(t[~t.case_id.isin(EXCLUDE)], f"{tag}:excl_induced_nonbio")
            # POST HOC: one row per PDB entry+ligand (8ros FMN appears in onaxis and h4)
            stats_rows += stats(t.drop_duplicates(["pdb_id", "ccd"]), f"{tag}:dedup_pdb_ccd (post hoc)")
            plot(t, OUT / "partner_gain_vs_share.png",
                 f"AF3 partner gain vs partner share, seed 0 top-ranked, n={len(t)} cases (exploratory)")
        else:
            plot(t, OUT / "partner_gain_vs_share_seeds012.png",
                 f"AF3 partner gain vs partner share, seeds 0-2 pooled, n={len(t)} cases (exploratory)")
    pd.DataFrame(stats_rows).to_csv(OUT / "stats.tsv", sep="\t", index=False, float_format="%.3f")
    print(pd.DataFrame(stats_rows).to_string())


if __name__ == "__main__":
    main()
