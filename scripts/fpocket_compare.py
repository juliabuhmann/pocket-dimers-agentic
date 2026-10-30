#!/usr/bin/env python
"""MBS vs SBS fpocket druggability at PDB ligand sites (EXPLORATORY, not pre-registered; plan declared in
docs/FPOCKET_DRUGGABILITY.md before this was run).

Inputs: data/interim/fpocket/pdbsite/part_*.parquet (scripts/fpocket_pdbsite.py), site labels and weights from
data/interim/pdbsite/sites_crop30.parquet, subsample from results/part1/fpocket/pdbsite_subsample_seed0.tsv,
stratum sizes from results/part1/pocket_descriptors/sample_strata_n12000_seed0.tsv.
Statistics follow scripts/pdbsite_compare.py (weighted Cliff's delta, cluster30 bootstrap, spline OLS, within-cluster
FE and pair delta); the helper functions are copied from there so this script does not depend on a file another
analysis is editing.
Outputs: results/part1/fpocket/pdbsite_*.tsv, data/interim/fpocket/pdbsite_sites.parquet
Status: first attempt, not yet verified.
"""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fpocket_common import WS, OUT, INTERIM, DRUG_THR, fpocket_version

SEED = 0
N_BOOT = 1000          # task spec / DESIGN.md §5
SPLINE_DF = 4          # same as scripts/pdbsite_common.py
MIN_COVER = 0.5        # sensitivity: chosen fpocket pocket covers >= 50% of ligand heavy atoms (our choice)
METRICS = {"fp_drug": "PRIMARY (exploratory): fpocket druggability score", "fp_score": "fpocket score",
           "fp_volume": "fpocket volume (A^3)", "fp_hydrophobicity": "fpocket hydrophobicity score",
           "fp_polarity": "fpocket polarity score", "drug_gt_thr": f"indicator fp_drug > {DRUG_THR}"}


def wcliff(x, y, wx=None, wy=None):
    """Weighted Cliff's delta (copied from scripts/pdbsite_compare.py)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    wx = np.ones(len(x)) if wx is None else np.asarray(wx, float)
    wy = np.ones(len(y)) if wy is None else np.asarray(wy, float)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    o = np.argsort(y); ys, wys = y[o], wy[o]
    cw = np.concatenate([[0], np.cumsum(wys)])
    lo = np.searchsorted(ys, x, "left"); hi = np.searchsorted(ys, x, "right")
    return float((wx * (cw[lo] - (cw[-1] - cw[hi]))).sum() / (wx.sum() * wys.sum()))


def within_cliff(d, col):
    """Pooled within-cluster Cliff's delta over MBS-SBS pairs in the same cluster30 (copied from pdbsite_compare.py)."""
    r = d.groupby("cluster30")[col].rank().values
    g1 = (d.g.values == 1).astype(float)
    a = pd.DataFrame({"c": d.cluster30.values, "one": 1.0, "g1": g1, "r1": r * g1}).groupby("c").sum()
    n1, n0 = a.g1.values, a.one.values - a.g1.values
    u = a.r1.values - n1 * (n1 + 1) / 2
    num, den = (2 * u - n1 * n0).sum(), (n1 * n0).sum()
    return num / den if den else np.nan


def boot(d, stat, seed=SEED, n_boot=None):
    n_boot = n_boot or N_BOOT
    rng = np.random.default_rng(seed)
    gi = d.groupby("cluster30").indices
    keys = list(gi); idx = [gi[k] for k in keys]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(keys), len(keys))
        dd = d.iloc[np.concatenate([idx[i] for i in pick])].copy()
        dd["cluster30"] = np.repeat(np.arange(len(pick)), [len(idx[i]) for i in pick])
        vals.append(stat(dd))
    return np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)


def wmean(v, w):
    return float(np.average(v, weights=w)) if len(v) else np.nan


def write(df, name, note=""):
    with open(OUT / name, "w") as fh:
        fh.write(f"# status=UNVERIFIED EXPLORATORY (not pre-registered) {fpocket_version()} {note}\n")
        df.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    print("wrote", name, len(df), flush=True)


def main():
    global N_BOOT, OUT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=N_BOOT, help="debug only; reported results use N_BOOT")
    ap.add_argument("--out-suffix", default="", help="debug only")
    a = ap.parse_args()
    N_BOOT = a.n_boot
    if a.out_suffix:
        OUT = OUT / a.out_suffix; OUT.mkdir(exist_ok=True)
    lab = pd.read_parquet(WS / "data/interim/pdbsite/sites_crop30.parquet",
                          columns=["site_id", "pdb_id", "ccd", "cluster30", "stratum", "stratum_lig", "y_mbs",
                                   "assembly_type", "is_on_axis", "n_heavy_site", "log_heavy"])
    sub = pd.read_csv(WS / "results/part1/fpocket/pdbsite_subsample_seed0.tsv", sep="\t", comment="#")
    fp = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(INTERIM / "pdbsite/part_*.parquet")))],
                   ignore_index=True).drop_duplicates("site_id")
    d = lab[lab.site_id.isin(sub.site_id)].merge(fp, on="site_id", how="left")
    strata = pd.read_csv(WS / "results/part1/pocket_descriptors/sample_strata_n12000_seed0.tsv", sep="\t", comment="#")
    avail = dict(zip(strata.stratum, strata.units_available))
    # ledger
    led = [("subsample (500 per class x MBS stratum, seed 0)", len(d), d.cluster30.nunique()),
           ("fpocket run returned", int(d.status.notna().sum()), d[d.status.notna()].cluster30.nunique())]
    for st, g in d.groupby(d.status.fillna("not_run")):
        led.append((f"status {st}", len(g), g.cluster30.nunique()))
    for e, n in d[d.status == "error"].error.str.split(":").str[0].value_counts().items():
        led.append((f"error: {e}", int(n), np.nan))
    ok = d[d.status == "ok"].copy()
    led.append((f"chosen pocket covers >= {MIN_COVER:.0%} of ligand heavy atoms", int((ok.fp_overlap_frac >= MIN_COVER).sum()),
                ok[ok.fp_overlap_frac >= MIN_COVER].cluster30.nunique()))
    write(pd.DataFrame(led, columns=["step", "n_sites", "n_clusters"]), "pdbsite_run_ledger.tsv")
    # weights: population units per analysed site in the stratum
    nst = ok.groupby("stratum").size()
    ok["w"] = ok.stratum.map(lambda s: avail[s] / nst[s])
    ok["drug_gt_thr"] = (ok.fp_drug > DRUG_THR).astype(float)
    ok["g"] = ok.y_mbs.astype(int)
    ok.to_parquet(INTERIM / "pdbsite_sites.parquet", index=False)

    # ---- Cliff's delta + fraction > threshold
    rows = []
    sets = {"all_weighted": (ok, True), "all_unweighted": (ok, False),
            "all_weighted_cover50": (ok[ok.fp_overlap_frac >= MIN_COVER], True),
            "drug_like": (ok[ok.stratum_lig == "drug_like"], False), "cofactor": (ok[ok.stratum_lig == "cofactor"], False),
            "other": (ok[ok.stratum_lig == "other"], False),
            "homomer": (ok[ok.assembly_type == "homomer"], True), "heteromer": (ok[ok.assembly_type == "heteromer"], True)}
    for m in METRICS:
        for name, (x, wt) in sets.items():
            x = x.dropna(subset=[m])
            def st(dd, m=m, wt=wt):
                a, b = dd[dd.g == 1], dd[dd.g == 0]
                return wcliff(a[m], b[m], a.w if wt else None, b.w if wt else None)
            def fr(dd, wt=wt):
                a, b = dd[dd.g == 1], dd[dd.g == 0]
                f = (lambda h: wmean(h.drug_gt_thr, h.w)) if wt else (lambda h: h.drug_gt_thr.mean())
                return f(a) - f(b)
            cols = list(dict.fromkeys(["cluster30", "g", "w", m, "drug_gt_thr"]))
            est = st(x); lo, hi = boot(x[cols], st)
            a, b = x[x.g == 1], x[x.g == 0]
            r = dict(metric=m, role=METRICS[m], stratum=name, weighted=wt, n_mbs=len(a), n_sbs=len(b),
                     n_clusters=x.cluster30.nunique(), median_mbs=a[m].median(), median_sbs=b[m].median(),
                     cliffs_delta=est, cd_lo=lo, cd_hi=hi)
            if m == "fp_drug":
                f1 = wmean(a.drug_gt_thr, a.w) if wt else a.drug_gt_thr.mean()
                f0 = wmean(b.drug_gt_thr, b.w) if wt else b.drug_gt_thr.mean()
                flo, fhi = boot(x[cols], fr)
                r.update(frac_gt_thr_mbs=f1, frac_gt_thr_sbs=f0, frac_diff=f1 - f0, frac_diff_lo=flo, frac_diff_hi=fhi)
            rows.append(r)
        print("delta", m, flush=True)
    write(pd.DataFrame(rows), "pdbsite_mbs_vs_sbs_cliffs_delta.tsv",
          note=f"MBS vs SBS; cluster30 bootstrap {N_BOOT} seed {SEED}; weights = units_available/n_analysed per stratum")

    # ---- regression (raw, percentile rank) + within-cluster FE + pair delta
    reg = []
    for m in METRICS:
        x = ok.dropna(subset=[m]).copy()
        x["y"] = x[m]; x["yp"] = x[m].rank(pct=True)
        for scale, yv in (("raw", "y"), ("percentile", "yp")):
            f = f"{yv} ~ g + bs(log_heavy, df={SPLINE_DF}) + C(stratum_lig)"
            fit = smf.ols(f, x).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(x.cluster30)[0]})
            ci = fit.conf_int().loc["g"]
            reg.append(dict(metric=m, model="pooled OLS spline(log heavy)+class, cluster-robust", scale=scale, n=len(x),
                            n_clusters=x.cluster30.nunique(), beta_mbs=fit.params["g"], lo=ci[0], hi=ci[1], p=fit.pvalues["g"]))
        both = x.groupby("cluster30").g.transform(lambda s: s.nunique() == 2)
        xb = x[both].copy()
        if xb.cluster30.nunique() >= 10:
            for scale, yv in (("raw", "y"), ("percentile", "yp")):
                f = f"{yv} ~ g + bs(log_heavy, df={SPLINE_DF}) + C(stratum_lig) + C(cluster30)"
                fit = smf.ols(f, xb).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(xb.cluster30)[0]})
                ci = fit.conf_int().loc["g"]
                reg.append(dict(metric=m, model="within-cluster FE (clusters with MBS and SBS)", scale=scale, n=len(xb),
                                n_clusters=xb.cluster30.nunique(), beta_mbs=fit.params["g"], lo=ci[0], hi=ci[1], p=fit.pvalues["g"]))
            st = lambda dd, m=m: within_cliff(dd, m)  # noqa: E731
            lo, hi = boot(xb[["cluster30", "g", m]], st)
            top5 = xb.groupby("cluster30").apply(lambda h: (h.g == 1).sum() * (h.g == 0).sum()).nlargest(5).sum()
            npairs = xb.groupby("cluster30").apply(lambda h: (h.g == 1).sum() * (h.g == 0).sum()).sum()
            reg.append(dict(metric=m, model="within-cluster pair Cliff's delta", scale="rank", n=len(xb),
                            n_clusters=xb.cluster30.nunique(), beta_mbs=st(xb), lo=lo, hi=hi, p=np.nan,
                            note=f"top-5 clusters supply {top5 / npairs:.0%} of pairs"))
    write(pd.DataFrame(reg), "pdbsite_mbs_vs_sbs_regression.tsv", note="g = MBS indicator")

    # ---- on-axis vs off-axis homomer MBS
    h = ok[(ok.assembly_type == "homomer") & (ok.g == 1) & ok.is_on_axis.notna()].copy()
    h["g"] = h.is_on_axis.astype(bool).astype(int)
    rows = []
    for m in METRICS:
        x = h.dropna(subset=[m])
        st = lambda dd, m=m: wcliff(dd[dd.g == 1][m], dd[dd.g == 0][m])  # noqa: E731
        lo, hi = boot(x[["cluster30", "g", m]], st)
        a, b = x[x.g == 1], x[x.g == 0]
        top2 = a.cluster30.value_counts().head(2)
        rows.append(dict(metric=m, n_on=len(a), n_off=len(b), clusters_on=a.cluster30.nunique(),
                         median_on=a[m].median(), median_off=b[m].median(), cliffs_delta=st(x), cd_lo=lo, cd_hi=hi,
                         frac_gt_thr_on=(a.fp_drug > DRUG_THR).mean(), frac_gt_thr_off=(b.fp_drug > DRUG_THR).mean(),
                         on_top2_clusters=";".join(f"{k}:{v}" for k, v in top2.items())))
    write(pd.DataFrame(rows), "pdbsite_onaxis_vs_offaxis_mbs.tsv", note="homomer MBS only, unweighted")


if __name__ == "__main__":
    main()
