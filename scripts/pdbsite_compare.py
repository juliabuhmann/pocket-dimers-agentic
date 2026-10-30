#!/usr/bin/env python
"""MBS vs SBS pocket descriptors at real PDB ligand sites (H2 secondary descriptor = burrow_score of the burrow pocket at
the ligand; everything else EXPLORATORY), and comparison with the predicted-dimer pattern of Q2.2.

Inputs: data/interim/pdbsite/sample_n*_seed0.parquet (scripts/pdbsite_sample.py),
        data/interim/pdbsite/describe/part_*.parquet (scripts/pdbsite_describe.py), results/part2/q22/q22_props.tsv.
Analyses (unit = one site per cluster30 x CCD unit; bootstrap unit = cluster30; 1,000 resamples):
 1. Cliff's delta MBS vs SBS: pooled with stratum weights (units available / sampled, so the pooled contrast reflects
    the population of organic oligomer units, not the equal-allocation sample), and per ligand stratum (unweighted).
 2. Adjusted: OLS  y ~ MBS + B-spline(log heavy atoms, df 4) + ligand stratum, cluster-robust SEs; y raw and as a
    mid-rank percentile (rank-based, scale-free). Within-cluster: same model with cluster30 fixed effects (demeaned),
    restricted to clusters with both MBS and SBS sites.
 3. Within-cluster Cliff's delta: all MBS-SBS pairs inside the same cluster30.
 4. On-axis vs off-axis MBS (homomer MBS sites, census is_on_axis of the sampled site).
 5. Side-by-side with the Q2.2 predicted-dimer deltas (interface vs single-chain pockets).
Outputs: results/part1/pocket_descriptors/*.tsv, data/part2/pdb_site_fingerprints.parquet, data/interim/pdbsite/sites.parquet
Run with .venvs/part1/bin/python.  Status: first attempt, not yet verified.
"""
import argparse, glob, sys
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrix
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdbsite_common import (BIOLIP_MIN_PAIRS, INTERIM, OUT, FP_OUT, Q22_PROPS, SEEDS, N_BOOT, SPLINE_DF, POWER_GATE_MIN_CLUSTERS,
                            write_tsv)

# metric -> (column, role). Lining chemistry / volume / burrow_score of the MERGED pocket containing the located pocket
# are the like-for-like counterparts of the Q2.2 atlas-style pockets; the located-member burrow_score is the
# pre-registered H2 secondary descriptor.
METRICS = {
    "burrow_score": ("burrow_score", "H2-secondary (pre-registered)"),
    "volume_A3": ("volume_A3", "exploratory (located member)"),
    "merged_burrow_score": ("merged_burrow_score", "exploratory (merged pocket, Q2.2 unit)"),
    "merged_volume_A3": ("merged_volume_A3", "exploratory (merged pocket, Q2.2 unit)"),
    "site_volume_A3": ("site_volume_A3", "exploratory"),
    "buriedness": ("merged_buriedness", "exploratory (merged pocket, Q2.2 unit)"),
    "n_lining": ("merged_n_lining", "exploratory (merged pocket, Q2.2 unit)"),
    "frac_hydrophobic": ("merged_frac_hydrophobic", "exploratory (merged pocket, Q2.2 unit)"),
    "frac_aromatic": ("merged_frac_aromatic", "exploratory (merged pocket, Q2.2 unit)"),
    "frac_charged": ("merged_frac_charged", "exploratory (merged pocket, Q2.2 unit)"),
    "member_frac_hydrophobic": ("frac_hydrophobic", "exploratory (located member)"),
    "member_frac_charged": ("frac_charged", "exploratory (located member)"),
    "member_frac_aromatic": ("frac_aromatic", "exploratory (located member)"),
    # census `contacts` = the 4 Å set (every residue with n_simple >= 1); BioLiP rule = n_vdw >= 2 (vdW+0.5 Å pairs)
    "contact4A_frac_hydrophobic": ("contact4A_frac_hydrophobic", "exploratory (census contacts, 4 Å set)"),
    "contact4A_frac_charged": ("contact4A_frac_charged", "exploratory (census contacts, 4 Å set)"),
    "contact4A_frac_aromatic": ("contact4A_frac_aromatic", "exploratory (census contacts, 4 Å set)"),
    "contactBioLiP_frac_hydrophobic": ("contactBioLiP_frac_hydrophobic", "exploratory (census contacts, vdW+0.5 Å >=2 pairs)"),
    "contactBioLiP_frac_charged": ("contactBioLiP_frac_charged", "exploratory (census contacts, vdW+0.5 Å >=2 pairs)"),
    "contactBioLiP_frac_aromatic": ("contactBioLiP_frac_aromatic", "exploratory (census contacts, vdW+0.5 Å >=2 pairs)"),
}
Q22_MAP = {"volume_A3": "merged_volume_A3", "burrow_score": "merged_burrow_score", "buriedness": "buriedness",
           "n_lining": "n_lining", "frac_hydrophobic": "frac_hydrophobic", "frac_aromatic": "frac_aromatic",
           "frac_charged": "frac_charged"}   # Q2.2 metric name -> our metric name
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
             "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V", "MSE": "M"}
from pdbsite_common import HYDROPHOBIC, AROMATIC, CHARGED


def wcliff(x, y, wx=None, wy=None):
    """Weighted Cliff's delta sum_ij wx_i wy_j sign(x_i - y_j) / (sum wx sum wy). O((n+m) log(n+m))."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    wx = np.ones(len(x)) if wx is None else np.asarray(wx, float)
    wy = np.ones(len(y)) if wy is None else np.asarray(wy, float)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    o = np.argsort(y); ys, wys = y[o], wy[o]
    cw = np.concatenate([[0], np.cumsum(wys)])
    lo = np.searchsorted(ys, x, "left"); hi = np.searchsorted(ys, x, "right")
    below, above = cw[lo], cw[-1] - cw[hi]
    return float((wx * (below - above)).sum() / (wx.sum() * wys.sum()))


def within_cliff(d, col):
    """Pooled within-cluster Cliff's delta over all MBS-SBS pairs in the same cluster30 (rank / Mann-Whitney U form)."""
    r = d.groupby("cluster30")[col].rank().values
    g1 = (d.y_mbs.values == 1).astype(float)
    a = pd.DataFrame({"c": d.cluster30.values, "one": 1.0, "g1": g1, "r1": r * g1}).groupby("c").sum()
    n1, n0 = a.g1.values, a.one.values - a.g1.values
    u = a.r1.values - n1 * (n1 + 1) / 2
    num, den = (2 * u - n1 * n0).sum(), (n1 * n0).sum()
    return num / den if den else np.nan


def boot(d, stat, seed, n_boot=N_BOOT):
    rng = np.random.default_rng(seed)
    gi = d.groupby("cluster30").indices
    keys = list(gi); idx = [gi[k] for k in keys]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(keys), len(keys))
        dd = d.iloc[np.concatenate([idx[i] for i in pick])].copy()
        # resampled clusters are distinct units for within-cluster statistics
        dd["cluster30"] = np.repeat(np.arange(len(pick)), [len(idx[i]) for i in pick])
        vals.append(stat(dd))
    return np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)


def contact_chem(contacts, biolip=False):
    """Lining fractions of the census contact residues: all protein entries (4 Å set, n_simple >= 1) or, with
    biolip=True, only residues with >= BIOLIP_MIN_PAIRS vdW+0.5 Å atom pairs (BioLiP2 rule, DESIGN §0 primary)."""
    if contacts is None:
        return np.nan, np.nan, np.nan
    codes = [THREE2ONE.get(c["comp"], "X") for c in contacts
             if c["kind"] == "protein" and (not biolip or c["n_vdw"] >= BIOLIP_MIN_PAIRS)]
    n = len(codes)
    if not n:
        return np.nan, np.nan, np.nan
    return (np.mean([c in HYDROPHOBIC for c in codes]), np.mean([c in AROMATIC for c in codes]),
            np.mean([c in CHARGED for c in codes]))


def load(sample, describe_dir):
    s = pd.read_parquet(sample)
    parts = sorted(glob.glob(str(Path(describe_dir) / "part_*.parquet")))
    r = pd.concat([pd.read_parquet(f) for f in parts], ignore_index=True).drop_duplicates("site_id")
    d = s.merge(r, on="site_id", how="inner")
    for tag, bl in (("contact4A", False), ("contactBioLiP", True)):
        ch = np.array([contact_chem(c, bl) for c in d.contacts], float)
        d[f"{tag}_frac_hydrophobic"], d[f"{tag}_frac_aromatic"], d[f"{tag}_frac_charged"] = ch.T
    d["n_heavy_site"] = d.n_heavy_atoms_instance.where(d.n_heavy_atoms_instance > 0, d.n_heavy)
    d["log_heavy"] = np.log(d.n_heavy_site.clip(lower=1))
    pop = s.groupby("stratum").size().rename("n_sampled").to_frame()
    return s, r, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=str(INTERIM / "sample_n12000_seed0.parquet"))
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--describe-dir", default=str(INTERIM / "describe_uncropped"))
    ap.add_argument("--variant", default="uncropped", help="'uncropped' = primary (results in OUT); else OUT/<variant>")
    a = ap.parse_args()
    s, r, d = load(a.sample, a.describe_dir)
    primary = a.variant == "uncropped"
    RES = OUT if primary else OUT / a.variant
    RES.mkdir(parents=True, exist_ok=True)
    strata = pd.read_csv(sorted(glob.glob(str(OUT / f"sample_strata_n{len(s)}_seed0.tsv")))[0], sep="\t", comment="#")
    w = dict(zip(strata.stratum, strata.units_available / strata.sampled))
    # ---- run ledger
    led = [("sampled sites", len(s), s.cluster30.nunique()),
           ("described (any status)", len(d), d.cluster30.nunique()),
           ("status ok", int((d.status == "ok").sum()), d[d.status == "ok"].cluster30.nunique())]
    ok = d[d.status == "ok"].copy()
    led.append(("ligand rebuilt within LIG_MATCH_A of census centroid", int(ok.lig_match_ok.sum()),
                ok[ok.lig_match_ok].cluster30.nunique()))
    ok = ok[ok.lig_match_ok.astype(bool)].copy()
    led.append(("D-mode match (coverage>=MIN or precision>=0.5)", int(ok.d_match.sum()), ok[ok.d_match].cluster30.nunique()))
    err = d[d.status != "ok"].error.str.split(":").str[0].value_counts()
    for k, v in err.items():
        led.append((f"error: {k}", int(v), np.nan))
    write_tsv(pd.DataFrame(led, columns=["step", "n_sites", "n_clusters"]), RES / "run_ledger.tsv")
    ok["w"] = ok.stratum.map(w)
    ok.drop(columns=["fingerprint", "contacts"]).to_parquet(INTERIM / ("sites.parquet" if primary else f"sites_{a.variant}.parquet"), index=False)
    # ---- fingerprints for Part 2 positive controls
    fp = ok[["site_id", "pdb_id", "ccd_key", "ccd", "cluster30", "cls", "stratum_lig", "y_mbs", "assembly_type",
             "is_on_axis", "burrow_score", "volume_A3", "merged_pocket_interface", "residues", "fingerprint"]].copy()
    fp["fingerprint"] = [np.asarray(v, np.float32) for v in fp.fingerprint]
    fp["fp_unit"] = "located burrow member (pocketcompare describe --ligand --locate burrow), not the merged pocket"
    if primary:
        FP_OUT.parent.mkdir(parents=True, exist_ok=True); fp.to_parquet(FP_OUT, index=False)
        print("fingerprints", len(fp), "->", FP_OUT)

    # ---- power gate
    gate = ok.groupby(["stratum_lig", "y_mbs"]).cluster30.nunique().rename("n_clusters").reset_index()
    gate["passes"] = gate.n_clusters >= POWER_GATE_MIN_CLUSTERS
    write_tsv(gate, RES / "power_gate.tsv")

    # ---- 1. Cliff's delta MBS vs SBS
    rows = []
    for m, (col, role) in METRICS.items():
        for stratum in ["all_weighted", "all_unweighted", "all_weighted_dmatch", "drug_like", "cofactor", "other",
                        "homomer", "heteromer", "pocketgeom_weighted"]:
            x = ok.dropna(subset=[col]).copy()
            x["g"] = x.y_mbs
            wt = stratum in ("all_weighted", "all_weighted_dmatch", "pocketgeom_weighted")
            if stratum in ("drug_like", "cofactor", "other"):
                x = x[x.stratum_lig == stratum]
            elif stratum in ("homomer", "heteromer"):
                x = x[x.assembly_type == stratum]
            elif stratum == "all_weighted_dmatch":          # sensitivity: located BURROW pocket matches the ligand
                x = x[x.d_match.astype(bool)]
            elif stratum == "pocketgeom_weighted":          # Q2.2 definition at PDB sites: interface pocket (>=3
                x["g"] = np.where(x.merged_pocket_interface, 1,  # lining residues on >=2 chains) vs single-chain
                                  np.where(x.merged_n_lining_chains == 1, 0, -1))   # pocket; marginal excluded
                x = x[x.g >= 0]
            def st(dd, col=col, wt=wt):
                g1, g0 = dd[dd.g == 1], dd[dd.g == 0]
                return wcliff(g1[col], g0[col], g1.w if wt else None, g0.w if wt else None)
            for seed in (SEEDS if stratum == "all_weighted" else SEEDS[:1]):
                est = st(x); lo, hi = boot(x[["cluster30", "g", "w", col]], st, seed, a.n_boot)
                g1, g0 = x[x.g == 1][col], x[x.g == 0][col]
                rows.append(dict(metric=m, column=col, role=role, stratum=stratum, seed=seed, n_g1=len(g1), n_g0=len(g0),
                                 n_clusters=x.cluster30.nunique(), median_g1=g1.median(), median_g0=g0.median(),
                                 cliffs_delta=est, cd_lo=lo, cd_hi=hi))
        print("delta", m, flush=True)
    cd = pd.DataFrame(rows)
    write_tsv(cd, RES / "mbs_vs_sbs_cliffs_delta.tsv", note="g1=MBS vs g0=SBS (pocketgeom: g1=interface pocket vs g0=single-chain pocket); cluster30 bootstrap")

    # ---- 2. adjusted regressions + within-cluster
    sp = dmatrix(f"bs(log_heavy, df={SPLINE_DF}, include_intercept=False) - 1", ok, return_type="dataframe")
    sp.columns = [f"sp{i}" for i in range(sp.shape[1])]
    X = pd.concat([ok.reset_index(drop=True), sp.reset_index(drop=True)], axis=1)
    both = X.groupby("cluster30").y_mbs.transform(lambda v: v.nunique() == 2)
    reg = []
    spl = " + ".join(sp.columns)
    for m, (col, role) in METRICS.items():
        Z = X.dropna(subset=[col]).copy()
        Z["y_raw"] = Z[col]
        Z["y_pct"] = Z[col].rank(pct=True)
        for yv in ("y_raw", "y_pct"):
            sd = Z[yv].std()
            f = smf.ols(f"{yv} ~ y_mbs + {spl} + C(stratum_lig)", Z).fit(cov_type="cluster",
                                                                          cov_kwds={"groups": Z.cluster30})
            reg.append(dict(metric=m, role=role, outcome=yv, model="adjusted: spline(log heavy)+ligand stratum",
                            n=len(Z), n_clusters=Z.cluster30.nunique(), beta_mbs=f.params["y_mbs"],
                            lo=f.conf_int().loc["y_mbs", 0], hi=f.conf_int().loc["y_mbs", 1], p=f.pvalues["y_mbs"],
                            beta_over_sd=f.params["y_mbs"] / sd))
            W = Z[both.loc[Z.index]].copy()
            cols = [yv, "y_mbs"] + list(sp.columns) + ["dl", "cf"]
            W["dl"] = (W.stratum_lig == "drug_like").astype(float); W["cf"] = (W.stratum_lig == "cofactor").astype(float)
            Wd = W[cols] - W.groupby("cluster30")[cols].transform("mean")
            Wd["cluster30"] = W.cluster30
            f2 = smf.ols(f"{yv} ~ y_mbs + {spl} + dl + cf - 1", Wd).fit(cov_type="cluster",
                                                                        cov_kwds={"groups": Wd.cluster30})
            reg.append(dict(metric=m, role=role, outcome=yv, model="within-cluster FE (clusters with MBS and SBS)",
                            n=len(W), n_clusters=W.cluster30.nunique(), beta_mbs=f2.params["y_mbs"],
                            lo=f2.conf_int().loc["y_mbs", 0], hi=f2.conf_int().loc["y_mbs", 1], p=f2.pvalues["y_mbs"],
                            beta_over_sd=f2.params["y_mbs"] / sd))
    write_tsv(pd.DataFrame(reg), RES / "mbs_vs_sbs_regression.tsv",
              note="beta of MBS; y_pct = mid-rank percentile; cluster-robust SEs (cluster30)")

    # ---- 3. within-cluster Cliff's delta
    wrows = []
    Wc = ok[ok.groupby("cluster30").y_mbs.transform(lambda v: v.nunique() == 2)]
    for m, (col, role) in METRICS.items():
        x = Wc.dropna(subset=[col])
        est = within_cliff(x, col); lo, hi = boot(x[["cluster30", "y_mbs", col]], lambda dd, col=col: within_cliff(dd, col), 0, a.n_boot)
        wrows.append(dict(metric=m, role=role, n=len(x), n_clusters=x.cluster30.nunique(), cliffs_delta_within=est,
                          lo=lo, hi=hi, n_boot=a.n_boot))
    write_tsv(pd.DataFrame(wrows), RES / "mbs_vs_sbs_within_cluster.tsv")

    # ---- 4. on-axis vs off-axis MBS (homomers)
    H = ok[(ok.y_mbs == 1) & (ok.assembly_type == "homomer")].copy()
    H["on"] = H.is_on_axis.fillna(False).astype(bool)
    orow = []
    for m, (col, role) in METRICS.items():
        for stratum in ["all", "drug_like", "cofactor", "other"]:
            x = H.dropna(subset=[col]) if stratum == "all" else H[H.stratum_lig == stratum].dropna(subset=[col])
            if x.on.sum() == 0:
                continue
            st = lambda dd, col=col: wcliff(dd[dd.on][col], dd[~dd.on][col])
            lo, hi = boot(x[["cluster30", "on", col]], st, 0, a.n_boot)
            orow.append(dict(metric=m, role=role, stratum=stratum, n_on=int(x.on.sum()), n_off=int((~x.on).sum()),
                             clusters_on=x[x.on].cluster30.nunique(), clusters_off=x[~x.on].cluster30.nunique(),
                             median_on=x[x.on][col].median(), median_off=x[~x.on][col].median(),
                             cliffs_delta=st(x), cd_lo=lo, cd_hi=hi))
    write_tsv(pd.DataFrame(orow), RES / "onaxis_vs_offaxis_mbs.tsv", note="on-axis(g1) vs off-axis(g0) homomer MBS")

    # ---- 5. pocket-chain composition and comparison with Q2.2
    comp = (ok.groupby(["y_mbs", "stratum_lig"]).agg(n=("site_id", "size"),
                                                    merged_pocket_interface=("merged_pocket_interface", "mean"),
                                                    member_pocket_interface=("pocket_interface", "mean"),
                                                    lining_ge2_chains=("merged_n_lining_chains", lambda v: (v >= 2).mean()),
                                                    median_frac_second_chain=("merged_frac_second_chain", "median"))
            .reset_index())
    write_tsv(comp, RES / "pocket_chain_composition.tsv")
    q = pd.read_csv(Q22_PROPS, sep="\t", comment="#")
    q = q[(q.stratum == "all") & (q.seed == 0) & q.contrast.str.startswith("interface")]
    cmp_ = []
    for contrast, stv in (("PDB MBS vs SBS ligand sites", "all_weighted"),
                          ("PDB interface vs single-chain pocket at ligand sites", "pocketgeom_weighted")):
        ours = cd[(cd.stratum == stv) & (cd.seed == 0)].set_index("metric")
        for qm, om in Q22_MAP.items():
            qq = q[q.metric == qm]
            if not len(qq) or om not in ours.index:
                continue
            cmp_.append(dict(pdb_contrast=contrast, metric=qm, pdb_metric=om, af_dimer_delta=qq.cliffs_delta.iloc[0],
                             af_lo=qq.cd_lo.iloc[0], af_hi=qq.cd_hi.iloc[0], pdb_delta=ours.loc[om, "cliffs_delta"],
                             pdb_lo=ours.loc[om, "cd_lo"], pdb_hi=ours.loc[om, "cd_hi"],
                             same_sign=np.sign(qq.cliffs_delta.iloc[0]) == np.sign(ours.loc[om, "cliffs_delta"]),
                             pdb_ci_excludes_0=(ours.loc[om, "cd_lo"] > 0) or (ours.loc[om, "cd_hi"] < 0)))
    write_tsv(pd.DataFrame(cmp_), RES / "pdb_vs_afdimer_pattern.tsv",
              note="AF: interface vs single-chain pockets (Q2.2); PDB: MBS vs SBS ligand sites, weighted pooled")
    print(pd.DataFrame(cmp_).to_string())


if __name__ == "__main__":
    main()
