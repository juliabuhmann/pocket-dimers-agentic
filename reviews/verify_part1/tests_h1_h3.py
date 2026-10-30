"""Verifier's own estimates of H1-H3 + H1 strict-MBS stress tests (does NOT import scripts/part1_analysis.py).

Input: reviews/verify_part1/out/units_verifier.parquet (from build_units.py).
Status: verifier's first attempt; estimates are compared against producer numbers in reviews/verify_part1_results.md.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.discrete.conditional_models import ConditionalLogit
from statsmodels.stats.contingency_tables import StratifiedTable

OUT = Path(__file__).resolve().parent / "out"
SEED = 0           # DESIGN §5 primary seed
N_BOOT = 1000      # DESIGN §5 bootstrap resamples
NONORGANIC = ["metal_ion", "metal_cluster", "inorganic"]  # DESIGN §0: metals never pooled with organic ligands

u = pd.read_parquet(OUT / "units_verifier.parquet")
org = u[~u.fine.isin(NONORGANIC)].copy()
org["cls"] = pd.Categorical(org.fine, categories=["drug_like"] + sorted(set(org.fine) - {"drug_like"}))
org["sg"] = org.is_structural_genomics.fillna(False).astype(int)
org["res"] = org.resolution.fillna(org.resolution.median())
org["yr"] = org.year.fillna(org.year.median())
rows = []


def h1(df, outcome, label):
    df = df.copy()
    df["y"] = df[outcome].astype(int)
    f = "y ~ C(cls) + bs(heavy, df=4) + res + C(method_class) + n_protein_chains + sg + yr"
    m = smf.glm(f, df, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": df.my_cluster30.values})
    t = [k for k in m.params.index if "cofactor" in k and "C(cls)" in k][0]
    b, se = m.params[t], m.bse[t]
    r = dict(test="H1", variant=label, n=len(df), n_clusters=df.my_cluster30.nunique(),
             est=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[t],
             frac_cof=df.loc[df.fine == "cofactor", "y"].mean(), frac_dl=df.loc[df.fine == "drug_like", "y"].mean())
    rows.append(r)
    print(r)
    return m


h1(org, "mbs", "primary")
h1(org, "mbs_any_entry", "MBS in any entry")
h1(org, "mbs_strict", "strict MBS")

# ---- H2: conditional logistic, drug-like (vs other organic) ~ MBS, strata = cluster30
def h2(df, outcome_mbs, label):
    df = df.copy()
    df["dl"] = (df.fine == "drug_like").astype(int)
    df["x"] = df[outcome_mbs].astype(int)
    g = df.groupby("my_cluster30")
    keep = g.x.transform("nunique").eq(2)
    df = df[keep]
    disc = df.groupby("my_cluster30").dl.transform("nunique").eq(2)
    n_disc = df.loc[disc, "my_cluster30"].nunique()
    dd = df[disc]
    # exact conditional likelihood for one binary predictor = Fisher noncentral hypergeometric per stratum
    agg = dd.groupby("my_cluster30").agg(N=("x", "size"), K=("x", "sum"), n=("dl", "sum"),
                                         a=("dl", lambda v: int((v.values & dd.loc[v.index, "x"].values).sum())))
    from scipy.optimize import minimize_scalar
    def nll(lb):
        return -stats.nchypergeom_fisher.logpmf(agg.a.values, agg.N.values, agg.K.values, agg.n.values, np.exp(lb)).sum()
    b = minimize_scalar(nll, bounds=(-5, 5), method="bounded").x
    h = 1e-3
    se = 1 / np.sqrt((nll(b + h) - 2 * nll(b) + nll(b - h)) / h ** 2)
    lr = 2 * (nll(0.0) - nll(b))
    class _M: pass
    m = _M(); m.pvalues = [stats.chi2.sf(lr, 1)]
    tabs = [np.array([[((s.x == 1) & (s.dl == 1)).sum(), ((s.x == 1) & (s.dl == 0)).sum()],
                      [((s.x == 0) & (s.dl == 1)).sum(), ((s.x == 0) & (s.dl == 0)).sum()]]) for _, s in dd.groupby("my_cluster30")]
    st = StratifiedTable(tabs)
    r = dict(test="H2", variant=label, n=len(df), n_clusters=df.my_cluster30.nunique(), n_discordant=n_disc,
             est=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[0], mh_or=st.oddsratio_pooled,
             frac_cof=df.loc[df.x == 1, "dl"].mean(), frac_dl=df.loc[df.x == 0, "dl"].mean())
    rows.append(r)
    print(r)


h2(org, "mbs", "primary (frac_cof=dl share at MBS, frac_dl=dl share at SBS)")
h2(org, "mbs_strict", "strict MBS")

# ---- H3
hm = org[(org.assembly_type == "homomer") & org.mbs & org.fine.isin(["drug_like", "cofactor"])].copy()
hm["dl"] = (hm.fine == "drug_like").astype(int)
hm["on"] = hm.on_axis.astype(int)


def h3(df, label, on="on"):
    a = ((df[on] == 1) & (df.dl == 1)).sum(); b_ = ((df[on] == 1) & (df.dl == 0)).sum()
    c = ((df[on] == 0) & (df.dl == 1)).sum(); d_ = ((df[on] == 0) & (df.dl == 0)).sum()
    orf, pf = stats.fisher_exact([[a, b_], [c, d_]])
    cells = {f"cl_{o}_{k}": df[(df[on] == oi) & (df.dl == ki)].my_cluster30.nunique()
             for o, oi in (("on", 1), ("off", 0)) for k, ki in (("dl", 1), ("cof", 0))}
    try:
        m = smf.glm(f"dl ~ {on}", df, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": df.my_cluster30.values})
        b, se = m.params[on], m.bse[on]
        lo, hi, p = np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), m.pvalues[on]
    except Exception as ex:  # noqa
        lo = hi = p = np.nan
    r = dict(test="H3", variant=label, n=len(df), n_clusters=df.my_cluster30.nunique(), a_on_dl=a, b_on_cof=b_,
             c_off_dl=c, d_off_cof=d_, est=orf, lo=lo, hi=hi, p=p, fisher_p=pf, **cells)
    rows.append(r)
    print(r)
    return r


h3(hm, "primary")
hm["on_geom"] = hm.on_axis_geom.astype(int)
h3(hm, "geometric on-axis", on="on_geom")
hs = org[(org.assembly_type == "homomer") & org.mbs_strict & org.fine.isin(["drug_like", "cofactor"])].copy()
hs["dl"] = (hs.fine == "drug_like").astype(int); hs["on"] = hs.on_axis.astype(int)
h3(hs, "strict MBS")
hm.to_parquet(OUT / "h3_units.parquet")

# ---- H1 strict-MBS stress: drop top cofactor ligands / top families
s = org.copy()
cof = s[s.fine == "cofactor"]
top_cof = cof.groupby("ccd_key").agg(n=("mbs_strict", "size"), n_strict=("mbs_strict", "sum"), name=("name", "first")).sort_values("n", ascending=False)
top_cof.head(25).to_csv(OUT / "h1_top_cofactors.tsv", sep="\t")
print(top_cof.head(15))
for k in (["NAD", "NAI", "NAP", "NDP"], ["FAD", "FMN"], ["HEM"], ["SAH", "SAM"], ["PLP", "PMP"], ["COA", "ACO"]):
    try:
        h1(s[~s.ccd_key.isin(k)], "mbs_strict", f"strict MBS, drop {'/'.join(k)}")
    except Exception as ex:
        print("fail", k, ex)
top10 = top_cof.index[:10].tolist()
h1(s[~s.ccd_key.isin(top10)], "mbs_strict", f"strict MBS, drop top-10 cofactor CCDs ({','.join(top10)})")
h1(s[~s.ccd_key.isin(top10)], "mbs", f"primary MBS, drop top-10 cofactor CCDs")
# per assembly type
for at in ("homomer", "heteromer"):
    h1(s[s.assembly_type == at], "mbs_strict", f"strict MBS, {at} only")
    h1(s[s.assembly_type == at], "mbs", f"primary MBS, {at} only")
# leave-one-cofactor-out influence on strict OR (top 20 by n)
infl = []
for k in top_cof.index[:20]:
    m = h1(s[s.ccd_key != k], "mbs_strict", f"strict LOO cofactor {k}")
    infl.append((k, rows[-1]["est"]))
rows[:] = [r for r in rows if not str(r["variant"]).startswith("strict LOO cofactor")]
pd.DataFrame(infl, columns=["dropped_ccd", "strict_or"]).to_csv(OUT / "h1_strict_loo_cofactor.tsv", sep="\t", index=False)
# leave-one-cluster-out: clusters contributing most strict-MBS cofactor or drug-like units
contrib = s[s.fine.isin(["cofactor", "drug_like"])].groupby("my_cluster30").size().sort_values(ascending=False)
infc = []
for c in contrib.index[:20]:
    h1(s[s.my_cluster30 != c], "mbs_strict", "tmp")
    infc.append((c, contrib[c], rows[-1]["est"]))
rows[:] = [r for r in rows if r["variant"] != "tmp"]
pd.DataFrame(infc, columns=["dropped_cluster30", "n_units_cof_or_dl", "strict_or"]).to_csv(OUT / "h1_strict_loo_cluster.tsv", sep="\t", index=False)
# strict-MBS cofactor vs drug-like fraction by heavy-atom size bin (is the reversal a size effect?)
sd = s[s.fine.isin(["cofactor", "drug_like"])].copy()
sd["hbin"] = pd.cut(sd.heavy, [0, 20, 30, 40, 60, 1000])
print(sd.groupby(["hbin", "fine"], observed=True).agg(n=("mbs", "size"), mbs=("mbs", "mean"), strict=("mbs_strict", "mean")).to_string())
sd.groupby(["hbin", "fine"], observed=True).agg(n=("mbs", "size"), mbs=("mbs", "mean"), strict=("mbs_strict", "mean")).to_csv(OUT / "h1_size_bins.tsv", sep="\t")

pd.DataFrame(rows).to_csv(OUT / "verifier_tests.tsv", sep="\t", index=False)
