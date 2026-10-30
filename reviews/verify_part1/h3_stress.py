"""H3 stress tests (verifier; first attempt, not yet reviewed by the user).
(a) clusters per on/off x drug-like/cofactor cell; (b) exclude viral / HIV-protease; (c) leave-one-cluster-out
for top-20 clusters; (d) cluster-level votes; (e) top-15 on-axis drug-like clusters + all on-axis cofactor units.
Viral call is a heuristic regex on the mmCIF source-organism name of the dominant entity (no taxonomy DB).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

OUT = Path(__file__).resolve().parent / "out"
VIRAL_RE = r"virus|viral|phage|HIV|HTLV|SIV\b|retrovir|virinae|viridae"  # heuristic, see docstring
HIV_PR_RE = r"protease|retropepsin|pol polyprotein|gag-pol|gag pol"       # within viral entries
GATE = 30  # DESIGN §1 power gate: >=30 clusters per stratum

hm = pd.read_parquet(OUT / "h3_units.parquet")
meta = pd.read_parquet(OUT / "entity_meta.parquet")
meta = meta.dropna(subset=["entity"]).rename(columns={"entity": "dominant_entity"})
hm = hm.merge(meta, on=["pdb_id", "dominant_entity"], how="left")
hm["viral"] = hm.organism.fillna("").str.contains(VIRAL_RE, case=False, regex=True)
hm["hiv_like_protease"] = hm.viral & hm.description.fillna("").str.contains(HIV_PR_RE, case=False, regex=True)
print("organism missing:", hm.organism.isna().sum(), "viral units:", hm.viral.sum(), "viral protease units:", hm.hiv_like_protease.sum())
lines = []


def fit(df, label):
    a = ((df.on == 1) & (df.dl == 1)).sum(); b = ((df.on == 1) & (df.dl == 0)).sum()
    c = ((df.on == 0) & (df.dl == 1)).sum(); d = ((df.on == 0) & (df.dl == 0)).sum()
    orf, pf = stats.fisher_exact([[a, b], [c, d]])
    try:
        m = smf.glm("dl ~ on", df, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": df.my_cluster30.values})
        lo, hi, p = np.exp(m.params.on - 1.96 * m.bse.on), np.exp(m.params.on + 1.96 * m.bse.on), m.pvalues.on
    except Exception:
        lo = hi = p = np.nan
    cells = {f"cl_{o}_{k}": df[(df.on == oi) & (df.dl == ki)].my_cluster30.nunique()
             for o, oi in (("on", 1), ("off", 0)) for k, ki in (("dl", 1), ("cof", 0))}
    r = dict(variant=label, n=len(df), clusters=df.my_cluster30.nunique(), on_dl=a, on_cof=b, off_dl=c, off_cof=d,
             fisher_or=orf, robust_lo=lo, robust_hi=hi, robust_p=p, fisher_p=pf, **cells,
             gate_per_cell=min(cells.values()) >= GATE)
    lines.append(r)
    return r


fit(hm, "primary")
fit(hm[~hm.viral], "exclude viral (organism regex)")
fit(hm[~hm.hiv_like_protease], "exclude viral proteases")

# (c) leave-one-cluster-out for top-20 clusters by on-axis drug-like units (then by total units)
on_dl = hm[(hm.on == 1) & (hm.dl == 1)].groupby("my_cluster30").size().sort_values(ascending=False)
loo = []
for c in on_dl.index[:20]:
    r = fit(hm[hm.my_cluster30 != c], f"LOO cluster {c}")
    loo.append(dict(cluster30=c, n_on_dl=on_dl[c], or_without=r["fisher_or"], lo=r["robust_lo"], hi=r["robust_hi"]))
lines[:] = [l for l in lines if not l["variant"].startswith("LOO")]
loo = pd.DataFrame(loo)
for k in (2, 5, 10, 20):
    fit(hm[~hm.my_cluster30.isin(on_dl.index[:k])], f"drop top-{k} on-axis-drug-like clusters jointly")

# (d) cluster-level votes
# d1: each cluster x axis-status cell votes 'drug-like' if majority of its units are drug-like
cv = hm.groupby(["my_cluster30", "on"]).dl.mean().reset_index()
cv["vote_dl"] = (cv.dl > 0.5).astype(int)
t = pd.crosstab(cv.on, cv.vote_dl)
orv, pv = stats.fisher_exact([[t.loc[1, 1], t.loc[1, 0]], [t.loc[0, 1], t.loc[0, 0]]])
cluster_vote = dict(on_dl=t.loc[1, 1], on_cof=t.loc[1, 0], off_dl=t.loc[0, 1], off_cof=t.loc[0, 0], fisher_or=orv, p=pv)
# d2: one unit per cluster x on-status x class (presence/absence), unit = cluster
pres = hm.drop_duplicates(["my_cluster30", "on", "dl"])
fit(pres, "presence-only (one unit per cluster x axis x class)")
# d3: within-cluster paired: clusters having both on- and off-axis units
both = hm.groupby("my_cluster30").on.transform("nunique").eq(2)
pc = hm[both].groupby(["my_cluster30", "on"]).dl.mean().unstack()
paired = dict(n_clusters=len(pc), on_more_dl=int((pc[1] > pc[0]).sum()), off_more_dl=int((pc[1] < pc[0]).sum()),
              tie=int((pc[1] == pc[0]).sum()))
nb = paired["on_more_dl"] + paired["off_more_dl"]
paired["sign_test_p"] = stats.binomtest(paired["on_more_dl"], nb).pvalue if nb else np.nan

# (e) top-15 on-axis drug-like clusters, example entries/ligands; all on-axis cofactor units
od = hm[(hm.on == 1) & (hm.dl == 1)]
top = (od.groupby("my_cluster30").agg(n_units=("ccd_key", "size"), protein=("description", lambda s: s.value_counts().index[0] if s.notna().any() else ""),
                                     organism=("organism", lambda s: s.value_counts().index[0] if s.notna().any() else ""),
                                     example_pdb=("pdb_id", lambda s: ",".join(s.iloc[:4])),
                                     example_ligands=("name", lambda s: "; ".join(s.dropna().astype(str).str[:30].iloc[:3])))
       .sort_values("n_units", ascending=False))
top["cum_frac"] = top.n_units.cumsum() / len(od)
oc = hm[(hm.on == 1) & (hm.dl == 0)][["my_cluster30", "pdb_id", "ccd_key", "name", "description", "organism"]]
# family breadth: distinct clusters / uniprot / organisms contributing on-axis drug-like
breadth = dict(on_dl_units=len(od), clusters=od.my_cluster30.nunique(), uniprots=od.dominant_uniprot.nunique(),
               viral_units=int(od.viral.sum()), viral_protease_units=int(od.hiv_like_protease.sum()),
               clusters_with_ge5=int((top.n_units >= 5).sum()), share_top2=top.n_units.iloc[:2].sum() / len(od),
               share_top10=top.n_units.iloc[:10].sum() / len(od))

pd.DataFrame(lines).to_csv(OUT / "h3_stress.tsv", sep="\t", index=False)
loo.to_csv(OUT / "h3_loo_top20.tsv", sep="\t", index=False)
top.to_csv(OUT / "h3_top_on_axis_druglike_clusters.tsv", sep="\t")
oc.to_csv(OUT / "h3_on_axis_cofactor_units.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(pd.DataFrame(lines).to_string())
print(loo.to_string())
print("cluster vote:", cluster_vote)
print("paired within-cluster:", paired)
print("breadth:", breadth)
print(top.head(15).to_string())
print(oc.to_string())
