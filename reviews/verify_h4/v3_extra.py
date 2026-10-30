"""Verifier extras: cluster sign-flip for raw MBS gain, family-level sign test, within-cluster permutation at region level,
design effect. Reads reviews/verify_h4/primary_with_regions.parquet (written by v1_stats.py). UNVERIFIED first attempt."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binomtest
OUT = Path(__file__).resolve().parent
N_PERM, SEED = 10000, 0
P = pd.read_parquet(OUT / "primary_with_regions.parquet")
rng = np.random.default_rng(SEED)

def signflip(df, n=N_PERM):
    cs = df.groupby("boot_cluster").gain.sum().values; obs = cs.sum() / len(df)
    fl = rng.choice([-1.0, 1.0], size=(n, len(cs)))
    null = fl @ cs / len(df)
    return obs, (np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (n + 1)
for lab, c in [("MBS", True), ("SBS", False)]:
    print(lab, "gain, cluster sign-flip p:", signflip(P[P.is_mbs == c]))

# family-level: per cluster MBS gain and SBS gain
cg = P.groupby(["boot_cluster", "is_mbs"]).gain.mean().unstack()
for lab, c in [("MBS", True), ("SBS", False)]:
    x = cg[c].dropna(); print(lab, "clusters", len(x), ">0:", int((x > 0).sum()), "<0:", int((x < 0).sum()),
                         "sign p:", binomtest(int((x > 0).sum()), int((x != 0).sum())).pvalue, "mean:", round(x.mean(), 4), "sd:", round(x.std(), 3))
mix = cg.dropna(); d = mix[True] - mix[False]
print("mixed clusters", len(mix), "DiD>0:", int((d > 0).sum()), "<0:", int((d < 0).sum()), "sign p:", binomtest(int((d > 0).sum()), int((d != 0).sum())).pvalue)

# design effect: ICC of gain within cluster -> effective n
def deff(df):
    g = df.groupby("boot_cluster").gain; m = g.size(); 
    msb = (m * (g.mean() - df.gain.mean()) ** 2).sum() / (len(m) - 1)
    msw = ((df.gain - g.transform("mean")) ** 2).sum() / (len(df) - len(m))
    m0 = (len(df) - (m ** 2).sum() / len(df)) / (len(m) - 1)
    icc = (msb - msw) / (msb + (m0 - 1) * msw)
    return icc, 1 + ((m ** 2).sum() / m.sum() - 1) * icc
for lab, c in [("MBS", True), ("SBS", False)]:
    icc, de = deff(P[P.is_mbs == c]); print(lab, "ICC of gain within cluster", round(icc, 3), "design effect", round(de, 1), "n_eff", round((P.is_mbs == c).sum() / de))

# within-cluster permutation on region-averaged units (one unit per cluster x region x class)
R = P.groupby("region_cls").agg(gain=("gain", "mean"), is_mbs=("is_mbs", "first"), boot_cluster=("boot_cluster", "first")).reset_index()
def wperm(df):
    df = df.sort_values("boot_cluster").reset_index(drop=True)
    inv = pd.factorize(df.boot_cluster)[0]; m = df.is_mbs.values; g = df.gain.values
    obs = g[m].mean() - g[~m].mean(); null = np.empty(N_PERM)
    for i in range(N_PERM):
        mp = m[np.argsort(inv + rng.random(len(inv)), kind="stable")]
        null[i] = g[mp].mean() - g[~mp].mean()
    mu = null.mean()
    return obs, mu, null.std(), (np.sum(np.abs(null - mu) >= abs(obs - mu)) + 1) / (N_PERM + 1)
print("region-level within-cluster perm (obs, null mean, null sd, centred p):", wperm(R))
# within-cluster perm restricted to mixed clusters only, null should be ~0
mixed = P.groupby("boot_cluster").is_mbs.transform(lambda x: x.any() and (~x).any())
print("site-level, mixed clusters only:", wperm(P[mixed]))
print("region-level, mixed clusters only:", wperm(R[R.groupby('boot_cluster').is_mbs.transform(lambda x: x.any() and (~x).any())]))
