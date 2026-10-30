"""Independent verifier: recompute H4 headline table, bootstrap, permutation and heterogeneity analyses
from the producer's site-level table (data/part2/h4/v2/benchmark_filtered.parquet). Does not import scripts/h4_*.py.
Status: verifier's first attempt, UNVERIFIED by anyone else."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[2]
SITE_TABLE = WS / "data/part2/h4/v2/benchmark_filtered.parquet"
OUT = WS / "reviews/verify_h4"
N_BOOT, SEED, N_PERM = 1000, 0, 10000
REGION_OC = 0.5   # overlap coefficient for merging binding-residue sets into one binding region (verifier choice)

d = pd.read_parquet(SITE_TABLE)
P = d[d.lig_set == "organic"].copy()          # primary set = organic, interface-filtered (table is already filtered)
assert P.iface_jaccard.min() >= 0.3

# --- own clusters: connected components over model and cluster30 (on primary set only; producer used the mapped set)
par = {}
def f(x):
    par.setdefault(x, x)
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for m, c in zip(P.model, P.cluster30):
    a, b = f("m:" + m), f("c:" + str(c)); par[a] = b
P["vclust"] = [f("m:" + m) for m in P.model]
print("own clusters on primary:", P.vclust.nunique(), "producer boot_cluster:", P.boot_cluster.nunique())
# do our components refine/equal producer's?
print("producer clusters split by ours:", (P.groupby("boot_cluster").vclust.nunique() > 1).sum(),
      "ours spanning >1 producer:", (P.groupby("vclust").boot_cluster.nunique() > 1).sum())

P["gain"] = P.dimer_hit.astype(float) - P.iso_hit.astype(float)
M = P.is_mbs.values

def point(df, w=None):
    w = np.ones(len(df)) if w is None else w
    m = df.is_mbs.values
    g = df.gain.values
    gm = np.sum(w[m] * g[m]) / w[m].sum(); gs = np.sum(w[~m] * g[~m]) / w[~m].sum()
    return gm - gs, gm, gs

rows = []
def boot(df, label, clcol="boot_cluster", stratify=False, weight=None, seed=SEED):
    """Cluster bootstrap: resample clusters with replacement; weight = None (site) or 'family' (equal weight per cluster per class)."""
    rng = np.random.default_rng(seed)
    df = df.reset_index(drop=True)
    cl, inv = np.unique(df[clcol].values, return_inverse=True)
    m = df.is_mbs.values; g = df.gain.values
    def stat(mult):  # mult per cluster
        w = mult[inv].astype(float)
        if weight == "family":
            # each cluster's MBS sites (and SBS sites) share total weight 1 per cluster draw
            nm = np.bincount(inv, weights=m.astype(float), minlength=len(cl))
            ns = np.bincount(inv, weights=(~m).astype(float), minlength=len(cl))
            w = np.where(m, w / np.maximum(nm[inv], 1), w / np.maximum(ns[inv], 1))
        a = np.sum(w[m] * g[m]); b = w[m].sum(); c = np.sum(w[~m] * g[~m]); e = w[~m].sum()
        gm = a / b if b > 0 else np.nan; gs = c / e if e > 0 else np.nan
        return gm - gs, gm, gs
    obs = stat(np.ones(len(cl)))
    kinds = df.groupby(clcol).kind.first().reindex(cl).values
    B = np.empty((N_BOOT, 3))
    for i in range(N_BOOT):
        if stratify:
            mult = np.zeros(len(cl))
            for k in np.unique(kinds):
                idx = np.where(kinds == k)[0]
                np.add.at(mult, rng.choice(idx, len(idx)), 1)
        else:
            mult = np.bincount(rng.integers(0, len(cl), len(cl)), minlength=len(cl))
        B[i] = stat(mult)
    lo, hi = np.nanpercentile(B[:, 0], [2.5, 97.5])
    p = min(1.0, 2 * min(np.mean(B[:, 0] <= 0), np.mean(B[:, 0] >= 0)))
    r = dict(analysis=label, n=len(df), n_mbs=int(m.sum()), n_clusters=len(cl), did=obs[0], gain_mbs=obs[1], gain_sbs=obs[2],
             did_lo=lo, did_hi=hi, p_boot_2sided=p,
             gain_mbs_lo=np.nanpercentile(B[:, 1], 2.5), gain_mbs_hi=np.nanpercentile(B[:, 1], 97.5),
             gain_sbs_lo=np.nanpercentile(B[:, 2], 2.5), gain_sbs_hi=np.nanpercentile(B[:, 2], 97.5))
    rows.append(r); print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    return r

# 1. headline
for kind_set, sub in [("homo+hetero", P), ("homo", P[P.kind == "homo"]), ("hetero", P[P.kind == "hetero"])]:
    for lab, col in [("MBS", True), ("SBS", False)]:
        s = sub[sub.is_mbs == col]
        print(kind_set, lab, len(s), "dimer", round(s.dimer_hit.mean(), 4), "iso", round(s.iso_hit.mean(), 4))
boot(P, "primary, unstratified cluster bootstrap (producer clusters)")
boot(P, "primary, kind-stratified cluster bootstrap (producer clusters)", stratify=True)
boot(P, "primary, unstratified, verifier clusters (primary-set components)", clcol="vclust")
for sd in (1, 2):
    boot(P, f"primary, unstratified, seed {sd}", seed=sd)

# discordant pairs
for lab, col in [("MBS", True), ("SBS", False)]:
    s = P[P.is_mbs == col]
    print(lab, "dimer-only", int((s.dimer_hit & ~s.iso_hit).sum()), "iso-only", int((~s.dimer_hit & s.iso_hit).sum()))

# 4. permutation: within-cluster MBS-label permutation, raw (uncentred) and centred p
rng = np.random.default_rng(SEED)
Pr = P.reset_index(drop=True)
cl, inv = np.unique(Pr.boot_cluster.values, return_inverse=True)
order = np.argsort(inv, kind="stable"); inv_s = inv[order]; m_s = Pr.is_mbs.values[order]; g_s = Pr.gain.values[order]
starts = np.r_[0, np.flatnonzero(np.diff(inv_s)) + 1, len(inv_s)]
obs = g_s[m_s].mean() - g_s[~m_s].mean()
null = np.empty(N_PERM)
for i in range(N_PERM):
    key = inv_s + rng.random(len(inv_s))
    mp = m_s[np.argsort(key, kind="stable")]
    null[i] = g_s[mp].mean() - g_s[~mp].mean()
mu = null.mean()
p_unc = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
p_cen = (np.sum(np.abs(null - mu) >= abs(obs - mu)) + 1) / (N_PERM + 1)
print(f"perm: obs {obs:.4f} null mean {mu:.4f} sd {null.std():.4f}; p uncentred {p_unc:.4g}, p centred {p_cen:.4g}")
# Within-cluster (stratified) DiD estimator: only mixed clusters, weighted by harmonic n
mix = []
for c, g in Pr.groupby("boot_cluster"):
    if g.is_mbs.any() and (~g.is_mbs).any():
        nm, ns = g.is_mbs.sum(), (~g.is_mbs).sum()
        mix.append((c, g[g.is_mbs].gain.mean() - g[~g.is_mbs].gain.mean(), nm * ns / (nm + ns), nm, ns))
mix = pd.DataFrame(mix, columns=["cluster", "did_c", "w", "nm", "ns"])
print("mixed clusters:", len(mix), "sites in mixed clusters:", int(mix.nm.sum() + mix.ns.sum()),
      "MH-type within-cluster DiD:", round(np.average(mix.did_c, weights=mix.w), 4),
      "unweighted mean of cluster DiDs:", round(mix.did_c.mean(), 4), "frac clusters DiD>0:", round((mix.did_c > 0).mean(), 3),
      "<0:", round((mix.did_c < 0).mean(), 3))
pure = Pr.groupby("boot_cluster").is_mbs.agg(["all", lambda x: (~x).all()])
print("MBS-only clusters", int(pure["all"].sum()), "SBS-only clusters", int(pure.iloc[:, 1].sum()))
mix.sort_values("w", ascending=False).to_csv(OUT / "within_cluster_did.tsv", sep="\t", index=False)

# family dominance: per cluster contribution
fam = Pr.groupby(["boot_cluster", "is_mbs"]).agg(n=("gain", "size"), gain=("gain", "mean"),
      accs=("model_acc_A", lambda x: ",".join(sorted(set(x))[:4]))).reset_index()
fam["contrib"] = fam.n * fam.gain / np.where(fam.is_mbs, M.sum(), (~M).sum()) * np.where(fam.is_mbs, 1, -1)
fam = fam.sort_values("contrib", key=np.abs, ascending=False)
fam.to_csv(OUT / "family_contributions.tsv", sep="\t", index=False)
print(fam.head(15).to_string())

# leave-one-cluster-out influence on DiD
loo = []
for c in Pr.boot_cluster.unique():
    s = Pr[Pr.boot_cluster != c]
    loo.append((c, point(s)[0]))
loo = pd.DataFrame(loo, columns=["cluster", "did_without"]).sort_values("did_without")
print("LOO DiD range:", loo.did_without.min().round(4), loo.did_without.max().round(4))
print(loo.head(5).to_string(), loo.tail(5).to_string())

# NAMPT
nam = Pr[(Pr.model_acc_A == "P43490") | (Pr.model_acc_B == "P43490")]
print("NAMPT sites:", len(nam), nam.groupby("is_mbs")[["dimer_hit", "iso_hit"]].agg(["mean", "size"]).to_string())
boot(Pr[Pr.boot_cluster != nam.boot_cluster.iloc[0]], "primary without NAMPT cluster")

# region dedup: binding region = union-find over sites in same cluster with overlap coeff >= REGION_OC on (acc, resnum) sets
def resset(r):
    return {(r.acc1, int(x)) for x in r.bind1} | {(r.acc2, int(x)) for x in r.bind2}
Pr["rs"] = [resset(r) for r in Pr.itertuples()]
reg = np.empty(len(Pr), dtype=object)
for c, idx in Pr.groupby("boot_cluster").indices.items():
    idx = list(idx); parent = list(range(len(idx)))
    def ff(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    sets = [Pr.rs.iloc[i] for i in idx]
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            a, b = sets[i], sets[j]
            if a and b and len(a & b) / min(len(a), len(b)) >= REGION_OC:
                parent[ff(i)] = ff(j)
    for k, i in enumerate(idx):
        reg[i] = f"{c}|{ff(k)}"
Pr["region"] = reg
# a region may contain both MBS and SBS sites; key on (cluster, region, class)
Pr["region_cls"] = Pr.region + "|" + Pr.is_mbs.astype(str)
print("regions (cluster x region x class):", Pr.region_cls.nunique(), "regions:", Pr.region.nunique())
rng2 = np.random.default_rng(SEED)
one = Pr.loc[Pr.groupby("region_cls").apply(lambda g: rng2.choice(g.index.values), include_groups=False).values]
boot(one, "one random site per (cluster, binding region, class), seed 0")
# region-averaged: each region weight 1
avg = Pr.groupby("region_cls").agg(gain=("gain", "mean"), is_mbs=("is_mbs", "first"), boot_cluster=("boot_cluster", "first"),
                                    kind=("kind", "first")).reset_index()
boot(avg, "region-averaged gain (each region weight 1)")
boot(Pr, "family-weighted (each cluster weight 1 per class)", weight="family")
P_nn = Pr[Pr.boot_cluster != nam.boot_cluster.iloc[0]]
boot(P_nn, "family-weighted without NAMPT cluster", weight="family")

# 5. collapse-rule sensitivities
for lab, a, b in [("no collapse", "dimer_nocollapse_hit", "iso_nocollapse_hit"), ("OC collapse", "dimer_occollapse_hit", "iso_occollapse_hit")]:
    s = Pr.copy(); s["gain"] = s[a].astype(float) - s[b].astype(float)
    for cls in (True, False):
        q = s[s.is_mbs == cls]; print(lab, "MBS" if cls else "SBS", round(q[a].mean(), 4), round(q[b].mean(), 4))
    boot(s, f"collapse sensitivity: {lab}")
    boot(s, f"collapse sensitivity: {lab}, family-weighted", weight="family")
# how many sites change hit status between rules
for a in ["dimer", "iso"]:
    for alt in ["nocollapse", "occollapse"]:
        ch = (Pr[f"{a}_hit"] != Pr[f"{a}_{alt}_hit"])
        print(a, alt, "changed:", int(ch.sum()), "MBS changed:", int((ch & Pr.is_mbs).sum()),
              "gained:", int((~Pr[f"{a}_hit"] & Pr[f"{a}_{alt}_hit"]).sum()), "lost:", int((Pr[f"{a}_hit"] & ~Pr[f"{a}_{alt}_hit"]).sum()))
pd.DataFrame(rows).to_csv(OUT / "stats_results.tsv", sep="\t", index=False)
Pr.drop(columns=["rs"]).to_parquet(OUT / "primary_with_regions.parquet")
