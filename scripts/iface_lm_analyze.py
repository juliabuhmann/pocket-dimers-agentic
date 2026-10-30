#!/usr/bin/env python
"""EXPLORATORY: lining chemistry of ligand-matched vs unmatched predicted interface pockets (H4 benchmark models).
Groups (atlas-pass, symmetric reps, no 1-2-residue 'marginal' pockets; q22_compare.py definitions):
  IM = interface pocket matched to a primary-set organic MBS site     IU = interface pocket matched to NO site (any row)
  SM = single-chain pocket matched to a primary-set organic SBS site  SU = single-chain pocket matched to NO site (reference)
Stats: cluster-bootstrap (paralog unit, seed 0, N_BOOT) CIs for medians and Cliff's delta; classifier check with cluster-split CV.
Run scripts/iface_lm_match.py first. Status: first attempt, not verified.
"""
import sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from q22_common import INTERIM, WS, N_BOOT, cliffs_delta, cluster_boot, write_tsv

OUT = WS / "results/part2/iface_ligand_matched"
MATCH = INTERIM.parent / "iface_ligand_matched/pocket_match.parquet"
METRICS = ["frac_charged", "frac_hydrophobic", "frac_aromatic", "volume_A3", "burrow_score", "buriedness"]
CONTRASTS = [("IM", "IU"), ("IM", "SM"), ("IU", "SM"), ("IU", "SU")]
SEED, K_FOLDS = 0, 5


def build(match_cols):
    d = pd.read_parquet(INTERIM / "dimer_pocket_features_units.parquet")
    m = pd.read_parquet(MATCH).set_index("pocket_id")
    d = d[d.pocket_id.isin(m.index)].join(m, on="pocket_id")
    S = d[d.atlas_pass & d.is_sym_rep]
    S = S[~(S.is_interface_ge1 & ~S.is_interface)].copy()
    mb, sb = match_cols
    S["grp"] = np.select([S.is_interface & S[mb], S.is_interface & ~S.match_any, ~S.is_interface_ge1 & S[sb],
                          ~S.is_interface_ge1 & ~S.match_any], ["IM", "IU", "SM", "SU"], "other")
    return S


def job_med(a):
    g, met, df = a
    return (g, met) + cluster_boot(df, "unit", lambda x: np.median(x[met]), N_BOOT, SEED)


def job_delta(a):
    g1, g0, met, df = a
    est, lo, hi = cluster_boot(df, "unit", lambda x: cliffs_delta(x.loc[x.grp == g1, met], x.loc[x.grp == g0, met]), N_BOOT, SEED)
    return dict(g1=g1, g0=g0, metric=met, n1=int((df.grp == g1).sum()), n0=int((df.grp == g0).sum()),
                clusters1=df.unit[df.grp == g1].nunique(), clusters0=df.unit[df.grp == g0].nunique(), delta=est, lo=lo, hi=hi)


def auc(y, s):
    from scipy.stats import rankdata
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def logit_fit(X, y, l2=1.0, it=50):
    X1 = np.c_[np.ones(len(X)), X]; w = np.zeros(X1.shape[1]); R = np.eye(len(w)) * l2; R[0, 0] = 0
    for _ in range(it):
        p = 1 / (1 + np.exp(-X1 @ w)); g = X1.T @ (p - y) + R @ w
        H = (X1 * (p * (1 - p))[:, None]).T @ X1 + R
        step = np.linalg.solve(H, g); w -= step
        if np.abs(step).max() < 1e-8: break
    return w


def cv_auc(D, feats):
    rng = np.random.default_rng(SEED)
    units = D.unit.unique(); fold_of = dict(zip(units, rng.permutation(len(units)) % K_FOLDS))
    f = D.unit.map(fold_of).to_numpy(); y = (D.grp == "IM").to_numpy().astype(float)
    X = D[feats].to_numpy(float); oof = np.zeros(len(D)); fa = []
    for k in range(K_FOLDS):
        tr, te = f != k, f == k
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-12
        w = logit_fit((X[tr] - mu) / sd, y[tr]); oof[te] = np.c_[np.ones(te.sum()), (X[te] - mu) / sd] @ w
        fa.append(auc(y[te], oof[te]))
    est, lo, hi = cluster_boot(D.assign(oof=oof, y=y), "unit", lambda x: auc(x.y.to_numpy(), x.oof.to_numpy()), N_BOOT, SEED)
    return dict(features="+".join(feats), n_pos=int(y.sum()), n_neg=int(len(y) - y.sum()), auc_oof=est, lo=lo, hi=hi,
                fold_auc_mean=float(np.mean(fa)), fold_auc_min=float(np.min(fa)), fold_auc_max=float(np.max(fa))), oof, y


def main():
    S = build(("match_mbs_prim", "match_sbs_prim"))
    cnt = S.grp.value_counts().rename_axis("group").reset_index(name="n_pockets")
    cnt["n_clusters"] = cnt.group.map(S.groupby("grp").unit.nunique()); cnt["n_models"] = cnt.group.map(S.groupby("grp").modelEntityId.nunique())
    cnt["n_interface_pockets_matched_only_to_SBS_or_other_excluded"] = ""
    I = S[S.is_interface]
    extra = pd.DataFrame([dict(group="interface_total", n_pockets=len(I), n_clusters=I.unit.nunique(), n_models=I.modelEntityId.nunique()),
                          dict(group="interface_matched_MBS_dcc_only", n_pockets=int((I.match_mbs_prim_dcc).sum())),
                          dict(group="interface_matched_MBS_oc_only_def", n_pockets=int((I.match_mbs_prim_oc).sum())),
                          dict(group="interface_matched_any_SBS", n_pockets=int((I.match_sbs_prim).sum())),
                          dict(group="interface_matched_MBS_and_SBS", n_pockets=int((I.match_sbs_prim & I.match_mbs_prim).sum()))])
    write_tsv(pd.concat([cnt.drop(columns=cnt.columns[-1]), extra]), OUT / "group_counts.tsv", "groups IM/IU/SM/SU")
    G = S[S.grp != "other"].copy()
    with Pool(4) as p:
        med = p.map(job_med, [(g, m, G[G.grp == g]) for g in ["IM", "IU", "SM", "SU"] for m in METRICS])
        dl = p.map(job_delta, [(a, b, m, G[G.grp.isin([a, b])]) for a, b in CONTRASTS for m in METRICS])
        # sensitivity: DCC-only matching (IU defined as before: no match under the main definition)
        S2 = build(("match_mbs_prim_dcc", "match_sbs_prim_dcc")); G2 = S2[S2.grp != "other"].copy()
        dl2 = p.map(job_delta, [(a, b, m, G2[G2.grp.isin([a, b])]) for a, b in [("IM", "IU"), ("IM", "SM")] for m in METRICS])
    write_tsv(pd.DataFrame(med, columns=["group", "metric", "median", "lo", "hi"]), OUT / "group_medians.tsv", "median, cluster-bootstrap 95% CI")
    write_tsv(pd.DataFrame(dl), OUT / "cliffs_delta.tsv", "delta = g1 vs g0 (positive: g1 larger); matching = DCC<=4A or OC>=0.5")
    write_tsv(pd.DataFrame(dl2), OUT / "cliffs_delta_dcc_only_sensitivity.tsv", "matching = DCC<=4A only")
    # classifier
    D = G[G.grp.isin(["IM", "IU"])].copy(); rows = []; keep = {}
    for feats in (["frac_charged"], ["burrow_score"], ["frac_charged", "burrow_score"], ["volume_A3"],
                  ["frac_charged", "burrow_score", "volume_A3"]):
        r, oof, y = cv_auc(D, feats); rows.append(r); keep["+".join(feats)] = (oof, y)
    write_tsv(pd.DataFrame(rows), OUT / "classifier_auc.tsv", "IM vs IU, logistic (L2=1), 5-fold CV split by paralog cluster, AUC CI = cluster bootstrap of OOF scores")
    pd.DataFrame({k + "_oof": v[0] for k, v in keep.items()}).assign(y=keep["frac_charged"][1], pocket_id=D.pocket_id.to_numpy()).to_parquet(INTERIM.parent / "iface_ligand_matched/oof.parquet")
    G.to_parquet(INTERIM.parent / "iface_ligand_matched/groups.parquet")
    print(pd.DataFrame(med, columns=["g", "m", "med", "lo", "hi"]).round(3).to_string()); 
    print(pd.DataFrame(dl).round(3).to_string()); print(pd.DataFrame(dl2).round(3).to_string()); print(pd.DataFrame(rows).round(3).to_string())
    print(cnt.to_string()); print(extra.to_string())


if __name__ == "__main__":
    main()
