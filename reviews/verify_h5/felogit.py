"""Verifier: fixed-effects (protein-stratified) logistic regression via Newton with Schur complement.
With large strata (hundreds of residues/protein) the unconditional FE-logit MLE approximates the conditional
(clogit) estimate; proteins without both outcomes carry no within-protein information and are dropped (as in clogit)."""
import numpy as np, pandas as pd

def fe_logit(y, X, g, maxit=50, tol=1e-8):
    y = np.asarray(y, float); X = np.asarray(X, float); g = np.asarray(g)
    codes, gi = np.unique(g, return_inverse=True)
    G = len(codes); k = X.shape[1]
    ev = np.bincount(gi, y, G); nn = np.bincount(gi, None, G)
    inf = (ev > 0) & (ev < nn)
    keep = inf[gi]
    y, X, gi = y[keep], X[keep], gi[keep]
    _, gi = np.unique(gi, return_inverse=True); G = gi.max() + 1
    beta = np.zeros(k); m = np.bincount(gi, y, G) / np.bincount(gi, None, G)
    alpha = np.log(m / (1 - m))
    ll_old = -np.inf
    for it in range(maxit):
        eta = alpha[gi] + X @ beta
        p = 1 / (1 + np.exp(-eta)); w = p * (1 - p); r = y - p
        ll = np.sum(y * eta - np.logaddexp(0, eta))
        ga = np.bincount(gi, r, G); gb = X.T @ r
        h = np.bincount(gi, w, G)
        Hab = np.zeros((G, k))
        for j in range(k): Hab[:, j] = np.bincount(gi, w * X[:, j], G)
        Hbb = (X * w[:, None]).T @ X
        S = Hbb - Hab.T @ (Hab / h[:, None])
        rhs = gb - Hab.T @ (ga / h)
        db = np.linalg.solve(S, rhs); da = (ga - Hab @ db) / h
        beta += db; alpha += da
        if abs(ll - ll_old) < tol and np.max(np.abs(db)) < 1e-6: break
        ll_old = ll
    cov = np.linalg.inv(S)
    return beta, np.sqrt(np.diag(cov)), dict(n=int(keep.sum()), n_strata=int(G), events=int(y.sum()), iters=it + 1, cov=cov)

def design(df, ref, levels, extra_num=("rsa_iso", "plddt"), aa=False, cat_col="cat"):
    cols = {}
    for l in levels:
        if l != ref: cols[f"{l}_vs_{ref}"] = (df[cat_col] == l).astype(float).values
    for c in extra_num:
        v = df[c].values.astype(float); cols[f"z_{c}"] = (v - np.nanmean(v)) / np.nanstd(v)
    if aa:
        for a in sorted(df.resname.unique()):
            if a != "ALA": cols[f"aa_{a}"] = (df.resname == a).astype(float).values
    return pd.DataFrame(cols)

def summarize(beta, se, names, info, label):
    out = []
    for b, s, n in zip(beta, se, names):
        out.append(dict(analysis=label, term=n, OR=np.exp(b), lo=np.exp(b - 1.96 * s), hi=np.exp(b + 1.96 * s),
                        p=float(2 * __import__("scipy.stats").stats.norm.sf(abs(b / s))), n=info["n"], strata=info["n_strata"], events=info["events"]))
    return out
