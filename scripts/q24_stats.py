#!/usr/bin/env python
"""Q2.4 (EXPLORATORY) summaries of data/interim/q24/q24_pocket_scores.parquet (scripts/q24_conservation.py).

Unit of observation: (human pocket, species). Conserved = score >= CONSERVED_OC on the chosen target set (primary `oc_atlas`).
CIs: cluster bootstrap over human protein clusters (Q2.2 units), N_BOOT resamples, seed 0; one resampling matrix is shared by
all statistics (same resamples for real and null and for both classes, so differences are paired).
Quintile-matched difference (as in Q2.2): sum_q w_q (p_a,q - p_b,q), w_q = share of class-a pockets in burrow_score quintile q.
Status: first attempt, not yet verified.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS / "scripts"))
from q22_common import write_tsv  # noqa: E402
from q24_conservation import CONSERVED_OC, N_BOOT, SEED, INTERIM, OUT, SPECIES  # noqa: E402

GROUPS = {"mouse": "mammals", "rat": "mammals", "zebrafish": "fish", "fly": "invertebrates", "worm": "invertebrates",
          "yeast": "fungi", "pombe": "fungi", "arabidopsis": "plant/amoeba", "dicty": "plant/amoeba",
          "ecoli": "bacteria", "mtb": "bacteria", "paeruginosa": "bacteria"}  # grouping by TimeTree divergence ties (q24_divergence.tsv)


class Boot:
    def __init__(self, units):
        self.codes, uniq = pd.factorize(units)
        rng = np.random.default_rng(SEED)
        n = len(uniq)
        self.W = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(N_BOOT)]).astype(np.float32)
        self.n = n

    def sums(self, mask, x):
        """per-bootstrap weighted sum of x over rows in mask -> (point, boot vector)."""
        u = np.bincount(self.codes[mask], weights=x[mask], minlength=self.n)
        return u.sum(), self.W @ u


def ci(v):
    return np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5)


def frac(B, mask, y):
    k, kb = B.sums(mask, y); n, nb = B.sums(mask, np.ones_like(y))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (k / n if n else np.nan), kb / nb, int(n)


def matched_diff(B, R, ma, mb, y):
    """quintile-matched difference p_a - p_b with weights from class a."""
    est, boot, w_n, w_b = 0.0, 0.0, 0.0, 0.0
    for q in sorted(R.score_q.dropna().unique()):
        mq = (R.score_q == q).values
        pa, pab, na = frac(B, ma & mq, y); pb, pbb, nb_ = frac(B, mb & mq, y)
        if na == 0 or nb_ == 0:
            continue
        _, wb = B.sums(ma & mq, np.ones_like(y))
        est += na * (pa - pb); w_n += na
        boot = boot + wb * (np.nan_to_num(pab) - np.nan_to_num(pbb)); w_b = w_b + wb
    return est / w_n, boot / w_b


def main():
    R = pd.read_parquet(INTERIM / "q24_pocket_scores.parquet")
    R["group"] = R.species.map(GROUPS)
    B = Boot(R.unit.values)
    out = []
    for col, tag in (("oc_atlas", "primary(atlas_pass targets)"), ("oc_all", "sensitivity(all targets)"), ("oc_atlas_itf", "interface target only")):
        y = (R[col].values >= CONSERVED_OC).astype(float)
        yn = (R["null_" + col].values >= CONSERVED_OC).astype(float)
        for scope, smask in [(s, (R.species == s).values) for s in SPECIES] + [(g, (R.group == g).values) for g in dict.fromkeys(GROUPS.values())] + [("ALL", np.ones(len(R), bool))] \
                + [(f"{g}_{k}", (((R.group == g) | (g == "ALL")) & (R.kind == k)).values) for g in list(dict.fromkeys(GROUPS.values())) + ["ALL"] for k in ("homo", "hetero")]:
            if smask.sum() == 0:
                continue
            res = {}
            for cl, cm in (("interface", R.cls == "interface"), ("single", R.cls == "single"), ("dimer_only", R["sub"] == "dimer_only"), ("preexisting", R["sub"] == "preexisting")):
                m = smask & cm.values
                p, pb, n = frac(B, m, y); q, qb, _ = frac(B, m, yn)
                res[cl] = (p, pb, q, qb, m)
                out.append(dict(target=tag, scope=scope, contrast=cl, n=n, n_units=len(np.unique(B.codes[m])), estimate=p, lo=ci(pb)[0], hi=ci(pb)[1],
                                null=q, null_lo=ci(qb)[0], null_hi=ci(qb)[1], excess=p - q, excess_lo=ci(pb - qb)[0], excess_hi=ci(pb - qb)[1]))
            for a, b in (("interface", "single"), ("dimer_only", "preexisting")):
                pa, pab, qa, qab, ma = res[a]; pb_, pbb, qb_, qbb, mb = res[b]
                d = pab - pbb; dn = (pab - qab) - (pbb - qbb)
                md, mdb = matched_diff(B, R, ma, mb, y)
                out.append(dict(target=tag, scope=scope, contrast=f"{a}-{b}", n=int(ma.sum() + mb.sum()), estimate=pa - pb_, lo=ci(d)[0], hi=ci(d)[1],
                                excess=(pa - qa) - (pb_ - qb_), excess_lo=ci(dn)[0], excess_hi=ci(dn)[1],
                                quintile_matched=md, qm_lo=ci(mdb)[0], qm_hi=ci(mdb)[1]))
    T = pd.DataFrame(out)
    T["divergence_mya"] = T.scope.map(R.groupby("species").divergence_mya.first())
    write_tsv(T, OUT / "q24_conservation.tsv", "conserved=score>=0.5; excess=real-null; quintile_matched=burrow_score quintile matched difference")

    # strata (pooled over species within each group, and ALL): class x burrow quintile, class x identity bin
    y = (R.oc_atlas.values >= CONSERVED_OC).astype(float); yn = (R.null_oc_atlas.values >= CONSERVED_OC).astype(float)
    st = []
    for scope, smask in [(g, (R.group == g).values) for g in dict.fromkeys(GROUPS.values())] + [("ALL", np.ones(len(R), bool))]:
        for strat in ("score_q", "id_bin"):
            for lev in [l for l in R[strat].cat.categories] if hasattr(R[strat], "cat") else sorted(R[strat].dropna().unique()):
                lm = (R[strat] == lev).values
                for cl in ("interface", "single", "dimer_only", "preexisting"):
                    cm = ((R.cls == cl) | (R["sub"] == cl)).values
                    m = smask & lm & cm
                    if m.sum() == 0:
                        continue
                    p, pb, n = frac(B, m, y); q, qb, _ = frac(B, m, yn)
                    st.append(dict(scope=scope, stratum=strat, level=str(lev), cls=cl, n=n, estimate=p, lo=ci(pb)[0], hi=ci(pb)[1], null=q, null_lo=ci(qb)[0], null_hi=ci(qb)[1]))
    write_tsv(pd.DataFrame(st), OUT / "q24_conservation_strata.tsv", "primary target set (atlas_pass); pooled over species in scope")

    # lining-residue sequence identity
    li = []
    ok = R.lining_identity.notna().values
    for scope, smask in [(s, (R.species == s).values) for s in SPECIES] + [(g, (R.group == g).values) for g in dict.fromkeys(GROUPS.values())] + [("ALL", np.ones(len(R), bool))]:
        res = {}
        for cl in ("interface", "single", "dimer_only", "preexisting"):
            m = smask & ok & ((R.cls == cl) | (R["sub"] == cl)).values
            if m.sum() == 0:
                continue
            for metric, x in (("lining_identity", R.lining_identity.values), ("lining_minus_chain", (R.lining_identity - R.chain_identity).values),
                              ("frac_aligned", R.frac_aligned.values)):
                s, sb = B.sums(m, np.nan_to_num(x)); n, nb = B.sums(m, np.ones(len(R)))
                res[(cl, metric)] = (s / n, sb / nb)
                li.append(dict(scope=scope, cls=cl, metric=metric, n=int(n), mean=s / n, lo=ci(sb / nb)[0], hi=ci(sb / nb)[1]))
        for a, b in (("interface", "single"), ("dimer_only", "preexisting")):
            for metric in ("lining_identity", "lining_minus_chain", "frac_aligned"):
                if (a, metric) in res and (b, metric) in res:
                    d = res[(a, metric)][1] - res[(b, metric)][1]
                    li.append(dict(scope=scope, cls=f"{a}-{b}", metric=metric, mean=res[(a, metric)][0] - res[(b, metric)][0], lo=ci(d)[0], hi=ci(d)[1]))
    write_tsv(pd.DataFrame(li), OUT / "q24_lining_identity.tsv", "mean per (pocket, species); identity over aligned lining residues")
    print(T[(T.target.str.startswith("primary"))].round(3).to_string())


if __name__ == "__main__":
    main()
