#!/usr/bin/env python
"""Q2.2 analysis 3 controls (EXPLORATORY): pocket-calling instability between isolated-chain and dimer runs.

Motivation (found in q22_compare.py output): 34% of atlas-pass isolated-chain pockets with NO lining residue at the
interface have no counterpart (overlap coefficient < 0.5) among the dimer pockets of the same coordinates, strongly
dependent on burrow_score. BURROW output is therefore context-dependent, so raw "dimer-only" and "destroyed" rates
include method noise. Controls:
 A. "New in dimer" rate of SINGLE-CHAIN dimer pockets (no residue on the partner chain) = null rate for dimer-only.
    Compared with interface pockets within burrow_score quintiles (quintile edges over the compared pockets);
    excess = weighted mean (by interface count) of within-quintile rate differences. Cluster bootstrap.
 B. Isolated-chain pocket loss with homodimer chain-swap matching (an iso pocket of chain A also matches the
    symmetric B copy in the dimer), partner-contact (>=3 lining interface residues) vs no-contact, within
    burrow_score quintiles; excess loss = contact - no-contact.
For homodimers the candidate sets in both directions include the A<->B swapped residue sets.
Overlap threshold: DIMER_ONLY_PRIMARY (0.5) with 0.2/0.8 sensitivities for A.
Outputs: results/part2/q22/q22_dimer_only_control.tsv, q22_destruction_control.tsv
Status: first attempt, not yet verified.
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import INTERIM, OUT, N_BOOT, DIMER_ONLY_THRESHOLDS, oc, write_tsv

N_Q = 5
SWAP = {"A": "B", "B": "A"}


def swap(s):
    return frozenset(SWAP[x[0]] + x[1:] for x in s)


def best_oc(query_sets, models, kinds, ref):
    out = []
    for s, m, k in zip(query_sets, models, kinds):
        cand = ref.get(m, [])
        best = max([oc(s, t) for t in cand] + [0.0])
        if k == "homo":
            ss = swap(s)
            best = max([best] + [oc(ss, t) for t in cand])
        out.append(best)
    return np.array(out)


def excess(df, flag, grp, seed=0):
    """weighted mean over burrow_score quintiles of rate(flag | grp) - rate(flag | ~grp); cluster bootstrap."""
    df = df.copy()
    df["q"] = pd.qcut(df.burrow_score.rank(method="first"), N_Q, labels=False)
    f, g, q = df[flag].to_numpy(float), df[grp].to_numpy(bool), df.q.to_numpy()
    def est(ix):
        diffs, w = [], []
        for k in range(N_Q):
            m = q[ix] == k; a = f[ix][m & g[ix]]; b = f[ix][m & ~g[ix]]
            if len(a) and len(b): diffs.append(a.mean() - b.mean()); w.append(len(a))
        return np.average(diffs, weights=w) if w else np.nan
    def raw(ix):
        return f[ix][g[ix]].mean() - f[ix][~g[ix]].mean()
    gi = df.groupby("unit").indices; keys = list(gi); rng = np.random.default_rng(seed)
    allix = np.arange(len(df)); e, r = est(allix), raw(allix)
    bs = []
    for _ in range(N_BOOT):
        ix = np.concatenate([gi[keys[i]] for i in rng.integers(0, len(keys), len(keys))]); bs.append((est(ix), raw(ix)))
    bs = np.array(bs)
    return dict(rate_g1=f[g].mean(), rate_g0=f[~g].mean(), n_g1=int(g.sum()), n_g0=int((~g).sum()),
                raw_diff=r, raw_lo=np.percentile(bs[:, 1], 2.5), raw_hi=np.percentile(bs[:, 1], 97.5),
                score_matched_diff=e, sm_lo=np.percentile(bs[:, 0], 2.5), sm_hi=np.percentile(bs[:, 0], 97.5))


def main():
    d = pd.read_parquet(INTERIM / "dimer_pocket_features_units.parquet")
    iso = pd.read_parquet(INTERIM / "iso_pocket_fate.parquet")  # atlas-pass, homo chain A only (for B)
    iso_all = pd.read_parquet(INTERIM / "iso_pocket_features.parquet")  # all iso pockets (reference sets for A)
    ref_iso = {}
    for m, r in zip(iso_all.modelEntityId, iso_all.residues):
        ref_iso.setdefault(m, []).append(frozenset(r.split()))
    ref_dim = {}
    for m, r in zip(d.modelEntityId, d.residues):
        ref_dim.setdefault(m, []).append(frozenset(r.split()))
    # ---- A
    S = d[d.atlas_pass & d.is_sym_rep & (d.is_interface | ~d.is_interface_ge1)].copy()
    S = S[S.modelEntityId.isin(ref_iso)]
    S["oc_iso_swap"] = best_oc([frozenset(r.split()) for r in S.residues], S.modelEntityId, S.kind, ref_iso)
    rows = []
    for t in DIMER_ONLY_THRESHOLDS:
        S["new"] = S.oc_iso_swap < t
        for label, sub in [("all", S), ("homo", S[S.kind == "homo"]), ("hetero", S[S.kind == "hetero"])]:
            r = excess(sub.reset_index(drop=True), "new", "is_interface")
            rows.append(dict(threshold_oc=t, stratum=label, contrast="new_in_dimer: interface(g1) vs single_chain(g0)", **r))
            print(t, label, {k: round(v, 3) for k, v in r.items()}, flush=True)
    write_tsv(pd.DataFrame(rows), OUT / "q22_dimer_only_control.tsv",
              "control A: new-in-dimer = max OC with iso pockets (homo: incl. chain swap) < threshold; score-matched = within burrow_score quintiles")
    S[["pocket_id", "oc_iso_swap"]].to_parquet(INTERIM / "dimer_oc_iso_swap.parquet")
    # ---- B
    iso = iso[iso.contact != "contact_1to2"].copy()
    iso["oc_dimer_swap"] = best_oc([frozenset(r.split()) for r in iso.residues], iso.modelEntityId, iso.kind, ref_dim)
    iso["lost_swap"] = iso.oc_dimer_swap < 0.5
    iso["grp"] = iso.contact == "partner_contact_ge3"
    rows = []
    for label, sub in [("all", iso), ("homo", iso[iso.kind == "homo"]), ("hetero", iso[iso.kind == "hetero"])]:
        r = excess(sub.reset_index(drop=True), "lost_swap", "grp")
        rows.append(dict(stratum=label, contrast="lost_in_dimer: partner_contact_ge3(g1) vs no_contact(g0)", **r))
        print("B", label, {k: round(v, 3) for k, v in r.items()}, flush=True)
    write_tsv(pd.DataFrame(rows), OUT / "q22_destruction_control.tsv",
              "control B: lost = max OC with dimer pockets (homo: incl. chain swap) < 0.5; score-matched within burrow_score quintiles")


if __name__ == "__main__":
    main()
