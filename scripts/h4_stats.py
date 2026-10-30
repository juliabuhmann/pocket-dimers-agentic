#!/usr/bin/env python
"""H4 / Q2.1 step 3-5 statistics (first attempt, NOT verified).

Input: data/part2/h4/<census>/recovery.parquet (h4_recovery.py).
Pipeline: mapping filter (SIFTS accession + >=90% residue-name identity) -> symmetric-copy collapse (census site group)
-> superposition ok and model has dimer pockets -> [interface filter Jaccard >= 0.3 for the primary estimand]
-> one representative entry per non-redundant unit (model x cluster30 x CCD x MBS label; best resolution, then id).
H4 primary contrast: DiD = mean(hit_dimer - hit_iso | MBS) - mean(hit_dimer - hit_iso | SBS), organic ligands,
cluster-stratified bootstrap CI (1,000, seed 0; strata = homo/hetero; resampling unit = connected component of
model and 30% cluster) and a within-cluster label-permutation p (N_PERM). Raw MBS gain: bootstrap CI + cluster
sign-flip p. Hetero enters the primary test only if the feasibility gate passes (gate.json); otherwise descriptive.
Outputs: results/part2/h4/<census>/*.tsv
"""
import argparse
import json

import numpy as np
import pandas as pd

from h4_common import (OUT_DATA, OUT_RES, IFACE_JACCARD_MIN, SEEDS, N_BOOT, N_PERM, POWER_GATE_CLUSTERS,
                       AF_TRAINING_CUTOFF, WS, write_tsv)
from h4_feasibility import mapping_filter, collapse_and_units, boot_clusters

HIT_CONDS = ["dimer", "iso", "dimer_nocollapse", "iso_nocollapse", "dimer_occollapse", "iso_occollapse", "dimer_atlas", "iso_atlas",
             "dimer_shuf", "iso_shuf"] + [f"{c}_null{s}" for c in ("dimer", "iso") for s in SEEDS]


def representatives(s):
    s = s.sort_values(["unit", "resolution", "pdb_id", "ligand_uid"])
    return s.drop_duplicates("unit").copy()


def site_region_dedup(rep):
    """EXPLORATORY sensitivity: per model and MBS label keep one site per binding region (greedy by resolution;
    a site is dropped if its chain-role-aware binding-residue set has Jaccard >= 0.5 with a kept site)."""
    keep = []
    for _, g in rep.sort_values(["resolution", "pdb_id", "ligand_uid"]).groupby(["model", "is_mbs"], sort=False):
        kept = []
        for i, b1, b2 in zip(g.index, g.bind1, g.bind2):
            L = {("1", p) for p in b1} | {("2", p) for p in b2}
            if g.loc[i, "kind"] == "homo":
                L2 = {("2", p) for p in b1} | {("1", p) for p in b2}
            else:
                L2 = L
            if all(max(len(L & K) / max(len(L | K), 1), len(L2 & K) / max(len(L2 | K), 1)) < 0.5 for K in kept):
                kept.append(L); keep.append(i)
    return rep.loc[keep]


def post_cutoff_pairs(models):
    """Pair is post-cutoff if no PDB entry released before the cutoff contains the pair as a complex:
    homo: >=2 chains of the accession in the AU (SIFTS chain table) or a census homomer assembly with the accession;
    hetero: both accessions in the same AU."""
    sc = pd.read_parquet(WS / "data/processed/sifts_chain_uniprot.parquet", columns=["pdb", "chain", "sp_primary"])
    ent = pd.read_parquet(WS / "data/part1/entries.parquet", columns=["pdb_id", "release_date", "assembly_type", "uniprot_accs"])
    rel = dict(zip(ent.pdb_id, ent.release_date))
    sc["pre"] = sc.pdb.map(rel).fillna("9999") < AF_TRAINING_CUTOFF
    sc = sc[sc.pre].drop_duplicates(["pdb", "chain", "sp_primary"])
    nch = sc.groupby(["pdb", "sp_primary"]).chain.nunique()
    homo_pre = set(nch[nch >= 2].index.get_level_values(1))
    eh = ent[(ent.assembly_type == "homomer") & (ent.release_date < AF_TRAINING_CUTOFF) & ent.uniprot_accs.notna()]
    for x in eh.uniprot_accs:
        homo_pre.update(x)
    accs_by_pdb = sc.groupby("pdb").sp_primary.agg(set)
    out = {}
    need = {(m, a, b) for m, a, b in models}
    het_needed = {frozenset((a, b)) for _, a, b in need if a != b}
    het_pre = set()
    for s in accs_by_pdb:
        if len(s) < 2:
            continue
        for p in het_needed:
            if p <= s:
                het_pre.add(p)
    for m, a, b in need:
        out[m] = (a not in homo_pre) if a == b else (frozenset((a, b)) not in het_pre)
    return out


# ------------------------------------------------------------------ resampling helpers
def _cluster_index(df):
    cl, inv = np.unique(df.boot_cluster.values, return_inverse=True)
    return cl, inv


def boot_weights(df, rng, n):
    """Cluster bootstrap stratified by kind: returns (n, n_clusters) multiplicity matrix."""
    cl, inv = _cluster_index(df)
    kind_of = pd.Series(df.kind.values).groupby(inv).first().values
    W = np.zeros((n, len(cl)))
    for k in np.unique(kind_of):
        idx = np.where(kind_of == k)[0]
        draws = rng.integers(0, len(idx), size=(n, len(idx)))
        for b in range(n):
            np.add.at(W[b], idx[draws[b]], 1)
    return W, inv


def did_stats(df, a="dimer_hit", b="iso_hit", seed=0, n_boot=N_BOOT, n_perm=N_PERM):
    g = df[a].astype(float).values - df[b].astype(float).values
    m = df.is_mbs.values.astype(bool)
    res = dict(n=len(df), n_mbs=int(m.sum()), n_sbs=int((~m).sum()), n_clusters=df.boot_cluster.nunique(),
               n_clusters_mbs=df[m].boot_cluster.nunique(), n_clusters_sbs=df[~m].boot_cluster.nunique(),
               rate_a_mbs=df[a][m].mean(), rate_b_mbs=df[b][m].mean(), rate_a_sbs=df[a][~m].mean(), rate_b_sbs=df[b][~m].mean())
    res["gain_mbs"] = g[m].mean() if m.any() else np.nan
    res["gain_sbs"] = g[~m].mean() if (~m).any() else np.nan
    res["did"] = res["gain_mbs"] - res["gain_sbs"]
    res["discordant_mbs_a_only"] = int(((g > 0) & m).sum()); res["discordant_mbs_b_only"] = int(((g < 0) & m).sum())
    res["discordant_sbs_a_only"] = int(((g > 0) & ~m).sum()); res["discordant_sbs_b_only"] = int(((g < 0) & ~m).sum())
    if not m.any() or m.all():
        return res
    rng = np.random.default_rng(seed)
    W, inv = boot_weights(df, rng, n_boot)
    nc = W.shape[1]
    sg_m = np.bincount(inv, weights=g * m, minlength=nc); n_m = np.bincount(inv, weights=m.astype(float), minlength=nc)
    sg_s = np.bincount(inv, weights=g * ~m, minlength=nc); n_s = np.bincount(inv, weights=(~m).astype(float), minlength=nc)
    with np.errstate(invalid="ignore", divide="ignore"):
        gm = W @ sg_m / (W @ n_m); gs = W @ sg_s / (W @ n_s)
    d = gm - gs
    res["did_lo"], res["did_hi"] = np.nanpercentile(d, [2.5, 97.5])
    res["gain_mbs_lo"], res["gain_mbs_hi"] = np.nanpercentile(gm, [2.5, 97.5])
    res["gain_sbs_lo"], res["gain_sbs_hi"] = np.nanpercentile(gs, [2.5, 97.5])
    # within-cluster permutation of the MBS label (keeps cluster composition)
    rng = np.random.default_rng(seed)
    order = np.argsort(inv, kind="stable")
    inv_s, g_s, m_s = inv[order], g[order], m[order]
    obs = res["did"]
    null = np.empty(n_perm)
    for i in range(n_perm):
        key = np.lexsort((rng.random(len(inv_s)), inv_s))
        mp = m_s[key]  # labels permuted within cluster blocks (inv_s sorted -> blocks preserved)
        null[i] = g_s[mp].mean() - g_s[~mp].mean()
    # the pooled statistic is not centred at 0 under within-cluster exchangeability (cluster composition differs),
    # so the two-sided p is computed around the permutation mean
    mu = null.mean()
    res["did_perm_null_mean"] = mu
    res["did_p_perm_within_cluster"] = (np.sum(np.abs(null - mu) >= abs(obs - mu) - 1e-12) + 1) / (n_perm + 1)
    res["did_p_boot"] = min(1.0, 2 * min(np.nanmean(d <= 0), np.nanmean(d >= 0)))
    # raw MBS gain: cluster sign-flip
    cs = np.bincount(inv[m], weights=g[m], minlength=nc)
    flips = rng.choice([-1.0, 1.0], size=(n_perm, nc))
    null = flips @ cs / m.sum()
    res["gain_mbs_p_signflip"] = (np.sum(np.abs(null) >= abs(res["gain_mbs"]) - 1e-12) + 1) / (n_perm + 1)
    cs = np.bincount(inv[~m], weights=g[~m], minlength=nc)
    null = flips @ cs / (~m).sum()
    res["gain_sbs_p_signflip"] = (np.sum(np.abs(null) >= abs(res["gain_sbs"]) - 1e-12) + 1) / (n_perm + 1)
    res["power_gate_pass"] = bool(res["n_clusters_mbs"] >= POWER_GATE_CLUSTERS and res["n_clusters_sbs"] >= POWER_GATE_CLUSTERS)
    return res


def rate_ci(df, col, seed=0, n_boot=N_BOOT):
    x = df[col].astype(float).values
    if not len(x):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    W, inv = boot_weights(df, rng, n_boot)
    nc = W.shape[1]
    s = np.bincount(inv, weights=x, minlength=nc); n = np.bincount(inv, minlength=nc).astype(float)
    with np.errstate(invalid="ignore"):
        r = W @ s / (W @ n)
    lo, hi = np.nanpercentile(r, [2.5, 97.5])
    return x.mean(), lo, hi


def prepare(census, ledger):
    r = pd.read_parquet(OUT_DATA / census / "recovery.parquet")
    for c in [c for c in r.columns if c.endswith("_hit")]:
        r[c] = r[c].fillna(False).astype(bool)
    r = mapping_filter(r, ledger)
    r = collapse_and_units(r)
    ledger.append(("after symmetric-copy collapse", len(r)))
    r = r[r.model_has_dimer_pockets & (r.n_assign_ok > 0)]
    ledger.append(("superposition ok and dimer pockets present", len(r)))
    pc = post_cutoff_pairs(set(zip(r.model, r.acc1, r.acc2)))
    r["post_cutoff_pair"] = r.model.map(pc)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    a = ap.parse_args()
    od = OUT_RES / a.census
    gate = json.load(open(od / "gate.json"))
    led = []
    r = prepare(a.census, led)
    sets = {}
    for filt in ("filtered", "unfiltered"):
        x = r if filt == "unfiltered" else r[r.iface_jaccard >= IFACE_JACCARD_MIN]
        rep = representatives(x)
        rep["boot_cluster"] = boot_clusters(rep)
        sets[filt] = rep
        for (kind, ls, mbs), g in rep.groupby(["kind", "lig_set", "is_mbs"]):
            led.append((f"{filt}: representatives {kind} {ls} {'MBS' if mbs else 'SBS'}", len(g)))
    # interface filter removal counts (units)
    u_all = representatives(r)
    u_f = sets["filtered"]
    rm = []
    for (kind, ls, mbs), g in u_all.groupby(["kind", "lig_set", "is_mbs"]):
        kept = u_f[(u_f.kind == kind) & (u_f.lig_set == ls) & (u_f.is_mbs == mbs)]
        rm.append(dict(kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", units_before=len(g), units_after=len(kept),
                       units_removed=len(g) - len(kept),
                       rows_jaccard_lt_min=int(((r.kind == kind) & (r.lig_set == ls) & (r.is_mbs == mbs) & (r.iface_jaccard < IFACE_JACCARD_MIN)).sum())))
    write_tsv(pd.DataFrame(rm), od / "interface_filter.tsv", {"census": a.census, "jaccard_min": IFACE_JACCARD_MIN})
    write_tsv(pd.DataFrame(led, columns=["step", "n"]), od / "analysis_ledger.tsv", {"census": a.census})
    for k, v in sets.items():
        v.drop(columns=["bind"]).to_parquet(OUT_DATA / a.census / f"benchmark_{k}.parquet")

    # ---- hit rates with cluster bootstrap CIs
    rows = []
    for filt, rep in sets.items():
        for (kind, ls, mbs), g in rep.groupby(["kind", "lig_set", "is_mbs"]):
            for c in HIT_CONDS:
                m, lo, hi = rate_ci(g, f"{c}_hit")
                rows.append(dict(filter=filt, kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", condition=c, n=len(g),
                                 n_clusters=g.boot_cluster.nunique(), hit_rate=m, ci_lo=lo, ci_hi=hi))
    rates = pd.DataFrame(rows)
    write_tsv(rates, od / "hit_rates.tsv", {"census": a.census, "dcc": 4.0, "topk": "n+2"})

    # ---- residue metrics (secondary) incl. AFDB monomer atlas (step 5)
    rows = []
    for filt, rep in sets.items():
        for (kind, ls, mbs), g in rep.groupby(["kind", "lig_set", "is_mbs"]):
            ga = g[g.atlas_has_model.fillna(False).astype(bool)]
            for c in ["dimer", "iso", "dimer_atlas", "iso_atlas", "atlas"]:
                src = ga if c == "atlas" else g
                for met in ["recall", "oc", "f1", "pksize"]:
                    col = f"{c}_{met}"
                    if col not in src:
                        continue
                    rows.append(dict(filter=filt, kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", condition=c, metric=met,
                                     n=src[col].notna().sum(), mean=src[col].mean(), median=src[col].median()))
            # paired on the subset with an atlas model
            for c in ["dimer", "iso"]:
                rows.append(dict(filter=filt, kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", condition=f"{c}[atlas subset]",
                                 metric="f1", n=len(ga), mean=ga[f"{c}_f1"].mean(), median=ga[f"{c}_f1"].median()))
            for sd in SEEDS:
                for c in ["dimer", "iso"]:
                    rows.append(dict(filter=filt, kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", condition=f"{c}_null{sd}",
                                     metric="f1", n=len(g), mean=g[f"{c}_null{sd}_f1"].mean(), median=g[f"{c}_null{sd}_f1"].median()))
    write_tsv(pd.DataFrame(rows), od / "residue_metrics.tsv", {"census": a.census})

    # ---- H4 tests
    kinds_primary = ["homo"] + (["hetero"] if gate["hetero_inferential"] else [])
    tests = []

    def add(label, df, role, **kw):
        if df.is_mbs.nunique() < 2:
            return
        res = did_stats(df, n_perm=a.n_perm, **kw)
        tests.append(dict(analysis=label, role=role, **res))

    F, U = sets["filtered"], sets["unfiltered"]
    org = lambda d: d[d.lig_set == "organic"]
    P = org(F)[org(F).kind.isin(kinds_primary)]
    add("H4 primary: organic, interface-filtered, kinds=" + "+".join(kinds_primary), P, "primary")
    for sd in SEEDS[1:]:
        add(f"H4 primary, bootstrap/permutation seed {sd}", P, "robustness_seed", seed=sd)
    for kind in ("homo", "hetero"):
        add(f"organic, filtered, {kind} only", org(F)[org(F).kind == kind],
            "secondary" if kind in kinds_primary else "descriptive (gate)")
    add("organic, UNFILTERED, kinds=" + "+".join(kinds_primary), org(U)[org(U).kind.isin(kinds_primary)], "reported (unfiltered)")
    add("sensitivity: no symmetric pocket collapse", P, "sensitivity", a="dimer_nocollapse_hit", b="iso_nocollapse_hit")
    add("sensitivity: overlap-coefficient collapse (annotate_pocket_chains.py rule)", P, "sensitivity", a="dimer_occollapse_hit", b="iso_occollapse_hit")
    add("sensitivity: atlas pocket filter (score>=0.03, pLDDT>=70)", P, "sensitivity", a="dimer_atlas_hit", b="iso_atlas_hit")
    add("EXPLORATORY sensitivity: one site per binding region per model", site_region_dedup(P), "sensitivity (exploratory)")
    add("negative control: shuffled ligand centroid", P, "control", a="dimer_shuf_hit", b="iso_shuf_hit")
    for sd in SEEDS:
        add(f"null: random surface patches seed {sd}", P, "control", a=f"dimer_null{sd}_hit", b=f"iso_null{sd}_hit")
    add("metal ions/clusters (separate), filtered", F[(F.lig_set == "metal") & F.kind.isin(kinds_primary)], "descriptive (metals)")
    # strata
    for v, lab in ((True, "post-cutoff pair"), (False, "pre-cutoff pair")):
        add(f"stratum: {lab} (cutoff {AF_TRAINING_CUTOFF})", P[P.post_cutoff_pair == v], "stratum")
    mods = P.drop_duplicates("model")
    cuts = np.quantile(mods.max_ipSAE, [1 / 3, 2 / 3])
    P = P.assign(ipsae_tertile=np.digitize(P.max_ipSAE, cuts) + 1)
    for t in (1, 2, 3):
        add(f"stratum: ipSAE tertile {t} (cuts {cuts[0]:.3f}, {cuts[1]:.3f})", P[P.ipsae_tertile == t], "stratum")
    tt = pd.DataFrame(tests)
    write_tsv(tt, od / "h4_tests.tsv", {"census": a.census, "n_boot": N_BOOT, "n_perm": a.n_perm,
                                        "af_cutoff": AF_TRAINING_CUTOFF, "kinds_primary": "+".join(kinds_primary)})
    pd.set_option("display.width", 250)
    print(pd.DataFrame(led).to_string())
    print(pd.DataFrame(rm).to_string())
    print(tt[["analysis", "n", "n_mbs", "n_clusters", "rate_a_mbs", "rate_b_mbs", "rate_a_sbs", "rate_b_sbs", "gain_mbs", "gain_sbs",
              "did", "did_lo", "did_hi", "did_p_boot", "did_perm_null_mean", "did_p_perm_within_cluster", "gain_mbs_p_signflip"]].to_string())


if __name__ == "__main__":
    main()
