#!/usr/bin/env python
"""Q2.2 analyses 1-3 (EXPLORATORY, DESIGN.md v3 §3; not pre-registered).

1. Interface vs single-chain pockets in the same dimer models (atlas filter; homodimer symmetric copies collapsed
   via is_sym_rep): medians, median difference and Cliff's delta with cluster-bootstrap 95% CIs; all / homo / hetero /
   ipSAE tertiles (tertile edges over models).
2. Expected vs observed interface-pocket rate vs the interface fraction of the solvent-accessible surface.
3. Dimer-only fraction among interface pockets (overlap thresholds 0.2/0.5/0.8), dimer-only vs pre-existing interface
   pocket properties, and isolated-chain pocket loss/shrinkage in the dimer (occlusion by the partner).
Bootstrap unit: paralog cluster (connected component of the mmseqs2 graph with >=30% identity AND >=80% coverage of
both sequences, scripts/q22_paralogs.sh) of chain A's protein. The coverage requirement avoids the single-linkage giant
component (2,172 proteins / 31% of models) that forms through shared domains without it; falls back to chain-A accession if the paralog file is missing (logged).
Outputs: results/part2/q22/q22_{props,props_dimer_only,expected_rate,dimer_only_frac,destruction,n_steps}.tsv
Status: first attempt, not yet verified.
"""
import sys
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import (INTERIM, OUT, SEEDS, N_BOOT, DIMER_ONLY_THRESHOLDS, DIMER_ONLY_PRIMARY, PARALOG_MIN_ID,
                        PARALOG_MAX_EVALUE, oc, write_tsv, ipsae_tertile, cliffs_delta)

UNIT_MIN_COV = 0.8  # coverage of both sequences for bootstrap-unit clustering (our choice, see docstring)
METRICS = ["volume_A3", "site_volume_A3", "burrow_score", "plddt_mean", "buriedness", "frac_hydrophobic",
           "frac_aromatic", "frac_charged", "n_lining"]
N_STEPS = []


def step(name, n, note=""):
    N_STEPS.append((name, n, note)); print(name, n, note, flush=True)


def paralog_clusters(accs):
    p = INTERIM / "paralogs_mmseqs.tsv"
    parent = {a: a for a in accs}
    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    if not p.exists():
        step("bootstrap_unit", 0, "FALLBACK: paralog file missing, unit = chain-A accession")
        return {a: a for a in accs}
    h = pd.read_csv(p, sep="\t", header=None, names=["q", "t", "fident", "alnlen", "evalue", "qcov", "tcov"])
    h = h[(h.q != h.t) & (h.fident >= PARALOG_MIN_ID) & (h.evalue <= PARALOG_MAX_EVALUE)
          & (h.qcov >= UNIT_MIN_COV) & (h.tcov >= UNIT_MIN_COV)]
    for q, t in zip(h.q, h.t):
        if q in parent and t in parent:
            a, b = f(q), f(t)
            if a != b: parent[a] = b
    cl = {a: f(a) for a in accs}
    sz = pd.Series(cl).value_counts()
    step("bootstrap_unit", len(sz), f"paralog components (>=30% id, >=80% cov both) over dimer accessions; largest={sz.iloc[0]}")
    return cl


def boot_compare(args):
    """Group contrast g1 (True) vs g0 (False) for all METRICS, cluster bootstrap."""
    df, label, seed, metrics = args
    rng = np.random.default_rng(seed)
    grp = df.groupby("unit").indices
    keys = list(grp); idxs = [grp[k] for k in keys]
    vals = {m: df[m].to_numpy(float) for m in metrics}
    g = df["grp"].to_numpy(bool)
    def stats(ix):
        out = {}
        gg = g[ix]
        for m in metrics:
            v = vals[m][ix]; a, b = v[gg], v[~gg]
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            out[m] = (np.median(a) - np.median(b) if len(a) and len(b) else np.nan, cliffs_delta(a, b))
        return out
    est = stats(np.arange(len(df)))
    bs = {m: [] for m in metrics}
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(keys), len(keys))
        s = stats(np.concatenate([idxs[i] for i in pick]))
        for m in metrics: bs[m].append(s[m])
    rows = []
    for m in metrics:
        arr = np.array(bs[m], float)
        a = df.loc[df.grp, m].dropna(); b = df.loc[~df.grp, m].dropna()
        rows.append(dict(stratum=label, seed=seed, metric=m, n_g1=len(a), n_g0=len(b), n_units=len(keys),
                         median_g1=a.median(), median_g0=b.median(),
                         median_diff=est[m][0], md_lo=np.nanpercentile(arr[:, 0], 2.5), md_hi=np.nanpercentile(arr[:, 0], 97.5),
                         cliffs_delta=est[m][1], cd_lo=np.nanpercentile(arr[:, 1], 2.5), cd_hi=np.nanpercentile(arr[:, 1], 97.5)))
    return rows


def strata_jobs(df, metrics, primary_only_seeds=True):
    jobs = []
    for label, sub in [("all", df), ("homo", df[df.kind == "homo"]), ("hetero", df[df.kind == "hetero"])] + \
                      [(f"ipSAE_{t}", df[df.tert == t]) for t in ["T1_low", "T2_mid", "T3_high"]]:
        if sub.grp.nunique() < 2:
            continue
        for s in (SEEDS if label == "all" else SEEDS[:1]):
            jobs.append((sub.reset_index(drop=True), label, s, metrics))
    return jobs


def boot_rate(df, num, den, unit="unit", seed=0):
    """pooled ratio sum(num)/sum(den) with cluster bootstrap"""
    rng = np.random.default_rng(seed)
    g = df.groupby(unit)[[num, den]].sum()
    n, d = g[num].to_numpy(float), g[den].to_numpy(float)
    est = n.sum() / d.sum()
    bs = []
    for _ in range(N_BOOT):
        pk = rng.integers(0, len(g), len(g)); bs.append(n[pk].sum() / d[pk].sum())
    return est, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def main():
    d = pd.read_parquet(INTERIM / "dimer_pocket_features.parquet")
    iso = pd.read_parquet(INTERIM / "iso_pocket_features.parquet")
    res = pd.read_parquet(INTERIM / "residues.parquet")
    step("dimer_pockets_all", len(d)); step("dimer_models", d.modelEntityId.nunique())
    accs = sorted(set(d.uniprot_A) | set(d.uniprot_B))
    cl = paralog_clusters(accs)
    d["unit"] = d.uniprot_A.map(cl)
    models = d.drop_duplicates("modelEntityId")[["modelEntityId", "max_ipSAE", "kind", "uniprot_A"]]
    _, edges = ipsae_tertile(models.max_ipSAE)
    step("ipSAE_tertile_edges", 0, f"{edges[0]:.4f},{edges[1]:.4f}")
    d["tert"], _ = ipsae_tertile(d.max_ipSAE, edges)

    # ---------- 1. interface vs single-chain ----------
    S = d[d.atlas_pass & d.is_sym_rep].copy()
    step("atlas_pass", int(d.atlas_pass.sum())); step("atlas_pass_sym_rep", len(S))
    S["single_chain"] = ~S.is_interface_ge1
    marginal = S.is_interface_ge1 & ~S.is_interface
    step("excluded_marginal_1to2_res_on_a_chain", int(marginal.sum()))
    A = S[~marginal].copy(); A["grp"] = A.is_interface
    step("interface_primary", int(A.grp.sum())); step("single_chain", int((~A.grp).sum()))
    with Pool(7) as p:
        rows = sum(p.map(boot_compare, strata_jobs(A, METRICS)), [])
    props = pd.DataFrame(rows); props.insert(0, "contrast", "interface(g1)_vs_single_chain(g0)")
    write_tsv(props, OUT / "q22_props.tsv", "analysis1")

    # ---------- 3a. dimer-only fractions ----------
    I = S[S.is_interface & S.max_overlap_iso.notna()].copy()
    step("interface_with_iso_overlap", len(I))
    fr = []
    for t in DIMER_ONLY_THRESHOLDS:
        I["do"] = (I.max_overlap_iso < t).astype(float); I["one"] = 1.0
        for label, sub in [("all", I), ("homo", I[I.kind == "homo"]), ("hetero", I[I.kind == "hetero"])] + \
                          [(f"ipSAE_{x}", I[I.tert == x]) for x in ["T1_low", "T2_mid", "T3_high"]]:
            e, lo, hi = boot_rate(sub, "do", "one")
            fr.append(dict(threshold_oc=t, primary=t == DIMER_ONLY_PRIMARY, stratum=label, n_interface=len(sub),
                           n_dimer_only=int(sub["do"].sum()), frac=e, lo=lo, hi=hi))
    write_tsv(pd.DataFrame(fr), OUT / "q22_dimer_only_frac.tsv", "analysis3: dimer-only = interface & max_overlap_iso < threshold")

    # ---------- 3b. dimer-only vs pre-existing interface pockets ----------
    I["grp"] = I.max_overlap_iso < DIMER_ONLY_PRIMARY
    with Pool(7) as p:
        rows = sum(p.map(boot_compare, strata_jobs(I, METRICS)), [])
    pdo = pd.DataFrame(rows); pdo.insert(0, "contrast", "dimer_only(g1)_vs_preexisting_interface(g0)")
    write_tsv(pdo, OUT / "q22_props_dimer_only.tsv", "analysis3")

    # ---------- 2. expected vs observed interface rate ----------
    res["iface_surf"] = res.is_interface & res.is_surface
    rs = res.groupby("modelEntityId").agg(n_surf=("is_surface", "sum"), n_iface_surf=("iface_surf", "sum"),
                                          sasa_iso=("sasa_iso", "sum"), sasa_cpx=("sasa_cpx", "sum"))
    rs["buried_sasa_frac"] = 1 - rs.sasa_cpx / rs.sasa_iso  # fraction of chain-alone SASA buried by the partner
    # residue-level: surface residues lining >=1 atlas-pass dimer pocket
    P = d[d.atlas_pass]
    tok = P[["modelEntityId", "residues"]].assign(tok=P.residues.str.split()).explode("tok")
    tok["chain"] = tok.tok.str[0]; tok["resnum"] = pd.to_numeric(tok.tok.str[2:], errors="coerce")
    tok = tok.dropna(subset=["resnum"])
    lin = set(zip(tok.modelEntityId, tok.chain, tok.resnum.astype(int)))
    res["lines_pocket"] = [(m, c, int(n)) in lin for m, c, n in zip(res.modelEntityId, res.chain, res.resnum)]
    sr = res[res.is_surface].assign(iface_lining=lambda x: x.lines_pocket & x.is_interface)
    rl = sr.groupby("modelEntityId").agg(n_surf_lining=("lines_pocket", "sum"), n_iface_lining=("iface_lining", "sum"))
    pk_all = d[d.atlas_pass].groupby("modelEntityId").agg(n_pk=("is_interface", "size"), n_pk_iface=("is_interface", "sum"))
    pk_rep = S.groupby("modelEntityId").agg(n_pk_rep=("is_interface", "size"), n_pk_iface_rep=("is_interface", "sum"))
    M = models.set_index("modelEntityId").join([rs, rl, pk_all, pk_rep]).fillna(0)
    M["unit"] = M.uniprot_A.map(cl); M["tert"], _ = ipsae_tertile(M.max_ipSAE, edges)
    M = M[M.n_surf > 0]
    step("models_expected_rate", len(M))
    er = []
    for label, sub in [("all", M), ("homo", M[M.kind == "homo"]), ("hetero", M[M.kind == "hetero"])] + \
                      [(f"ipSAE_{x}", M[M.tert == x]) for x in ["T1_low", "T2_mid", "T3_high"]]:
        sub = sub.reset_index()
        for name, num, den in [("expected_iface_frac_of_surface_residues", "n_iface_surf", "n_surf"),
                               ("observed_iface_frac_of_pockets_all_copies", "n_pk_iface", "n_pk"),
                               ("observed_iface_frac_of_pockets_sym_collapsed", "n_pk_iface_rep", "n_pk_rep"),
                               ("observed_iface_frac_of_pocket_lining_surface_residues", "n_iface_lining", "n_surf_lining")]:
            e, lo, hi = boot_rate(sub, num, den)
            er.append(dict(stratum=label, quantity=name, n_models=len(sub), estimate=e, lo=lo, hi=hi))
        # enrichment ratios (bootstrap of ratio of pooled ratios)
        rng = np.random.default_rng(0); g = sub.groupby("unit")[["n_iface_surf", "n_surf", "n_pk_iface", "n_pk",
                                                                  "n_iface_lining", "n_surf_lining"]].sum().to_numpy(float)
        def ratios(x):
            exp = x[:, 0].sum() / x[:, 1].sum()
            return x[:, 2].sum() / x[:, 3].sum() / exp, x[:, 4].sum() / x[:, 5].sum() / exp
        est = ratios(g); bs = np.array([ratios(g[rng.integers(0, len(g), len(g))]) for _ in range(N_BOOT)])
        for k, nm in enumerate(["enrichment_pocket_level(obs_all_copies/expected)", "enrichment_residue_level(lining/expected)"]):
            er.append(dict(stratum=label, quantity=nm, n_models=len(sub), estimate=est[k],
                           lo=np.percentile(bs[:, k], 2.5), hi=np.percentile(bs[:, k], 97.5)))
        # residue-level odds ratio: surface residue lines a pocket ~ interface
        a = sub.n_iface_lining.sum(); b = (sub.n_iface_surf - sub.n_iface_lining).sum()
        c = (sub.n_surf_lining - sub.n_iface_lining).sum(); dd = (sub.n_surf - sub.n_iface_surf - (sub.n_surf_lining - sub.n_iface_lining)).sum()
        er.append(dict(stratum=label, quantity="rate_lining_iface_surface_res", n_models=len(sub), estimate=a / (a + b), lo=np.nan, hi=np.nan))
        er.append(dict(stratum=label, quantity="rate_lining_noniface_surface_res", n_models=len(sub), estimate=c / (c + dd), lo=np.nan, hi=np.nan))
    write_tsv(pd.DataFrame(er), OUT / "q22_expected_rate.tsv", "analysis2: surface = rel SASA(chain alone) >= 0.25; interface = heavy atom <= 5 A of partner")
    M.reset_index().to_parquet(INTERIM / "model_surface_rates.parquet")

    # ---------- 3c. isolated-chain pocket loss in the dimer ----------
    iso = iso[iso.atlas_pass & ~((iso.kind == "homo") & (iso.chain != "A"))].copy()  # homodimers: chain A only (copies)
    step("iso_pockets_atlas_pass_(homo_chainA_only)", len(iso))
    dsets = {}
    for m, rsd, vol in zip(d.modelEntityId, d.residues, d.volume_A3):
        dsets.setdefault(m, []).append((frozenset(rsd.split()), vol))
    best_oc, vr = [], []
    for m, rsd, vol in zip(iso.modelEntityId, iso.residues, iso.volume_A3):
        s = frozenset(rsd.split()); cand = dsets.get(m, [])
        if not cand:
            best_oc.append(0.0); vr.append(0.0); continue
        ocs = [oc(s, t) for t, _ in cand]; k = int(np.argmax(ocs))
        best_oc.append(ocs[k]); vr.append(cand[k][1] / vol if vol > 0 else np.nan)
    iso["best_oc_dimer"] = best_oc; iso["vol_ratio_dimer"] = vr
    iso["fate"] = np.where(iso.best_oc_dimer < 0.5, "lost", np.where(iso.vol_ratio_dimer < 0.5, "shrunk", "retained"))
    iso["contact"] = np.select([iso.n_lining_interface >= 3, iso.n_lining_interface == 0], ["partner_contact_ge3", "no_contact"], "contact_1to2")
    iso["unit"] = iso.modelEntityId.map(d.drop_duplicates("modelEntityId").set_index("modelEntityId").uniprot_A).map(cl)
    iso["tert"] = iso.modelEntityId.map(models.set_index("modelEntityId").max_ipSAE).pipe(lambda x: ipsae_tertile(x, edges)[0])
    de = []
    for label, sub in [("all", iso), ("homo", iso[iso.kind == "homo"]), ("hetero", iso[iso.kind == "hetero"])]:
        for ct in ["partner_contact_ge3", "contact_1to2", "no_contact"]:
            s2 = sub[sub.contact == ct].copy(); s2["one"] = 1.0
            for fate in ["lost", "shrunk", "retained"]:
                s2["x"] = (s2.fate == fate).astype(float)
                e, lo, hi = boot_rate(s2, "x", "one")
                de.append(dict(stratum=label, contact=ct, fate=fate, n_iso_pockets=len(s2), n=int(s2.x.sum()), frac=e, lo=lo, hi=hi))
    write_tsv(pd.DataFrame(de), OUT / "q22_destruction.tsv",
              "analysis3c: lost = best overlap coeff with any dimer pocket < 0.5; shrunk = OC>=0.5 & volume ratio < 0.5")
    iso.to_parquet(INTERIM / "iso_pocket_fate.parquet")
    d.to_parquet(INTERIM / "dimer_pocket_features_units.parquet")
    write_tsv(pd.DataFrame(N_STEPS, columns=["step", "n", "note"]), OUT / "q22_n_steps.tsv", "q22_compare")


if __name__ == "__main__":
    main()
