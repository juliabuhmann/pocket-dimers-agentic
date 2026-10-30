#!/usr/bin/env python
"""H5 step 3-5: power gate, mixed models (via scripts/h5_glmm.R), CMH robustness, descriptive rates with protein
bootstrap CIs, known-answer (TTR) and shuffled-label negative control. First attempt, not independently verified.

Usage: .venvs/part1/bin/python scripts/h5_analysis.py --stage fit --workers 8
       .venvs/part1/bin/python scripts/h5_analysis.py --stage summarize
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS / "scripts"))
from h5_build_table import provenance, CATS, AA3  # noqa: E402

TAB = WS / "data/part2/h5/tables"
RES = WS / "results/part2/h5"
FITS = RES / "fits"
RSCRIPT = WS / ".venvs/h5_r/bin/Rscript"
SEEDS = [0, 1, 2]                 # DESIGN.md §5; seed 0 primary
N_BOOT = 1000                     # DESIGN.md §5 bootstrap resamples over clusters (proteins)
POWER_MIN_PROTEINS = 30           # DESIGN.md §1 power gate
REF = "interface"
PRIMARY = "dimer_only"
TTR_ACC, TTR_MATURE_POS, TTR_WT, TTR_MUT = "P02766", 30, "V", "M"   # task spec known-answer: TTR V30M (mature numbering)

# (name, table, outcome, catcol, shuffle_seed, label)
JOBS = [("primary_plp", "primary", "plp", "cat", -1, "PRIMARY"),
        ("primary_am", "primary", "am", "cat", -1, "SECONDARY"),
        ("sens_bestmodel_plp", "sens_bestmodel", "plp", "cat", -1, "SENSITIVITY"),
        ("homo_plp", "homo", "plp", "cat", -1, "SECONDARY"),
        ("hetero_plp", "hetero", "plp", "cat", -1, "SECONDARY"),
        ("ipsae_T1_plp", "ipsae_T1", "plp", "cat", -1, "SECONDARY"),
        ("ipsae_T2_plp", "ipsae_T2", "plp", "cat", -1, "SECONDARY"),
        ("ipsae_T3_plp", "ipsae_T3", "plp", "cat", -1, "SECONDARY"),
        ("no_pdb_complex_plp", "no_pdb_complex", "plp", "cat", -1, "EXPLORATORY"),
        ("dimeronly_ov0.2_plp", "primary", "plp", "cat_s02", -1, "SENSITIVITY"),
        ("dimeronly_ov0.8_plp", "primary", "plp", "cat_s08", -1, "SENSITIVITY"),
        # PI-requested after the primary fit was seen (labelled, not pre-registered):
        ("dimeronly_burrowtop2quintiles_plp", "primary", "plp", "cat_q40", -1, "SENSITIVITY_PI"),
        ("split_preexisting_ifpocket_plp", "primary", "plp", "cat_split", -1, "COMPARATOR_PI")] + \
       [(f"shuffle_seed{s}_plp", "primary", "plp", "cat", s, "NEG_CONTROL") for s in SEEDS]


def fit(workers):
    FITS.mkdir(parents=True, exist_ok=True)

    def run(j):
        name, tab, outc, catcol, seed, _ = j
        out = FITS / f"{name}.tsv"
        if out.exists():
            return name, "cached"
        log = open(WS / f"logs/h5/fit_{name}.log", "w")
        r = subprocess.run([str(RSCRIPT), str(WS / "scripts/h5_glmm.R"), str(TAB / f"protein_residues_{tab}.parquet"),
                            str(out), outc, catcol, str(seed), "1"], stdout=log, stderr=subprocess.STDOUT,
                           env={**__import__("os").environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
        return name, r.returncode

    with ThreadPoolExecutor(workers) as ex:
        for name, rc in ex.map(run, JOBS):
            print(name, rc, flush=True)


def power_gate(g, catcol="cat"):
    rows = []
    for c in [c for c in g[catcol].cat.categories if c not in ("buried", "excluded")]:
        x = g[g[catcol] == c]
        rows.append((c, len(x), int(x.y.sum()), x.accession.nunique(), x[x.y == 1].accession.nunique()))
    return pd.DataFrame(rows, columns=["category", "n_residues", "n_plp_residues", "n_proteins", "n_proteins_with_plp"])


def boot_rates(g, seed=0, n=N_BOOT):
    """P/LP residue rate per category; 95% percentile CI from resampling proteins (clusters) with replacement."""
    g = g[g.cat != "buried"]
    acc = g.accession.astype("category")
    codes = acc.cat.codes.values; P = len(acc.cat.categories)
    out = []
    rng = np.random.default_rng(seed)
    W = rng.multinomial(P, np.full(P, 1 / P), size=n)  # bootstrap protein weights
    for c in CATS[:-1]:
        m = (g.cat == c).values
        num = np.bincount(codes[m], weights=g.y.values[m], minlength=P)
        den = np.bincount(codes[m], minlength=P).astype(float)
        est = num.sum() / den.sum()
        b = (W @ num) / (W @ den)
        am = g.am_mean.values[m]
        out.append((c, int(den.sum()), int(num.sum()), est, *np.nanpercentile(b, [2.5, 97.5]), np.nanmean(am)))
    return pd.DataFrame(out, columns=["category", "n_residues", "n_plp", "plp_rate", "ci_lo", "ci_hi", "am_mean"])


def cmh(g, a=PRIMARY, b=REF):
    """Mantel-Haenszel OR stratified by protein (unadjusted robustness check, no covariates)."""
    from statsmodels.stats.contingency_tables import StratifiedTable
    x = g[g.cat.isin([a, b])]
    tabs = []
    for _, s in x.groupby("accession"):
        t = np.array([[((s.cat == a) & (s.y == 1)).sum(), ((s.cat == a) & (s.y == 0)).sum()],
                      [((s.cat == b) & (s.y == 1)).sum(), ((s.cat == b) & (s.y == 0)).sum()]], float)
        if t.sum(1).min() > 0 and t[:, 0].sum() > 0:
            tabs.append(t)
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint()
    return dict(contrast=f"{a}_vs_{b}", n_strata=len(tabs), OR_MH=st.oddsratio_pooled, lo=lo, hi=hi,
                p=st.test_null_odds(correction=False).pvalue)


def ttr_check():
    """Known-answer: TTR V30M. Mature-numbering offset = UniProt signal-peptide length (fetched from UniProt REST)."""
    import requests
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{TTR_ACC}.json", timeout=60).json()
    sp = [f for f in r.get("features", []) if f["type"] == "Signal"]
    off = sp[0]["location"]["end"]["value"] if sp else 0
    pos = TTR_MATURE_POS + off
    seq = r["sequence"]["value"]
    m = pd.read_parquet(TAB / "model_residues.parquet", filters=[("accession", "==", TTR_ACC)])
    g = pd.read_parquet(TAB / "protein_residues_primary.parquet", filters=[("accession", "==", TTR_ACC)])
    v = pd.read_parquet(WS / "data/processed/variants_humsavar.parquet")
    v = v[(v.accession == TTR_ACC) & (v.position == pos)]
    rows = m[m.position == pos][["modelEntityId", "chain", "kind", "max_ipSAE", "resname", "plddt", "rsa_dimer", "rsa_iso",
                                 "min_dist_partner", "interface", "cat"]].copy()
    rows["cat"] = rows.cat.map(lambda c: CATS[c - 1])
    info = dict(signal_peptide_len=off, uniprot_pos=pos, uniprot_wt=seq[pos - 1], expected_wt=TTR_WT,
                humsavar=";".join(f"{a.wt}{a.position}{a.mut}:{a.category}" for a in v.itertuples()),
                aggregated_cat=str(g[g.position == pos].cat.iloc[0]) if (g.position == pos).any() else "absent",
                ttr_plp_residues=int(g.y.sum()), ttr_plp_by_cat=g.groupby("cat", observed=True).y.sum().to_dict())
    return info, rows


def summarize():
    prov = provenance()
    hdr = f"# design_sha256={prov['design_sha256']} git_commit={prov['git_commit']}\n"

    def w(df, name):
        with open(RES / name, "w") as fh:
            fh.write(hdr); df.to_csv(fh, sep="\t", index=False)

    g = pd.read_parquet(TAB / "protein_residues_primary.parquet")
    # power gate for every table/contrast
    pg = []
    for tab in ["primary", "sens_bestmodel", "homo", "hetero", "ipsae_T1", "ipsae_T2", "ipsae_T3", "no_pdb_complex"]:
        x = pd.read_parquet(TAB / f"protein_residues_{tab}.parquet")
        cols = ["cat"] + (["cat_s02", "cat_s08", "cat_q40", "cat_split"] if tab == "primary" else [])
        for cc in cols:
            p = power_gate(x, cc); p.insert(0, "catcol", cc); p.insert(0, "table", tab); pg.append(p)
    pg = pd.concat(pg)
    w(pg, "h5_power_gate.tsv")
    # model results
    fr = []
    for name, tab, outc, catcol, seed, label in JOBS:
        f = FITS / f"{name}.tsv"
        if f.exists():
            d = pd.read_csv(f, sep="\t"); d.insert(0, "analysis", name); d.insert(1, "label", label); fr.append(d)
    fr = pd.concat(fr)
    gate = pg.set_index(["table", "catcol", "category"]).n_proteins_with_plp
    def gate_ok(r):
        if not r.term.startswith("cat"):
            return np.nan
        tab = dict((j[0], j[1]) for j in JOBS)[r.analysis]
        c = r.term[3:]
        return bool(min(gate.get((tab, r.catcol, c), 0), gate.get((tab, r.catcol, REF), 0)) >= POWER_MIN_PROTEINS)
    fr["power_gate_pass"] = fr.apply(gate_ok, axis=1)
    w(fr, "h5_model_results.tsv")
    rates = boot_rates(g, seed=SEEDS[0])
    rates_r = [boot_rates(g, seed=s).assign(seed=s) for s in SEEDS]
    w(rates, "h5_category_rates.tsv"); w(pd.concat(rates_r), "h5_category_rates_seeds.tsv")
    # rates for the split tables (for the forest/descriptive)
    w(pd.DataFrame([cmh(g)] + [cmh(g, c, REF) for c in ["sc_pocket", "surface"]] + [cmh(g, PRIMARY, "sc_pocket")]), "h5_cmh.tsv")
    # single-chain pocket split (exploratory, descriptive)
    sc = g[g.cat == "sc_pocket"].groupby("sc_from_ifpocket").agg(n=("y", "size"), plp=("y", "sum"), rate=("y", "mean"), am=("am_mean", "mean")).reset_index()
    w(sc, "h5_sc_pocket_split_exploratory.tsv")
    info, rows = ttr_check()
    w(rows, "h5_ttr_known_answer.tsv")
    with open(RES / "h5_ttr_known_answer.json", "w") as fh:
        json.dump({**info, **prov}, fh, indent=1, default=str)
    print(fr[fr.term.str.startswith("cat")][["analysis", "term", "estimate", "se", "p", "power_gate_pass"]].to_string())
    print(rates.to_string()); print(info)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", required=True); ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    fit(a.workers) if a.stage == "fit" else summarize()
