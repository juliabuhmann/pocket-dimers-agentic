#!/usr/bin/env python
"""Part 1 analysis: descriptives, H1-H3 (DESIGN.md v3 §1 family A, §2, §5) + sensitivities.

Status: first attempt, UNVERIFIED until independent re-derivation. Numbers are produced by this script only.

Rerun (after census and scripts/build_druglike_flags.py):
    source scripts/env_threads.sh
    .venvs/part1/bin/python scripts/part1_analysis.py           # -> results/part1/analysis/*.tsv
    .venvs/part1/bin/python scripts/part1_figures.py            # -> results/part1/figures/*.png

Inputs: data/part1/{ligand_sites,entries,ligand_classes}.parquet, data/processed/ccd_druglike.parquet.

Operationalisation of DESIGN §0/§1 on census columns (choices are ours where the design is silent; see docs/PART1_RESULTS.md):
- Binding residue = >=2 atom contacts (DESIGN §0 primary; census `*_biolip2` columns). Bound = >=2 binding residues
  (`is_bound_biolip2`), MBS = `is_mbs_biolip2`. The census' >=1-pair rule is a reported sensitivity ("marsh_ge1pair").
- Primary set = not duplicate copy & bound & not glycan & not additive (`is_additive`: BioLiP2 label, FALLBACK for
  entries newer than BioLiP2) & no contact with an antibody chain (DESIGN §0: antibody chains excluded).
- Oligomers only (assembly_type homomer/heteromer), cluster30 present.
- Symmetric copies collapsed by (pdb_id, sym_group); site MBS = any copy MBS.
- Unit = (cluster30 x ccd_key). Representative entry = best resolution, ties earliest pdb_id. Unit MBS = any
  (collapsed) site of that ligand in the representative entry is MBS. `mbs_any_entry` = MBS in any entry (H1 sensitivity).
- Fine class: census ligand_class, with drug_like/other_organic re-assigned from ccd_druglike.parquet (`druglike_primary`).
- Metals (metal_ion, metal_cluster, inorganic) are never pooled with organic ligands in H1-H3 (DESIGN §0).
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/part1/analysis"
SEEDS = [0, 1, 2]                 # DESIGN §5; seed 0 primary
N_BOOT = 1000                     # DESIGN §5
POWER_GATE_MIN_CLUSTERS = 30      # DESIGN §1 power gate
ALPHA = 0.05                      # DESIGN §1
STRICT_MIN_RES_PER_CHAIN = 3      # DESIGN §0 strict MBS: >=3 residues on each of >=2 chains
C2_AXIS_MAX_DIST = 3.0            # DESIGN §1 H3 geometric sensitivity (A)
N_LEAVE_OUT = 10                  # DESIGN §2 leave-out top-10 clusters / ligands
SPLINE_DF = 4                     # our choice: B-spline df for heavy atoms in H1
REF_FRACTIONS = {"Q-BioLiP (Wei et al. 2024)": 0.288, "PLINDER (Durairaj et al. 2024)": 0.25}  # LITERATURE.md
METAL_CLASSES = {"metal_ion", "metal_cluster", "inorganic"}
OLIGO = ("homomer", "heteromer")

# Census column map. v2 (census_version column present): primary residue rule = >=2 atom pairs in the unsuffixed
# columns, >=1 pair in *_1atom, symmetric collapse by Jaccard `site_group`. v1: >=2 pairs in *_biolip2, `sym_group`.
COLMAP = {
    "v2": dict(bound="is_bound", mbs="is_mbs", cc="chain_counts", bound1="is_bound_1atom", mbs1="is_mbs_1atom",
               sym="site_group"),
    "v1": dict(bound="is_bound_biolip2", mbs="is_mbs_biolip2", cc="chain_counts_biolip2", bound1="is_bound",
               mbs1="is_mbs", sym="sym_group"),
}
BASE_COLS = ["in_primary_set", "pdb_id", "ligand_uid", "ccd", "ccd_key", "assembly_type", "n_protein_chains", "method_class",
             "resolution", "release_date", "is_dup_copy", "is_glycan", "is_additive", "plinder_artifact",
             "biolip_list_geom_artifact", "contacts_antibody", "cluster30", "marsh_class", "ligand_class",
             "n_heavy_atoms_instance", "n_heavy", "is_bound_4A", "is_mbs_4A", "is_on_axis", "on_axis_jaccard",
             "c2_axis_dist", "chainpair_rot_angle", "mbs_contrib_symmetry_copy", "biolip_status"]
OPTIONAL_COLS = ["census_version", "biolip_covered", "is_photosynthetic", "is_on_axis_geom", "additive_label_source"]
CENSUS_DIR = ROOT / "data/part1"
CENSUS = None   # "v1" | "v2", set in load()


# ------------------------------------------------------------------ provenance
def provenance(census_dir):
    sha = (ROOT / "results/DESIGN_FROZEN.sha256").read_text().split()[0]
    cur = hashlib.sha256((ROOT / "DESIGN.md").read_bytes()).hexdigest()
    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no", "scripts/part1_analysis.py", "scripts/build_druglike_flags.py", "src/part1/"],
                                capture_output=True, text=True).stdout.strip())
    dl = json.loads((ROOT / "data/processed/ccd_druglike.json").read_text())
    ls = census_dir / "ligand_sites.parquet"
    st = ls.stat()
    return dict(design_sha256_frozen=sha, design_sha256_current=cur, design_matches=sha == cur, git_commit=git,
                code_dirty=dirty, census_dir=str(census_dir.relative_to(ROOT)), census_version=CENSUS,
                ligand_sites_bytes=st.st_size,
                ligand_sites_mtime=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
                druglike_fallback=dl["fallback"], druglike_sentinels_pass=dl.get("sentinels_pass"),
                chembl_release=None if dl["fallback"] else dl["chembl_release_expected"],
                generated=time.strftime("%Y-%m-%dT%H:%M:%S%z"), status="UNVERIFIED_until_independent_re-derivation")


PROV = None


def write(df, name):
    p = OUT / f"{name}.tsv"
    with open(p, "w") as fh:
        fh.write("# " + " ".join(f"{k}={v}" for k, v in PROV.items()) + "\n")
        df.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    print("wrote", p, len(df), flush=True)


# ------------------------------------------------------------------ load
def load(census_dir):
    global CENSUS
    import pyarrow.parquet as pq
    names = set(pq.read_schema(census_dir / "ligand_sites.parquet").names)
    CENSUS = "v2" if "census_version" in names else "v1"
    cm = COLMAP[CENSUS]
    cols = [c for c in BASE_COLS + OPTIONAL_COLS if c in names] + [cm[k] for k in ("bound", "mbs", "cc", "bound1", "mbs1", "sym")]
    s = pd.read_parquet(census_dir / "ligand_sites.parquet", columns=list(dict.fromkeys(cols)))
    e = pd.read_parquet(census_dir / "entries.parquet", columns=["pdb_id", "is_structural_genomics", "pisa_agrees"])
    dl = pd.read_parquet(ROOT / "data/processed/ccd_druglike.parquet",
                         columns=["ccd", "druglike_primary", "druglike_ro5qed"])
    s = s.merge(e, on="pdb_id", how="left").merge(dl, on="ccd", how="left")
    # standardised names
    s["bound_p"], s["mbs_p"] = s[cm["bound"]], s[cm["mbs"]]
    s["bound_1atom"], s["mbs_1atom"] = s[cm["bound1"]], s[cm["mbs1"]]
    s["symkey"] = s[cm["sym"]].astype(str)
    if "biolip_covered" not in s:
        s["biolip_covered"] = s.biolip_status != "not_covered"
    if "is_on_axis_geom" not in s:   # v1: geometric rule from stored axis distance (C2 detected only if 180+-15 deg)
        s["is_on_axis_geom"] = s.c2_axis_dist.le(C2_AXIS_MAX_DIST)
    # fine class (primary drug-like = ChEMBL rule); secondary = Ro5+QED. Only drug_like/other_organic are re-assigned.
    # drug_like census label is discarded; drug-like flags may override the classes listed in ccd_druglike.json
    over = json.loads((ROOT / "data/processed/ccd_druglike.json").read_text())["druglike_overridable_classes"]
    base = s.ligand_class.fillna("other_organic").replace({"drug_like": "other_organic", "unknown_ccd": "other_organic"})
    reassign = base.isin(over)
    s["cls"] = np.where(reassign & s.druglike_primary.fillna(False).astype(bool), "drug_like", base)
    s["cls_ro5qed"] = np.where(reassign & s.druglike_ro5qed.fillna(False).astype(bool), "drug_like", base)
    s["year"] = pd.to_datetime(s.release_date, errors="coerce").dt.year
    s["sg"] = s.is_structural_genomics.fillna(False).astype(int)
    s["is_additive"] = s.is_additive.astype(bool)
    cc = s[cm["cc"]]
    s["is_mbs_strict"] = [sum(c["n_res"] >= STRICT_MIN_RES_PER_CHAIN for c in x) >= 2 if x is not None else False for x in cc]
    s["_chk_mbs_from_counts"] = [sum(c["n_res"] >= 1 for c in x) >= 2 if x is not None else False for x in cc]
    s = s.drop(columns=[cm["cc"]])
    s["n_heavy_lig"] = s.n_heavy_atoms_instance.where(s.n_heavy_atoms_instance > 0, s.n_heavy)
    return s


# ------------------------------------------------------------------ variants
PRIMARY = dict(bound="bound_p", mbs="mbs_p", additive="is_additive", cls="cls", pisa_only=False, covered=True,
               excl_photo=False, onaxis="is_on_axis")
VARIANTS = {
    "primary": {},
    "strict_mbs": dict(mbs="is_mbs_strict"),
    "contacts_4A": dict(bound="is_bound_4A", mbs="is_mbs_4A"),
    "plinder_additives": dict(additive="plinder_artifact"),
    "biolip_geom_additives": dict(additive="biolip_list_geom_artifact"),
    "assembly_author_software_agree": dict(pisa_only=True),
    "exclude_photosynthetic": dict(excl_photo=True),                       # PI request 2026-10-29 (audit M-1)
    "post_biolip_fallback_set": dict(covered=False),                      # PI request: entries newer than BioLiP2
    "druglike_ro5qed_secondary": dict(cls="cls_ro5qed"),                  # DESIGN §0 secondary definition
    "marsh_ge1pair_residue": dict(bound="bound_1atom", mbs="mbs_1atom"),  # EXPLORATORY (not in DESIGN §2)
}


def site_filter(s, v, ledger=None, tag="", apply_additive=True):
    """Return collapsed oligomer sites for variant v; append cascade to ledger."""
    steps = [("all ligand instances (assembly 1)", np.ones(len(s), bool)),
             ("not duplicate copy", ~s.is_dup_copy.values),
             (f"bound ({v['bound']})", s[v["bound"]].fillna(False).values.astype(bool)),
             ("not glycan", ~s.is_glycan.values)]
    if v["covered"] is True:
        steps.append(("entry covered by BioLiP2 snapshot", s.biolip_covered.values.astype(bool)))
    elif v["covered"] is False:
        steps.append(("entry NOT covered by BioLiP2 (FALLBACK additive label)", ~s.biolip_covered.values.astype(bool)))
    if apply_additive:
        steps.append((f"not additive ({v['additive']})", ~s[v["additive"]].fillna(False).astype(bool).values))
    steps += [("no antibody-chain contact", ~s.contacts_antibody.fillna(False).values.astype(bool)),
              ("oligomer assembly (homomer/heteromer)", s.assembly_type.isin(OLIGO).values),
              ("cluster30 assigned", s.cluster30.notna().values)]
    if v["pisa_only"]:
        steps.append(("assembly author/software agree", s.pisa_agrees.fillna(False).astype(bool).values))
    if v["excl_photo"]:
        steps.append(("not photosynthetic entry", ~s.is_photosynthetic.fillna(False).astype(bool).values))
    m = np.ones(len(s), bool)
    for name, mask in steps:
        m &= mask
        if ledger is not None:
            ledger.append(dict(variant=tag, level="ligand_instance", step=name, n=int(m.sum()),
                               n_clusters=int(s.cluster30[m].nunique())))
    x = s[m].copy()
    x["y_mbs"] = x[v["mbs"]].fillna(False).astype(bool)
    x["c"] = x[v["cls"]]
    # symmetric-copy collapse (site MBS / on-axis = any copy)
    x = x.sort_values("ligand_uid")
    g = x.groupby(["pdb_id", "symkey"], sort=False)
    col = g.first()
    col["y_mbs"] = g.y_mbs.max()
    col["onax"] = g.is_on_axis.max()
    col["onax_geom"] = g.is_on_axis_geom.max()
    col["c2d"] = g.c2_axis_dist.min()
    col = col.reset_index()
    if ledger is not None:
        ledger.append(dict(variant=tag, level="site (symmetric copies collapsed)", step=f"collapse by {COLMAP[CENSUS]['sym']}",
                           n=len(col), n_clusters=int(col.cluster30.nunique())))
    return col


def build_units(sites, ledger=None, tag=""):
    """Non-redundant units (cluster30 x ccd_key) from collapsed sites."""
    rep = (sites.sort_values(["resolution", "pdb_id"]).drop_duplicates(["cluster30", "ccd_key"])
           [["cluster30", "ccd_key", "pdb_id"]].rename(columns={"pdb_id": "rep_pdb"}))
    x = sites.merge(rep, on=["cluster30", "ccd_key"])
    anyent = x.groupby(["cluster30", "ccd_key"]).y_mbs.max().rename("mbs_any_entry")
    r = x[x.pdb_id == x.rep_pdb].copy()
    r["onax_mbs"] = r.onax.fillna(False).astype(bool) & r.y_mbs
    r["onaxg_mbs"] = r.onax_geom.fillna(False).astype(bool) & r.y_mbs
    r["c2_mbs"] = r.c2d.where(r.y_mbs)
    g = r.groupby(["cluster30", "ccd_key"])
    u = g.agg(pdb_id=("pdb_id", "first"), ccd=("ccd", "first"), cls=("c", "first"), marsh_class=("marsh_class", "first"),
              assembly_type=("assembly_type", "first"), y_mbs=("y_mbs", "max"), n_sites=("y_mbs", "size"),
              n_heavy=("n_heavy_lig", "max"), resolution=("resolution", "first"), method=("method_class", "first"),
              n_chains=("n_protein_chains", "first"), sg=("sg", "first"), year=("year", "first"),
              onax=("onax_mbs", "max"), onax_geom=("onaxg_mbs", "max"), c2d=("c2_mbs", "min"),
              symcopy=("mbs_contrib_symmetry_copy", "max")).reset_index()
    u = u.join(anyent, on=["cluster30", "ccd_key"])
    u["y_mbs"] = u.y_mbs.astype(int)
    u["mbs_any_entry"] = u.mbs_any_entry.astype(int)
    u["is_metal"] = u.cls.isin(METAL_CLASSES)
    if ledger is not None:
        ledger.append(dict(variant=tag, level="unit (cluster30 x CCD)", step="representative entry",
                           n=len(u), n_clusters=int(u.cluster30.nunique())))
    return u


# ------------------------------------------------------------------ bootstrap
def cluster_boot_frac(df, strata_cols, seed, ycol="y_mbs", n_boot=N_BOOT):
    """MBS fraction per stratum with percentile CIs from resampling cluster30 (with replacement)."""
    d = df.assign(_y=df[ycol].astype(float))
    clusters, cidx = np.unique(d.cluster30.values, return_inverse=True)
    d["_ci"] = cidx
    if strata_cols:
        key = d[strata_cols[0]].astype(str)
        for c in strata_cols[1:]:
            key = key + "|" + d[c].astype(str)
    else:
        key = pd.Series("all", index=d.index)
    skeys, sidx = np.unique(key.values, return_inverse=True)
    C, S = len(clusters), len(skeys)
    N = np.zeros((C, S)); K = np.zeros((C, S))
    np.add.at(N, (cidx, sidx), 1.0); np.add.at(K, (cidx, sidx), d._y.values)
    rng = np.random.default_rng(seed)
    fr = np.empty((n_boot, S))
    for b0 in range(0, n_boot, 100):
        nb = min(100, n_boot - b0)
        W = np.stack([np.bincount(rng.integers(0, C, C), minlength=C) for _ in range(nb)]).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            fr[b0:b0 + nb] = (W @ K) / (W @ N)
    est = K.sum(0) / N.sum(0)
    lo, hi = np.nanpercentile(fr, [2.5, 97.5], axis=0)
    ncl = (N > 0).sum(0)
    out = pd.DataFrame({"stratum": skeys, "n": N.sum(0).astype(int), "n_mbs": K.sum(0).astype(int),
                        "n_clusters": ncl, "mbs_frac": est, "ci_lo": lo, "ci_hi": hi, "seed": seed})
    if strata_cols:
        out[strata_cols] = out.stratum.str.split("|", expand=True)
    return out.drop(columns="stratum")


def boot_diff(df, ga, gb, seed, ycol="y_mbs", col="cls", n_boot=N_BOOT):
    """Cluster-bootstrap CI for frac(ga) - frac(gb)."""
    d = df[df[col].isin([ga, gb])]
    clusters, cidx = np.unique(d.cluster30.values, return_inverse=True)
    isa = (d[col] == ga).values
    y = d[ycol].astype(float).values
    C = len(clusters)
    Na, Ka = np.bincount(cidx, isa, C), np.bincount(cidx, isa * y, C)
    Nb, Kb = np.bincount(cidx, ~isa, C), np.bincount(cidx, (~isa) * y, C)
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, C, C), minlength=C)
        diffs.append((w @ Ka) / (w @ Na) - (w @ Kb) / (w @ Nb))
    est = Ka.sum() / Na.sum() - Kb.sum() / Nb.sum()
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    return est, lo, hi


# ------------------------------------------------------------------ hypotheses
def h1(u, ycol="y_mbs"):
    d = u[~u.is_metal].dropna(subset=["n_heavy", "resolution", "year", "n_chains"]).copy()
    d["y"] = d[ycol]
    gate = d.groupby("cls").cluster30.nunique()
    ok = all(gate.get(c, 0) >= POWER_GATE_MIN_CLUSTERS for c in ("cofactor", "drug_like"))
    res = dict(test="H1", n_units=len(d), n_clusters=int(d.cluster30.nunique()),
               gate_clusters_cofactor=int(gate.get("cofactor", 0)), gate_clusters_druglike=int(gate.get("drug_like", 0)),
               power_gate_pass=ok, mbs_frac_cofactor=d.y[d.cls == "cofactor"].mean(),
               mbs_frac_druglike=d.y[d.cls == "drug_like"].mean())
    if d.method.nunique() > 1:
        meth = " + C(method)"
    else:
        meth = ""
    f = (f"y ~ C(cls, Treatment('drug_like')) + bs(n_heavy, df={SPLINE_DF}) + resolution{meth} + n_chains + sg + year")
    try:
        m = smf.glm(f, d, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": d.cluster30.values})
        k = "C(cls, Treatment('drug_like'))[T.cofactor]"
        ci = m.conf_int().loc[k]
        res.update(estimand="OR cofactor vs drug-like (adjusted)", estimate=np.exp(m.params[k]), ci_lo=np.exp(ci[0]),
                   ci_hi=np.exp(ci[1]), p=m.pvalues[k], converged=bool(m.converged))
        coefs = pd.DataFrame({"term": m.params.index, "log_or": m.params.values, "se_cluster": m.bse.values,
                              "or": np.exp(m.params.values), "ci_lo": np.exp(m.conf_int()[0].values),
                              "ci_hi": np.exp(m.conf_int()[1].values), "p": m.pvalues.values})
    except Exception as ex:  # noqa: BLE001
        res.update(estimand="OR cofactor vs drug-like (adjusted)", estimate=np.nan, error=repr(ex))
        coefs = None
    return res, coefs


def clr_binary(d, xcol, ycol, gcol):
    """Conditional logistic regression with one binary predictor, exact: the conditional likelihood per stratum is
    Fisher's noncentral hypergeometric of the count (x=1 & y=1) given the stratum margins (identical to the CLR
    likelihood; statsmodels ConditionalLogit hit a recursion limit on large strata). Wald CI from the observed
    information (numerical), p from the likelihood-ratio test vs OR = 1."""
    from scipy.optimize import minimize_scalar
    t = d.assign(_xy=d[xcol] * d[ycol]).groupby(gcol).agg(M=(xcol, "size"), n=(xcol, "sum"), N=(ycol, "sum"),
                                                          a=("_xy", "sum"))
    t = t[(t.n > 0) & (t.n < t.M) & (t.N > 0) & (t.N < t.M)]
    M, n, N, a = (t[c].values.astype(int) for c in ("M", "n", "N", "a"))

    def nll(b):
        return -stats.nchypergeom_fisher.logpmf(a, M, n, N, np.exp(b)).sum()
    opt = minimize_scalar(nll, bounds=(-10, 10), method="bounded", options={"xatol": 1e-8})
    b, h = opt.x, 1e-3
    info = (nll(b + h) - 2 * nll(b) + nll(b - h)) / h ** 2
    se = 1 / np.sqrt(info)
    lr = 2 * (nll(0.0) - opt.fun)
    return dict(estimate=float(np.exp(b)), ci_lo=float(np.exp(b - 1.959964 * se)), ci_hi=float(np.exp(b + 1.959964 * se)),
                p=float(stats.chi2.sf(max(lr, 0), 1)), clr_strata=len(t))


def h2(u, ycol="y_mbs"):
    d = u[~u.is_metal].copy()
    d["dl"] = (d.cls == "drug_like").astype(int)
    d["x"] = d[ycol].astype(int)
    g = d.groupby("cluster30").agg(xmin=("x", "min"), xmax=("x", "max"), dmin=("dl", "min"), dmax=("dl", "max"))
    both = g[(g.xmin == 0) & (g.xmax == 1)].index
    disc = g[(g.xmin == 0) & (g.xmax == 1) & (g.dmin == 0) & (g.dmax == 1)].index
    dd = d[d.cluster30.isin(both)]
    res = dict(test="H2", n_units=len(dd), n_clusters=len(both), n_discordant_clusters=len(disc),
               power_gate_pass=len(disc) >= POWER_GATE_MIN_CLUSTERS,
               druglike_frac_mbs=dd.dl[dd.x == 1].mean(), druglike_frac_sbs=dd.dl[dd.x == 0].mean(),
               estimand="conditional OR drug-like occupancy, MBS vs SBS (within cluster30)")
    di = dd[dd.cluster30.isin(disc)]
    try:
        res.update(clr_binary(di, "x", "dl", "cluster30"))
    except Exception as ex:  # noqa: BLE001
        res.update(estimate=np.nan, error=repr(ex))
    # Mantel-Haenszel cross-check over the same strata
    tabs = []
    for _, gg in di.groupby("cluster30"):
        t = np.array([[((gg.x == 1) & (gg.dl == 1)).sum(), ((gg.x == 1) & (gg.dl == 0)).sum()],
                      [((gg.x == 0) & (gg.dl == 1)).sum(), ((gg.x == 0) & (gg.dl == 0)).sum()]])
        tabs.append(t)
    if tabs:
        st = StratifiedTable(tabs)
        res.update(mh_or=st.oddsratio_pooled, mh_ci_lo=st.oddsratio_pooled_confint()[0],
                   mh_ci_hi=st.oddsratio_pooled_confint()[1], mh_p=st.test_null_odds().pvalue)
    return res


def h3(u, axis="jaccard"):
    d = u[(u.assembly_type == "homomer") & (u.y_mbs == 1) & u.cls.isin(["drug_like", "cofactor"])].copy()
    if axis == "jaccard":
        d["on"] = d.onax.fillna(False).astype(int)
    elif axis == "union":  # EXPLORATORY (audit M-3): Jaccard OR geometric
        d["on"] = (d.onax.fillna(False).astype(bool) | d.onax_geom.fillna(False).astype(bool)).astype(int)
    else:  # geometric: ligand centroid <= C2_AXIS_MAX_DIST from the C2 axis (census is_on_axis_geom; no C2 -> off-axis)
        d["on"] = d.onax_geom.fillna(False).astype(int)
    d["dl"] = (d.cls == "drug_like").astype(int)
    tab = pd.crosstab(d.on, d.dl).reindex(index=[1, 0], columns=[1, 0], fill_value=0)
    ncl = d.groupby("on").cluster30.nunique()
    cells = d.groupby(["on", "dl"]).cluster30.nunique()
    ok = all(ncl.get(k, 0) >= POWER_GATE_MIN_CLUSTERS for k in (0, 1))
    res = dict(test={"jaccard": "H3", "geom": "H3_geom_c2", "union": "H3_union_jaccard_or_geom"}[axis], n_units=len(d), n_clusters=int(d.cluster30.nunique()),
               n_on_druglike=int(tab.loc[1, 1]), n_on_cofactor=int(tab.loc[1, 0]), n_off_druglike=int(tab.loc[0, 1]),
               n_off_cofactor=int(tab.loc[0, 0]), gate_clusters_on=int(ncl.get(1, 0)), gate_clusters_off=int(ncl.get(0, 0)),
               min_cell_clusters=int(min(cells.get((a, b), 0) for a in (0, 1) for b in (0, 1))),
               power_gate_pass=ok, estimand="OR drug-like (vs cofactor), on-axis vs off-axis homomer MBS")
    o, pf = stats.fisher_exact(tab.values)
    res.update(fisher_or=o, fisher_p=pf)
    try:
        m = smf.glm("dl ~ on", d, family=sm.families.Binomial()).fit(cov_type="cluster",
                                                                     cov_kwds={"groups": d.cluster30.values})
        ci = m.conf_int().loc["on"]
        res.update(estimate=np.exp(m.params["on"]), ci_lo=np.exp(ci[0]), ci_hi=np.exp(ci[1]), p=m.pvalues["on"])
    except Exception as ex:  # noqa: BLE001
        res.update(estimate=np.nan, error=repr(ex))
    return res


def run_tests(u, variant):
    rows = []
    r1, coefs = h1(u)
    rows.append(r1)
    if variant == "primary":
        r = h1(u, ycol="mbs_any_entry")[0]
        r["test"] = "H1"; r["sens_extra"] = "mbs_any_entry"
        rows.append(r)
    rows.append(h2(u))
    rows.append(h3(u, "jaccard"))
    rows.append(h3(u, "geom"))
    rows.append(h3(u, "union"))
    for r in rows:
        r["variant"] = variant + ("/" + r.pop("sens_extra") if "sens_extra" in r else "")
    return rows, coefs


# ------------------------------------------------------------------ main
def main():
    global PROV
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--census-dir", default=str(CENSUS_DIR))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    globals()["OUT"] = Path(a.out).resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    cdir = Path(a.census_dir).resolve()
    t0 = time.time()
    s = load(cdir)
    PROV = provenance(cdir)
    (OUT / "provenance.json").write_text(json.dumps(PROV, indent=1))
    print(f"loaded {len(s)} rows in {time.time() - t0:.0f}s")
    chk = pd.DataFrame([dict(check="primary MBS column == (>=2 chains with >=1 residue in primary chain_counts)",
                             n_rows=len(s), n_disagree=int((s.mbs_p.fillna(False).astype(bool) != s._chk_mbs_from_counts).sum()))])
    if "in_primary_set" in s:   # our primary filter (before antibody/oligomer/cluster steps) vs census in_primary_set
        mine = (~s.is_dup_copy & s.bound_p.fillna(False).astype(bool) & ~s.is_glycan & s.biolip_covered.astype(bool)
                & ~s.is_additive & s.assembly_type.isin(OLIGO))
        cen = s.in_primary_set.fillna(False).astype(bool) & s.assembly_type.isin(OLIGO)
        chk = pd.concat([chk, pd.DataFrame([dict(check="oligomer primary filter == census in_primary_set (oligomers)",
                                                 n_rows=int(cen.sum()), n_disagree=int((mine != cen).sum()))])])
    write(chk, "consistency_checks")

    ledger, test_rows, desc = [], [], []
    units_by_variant = {}
    for name, ov in VARIANTS.items():
        v = {**PRIMARY, **ov}
        if v["excl_photo"] and "is_photosynthetic" not in s:
            print(f"skip {name}: is_photosynthetic not in census {CENSUS}")
            continue
        sites = site_filter(s, v, ledger, name)
        u = build_units(sites, ledger, name)
        units_by_variant[name] = u
        rows, coefs = run_tests(u, name)
        test_rows += rows
        if name == "primary":
            write(coefs, "h1_primary_coefficients")
            prim_sites = sites
        print(f"{name}: {len(u)} units, {time.time() - t0:.0f}s")

    u = units_by_variant["primary"]
    # leave-out top-10 clusters / ligands (by number of primary units)
    top_cl = u.cluster30.value_counts().index[:N_LEAVE_OUT]
    top_lig = u.ccd_key.value_counts().index[:N_LEAVE_OUT]
    for name, uu in [("leave_out_top10_clusters", u[~u.cluster30.isin(top_cl)]),
                     ("leave_out_top10_ligands", u[~u.ccd_key.isin(top_lig)])]:
        rows, _ = run_tests(uu, name)
        test_rows += [r for r in rows if r["variant"] == name]
        ledger.append(dict(variant=name, level="unit (cluster30 x CCD)", step="leave-out", n=len(uu),
                           n_clusters=int(uu.cluster30.nunique())))
    write(pd.DataFrame({"left_out_cluster30": top_cl.astype(int),
                        "n_units": u.cluster30.value_counts().values[:N_LEAVE_OUT]}), "leave_out_top10_clusters")
    write(pd.DataFrame({"left_out_ccd_key": top_lig, "n_units": u.ccd_key.value_counts().values[:N_LEAVE_OUT]}),
          "leave_out_top10_ligands")

    tests = pd.DataFrame(test_rows)
    # BH family A over primary H1-H3 that pass the power gate (DESIGN §1)
    prim = tests[(tests.variant == "primary") & tests.test.isin(["H1", "H2", "H3"]) & tests.power_gate_pass]
    tests["p_bh_familyA"] = np.nan
    if len(prim):
        tests.loc[prim.index, "p_bh_familyA"] = multipletests(prim.p.values, method="fdr_bh")[1]
    tests["role"] = np.where((tests.variant == "primary") & tests.test.isin(["H1", "H2", "H3"]), "primary",
                             np.where(tests.variant.str.startswith("marsh_ge1pair") | tests.test.str.contains("union"),
                                      "EXPLORATORY", "sensitivity"))
    tests["inference"] = np.where(tests.power_gate_pass, "inferential", "DESCRIPTIVE (power gate failed)")
    write(tests, "hypothesis_tests")
    write(pd.DataFrame(ledger), "analysis_ledger")

    # descriptives ---------------------------------------------------
    for seed in SEEDS:
        tabs = []
        for filt in ("before_additive_filter", "after_additive_filter"):
            ss = prim_sites if filt == "after_additive_filter" else site_filter(s, PRIMARY, None, "", apply_additive=False)
            for lvl, dd in (("site_symcollapsed", ss), ("unit_nonredundant", build_units(ss))):
                for strata in (["assembly_type", "marsh_class"], ["assembly_type", "cls"], ["assembly_type"], []):
                    t = cluster_boot_frac(dd, strata, seed, n_boot=a.n_boot)
                    t["level"], t["filter"], t["grouping"] = lvl, filt, "+".join(strata) or "all_oligomers"
                    tabs.append(t)
                # organic-only overall (metals separate)
                t = cluster_boot_frac(dd[~dd.cls.isin(METAL_CLASSES)],
                                      ["assembly_type"], seed, n_boot=a.n_boot)
                t["level"], t["filter"], t["grouping"] = lvl, filt, "assembly_type|organic_only"
                tabs.append(t)
        desc_t = pd.concat(tabs, ignore_index=True)
        cols = ["seed", "level", "filter", "grouping", "assembly_type", "marsh_class", "cls", "n", "n_mbs",
                "n_clusters", "mbs_frac", "ci_lo", "ci_hi"]
        desc_t = desc_t.reindex(columns=cols)
        write(desc_t, f"mbs_fraction_descriptives_seed{seed}")
        if seed == 0:
            desc0 = desc_t
    # cofactor - drug-like MBS fraction difference (unit level, organic), cluster bootstrap, all seeds
    diffs = []
    for seed in SEEDS:
        for at in ("all", "homomer", "heteromer"):
            uu = u if at == "all" else u[u.assembly_type == at]
            est, lo, hi = boot_diff(uu, "cofactor", "drug_like", seed, n_boot=a.n_boot)
            diffs.append(dict(seed=seed, assembly_type=at, contrast="MBS frac cofactor - drug_like (units)",
                              estimate=est, ci_lo=lo, ci_hi=hi))
    write(pd.DataFrame(diffs), "h1_effect_size_bootstrap")
    # comparison to published fractions
    cmp_rows = []
    for (lvl, filt, grp), g in desc0.groupby(["level", "filter", "grouping"]):
        if grp != "all_oligomers":
            continue
        r = g.iloc[0]
        for ref, val in REF_FRACTIONS.items():
            cmp_rows.append(dict(level=lvl, filter=filt, ours=r.mbs_frac, ci_lo=r.ci_lo, ci_hi=r.ci_hi, reference=ref,
                                 ref_value=val, ref_in_ci=bool(r.ci_lo <= val <= r.ci_hi)))
    write(pd.DataFrame(cmp_rows), "comparison_published_fractions")
    (ROOT / "data/part1/analysis").mkdir(parents=True, exist_ok=True)   # data/ is git-ignored
    u.to_parquet(ROOT / f"data/part1/analysis/units_primary_{CENSUS}.parquet", index=False)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
