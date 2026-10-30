"""Verifier fits for H5 on the producer's primary table (and the verifier-rebuilt subset). Not verified itself."""
from pathlib import Path
import sys, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
from felogit import fe_logit, design, summarize
ROOT = str(Path(__file__).resolve().parents[2]); OUT = f"{ROOT}/reviews/verify_h5/out"
LEVELS = ["dimer_only", "sc_pocket", "interface", "surface"]
N_BOOT, SEED, TOP_GENES = 200, 0, 20
AROM = {"PHE", "TYR", "TRP", "HIS"}

def fit(df, label, ref="interface", levels=LEVELS, cat_col="cat", aa=False, extra=("rsa_iso", "plddt"), bins=False):
    d = df.dropna(subset=list(extra)) if extra else df
    X = design(d, ref, levels, extra_num=() if bins else extra, aa=aa, cat_col=cat_col)
    if bins:  # flexible covariates: decile dummies of rsa_iso and pLDDT
        for c in extra:
            q = pd.qcut(d[c], 10, labels=False, duplicates="drop").values
            for b in range(1, q.max() + 1): X[f"{c}_q{b}"] = (q == b).astype(float)
    b, se, info = fe_logit(d.y.values, X.values, d.accession.values)
    return summarize(b, se, list(X.columns), info, label)

def boot_one(args):
    seed, = args
    rng = np.random.default_rng(seed)
    pr = np.array(BOOT_PROTS); pick = rng.choice(pr, len(pr), replace=True)
    parts = [BOOT_GROUPS[p].assign(accession=f"{p}#{i}") for i, p in enumerate(pick)]
    d = pd.concat(parts)
    X = design(d, "interface", LEVELS)
    b, se, info = fe_logit(d.y.values, X.values, d.accession.values)
    return b[list(X.columns).index("dimer_only_vs_interface")]

def mh(df, a, b, strata):
    d = df[df.cat.isin([a, b])]
    t = d.groupby(strata + ["cat"]).y.agg(["sum", "count"]).unstack("cat").fillna(0)
    A = t[("sum", a)]; B = t[("count", a)] - A; C = t[("sum", b)]; D = t[("count", b)] - C; N = A + B + C + D
    R = A * D / N; S = B * C / N
    orr = R.sum() / S.sum()
    # Robins-Breslow-Greenland variance
    P = (A + D) / N; Q = (B + C) / N
    v = (P * R).sum() / (2 * R.sum() ** 2) + ((P * S + Q * R).sum()) / (2 * R.sum() * S.sum()) + (Q * S).sum() / (2 * S.sum() ** 2)
    se = np.sqrt(v)
    return dict(analysis=f"MH {a} vs {b} | " + "x".join(strata), term=f"{a}_vs_{b}", OR=orr, lo=orr * np.exp(-1.96 * se), hi=orr * np.exp(1.96 * se),
                p=2 * stats.norm.sf(abs(np.log(orr) / se)), n=len(d), strata=int(((t[("sum", a)] + t[("sum", b)]) > 0).sum()), events=int(d.y.sum()))

def main():
    global BOOT_PROTS, BOOT_GROUPS
    P = pd.read_parquet(f"{ROOT}/data/part2/h5/tables/protein_residues_primary.parquet")
    for c in ["cat", "cat_split"]: P[c] = P[c].astype(str)
    D = P[P.cat != "buried"].copy()
    rows = []
    rows += fit(D, "A primary-equivalent FE-logit (cat + z rsa_iso + z pLDDT)")
    rows += fit(D, "B + amino-acid identity", aa=True)
    rows += fit(D, "C decile-binned rsa_iso & pLDDT + aa", aa=True, bins=True)
    rows += fit(D, "C2 FE-logit, no covariates", extra=())
    # gene burden: leave out top-20 genes by P/LP residue count (all residues incl. buried)
    burden = P.groupby("accession").y.sum().sort_values(ascending=False)
    top = list(burden.index[:TOP_GENES]); pd.Series(burden.head(TOP_GENES)).to_csv(f"{OUT}/top20_burden.tsv", sep="\t")
    rows += fit(D[~D.accession.isin(top)], f"D leave-out top-{TOP_GENES} proteins by P/LP")
    rows += fit(D[~D.accession.isin(top)], f"D2 leave-out top-{TOP_GENES} + aa", aa=True)
    # pLDDT restriction
    rows += fit(D[D.plddt >= 90], "E residues with pLDDT >= 90 only")
    # specificity
    rows += fit(D, "F ref = dimer_only (vs sc_pocket etc.)", ref="dimer_only")
    rows += fit(D, "F2 ref = dimer_only + aa", ref="dimer_only", aa=True)
    D["cs"] = D.cat_split
    rows += fit(D, "G cat_split ref = dimer_only", ref="dimer_only", levels=["dimer_only", "ifpocket_preexisting", "sc_pocket", "interface", "surface"], cat_col="cs")
    rows += fit(D, "G2 cat_split ref = interface", ref="interface", levels=["dimer_only", "ifpocket_preexisting", "sc_pocket", "interface", "surface"], cat_col="cs")
    # MH adjusted by strata
    D["rq"] = pd.qcut(D.rsa_iso, 5, labels=False, duplicates="drop"); D["pb"] = pd.cut(D.plddt, [0, 70, 80, 90, 101], labels=False)
    for a, b in [("dimer_only", "interface"), ("dimer_only", "sc_pocket"), ("sc_pocket", "interface")]:
        rows.append(mh(D, a, b, ["accession"])); rows.append(mh(D, a, b, ["accession", "rq", "pb"]))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/fits_full.tsv", sep="\t", index=False)
    print(R[R.term.str.contains("_vs_")].to_string())
    # composition / residue type
    comp = D.assign(G=D.resname == "GLY", Pr=D.resname == "PRO", C=D.resname == "CYS", Ar=D.resname.isin(AROM)).groupby("cat")[["G", "Pr", "C", "Ar", "y"]].mean()
    comp.to_csv(f"{OUT}/composition.tsv", sep="\t"); print(comp)
    aa_rate = D.groupby("resname").y.agg(["mean", "sum", "count"]); aa_rate.to_csv(f"{OUT}/plp_rate_by_aa.tsv", sep="\t")
    # burden of dimer-only P/LP
    do = D[(D.cat == "dimer_only") & (D.y == 1)].accession.value_counts()
    print("dimer-only P/LP residues", do.sum(), "proteins", len(do), "top5 share", do.head(5).sum() / do.sum(), do.head(10).to_dict())
    print("top20 burden proteins", burden.head(TOP_GENES).to_dict())
    do.to_csv(f"{OUT}/dimer_only_plp_by_protein.tsv", sep="\t")
    # protein bootstrap CI of primary-equivalent FE-logit
    inf = D.groupby("accession").y.agg(["sum", "count"]); inf = inf[(inf["sum"] > 0) & (inf["sum"] < inf["count"])].index
    Dz = D[D.accession.isin(inf)].copy()
    for c in ["rsa_iso", "plddt"]: Dz[c] = Dz[c]  # z-scoring happens inside design() per resample
    BOOT_PROTS = list(inf); BOOT_GROUPS = {k: g[["cat", "rsa_iso", "plddt", "y"]] for k, g in Dz.groupby("accession")}
    with ProcessPoolExecutor(6) as ex:
        bs = np.array(list(ex.map(boot_one, [(SEED * 10000 + i,) for i in range(N_BOOT)])))
    np.savetxt(f"{OUT}/boot_dimer_only.txt", bs)
    print("bootstrap OR CI", np.exp(np.percentile(bs, [2.5, 50, 97.5])))

if __name__ == "__main__":
    main()
