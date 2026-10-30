"""Verifier: FE-logit on the rebuilt 300-protein subset (own categories + own humsavar y) vs the producer table on the same proteins."""
from pathlib import Path
import sys, numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from felogit import fe_logit, design, summarize
from statsmodels.discrete.conditional_models import ConditionalLogit
ROOT = str(Path(__file__).resolve().parents[2]); OUT = f"{ROOT}/reviews/verify_h5/out"
LEVELS = ["dimer_only", "sc_pocket", "interface", "surface"]
plp = pd.read_parquet(f"{OUT}/plp_in_models.parquet")
plp = plp[plp.model_aa == plp.wt][["accession", "position"]].drop_duplicates().assign(y_ver=1)
prod = pd.read_parquet(f"{ROOT}/data/part2/h5/tables/protein_residues_primary.parquet"); prod["cat"] = prod.cat.astype(str)
rows = []
for v in ["collapse", "nocollapse"]:
    S = pd.read_parquet(f"{OUT}/sample_protein_residues_{v}.parquet").merge(plp, on=["accession", "position"], how="left")
    S["y"] = S.y_ver.fillna(0).astype(int)
    if v == "collapse":
        chk = S.merge(prod[["accession", "position", "y"]], on=["accession", "position"], suffixes=("", "_prod"))
        print("y agreement with producer on sample:", (chk.y == chk.y_prod).mean(), "P/LP residues ver/prod", chk.y.sum(), chk.y_prod.sum())
    D = S[S.cat != "buried"]
    print(v, D.groupby("cat").y.agg(["sum", "count"]).to_dict())
    X = design(D, "interface", LEVELS); b, se, info = fe_logit(D.y.values, X.values, D.accession.values)
    rows += summarize(b, se, list(X.columns), info, f"rebuilt subset ({v}) FE-logit")
    if v == "collapse":
        inf = D.groupby("accession").y.transform(lambda s: 0 < s.sum() < len(s))
        Di = D[inf.astype(bool)]; Xi = design(Di, "interface", LEVELS)
        try:
            r = ConditionalLogit(Di.y.values, Xi.values, groups=Di.accession.values).fit(disp=0)
            for n, bb, ss in zip(Xi.columns, r.params, r.bse):
                rows.append(dict(analysis="rebuilt subset (collapse) exact clogit statsmodels", term=n, OR=np.exp(bb), lo=np.exp(bb - 1.96 * ss), hi=np.exp(bb + 1.96 * ss), p=float(r.pvalues[list(Xi.columns).index(n)]), n=len(Di), strata=Di.accession.nunique(), events=int(Di.y.sum())))
        except Exception as e: print("clogit failed", e)
samp = set(pd.read_csv(f"{OUT}/sample_accessions.tsv", sep="\t").accession)
Pp = prod[prod.accession.isin(samp) & (prod.cat != "buried")]
X = design(Pp, "interface", LEVELS); b, se, info = fe_logit(Pp.y.values, X.values, Pp.accession.values)
rows += summarize(b, se, list(X.columns), info, "producer table, same 300 proteins, FE-logit")
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/fits_subset.tsv", sep="\t", index=False); print(R[R.term.str.contains("_vs_")].to_string())
