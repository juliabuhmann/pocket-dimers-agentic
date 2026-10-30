#!/usr/bin/env python
"""fpocket on the AF dimer models of the curated dimer-only list (EXPLORATORY; docs/FPOCKET_DRUGGABILITY.md).

Pockets compared (all BURROW pockets from data/part2/human_dimer_pockets_annotated.parquet, the Q2.2/Q2.3a source):
 - dimer-only: every row of results/part2/q23a_curated.tsv (curated dimer-only interface pockets);
 - single_chain_matched (primary control): for each dimer-only pocket, one single-chain pocket (all lining residues on
   one chain) of the SAME dimer model, drawn at random (SEED) without replacement among atlas-pass symmetry
   representatives (burrow_score >= 0.03, lining pLDDT >= 70, DESIGN.md §0);
 - single_chain_strict (sensitivity): the same draw restricted to pockets that also pass the curated list's
   pocket filters (volume_A3 in DRUGLIKE_VOL, burrow_score >= score_thr read from the curated TSV header). Few models
   have one (logged).
fpocket (default parameters) is run once on each dimer model PDB. Each BURROW pocket is mapped to the fpocket pocket
with the highest residue Jaccard |B & F| / |B | F| (ties: overlap coefficient, then better fpocket rank);
"mapped" = Jaccard >= MAP_MIN_JAC. Best-case sensitivity: max fpocket druggability over all fpocket pockets with
overlap coefficient |B & F| / min(|B|, |F|) >= BEST_MIN_OC (fp_drug_best).
Outputs: results/part1/fpocket/afdimer_pockets.tsv, afdimer_summary.tsv
Run: source scripts/env_threads.sh; .venvs/part1/bin/python scripts/fpocket_afdimer.py --workers 1
Status: first attempt, not yet verified.
"""
import argparse, re, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fpocket_common import WS, OUT, INTERIM, DRUG_THR, run_fpocket, fpocket_version
from q22_common import DIMER_ROOT, KINDS, ANNOT, DRUGLIKE_VOL, write_tsv

CURATED = WS / "results/part2/q23a_curated.tsv"
SEED = 0
MAP_MIN_JAC = 0.2    # our choice (residue Jaccard of the chosen fpocket pocket)
BEST_MIN_OC = 0.5    # our choice: half of the smaller residue set shared


def model_path(mid, kind):
    p = DIMER_ROOT / KINDS[kind] / f"{mid}.pdb"
    if not p.exists():
        raise FileNotFoundError(p)
    return p


def fp_model(args):
    mid, kind = args
    try:
        pk, _, res = run_fpocket(model_path(mid, kind), keep_residues=True)
        return mid, pk, res, ""
    except Exception as e:  # noqa: BLE001
        return mid, {}, {}, f"{type(e).__name__}: {e}"[:300]


def best_map(bres, res):
    best = (0.0, 0.0, None)
    for k in sorted(res):
        inter = len(bres & res[k])
        if not inter:
            continue
        oc = inter / min(len(bres), len(res[k])); jac = inter / len(bres | res[k])
        if (jac, oc) > best[:2]:
            best = (jac, oc, k)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    assert a.workers <= 10
    print(fpocket_version(), flush=True)
    hdr = open(CURATED).readline()
    thr = float(re.search(r"score_thr=([0-9.]+)", hdr).group(1))
    q = pd.read_csv(CURATED, sep="\t", comment="#")
    an = pd.read_parquet(ANNOT)
    an = an[an.modelEntityId.isin(q.modelEntityId)]
    single = an[(~an.is_interface) & ((an.n_res_A == 0) | (an.n_res_B == 0)) & an.atlas_pass & an.is_sym_rep]
    strict = single[single.volume_A3.between(*DRUGLIKE_VOL) & (single.burrow_score >= thr)]
    rows, log = [], []
    for grp, pool in (("single_chain_matched", single), ("single_chain_strict", strict)):
        rng = np.random.default_rng(SEED)
        for mid, g in q.sort_values("rank").groupby("modelEntityId", sort=False):
            cand = pool[pool.modelEntityId == mid].sort_values("pocket_id")
            pick = cand.iloc[rng.permutation(len(cand))[:len(g)]] if len(cand) else cand
            if len(pick) < len(g):
                log.append((grp, mid, len(g), len(pick)))
            for _, r in pick.iterrows():
                rows.append(dict(group=grp, pocket_id=r.pocket_id, modelEntityId=mid, kind=r.kind,
                                 burrow_score=r.burrow_score, volume_A3=r.volume_A3, residues=r.residues))
    for _, r in q.iterrows():
        rows.append(dict(group="dimer_only_interface", pocket_id=r.pocket_id, modelEntityId=r.modelEntityId, kind=r.kind,
                         burrow_score=r.burrow_score, volume_A3=r.volume_A3, residues=r.residues))
    P = pd.DataFrame(rows)
    print(P.group.value_counts().to_dict(), flush=True)
    print("models short of single-chain matches:", pd.Series([l[0] for l in log]).value_counts().to_dict(), flush=True)
    models = P[["modelEntityId", "kind"]].drop_duplicates().itertuples(index=False, name=None)
    from concurrent.futures import ProcessPoolExecutor
    cache = INTERIM / "afdimer_fpocket_raw.pkl"   # raw fpocket output per model (resumable / re-mappable)
    import pickle
    fp = pickle.load(open(cache, "rb")) if cache.exists() else {}
    todo = [m for m in models if m[0] not in fp]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        fp.update({mid: (pk, res, err) for mid, pk, res, err in ex.map(fp_model, todo)})
    pickle.dump(fp, open(cache, "wb"))
    out = []
    for _, r in P.iterrows():
        pk, res, err = fp[r.modelEntityId]
        rec = r.to_dict()
        rec.update(fp_error=err, fp_n_pockets=len(pk))
        if not err:
            B = set(r.residues.split())
            jac, oc, k = best_map(B, res)
            rec.update(map_oc=oc, map_jaccard=jac, mapped=jac >= MAP_MIN_JAC)
            if k is not None:
                rec.update(pk[k]); rec["fp_rank"] = k
            ds = [pk[j]["fp_drug"] for j in res if res[j] and len(B & res[j]) / min(len(B), len(res[j])) >= BEST_MIN_OC]
            rec["fp_drug_best"] = max(ds) if ds else np.nan
        out.append(rec)
    O = pd.DataFrame(out)
    O.to_parquet(INTERIM / "afdimer_pockets.parquet", index=False)
    write_tsv(O.drop(columns=["residues"]), OUT / "afdimer_pockets.tsv",
              note=f"fpocket on AF dimer models; {fpocket_version()}; map = max residue Jaccard, mapped if >= {MAP_MIN_JAC}; fp_drug_best = max over fpocket pockets with overlap coef >= {BEST_MIN_OC}")
    S = []
    for grp, g in O.groupby("group"):
        for sub, h, col in (("max-Jaccard map, any overlap", g[g.fp_drug.notna()], "fp_drug"),
                            (f"max-Jaccard map, Jaccard>={MAP_MIN_JAC}", g[g.mapped == True], "fp_drug"),  # noqa: E712
                            (f"best-case: max drug over pockets with oc>={BEST_MIN_OC}", g[g.fp_drug_best.notna()], "fp_drug_best")):
            S.append(dict(group=grp, subset=sub, n=len(h), n_models=h.modelEntityId.nunique(),
                          median_drug=h[col].median(), q25_drug=h[col].quantile(.25), q75_drug=h[col].quantile(.75),
                          frac_drug_gt_thr=(h[col] > DRUG_THR).mean(), median_fp_score=h.fp_score.median(),
                          median_fp_volume=h.fp_volume.median(), median_map_jaccard=h.map_jaccard.median(),
                          median_burrow_volume=h.volume_A3.median(), median_burrow_score=h.burrow_score.median()))
    S = pd.DataFrame(S)
    write_tsv(S, OUT / "afdimer_summary.tsv", note=f"druggability threshold {DRUG_THR}")
    print(S.to_string(), flush=True)
    # paired within model: dimer-only vs its matched single-chain pocket(s) (model-level means)
    pr = []
    for grp in ("single_chain_matched", "single_chain_strict"):
        a = O[O.group == "dimer_only_interface"].groupby("modelEntityId").fp_drug.mean()
        b = O[O.group == grp].groupby("modelEntityId").fp_drug.mean()
        j = pd.concat([a.rename("dimer_only"), b.rename("single")], axis=1).dropna()
        if len(j):
            from scipy.stats import wilcoxon
            p = wilcoxon(j.dimer_only, j.single).pvalue if (j.dimer_only != j.single).any() else np.nan
            pr.append(dict(control=grp, n_models=len(j), median_dimer_only=j.dimer_only.median(),
                           median_single=j.single.median(), frac_dimer_only_higher=(j.dimer_only > j.single).mean(),
                           frac_tied=(j.dimer_only == j.single).mean(), wilcoxon_p=p))
    write_tsv(pd.DataFrame(pr), OUT / "afdimer_paired.tsv", note="per-model means of max-Jaccard-mapped fp_drug")
    print(pd.DataFrame(pr).to_string(), flush=True)


if __name__ == "__main__":
    main()
