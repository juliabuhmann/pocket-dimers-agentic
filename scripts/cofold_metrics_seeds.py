#!/usr/bin/env python3
"""Per-seed metrics and success tables for the replicate seeds of the on-axis and dimer-only co-folding experiments.
Status: first attempt, UNVERIFIED.

Metric code is the existing code: dimer-only cases -> cofold_metrics_dimeronly.evaluate (symmetry-corrected spyrmsd
ligand RMSD, success = lig_rmsd_sym <= 2 A after global CA superposition). On-axis cases -> cofold_metrics_onaxis.evaluate
(name-matched RMSD `lig_rmsd`, the ORIGINAL success) plus the SAME symmetry-corrected RMSD as the dimer-only
experiment (DM.sym_rmsd, min over chain mappings and predicted ligand copies, CCD bond graphs), reported as
lig_rmsd_sym / success_sym. Success is judged on the top-ranked sample of each seed's 5 samples.

Outputs (results/cofold/<exp>/): seeds_metrics.tsv (top-ranked sample per case x condition x seed),
seeds_summary.tsv (table = per_seed | pooled | per_case). Per-sample metrics cached in data/cofold/seeds/metric_cache.
usage: cofold_metrics_seeds.py [onaxis] [dimeronly]"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cofold_metrics_dimeronly as DM  # noqa: E402
import cofold_metrics_onaxis as OM  # noqa: E402

SEEDS = [0, 1, 2]                      # project's fixed seed list
CONDS = ("dimer", "monomer")
SUCCESS = 2.0                          # A, same threshold as the original analyses
CACHE = ROOT / "data/cofold/seeds/metric_cache"
NEW_AF3 = ROOT / "data/cofold/seeds/af3"
EXP = {
    "onaxis": dict(cases=ROOT / "results/cofold/onaxis/cases.tsv", strata={"on_axis_druglike", "offaxis_cofactor"},
                   old=ROOT / "data/cofold/onaxis/af3", group="stratum"),
    "dimeronly": dict(cases=ROOT / "results/cofold/dimeronly/cases.tsv", strata=None,
                      old=ROOT / "data/cofold/dimeronly/af3", group=None),
}


def run_dir(exp, case_id, cond, seed):
    """AF3 output dir of one (case, condition, seed): original seed 0, old on-axis seed-1 dimers, new seeds."""
    old = EXP[exp]["old"]
    cands = [old / f"{case_id}_{cond}"] if seed == 0 else [old / f"{case_id}_{cond}_seed{seed}", NEW_AF3 / f"{case_id}_{cond}_seed{seed}"]
    for d in cands:
        if (d / f"{d.name}_ranking_scores.csv").exists():
            return d
    return None


def onaxis_sym(case, cond, model_path):
    """symmetry-corrected ligand RMSD for an on-axis case: same mapping/copy minimisation as OM.evaluate, DM.sym_rmsd."""
    rchains, rprot, rligs = OM.ref_structure(case)
    pprot, pligs, finite = OM.pred_structure(model_path, case.ccd)
    if not finite or not pligs or not rligs:
        return np.nan
    pch = sorted(pprot)
    maps = [dict(zip(pch, rchains)), dict(zip(pch, rchains[::-1]))] if cond == "dimer" else [{pch[0]: rc} for rc in rchains]
    best = np.inf
    for m in maps:
        X = np.vstack([OM.ca_pairs(pprot[p], rprot[r])[0] for p, r in m.items()])
        Y = np.vstack([OM.ca_pairs(pprot[p], rprot[r])[1] for p, r in m.items()])
        R, t = OM.kabsch(X, Y)
        for pl in pligs:
            best = min(best, DM.sym_rmsd(DM.tr(pl, R, t), rligs[0], case.ccd))
    return float(best)


def samples(exp, c, cond, seed, ctx):
    d = run_dir(exp, c.case_id, cond, seed)
    if d is None:
        return None
    cf = CACHE / f"{exp}_{d.name}.json"
    if cf.exists():
        return json.loads(cf.read_text())
    rows = []
    r = pd.read_csv(d / f"{d.name}_ranking_scores.csv").sort_values("ranking_score", ascending=False).reset_index(drop=True)
    for rank, x in r.iterrows():
        sd = d / f"seed-{int(x.seed)}_sample-{int(x['sample'])}"
        pre = f"{d.name}_seed-{int(x.seed)}_sample-{int(x['sample'])}"
        conf = json.loads((sd / f"{pre}_summary_confidences.json").read_text())
        try:
            if exp == "onaxis":
                m = OM.evaluate(c, cond, sd / f"{pre}_model.cif", conf)
                m["lig_rmsd_sym"] = onaxis_sym(c, cond, sd / f"{pre}_model.cif")
            else:
                m = DM.evaluate(c, cond, sd / f"{pre}_model.cif", conf, *ctx)
            st = "ok"
        except Exception as e:  # noqa: BLE001
            m, st = {}, f"error:{e}"
        rows.append(dict(rank=rank, sample=int(x["sample"]), status=st, **{k: (v.item() if hasattr(v, "item") else v) for k, v in m.items()}))
    CACHE.mkdir(parents=True, exist_ok=True)
    cf.write_text(json.dumps(rows))
    return rows


def main(exp):
    cfg = EXP[exp]
    cases = pd.read_csv(cfg["cases"], sep="\t")
    if cfg["strata"]:
        cases = cases[cases.stratum.isin(cfg["strata"])]
    chains = json.loads((ROOT / "data/cofold/dimeronly/case_chains.json").read_text()) if exp == "dimeronly" else None
    pockets = json.loads((ROOT / "data/cofold/dimeronly/case_pockets.json").read_text()) if exp == "dimeronly" else None
    top = []
    for _, c in cases.iterrows():
        ctx = (chains[c.case_id], pockets[c.case_id]) if exp == "dimeronly" else None
        for cond in CONDS:
            for seed in SEEDS:
                rows = samples(exp, c, cond, seed, ctx)
                if not rows:
                    continue
                ok = [r for r in rows if r["status"] == "ok" and r.get("lig_rmsd_sym", r.get("lig_rmsd")) is not None]
                if not ok:
                    continue
                t = next(r for r in ok if r["rank"] == 0)
                key = "lig_rmsd_sym"
                rec = dict(case_id=c.case_id, group=c[cfg["group"]] if cfg["group"] else "all", post_cutoff=c.post_cutoff,
                           cond=cond, seed=seed, n_samples=len(ok), lig_rmsd_sym=t[key], success_sym=bool(t[key] <= SUCCESS),
                           n_success5_sym=int(sum(r[key] <= SUCCESS for r in ok)), min5_sym=min(r[key] for r in ok),
                           ranking_score=t.get("ranking_score"), contacts_both=t.get("contacts_both"), ca_rmsd=t.get("ca_rmsd"))
                if exp == "onaxis":
                    rec.update(lig_rmsd_orig=t["lig_rmsd"], success_orig=bool(t["lig_rmsd"] <= SUCCESS),
                               lig_rmsd_hung=t.get("lig_rmsd_hung"))
                top.append(rec)
        print(exp, c.case_id, flush=True)
    T = pd.DataFrame(top)
    OUT = ROOT / f"results/cofold/{exp}"
    T.to_csv(OUT / "seeds_metrics.tsv", sep="\t", index=False, float_format="%.3f")
    defs = [("sym", "success_sym")] + ([("orig", "success_orig")] if exp == "onaxis" else [])
    rows = []
    for metric, col in defs:
        for (grp, cond), g in T.groupby(["group", "cond"]):
            for seed, gs in g.groupby("seed"):          # per seed
                rows.append(dict(table="per_seed", metric=metric, group=grp, cond=cond, seed=str(seed), n=len(gs),
                                 n_success=int(gs[col].sum()), rate=gs[col].mean()))
            # pooled over seeds, restricted to cases that have every seed that is present for >= 1 case of the group
            full = g.groupby("case_id").seed.nunique()
            ncomp = full.max()
            comp = g[g.case_id.isin(full[full == ncomp].index)]
            rows.append(dict(table="pooled", metric=metric, group=grp, cond=cond, seed=f"0-{ncomp - 1}_complete_cases",
                             n=len(comp), n_success=int(comp[col].sum()), rate=comp[col].mean()))
            rows.append(dict(table="pooled_all_runs", metric=metric, group=grp, cond=cond, seed="all_available",
                             n=len(g), n_success=int(g[col].sum()), rate=g[col].mean()))
            for cid, gc in g.groupby("case_id"):       # per-case consistency
                k = int(gc[col].sum())
                rows.append(dict(table="per_case", metric=metric, group=grp, cond=cond, seed=",".join(map(str, gc.seed)),
                                 case_id=cid, n=len(gc), n_success=k, rate=k / len(gc),
                                 consistent=bool(k in (0, len(gc))),
                                 lig_rmsd_by_seed=";".join(f"{v:.2f}" for v in gc[{"sym": "lig_rmsd_sym", "orig": "lig_rmsd_orig"}[metric]])))
    pd.DataFrame(rows).to_csv(OUT / "seeds_summary.tsv", sep="\t", index=False, float_format="%.3f")


if __name__ == "__main__":
    for e in (sys.argv[1:] or list(EXP)):
        main(e)
