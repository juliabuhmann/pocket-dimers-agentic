#!/usr/bin/env python3
"""Metrics + paired statistics for the EXPLORATORY heterodimer MBS extension (post-hoc; docs/COFOLD_HETERO_EXT.md).

Status: first attempt, UNVERIFIED.

Rerun:  source scripts/env_threads.sh; .venvs/part1/bin/python scripts/cofold_metrics_hetero_ext.py
        -> results/cofold/hetero_ext/{metrics,summary,stats}.tsv

Scores the top-ranked seed-0 sample of every finished run of
  * the extension (results/cofold/hetero_ext/cases.tsv, AF3 outputs data/cofold/hetero_ext/af3) and
  * the 6 original heterodimer cases (results/cofold/h4/cases.tsv, stratum hetero; data/cofold/h4/af3),
with two independent code paths (both imported unmodified):
  * verifier (PRIMARY): reviews/verify_cofold/verify_rmsd.py score_job -> symmetry-aware ligand RMSD (CCD graph
    automorphisms, min over chain assignments / ligand copies), clash flag; site_share -> partner share
    (BioLiP-like binding residues on the minor chain / all binding residues of the pair, PDB reference copy);
  * producer: scripts/cofold_metrics_h4.py evaluate -> lig_rmsd (symmetry-corrected), dcc, contacts_both,
    lig_prot_pair_iptm (ligand ipTM = max ligand-protein chain_pair_iptm), ranking_score.
Success = RMSD <= SUCCESS_RMSD. Paired A (heterodimer) vs B (single chain with most contacts) statistics:
exact McNemar (two-sided binomial on discordant pairs), sign test on per-case RMSD (dimer < single), Spearman of
partner share vs partner gain (log10 RMSD_single / RMSD_dimer). Splits: origin (original / extension), ligand tier,
release vs AF3 cutoff, pair complex released before the cutoff.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "reviews/verify_cofold"))
import cofold_metrics_h4 as prod   # noqa: E402  producer metric
import verify_rmsd as ver          # noqa: E402  verifier metric

OUT = ROOT / "results/cofold/hetero_ext"
SETS = {   # origin -> (cases table, AF3 output dir, row filter)
    "original": (ROOT / "results/cofold/h4/cases.tsv", ROOT / "data/cofold/h4/af3", lambda d: d.stratum == "hetero"),
    "extension": (OUT / "cases.tsv", ROOT / "data/cofold/hetero_ext/af3", lambda d: d.stratum == "hetero_ext"),
}
TIER1 = ("drug_like", "cofactor")   # ligand classes of the original deep dive
SUCCESS_RMSD = prod.SUCCESS_RMSD    # 2 A
SEED = 0


def score(origin, case, af3dir, cond):
    job = f"{case.case_id}_{cond}"
    jd = af3dir / job
    rec = dict(origin=origin, case_id=case.case_id, cond=cond)
    if not (jd / f"{job}_ranking_scores.csv").exists():
        rec["status"] = "missing"
        return rec
    v = ver.score_job(("h4", case, cond, SEED, jd, job))   # batch "h4" -> chains chainA_pdb/chainB_pdb
    rec.update(status=v.get("status"), rmsd_ver=v.get("lig_rmsd"), ca_rmsd_ver=v.get("ca_rmsd"),
               clash_ver=v.get("clash"), min_lig_prot_dist=v.get("min_lig_prot_dist"), model=v.get("model"))
    if rec["status"] != "ok":
        return rec
    pre = Path(v["model"]).name.replace("_model.cif", "")
    sd = ROOT / Path(v["model"]).parent
    conf = json.loads((sd / f"{pre}_summary_confidences.json").read_text())
    p = prod.evaluate(case, cond, ROOT / v["model"], conf)
    rec.update(rmsd_prod=p.get("lig_rmsd"), dcc=p.get("dcc"), contacts_both=p.get("contacts_both"),
               n_res_A=p.get("n_res_A"), n_res_B=p.get("n_res_B"), lig_iptm=p.get("lig_prot_pair_iptm"),
               lig_chain_iptm=p.get("lig_chain_iptm"), ranking_score=p.get("ranking_score"), iptm=p.get("iptm"))
    return rec


def mcnemar_exact(b, c):
    n = b + c
    return 1.0 if n == 0 else float(min(1.0, 2 * stats.binom.cdf(min(b, c), n, 0.5)))


def paired_stats(w, label):
    d = w.dropna(subset=["rmsd_ver_dimer", "rmsd_ver_monomer"])
    sd, sm = d.rmsd_ver_dimer <= SUCCESS_RMSD, d.rmsd_ver_monomer <= SUCCESS_RMSD
    b, c = int((sd & ~sm).sum()), int((~sd & sm).sum())
    better, worse = int((d.rmsd_ver_dimer < d.rmsd_ver_monomer).sum()), int((d.rmsd_ver_dimer > d.rmsd_ver_monomer).sum())
    pd_, pm = d.rmsd_prod_dimer <= SUCCESS_RMSD, d.rmsd_prod_monomer <= SUCCESS_RMSD
    gain = np.log10(d.rmsd_ver_monomer / d.rmsd_ver_dimer)
    rho, prho = stats.spearmanr(d.partner_share, gain) if len(d) >= 4 else (np.nan, np.nan)
    return dict(subset=label, n=len(d), dimer_success=int(sd.sum()), single_success=int(sm.sum()),
                both=int((sd & sm).sum()), dimer_only=b, single_only=c, neither=int((~sd & ~sm).sum()),
                mcnemar_exact_p=mcnemar_exact(b, c),
                dimer_success_prod=int(pd_.sum()), single_success_prod=int(pm.sum()),
                dimer_only_prod=int((pd_ & ~pm).sum()), single_only_prod=int((~pd_ & pm).sum()),
                mcnemar_exact_p_prod=mcnemar_exact(int((pd_ & ~pm).sum()), int((~pd_ & pm).sum())),
                rmsd_dimer_lower=better, rmsd_dimer_higher=worse,
                sign_test_p=float(stats.binomtest(better, better + worse, 0.5).pvalue) if better + worse else 1.0,
                median_rmsd_dimer=d.rmsd_ver_dimer.median(), median_rmsd_single=d.rmsd_ver_monomer.median(),
                median_partner_share=d.partner_share.median(),
                spearman_share_vs_gain=rho, spearman_p=prho,
                dimer_contacts_both=int(d.contacts_both_dimer.fillna(False).astype(bool).sum()),
                median_lig_iptm_dimer=d.lig_iptm_dimer.median(), median_lig_iptm_single=d.lig_iptm_monomer.median())


def main():
    ccds, rows, shares, meta = set(), [], [], []
    for origin, (ct, af3, filt) in SETS.items():
        cases = pd.read_csv(ct, sep="\t")
        cases = cases[filt(cases)]
        for c in cases.ccd.unique():
            ver.ccd_graph(c)
        for _, c in cases.iterrows():
            for cond in ("dimer", "monomer"):
                rows.append(score(origin, c, af3, cond))
            s = ver.site_share(("h4", c))
            shares.append(dict(case_id=c.case_id, partner_share=s["partner_share"], n_site=s["n_site"],
                               n_minor=s["n_minor"], site_chain_counts=s["site_chain_counts"]))
            meta.append(dict(origin=origin, case_id=c.case_id, pdb_id=c.pdb_id, ccd=c.ccd, ligand_class=c.ligand_class,
                             tier=1 if c.ligand_class in TIER1 else 2, release_date=c.release_date,
                             post_cutoff=bool(c.post_cutoff), pair_complex_pre_cutoff=not bool(c.post_cutoff_pair),
                             n_lig_dimer_run=c.n_lig_dimer_run, tokens_dimer=c.tokens_dimer))
            print(origin, c.case_id, flush=True)
    m = pd.DataFrame(rows)
    meta = pd.DataFrame(meta).merge(pd.DataFrame(shares), on="case_id")
    long = meta.merge(m, on=["origin", "case_id"], how="left")
    long["success_ver"] = long.rmsd_ver <= SUCCESS_RMSD
    long["success_prod"] = long.rmsd_prod <= SUCCESS_RMSD
    long.to_csv(OUT / "metrics.tsv", sep="\t", index=False, float_format="%.3f")

    ok = m[m.status == "ok"]
    vals = ["rmsd_ver", "rmsd_prod", "dcc", "contacts_both", "lig_iptm", "ranking_score", "clash_ver"]
    w = ok.pivot(index="case_id", columns="cond", values=vals)
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    for col in w.columns:
        if not col.startswith(("contacts_both", "clash")):
            w[col] = pd.to_numeric(w[col])
    w = meta.merge(w.reset_index(), on="case_id", how="inner")
    w["partner_gain_log10"] = np.log10(w.rmsd_ver_monomer / w.rmsd_ver_dimer)
    w.to_csv(OUT / "paired.tsv", sep="\t", index=False, float_format="%.3f")

    subsets = [("all", w), ("original_6", w[w.origin == "original"]), ("extension", w[w.origin == "extension"]),
               ("tier1_ligands", w[w.tier == 1]), ("tier2_ligands", w[w.tier == 2]),
               ("release_post_cutoff", w[w.post_cutoff]), ("release_pre_cutoff", w[~w.post_cutoff]),
               ("pair_complex_pre_cutoff", w[w.pair_complex_pre_cutoff]),
               ("pair_new_post_cutoff", w[~w.pair_complex_pre_cutoff]),
               ("extension_post_cutoff", w[(w.origin == "extension") & w.post_cutoff])]
    st = pd.DataFrame([paired_stats(g, lab) for lab, g in subsets if len(g)])
    st.to_csv(OUT / "summary.tsv", sep="\t", index=False, float_format="%.4g")
    agree = ok.dropna(subset=["rmsd_ver", "rmsd_prod"])
    print("verifier vs producer success agreement:",
          int(((agree.rmsd_ver <= SUCCESS_RMSD) == (agree.rmsd_prod <= SUCCESS_RMSD)).sum()), "/", len(agree))
    print(st.T.to_string())


if __name__ == "__main__":
    main()
