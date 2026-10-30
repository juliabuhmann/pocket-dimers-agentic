"""EXPLORATORY / UNVERIFIED. Recall tables (overall, multi-chain vs single-chain sites, ion-only, pLDDT) from ligysis_site_hits.tsv.
Multi-chain annotation: from LIGYSIS_human_chains_per_lig_MAY_2024.pkl (chains in contact with each ligand instance):
site_multichain_any = at least one supporting ligand contacts >1 chain; site_multichain_all = every ligand does.
Bootstrap CIs are over proteins (ACC clusters), B resamples, seed fixed."""
import pickle, numpy as np, pandas as pd
R = 'data/raw/ligysis/'; OUT = 'results/benchmarks/ligysis/'; B = 2000; SEED = 0
HITS = ['hit_dcc12_topN2', 'hit_dcc4_topN2', 'hit_ovl_topN2', 'hit_either_topN2', 'hit_dcc12_max', 'hit_ovl_max']

def ann():
    c = pickle.load(open(R + 'LIGYSIS_human_chains_per_lig_MAY_2024.pkl', 'rb'))
    ion = pickle.load(open(R + 'LIGYSIS_human_sites_ion_prop_APRIL_2024.pkl', 'rb'))
    rows = []
    for seg, v in c.items():
        for site, l in v.items():
            nch = [len(x) for x in l.values()]
            acc, s = seg.rsplit('_', 1)
            rows.append(dict(SITE_NAME=f'{seg}_{site}', n_lig=len(nch), site_multichain_any=max(nch) > 1, site_multichain_all=min(nch) > 1,
                             ion_prop=ion.get((acc, int(s), site), np.nan)))
    return pd.DataFrame(rows)

def boot(df, col, rng):
    g = df.groupby('ACC')[col].agg(['sum', 'count']); s, n = g['sum'].to_numpy(), g['count'].to_numpy()
    idx = rng.integers(0, len(g), (B, len(g)))
    r = s[idx].sum(1) / n[idx].sum(1)
    return np.percentile(r, [2.5, 97.5])

def summarise(df, label, rng):
    out = []
    for col in HITS:
        if col not in df: continue
        x = df.dropna(subset=[col]) if False else df
        v = x[col].fillna(False).astype(bool)
        d = x.assign(_h=v.astype(float)); lo, hi = boot(d, '_h', rng) if len(d) else (np.nan, np.nan)
        out.append(dict(subset=label, criterion=col, n_sites=len(d), n_proteins=d.ACC.nunique(), recall=v.mean(), ci_lo=lo, ci_hi=hi))
    return out

if __name__ == '__main__':
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(OUT + 'ligysis_site_hits.tsv', sep='\t').merge(ann(), on='SITE_NAME', how='left')
    allr = []
    for level in ('acc', 'seg'):
        d = df[df.level == level]
        # denominators: (A) every site of a protein with an AFDB structure row; (B) only proteins with >=1 atlas pocket; (C) B with a computable site centre
        dB = d[d.n_pockets_atlas > 0]; dC = dB[dB.has_centre.fillna(False)]
        subs = {'A_all_sites_with_AF_structure': d, 'B_proteins_with_atlas_pockets': dB, 'C_B_and_centre_found': dC}
        for nm, x in subs.items():
            allr += summarise(x, f'{level}|{nm}|all', rng)
        # splits on B (primary denominator)
        splits = {'single_chain_only(all ligands 1 chain)': ~dB.site_multichain_any, 'multichain_any': dB.site_multichain_any,
                  'multichain_all': dB.site_multichain_all, 'not_ion_only(ion_prop<1)': dB.ion_prop < 1,
                  'ion_only(ion_prop==1)': dB.ion_prop == 1, 'site_plddt>=70': dB.site_plddt >= 70, 'site_plddt<70': dB.site_plddt < 70,
                  'n_site_res<=5': dB.n_site_res <= 5, 'n_site_res>=16': dB.n_site_res >= 16}
        for nm, m in splits.items():
            allr += summarise(dB[m.fillna(False)], f'{level}|B|{nm}', rng)
    r = pd.DataFrame(allr); r.to_csv(OUT + 'ligysis_recall_summary.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
    print(r[r.criterion.isin(['hit_dcc12_topN2', 'hit_dcc4_topN2', 'hit_ovl_topN2', 'hit_dcc12_max'])].round(3).to_string())
