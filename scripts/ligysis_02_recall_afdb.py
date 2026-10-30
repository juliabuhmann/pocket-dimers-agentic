"""EXPLORATORY / UNVERIFIED. LIGYSIS-style top-(N+2) recall of atlas (AFDB monomer) BURROW pockets.

Site centre PROXY: LIGYSIS ligand centroids / site centres live in a PDBe superposition frame that is not available
locally, so the site centre is approximated by the centroid of the CA atoms of the site's UniProt residues, taken
directly in the AFDB model (no superposition). DCC = distance atlas-pocket centre -> that proxy centre.
Hit criteria (per site): DCC_proxy <= 12 A (LIGYSIS headline threshold), <= 4 A, residue overlap coefficient >= 0.5
(UniProt numbering). Pockets considered: top-(N+2) by burrow_score (N = LIGYSIS sites of the protein, ACC-level;
sensitivity: N = sites of the structural segment), and all atlas pockets ("max").
"""
import gzip, pickle, sys, numpy as np, pandas as pd, gemmi
from concurrent.futures import ProcessPoolExecutor
I = 'data/interim/ligysis/'; R = 'data/raw/ligysis/'; OUT = 'results/benchmarks/ligysis/'
WORKERS = 6; DCC_THR = (12.0, 4.0); OVL_THR = 0.5; EXTRA = 2

def load_model(acc, frag):
    """-> (offset, {local_resnum: (CA xyz, plddt)}) or None. UniProt = local + offset (from DBREF)."""
    try:
        txt = gzip.open(f'{I}af_models/AF-{acc}-F{frag}-model_v6.pdb.gz', 'rt').read().splitlines()
    except FileNotFoundError:
        return None
    off = 0
    ca = {}
    for l in txt:
        if l.startswith('DBREF'):
            w = l.split(); off = int(w[8]) - int(w[3])
        elif l.startswith('ATOM') and l[12:16] == ' CA ':
            ca[int(l[22:26])] = (np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66]))
    return off, ca

def work(args):
    acc, sites, pk = args   # sites: list of (site_name, seg, up_aas); pk: atlas pockets DataFrame for acc
    rows = []
    models = {}
    # fragment of each site residue: AF fragment k covers UniProt (k-1)*200+1 .. (k-1)*200+1400
    nfrag_needed = {}
    for name, seg, aas in sites:
        cnt = {}
        for a in aas:
            for k in range(max(1, (a - 1400) // 200 + 1), (a - 1) // 200 + 2):
                if (k - 1) * 200 + 1 <= a <= (k - 1) * 200 + 1400:
                    cnt[k] = cnt.get(k, 0) + 1
        nfrag_needed[name] = cnt
    for name, seg, aas in sites:
        cnt = nfrag_needed[name]
        best = max(cnt, key=lambda k: (cnt[k], -k)) if cnt else None
        centre = None; plddt = np.nan; n_in = 0
        if best is not None:
            if best not in models: models[best] = load_model(acc, best)
            m = models[best]
            if m is not None:
                off, ca = m
                pts = [ca[a - off] for a in aas if (a - off) in ca]
                n_in = len(pts)
                if n_in:
                    centre = np.mean([p[0] for p in pts], 0); plddt = np.mean([p[1] for p in pts])
        rows.append(dict(SITE_NAME=name, SEG_ID=f'{acc}_{seg}', best_frag=best, n_res_found=n_in, site_plddt=plddt,
                         centre=centre, aas=aas))
    return acc, rows

def main():
    h = pd.read_pickle(R + 'LIGYSIS_human_sites_APRIL_2024.pkl')
    pk = pd.read_parquet(I + 'atlas_human_pockets.parquet')
    st = pd.read_parquet(I + 'atlas_human_structures.parquet')
    jobs = []
    for acc, g in h.groupby('ACC'):
        if acc not in set(st.accession): continue
        jobs.append((acc, [(r.SITE_NAME, r.SEG, list(r.up_aas)) for r in g.itertuples()], pk[pk.accession == acc]))
    print('proteins with atlas structure row:', len(jobs), file=sys.stderr)
    with ProcessPoolExecutor(WORKERS) as ex:
        res = dict(ex.map(work, jobs, chunksize=20))
    # ---- evaluation
    pkg = {a: g for a, g in pk.groupby('accession')}
    out = []
    for acc, rows in res.items():
        g = pkg.get(acc)
        if g is not None:
            g = g.sort_values('burrow_score', ascending=False).reset_index(drop=True)
        n_acc = len(rows)
        nseg = pd.Series([r['SEG_ID'] for r in rows]).value_counts().to_dict()
        # UniProt-numbered residue sets per pocket
        pres = []
        # pocket residues in UniProt numbering: need fragment offset
        offs = {}
        if g is not None:
            for fr in g.fragment.unique():
                m = load_model(acc, int(fr)); offs[int(fr)] = m[0] if m else 0
            pres = [set(int(x.split(':')[1]) + offs[int(f)] for x in s.split()) for s, f in zip(g.residues, g.fragment)]
        for r in rows:
            for label, N in (('acc', n_acc), ('seg', nseg[r['SEG_ID']])):
                k = N + EXTRA
                d = dict(SITE_NAME=r['SITE_NAME'], SEG_ID=r['SEG_ID'], ACC=acc, level=label, N=N, n_pockets_atlas=0 if g is None else len(g),
                         n_site_res=len(r['aas']), n_res_found=r['n_res_found'], site_plddt=r['site_plddt'])
                if g is None or r['centre'] is None:
                    d.update(has_centre=r['centre'] is not None)
                    out.append(d); continue
                S = set(r['aas'])
                dcc = np.array([np.linalg.norm(np.array([p.center_x, p.center_y, p.center_z]) - r['centre']) if int(p.fragment) == r['best_frag'] else np.inf
                                for p in g.itertuples()])
                ovl = np.array([len(S & P) / min(len(S), len(P)) if P else 0 for P in pres])
                d['has_centre'] = True
                top = slice(0, k)
                for t in DCC_THR:
                    d[f'hit_dcc{int(t)}_topN2'] = bool((dcc[top] <= t).any()); d[f'hit_dcc{int(t)}_max'] = bool((dcc <= t).any())
                d['hit_ovl_topN2'] = bool((ovl[top] >= OVL_THR).any()); d['hit_ovl_max'] = bool((ovl >= OVL_THR).any())
                d['hit_either_topN2'] = d['hit_dcc4_topN2'] or d['hit_ovl_topN2']
                d['min_dcc_topN2'] = float(dcc[top].min()); d['max_ovl_topN2'] = float(ovl[top].max())
                out.append(d)
    df = pd.DataFrame(out)
    df.to_csv(OUT + 'ligysis_site_hits.tsv', sep='\t', index=False)
    print(df.groupby('level').size())

if __name__ == '__main__':
    main()
