"""EXPLORATORY / UNVERIFIED. Like-for-like check: BURROW (pocketcompare describe, merged pockets, NOT AF mode) on the
LIGYSIS representative EXPERIMENTAL chains (sample), top-(N+2) recall with the same site-centre PROXY as the AFDB run
(centroid of CA of site residues in the experimental chain; the LIGYSIS site-centre frame is not available locally).
Sample: SAMPLE_N rep chains (seed SEED) from LIGYSIS segments that have atlas pockets and whose DEF_TRANS site count
equals the human-site count. Chain is cut out (protein only) and residues renumbered to label_seq_id (= LIGYSIS 'aas').
Pockets filtered like the atlas: burrow_score >= SCORE_MIN (pocketcompare itself applies no score cut).
"""
import os, subprocess, sys, numpy as np, pandas as pd, gemmi
SAMPLE_N = 200; SEED = 0; SCORE_MIN = 0.03; EXTRA = 2; WORKERS = 6
PDB = '/path/to/pdb_mirror/'; W = 'data/interim/ligysis/exp/'; OUT = 'results/benchmarks/ligysis/'
os.makedirs(W + 'chains', exist_ok=True)

def main():
    t = pd.read_pickle('data/raw/ligysis/LIGYSIS_sites_DEF_TRANS.pkl')
    h = pd.read_pickle('data/raw/ligysis/LIGYSIS_human_sites_APRIL_2024.pkl')
    atl = set(pd.read_parquet('data/interim/ligysis/atlas_human_pockets.parquet').accession)
    nt = t.groupby('SEG_ID').size(); nh = h.groupby('SEG_ID').size()
    ok = [s for s in nt.index if nh.get(s) == nt[s]]
    reps = t[t.SEG_ID.isin(ok) & t.ACC.isin(atl)].drop_duplicates('SEG_ID')
    reps = reps.sample(SAMPLE_N, random_state=SEED)
    # cut chains
    keep = []
    for r in reps.itertuples():
        f = f'{W}chains/{r.rep_chain}.pdb'
        if not os.path.exists(f):
            st = gemmi.read_structure(f'{PDB}{r.pdb_id[1:3]}/{r.pdb_id}.cif.gz'); st.setup_entities()
            ch = st[0][r.auth_asym_id]
            ns = gemmi.Structure(); ns.cell = st.cell; ns.spacegroup_hm = st.spacegroup_hm
            m = gemmi.Model('1'); c = gemmi.Chain('A')
            for res in ch:
                if res.het_flag == 'H' or res.label_seq is None: continue
                r2 = res.clone(); r2.seqid = gemmi.SeqId(int(res.label_seq), ' '); c.add_residue(r2)
            m.add_chain(c); ns.add_model(m); ns.write_pdb(f)
        keep.append(f)
    open(W + 'list.txt', 'w').write('\n'.join(keep) + '\n')
    pq = W + 'pockets.parquet'
    if not os.path.exists(pq):
        subprocess.run(['pocketcompare', 'describe', W + 'chains', '--workers', str(WORKERS), '--backend', 'fast', '-o', pq],
                       check=True, stdout=open(W + 'describe.log', 'w'), stderr=subprocess.STDOUT)
    import pocketcompare.fpindex as fi
    meta, _, _ = fi.read_table(pq)
    meta = meta[meta.burrow_score >= SCORE_MIN]
    out = []
    for r in reps.itertuples():
        stem = r.rep_chain
        sites = t[t.SEG_ID == r.SEG_ID]
        ca = {res.seqid.num: np.array(a.pos.tolist()) for res in gemmi.read_structure(f'{W}chains/{stem}.pdb')[0]['A']
              for a in res if a.name == 'CA'}
        g = meta[meta.pocket_id.str.startswith(stem + ':')].sort_values('burrow_score', ascending=False)
        pres = [set(int(x.split(':')[1]) for x in s.split()) for s in g.residues]
        cen = g[['center_x', 'center_y', 'center_z']].to_numpy()
        N = len(sites)
        for s in sites.itertuples():
            aas = [a for a in s.aas if str(a).lstrip('-').isdigit()]  # drops the few insertion-code labels (e.g. '66A')
            pts = [ca[int(a)] for a in aas if int(a) in ca]
            d = dict(SITE_NAME=s.SITE_NAME, SEG_ID=s.SEG_ID, rep_chain=stem, N=N, n_pockets=len(g), n_res_found=len(pts))
            if pts and len(g):
                c = np.mean(pts, 0); dcc = np.linalg.norm(cen - c, axis=1)
                S = set(int(a) for a in aas); ovl = np.array([len(S & P) / min(len(S), len(P)) for P in pres])
                k = N + EXTRA
                for th in (12, 4): d[f'hit_dcc{th}_topN2'] = bool((dcc[:k] <= th).any()); d[f'hit_dcc{th}_max'] = bool((dcc <= th).any())
                d['hit_ovl_topN2'] = bool((ovl[:k] >= 0.5).any()); d['hit_ovl_max'] = bool((ovl >= 0.5).any())
                d['min_dcc_topN2'] = float(dcc[:k].min())
            out.append(d)
    df = pd.DataFrame(out); df.to_csv(OUT + 'ligysis_exp_site_hits.tsv', sep='\t', index=False)
    print(len(df), df.filter(like='hit').mean())

if __name__ == '__main__':
    main()
