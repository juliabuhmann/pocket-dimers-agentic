"""EXPLORATORY. Extract human (taxid 9606) atlas pockets + structure rows from the pocketcompare AFDB snapshot."""
import glob, pandas as pd, pyarrow.dataset as ds, pyarrow.compute as pc
P = 'data/raw/pocketcompare_afdb_snapshot/'
OUT = 'data/interim/ligysis/'
TAXID = 9606
PK = "pocket_id,accession,fragment,burrow_rank,burrow_score,center_x,center_y,center_z,residues,plddt_mean,volume_A3,n_points".split(',')
ST = "accession,fragment,n_res,plddt_mean,n_pockets_all,n_pockets_kept,status".split(',')

def load(kind, cols):
    out = []
    for f in sorted(glob.glob(P + kind + '/tier*/*.parquet')):
        t = ds.dataset(f).to_table(columns=cols + ['taxid'], filter=pc.field('taxid') == TAXID)
        if t.num_rows: out.append(t.to_pandas())
    return pd.concat(out, ignore_index=True)

if __name__ == '__main__':
    d = load('pockets', PK); s = load('structures', ST)
    print(d.shape, d.accession.nunique(), s.shape, s.status.value_counts().to_dict())
    d.to_parquet(OUT + 'atlas_human_pockets.parquet'); s.to_parquet(OUT + 'atlas_human_structures.parquet')
    t = pd.read_pickle('data/raw/ligysis/LIGYSIS_sites_DEF_TRANS.pkl')
    h = pd.read_pickle('data/raw/ligysis/LIGYSIS_human_sites_APRIL_2024.pkl')
    for n, x in (('DEF_TRANS', t), ('human_sites', h)):
        a = set(x.ACC); print(n, len(a), 'with pockets', len(a & set(d.accession)), 'with structure row', len(a & set(s.accession)))
    print(s.fragment.value_counts().head(3).to_dict())
