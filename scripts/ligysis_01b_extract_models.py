"""EXPLORATORY. Extract AFDB v6 human models (pdb.gz) for every LIGYSIS human-site UniProt (all fragments in the atlas) from the proteome tar."""
import subprocess, os, pandas as pd
I = 'data/interim/ligysis/'; TAR = '/cluster/work/beltrao/shared/alphafold/v6/UP000005640_9606_HUMAN_v6.tar'
h = pd.read_pickle('data/raw/ligysis/LIGYSIS_human_sites_APRIL_2024.pkl')
s = pd.read_parquet(I + 'atlas_human_structures.parquet'); s = s[s.accession.isin(set(h.ACC))]
names = [f'AF-{a}-F{f}-model_v6.pdb.gz' for a, f in zip(s.accession, s.fragment)]
names = [n for n in names if not os.path.exists(I + 'af_models/' + n)]
open(I + 'needed_models_2.txt', 'w').write('\n'.join(names) + '\n'); print(len(names), 'to extract')
subprocess.run(['tar', 'xf', TAR, '-C', I + 'af_models', '-T', I + 'needed_models_2.txt'])  # a few names are absent from the tar (reported by tar)
