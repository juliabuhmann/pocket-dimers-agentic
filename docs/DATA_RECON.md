# Data reconnaissance (2026-10-29)

All mounts read-only. Counts below were obtained by the commands shown; nothing was copied.
Scratch (`/cluster/scratch/...`) is purged after ~15 days: copy needed indices to project/work storage early.

## 1. Predicted dimers: `/cluster/scratch/gmagnusson/dimers`

Layout: `afdb_<species>_dimers/{homodimers,heterodimers}/` with `AF-<modelEntityId>.pdb` plus two CSVs per folder.
Species (13): arabidopsis, dicty, ecoli, fly, human, mouse, mtb, paeruginosa, pombe, rat, worm, yeast, zebrafish.
These are **AFDB quaternary-structure (complex) models** downloaded from `https://alphafold.ebi.ac.uk/files/AF-<id>-model_v1.pdb`
(header: "AlphaFold DB expands to proteome-scale quaternary structures", MODEL_V1, dated 21-FEB-26).
Only 2-chain models (chains `A`, `B`), already filtered (candidates have max_ipSAE >= 0.6 and max_pDockQ2 >= 0.23; observed minima 0.600 / 0.230).

Models per species (files = rows in candidates.csv):

| species | homo | hetero |
|---|---|---|
| arabidopsis | 3611 | 6440 |
| dicty | 878 | 1045 |
| ecoli | 1208 | 1449 |
| fly | 1621 | 2060 |
| **human** | **4821** (4819 in CSV) | **5793** (5791 in CSV) |
| mouse | 3223 | 4926 |
| mtb | 1072 | 655 |
| paeruginosa | 1619 | 1217 |
| pombe | 583 | 997 |
| rat | 3726 | 4495 |
| worm | 1717 | 1825 |
| yeast | 650 | 1252 |
| zebrafish | 3157 | 4224 |

(Counts from `find <dir> -type f | wc -l` include the 2 CSVs; the CSV row counts are authoritative: human 4819 homo, 5791 hetero; 4819 / 5754 unique accession pairs
[hetero: some pairs have >1 model; check `uniprotAccession_*` duplicates]. Species totals above for non-human are file counts incl. 2 CSVs; subtract 2.)

### Homo vs hetero
Directory name (`homodimers`/`heterodimers`) and `uniprotAccession_1 == uniprotAccession_2` in the CSV. Both are consistent.

### Chain -> UniProt mapping
**Do not use the PDB header.** TITLE/COMPND/SOURCE describe only the first protein ("MOL_ID", both chains A and B listed under one molecule,
TITLE `NAME (ACCESSION)`), even for heterodimers. Use `candidates.csv`:
`uniprotAccession_1` = chain A, `uniprotAccession_2` = chain B (verified by chain lengths on two heterodimers: O00203/Q92572 -> A 1094 + B 193 = n_residues 1287;
Q12962/Q01658 -> A 218, B 176). SEQRES block order may list chain B first; ATOM order is A then B. Residue numbers restart at 1 per chain
and are UniProt positions of the full-length (unfragmented) sequence. For chains > 999 residues, columns 23-26 overflow in the PDB format (e.g. `A1000`);
use gemmi/Biopython (not whitespace `split`) to parse.
No DBREF records and no CIF files.

### Confidence metrics
- Per residue pLDDT in the B-factor column (cols 61-66), both chains.
- Per model, only in the CSVs (not in PDB headers; no PAE JSON present): `ipTM`, `max_ipSAE`, `max_pDockQ2`, `n_residues` (total of both chains), `gene_1`, `gene_2`.
- Human ranges (min / median / max): homo ipTM 0.26/0.76/0.96, ipSAE 0.60/0.77/0.95, pDockQ2 0.23/0.63/0.96, n_res 76/732/2524;
  hetero ipTM 0.16/0.80/0.96, ipSAE 0.60/0.75/0.93, pDockQ2 0.23/0.51/0.96, n_res 162/674/2564.
- `download_manifest.csv` = candidates + `status` (`downloaded` / `present` / `download_error`), `detail`, `source_url`, `output_file`
  (relative path from `/cluster/scratch/gmagnusson/`). Each model appears twice (first run `downloaded`, rerun `present`), plus a few `download_error` rows
  (human: 18 homo, 42 hetero rows); dedupe by `modelEntityId`; verify with file existence. Human sizes: homodimers 2.4 GB.

```python
import pandas as pd
R="/cluster/scratch/gmagnusson/dimers/afdb_human_dimers"
c=pd.concat([pd.read_csv(f"{R}/{t}/candidates.csv").assign(kind=t) for t in ("homodimers","heterodimers")])
c["path"]=R+"/"+c.kind+"/"+c.modelEntityId+".pdb"
import gemmi; st=gemmi.read_structure(c.path.iloc[0]); [ (ch.name,len(ch)) for ch in st[0] ]
```

## 2. AFDB pocket atlas: `/cluster/scratch/gmagnusson/pocketcompare_afdb`
Documentation: `POCKETOME.md` (read in full; key points below). Pipeline `afdb/worker.py` from github.com/zoidpuff/pocketcompare (not mounted here).

- **Monomers only** (single-chain AFDB models; interface pockets absent; no ligands). Dimers were NOT processed.
- Detector: BURROW (AF mode prunes pLDDT < 50), fingerprints: pocketcompare `A4b:n1_linear`, 384-d unit vectors (float16), similarity = dot product.
  Comparable only within the same method record (Parquet schema metadata; config_sha256 b3ea8fe9...).
- Done here: tier 1 (46 reference proteomes, v6) = 555,846 structures, 1,788,750 pockets (162 units); tier 2 (Swiss-Prot v6 not in tier 1) = 436,886 structures,
  1,442,678 pockets (110 units). Tier 3 (rest of AFDB v4/v3, ~214 M structures) is planned (`worklist.parquet`, 48,341 units) but **not done**: `pockets/` has only tier1, tier2
  (`state/done` has 272 files). `foldcomp/afdb_uniprot_v4*` (~938 GB, mtime still changing 2026-10-29) is an in-progress download for tier 3; do not rely on it.
- Pocket table `pockets/tier{1,2}/<unit:06d>.parquet` (26 columns): `pocket_id` (`AF-<acc>-F<frag>-model_v6:m<i>`), `accession`, `fragment`, `taxid`, `organism`,
  `burrow_score` (>=0.03), `n_points`, `volume_A3`, `site_n_points`, `site_volume_A3`, `center_x/y/z`, `residues` (`A:53 A:79 ...`, UniProt numbering for fragment 1),
  `plddt_mean/min/frac_lt70` (lining residues, mean >= 70), `fingerprint` (fixed_size_list<halffloat>[384]).
  **No druggability score** other than `burrow_score`/volume; no ligand predictions.
- Structure table `structures/tier*/*.parquet`: per-structure `n_res`, plddt stats, `n_pockets_all`, `n_pockets_kept`, `status`.
- Other: `worklist.{parquet,tsv,json}`, `skip_accessions.txt` (989,696 acc -> tier), `ncbi/{rankedlineage,merged}.dmp`, `tar_index/` (94 files), `gcs_listing.tsv.gz`, `logs/`, `state/`.
- Join key to dimers: `accession` (UniProt). Pockets are of the monomer model; to ask "does the pocket overlap the dimer interface" compute residue overlap with the dimer
  model (residues are `A:<UniProt pos>`), noting fragments > 1 for long human proteins and the pLDDT >= 70 / score >= 0.03 filters.

```python
import pyarrow.dataset as ds, pyarrow.compute as pc
d=ds.dataset("/cluster/scratch/gmagnusson/pocketcompare_afdb/pockets",format="parquet")
t=d.to_table(columns=["pocket_id","accession","residues","burrow_score","volume_A3"],filter=pc.field("taxid")==9606).to_pandas()
```
(venv with pyarrow/pandas: `.venvs/recon`, created with `uv venv .venvs/recon; uv pip install pyarrow pandas`.)

## 3. PDB: `/cluster/scratch/jjaenes/pdb`
- mmCIF, gzipped: `<id[1:3]>/<id>.cif.gz` (e.g. `01/101m.cif.gz`, `1a/11aa.cif.gz`; 1294 shard dirs, characters 2-3 of the lowercase ID). 261,139 files, no other files.
  Newest files dated 2026-10-09 (~current PDB). No SIFTS, no manifests, no components.cif, no pre-built assembly files.
- Asymmetric-unit coordinates, but biological-assembly *definitions* are in each file: `_pdbx_struct_assembly` (oligomeric_details e.g. dimeric/monomeric, oligomeric_count),
  `_pdbx_struct_assembly_gen` (asym_id_list, oper_expression), `_pdbx_struct_oper_list`. Generate assemblies with gemmi (`gemmi convert --assembly=1`
  / `st.make_assembly`) or biotite.
- UniProt per entity: `_struct_ref` (db_name `UNP`, `pdbx_db_accession`) and `_struct_ref_seq` (chain `pdbx_strand_id`, seq_align_beg/end, db_align_beg/end). Gives author-level UniProt mapping;
  SIFTS (residue-level) must be downloaded (section 8).
- Count example: `find /cluster/scratch/jjaenes/pdb -type f | wc -l` (took ~1 min).

## 4. AFDB copy: `/cluster/work/beltrao/shared/alphafold/v6`
- Model version **v6**, `*_v6.tar` per proteome: 46 reference proteomes (`UP000005640_9606_HUMAN_v6.tar`, ...YEAST, MOUSE, ECOLI, ... ; name = UPID_taxid_CODE) + `swissprot_cif_v6.tar` (40 GB) and `swissprot_pdb_v6.tar` (29 GB)
  (48 tars incl. swissprot, ~143 GB total) + `msa_depths.csv` (241,070,489 rows, `accession,msa_depth`; not checksummed).
- Inside tars: `AF-<acc>-F<n>-model_v6.cif.gz` and `.pdb.gz` per fragment (yeast tar: 6055 of each). Read with `tarfile` / `tar xOf <tar> <member>`; `pocketcompare_afdb/tar_index/*.parquet` has member offsets.
- Monomers, UniProt accession in file name.

## 5. AlphaFold3 MSAs: `/cluster/project/beltrao/shared/`
- `25.05_alphafold3_msas_human`: 20,293 files `<acc_lowercase>_data.json.gz`; `25.06_alphafold3_msas_yeast`: 6,039 files. No other files.
- Keys: `dialect`("alphafold3"), `version`(2), `name`(accession lowercase), `sequences`, `modelSeeds` ([4]), `bondedAtomPairs`, `userCCD`.
  `sequences[0].protein` has `id`, `sequence`, `modifications`, `unpairedMsa` (string, A3M), `pairedMsa`, `templates` (4).
  Example human a0a024r1r8: seq len 64, unpairedMsa 261 kB. Yeast a0a023pyf4: len 145. Single-chain MSAs (monomers); for dimers, build a 2-protein input
  by combining two files (see `scripts/make_af3_input.py` in this workspace).

## 6. AlphaFold3 model weights: `/cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst` (1.02 GB, zstd-compressed). Weights licence: AF3 terms of use (non-commercial, not redistributable).

## 7. Conda envs: `/cluster/project/beltrao/jbuhmann/conda/envs`
bioagent_probe, hf, intellifold-host (contains `alphafold3` 3.0.2, jax 0.9.1, rdkit), opendde (biopython, biotite, gemmi, torch, rdkit), openfold3 (openfold3 0.4.1, torch 2.10, gemmi, biotite),
paper-reading, pooled-ppi-2026-08 (biopython, pandas, pyarrow), ppi-analysis (biopython, pandas, pyarrow).
No env named burrow/pocketcompare/adock and no foldseek/mmseqs/fpocket/vina binaries found in any env `bin/`. Workspace `tools/` has `adock`, `alphafold3`, `burrow` (with own README/DATA.md).

## 8. Internet access (checked 2026-10-29 with curl from the container; HTTP status of final response)
200: rest.uniprot.org (`/uniprotkb/P69905.fasta`), files.rcsb.org, data.rcsb.org, www.ebi.ac.uk/pdbe/api (`/mappings/uniprot/1cbs`), alphafold.ebi.ac.uk/api,
files.wwpdb.org (`components.cif.gz` = ligand dictionary), ftp.ebi.ac.uk SIFTS (`.../msd/sifts/flatfiles/tsv/uniprot_pdb.tsv.gz`), huggingface.co.
403: `https://alphafold.ebi.ac.uk/files/AF-0000000065760335-model_v1.pdb` (complex model file; the dimers were downloaded earlier by the owner, so just use the scratch copies).
=> SIFTS and the PDB Chemical Component Dictionary can be downloaded; do so into `data/raw/` with a script (record URL and date).

## Reproducibility
`scripts/data_checksums.sh` writes `data/raw/CHECKSUMS.sha256` (sha256 of candidates/manifests, worklist, skip list, listing, taxdump, POCKETOME.md).
