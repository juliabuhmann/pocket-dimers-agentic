# Review r2 — data sources & algorithms (DESIGN v2)

Checked 2026-10-29 from the cluster with curl HEAD/GET only; nothing downloaded. Sizes are Content-Length.

| Source | URL | Size | Format | Status |
|---|---|---|---|---|
| BioLiP2 annotation | zhanggroup.org/BioLiP/download/BioLiP.txt.gz | 78.6 MB | TSV | OK. Contains only ligands passing its relevance filter. "Annotated relevant" = row present, so additives = PDB ligands absent from it. Records weekly release, so pin sha256 |
| Q-BioLiP | yanglab.qd.sdu.edu.cn/Q-BioLiP/Download | page 40 KB | HTML (links) | Reachable, file names not checked; use BioLiP2 as primary |
| UniChem PDBe(3)->ChEMBL(1) | ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src3.txt.gz | 144 KB | TSV | OK (dated 2026-10-24) |
| ChEMBL max phase | .../chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz (or chemreps 279 MB) | 5.4 GB | SQLite | Problem: chemreps has no max_phase/activities. Use the SQLite (one-off, off-scratch) or the REST API for the few thousand mapped IDs |
| SIFTS chain-level | .../sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz | 6 MB | TSV | OK, chain-level only |
| SIFTS per-residue (a) | ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/<id>.xml.gz (1cbs 18 KB) | ~18 KB per entry | XML | OK, ~100k entries needed, parallel fetch is ~2 GB total |
| SIFTS per-residue (b) | www.ebi.ac.uk/pdbe/entry-files/download/<id>_updated.cif | ~190 KB/entry | mmCIF with _atom_site.pdbx_sifts_xref_db_{name,acc,num,res} | OK. **The local wwPDB mmCIFs have no _pdbx_sifts columns (checked 1cbs)** |
| SIFTS (c) PDBe API | www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/<id> | 0.3 KB | JSON | Segments only, not per-residue. Use only as QC |
| AlphaMissense | storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz | 1.21 GB | TSV (UniProt, protein-level) | OK. The hg38 name in the notes returned 404. Use aa_substitutions; no coordinate liftover needed |
| ClinVar | ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz; tab_delimited/variant_summary.txt.gz | 201 MB; 457 MB | VCF / TSV | OK, but the protein position needs HGVS p. parsing. Prefer the UniProt variation API (stated in DESIGN) with the ClinVar snapshot date recorded |
| LIGYSIS | zenodo.org/api/records/14645504 | 24 KB metadata; LIGYSIS_human_sites_APRIL_2024.pkl 0.9 MB, ..._ress_binary_labels.pkl 97 MB | pickle | OK, small. Avoid the 3.3 GB MASTER_POCKET_SHAPE and 350 MB COMPARISONS files. Pickles need a compatible pandas, and loading them runs untrusted code (the record is public, pin the md5) |
| PLINDER annotation_table | storage.googleapis.com/plinder/2024-06/v2/index/annotation_table.parquet | unknown | parquet | **Problem: anonymous HTTPS gives 403 and bucket listing 401.** Not verified as fetchable. Try the `plinder` package/gcsfs with token="anon", else drop (it is already a "cut" item) |
| ECOD | prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt | 663 MB (redirect, Content-Length as reported) | TXT | OK but large. For the kNN split, a Pfam-clan split from pdb_chain_pfam.tsv.gz (4 MB) or pdb_chain_scop_uniprot is cheaper |
| CCD | files.wwpdb.org components.cif.gz | already noted in DATA_RECON | CIF | OK |

## Algorithm notes
- Per-residue SIFTS: best is (a) the XML, or (b) _updated.cif. Option (b) also gives UniProt residue names for the 90% identity check. Do not use struct_ref_seq offsets.
- Obsolete or superseded entries: the local snapshot is 2026-10-09, SIFTS is 2026-10-25. Record the mismatch and drop entries with no SIFTS mapping, in the ledger.
- BioLiP2 join key: PDB id + ligand chain + residue number. Its chain ids are author ids, and its assembly handling differs from "assembly 1", so check on sentinels that the join does not drop on-axis homodimer sites.

## Proposed changes
MUST
1. Replace "SIFTS per-residue download" with a scripted per-entry fetch of sifts/xml (or _updated.cif) restricted to entries in the census that have a ligand and are in a human pair. Record the date.
2. AlphaMissense: use AlphaMissense_aa_substitutions.tsv.gz (1.2 GB). The ChEMBL source must be stated as SQLite (5.4 GB) or API, and the ChEMBL release (37) recorded.
SHOULD
3. PLINDER: test anonymous access with gcsfs first; if not working today, drop the cross-check and the PLINDER-additive sensitivity (it needs the same table).
4. Swap ECOD for a Pfam-clan split unless ECOD is already on disk.
