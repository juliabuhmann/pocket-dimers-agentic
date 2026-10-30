# Independent verification: PDB site descriptors, Q2.4, iface-ligand-matched, q23a curated

Verifier code (own, no producer scripts imported): `reviews/verify_misc/v*.py`; outputs `reviews/verify_misc/*.tsv`.
Status of this review: single pass by one agent; numbers below were recomputed, not eyeballed.

## 1. docs/PDB_SITE_DESCRIPTORS.md
Source: `data/interim/pdbsite/sites.parquet` (11,983 sites, per-site table incl. weight `w`). Own weighted Cliff's delta (pairwise, weights w_i*w_j), 1,000-resample cluster30 bootstrap (own seed), `v1a_cliffs.py`.

| number | doc | mine [95% CI] | verdict |
|---|---|---|---|
| volume (merged_volume_A3) | +0.11 | +0.107 [+0.050, +0.164] | VERIFIED |
| charged (merged_frac_charged) | 0.00 | -0.0002 [-0.059, +0.056] | VERIFIED |
| hydrophobic (merged_frac_hydrophobic) | +0.03 | +0.033 [-0.026, +0.087] | VERIFIED |
| burrow_score (located member) | -0.03 | -0.031 [-0.094, +0.032] | VERIFIED |

Note: the "volume" row is the merged-pocket volume; the located-member volume_A3 is +0.078, so the +0.11 depends on which pocket unit is used (the doc says "merged pocket", fine).

**Lining fractions of ligand-contact residues from raw mmCIF** (`v1b_contacts.py`; 20 random sites, seed 0; first biological assembly via gemmi, ligand re-found by CCD + nearest to stored centroid, protein residues with any heavy atom <= 4.0 A): hydrophobic and charged fractions identical to the stored `contact_frac_*` for 18/20 sites (n residues identical). The two others (1htg, ligand match 2.1 A off in my re-finding; 8d88, 0.09 A) differ by <=0.05 / <=0.02 and are plausibly my ligand-instance choice; overall corr 0.998 / 0.999. VERIFIED.
**DISCREPANT (documentation only):** the stored "contacts" list is the 4 A set (all residues with n_simple >= 1), not the vdW+0.5 A BioLiP set that DESIGN.md names primary. A vdW+0.5 rule gives fewer residues (e.g. 17 vs 18, 21 vs 33) and different fractions (differences up to 0.25 in a few sites). The doc calls them "census contacts"; the definition used should be stated. The MBS-vs-SBS contrast was not recomputed with the vdW rule.

**30 A crop** (`v1c_prep.py`, CLI `pocketcompare describe --ligand --locate burrow --backend fast`, 10 random sites among those where the crop removed chains, assembly <=24 chains):
- 2,500/11,983 sites (21%) had chains cropped; cropped fraction is higher for MBS (24.4%) than SBS (17.8%), so the crop is differential across the contrast.
- Control: my own 30 A crop through the CLI reproduces the stored burrow_score/volume/d_match exactly for 10/10, so the route is equivalent.
- Uncropped vs stored: located-pocket burrow_score changes substantially (median |delta| 0.09, Spearman between uncropped and cropped only 0.13; e.g. 3pci 0.13 vs 0.59, 6vaw 0.77 vs 0.11, 6hlr 0.74 vs 0.41), volume changes in 8/10, D-match fails in 2/10 uncropped (both cropped runs matched), residue-set Jaccard median ~0.68, 0.0 where the pocket changed. Details: `v1c_crop_compare.tsv`.
- DISCREPANT with the doc statement "BURROW features are local ... not tested": it is now tested on 10 sites and the located pocket is NOT robust to the crop where chains are removed. Sample is small and selected (cropped sites only); effect on pooled deltas is unknown, but because the burrow_score H2-secondary and the pocket-volume contrasts rely on those sites, an uncropped rerun of a stratified subsample (or all 2,500 cropped sites) is recommended before publishing. Chemistry fractions were not recomputed uncropped.

## 2. docs/Q24_RESULTS.md
`data/interim/q24/q24_pocket_scores.parquet` (80,557 rows, 2,341 units). conserved = oc_atlas >= 0.5; own burrow_score quintiles over unique pockets (agree with stored score_q up to boundary ties); matched difference weighted by interface quintile shares; bootstrap over `unit` (500 resamples), `v2_q24.py`.

| number | doc | mine | verdict |
|---|---|---|---|
| ALL matched interface - single | +0.030 [0.019, 0.041] | +0.030 [0.019, 0.042] | VERIFIED |
| homo | +0.006 | +0.006 [-0.011, 0.023] | VERIFIED |
| hetero | +0.047 | +0.047 [0.032, 0.063] | VERIFIED |
| raw ALL (interface 0.75 vs single 0.69) | 0.75/0.69 | 0.751/0.692 | VERIFIED |

RBH check: my own mmseqs-hit RBH (e <= 1e-3, mutual best bitscore) reproduces the per-species counts in the doc (mouse 4,290, rat 3,906, ... all 12). 5 random pairs (seed 0; `v2b_rbh_5.tsv`, UniProt gene names by REST): mouse Q9NWW5-Q3U466 CLN6/Cln6, rat Q9H444-D4A9Z8 CHMP4B/Chmp4b, mouse EGFR/Egfr, zebrafish CRYZL1/cryzl1, mouse SAMD13/Samd13: 5/5 same gene name, identity 0.60-0.97. Extra 70 eukaryote pairs: mouse 93%, rat 78%, zebrafish 75% same gene name; fly/worm/yeast/arabidopsis low (names differ by nomenclature, not tested by sequence). VERIFIED for plausibility; non-vertebrates not assessed.

## 3. docs/IFACE_LIGAND_MATCHED.md
Own rebuild from `pocket_match.parquet` + `dimer_pocket_features_units.parquet`: atlas_pass & is_sym_rep, interface = >=3 residues per chain, IM = interface & match_mbs_prim (173 pockets, 126 clusters), IU = interface & ~match_any (1,193; 501) (group sizes identical to doc). Own L2 logistic (L-BFGS), 5-fold CV split by paralog unit, 3 different fold seeds, 500-resample cluster bootstrap (`v3_auc.py`).

| feature | doc AUC | mine (3 fold seeds) [CI] | verdict |
|---|---|---|---|
| burrow_score | 0.87 [0.84, 0.90] | 0.871 / 0.864 / 0.861 [0.835, 0.902] | VERIFIED |
| charged | 0.51 | 0.49 [0.44, 0.53] | VERIFIED (no information) |
| charged + burrow | 0.87 | 0.869 / 0.865 / 0.860 | VERIFIED |
| volume | 0.77 | 0.778 / 0.770 / 0.760 | VERIFIED |

Caveat (not a discrepancy): the 0.87 is largely by construction. Matching uses a lenient overlap rule that favours large, well-scored pockets, so burrow_score separates IM from IU. IM/IU definition itself was taken from `match_*` flags, not re-derived from the H4 data.

## 4. results/part2/q23a_curated.tsv (114 rows)
5 random rows (seed 0): AF-0000000210462152:m4, AF-0000000066625148:m0, AF-0000000211655060:m4, AF-0000000066611529:m0, AF-0000000065828223:m4. Each re-evaluated from `dimer_pocket_features_units.parquet` and the raw model PDBs (`v4_helix.py`, pydssp numpy backend, own parsing):
- F1 score >= p60: own p60 of the 69,857 atlas-pass pockets = 0.11269 (matches 0.1127); all 5 above (0.157-0.545). PASS
- F2 lining pLDDT >= 80 both chains (min 83.5), F3 >= 5 residues per chain (min 7), F5 volume 300-1500 (401-971), F6 max ipSAE >= 0.7 (0.741-0.938). PASS
- F7a max_overlap_iso < 0.5 (max 0.455), is_sym_rep True, atlas_pass True, dimer_only True; F7b atlas overlap < 0.5 on both chains or no atlas protein (one row atlas_has_protein=False, vacuous) (max 0.42). PASS
- F4 helix fraction from raw model: 0.323, 0.326, 0.269, 0.343, 0.021, identical to stored; all < 0.7. PASS
- Dedup and the PDB-complex/ChEMBL annotation columns were not re-derived. VERIFIED for the 5 rows and the cascade's first-step threshold; the full 114-row set was not re-evaluated.

## Summary
All headline numbers VERIFIED. Two DISCREPANT/doc-level issues: (a) the ligand-contact definition is the 4 A set, not vdW+0.5 A; (b) the 30 A crop is not innocuous for the located BURROW pocket (10-site test) and is applied more often to MBS sites. 
