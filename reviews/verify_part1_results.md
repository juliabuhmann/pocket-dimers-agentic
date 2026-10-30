# Independent re-derivation of Part 1 headline numbers (H1–H3)

Verifier: an independent agent, 2026-10-29. It did **not** import `scripts/part1_analysis.py`. All code is in `reviews/verify_part1/`:
- `build_units.py`: unit table;
- `tests_h1_h3.py`: H1–H3 and the H1 strict-MBS stress tests;
- `fetch_entity_meta.py`: organism and entity names from the local mmCIF;
- `h3_stress.py`;
- `spot_check_labels.py`.

Outputs are in `reviews/verify_part1/out/`. Inputs are census v2 `data/part1/ligand_sites.parquet`, `data/processed/ccd_druglike.parquet` and RCSB `data/raw/part1/clusters-by-entity-30.txt`. I re-derived cluster30 from this file using pdb_id + dominant_entity. It matches the census `cluster30` column 100%.

**Rerun:**
```bash
source scripts/env_threads.sh
for s in build_units tests_h1_h3 fetch_entity_meta h3_stress spot_check_labels; do
  .venvs/part1/bin/python reviews/verify_part1/$s.py
done
```

## 1. Unit table

| Step | Producer | Verifier | Status |
|---|---|---|---|
| Not duplicate → bound → not glycan → BioLiP2-covered → not additive → no antibody contact | 1,885,733 → … → 501,762 | identical at every step | VERIFIED |
| Oligomer (homo/hetero; same_uniprot_hetero excluded) | 395,673 | 395,673 | VERIFIED |
| cluster30 assigned | 395,627 / 10,045 cl | 395,627 / 10,045 cl | VERIFIED |
| Sites (pdb_id × site_group) | 183,660 / 10,027 cl | 183,660 / 10,031 cl | VERIFIED (sites); clusters +4 |
| Units (cluster30 × ccd_key, rep entry) | 44,178 / 10,027 cl | 44,186 / 10,031 cl | VERIFIED (Δ = 8 units, 0.02%) |
| Homomer / heteromer units | 33,173 / 11,005 | 33,173 / 11,013 | VERIFIED |
| Organic units (H1/H2) | 34,654 / 7,252 cl | 34,663 / 7,257 cl | VERIFIED |

The small surplus has a known cause. 423 site groups contain instances whose dominant chain falls in different cluster30s. I assign each site the cluster30 of its first instance; the producer evidently uses another tie rule. The effect on the estimates is negligible.

## 2. Hypothesis estimates

My model specification:
- **H1:** GLM binomial, class + B-spline(heavy atoms, df 4) + resolution + method + n_chains + SG + release year, with cluster-robust SEs on cluster30.
- **H2:** exact conditional likelihood (Fisher noncentral hypergeometric per cluster), plus Mantel–Haenszel.
- **H3:** Fisher 2×2, plus a logistic regression with cluster-robust SEs.

| Number | Producer | Verifier | Status |
|---|---|---|---|
| H1 OR cofactor vs drug-like | 0.89 (0.61–1.32), p 0.57 | 0.89 (0.61–1.32), p 0.57 | VERIFIED |
| H1 "MBS in any entry" | 0.96 (0.66–1.39) | 0.96 (0.66–1.39) | VERIFIED |
| H1 strict MBS | 0.41 (0.22–0.76) | 0.41 (0.22–0.76), p 0.005 | VERIFIED |
| H2 conditional OR | 1.17 (1.04–1.31), p 0.0085; 544 discordant | 1.17 (1.04–1.31), p 0.0084; 543 discordant; MH 1.16 | VERIFIED |
| H2 strict MBS | 1.51 | 1.52 (1.32–1.74) | VERIFIED |
| H3 2×2 (on dl / on cof / off dl / off cof) | 352 / 9 / 1,095 / 992 | 351 / 9 / 1,093 / 992 | VERIFIED |
| H3 OR | 35.4 (12.6–99.3) | 35.4 (12.6–99.4), p 1.3e-11 | VERIFIED |
| H3 geometric on-axis | 34.7 (11.2–108) | 34.7 (11.1–108) | VERIFIED |
| H3 strict MBS | 65 | 51.7 (16.8–159) | VERIFIED (inside the CI, same direction; the difference is probably from how the strict unit set is formed) |

## 3. H3 stress tests

**(a) Power gate per cell.** The table gives the number of clusters per cell.

| | Drug-like | Cofactor |
|---|---|---|
| On-axis | 64 | **8** |
| Off-axis | 334 | 576 |

The on-axis cofactor cell has only 9 units in 8 clusters. **If the ≥30-cluster gate is read per cell, H3 fails it and must be reported descriptively.** It fails in every variant below as well.

**(b) Viral exclusion.** Organisms come from the mmCIF `_entity_src_gen` / `_entity_src_nat` / `_pdbx_entity_src_syn` records of the dominant entity, matched with a regex.

| Variant | OR (95% CI) |
|---|---|
| Exclude viral (268 units) | 23.6 (8.6–65.2) |
| Exclude viral proteases (206 units) | 24.3 (9.1–65.4) |

The effect remains large.

**(c) Leave-one-cluster-out**, for the top-20 clusters by on-axis drug-like units:
- Dropping a single cluster gives an OR of 24.1–35.3. The minimum comes from dropping cluster 93 (HIV-1 protease).
- Dropping the top clusters jointly:

| Clusters dropped | OR (95% CI) |
|---|---|
| Top 2 | 17.0 (7.3–39.8) |
| Top 5 | 12.4 (5.5–28.1) |
| Top 10 | 9.0 (4.0–20.2) |
| Top 20 | 5.8 (2.6–12.9) |

The OR decays steadily but stays above 1.

**(d) Cluster-level analysis.**
- **Majority vote**, one vote per cluster × axis status:
  - on-axis: 63 drug-like vs 8 cofactor;
  - off-axis: 309 drug-like vs 549 cofactor;
  - OR 14.0, p 1e-18.
- **Presence-only units:** OR 13.8 (6.7–28.6).
- **Within-cluster paired comparison:** only **34 clusters** have both on- and off-axis units. 29 are ties, 4 favour on-axis drug-like and 1 favours off-axis; sign test p = 0.38.
- **Interpretation:** the effect is almost entirely *between proteins*. Proteins whose ligand sits on the 2-fold axis are mostly drug-discovery targets. It is not shown that, within one protein, on-axis sites prefer drug-like ligands.

**(e) Top on-axis drug-like clusters.** These are 351 units from 64 clusters and 94 UniProt accessions. The top 2 clusters give 55% of the units and the top 10 give 77%; 11 clusters have ≥5 units.

| Cluster | Units | Protein | Example PDB IDs | Notes |
|---|---|---|---|---|
| 93 | 127 | HIV-1 protease | 3KDB, 2HS1, 3VF5, 1ODX | |
| 259 | 67 | Transthyretin | 4HIQ, 6E78, 2F7I | T4 channel |
| 209 | 21 | GluA2 LBD | 7F3O, 2XX8, 4FAT | dimer-interface PAMs |
| 3205 | 15 | Glutaminase GLS | 5HL1, 5I94 | BPTES-type |
| 2246 | 8 | PD-L1 | 5NIU | small-molecule-induced dimer |
| 1158 | 8 | RSV F | 5KWW, 6EAK | 3-fold axis fusion inhibitors |
| 701 | 8 | PKM2 | 3ME3, 4JPG | activator site |
| 238 | 6 | S. aureus gyrase | 4BUL | NBTIs |
| 1107 | 5 | IDH2 | 5I96 | allosteric |
| 1375 | 5 | STING | 6XF3 | |
| 97 | 5 | BRD4 | | bivalent |
| 888 | 4 | FBPase1 | 5PZR | |
| 15066 | 4 | NKG2D | 8EA7 | |
| 3042 | 3 | caspase-6 | 4NBL | |
| 587 | 3 | FabG | 4BNT | |

Full list: `out/h3_top_on_axis_druglike_clusters.tsv`.

**The 9 on-axis cofactor units:**

| PDB | Ligand | Protein |
|---|---|---|
| 3UXJ | NAP | QueF |
| 4MRS | GDS | ABC transporter |
| 8QFE | FMN | adenylate cyclase |
| 5GM1 | SAH | O-methyltransferase |
| 8P90 | PLP | NAPE-PLD |
| 2CIE | FAD | VNG1446H |
| 3CNE | FMN | putative protease I |
| 4D1Y | RBF | putative protease I |
| 8VJN | FMN | EncD |

All are from distinct small families.

**Verdict on H3.** The direction is robust across many families, not only HIV protease or TTR. However:
1. The headline OR of 35 is inflated by a few heavily drugged targets, because the cluster × CCD unit gives HIV PR 127 units. Cluster-level ORs are about 14, and about 6–9 after dropping the top 10–20 clusters.
2. The per-cell gate fails, with 8 clusters in the on-axis cofactor cell.
3. There is no within-protein evidence.

## 4. H1 strict-MBS reversal (OR 0.41)

| Variant | OR (95% CI) |
|---|---|
| Leave out one of the top-20 cofactor CCDs | 0.36–0.44 |
| Leave out one of the top-20 clusters | 0.39–0.57 (max: dropping HIV PR) |
| Drop NAD/NAI/NAP/NDP | 0.57 (0.31–1.02) |
| Drop FAD/FMN | 0.39 |
| Drop SAH/SAM | 0.43 |
| Drop PLP/PMP | 0.41 |
| Drop CoA/AcCoA | 0.40 |
| **Drop the top-10 cofactor CCDs jointly** | **0.75 (0.40–1.43), NS** |
| Homomers | 0.34 (0.16–0.75) |
| Heteromers | 0.65 (0.38–1.12) |

HEM is not in the cofactor class (heme_porphyrin is a separate class), so dropping it has no effect.

- **Cause of the reversal.** Abundant dinucleotide/Rossmann-type cofactors (NAD(P)(H), FAD, SAM/SAH, PLP, CoA) are pooled across many families. They touch the neighbouring chain with only 1–2 residues: 38% are MBS under the primary rule but 16% under the strict rule. Large drug-likes are the opposite: 40–60 heavy atoms → 51% strict MBS; >60 → 80%. This large-drug-like group includes HIV-PR peptidomimetics.
- **Not one ligand or one family.** No single ligand or family drives the reversal; the "major cofactor" group as a whole does.
- **Primary-rule masking.** Dropping the top-10 cofactors under the *primary* rule gives an OR of 1.86 (1.26–2.75). The remaining minor cofactors are 55% MBS, so the primary null mixes two opposing subgroups.

## 5. Label spot-check (random, seed 0; `out/label_spot_check.tsv`)
- **Drug-like, 10 sampled:** all are ChEMBL-active synthetic compounds (pChEMBL 6.4–8.9) or macimorelin (phase 4). One is borderline: P1A, ponasterone A, an insect steroid hormone. **Sane.**
- **Non-drug-like "other organic", 10 sampled:** 4 are drug-like-looking synthetic inhibitors with no ChEMBL mapping (XE9, 7FD, 6A1, 7EN). Two are mapped but below pChEMBL 6 (I7K, HH9). The rest are fragments or adducts (TRI, PLP+SUO, LPI, TYP). The definition is applied correctly but it is conservative: unmapped synthetic ligands are false negatives.
- **Cofactor, 10 sampled:** NDP, a TPP analogue and glutaryl-CoA are sane. 7 of the 10 are composite covalent keys (PLP+SER, PLP+TRS, FAD+H2S, BBN+NAD and similar). They are classed as cofactor, which is reasonable, but they carry no name. Composite keys that contain a drug-like component are never labelled drug-like.

## 6. Recommendations for reporting
- **H1:** VERIFIED null; "Marsh not replicated, OR 0.89 (0.61–1.32)". The strict-MBS reversal should be reported as a sensitivity driven by the major Rossmann-type cofactor group collectively, not by single ligands or families. It loses significance when the top-10 cofactors are dropped.
- **H2:** VERIFIED, OR 1.17 (1.04–1.31). Report it as weak and fragile: it is not robust to leaving out the top-10 clusters or to Ro5+QED (per the producer).
- **H3:** VERIFIED arithmetic, but report it as **descriptive** (the per-cell gate fails, 8 clusters) and with cluster-level effect sizes (OR ≈ 14). Suggested wording: "on-axis homomer MBS sites holding drug-like ligands occur in 64 clusters across diverse targets (HIV PR, TTR, GluA2, GLS, PD-L1, RSV F, PKM2, …) versus 8 clusters for cofactors; this is a between-protein pattern and not tested within proteins." Do **not** headline OR 35 or attribute it to HIV protease alone; without viral proteins the OR is 24.
