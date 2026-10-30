# Part 1 results: SBS vs MBS ligand sites in the PDB (family A, H1–H3)

**Status: UNVERIFIED until independent re-derivation** (DESIGN §5). All numbers come from one run of the scripts below
on census **v2** (`data/part1/`, census_version = v2). Nobody has checked them by hand yet. Every results TSV starts with a
comment line holding: the frozen design sha256 (`e3f9b550…eee28`, which matches the current DESIGN.md), the git commit,
`code_dirty`, the census version, the size and mtime of the ligand_sites file, the ChEMBL release and whether a drug-like
FALLBACK was used. The same fields are in `results/part1/analysis/provenance.json`.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/ext/bin/python scripts/fetch_external.py --only chebi_relation chebi_compounds chebi_dbacc chebi_source chebi_relation_type
scripts/fetch_chembl_sqlite.sh                       # ChEMBL 37 SQLite (tarball sha256 33c20374…d281)
cp data/raw/chembl_sqlite/chembl_37/chembl_37_sqlite/chembl_37.db $LOCAL/   # optional: random reads on /cluster/work are ~100x slower
.venvs/part1/bin/python scripts/build_druglike_flags.py --classes data/part1/ligand_classes.parquet [--db $LOCAL/chembl_37.db]  # ~1 min
.venvs/part1/bin/python scripts/part1_analysis.py    # ~4 min -> results/part1/analysis/*.tsv, data/part1/analysis/units_primary_v2.parquet
.venvs/part1/bin/python scripts/part1_figures.py     # -> results/part1/figures/*.png
```
`part1_analysis.py --census-dir data/part1/v1` also runs on census v1. The column map is chosen automatically.

## 1. Drug-like flags (`data/processed/ccd_druglike.parquet`, + `.json`)
- **Source.** ChEMBL **37** SQLite (release 2026-05-01), with CCDs mapped to ChEMBL through UniChem (`ccd_to_chembl.parquet`). There is **no FALLBACK**: all 20,683 mapped ChEMBL IDs were found in the SQLite.
- **Primary definition.** A CCD is drug-like when all of these hold:
  - it maps to ChEMBL;
  - it has max_phase ≥ 1 OR a pChEMBL value ≥ 6 against any target (molecule or its parent, max taken);
  - it has ≥ 12 heavy atoms (census `DRUGLIKE_MIN_HEAVY_ATOMS`);
  - it is not a cofactor, metabolite or buffer.
- **Unmapped CCDs** are never drug-like. Ro5+QED is only the secondary definition, kept in its own column.
- **Cofactor, metabolite and buffer exclusion:**
  - census class precedence: the metal, inorganic, heme, PDBe-CoFactor, nucleotide-like, sugar and lipid classes win over drug-like;
  - BioLiP artifact list or PLINDER badlist;
  - **ChEBI roles**. CCD→ChEBI comes from the PDBeChem xrefs plus ChEBI's ChEMBL xrefs. Roles are expanded over conjugate acid/base and tautomer links, then over is_a ancestors. A CCD is excluded if it has a role under "mammalian / *E. coli* / *S. cerevisiae* / fundamental metabolite" and not under "xenobiotic / drug metabolite", or a role under "buffer".
- **Deviation (logged).** The census class `amino_acid` (CCD type "peptide-linking") may be overridden by drug-like. Without this, ritonavir (RIT) was classed `amino_acid`. Free amino acids are still removed by the ChEBI metabolite rule.
- **Counts:**
  - CCDs meeting the ChEMBL rule: 12,986 (phase ≥ 1: 2,027; pChEMBL ≥ 6: 13,185 before the size floor);
  - removed: 581 cofactor-class, 13 buffer/solvent-list, 116 ChEBI-metabolite;
  - **12,276 drug-like CCDs**; 23,475 Ro5+QED-secondary CCDs.
- REST-era max_phase agrees with the SQLite value for 93.5% of CCDs.
- **Sentinels pass** (the script exits non-zero otherwise):
  - drug-like: STI (imatinib), RIT (ritonavir), MTX (methotrexate), 3MI (tafamidis);
  - not drug-like: MES, EPE (HEPES), CIT (citrate), CLR (cholesterol).
- **Spot check.** I read the 40 most frequent drug-like CCDs in the primary set (by number of clusters). Most are real ligands or inhibitors: staurosporine, MTX, imatinib, indomethacin, bortezomib BO2. A few are borderline:
  - EDTA (phase 4), pyridine-dicarboxylates and zeaxanthin (phase only);
  - flavonoids (quercetin, naringenin) and nicotine.

  The definition is applied as written; these are not hand-curated out.

## 2. Analysis unit table (n at each step; `analysis_ledger.tsv`)
| Step (primary) | Ligand instances | Clusters (30%) |
|---|---|---|
| All ligand instances, assembly 1 (incl. alt conformers) | 1,932,129 | |
| Not a duplicate copy | 1,885,733 | |
| Bound (≥ 2 residues with ≥ 2 atom contacts) | 1,158,667 | |
| Not glycan | 1,083,502 | |
| Entry covered by the BioLiP2 snapshot | 972,538 | |
| Not additive (BioLiP2-relevant OR SoI) | 502,621 | 14,330 |
| No antibody-chain contact | 501,762 | |
| Oligomer (homomer/heteromer; same_uniprot_hetero excluded) | 395,673 | 10,045 |
| cluster30 assigned | 395,627 | 10,045 |
| Symmetric copies collapsed (census `site_group`, Jaccard) | 183,660 sites | 10,027 |
| **Units (cluster30 × ccd_key), representative entry** | **44,178** (33,173 homomer, 11,005 heteromer) | 10,027 |
| Organic units (metals excluded; H1/H2 set) | 34,654 | 7,252 |

**Consistency checks** (`consistency_checks.tsv`):
- The primary MBS column equals the rule recomputed from `chain_counts`: 0 disagreements.
- My oligomer primary filter equals the census `in_primary_set`: 0 disagreements in 396,465 rows.

**Choices where the design is silent** (ours):
- The representative entry is the best resolution, with ties broken by the earliest PDB ID.
- Unit MBS = any collapsed site of that ligand in the representative entry is MBS.
- Sites contacting an antibody chain are dropped. DESIGN §0 excludes antibody chains.

## 3. Descriptives (`mbs_fraction_descriptives_seed{0,1,2}.tsv`; 1,000 cluster-bootstrap resamples, seed 0)
Overall MBS fraction among oligomer ligands, all classes:

| Level | Before additive filter | After additive filter |
|---|---|---|
| Sites, symmetric-collapsed | 26.7% (24.9–28.5) | **26.3% (23.9–28.6)** |
| Non-redundant units | 32.9% (31.2–34.4) | **30.6% (28.1–33.1)** |

- **Comparison with published values** (`comparison_published_fractions.tsv`):
  - PLINDER 25% lies inside the site-level CI; Q-BioLiP 28.8% lies just above it.
  - At unit level, Q-BioLiP lies inside the CI and PLINDER lies below it.
  - The published figures are instance-level with different definitions, so this is a plausibility check only.
- **By assembly type** (units, after the filter):
  - homomer 29.2% (26.0–32.2), heteromer 34.9% (32.1–37.7);
  - organic only: homomer 32.3% (28.2–36.1), heteromer 42.1% (38.1–46.0).
- **Marsh 3-class** (units, after the filter):

  | | Cofactor | Small molecule | Metal |
  |---|---|---|---|
  | Homomer | 36.1% (33.6–38.9) | 31.2% (27.0–35.3) | 14.7% (13.6–16.0) |
  | Heteromer | 45.7% (42.2–49.2) | 40.2% (36.1–44.5) | 10.5% (9.2–11.9) |
- **Fine classes** (units, after the filter; homomer / heteromer):
  - cofactor 38.0% / 41.0%;
  - drug-like 35.0% (28.2–43.6) / 30.6% (21.8–39.8);
  - heme 39.0% / 64.7%;
  - lipid 42.6% / 70.3%;
  - nucleotide-like 34.7% / 30.4%;
  - metal ion 14.7% / 10.5%.
- **Unadjusted cofactor − drug-like MBS difference** (organic units): +4.4 percentage points (−2.8 to +10.9). Homomer +3.0 (−5.3 to +10.2); heteromer +10.4 (−0.8 to +21.1).
- **Seeds 1 and 2.** The CI bounds move by at most 1.1 and 1.9 percentage points.

## 4. Hypotheses (`hypothesis_tests.tsv`; BH over the primary H1–H3 that passed the power gate)
Power gate: ≥ 30 clusters per compared stratum. It **passed for H1, H2 and H3**.

| Test | Estimate (95% CI) | p | p (BH, family A) | n |
|---|---|---|---|---|
| **H1**: adjusted OR of MBS, cofactor vs drug-like | **0.89 (0.61–1.32)** | 0.57 | 0.57 | 34,654 units, 7,252 clusters (cofactor 1,783 / drug-like 1,211 clusters) |
| **H2**: conditional OR of drug-like occupancy, MBS vs SBS (within cluster30) | **1.17 (1.04–1.31)** | 0.0085 | 0.013 | 20,155 units in 1,733 clusters with both MBS and SBS; 544 discordant clusters |
| **H3**: OR of drug-like (vs cofactor), on-axis vs off-axis homomer MBS | **35.4 (12.6–99.3)** | 1.2e-11 | 3.5e-11 | 2,448 units, 895 clusters |

**H1.** Marsh is not replicated under the pre-registered model.
- Unadjusted MBS fraction: cofactor 38.4% vs drug-like 34.0%.
- Model: logistic regression with class + B-spline(heavy atoms, df 4) + resolution + method + n_chains + SG + release year, cluster-robust SEs on cluster30.
- Metals are excluded from H1 (DESIGN §0 forbids pooling them).
- The only class with an OR significantly different from drug-like is lipid/fatty acid: 1.71 (1.26–2.33).
- "MBS in any entry" sensitivity: 0.96 (0.66–1.39).
- **Strict MBS** (≥ 3 residues per chain) **reverses the direction**: 0.41 (0.22–0.76). Drug-like ligands are more often deep, two-sided MBS.

**H2.**
- Drug-like share is 18.8% at MBS and 17.0% at SBS.
- The test is an exact conditional logistic regression with one binary predictor, i.e. Fisher's noncentral hypergeometric likelihood per cluster, with an LR p-value. Mantel–Haenszel cross-check: 1.16 (1.04–1.30).
- The effect is **not robust**:
  - leave-out top-10 clusters: 1.08 (0.95–1.22), p = 0.26;
  - post-BioLiP fallback set: 0.92 (0.62–1.34);
  - Ro5+QED drug-like: 1.05 (0.95–1.16);
  - strict MBS strengthens it: 1.51 (1.32–1.74).

**H3.**
- 2×2 table (units):

  | | Drug-like | Cofactor |
  |---|---|---|
  | On-axis | 352 | 9 |
  | Off-axis | 1,095 | 992 |
- Fisher OR 35.4. The model is a logistic regression with cluster-robust SEs.
- **Caveat on the power gate.** The comparison strata pass (on-axis 72 clusters, off-axis 858). But the on-axis cofactor cell has only 9 units in 8 clusters. If the gate is read per 2×2 cell, H3 would be descriptive.
- **Concentration.** 194 of the 352 on-axis drug-like units come from two cluster30 groups: cluster 93, representative 3KDB, 127 units, and cluster 259, representative 4HIQ, 67 units. These look like HIV-1 protease-type C2 sites (not checked per entry). Leaving out the top-10 clusters gives an OR of 18.8 (8.1–43.7).
- Geometric on-axis rule (centroid ≤ 3 Å from the C2 axis): 34.7 (11.2–108). Jaccard OR geometric (EXPLORATORY): 35.3.

## 5. Sensitivities (DESIGN §2 + PI additions; all in `hypothesis_tests.tsv`, forest plot `figures/forest_h1_h3.png`)
H1, H2 and H3 under each variant:

| Variant | H1 OR | H2 OR | H3 OR |
|---|---|---|---|
| strict MBS | 0.41 | 1.51 | 65 |
| 4 Å contacts | 0.98 | 1.12 | 35 |
| PLINDER additives | 0.90 | 1.16 | 32 |
| BioLiP geometric additives | 0.90 | 1.13 | 32 |
| author/software assembly agree | 0.90 | 1.18 | 38 |
| exclude photosynthetic | 0.89 | 1.19 | 35 |
| post-BioLiP FALLBACK set | 0.85 | 0.92 | 3.6 (gate fails: on-axis 12 clusters → descriptive) |
| Ro5+QED drug-like | 1.01 | 1.05 | 21 |
| leave-out top-10 clusters | 1.10 | 1.08 | 19 |
| leave-out top-10 ligands | 0.97 | 1.18 | 34 |
| ≥ 1-atom residue rule (EXPLORATORY) | 0.94 | 1.17 | 35 |

The left-out clusters and ligands are listed in `leave_out_top10_{clusters,ligands}.tsv`.

## 6. Not done / limitations
- **The H2 secondary (burrow_score of the pocket at the ligand) was not run** because of wall-clock budget. No pocket descriptors were computed, and fpocket was not attempted.
- **Not pre-registered (our choices):**
  - n_chains enters H1 as a linear term;
  - the spline has df = 4;
  - the H3 gate uses comparison strata rather than cells;
  - drug-like occupancy in H2 is drug-like vs all other organic ligands.
- **Census v2 limits carry over.** On-axis is computed on the top-2 chains only. Homomer typing follows the coverage deviation accepted by the PI. The additive labels carry BioLiP2 curation noise.
- **Figures** are PNG at dpi 150:
  - `mbs_fraction_by_class_units.png`, `…_units_before_additive.png`, `…_sites.png`;
  - `forest_h1_h3.png`.
