# Pocket descriptors at real PDB ligand sites: MBS vs SBS

**Status: UNVERIFIED (v2: uncropped primary).** The first version ran on assemblies cropped to 30 Å, and that version was partially checked in `reviews/verify_misc.md`. The review found the crop changes the located pocket, so this version reruns everything on the whole assembly 1 (PI request, 2026-10-29). The cropped run is kept as a sensitivity. The numbers below have not been re-derived independently.
- **H2-secondary (pre-registered descriptor, DESIGN.md §1 H2):** burrow_score of the BURROW pocket at the ligand.
- **EXPLORATORY:** everything else (volume, buriedness, lining physicochemistry, on/off-axis, comparison with Q2.2).

Provenance:
- DESIGN.md sha256 is `e3f9b550…eee28`, which equals the frozen tag.
- Each TSV starts with a line holding the design hash and git HEAD.
- Census v2. Units come from `data/part1/analysis/units_primary_v2.parquet`.
- Bootstrap: 1,000 resamples over cluster30. Seed 0 is primary; seeds 1 and 2 were run for the pooled contrasts and their CI bounds move by ≤0.008.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/pdbsite_sample.py --n 12000            # ~5 min -> data/interim/pdbsite/sample_n12000_seed0.parquet
source scripts/env_tools.sh pockets
python scripts/pdbsite_describe.py --sample data/interim/pdbsite/sample_n12000_seed0.parquet --workers 12 \
    --out data/interim/pdbsite/describe_uncropped                       # primary, whole assembly 1: 31 min, resumable
python scripts/pdbsite_describe.py --sample data/interim/pdbsite/sample_n12000_seed0.parquet --workers 12 \
    --crop 30 --out data/interim/pdbsite/describe_crop30                # sensitivity: 22 min
.venvs/part1/bin/python scripts/pdbsite_compare.py                    # primary -> results/part1/pocket_descriptors/*.tsv
.venvs/part1/bin/python scripts/pdbsite_compare.py --describe-dir data/interim/pdbsite/describe_crop30 --variant crop30
.venvs/part1/bin/python scripts/pdbsite_crop_sensitivity.py           # crop_vs_uncropped_*.tsv
.venvs/part1/bin/python scripts/pdbsite_figures.py                    # -> results/part1/figures/pdbsite_*.png (primary)
```
The crop30 run in `describe_crop30/` was produced by the first version of `pdbsite_describe.py`, where the crop was hard-wired. That code path is the same as `--crop 30` now, and the crop30 comparison reproduces the v1 numbers exactly.

## Method
1. **Sample** (`pdbsite_sample.py`).
   - Units: the organic oligomer units of the Part 1 primary set (34,654 units / 7,252 clusters).
   - One site per (cluster30 × CCD) unit. It is drawn at random (seed 0) from the symmetric-collapsed sites in the unit's representative entry, with the same ccd_key and the same MBS label as the unit. Every unit had a matching site.
   - Strata: MBS/SBS × ligand class. Ligand classes are drug-like (primary ChEMBL rule), cofactor (Marsh cofactor class: PDBe cofactor + heme) and other organic.
   - Allocation is equal across strata, ≤12,000 sites in total (`sample_strata_n12000_seed0.tsv`). Sizes: MBS-cofactor 1,577 (all available units), MBS-drug-like 1,815 (all), and 2,152 in each other stratum. That is 12,000 sites in 4,257 clusters.
   - The pilot of 100 sites (cropped) took 0.35 s/site aggregate with 12 workers. The uncropped full run took 0.16 s/site aggregate, 31 min; the median per site was 0.6 s and the maximum 52 s.
2. **Describe** (`pdbsite_describe.py`, pocketcompare CPU fast backend, 12 workers).
   - Structure: gemmi builds the first assembly (as the census does), keeping the first altloc only and dropping H.
   - Ligand: the census instance is re-found among het residues with the CCD component names, nearest to the census centroid. 11,983 of 11,999 lie within 1.5 Å.
   - **Primary: all polymer chains of assembly 1**, up to 60. Sensitivity: chains with no heavy atom within 30 Å of the ligand are dropped.
   - All ligands and waters are removed. The site is then described with the ligand heavy atoms as the site and `--locate burrow`, the same code path as `describe --ligand lig.sdf --locate burrow`, with the coordinates passed in-process. Two pocket units are stored:
     - **Located member** (the pre-registered "BURROW pocket at the ligand"): burrow_score, volume_A3, fingerprint.
     - **Merged pocket** that contains it: the unit of `describe` without a site, as in the AFDB atlas and Q2.2 dimer pockets. It gets merged_burrow_score (max member score), merged_volume_A3 and merged lining.
   - Lining chemistry uses exactly the classes of `scripts/q22_common.py`: hydrophobic AVLIMFWC, aromatic FWYH, charged DEKR. MSE is counted as M.
   - Buriedness is mean(1 − min(rel. SASA, 1)) over lining residues in the described complex (FreeSASA, as in `q22_residue_features.py`). It is missing at 75 sites.
   - Errors: 1 (a gemmi "no subchain" error).
   - D-mode match (the located pocket covers the ligand or is precise): 10,530 of 11,983 (88%); cropped run: 10,750.
3. **Ligand-contact chemistry (BURROW-independent).** It is computed from the census `contacts` list in two versions:
   - **contact4A_*** uses every protein entry of the stored list. The stored list is the **4 Å set** (all residues with n_simple ≥ 1), **not** the vdW+0.5 Å BioLiP rule that DESIGN §0 names primary (finding of `reviews/verify_misc.md`).
   - **contactBioLiP_*** keeps only residues with ≥2 vdW+0.5 Å atom pairs (`n_vdw ≥ 2`, src/part1/config.py `MIN_ATOM_CONTACTS_PRIMARY`), which is the BioLiP2 rule.
   - Both are reported. Neither depends on the crop.
4. **Compare** (`pdbsite_compare.py`). The code and weights are unchanged between the cropped and uncropped runs.
   - Cliff's δ MBS vs SBS. The pooled contrast is **weighted** by units-available/sampled, so it estimates the population of organic oligomer units rather than the equal-allocation sample. Unweighted, per class and per assembly type are also reported.
   - OLS on the raw value and on its percentile rank: y ~ MBS + B-spline(log ligand heavy atoms, df 4) + class, with cluster-robust SEs. A within-cluster fixed-effects version uses the 682 clusters (4,980 sites) that have both MBS and SBS sites.
   - Pair-pooled within-cluster δ. The top 5 clusters supply 66% of the pairs.
   - On-axis vs off-axis homomer MBS.
   - The Q2.2 contrast (interface pocket ≥3 lining residues on ≥2 chains vs single-chain pocket) applied to the PDB ligand sites.
   - Power gate: all six strata pass, with ≥490 clusters each.

## Results (primary: uncropped, whole assembly 1)

### H2-secondary: burrow_score of the pocket at the ligand
No MBS/SBS difference overall. Slightly lower for MBS within cluster.

| Contrast | Result |
|---|---|
| Medians | 0.320 (MBS) vs 0.343 (SBS) |
| Weighted Cliff's δ | **−0.01 [−0.07, +0.06]** |
| Adjusted β (percentile scale) | +0.003 [−0.034, +0.040], p = 0.86 |
| Within-cluster fixed effects | **−0.030 [−0.053, −0.006], p = 0.014** |
| Within-cluster pair δ | −0.10 [−0.20, −0.01] |
| Per class: cofactor | +0.06 [+0.01, +0.11] |
| Per class: drug-like | −0.03 [−0.21, +0.17] |
| Per class: other | −0.02 [−0.08, +0.03] |
| Restricted to D-matched sites | −0.04 [−0.11, +0.04] |

Reading: MBS ligand sites are not better-scored BURROW pockets. Within a protein cluster they score slightly lower, a small effect driven by a few large clusters.

### EXPLORATORY: is the predicted-dimer pattern reproduced at real sites?
`pdb_vs_afdimer_pattern.tsv`, figure `pdbsite_vs_afdimer_forest.png`. All values are Cliff's δ.

| Metric | AF dimers: interface vs single-chain pocket (Q2.2) | PDB: MBS vs SBS (weighted) | PDB: interface vs single-chain pocket at ligand sites |
|---|---|---|---|
| Volume (merged pocket) | +0.23 | **+0.07 [+0.01, +0.13]** | +0.49 [+0.43, +0.54] |
| Volume (located member) | — | +0.05 [−0.01, +0.11] | +0.36 |
| # lining residues | +0.31 | **+0.08 [+0.02, +0.14]** | +0.50 |
| burrow_score (merged pocket) | +0.09 | +0.02 [−0.03, +0.09] | +0.08 [+0.02, +0.15] |
| Buriedness | +0.10 | −0.02 [−0.09, +0.04] | −0.11 [−0.16, −0.06] |
| Hydrophobic fraction | −0.15 | +0.03 [−0.03, +0.08] | +0.02 [−0.04, +0.08] |
| Aromatic fraction | −0.12 | −0.05 [−0.12, +0.02] | −0.10 [−0.17, −0.04] |
| Charged fraction (DEKR) | +0.20 | −0.01 [−0.07, +0.05] | +0.06 [+0.00, +0.12] |
| Contact residues (4 Å): charged / hydrophobic | — | −0.02 / +0.06 [+0.00, +0.11] | — |
| Contact residues (BioLiP rule): charged / hydrophobic | — | −0.02 [−0.07, +0.03] / +0.04 [−0.01, +0.09] | — |

- **Size: a weak version of the AF pattern.**
  - MBS merged pockets are slightly larger and have more lining residues.
  - Adjusted for ligand size: merged volume β = +0.024 percentile, p = 0.03; lining β = +0.032, p = 0.004.
  - Within cluster: merged volume FE β = +0.040, p = 4e-4.
  - The effect is about a third of the AF effect, and smaller than in the cropped run (+0.11).
- **Chemistry: not reproduced.** Pooled charged δ is −0.01 and hydrophobic δ is +0.03, against +0.20 and −0.15 in the AF dimers. The BURROW lining and both contact definitions agree on this.
- **Small charged signals in subsets:**
  - within-cluster FE, charged: β = +0.028, p = 0.009; pair δ +0.05 [−0.03, +0.14];
  - homomer stratum: +0.08 [+0.02, +0.14];
  - cofactor stratum: +0.06 [+0.02, +0.11].

  Heteromer sites go the other way: −0.08 [−0.16, −0.01] (BioLiP contacts −0.12 [−0.17, −0.06]).
- **Q2.2 definition applied to PDB ligand sites:** interface pockets are much larger, less aromatic and slightly more charged (+0.06, lower CI at 0.00). They are less buried, which is the opposite of the AF result.
  - So ligand-occupied interface pockets carry only a faint version of the charged/polar signature of predicted interface pockets.
  - One possible reason: occupancy selects the less polar subset. This is an interpretation and has not been tested.
- **Sanity check** (`pocket_chain_composition.tsv`): 69–74% of MBS sites sit in a merged pocket with ≥3 lining residues on ≥2 chains, against 15–22% of SBS sites.

### EXPLORATORY: on-axis vs off-axis homomer MBS
`onaxis_vs_offaxis_mbs.tsv`, figure `pdbsite_onaxis_vs_offaxis.png`. The groups are 517 on-axis sites (125 clusters) and 3,509 off-axis sites.

- No burrow_score or size difference: burrow_score δ +0.01 [−0.37, +0.35].
- On-axis lining is much less aromatic: δ −0.46 [−0.63, −0.03].
- It is more hydrophobic: δ +0.23 [−0.03, +0.42]; 4 Å contacts +0.24 [+0.09, +0.35].
- **Concentration:** 47% of on-axis sites come from two clusters, 93 (rep. 1ZPA) and 259 (9JIS). These are HIV-protease-like C2 sites (not checked per entry), so this contrast largely describes those two families.

## Sensitivity: 30 Å crop vs whole assembly
`crop_vs_uncropped_sites.tsv` (per site) and `crop_vs_uncropped_deltas.tsv` (pooled).

**Per site.**
- The crop removed chains at 2,500 sites: 1,352 of 5,532 MBS (24.4%) and 1,148 of 6,451 SBS (17.8%). So it acted differentially across the contrast, as the review found.
- At the 9,483 sites where it removed nothing, all outputs are identical.

| Sites where chains were cropped | SBS (n = 1,148) | MBS (n = 1,352) |
|---|---|---|
| Spearman burrow_score, cropped vs uncropped | 0.47 | 0.40 |
| burrow_score changed by >0.05 | 60% | 61% |
| Median abs. Δ burrow_score | 0.071 | 0.083 |
| Located-pocket volume unchanged | 3.7% | 2.8% |
| D-match lost / gained when uncropped | 9.7% / 3.0% | 15.2% / 4.6% |
| Median Jaccard of located-pocket lining | 0.40 | 0.33 |
| Spearman merged charged fraction | 0.54 | 0.49 |

**Pooled MBS vs SBS (weighted δ, seed 0).**

| Metric | Uncropped (primary) | Cropped 30 Å |
|---|---|---|
| burrow_score (H2-sec.) | −0.009 [−0.071, +0.056] | −0.031 [−0.093, +0.038] |
| Volume (located member) | +0.046 [−0.015, +0.109] | +0.078 [+0.020, +0.141] |
| Merged volume | +0.071 [+0.010, +0.132] | +0.107 [+0.048, +0.166] |
| # lining (merged) | +0.081 [+0.020, +0.142] | +0.115 [+0.060, +0.173] |
| burrow_score (merged) | +0.025 | +0.019 |
| Buriedness | −0.023 | −0.028 |
| Hydrophobic fraction | +0.027 [−0.030, +0.078] | +0.033 [−0.023, +0.085] |
| Aromatic fraction | −0.052 | −0.045 |
| Charged fraction | −0.009 [−0.068, +0.047] | −0.000 [−0.059, +0.057] |

The crop does not change the conclusions. It did inflate the MBS size effect by about 0.03–0.04 δ, because cropped MBS sites lost partner-chain context more often than SBS sites. The burrow_score H2 secondary is null in both runs. The within-cluster burrow_score effect reaches p < 0.05 only uncropped: FE p = 0.014, against 0.058 cropped.

## Outputs
- Primary tables in `results/part1/pocket_descriptors/`:
  - `sample_strata_n12000_seed0.tsv`, `run_ledger.tsv`, `power_gate.tsv`;
  - `mbs_vs_sbs_cliffs_delta.tsv`, `mbs_vs_sbs_regression.tsv`, `mbs_vs_sbs_within_cluster.tsv`;
  - `onaxis_vs_offaxis_mbs.tsv`, `pocket_chain_composition.tsv`, `pdb_vs_afdimer_pattern.tsv`;
  - `crop_vs_uncropped_{sites,deltas}.tsv`.
- Sensitivity: the same tables in `results/part1/pocket_descriptors/crop30/`.
- Figures (primary) in `results/part1/figures/`: `pdbsite_vs_afdimer_forest.png`, `pdbsite_mbs_vs_sbs_by_stratum.png`, `pdbsite_burrow_score_ecdf.png`, `pdbsite_onaxis_vs_offaxis.png`.
- `data/part2/pdb_site_fingerprints.parquet` (uncropped; not in git because data/ is ignored):
  - 11,983 sites with a 384-d fingerprint plus labels.
  - **The fingerprint is of the located BURROW member, not the merged pocket.** The atlas and Q2.2 fingerprints are of merged pockets. In 69% of SBS and 65% of MBS sites the merged pocket has a single member, where the two are identical.
- Per-site tables: `data/interim/pdbsite/sites.parquet` (uncropped) and `sites_crop30.parquet`. Raw describe parts: `data/interim/pdbsite/describe_uncropped/` and `describe_crop30/`.

## Caveats and choices (ours)
- **Unit definitions differ between the predicted and PDB contrasts.**
  - Q2.2 compares all atlas-pass pockets split by geometry.
  - Here every pocket is selected by a bound ligand, and the main PDB contrast is by ligand contacts (MBS).
  - The "interface vs single-chain pocket" column is the closest like-for-like comparison, but it is still conditioned on occupancy.
- **No atlas-style burrow_score filter** is applied here, because the site is given by the ligand.
- **Regression models are unweighted.** The equal-allocation sample is adjusted for class, but not reweighted.
- **Covariate:** heavy-atom count is that of the sampled instance. Resolution and method are not in the models.
- **BH:** the H2 secondary is a single descriptor, so no BH is applied. Family A H1–H3 are already tested in docs/PART1_RESULTS.md.
- **The pair-pooled within-cluster δ is dominated by a few large clusters.** Prefer the FE regression.
- **Bug fixed in v1:** an early compare run had the "all_unweighted" contrast accidentally weighted. All reported numbers come from the corrected code.
