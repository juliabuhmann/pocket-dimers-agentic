# fpocket druggability at MBS vs SBS ligand sites and at predicted dimer-only pockets

**Status: UNVERIFIED. EXPLORATORY (not pre-registered).** Added after REPORT v2 to fill the gap noted in REPORT.md §4.3
("there is no druggability predictor"). Nothing here is part of the frozen DESIGN.md hypotheses.

## Analysis plan (declared before any MBS/SBS comparison was computed)

Written and committed after the fpocket install, the 1-structure test and the 50-site timing pilot, and before the
comparison script was run.

- **Primary metric for this exploratory analysis:** fpocket *Druggability Score* (`fp_drug`, 0–1) of the fpocket
  pocket that best overlaps the bound ligand.
- **Primary contrast:** MBS vs SBS, weighted Cliff's δ (weights = units available / units run per stratum, so the
  estimate refers to the population of organic oligomer units, as in docs/PDB_SITE_DESCRIPTORS.md), cluster30
  bootstrap with 1,000 resamples, seed 0.
- **Secondary (all exploratory):**
  - fraction of sites with `fp_drug > 0.5` (fpocket's conventional threshold, Schmidtke & Barril 2010), MBS vs SBS;
  - OLS `fp_drug ~ MBS + B-spline(log ligand heavy atoms, df 4) + class`, on the raw value and percentile rank,
    cluster-robust SEs; within-cluster fixed effects on clusters that have both MBS and SBS sites; pooled
    within-cluster pair δ;
  - by ligand class (drug-like / cofactor / other) and by assembly type;
  - on-axis vs off-axis homomer MBS;
  - other fpocket descriptors: fpocket score, volume, hydrophobicity score, polarity score;
  - sensitivity: sites where the chosen fpocket pocket covers ≥50% of the ligand heavy atoms.
- **Predicted dimers:** fpocket druggability of the curated dimer-only interface pockets (results/part2/q23a_curated.tsv)
  vs a matched random sample of single-chain pockets from the same AF dimer models; descriptive only.
  *Disclosure:* a first AF-dimer run was made before this plan was written. Its pocket mapping (max residue overlap
  coefficient) and its single-chain matching (only 19 matches) were judged inadequate and are revised; see Deviations.
- No multiple-testing correction: this is a single exploratory family and is reported as such.

<!-- RESULTS BELOW ARE ADDED AFTER THE RUN -->

## Rerun
```bash
conda create -p conda_envs/fpocket -c conda-forge -c bioconda fpocket      # fpocket 4.2.3 (build h5fd1fdb_2)
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/fpocket_pdbsite.py --limit 50 --workers 10 --out data/interim/fpocket/pilot50   # timing pilot
.venvs/part1/bin/python scripts/fpocket_subsample.py                    # -> results/part1/fpocket/pdbsite_subsample_seed0.tsv
.venvs/part1/bin/python scripts/fpocket_pdbsite.py --subsample results/part1/fpocket/pdbsite_subsample_seed0.tsv --workers 9   # resumable
.venvs/part1/bin/python scripts/fpocket_compare.py                      # -> results/part1/fpocket/pdbsite_*.tsv
.venvs/part1/bin/python scripts/fpocket_afdimer.py --workers 1          # -> results/part1/fpocket/afdimer_*.tsv
.venvs/part1/bin/python scripts/fpocket_figures.py                      # -> results/part1/fpocket/figures/*.png
```

## Method

### Tool
- fpocket 4.2.3 from bioconda (`conda_envs/fpocket`), default parameters (`fpocket -f x.pdb`).
- Tested on one structure (site 5b0w_05103, 3 chains, 46 pockets, ~2.4 s) before the pilot.
- Descriptors kept per pocket (from `<name>_info.txt`): Score, Druggability Score, Volume, Hydrophobicity score,
  Polarity score, plus SASA, charge score, alpha-sphere counts (full list in `scripts/fpocket_common.py`).
- The druggability score is fpocket's built-in logistic model trained on the DCC druggable/non-druggable set
  (Schmidtke & Barril 2010). 0.5 is its conventional cut-off.

### PDB sites (`scripts/fpocket_pdbsite.py`)
- Same sites as docs/PDB_SITE_DESCRIPTORS.md: the 11,983 described sites, read from
  `data/interim/pdbsite/sites_crop30.parquet` (same site_id, ligand instance, labels and weights; read-only).
- Structure: gemmi, first altloc, no H, assembly 1 (first listed), **whole assembly, no crop**.
- The ligand instance is re-found with the same rule as `pdbsite_describe.py` (CCD components nearest the census
  centroid).
- Only polymer residues are written (gemmi entity type Polymer, so MSE and other modified residues stay). All
  ligands, waters and additives are removed. Chains are renamed to single characters for PDB format.
- Pocket choice: the fpocket pocket with the most ligand heavy atoms within 3 Å of its alpha-sphere centres
  (ties: better fpocket rank). The count and the fraction of ligand atoms covered are stored
  (`fp_overlap_n`, `fp_overlap_frac`).
- Per-site fpocket timeout 900 s (our choice; pilot max 331 s). Timeouts are logged as errors.

### Subsample (deviation from "all 11,983")
- Pilot: 50 sites (seed-0 shuffle), all ok, 13.7 CPU-s/site mean (median 2.5 s; max 331 s for a 66,526-atom assembly).
  All 11,983 sites would take ~4.6 h on 10 workers, beyond the ~1.5 h budget.
- So `scripts/fpocket_subsample.py` draws 500 sites per ligand class × MBS/SBS stratum (seed 0): 3,000 sites,
  1,752 clusters. Pooled contrasts are re-weighted by units available / sites analysed per stratum
  (`results/part1/pocket_descriptors/sample_strata_n12000_seed0.tsv`).
- Equal allocation means MBS drug-like (1,815 available) and MBS cofactor (1,577) are sampled at a higher rate than
  the others; per-class contrasts are unweighted.

### Statistics (`scripts/fpocket_compare.py`)
- As declared above; the code mirrors `scripts/pdbsite_compare.py` (helper functions copied).
- Analysis set: status ok (fpocket ran and at least one pocket has ≥1 ligand atom within 3 Å).

### Predicted dimers (`scripts/fpocket_afdimer.py`)
- Models: the 112 AF dimer models (79 homo / 35 hetero rows) of the 114 curated dimer-only pockets in
  `results/part2/q23a_curated.tsv`, from `/path/to/afdb_dimers/afdb_human_dimers` (as `q22_common.py`).
- fpocket is run once on each full dimer model.
- BURROW pockets come from `data/part2/human_dimer_pockets_annotated.parquet`. Each is mapped to the fpocket pocket
  with the highest residue Jaccard (chain:resnum sets). "Mapped" = Jaccard ≥ 0.2 (our choice).
  Best-case sensitivity: max druggability over all fpocket pockets with overlap coefficient ≥ 0.5.
- Controls, drawn at random (seed 0) without replacement from the same model as each dimer-only pocket:
  - `single_chain_matched` (primary): any atlas-pass, symmetry-representative single-chain pocket;
  - `single_chain_strict`: additionally volume 300–1,500 Å³ and burrow_score ≥ the curated list's threshold (0.1127,
    read from the TSV header), the same pocket filters as the curated list.

## Results (UNVERIFIED: one run, not independently re-derived)

### Run ledger (`pdbsite_run_ledger.tsv`)
- 3,000 subsampled sites (1,752 clusters). 2,901 ok (1,690 clusters).
- 81 sites had no fpocket pocket within 3 Å of any ligand atom (`no_overlap`).
- 18 errors: 11 assemblies over the PDB-format atom limit (116k–178k atoms) and 7 fpocket timeouts (>900 s).
  **All 18 are very large assemblies, so they are excluded non-randomly.**
- Chosen pocket: median 96% of ligand heavy atoms within 3 Å. 2,473 sites (85%) are ≥50% covered. The chosen
  pocket is fpocket's top-ranked pocket in 15% of sites (median rank 12). So fpocket's own ranking often does not put
  the ligand pocket first in whole assemblies.

### Primary (exploratory): fpocket druggability, MBS vs SBS: no difference
`pdbsite_mbs_vs_sbs_cliffs_delta.tsv`, `pdbsite_mbs_vs_sbs_regression.tsv`, figures `fpocket_pdbsite_druggability_ecdf.png`, `fpocket_pdbsite_forest.png`.

| Contrast | Result |
|---|---|
| Medians | 0.33 (MBS) vs 0.35 (SBS) |
| **Weighted Cliff's δ** | **+0.01 [−0.08, +0.09]** |
| Unweighted δ | −0.04 [−0.11, +0.03] |
| Pocket covers ≥50% of ligand (weighted) | +0.02 [−0.06, +0.11] |
| Adjusted β (raw / percentile) | −0.019 [−0.065, +0.026] / −0.016 [−0.050, +0.017] |
| Within-cluster FE (116 clusters, 588 sites) | +0.039 [−0.049, +0.126] |
| Within-cluster pair δ | +0.04 [−0.13, +0.34] (top 5 clusters supply 72% of pairs) |
| **Fraction > 0.5 (weighted)** | **40.3% MBS vs 39.3% SBS**, difference +1.0 pp [−5.8, +7.1] |

By stratum (δ; fraction >0.5 MBS vs SBS):
- drug-like: −0.13 [−0.26, +0.01]; 41.5% vs 50.5%. MBS drug-like sites lean *less* druggable, but the CI touches 0;
- cofactor: −0.04 [−0.11, +0.05]; 43% vs 44%;
- other organic: +0.04 [−0.06, +0.13]; 40% vs 36%;
- homomer: +0.05 [−0.05, +0.15]; heteromer: −0.10 [−0.24, +0.05].

### Other fpocket descriptors (weighted δ, MBS vs SBS)
- fpocket volume: +0.08 [+0.01, +0.14]. MBS pockets are slightly larger, as with BURROW (docs/PDB_SITE_DESCRIPTORS.md).
  Adjusted for ligand size: +143 Å³ [+49, +237] raw, but percentile β +0.01 (n.s.); within-cluster pair δ +0.25 [+0.11, +0.41].
- fpocket score: −0.09 [−0.16, −0.02]. MBS pockets have a slightly lower fpocket *ranking* score.
- hydrophobicity score +0.02 [−0.06, +0.09]; polarity score 0.00 [−0.08, +0.08].
- These secondary descriptors are not corrected for multiplicity.

### On-axis vs off-axis homomer MBS (`pdbsite_onaxis_vs_offaxis_mbs.tsv`)
- 135 on-axis (53 clusters) vs 911 off-axis sites. Druggability δ −0.03 [−0.19, +0.25]; >0.5: 38% vs 43%.
- On-axis pockets have a lower polarity score: δ −0.48 [−0.64, −0.28]. This is consistent with the lower aromatic /
  higher hydrophobic lining reported for on-axis sites from BURROW.
- 61 of the 135 on-axis sites come from clusters 93 and 259 (the HIV-protease-like families noted in
  docs/PDB_SITE_DESCRIPTORS.md), so this mostly describes those two families.

### Predicted dimers (`afdimer_summary.tsv`, `afdimer_paired.tsv`, figure `fpocket_afdimer_druggability.png`)
- fpocket ran on all 112 models. 114 of 114 dimer-only pockets overlap an fpocket pocket; 107 map with Jaccard ≥0.2
  (median Jaccard 0.36).
- **Dimer-only interface pockets: median druggability 0.010 (IQR 0.001–0.085); 9.6% > 0.5.**
- Matched single-chain pockets, same models (96 with an fpocket hit, 94 models): median 0.006; 8.3% > 0.5.
  - Paired per model: dimer-only higher in 45% of 94 models; Wilcoxon p = 0.69.
  - These controls are much smaller BURROW pockets (median 97.5 vs 465 Å³), so they are not size-matched.
- Strict single-chain control (same pocket filters as the curated list): only 19 models have one. Median 0.005;
  21% > 0.5; paired p = 0.07 (dimer-only lower).
- Best-case mapping (max druggability over fpocket pockets with overlap coefficient ≥0.5): dimer-only median 0.020,
  10.7% > 0.5.
- **Reading:** by fpocket, the curated dimer-only pockets are no more druggable than single-chain pockets in the same
  models. Only ~10% pass the conventional 0.5 cut-off, against ~40% of real ligand-bound PDB sites (whose pocket was
  chosen by the ligand, so the two numbers are not like-for-like).

## Interpretation (ours, provisional)
- fpocket druggability does not separate MBS from SBS ligand sites, which matches the pre-registered burrow_score
  null (H2 secondary). The modest size excess of MBS pockets is seen again with an independent pocket finder.
- The q23a dimer-only list is not enriched for fpocket-druggable pockets. This weakens any reading of that list as
  "druggable" beyond BURROW's own score.

## Deviations and caveats
- **Subsample:** 3,000 of 11,983 sites (500 per stratum, seed 0), forced by run time (see Method). The within-cluster
  analyses are therefore thin (116 clusters with both MBS and SBS).
- **Large assemblies excluded:** 18 sites (PDB-format limit or timeout). An mmCIF input path or a crop would recover
  them; neither was done.
- **Whole-assembly fpocket** differs from the BURROW run here: the current PDB-site BURROW numbers are cropped
  (uncropped rerun in progress by another agent). fpocket was not cropped.
- **Polymer filter differs slightly** from `pdbsite_describe.py`, which keeps non-HETATM residues: here modified
  residues (MSE etc.) are kept as polymer.
- **Pocket choice by ligand overlap** guarantees a ligand-bound pocket for PDB sites. For AF dimers the choice is by
  residue overlap with a BURROW pocket. The two selection rules differ.
- **AF-dimer mapping and controls were revised after a first look** (disclosed above): the first run used the max
  overlap coefficient (which picked small sub-pockets inside the BURROW pocket) and only the strict control (19
  matches). The druggability medians were low (≤0.02) under every mapping tried.
- fpocket's druggability model was trained on apo/holo crystal pockets, not on AF models. Its calibration on predicted
  interface pockets is unknown.
- Statistics: same estimators as docs/PDB_SITE_DESCRIPTORS.md; regression is unweighted.

## Outputs
- `results/part1/fpocket/`: `pdbsite_subsample_seed0.tsv`, `pdbsite_run_ledger.tsv`, `pdbsite_mbs_vs_sbs_cliffs_delta.tsv`,
  `pdbsite_mbs_vs_sbs_regression.tsv`, `pdbsite_onaxis_vs_offaxis_mbs.tsv`, `afdimer_pockets.tsv`, `afdimer_summary.tsv`,
  `afdimer_paired.tsv`.
- Figures: `results/part1/fpocket/figures/fpocket_pdbsite_druggability_ecdf.png`, `fpocket_pdbsite_forest.png`,
  `fpocket_afdimer_druggability.png`.
- Per-site data (not in git): `data/interim/fpocket/pdbsite_sites.parquet`, `data/interim/fpocket/pdbsite/part_*.parquet`,
  `data/interim/fpocket/afdimer_pockets.parquet`, `data/interim/fpocket/afdimer_fpocket_raw.pkl`.
- Logs: `logs/fpocket/{pilot50,pdbsite_sub3000,afdimer,compare}.log`.
