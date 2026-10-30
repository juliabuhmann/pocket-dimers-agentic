# H4 / Q2.1: do predicted human dimers recover PDB ligand sites, and multi-chain sites better than isolated chains?

**Status: UNVERIFIED.** This is a first run, done 2026-10-29. Nobody has re-derived the headline numbers independently (DESIGN §5 requires this). Only one site was checked by hand (the sentinel below). Every number here comes from the TSVs and is not a confirmed result.

## Provenance

| Item | Value |
|---|---|
| DESIGN.md sha256 | `e3f9b550eb29cfe6eb444ebd884920b847687518470834c0bcee3439998eee28` (equals `results/DESIGN_FROZEN.sha256`) |
| Census | **v2** (`data/part1/ligand_sites.parquet`, `census_version == "v2"`) |
| Scripts | `scripts/h4_*.py`, committed in f74eec9; the figure-label fix was committed with this document |
| Git HEAD stamped in the outputs | 41b1b2d / 55df208. Other agents committed to the branch during the run, so the stamp is HEAD at run time. |

### Rerun

```bash
source scripts/env_threads.sh
P=.venvs/part1/bin/python
$P scripts/h4_select.py --census v2 --workers 6      # step 1: candidates, SIFTS XML fetch (<=4 threads), assembly, interface, mapping
$P scripts/h4_feasibility.py --census v2             # feasibility gate (run BEFORE recovery) -> feasibility.tsv, gate.json
$P scripts/h4_recovery.py --census v2 --workers 6    # steps 2, 3, 5: superposition, DCC, residue metrics, null, shuffled control, atlas
$P scripts/h4_sentinels.py --census v2               # step 6: hand-check sentinel + negative control (exit 1 on failure)
$P scripts/h4_stats.py --census v2                   # step 4: H4 tests, strata, sensitivities
$P scripts/h4_figures.py --census v2
```

The pipeline was built and tested on census v1, on a 30-entry subset only. The v1 outputs were test runs and were deleted. Every number below is from v2.

- Outputs: `results/part2/h4/v2/*.tsv` and `gate.json`. Each TSV header carries the design sha256 and the git commit.
- Figures: `results/part2/figures/h4_{hit_rates,did_forest,dcc_ecdf,residue_f1}_v2.png`.
- Per-site data: `data/part2/h4/v2/` (git-ignored).

## Methods as implemented, with the choices DESIGN leaves open

### Benchmark sites
- **Sites.** Census v2 `in_primary_set` sites: BioLiP2-covered, SoI arm, ≥2-atom-contact binding residues, not duplicates, not glycans, not additives.
- **Assemblies.** Only `homomer` and `heteromer` assemblies are used. `same_uniprot_hetero` is therefore excluded.
- **Ligand sets.** The primary set is organic ligands (`ligand_class` not in metal_ion, metal_cluster, inorganic). Metal ions and clusters are reported separately.

### Pair assignment
- **MBS.** The two chains with the most binding residues must map to the UniProt pair of a dimer model.
- **SBS.** The binding chain's protein must be in a dimer model of a pair present in the entry. The partner is the chain of the partner accession with the largest interface (5 Å) to the binding chain. "Partner accession" means the same accession for homodimers and the model partner for heterodimers.

### Residue mapping and redundancy
- **Mapping.** Per-residue SIFTS XML, keyed by (auth chain, label_seq). The XML was fetched only for the 5,462 selected entries, with 0 failures.
- **Filter.** Both chains must map to the model accessions, and the binding residues must have ≥90% residue-name identity with UniProt.
- **Non-redundant unit.** Model × cluster30 × CCD × MBS label. Symmetric copies are collapsed first by the census `site_group` within each entry.
- **Representative.** For each unit, the best-resolution entry is used (ties go to the earliest ID), chosen from up to 3 candidate entries.

### Interface filter
- PDB interface residues (≤5 Å heavy atom) are compared with model interface residues (≤5 Å), in UniProt numbering and aware of chain role.
- The score is a Jaccard index. For homodimers the maximum over the two chain assignments is used.
- A site is kept if Jaccard ≥0.3.

### Superposition (choice)
- **Local frame.** Mapped CA pairs of both chains within 15 Å of the ligand centroid are used. If there are fewer than 10, the whole chain pair is used instead; this happened for 1 of 4,013 sites.
- **Fit quality.** Median local CA RMSD is 0.51 Å.
- **Shared transform.** The same transform serves the dimer and the isolated chains, because the isolated chains have the same coordinates.

### Hit and top-k
- A hit is DCC ≤4 Å to any of the top-(n+2) pockets by burrow_score.
- n = the number of collapsed site groups of the same ligand set on that chain pair in that entry.
- Isolated chains: the pooled pockets of both chains.
- Homodimers: both chain assignments are scored and the maximum is taken.

### Homodimer pocket symmetric collapse (deviation, see caveat 1)
- Pockets are copies if Jaccard(swap_AB(residues), residues) ≥0.5. Pockets are merged by union-find, and the highest-scoring pocket in each group is kept.
- The same rule is applied to the dimer pockets and to the pooled isolated pockets.
- Sensitivities: no collapse, and the overlap-coefficient rule from `annotate_pocket_chains.py`.

### Other metrics and controls
- **Residue metrics.** Recall, overlap coefficient and F1 are taken from the best-F1 pocket among the top-(n+2), and the size of that pocket is reported.
- **Null.** For each top-k pocket, a random surface patch of the same residue count. The seed is a random residue with fewer than 20 CB neighbours within 10 Å; the patch is the seed plus its nearest residues, kept on the same chain for isolated chains. The patch centre is the mean CB. Seeds 0–2 are used.
- **Negative control.** The transformed ligand centroid is randomly rotated about the model's centre of mass (seed 0).

### Statistics
- **Contrast.** DiD of the paired per-site gain (hit_dimer − hit_iso), MBS minus SBS.
- **Bootstrap.** 1,000 resamples, seed 0, stratified by homo/hetero. The resampling unit is the connected component of models and cluster30 IDs, which gives 470 clusters in the primary set.
- **Permutation.** The MBS label is permuted within clusters, 10,000 times. The p-value is centred on the permutation mean because the pooled statistic is not centred at 0 when cluster composition differs.
- **Raw MBS gain.** Bootstrap CI and a cluster sign-flip p.
- **Multiple testing.** Family-B BH correction has to wait for the H5 p-value; the p-values here are raw.

### Strata
- **AF training cutoff.** **2021-09-30**. LITERATURE.md gives this date only for the post-cutoff split of AF3 and Runs N' Poses. That the AFDB-complex (AF-Multimer) models used the same PDB cutoff is an assumption [UNVERIFIED].
- **Post-cutoff pair.** No pre-cutoff entry contains the pair as a complex:
  - homodimers: no pre-cutoff entry has ≥2 chains of the accession in its asymmetric unit, and no pre-cutoff census homomer assembly contains it;
  - heterodimers: no pre-cutoff entry has both accessions in the same asymmetric unit.
- **ipSAE tertiles.** Taken over the models in the primary set; the cuts are 0.802 and 0.878.

## Counts per step (v2)

| Step | n |
|---|---|
| Primary-set homomer/heteromer sites in entries containing a model protein | 44,370 (10,800 entries) |
| Dropped: SBS without a model pair in the entry / MBS pair not modelled / chain without UniProt | 13,424 / 2,495 / 242 |
| Site × model candidate rows | 29,436 (6,592 entries) |
| Non-redundant units before mapping | 5,697 |
| Entries processed (≤3 per unit) | 5,462 (0 errors) |
| Chain pair resolved in assembly 1 | 24,174 (237 rows had no partner chain) |
| SIFTS accession match / ≥90% residue identity | 24,168 / 23,851 |
| After symmetric-copy collapse | 11,601 |
| Superposition ok and model has pockets | 11,597 |

### Feasibility gate

The gate was logged before recovery was computed (`feasibility.tsv`, `gate.json`). Counts are non-redundant organic units after mapping; numbers in brackets are bootstrap clusters.

| | MBS | SBS |
|---|---|---|
| Homodimers | 1,346 (172) | 2,417 (317) |
| Heterodimers | 135 (63) | 668 (123) |

- Heterodimers have 135 MBS units, at least 100, so they **pass** the gate and are included in the primary test.
- Metals: homo MBS 42 / SBS 356 units; hetero MBS 19 / SBS 295.

### Interface filter (`interface_filter.tsv`, units removed by Jaccard <0.3)

| Organic units | Before | Removed |
|---|---|---|
| Homo MBS | 1,346 | **327** |
| Homo SBS | 2,417 | 157 |
| Hetero MBS | 135 | 16 |
| Hetero SBS | 668 | 53 |

The primary set after filtering is 4,013 organic sites: 1,138 MBS and 2,875 SBS, from 607 models and 3,509 entries.

## Results

### Hit rates

Organic ligands, interface-filtered, DCC ≤4 Å within the top-(n+2) pockets (`hit_rates.tsv`, with cluster-bootstrap CIs).

| | n | Dimer | Isolated | Random patches (seed 0) | Shuffled centroid (dimer) |
|---|---|---|---|---|---|
| Homo MBS | 1,019 | **0.401** | 0.278 | 0.005 | 0.045 |
| Homo SBS | 2,260 | 0.432 | 0.469 | 0.005 | 0.025 |
| Hetero MBS | 119 | **0.227** | 0.134 | 0.000 | 0.034 |
| Hetero SBS | 615 | 0.283 | 0.350 | 0.007 | 0.020 |
| Metal, homo MBS / SBS | 33 / 323 | 0.121 / 0.176 | 0.030 / 0.176 | | |
| Metal, hetero MBS / SBS | 18 / 275 | 0.111 / 0.178 | 0.389 / 0.215 | | |

### H4, primary contrast

Homo + hetero, organic, interface-filtered (`h4_tests.tsv`, row 1).

| Quantity | Estimate | 95% CI / p |
|---|---|---|
| Dimer gain at MBS (dimer − isolated hit rate) | **+0.120** (0.383 vs 0.263) | CI −0.048 to 0.274; sign-flip p = 0.17 |
| Dimer gain at SBS | −0.043 (0.400 vs 0.443) | CI −0.119 to 0.025; sign-flip p = 0.24 |
| **DiD (MBS gain − SBS gain)** | **+0.164** | cluster-bootstrap CI **−0.015 to 0.332**; bootstrap p = 0.066; within-cluster permutation p = 1e-4 (centred) |

- Discordant pairs at MBS: 310 sites hit only in the dimer and 173 only in the isolated chains. At SBS: 369 and 494.
- Power gate passed: 206 MBS clusters and 412 SBS clusters.
- Seeds 1 and 2 give DiD CIs of (−0.004, 0.330) and (−0.006, 0.336).

**How to read this.** The point estimate is in the hypothesised direction. Dimer models recover multi-chain sites better than isolated chains (+12 points), while single-chain sites are not helped and are slightly worse (−4 points). The two pre-registered inference routes disagree:
- The cluster bootstrap, which treats protein families as the sampling unit, gives a 95% CI that just includes 0.
- The within-cluster label permutation, which holds the observed families fixed, is highly significant.

The difference comes from large between-family heterogeneity in the gain. For example, NAMPT (P43490) has an MBS gain of −0.88 over 56 inhibitor sites in one family. At α = 0.05 under the family-level bootstrap, **H4 is not established**. The result is directionally supported and robust in sign across every sensitivity analysis and stratum.

### Sensitivities and strata (DiD, 95% bootstrap CI)

| Analysis | DiD | 95% CI |
|---|---|---|
| Homo only (n = 3,279) | 0.161 | −0.035 to 0.356 |
| Hetero only (n = 734) | 0.159 | −0.053 to 0.357 |
| Unfiltered: no interface filter (n = 4,566) | 0.136 | −0.006 to 0.292 |
| No symmetric pocket collapse | 0.104 | −0.057 to 0.254 |
| Overlap-coefficient collapse (annotation rule) | 0.175 | 0.014 to 0.326 |
| Atlas pocket filter | 0.164 | −0.015 to 0.333 |
| *EXPLORATORY*: one site per binding region per model (n = 1,708) | 0.145 | 0.071 to 0.214 |
| Post-cutoff pairs (n = 224; 84 models) | 0.194 | 0.030 to 0.349 |
| Pre-cutoff pairs (n = 3,789) | 0.162 | −0.019 to 0.343 |
| ipSAE tertile 1 / 2 / 3 | 0.246 / 0.106 / 0.153 | All three CIs include 0 |
| Metals (descriptive; n = 649, 51 MBS) | −0.022 | −0.182 to 0.174 |
| Shuffled-centroid control | 0.019 | −0.001 to 0.036 |
| Random-patch null, seed 0 | 0.001 | |

- In the unfiltered set, the MBS hit rate falls from 0.383 to 0.305 in the dimer and from 0.263 to 0.205 in the isolated chains. Sites whose PDB interface the model does not reproduce are recovered less often by both.
- The region-deduplicated sensitivity is the only one whose CI clearly excludes 0. It is exploratory and was not pre-registered: it removes the many-ligands-at-one-site redundancy that cluster30 × CCD keeps.

### Secondary: residue overlap and the AFDB monomer atlas (step 5)

Mean F1 of the best pocket among the top-(n+2), organic, filtered (`residue_metrics.tsv`). The atlas columns use sites whose proteins have tier-1 atlas pockets: homo 1,019 / 2,255, hetero 116 / 614.

| | Dimer | Isolated | AFDB monomer atlas | Random patches |
|---|---|---|---|---|
| Homo MBS | **0.371** | 0.311 | 0.296 | 0.084 |
| Homo SBS | 0.375 | 0.391 | 0.372 | 0.063 |
| Hetero MBS | **0.265** | 0.200 | 0.196 | 0.033 |
| Hetero SBS | 0.342 | 0.336 | 0.322 | 0.037 |

Median residue recall at MBS:
- Homodimers: dimer 0.67, isolated 0.40, atlas 0.40.
- Heterodimers: dimer 0.53, isolated 0.27, atlas 0.28.

Dimer pockets at MBS are larger (median 17 vs 12 residues in homodimers).

The atlas comparison uses residues only. No atlas monomer structures are available locally, so atlas DCC was not computed. Note also that the atlas stores only filtered ("kept") pockets, about 4.3 per protein. Isolated chains cut from the dimer and the AFDB monomer models perform similarly at MBS.

### Known-answer sentinels (`sentinels.tsv`, all PASS)

**Hand check: GSTP1 (P09211) + glutathione, PDB 11gs.**
- This is an MBS site in a C2 homodimer, interface Jaccard 0.97.
- Pipeline: dimer DCC 2.80 Å, a hit (third-ranked pocket of the top 3); isolated-chain DCC 5.07 Å, a miss.
- Independent recomputation: plain gemmi, global CA superposition (RMSD 0.39 Å over 416 CA) gives DCC 3.02 Å to the third-ranked pocket. This agrees with the pipeline's local fit (RMSD 0.37 Å).

**Negative control: shuffled centroid.**
- Hit rate: dimer 0.030, isolated 0.027.
- Random-patch null: 0.005 / 0.004.
- Real hit rates: 0.395 / 0.392.
- This is within the pre-set margin of +0.05 over the null, but it is not "≈ null". A rotation about the centre of mass keeps the centroid's distance from the centre, so interface-proximal centroids stay near the central interface pockets.

**Order of checks.** The sentinels run after the recovery script, because they need its output, and they gate the statistics. They did not abort a pipeline start as DESIGN §5 describes.

### Not done
- **Step 7, LIGYSIS calibration (optional).** Not done. The published 60.4% uses ligand-centroid DCC in monomer structures, and no AFDB monomer coordinates are available locally; a residue-only comparison would not be like-for-like. This is the first item on the DESIGN §7 cut list.
- **Other required items still open:** the independent re-derivation of the headline numbers and the PI spot-check of ≥10 sites per class.

## Caveats
1. **Symmetric pocket collapse rule (deviation, decided during analysis).**
   - `annotate_pocket_chains.py` (`is_sym_rep`) uses the overlap coefficient ≥0.5 with union-find. This lets a large interface pocket absorb distinct single-chain pockets: in the GSTP1 dimer, the 58-residue m0 absorbed the chain-A G-site pockets m1 and m2. This changed H4 test-subset hits.
   - H4 therefore uses swapped-residue **Jaccard** ≥0.5, the rule from the burrow dimer benchmark (TOOLS.md §2), and reports both other rules.
   - The annotation's `is_sym_rep` / `dimer_only` columns used by other analyses may be affected. The PI should review this.
2. **DCC ≤4 Å to burrow pocket centres is strict.** Pocket centres are cavity centres, and large merged pockets have centres far from the ligand. This penalises dimer interface pockets, which are larger. Hit rates of about 0.4 are therefore not comparable with LIGYSIS-style benchmarks.
3. **Redundancy.** The pre-registered unit, cluster30 × CCD, lets one binding site with many ligands (inhibitor series) count many times. This inflates between-family variance. See the exploratory region-dedup row above.
4. **Training overlap.** 3,789 of 4,013 primary sites are on pairs with a pre-cutoff complex in the PDB. The post-cutoff stratum is small (n = 224, 84 models), and its pair definition is AU/assembly based (approximate).
5. **AF cutoff date.** Using 2021-09-30 for the AFDB-complex models is an assumption [UNVERIFIED].
6. **Partner-chain choice.** For SBS sites the partner is the chain with the largest interface. MBS sites with ≥3 contributing chains are reduced to the top-2 chains.
7. **Confident models only.** Only confident models (ipSAE ≥0.6, pDockQ2 ≥0.23) were available, so the results apply to confident interfaces only.
8. **Isolated-chain pockets** are computed on chains cut from the dimer model, not on independent monomer predictions. The atlas comparison covers that case at residue level only.
