# H5 — Disease variants at dimer-only pockets (Part 2, Q2.4)

**Status: UNVERIFIED.** This is a first pass. It has not been re-derived by an independent agent and nobody has inspected the residue-level output by hand.
Provenance: DESIGN.md sha256 `e3f9b550eb29cfe6eb444ebd884920b847687518470834c0bcee3439998eee28` (v3 FROZEN). The git commit is written into every results TSV header. Results files start with `# design_sha256=… git_commit=…`.

## Pipeline (rerun from scratch)
```
source scripts/env_threads.sh
python3 -m venv .venvs/h5 && .venvs/h5/bin/pip install -r envs/h5_requirements.txt         # freesasa, gemmi
mamba env create -p .venvs/h5_r -f envs/h5_r.yml                                             # R 4.4 + lme4 2.0.6 (lock: envs/h5_r_lock.txt)
.venvs/h5/bin/python   scripts/h5_residue_features.py --workers 8   # SASA, interface, pLDDT per model residue (~13 min)
.venvs/part1/bin/python scripts/h5_build_table.py                    # categories, aggregation, outcomes, ledger
.venvs/part1/bin/python scripts/h5_analysis.py --stage fit --workers 8   # lme4 fits via scripts/h5_glmm.R (~35 min/fit, full set)
.venvs/part1/bin/python scripts/h5_analysis.py --stage summarize
.venvs/part1/bin/python scripts/h5_figures.py
```
Note: the two PI-requested fits and the overlap 0.2/0.8 fits were launched by hand with the same `Rscript scripts/h5_glmm.R` call that `--stage fit` uses. They are listed in `JOBS`, so `--stage fit` reproduces them.

## Definitions as implemented
- **Dimer set**: every model in `afdb_human_dimers/{homo,hetero}dimers/candidates.csv` (10,610 models: 4,820 homo, 5,790 hetero after dedup on modelEntityId; 0 feature failures). Residue numbers in the models are UniProt positions. Check: 10,519 humsavar P/LP rows fall in model proteins, and only 1 has a wild-type mismatch against the model residue; that variant is dropped.
- **Pockets**: rows of `human_dimer_pockets_annotated.parquet` with `atlas_pass` (burrow_score ≥0.03 and lining pLDDT mean ≥70). For homodimers only symmetry representatives (`is_sym_rep`) are used, following the §0 collapse rule. That leaves 54,876 pockets. A lining residue counts only if its own pLDDT is ≥70; 44,654 of 867,552 model-chain lining rows were dropped this way.
- **Categories** (mutually exclusive, priority order):
  1. **dimer-only pocket**: `is_interface` (≥3 lining residues per chain) and `max_overlap_iso` <0.5.
  2. **single-chain pocket**: any other atlas-pass pocket. This includes interface pockets that also exist in the isolated chain, so that category 3 is "interface residues outside any pocket" as H5 states.
  3. **interface**: any heavy atom within **5.0 Å** of a heavy atom of the partner chain. DESIGN.md does not define this; 5.0 Å comes from the task specification and was fixed before any outcome was seen.
  4. **other surface**: relative SASA ≥0.2 in the dimer context.
  5. **buried**: everything else; excluded from the models.
- **SASA**: freesasa (Lee–Richards, ProtOr radii, probe 1.4 Å). Relative SASA is freesasa `relativeTotal`. The covariate `rsa_iso` is the isolated-chain relative SASA on the same coordinates (freesasa `separate-chains`), so interface burial does not count as burial.
- **Aggregation (primary)**: for each protein × UniProt position, the highest-priority category across every model and both chains. Covariates are the mean over the instances that carry that category. **Sensitivity**: one model per protein, the one with the highest max_ipSAE.
- **Outcomes**: (primary) the residue has ≥1 humsavar **P/LP** missense variant (UniProt humsavar 2026_03); (secondary) AlphaMissense mean pathogenicity per residue.
- **Data-source deviation (declared, decided before any results)**: DESIGN.md names ClinVar P/LP. The variant table in this project is UniProt humsavar (`docs/EXTERNAL_DATA.md`), which merges P with LP. ClinVar itself was not fetched.

## Model
`glmer(y ~ category + z(rsa_iso) + z(pLDDT) + (1|protein), binomial, Laplace nAGQ=1, bobyqa)` in lme4 2.0.6. The reference category is interface non-pocket and buried residues are excluded. A true GLMM fit in about 33 min on the full set, so no fallback was needed. CIs and p-values are Wald. **Conservation is not included as a covariate.** No conservation score was cheaply available in the workspace: there are no MSAs or per-residue conservation tables. AlphaMissense is an outcome, so it cannot be used as a covariate.
AlphaMissense: `lmer(am_mean ~ same, ML)`. AlphaMissense covers only 5,028 of the 8,279 proteins; the missing ones are mainly TrEMBL accessions that are absent from the aa_substitutions file.
The random-intercept SD is large (about 6 on the logit scale) because most proteins have no P/LP variant at all. Fixed effects are identified from within-protein contrasts. The Mantel–Haenszel OR stratified by protein (no covariates) is reported as a robustness check.

## n per step (`results/part2/h5/h5_ledger.tsv`, `h5_power_gate.tsv`)
| | residues | P/LP residues | proteins | proteins with ≥1 P/LP in category |
|---|---|---|---|---|
| dimer-only pocket | 109,511 | 547 | 5,740 | 254 |
| single-chain pocket | 401,695 | 2,061 | 7,634 | 595 |
| interface (ref) | 349,916 | 954 | 8,271 | 420 |
| other surface | 1,722,273 | 1,821 | 8,279 | 593 |
| buried (excluded) | 823,126 | 3,541 | — | 739 |
| model input | 2,583,395 | 5,383 | 8,279 | |

Power gate (≥30 proteins with P/LP in each compared category) **passes** for the primary contrast and for every secondary or sensitivity split, including the smallest (burrow top-2-quintile dimer-only: 116 proteins).

## Results
**Primary (pre-registered): dimer-only pocket vs interface non-pocket, P/LP, adjusted GLMM. OR 1.20 (95% CI 1.07–1.35), p = 0.0024.** Applying BH within family B together with H4 gives an adjusted p ≤0.0048, whichever H4 p is used (bootstrap 0.066 or permutation 1e-4).

Unadjusted P/LP rates with 95% CIs from a protein bootstrap (1,000 resamples, seed 0; `h5_category_rates.tsv`):

| category | rate |
|---|---|
| dimer-only | 0.50% (0.42–0.58) |
| single-chain pocket | 0.51% (0.45–0.58) |
| interface | 0.27% (0.23–0.32) |
| surface | 0.11% (0.09–0.12) |

Seeds 1 and 2 give the same CIs to within ±0.01 percentage points (`h5_category_rates_seeds.tsv`).

| analysis (`h5_model_results.tsv`) | label | dimer-only vs interface OR [95% CI], p |
|---|---|---|
| primary (aggregate per protein) | PRIMARY | **1.20 [1.07, 1.35], 0.0024** |
| one model per protein (best ipSAE) | SENSITIVITY | 1.18 [1.03, 1.36], 0.018 |
| dimer-only overlap <0.2 | SENSITIVITY (§0) | 1.18 [1.02, 1.36], 0.027 (174 proteins with P/LP) |
| dimer-only overlap <0.8 | SENSITIVITY (§0) | 1.26 [1.13, 1.40], 3.5e-5 |
| homodimer models only | SECONDARY | 1.22 [1.01, 1.47], 0.035 |
| heterodimer models only | SECONDARY | 1.16 [1.00, 1.34], 0.043 |
| ipSAE tertile 1 (<0.714) | SECONDARY | 1.04 [0.84, 1.31], 0.70 |
| ipSAE tertile 2 (0.714–0.806) | SECONDARY | 1.27 [1.05, 1.54], 0.012 |
| ipSAE tertile 3 (≥0.806) | SECONDARY | 1.17 [0.98, 1.39], 0.086 |
| models without a PDB complex of the pair | EXPLORATORY | 1.19 [1.01, 1.39], 0.033 |
| dimer-only restricted to burrow_score top-2 quintiles (≥0.110) | SENSITIVITY_PI (post hoc) | 1.24 [1.05, 1.47], 0.013 (44,178 residues, 116 proteins with P/LP; 65,333 lower-quintile dimer-only residues excluded) |
| shuffled categories within protein, seeds 0/1/2 | NEG_CONTROL | 1.06 [0.92, 1.22]; 1.08 [0.94, 1.25]; 0.96 [0.83, 1.11] |

**Other contrasts in the primary model (secondary)**:
- single-chain pocket vs interface: OR 1.17 [1.06, 1.28], p = 9e-4.
- other surface vs interface: OR 0.53 [0.48, 0.57], p = 2e-45.
- Covariates: rsa_iso per SD OR 0.69, pLDDT per SD OR 1.34.

**Pocket effect vs dimer-only effect (PI request, post hoc COMPARATOR)**: pre-existing interface pockets are interface pockets that also exist in the isolated chain (overlap ≥0.5). In a 5-level model (`cat_split`):
- pre-existing interface pocket vs interface: **1.26 [1.12, 1.42], p = 1e-4** (256 proteins with P/LP).
- dimer-only vs interface: 1.21 [1.07, 1.36], p = 0.0021.
- other single-chain pocket vs interface: 1.12 [1.02, 1.24], p = 0.021.

Pre-existing interface pockets show at least as large an enrichment as dimer-only pockets, which supports a pocket (interface-pocket) effect rather than a dimer-only-specific one.

Mantel–Haenszel OR stratified by protein (unadjusted, `h5_cmh.tsv`):
- dimer-only vs interface: 1.39 [1.23, 1.58] (456 strata).
- single-chain pocket vs interface: 1.58 [1.44, 1.73].
- **dimer-only vs single-chain pocket: 0.90 [0.80, 1.01], p = 0.07** (exploratory).

**AlphaMissense (secondary, LMM, 1.63 M residues, 5,028 proteins)**:
- dimer-only vs interface: +0.025 [0.023, 0.027], p ≈ 1e-154.
- single-chain pocket vs interface: +0.011.
- surface vs interface: −0.137.
- Unadjusted means: dimer-only 0.700, single-chain pocket 0.699, interface 0.634, surface 0.414.

## Sanity checks and controls
- **Negative control**: with category labels shuffled within protein, the dimer-only OR is 0.96–1.08 with every CI spanning 1 (seeds 0, 1, 2). This passes.
- **Known answer (TTR V30M)**: the UniProt signal peptide is 20 aa, fetched from the UniProt REST API. Mature V30 is therefore **UniProt V50**; the residue in the model is V and humsavar lists V50M/A/G/L as P/LP. In the AFDB TTR homodimer model (AF-0000000066262578, ipSAE 0.83), V50 is **buried** (rsa 0.001 in both dimer and isolated context, pLDDT 99) and **not at the dimer interface**: no heavy atom lies within 5 Å of the partner chain. The task premise that V30M sits at the dimer interface was not confirmed. V30 sits in the β-sheet core of the monomer and is thought to destabilise the tetramer indirectly. Separately, the TTR tetramer's T4 channel lies at the dimer–dimer interface, which a dimer model does not contain. What this check does confirm is the mapping and numbering (mature → UniProt offset, wild type matches). TTR P/LP residues by category: dimer-only 6, single-chain pocket 1, interface 4, surface 18, buried 19.
- **Expected direction**: P/LP rates are buried > interface > surface (buried 0.43% vs surface 0.11%), which matches the known enrichment of disease variants in protein cores.

## Interpretation and caveats
1. The pre-registered contrast is positive (OR ≈1.2) and robust to the one-model-per-protein sensitivity, the overlap thresholds, the homo/hetero split and PDB-complex exclusion. It weakens in the lowest ipSAE tertile.
2. **The enrichment does not look specific to dimer-only pockets.** Single-chain pockets show the same OR (1.17), and the unadjusted MH comparison of dimer-only vs single-chain pockets is 0.90 (p = 0.07). The most parsimonious reading is "pocket-lining residues carry more P/LP variants than non-pocket interface residues", not a dimer-specific effect. The PI's comparator rows above address this directly.
3. Q2.2 found that burrow output differs between dimer and isolated runs even away from the interface. Some "dimer-only" pockets are therefore probably detection noise, which would dilute any dimer-specific effect. The burrow top-2-quintile sensitivity is post hoc and requested by the PI; it is not pre-registered.
4. There is no conservation covariate. Pocket residues are typically more conserved, so part of the pocket effect may be conservation.
5. humsavar P/LP is used rather than ClinVar. It is curated, smaller, and biased toward well-studied disease genes.
6. Aggregating across all models of a protein gives proteins with many models more chances to be labelled dimer-only. The best-model sensitivity gives a similar OR.
7. The PDB-complex flag is approximate: it uses SIFTS chain-level mapping and asymmetric-unit co-occurrence, not assemblies.
8. The PI-requested analyses (burrow quintile, pre-existing interface pocket comparator) were added after the primary GLMM had been seen. They are labelled and do not replace the primary.

## Files
- Scripts: `scripts/h5_residue_features.py`, `scripts/h5_build_table.py`, `scripts/h5_glmm.R`, `scripts/h5_analysis.py`, `scripts/h5_figures.py`
- Tables: `results/part2/h5/`. `h5_model_results.tsv` (all fits, power-gate flag), `fits/*.tsv` (raw), `h5_power_gate.tsv`, `h5_category_rates*.tsv`, `h5_cmh.tsv`, `h5_sc_pocket_split_exploratory.tsv`, `h5_ttr_known_answer.{tsv,json}`, `h5_ledger.tsv`
- Figures: `results/part2/figures/h5_plp_rate_by_category.png`, `results/part2/figures/h5_forest.png`
- Data (not committed): `data/part2/h5/residue_features/`, `data/part2/h5/tables/`
