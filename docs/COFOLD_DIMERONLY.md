# Deep dive #3: does AF3 need the partner chain for ligands in curated dimer-only pockets?

**Status: UNVERIFIED, EXPLORATORY.** All numbers come from one run of the scripts below: AF3, seed 0, 5 diffusion samples,
n = 1 run per condition, 8 cases. Nobody has inspected the structures by eye. The metric code is a first attempt (one
check of the symmetry RMSD is described under Sanity checks).

## Question
The curated dimer-only pockets of human AF dimer models (`results/part2/q23a_curated.tsv`; `docs/Q22_RESULTS.md`
section 6) are mostly canonical active sites of obligate-homodimer enzymes. These pockets are missing from the
isolated-chain prediction and from the monomer pocket atlas. When the natural ligand, cofactor or a known drug is
co-folded into these sites, does AF3 need the partner chain? Ground truth is the PDB complex.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/cofold_select_dimeronly.py    # cases.tsv + AF3 JSONs (~5 min, census read + UniProt REST)
scripts/cofold_run_dimeronly.sh 2 3                           # AF3 on GPUs 2 and 3 (never 0); ~1 h 20 min incl. ~9 min start-up
.venvs/part1/bin/python scripts/cofold_metrics_dimeronly.py   # metrics.tsv, metrics_samples.tsv, summary.tsv
.venvs/part1/bin/python scripts/cofold_figures_dimeronly.py   # figures/*.png
```
- spyrmsd 0.9.0 and networkx 3.7 were added to `.venvs/part1` with `uv pip install --no-deps spyrmsd networkx`.
- Big outputs are in `data/cofold/dimeronly/` (gitignored): AF3 outputs, inputs, caches, `rejects.tsv`, `candidates_all.tsv`.
- Logs are in `logs/cofold_dimeronly/`; the run is recorded in `logs/compute_ledger.txt`.
- Reused code:
  - `cofold_select_onaxis.py` is imported for the MSA helpers;
  - `cofold_metrics_onaxis.py` is imported for the geometry helpers;
  - the run and split scripts are copies of the on-axis ones.

## Case selection (`scripts/cofold_select_dimeronly.py`; every rule is a constant in the code)
- **Curated pocket → canonical UniProt.**
  - Isoform suffixes are stripped.
  - A TrEMBL accession is replaced by the unique reviewed human accession with the same gene name (UniProt REST).
  - Pocket residues are mapped from model numbering to UniProt positions by aligning each model chain's sequence.
- **PDB sites.** Census v2 primary set, restricted to:
  - a single CCD component, not covalent;
  - not a metal ion, inorganic or metal cluster;
  - contacting both proteins of the pair: ≥ 2 chains for a homodimer, both accessions for a heterodimer.
- **Overlap.** The overlap coefficient between the ligand's contact residues and the curated pocket residues must be ≥ 0.3. Positions are (accession, UniProt position); for homodimers the two chains are pooled.
- **One case per family** (census cluster30) and per protein pair.
  - Families are taken in curated-rank order.
  - Within a family the order is: post-2021-09-30 first, then overlap, then resolution.
- **Stoichiometry.** Bound copies of the CCD on the reference chain pair (assembly 1). Run (B) gets 1 copy.
- **Sequences and limits.**
  - PDB SEQRES trimmed to the UniProt-aligned span.
  - Cached AF3 MSA remapped onto it.
  - No templates.
  - ≤ 1,280 tokens.
- **Funnel:**
  - 114 curated pockets; 137/139 accessions mapped to a canonical entry;
  - 843 census sites on these proteins; 928 pocket × site pairs;
  - 285 pairs with overlap ≥ 0.3, of which 248 also contact both proteins of the pair;
  - 12 curated pockets have at least one such site.
- **The 8 cases** (one per family):

  | case | rank | ligand | PDB | released | dimer tokens | copies |
  |---|---|---|---|---|---|---|
  | STING1 | 3 | A1A4R (agonist) | 9CUC | 2025 | 490 | 2 |
  | SGPL1 | 4 | 30J (inhibitor) | 4Q6R | 2014 | 1086 | 2 |
  | IVD | 5 | FAD | 8SGR (EM) | 2024 | 958 | 2 |
  | ACADVL | 6 | FAD | 3B96 | 2008 | 1280 | 2 |
  | MAT2A | 12 | U8E (allosteric inhibitor) | 8QE1 | 2024 | 852 | 2 |
  | PCCB | 74 | BTI (biotin) | 7YBU (EM) | 2023 | 1108 | 2 |
  | SOD1 | 94 | 9JT | 5O40 | 2018 | 338 | 2 |
  | IDH3G–IDH3A (hetero) | 110 | CIT (citrate) | 8GRH | 2022 | 706 | 1 |

- **Lost families.** These had candidate sites but produced no case (`data/cofold/dimeronly/rejects.tsv`):
  - GCDH: no bound copy passed the majority-contact rule in 8 entries. This may be a selection artefact; see Limitations.
  - PADI2 (1,476 tokens) and PFKP (1,562 tokens): over the token limit.
  - VEGFA: the cached MSA differs from UniProt.
  - The ACADVL post-cutoff entries were over 1,280 tokens, so the pre-cutoff 3B96 was used.
- **Not represented.** GCAT, DPYD and OGDHL have no qualifying PDB site in the census: no human PDB ligand complex, or the protein is too large.
- **Release dates.** 5 of the 8 cases are post-cutoff.

## Results (seed 0, top-ranked sample; `results/cofold/dimeronly/{metrics,summary}.tsv`)
Success means a symmetry-corrected ligand RMSD ≤ 2 Å after global CA superposition.

| | Dimer + ligand(s) | Single chain + 1 ligand |
|---|---|---|
| Success, top-ranked | **6/8** | **5/8** |
| Success in any of 5 samples | 7/8 | 5/8 |
| Median symmetry-corrected RMSD | 0.51 Å | 1.04 Å |
| Median DCC after pocket superposition | 0.26 Å | 0.39 Å |
| Ligand contacts both chains | 8/8 | n/a |
| Post-cutoff successes (n = 5) | 4/5 | 3/5 |

Per case (dimer → monomer, symmetry-corrected RMSD in Å):

| Case | Dimer → monomer | Outcome |
|---|---|---|
| STING1·A1A4R | 0.49 → **5.0** | The only clear case where the partner is needed. The site is split 7/7 residues across the two chains. The monomer keeps the ligand at its half-site (it contacts all 7 reference site residues of that chain) but in the wrong pose. |
| IDH3G–IDH3A·CIT | 0.52 → 1.40 | Both succeed. The monomer is worse. |
| SGPL1·30J | 0.40 → 0.60 | Both succeed. |
| IVD·FAD | 0.51 → 0.67 | Both succeed. |
| ACADVL·FAD | 0.21 → 0.26 | Both succeed. Pre-cutoff. |
| MAT2A·U8E | 0.50 → 0.68 | Both succeed. |
| SOD1·9JT | 3.5 → 2.1 | Both fail. The site is 5 + 1 residues; the best dimer sample is 0.95 Å. |
| PCCB·BTI | 4.8 → 7.6 | Both fail. Biotin is small, the site has 3 + 3 residues, and it is cryo-EM. |

Observations (descriptive; n = 8):
- **Dimer input.** All dimer runs place a ligand copy contacting both chains.
  - Where the dimer succeeds, every reference copy is recovered (2/2; 1/1 for IDH3). The exception is SOD1, where the dimer fails.
  - The inter-chain arrangement is correct: iRMSD 0.26–0.80 Å, except STING 2.3 Å.
- **Single-chain input.**
  - AF3 still puts the ligand into the half-site of the one chain. Site fraction on that chain is 0.83–1.0, except PCCB 0.33.
  - For the large cofactor sites (FAD in IVD/ACADVL) and drug sites where the partner contributes 31–50% of the reference site residues (SGPL1, MAT2A, IDH3), the pose stays sub-Å. For FAD, most contacts and the dominant chain sit on one chain.
  - So for 5 of 6 dimer successes, the single chain alone was enough for a correct pose relative to that chain.
- **Distance to the curated pocket centre** (`dcc_curated`).
  - Median 4.2 Å for dimer runs and 5.9 Å for monomer runs.
  - SGPL1 is 12–14 Å: its curated pocket overlaps the 30J site only partly (overlap 0.36). The curated pocket is a larger dimer cleft, and the ligand occupies part of it.
- **Answer to the question (provisional).** For these PDB-backed dimer-only pockets, AF3 mostly does **not** need the partner chain to place the ligand correctly relative to the chain it binds most. The partner chain:
  - completes the site (contacts on both chains);
  - is required for the one symmetric two-half-site drug pocket (STING);
  - gives modestly lower RMSDs (median 0.51 vs 1.04 Å).

  This contrasts with the "dimer-only" label from pocket detection: a pocket can be absent from the isolated-chain *apo* prediction while AF3 co-folding still builds the half-site around the ligand.

## Sanity checks
- **Runs.** All 16 runs finished (EXIT 0). All 80 samples have finite coordinates, with no AF3 clash flags in the top-ranked samples. Only GPUs 2 and 3 were used.
- **Symmetry RMSD.** Swapping the coordinates of the equivalent FAD phosphate oxygens O1A/O2A gives a name-matched RMSD of 0.49 Å and a spyrmsd RMSD of 0.0 Å. So the symmetry correction is active (checked on 3B96 and 8SGR).
  - For most top-ranked poses the name-matched and symmetry-corrected RMSDs agree.
  - They differ for IDH3·CIT (4.05 → 0.52 Å; citrate's two equivalent carboxymethyl arms), STING (0.90 → 0.49 Å) and MAT2A (1.04 → 0.50 Å).
- **Reference copies.** Bound reference copies on the chain pair use the same majority-contact rule in selection and metrics. This gives 2 copies per homodimer case and 1 for IDH3.
- **No known-answer re-run.** No independent known-answer test was done beyond the pipeline already exercised in COFOLD_ONAXIS.md.

## Limitations
- **Small n and selection bias.**
  - 8 cases, seed 0 only, no CIs; this is a descriptive case series.
  - Only pockets with human PDB ligand complexes ≤ 1,280 tokens qualify, so the top curated hits GCAT, DPYD and OGDHL are not tested.
- **Memorisation.**
  - 3 of 8 cases are pre-cutoff (SGPL1, ACADVL, SOD1).
  - Post-cutoff entries may have pre-cutoff homologue or re-determined complexes; for example, FAD in acyl-CoA dehydrogenases and STING agonists are well represented before 2021.
  - Templates were off, but the MSA and weights still carry information.
- **Monomer RMSD.** The monomer RMSD is computed against the full reference pose after superposition on one chain. A half-site ligand can pass even though the partner residues are absent. `ref_site_minor_frac` (0.17–0.50) gives the share of the site that a monomer cannot reproduce.
- **Curated-pocket DCC.**
  - It uses curated residues mapped from model to canonical to construct numbering, and only the residues present in the run.
  - For a monomer run it uses only half of the pocket, so dimer and monomer values are not strictly comparable.
- **Selection artefact (suspected).** The selection script's chain-contact helper uses `gemmi.NeighborSearch` with the unit cell, so crystal-symmetry images can be counted.
  - The `ref_contacts` column in `cases.tsv` (e.g. ACADVL 30,30) is therefore unreliable; the metrics script recomputes site sizes without symmetry images.
  - This may also explain the GCDH `no_bound_copy` rejections. Not investigated.
- **Shared environment.** spyrmsd was installed into the shared `.venvs/part1` environment.
- **No DockQ.** iRMSD only.

## Files
- **Scripts:** `scripts/cofold_select_dimeronly.py`, `scripts/cofold_split_jobs_dimeronly.py`, `scripts/cofold_run_dimeronly.sh`, `scripts/cofold_metrics_dimeronly.py`, `scripts/cofold_figures_dimeronly.py`.
- **Tables:** `results/cofold/dimeronly/{cases,metrics,metrics_samples,summary}.tsv`.
- **Figures:** `results/cofold/dimeronly/figures/{lig_rmsd_by_case,pocket_placement}.png`.

## Seeds 1–2 replicate (UNVERIFIED)
Status: first attempt, not inspected by eye. AF3 with the same inputs and settings as seed 0, only `modelSeeds` changed (5 samples per seed; project seed list [0,1,2]). Dimer and monomer for all 8 cases. Code: `scripts/cofold_seeds_jobs.py`, `scripts/cofold_run_seeds.sh`, `scripts/cofold_metrics_seeds.py` (calls `cofold_metrics_dimeronly.evaluate` unchanged; seed-0 values reproduce `metrics.tsv` exactly). Tables: `results/cofold/dimeronly/seeds_{metrics,summary}.tsv`.

**Coverage (complete, UNVERIFIED).** Seeds 1 and 2 are complete for all 8 cases, both conditions (48/48 runs incl. seed 0). The missing runs were completed on 2026-10-30 with AF3 from a node-local /tmp copy (see the on-axis doc). Previous partial: seed 2 had only 3 of 8 cases (dimer 2/3, monomer 1/3; pooled 6/9 vs 3/9; all 19 runs 14/19 vs 11/19).

| Condition | Seed 0 | Seed 1 | Seed 2 | Pooled, 3 seeds |
|---|---|---|---|---|
| Dimer, top-ranked success | 6/8 | 6/8 | 6/8 | 18/24 |
| Monomer, top-ranked success | 5/8 | 5/8 | 5/8 | 15/24 |

**Per-case consistency.** All 16 case x condition pairs have the same outcome at all three seeds.
- Dimer fails at every seed for SOD1 (3.5, 5.0, 3.7 A) and PCCB (4.8, 10.5, 5.2 A); all other dimers succeed (RMSD 0.2-0.6 A).
- Monomer fails at every seed for STING1 (5.0, 5.1, 5.1 A), SOD1 (2.1, 2.2, 2.3 A) and PCCB (7.6, 9.8, 12.2 A); IDH3 succeeds at all seeds (1.4, 1.6, 1.6 A).
- Reading: the seed-0 picture is reproduced by seeds 1 and 2 (STING1 needs the partner chain; the other cases do not). Not inspected by eye.
