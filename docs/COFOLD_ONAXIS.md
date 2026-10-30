# Deep dive #1: AF3 co-folding of on-axis drug sites in human homodimers

**Status: UNVERIFIED.** All numbers come from one run of the scripts below: AF3, seed 0, 5 diffusion samples, n = 1 run
per condition. Nobody has inspected the structures by eye. The metric code is a first attempt.

## Question
Part 1 (H3, descriptive) found that ligands sitting on the twofold axis of homomers are mostly drug-like, and that
cofactors rarely sit there. This deep dive asks three things:
- Does AlphaFold3 put these drugs on the axis when it is given the homodimer?
- Does it fail when given only one chain, where the site does not exist?
- How does it behave for the two controls: off-axis MBS cofactor sites and SBS drug-like sites?

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/cofold_select_onaxis.py      # cases.tsv + AF3 JSONs in data/cofold/onaxis/inputs (~6 min)
scripts/cofold_run_onaxis.sh 2 3                             # AF3 on GPUs 2 and 3 (never 0), ~1.5-2 h; skips finished jobs
.venvs/part1/bin/python scripts/cofold_metrics_onaxis.py     # metrics.tsv, metrics_samples.tsv, summary.tsv
.venvs/part1/bin/python scripts/cofold_figures_onaxis.py     # figures/*.png
```
Large AF3 outputs are in `data/cofold/onaxis/af3/` (gitignored). The logs are in `logs/cofold_onaxis/`.

## Case selection (`scripts/cofold_select_onaxis.py`; every rule is a constant in the code)
- **Sites.** Part 1 primary collapsed sites (`part1_analysis.load` + `site_filter`, primary variant), restricted to:
  - homomer, with exactly 2 protein chains in assembly 1;
  - one CCD component, not covalent;
  - accession organism Homo sapiens (taxon 9606), from the UniProt REST API.
- **Strata:**

  | Stratum | Rule |
  |---|---|
  | on-axis drug-like | MBS, census `is_on_axis` (Jaccard), drug-like |
  | off-axis cofactor | MBS, not on-axis, cofactor class |
  | SBS drug-like | not MBS, drug-like |

- **One case per cluster30.** Within a cluster, a post-2021-09-30 release is preferred, then the best resolution.
- **Order across clusters:** post-cutoff first, then size. On-axis cases are ordered by fewest tokens. Controls are ordered by how close they are to the median on-axis token count.
- **Size limit:** tokens ≤ 1,024 for the dimer with its ligands. This is our cap for the 4 h budget; COFOLDING.md allows ≤ 1,280.
- **Ligand copies.**
  - On-axis cases: 1 copy, as the task specifies.
  - Controls: the number of bound, deduplicated copies of the CCD in assembly 1. It must equal 2, so that it matches the biological stoichiometry.
- **AF3 sequence.**
  - The PDB entity SEQRES, trimmed to the first and last residue that aligns to UniProt. This drops expression tags.
  - MSA columns come from the cached full-length AF3 MSA (`/cluster/project/beltrao/shared/25.05_alphafold3_msas_human`). The case is used only if the cached query equals the UniProt sequence fetched by the script.
  - Unaligned positions become gap columns.
- **No templates are used.** The cached templates go up to 2025 and would include the holo entries. This is a deviation from COFOLDING.md, which used templates.
- **Funnel:**
  - 6,513 candidate collapsed sites in 2-chain homomers;
  - 2,125 after the human, single-component and non-covalent filters;
  - by cluster: 18 on-axis, 51 off-axis cofactor and 265 SBS drug-like.
  - Rejection reasons are in `data/cofold/onaxis/rejects.tsv`. Four on-axis candidates were rejected because the cached MSA query differs from the UniProt sequence.

## Cases (`results/cofold/onaxis/cases.tsv`)
- **24 cases:** 12 on-axis drug-like, 6 off-axis cofactor, 6 SBS drug-like.
- **Release date.** 20 cases were released after the AF3 training cutoff (2021-09-30). Four on-axis cases are pre-cutoff: 6RV3 KCNK3, 4JA8 IDH2, 4H30 MMP12, 4MTZ XIAP. All 12 controls are post-cutoff.
- **On-axis targets:**
  - IL-17A (7AMA), PD-L1 (9I0U) and caspase-6 (8F78);
  - KCNK6 (9MEL), STING (8STI) and NAPE-PLD (8PC4);
  - BRDT and BRD4 with the bivalent ligand ZN1 (8CZA, 7MR7);
  - KCNK3 (6RV3), IDH2 (4JA8), MMP12 (4H30) and XIAP BIR3 (4MTZ).
- **Off-axis cofactor controls:** GSTP1·GSH, PNPO·FMN, MAT2A·SAM, GOT1·PLP, SMYD2·SAH and NADK2·NADP.
- **SBS drug-like controls:** BLVRB, TTR, PRDX1, FABP7, galectin-1 and EBP.
- **Size.** Dimer tokens are 253–918. On-axis median is about 390; controls 322–918.
- **Caveat on the controls.** "Homodimer" means a 2-chain assembly 1. For some SBS controls this is probably not the biological unit:
  - FABP7 and BLVRB are monomeric proteins;
  - TTR is a tetramer.

  AF3 gives these low ranking scores, and they inflate the dimer-condition failures in that stratum.

## AF3 runs
- **Runs:** 24 cases × (A) dimer + ligand and (B) monomer + 1 ligand. Seed 0, 5 diffusion samples, no templates, GPUs 2 and 3.
  - All 48 runs finished.
  - All 240 samples have finite coordinates (`metrics_samples.tsv` column `finite`).
  - Wall time was about 1 h 25 min.
- **Seed 1** was added for the 12 on-axis dimer runs only, because of the time budget.
- **Ligand copies in (A):** on-axis cases get 1 copy; controls get 2.

## Results (seed 0, top-ranked sample; `metrics.tsv`, `summary.tsv`)
Success means a ligand RMSD ≤ 2 Å after superposing CA atoms. Both chain mappings are tried, and atoms are matched by CCD name.

| Stratum | Dimer: success | Monomer: success | Dimer: ligand contacts both chains | Dimer: ligand ≤ 3 Å from predicted C2 axis |
|---|---|---|---|---|
| On-axis drug-like (n = 12) | **7/12** | **3/12** | 11/12 | 11/12 |
| Off-axis MBS cofactor (n = 6) | 5/6 | 6/6 | 5/6 | 0/6 |
| SBS drug-like (n = 6) | 2/6 | 1/6 | 0/6 | 0/6 |

**Split by the AF3 training cutoff (on-axis only):**
- Pre-cutoff: dimer 3/4, monomer 1/4.
- Post-cutoff: dimer 4/8, monomer 2/8.
- All controls are post-cutoff, so their rows above are post-cutoff as well.

**Seed 1 (on-axis dimers only).** Seed 1 gives the same result as seed 0.
- 7/12 successes, in the same 7 cases.
- By cutoff: pre 3/4, post 4/8.
- 12/12 ligands contact both chains; median axis distance 0.62 Å.
- The largest changes are BRD4·ZN1 (13.4 → 8.5 Å, still a failure) and KCNK6 (best sample now 1.9 Å, top-ranked 4.4 Å).

**What the case-level data show** (n is small; descriptive only):
- **Dimer input, on-axis.** Whether or not the pose is right, AF3 almost always puts the drug on the twofold axis and in contact with both chains:
  - median distance from the predicted axis 0.49 Å (11/12 within 3 Å);
  - in contrast, off-axis and SBS ligands sit 10–32 Å from the axis.
- **Monomer input, on-axis.** Ligand RMSD is worse than with the dimer in 11 of the 12 cases. Four cases succeed with the dimer but fail with the monomer: PD-L1, caspase-6, MMP12 and XIAP. No case succeeds only as a monomer.
- **Monomer successes.** In the three on-axis cases where the monomer succeeds (IL-17A, STING, KCNK3), a single chain contributes most of the contacts. With the dimer the RMSD is still lower or similar.
- **Where the monomer ligand goes.** In most on-axis monomer runs the ligand stays at the half-site of the single chain: a median 87% of that chain's true site residues are contacted, and 9/12 cases contact ≥ 50%. It is placed elsewhere in caspase-6, IDH2 and NAPE-PLD: site_frac 0, ligand 19–35 Å away.
- **Off-axis cofactors** are predicted well with or without the partner chain (monomer 6/6). The cofactor sits mainly in one chain's pocket.
  - The one dimer failure, SMYD2·SAH, is a wrong inter-chain arrangement (CA RMSD 22.7 Å).
  - Its pocket-superposed ligand RMSD is 4.4 Å. In the separate monomer run the SAH pose is 0.5 Å.
- **On-axis dimer failures:**
  - **BRDT/BRD4 + bivalent ZN1.** This dimer is probably induced by the ligand. AF3 builds a different symmetric BD dimer (CA RMSD 15–19 Å) but still puts ZN1 on its axis.
  - **NAPE-PLD.** Wrong dimer (CA RMSD 23 Å); the ligand contacts only one chain.
  - **IDH2 + allosteric 1K9.** The dimer is correct (CA RMSD 1.1 Å) and the ligand is on the axis (0.8 Å) and in contact with both chains, but the pose is wrong (7.5 Å). Pre-cutoff.
  - **KCNK6.** The channel dimer is correct; the pose is 5.0 Å off, with a best of 5 samples of 2.8 Å.
- **Confidence.** Ranking scores are high in most successes (0.78–0.98). But high confidence also occurs in wrong poses, e.g. IDH2 dimer 0.82 and NAPE-PLD monomer 0.99 with RMSD 19 Å. ipTM is not a reliable success filter here (n too small to quantify).

## Sanity checks
- **GPU outputs.** All 240 seed-0 samples and all 60 seed-1 samples have finite coordinates, and every run has a ranking_scores file. Only GPUs 2 and 3 were used.
- **Known-answer re-run.** COFOLDING.md has no case that was compared against an experimental structure. Its TPI1+PGA test was checked only by confidence, so there was no known answer to reproduce.
  - The closest substitute is a post-cutoff on-axis case that reproduces with sub-Å accuracy: IL-17A·RMK (7AMA, released 2022-04), with dimer ligand RMSD 0.75 Å, CA RMSD 0.59 Å and all 5 samples < 1 Å.
  - This exercises the whole pipeline: input building, MSA remapping, residue mapping and the metrics. It is not an independent known-answer test.
- **Bug found and fixed during analysis.** The first metrics pass had an off-by-one residue mapping, which gave a systematic CA RMSD of about 3.8 Å. The fix is in `_ref_structure`: label_seq is 1-based. After the fix the correctly folded cases have CA RMSD 0.4–1.5 Å.

## Limitations
- **Small n:**
  - 12 / 6 / 6 cases, one run per condition, no CIs;
  - on-axis pre-cutoff n = 4;
  - all controls are post-cutoff.
  This is a descriptive case series, not a test.
- **Memorisation.** Templates are off, but MSAs and weights may still carry pre-cutoff information. Post-cutoff entries can also be close homologues or re-determinations of pre-cutoff complexes; this was not checked per case (e.g. other PD-L1 or STING ligand complexes from before 2021).
- **Monomer RMSD.** In monomer runs the RMSD is computed against the full on-axis ligand pose. A ligand that correctly sits in the half-site can still pass when the other chain adds few contacts. For that reason `site_frac` is reported as well.
- **Ligand-induced dimers.** These include BRD + ZN1 and XIAP BIR3 + a bivalent Smac mimetic. They are a special case: the dimer exists only because of the ligand.
- **Symmetry handling.** Ligand atom matching is by CCD name. `lig_rmsd_hung` (a Hungarian match within each element) is a lower bound on the symmetry-corrected RMSD. Under it, 3 of 48 seed-0 top-ranked calls flip to success; all three were 2.1–2.8 Å by name match:
  - MMP12 monomer;
  - TTR SBS dimer;
  - EBP monomer.

  With the Hungarian match the totals become: on-axis monomer 4/12, SBS dimer 3/6, SBS monomer 2/6. RDKit or spyrmsd symmetry-corrected RMSD was not used, because spyrmsd is not installed.
- **No DockQ.** It is not installed; interface backbone iRMSD (10 Å CAPRI interface) is reported instead.
- **Off-axis controls are matched only roughly by size.** They are larger: median 851 tokens vs 405 for on-axis (SBS 399).

## Files
- **Scripts:** `scripts/cofold_select_onaxis.py`, `scripts/cofold_split_jobs.py`, `scripts/cofold_run_onaxis.sh`, `scripts/cofold_seed1_onaxis.py`, `scripts/cofold_metrics_onaxis.py`, `scripts/cofold_figures_onaxis.py`.
- **Tables:** `results/cofold/onaxis/{cases,metrics,metrics_samples,summary}.tsv`.
- **Figures:** `results/cofold/onaxis/figures/{lig_rmsd_by_case,success_rates}.png`.

## Seeds 1–2 replicate (UNVERIFIED)
Status: first attempt, not inspected by eye. AF3 with the same inputs and settings as seed 0, only `modelSeeds` changed (5 samples per seed, no templates, GPUs 2 and 3). The project's seed list is [0,1,2]. Scope: 12 on-axis + 6 off-axis cofactor cases, dimer and monomer; SBS controls not rerun.

**Code.** `scripts/cofold_seeds_jobs.py` (inputs + job split, `--seeds`), `scripts/cofold_run_seeds.sh` (runner with restart loop), `scripts/cofold_metrics_seeds.py` (metrics and tables). Original scripts are unchanged. Order was seed-major (all of seed 1 before any seed 2); within a seed, cases by dimer token count ascending, both experiments pooled.

**Coverage (complete, UNVERIFIED).** Seeds 1 and 2 are now complete for all 18 cases, dimer and monomer. The 21 missing jobs were run on 2026-10-30 with AF3 from a node-local /tmp copy of the env and weights (`scripts/stage_af3_local.sh`, `scripts/af3_local_env.sh`; Lustre import stalls, LESSONS rule 7). Previous partial: seed 2 had 12/18 dimers and 13/18 monomers (on-axis dimer seed 2 7/10, monomer 2/11, pooled dimer 21/30, monomer 8/33).

**Metrics.** Same metric code as the other experiment. The original on-axis RMSD is name-matched (`lig_rmsd_orig`). The symmetry-corrected RMSD (spyrmsd, as in the dimer-only experiment) is `lig_rmsd_sym`. Success is RMSD <= 2 A for the top-ranked sample of each seed. Seed 0 reproduces `metrics.tsv` for the original RMSD, and the dimer-only seed-0 values reproduce exactly. Tables: `results/cofold/onaxis/seeds_{metrics,summary}.tsv` (summary has per_seed, pooled = cases with all seeds present, pooled_all_runs, per_case).

| Stratum, condition | Seed 0 | Seed 1 | Seed 2 | Pooled, 3 seeds |
|---|---|---|---|---|
| On-axis dimer, symmetry-corrected | 7/12 | 7/12 | 7/12 | 21/36 (58%) |
| On-axis dimer, original RMSD | 7/12 | 7/12 | 6/12 | 20/36 |
| On-axis monomer, symmetry-corrected | 3/12 | 3/12 | 2/12 | 8/36 (22%) |
| Off-axis cofactor dimer (6 cases) | 5/6 | 5/6 | 5/6 | 15/18 |
| Off-axis cofactor monomer (6 cases) | 6/6 | 6/6 | 6/6 | 18/18 |

The pooled rows are per-run rates, not independent samples.

**Per-case consistency (all 3 seeds, 12 on-axis cases)**:
- On-axis dimer: 10 of 12 consistent (4 fail at every seed, 6 succeed at every seed). Inconsistent: XIAP BIR3 (4MTZ) 1.14, 0.97, 3.06 A (original RMSD 1.7, 1.6, 3.3) and KCNK6 (9MEL) 3.15, 2.12, 1.40 A (1/3).
- On-axis monomer: 11 of 12 consistent. Inconsistent: KCNK3 (6RV3) 1.54, 1.43, 6.78 A.
- Off-axis cofactors: all 12 case x condition pairs consistent; SMYD2/SAH (9CKF) fails dimer at all seeds (about 21 A).
- Reading: the seed-0 result (7/12 dimer vs 3/12 monomer) is unchanged by the complete seeds 1-2 (pooled 21/36 vs 8/36). Not inspected by eye.
