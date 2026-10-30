# Co-folding synthesis: does the partner chain help more when the PDB site is split across chains?

**Status: EXPLORATORY, UNVERIFIED by a second party.** The analysis plan below was written and committed before any
of the synthesis numbers were computed. The results section was added afterwards.

## Pre-declared analysis (written before computing)
- **Cases.** Every co-folding case in the three deep dives that has both a dimer run and a single-chain run:
  - on-axis batch (24 cases; `docs/COFOLD_ONAXIS.md`);
  - dimer-only batch (8 cases; `docs/COFOLD_DIMERONLY.md`);
  - H4 batch (24 cases; `docs/COFOLD_H4.md`).

  That is 56 cases in total. Primary analysis: seed 0, top-ranked sample. Secondary: seeds 0–2 pooled per case (mean RMSD over the available seeds, success rate over seeds), for the cases whose seed-1/2 runs exist when this analysis is finished.
- **Pose metric.** The ligand RMSD is recomputed independently (`reviews/verify_cofold/verify_rmsd.py`), not taken from the producers' tables:
  - global CA superposition;
  - symmetry-aware matching through automorphisms of the CCD heavy-atom graph;
  - minimum over symmetry-equivalent protein chain assignments and over pairs of predicted and reference ligand copies.

  Success means RMSD ≤ 2 Å.
- **Partner share (predictor).** For the case's PDB site:
  - Use the reference ligand copy nearest the case centroid in assembly 1.
  - A residue on either of the two reference protein chains is a binding residue if it has ≥ 2 heavy-atom pairs with distance ≤ vdW_i + vdW_j + 0.5 Å. This is the BioLiP rule; gemmi vdW radii are used.
  - partner_share = (binding residues on the chain with fewer binding residues) / (all binding residues). Its range is 0–0.5.
- **Outcomes.**
  1. partner_gain_rmsd = single-chain RMSD − dimer RMSD (Å). It is positive when the dimer is better. Because RMSDs of failed poses are large and noisy, the gain is also computed on log10 RMSD.
  2. success_diff = dimer success − single success, which takes the values −1, 0 or +1.
  3. dimer_only = dimer success AND single-chain failure (binary).
- **Tests.**
  - Spearman ρ of partner_gain_rmsd vs partner_share, and of success_diff vs partner_share.
  - Logistic regression: dimer_only ~ partner_share.
  - 95% CIs from a percentile bootstrap over cases (2,000 resamples). Resamples in which the logistic fit fails to converge or separates the data are counted and reported.
  - One-sided reading per the hypothesis: the gain increases with partner_share.
- **Hypothesis (from the deep dives).** Cofactor half-sites (low partner share) do not need the partner. Even-split sites (on-axis drugs, STING1) do.
- **Sensitivity.**
  - The analysis is repeated per batch, and with the ligand-induced or non-biological dimers removed: BRD ZN1 (8cza, 7mr7), the FABP7, BLVRB and TTR SBS controls, and XIAP 4mtz.
  - The H4 batch uses 1 ligand copy for NAMPT 8y55. This is kept as is.
- **Figure.** partner_gain_rmsd vs partner_share, coloured by batch/class, saved as a PNG at dpi 150.

No other analyses are planned. Anything added later is labelled post hoc.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python reviews/verify_cofold/verify_rmsd.py --seeds 0 --workers 4 --sites
.venvs/part1/bin/python reviews/verify_cofold/verify_rmsd.py --seeds 1,2 --workers 4 --out reviews/verify_cofold/rescored_seeds12.tsv
.venvs/part1/bin/python reviews/verify_cofold/synthesis.py     # ~5 min (bootstraps)
```
- **Outputs.**
  - Cases: `results/cofold/synthesis/cases_seed0.tsv`, `cases_seeds012_pooled.tsv`.
  - Tests: `results/cofold/synthesis/stats.tsv`.
  - Figures: `results/cofold/synthesis/partner_gain_vs_share.png` and `partner_gain_vs_share_seeds012.png`.
- **How the RMSDs were checked.** The RMSDs are verified against the producers' calls in `reviews/verify_cofold.md` (109/112 agree). The synthesis itself has had no second-party check.

## Results (added after computing)

### Partner share by class
Median partner share:

| Class | Median partner share |
|---|---|
| SBS drug | 0.00 |
| Off-axis cofactor | 0.16 |
| H4 MBS | 0.28 |
| Dimer-only pocket | 0.36 |
| On-axis drug | 0.48 |

The BioLiP rule is stricter than the producers' contact rule. For the dimer-only cases it gives smaller sites, but the ordering is similar. The one exception is SOD1 5o40, which gives 4 + 0 residues here vs 5 + 1 in the producers' count. 7 cases have fewer than 5 site residues, so their share is coarse.

### Seed 0, all 56 cases (pre-declared primary)
Totals: dimer success 32/56, single-chain success 27/56. There are 9 dimer-only successes and 4 single-only successes.

| Test | Estimate | 95% bootstrap CI | Bootstrap P(≤ 0) |
|---|---|---|---|
| Spearman ρ, RMSD gain vs share | 0.34 | 0.10 to 0.55 | 0.004 |
| Spearman ρ, log10 RMSD gain vs share | 0.32 | 0.07 to 0.53 | 0.009 |
| Spearman ρ, success difference vs share | 0.17 | −0.12 to 0.41 | 0.12 |
| Logit slope, dimer-only success ~ share | 2.7 per unit share | −1.6 to 10.6 | 0.11 (1 failed fit) |

### Sensitivity analyses
- **Without the ligand-induced and non-biological dimers (n = 50).**
  - ρ = 0.30 (0.02 to 0.53) for the RMSD gain.
  - The logit slope is 4.5 (0.2 to 12.4), so its CI now excludes 0.
- **Post hoc: one row per PDB entry and ligand.** This removes the duplicate 8ros FMN case (n = 55). The results are essentially unchanged: ρ = 0.34, logit 2.7 (−1.7 to 10.0).
- **Per batch.** All point estimates for the RMSD gain are positive, but the CIs are wide:
  - on-axis (n = 24): ρ 0.39 (0.03 to 0.66);
  - H4 (n = 24): ρ 0.34 (−0.14 to 0.70);
  - dimer-only (n = 8): ρ 0.19 (−0.64 to 0.88).
- **Seeds 0–2 pooled.**
  - Only seed 1 had finished, for 14 cases in both conditions (10 on-axis drugs, GSH and FMN off-axis cofactors, STING1, SOD1), when this analysis was completed. Seed 2 was not yet available.
  - ρ = 0.27 (−0.45 to 0.78) for the RMSD gain; logit slope 6.5 (−4.4 to 21.2).
  - The point estimates go in the same direction, but there is too little data for this to be informative.

### Descriptive reading (post hoc split at share 0.3)
| Group | Dimer-only successes | Single-only successes | Median log10 gain |
|---|---|---|---|
| Share < 0.3 (n = 26) | 3 | 2 | 0.03 |
| Share ≥ 0.3 (n = 30) | 6 | 2 | 0.20 |

- **Low share.** Of the 3 dimer-only successes at low share, one is TTR 8ii3, which is borderline (1.97 vs 2.07 Å). The other two are H4 hetero 6wih EDT and homodimer 8pyv 48T, both with share 0.25.
- **High share without a partner gain.** Several even-split sites need no partner: the MAT2A·SAM cofactor (0.47), IVD·FAD (0.44), SGPL1 (0.50), IL-17A 7ama (0.50) and STING 8sti (0.50).
- **Single-only successes.** There are four. Three of them have a visible dimer-run failure mode:
  - SOD1: duplicate overlapping copies;
  - 9ckf: wrong dimer;
  - 8tcz: wrong top-ranked sample, although 3/5 dimer samples succeed.

  7p0o 4OY has no explanation.

## Interpretation (exploratory)
- **Direction.** The direction matches the hypothesis. The more of the PDB site lies on the partner chain, the more AF3's ligand RMSD improves when the partner is present (ρ ≈ 0.3, CI excludes 0).
- **Success calls.** The effect on binary success is weaker and its CI includes 0 in the primary analysis. It is driven by the on-axis drugs and the heterodimer interface ligands (8f78 caspase-6, 7msr, STING1, 4h30, 4mtz, 9i0u).
- **Not "cofactor vs drug".** Cofactors with even-split sites (MAT2A·SAM, IVD/ACADVL·FAD) are predicted well from one chain. So a high partner share is roughly necessary for a partner gain, but not sufficient.
- **Caveats.**
  - These are 56 cases from three convenience-selected batches.
  - The predictor is computed from the same PDB structures that define success.
  - Memorisation is not controlled.
  - One seed in the primary analysis.
  - The p-values are bootstrap tail fractions and are not corrected for the several tests and subsets.
