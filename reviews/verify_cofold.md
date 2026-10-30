# Independent verification of the AF3 co-folding ligand RMSDs (on-axis, dimer-only and H4 batches)

**Verifier:** a separate agent. It did not produce any of the co-folding runs or metrics, and it does not import the producers' `cofold_*` code.

**Scope:** every seed-0 top-ranked prediction in the three batches, plus the seed-1 runs that had finished by the time of checking.

**Status:** this is a numerical re-computation. Nobody looked at the structures by eye.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python reviews/verify_cofold/verify_rmsd.py --seeds 0 --workers 4 --sites          # ~1 min
.venvs/part1/bin/python reviews/verify_cofold/verify_rmsd.py --seeds 1,2 --workers 4 \
    --out reviews/verify_cofold/rescored_seeds12.tsv
```
- **Outputs.**
  - `reviews/verify_cofold/rescored_seed0.tsv` and `rescored_seeds12.tsv`: per job;
  - `compare_seed0.tsv`: joined to the producers' `metrics.tsv`;
  - `site_share.tsv`.
- **CCD definitions.** These are fetched from `files.rcsb.org/ligands/download/<CCD>.cif` into `reviews/verify_cofold/ccd_cache/`. That folder is gitignored and can be re-fetched.

## Method (mine, stated because it differs in places from the producers')
- **Top-ranked sample.** The sample with the highest `ranking_score` for the requested seed in `<job>_ranking_scores.csv`.
- **Reference structure.**
  - Built from PDB assembly 1, using gemmi `make_assembly` with AddNumber chain naming. Hydrogens are removed and only the first altloc is kept.
  - The two reference protein chains come from the case tables' chain columns (`pchains`, `chainA_pdb`/`chainB_pdb`). These are inputs to the runs, not results.
- **Residue mapping.** Each predicted chain is aligned globally (gemmi) to the *observed* residues of each reference chain. Only identical aligned residues are paired.
- **Superposition.** Kabsch on all paired CA atoms of all predicted protein chains, i.e. global CA, the same as the producers. Pocket superposition is not used.
- **Chain assignment.**
  - Every assignment of predicted chains to reference chains with ≥ 90% identity is tried. For a single-chain run, that means either copy of a homodimer.
  - The reported RMSD is the **minimum over assignments**, on the grounds that identical copies are symmetry-equivalent.
  - The RMSD under the assignment with the best CA RMSD is also kept (`lig_rmsd_cabest`).
- **Ligand symmetry.**
  - Uses all element-preserving automorphisms of the CCD heavy-atom bond graph (networkx VF2). There are 1–93,312 per ligand (4mtz 0BU); none were truncated.
  - The RMSD is taken over the atoms present in both copies.
  - It is the minimum over (predicted copy, reference copy) pairs. Reference copies are residues of that CCD within 5 Å of either reference chain.
- **Success** means a ligand RMSD ≤ 2 Å.

## A1. Agreement with the producers' success calls (seed 0, top-ranked, 112 structures)

| Batch | Structures | Calls agree | Producer RMSD column compared |
|---|---|---|---|
| On-axis | 48 | 46 | `lig_rmsd`, name-matched (the documented primary metric) |
| Dimer-only | 16 | 15 | `lig_rmsd_sym` |
| H4 | 48 | **48** | `lig_rmsd`, symmetry-corrected |
| **All** | **112** | **109 (97%)** | |

- **Other comparisons.**
  - My on-axis calls agree 36/36 with the producers' own later symmetry-corrected re-score of seed 0 (`results/cofold/onaxis/seeds_metrics.tsv`, `success_sym`).
  - Seed 1: 31 runs were re-scored. My calls agree 12/12 with the on-axis `seeds_metrics.tsv` rows that overlap.
- **Size of the RMSD differences.** The median difference (mine − producer) is −0.004 Å. Where the difference is > 0.3 Å (24 rows), mine is always the lower value (maximum difference +0.10 Å). The causes are listed below.

**The three disagreements** (all within 0.3 Å of the cutoff):

| Case | Producer | Mine | Cause |
|---|---|---|---|
| on-axis `sbs_8ii3_PKK` dimer (TTR) | 2.80 (name-matched) → fail | 1.97 → success | Ligand symmetry: the producers' Hungarian lower bound is 1.77. The call is also sensitive to chain assignment: under the CA-best assignment I get 2.01, a failure. Treat it as a coin flip at the cutoff. |
| on-axis `sbs_8w0r_A1AEU` monomer (EBP) | 2.11 (name-matched) → fail | 1.80 → success | Ligand symmetry, which the producers' doc already notes as a Hungarian flip. |
| dimer-only `SOD1_5o40_9JT` monomer | 2.08 → fail | 0.85 → success | Chain choice for the single chain. The producers score the single chain only against `mono_chain` A1, which gives 2.08 (I reproduce this as `lig_rmsd_cabest`). Matching the single chain to F1 instead, against F1's own 9JT copy, gives 0.85. The two 9JT copies in 5o40 are not symmetric (centroids 4.4 Å apart), so this is a definitional choice. |

- **Effect on the headline numbers.**
  - Under my definitions: on-axis SBS dimer 3/6 and SBS monomer 2/6, matching the producers' Hungarian-based totals; dimer-only single chain 6/8 instead of 5/8.
  - No dimer-only or H4 dimer call changes, and no on-axis drug call changes.
- **Large RMSD differences that do not change a call.** All are failures in both scorings, at 3–36 Å:
  - name-matched vs symmetric RMSD in on-axis cofactors (GSH 1.01 → 0.21 Å; still successes);
  - different chain or copy choice on wrongly assembled dimers (9ckf SAH 20.7 vs 7.7 Å, CA RMSD 22.7 Å);
  - 9mel 1II (5.0 vs 3.1 Å);
  - 8x71 YCC.
- **CA RMSD.** It differs from the producers' value in a few cases: STING1 4.1 vs 2.1 Å, and 8f78 4.7 Å under the min-ligand assignment vs 1.0 Å under the CA-best one. This is because I pair all identical observed residues, while the producers use their own residue-mapping scheme. The ligand calls are unaffected.

## A2. Sanity checks
- **Finite coordinates.** All 112 seed-0 and all 31 seed-1 top-ranked models are finite.
- **Ligand–protein clashes (min heavy-atom distance < 1.5 Å).** 4/112, all H4 dimer runs, all already failures:
  - 8r5q Y5N 0.97 Å;
  - 8y55 A1D52 (NAMPT) 1.14 Å;
  - 8g9v YYC 1.33 Å;
  - 8dyg U5Q 1.38 Å.

  AF3's own `has_clash` is 0 for all four, because that flag concerns protein chains. No clashes among the seed-1 models. The median minimum distance is 2.6 Å.
- **Duplicate ligands.** 1/112: in the dimer-only `SOD1_5o40_9JT` dimer, the two predicted 9JT copies sit on top of each other (centroids 1.9 Å apart, closest atoms 1.0 Å). The reference copies are 4.4 Å apart, with closest atoms 3.1 Å. This is a physically implausible prediction and a failure (3.0 Å). It should be noted in COFOLD_DIMERONLY.md.
- **Partial reference ligand.** The PDB copy of 8pyv 48T (H4) is modelled with 21/35 heavy atoms, so its RMSD covers 21 atoms only. Every other case uses all CCD heavy atoms.
- **Reference copy for the site.** The copy nearest each case's centroid is ≤ 0.01 Å from it in all 56 cases.

## Verdict
- The success calls are reproduced: 109/112 at 2 Å.
- All three disagreements sit at the cutoff. Two come from the on-axis doc's name-matched primary metric, which that doc already caveats. One comes from a defensible but different single-chain chain choice (SOD1).
- No disagreement changes a qualitative conclusion of the three deep dives.
- Two items should be added to the producers' docs:
  - the SOD1 duplicate-ligand dimer prediction;
  - four clashing ligands in H4 dimer failures.

## Limitations of this check
- One verifier and one code path. The structures were not inspected by eye.
- Both scorings use global CA superposition, so neither separates a wrong pose from a wrong dimer arrangement.
- The "min over chain assignments" rule is more lenient than a fixed-chain rule.
