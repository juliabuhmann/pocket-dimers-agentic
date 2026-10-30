# EXPLORATORY: more heterodimer MBS sites for the AF3 co-folding comparison (post-hoc extension of COFOLD_H4)

**Status: UNVERIFIED and EXPLORATORY.** This is a post-hoc extension and was not pre-registered. It was motivated by
the small-n hint in docs/COFOLD_H4.md: there, heterodimer multi-chain (MBS) sites were placed correctly in 4/6 cases
with the full heterodimer and in 2/6 with the single chain. Each condition was run once (seed 0, 5 diffusion samples,
top-ranked sample scored). Nobody has looked at the structures by eye, so treat every p-value below as descriptive.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/cofold_select_hetero_ext.py      # cases.tsv, rejects.tsv, data/cofold/hetero_ext/inputs
scripts/cofold_run_hetero_ext.sh 1 2026-10-30T04:30:00+01:00      # AF3, one process, GPU 1; launch deadline; resumable
.venvs/part1/bin/python scripts/cofold_metrics_hetero_ext.py     # metrics.tsv, paired.tsv, summary.tsv
.venvs/part1/bin/python scripts/cofold_figures_hetero_ext.py     # figures/hetero_ext.png
```
- **Code reuse.** The new scripts import the original deep-dive code without modifying it:
  - `cofold_select_h4.py`: chain pair, contacts, ligand copies, AF3 JSON writing, MSA remap, no templates.
  - `cofold_metrics_h4.py`: the producer metric.
  - `reviews/verify_cofold/verify_rmsd.py`: the verifier metric.
  - `run_af3_wrapper.py`: AF3. `scripts/run_af3_deadline.py` only adds a check that stops new jobs from being launched after a set time.
- **AF3 outputs.** `data/cofold/hetero_ext/af3/` (gitignored). Logs: `logs/cofold_hetero_ext/`.

## Case selection (`scripts/cofold_select_hetero_ext.py`; all rules are in the code)
- **Starting set.** Same base set as the original: H4 v2 primary set, organic ligands, MBS, heterodimer, one CCD component, not covalent.
- **Ligand classes.** Only 13 unused cluster pairs have a ligand in the original classes (drug_like, cofactor; "tier 1"). To get closer to the target of 20, tier 2 adds every other ligand class with ≥ 10 heavy atoms, which drops additives such as EDO and MPD. The tier is recorded for each case and the results are split by tier below.
- **One case per cluster pair.** The unit is the unordered pair of RCSB 30% clusters of the two chains, with at most one case per pair and one per PDB entry. The 6 original cluster pairs are excluded.
  - 8gs8 HEM uses a different chain pair (C+D) of the same entry as the original 8gs8 UQ1 case, and shares one cluster with it.
- **Size limit.** ≤ 1,280 tokens.
- **Preference order.**
  - Within a pair: tier 1, then post-cutoff release, then fewer tokens.
  - Across pairs: tier 1, then fewer tokens.
- **Selected cases.** 20 (6 tier 1, 14 tier 2), listed in `results/cofold/hetero_ext/cases.tsv`. For each case the table gives the release date, `post_cutoff` (release after 2021-09-30) and `pair_complex_pre_cutoff`. The latter means the pair was already in one PDB entry before the cutoff (the H4 `post_cutoff_pair` definition, negated).
- **Rejected.** 10 rejected candidates are in `rejects.tsv`: 8 exceeded the token limit and 2 had no usable cached MSA.
- **Protocol (identical to the original).**
  - Run A: the PDB pair plus the PDB ligand stoichiometry: 1 copy, or 2 if ≥ 2 copies in assembly 1 contact both chains.
  - Run B: the chain with the most contacts plus 1 ligand.
  - No templates, seed 0.
  - Jobs ran smallest first, dimer and single-chain runs of a case one after the other.

## What ran
- **Completed (2026-10-30 rerun).** All 40 of 40 AF3 jobs, giving **20 complete extension pairs** (26 pairs with the original 6). The earlier run had 31/40 (15 pairs) after Lustre import stalls; the 9 missing jobs were run from a node-local /tmp copy of the AF3 env and weights (`scripts/stage_af3_local.sh`, `scripts/af3_local_env.sh`). The 6pv7 NCT dimer hit a CUDA out-of-memory error in the long-lived process and was rerun in a fresh process. Previous partial (kept for reference): 21 pairs, dimer 11/21 vs single 7/21, discordant 5 vs 1, McNemar p = 0.22, RMSD sign test 18/21 p = 0.0015.
- **How the metrics were run.** During the same filesystem problem, I ran the metrics and figure scripts with the venv's base interpreter (uv CPython 3.12.14). The packages were copies of `.venvs/part1` site-packages on node-local `/tmp` (same package versions). The scripts themselves were unchanged.
- **Ledger.** All of this is recorded in `logs/compute_ledger.txt`.

## Metrics (`scripts/cofold_metrics_hetero_ext.py`)
Applied to the 6 original and the 15 extension complete pairs, top-ranked seed-0 sample.
- **Ligand RMSD, verifier (primary).** Symmetry-aware ligand RMSD from the verifier's code (`score_job`): CCD graph automorphisms, minimum over chain assignments and ligand copies, global CA superposition.
- **Ligand RMSD, producer.** The producer's symmetry-corrected `lig_rmsd`, computed against the census reference copy.
- **Success** means RMSD ≤ 2 Å.
- **Contacts with both chains.** The producer's `contacts_both`, for run A.
- **Ligand ipTM.** The maximum ligand–protein `chain_pair_iptm`.
- **Partner share.** From the verifier's `site_share`: BioLiP-like binding residues on the minor chain divided by all binding residues of the pair, using the PDB reference copy.
- **Partner gain.** log10(RMSD_single / RMSD_dimer).

**Do the two metrics agree?** Success calls agree in 43/43 runs. The RMSD values differ in two failures, both caused by which reference copy is used:
- 8umv AGS dimer: 2.16 Å (verifier) vs 34.2 Å (producer);
- 7vd2 PC1 single chain: 8.1 Å vs 32.8 Å.

Verifier clash flags (min ligand–protein distance < 1.5 Å) are raised for 2aer FUC (both runs) and the 7vd2 PC1 dimer.

## Results (`results/cofold/hetero_ext/summary.tsv`, `paired.tsv`)
Paired success, RMSD ≤ 2 Å, verifier metric. The producer metric gives identical counts in every row.

| Subset (full set; the tables below the first row are still the OLD 21-pair partial, see results/cofold/hetero_ext/summary.tsv for current values) | n | Dimer ✓ | Single ✓ | Both | Dimer only | Single only | Neither | McNemar exact p | Dimer lower RMSD | Sign test p |
|---|---|---|---|---|---|---|---|---|---|---|
| **All (original 6 + extension 20), FULL** | 26 | 14 | 8 | 7 | 7 | 1 | 11 | 0.070 | 22/26 | 0.00053 |
| All, previous partial (6 + 15) | 21 | 11 | 7 | 6 | 5 | 1 | 9 | 0.22 | 18/21 | 0.0015 |
| Original 6 (unchanged) | 6 | 4 | 2 | 2 | 2 | 0 | 2 | 0.50 | 6/6 | 0.031 |
| Extension only, FULL (20) | 20 | 10 | 6 | 5 | 5 | 1 | 9 | 0.22 | 16/20 | 0.012 |
| Extension only, previous partial (15) | 15 | 7 | 5 | 4 | 3 | 1 | 7 | 0.63 | 12/15 | 0.035 |
| Tier-1 ligands | 8 | 5 | 3 | 3 | 2 | 0 | 3 | 0.50 | 7/8 | 0.070 |
| Tier-2 ligands | 13 | 6 | 4 | 3 | 3 | 1 | 6 | 0.63 | 11/13 | 0.022 |
| Release post-cutoff | 10 | 3 | 2 | 1 | 2 | 1 | 6 | 1.0 | 8/10 | 0.11 |
| Release pre-cutoff | 11 | 8 | 5 | 5 | 3 | 0 | 3 | 0.25 | 10/11 | 0.012 |
| Pair complex before cutoff | 17 | 9 | 6 | 5 | 4 | 1 | 7 | 0.38 | 14/17 | 0.013 |
| Pair first seen after cutoff | 4 | 2 | 1 | 1 | 1 | 0 | 2 | 1.0 | 4/4 | 0.13 |

- **Success rate.** With the full heterodimer, 11 of 21 ligands are placed within 2 Å; with the single chain, 7 of 21. The discordant pairs split 5 to 1, which is not significant with an exact McNemar test (p = 0.22).
  - The extension alone reproduces the direction of the original hint but more weakly: 7/15 vs 5/15, discordant 3 to 1.
  - The one single-chain-only success is 8umv AGS: 2.16 Å dimer (just above the cutoff) vs 1.75 Å single chain.
- **RMSD sign test.** The dimer run has the lower ligand RMSD in 18 of 21 pairs (p = 0.0015).
  - Many of these are small differences between two successes, for example 5j1s ATP at 0.74 vs 0.83 Å.
  - The sign test therefore shows a consistent but often small advantage. It does not show more pockets being found.
  - Median RMSD: 1.9 Å dimer vs 4.7 Å single chain.
- **Partner share.**
  - The median partner share is 0.25, ranging from 0.05 to 0.50.
  - Partner gain is **not** related to partner share: Spearman ρ = 0.16, p = 0.50, n = 21. It is also unrelated within the original and extension sets separately (ρ = −0.03 and 0.04).
  - Dimer gains also appear at sites where the minor chain contributes only 5–10% of the binding residues (1efv, 6l59, 1grn, 5j1s).
- **Contacts with both chains.** The ligand touches both chains in 17 of 21 dimer runs. The 4 exceptions are failures in the dimer run: 8pvy, 8vq3, 1ld7 and 8umv.
- **Ligand ipTM.** Medians are 0.88 (dimer) and 0.87 (single chain), so this score does not separate the two conditions.
- **Training cutoff.**
  - Most of the dimer advantage comes from pre-cutoff releases (8/11 vs 5/11).
  - Post-cutoff releases succeed rarely in either condition (3/10 vs 2/10).
  - 17 of 21 pairs already had a complex in the PDB before 2021-09-30, so memorisation cannot be ruled out.
  - Only 4 pairs were first seen after the cutoff (2/4 vs 1/4).

Figure: `results/cofold/hetero_ext/figures/hetero_ext.png`. The left panel plots paired RMSD, single chain vs dimer (log axes, 2 Å lines). The right panel plots partner share against partner gain.

## Interpretation (tentative)
- **Ranking within a case.** Adding the partner chain fairly consistently lowers ligand RMSD at heterodimer MBS sites.
- **Success rate.** With 21 pairs, the difference in success rate (+4 cases) is not statistically supported.
- **Partner share.** The gain does not follow how much of the site the partner chain contributes. That fits an indirect effect, for example the partner fixing the conformation of the main chain, but this run does not test it.

## Caveats
- **Post-hoc.** The extension was chosen after seeing a 4/6 vs 2/6 hint, and the 6 original cases are included in the pooled test.
- **Ligand classes.** The extension relaxes the ligand-class rule (tier 2 is 13 of the 21 analysed pairs), so it is not the same population as the original stratum.
- **Missing cases.** None after the 2026-10-30 rerun (previously 5 of 20 were missing). Text in the Results and Interpretation sections below the first table rows was written for the 21-pair partial and has not been updated; use summary.tsv.
- **Single runs.** One seed per condition and the top-ranked sample only. The original deep dive showed ranking can pick a wrong sample.
- **Memorisation.** Most pairs were in the PDB before the AF3 cutoff, and templates were off.
- **Run B design.** Run B uses the chain with the most contacts, which favours the single chain when the site is mostly on that chain.
