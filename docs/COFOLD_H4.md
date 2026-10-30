# Deep dive #2: AF3 co-folding at H4 multi-chain (MBS) ligand sites

**Status: UNVERIFIED.** One AF3 run per condition (seed 0, 5 diffusion samples). Nobody has inspected the structures
by eye; one case (8qdy, MAT2A) was re-scored with an independent gemmi superposition (ligand RMSD 0.49 vs 0.48 Å
dimer, 1.37 vs 1.37 Å single chain, name-matched). n is small (≤ 24 cases), so every rate below is descriptive.

## Question
H4 (docs/H4_RESULTS.md, reviews/verify_h4.md) found that ligand-free AFDB dimer models recover PDB MBS sites better
than isolated chains, and that the dimer misses are mostly oversized merged interface pockets or apo tunnels that are
too narrow (NAMPT). Here: does AF3 co-folding of the dimer WITH the ligand place it correctly at MBS sites, including
sites the apo pocket finder missed, and does it beat co-folding with a single chain?

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/cofold_select_h4.py    # cases.tsv + AF3 JSONs in data/cofold/h4/inputs (~3 min)
scripts/cofold_run_h4.sh 1                             # AF3, one process on GPU 1; skips finished jobs
.venvs/part1/bin/python scripts/cofold_metrics_h4.py   # metrics.tsv, metrics_samples.tsv, summary.tsv
.venvs/part1/bin/python scripts/cofold_figures_h4.py   # results/cofold/h4/figures/*.png
```
AF3 outputs: `data/cofold/h4/af3/` (gitignored); logs `logs/cofold_h4/`.

## Case selection (`scripts/cofold_select_h4.py`, all rules are constants in the code)
- H4 v2 primary set (`data/part2/h4/v2/benchmark_filtered.parquet`, organic, interface-filtered), MBS,
  `ligand_class` drug_like or cofactor, single CCD component, not covalent (census). 568 sites / 119 clusters, 499
  sites after the component/covalent filter.
- Strata: **homo_miss** (homodimer, H4 apo `dimer_hit` False), **homo_hit**, **hetero** (heterodimer, either).
- NAMPT and MAT2A are not hard-coded: the two MBS families with the most negative contribution to the H4 DiD
  (`reviews/verify_h4/family_contributions.tsv`) are forced first into homo_miss.
- One case per cluster30, no cluster reused across strata; post-cutoff release (> 2021-09-30) preferred, then
  fewer tokens. Size: ≤ 768 tokens (budget cap), forced cases ≤ 1,280.
- AF3 input: PDB SEQRES trimmed to the UniProt-aligned span, MSA columns remapped from the cached AF3 MSA (cached
  query must equal the UniProt REST sequence), **no templates** (cached templates include the holo entries).
- Ligand copies (run A): copies of the CCD in assembly 1 contacting both chains; 2 if ≥ 2 ("symmetric-paired").
  **Deviation:** the 1,280-token bucket ran out of GPU memory for NAMPT (2 copies, 1,034 tokens), so cases whose
  2-copy input exceeds 1,024 tokens get 1 copy (`lig_copies_reduced`; only 8y55). Run B: the chain with the most
  contacts + 1 copy.
- Rejects: `data/cofold/h4/rejects.tsv`.

## Metrics (`scripts/cofold_metrics_h4.py`)
Reuses the deep-dive #1 helpers (imported, unmodified). Per sample and for the top-ranked sample:
- **lig_rmsd**: ligand heavy-atom RMSD after CA superposition of the protein(s) onto the PDB pair; chain mapping
  symmetry-aware (both A/B assignments for homodimers; the single chain is also tried on the other copy); minimum
  over predicted ligand copies; **symmetry-corrected** over ligand graph automorphisms (bonds inferred from the PDB
  geometry). Name-matched and Hungarian lower-bound RMSDs are also in the table. Success = lig_rmsd ≤ 2 Å.
  Symmetry correction mattered for citrate (6l59: 4.0 → 0.54 Å dimer, 4.3 → 1.5 Å single chain; both calls flip
  to success) and changed no other success call.
- **dcc**: predicted vs PDB ligand centroid (pocket level), hit ≤ 4 Å (the H4 threshold).
- **contacts_both**: the matched predicted copy has a heavy atom within 4 Å of both chains (run A only).
- AF3 ranking_score, iptm, max ligand–protein chain_pair_iptm (`lig_prot_pair_iptm`), has_clash.

## Results (top-ranked sample, seed 0; `results/cofold/h4/{metrics,summary}.tsv`)
24 cases ran in both conditions (48 AF3 jobs, all finite, no clashes flagged). 21 cases are post-cutoff releases;
3 are pre-cutoff (all heterodimers). Only 3 cases are post-cutoff at the *pair* level (H4 definition), so almost all
pairs (and in several families, close analogues of the ligand) were in AF3's training PDB.

| Stratum | n | A dimer: RMSD ≤ 2 Å | B single chain: RMSD ≤ 2 Å | A: DCC ≤ 4 Å | B: DCC ≤ 4 Å | A: contacts both chains |
|---|---|---|---|---|---|---|
| Homodimer, apo-dimer MISS | 10 | 4 | 5 | 8 | 7 | 10 |
| Homodimer, apo-dimer HIT | 8 | 3 | 3 | 4 | 4 | 6 |
| Heterodimer | 6 | 4 | 2 | 5 | 4 | 6 |
| **All** | 24 | **11** | **10** | 17 | 15 | 22 |
| Post-cutoff release only | 21 | 8 | 8 | 14 | 12 | 19 |
| Pre-cutoff release (hetero) | 3 | 3 | 2 | 3 | 3 | 3 |

- **A vs B (paired, RMSD success):** both succeed 8, dimer only 3 (8pyv 48T, 6wih EDT, 7msr ZNS), single chain only 2
  (7p0o 4OY, 8tcz ZRI), both fail 11. Median RMSD 3.1 vs 4.3 Å; median DCC 1.8 vs 3.1 Å. On these numbers the dimer
  input does **not** clearly beat the single chain for homodimers; the advantage is concentrated in heterodimers
  (4/6 vs 2/6), where the single chain often places the ligand elsewhere (7msr: 0.70 vs 40 Å; 8gs8: 4.7 vs 13 Å).
- **Apo misses:** co-folding the dimer places the ligand within 2 Å at 4/10 apo-dimer misses and within 4 Å
  (centroid) at 8/10. It recovers the oversized-pocket misses MAT2A 8qdy (0.48 Å, 2 copies, both chains) and
  9hjc VIB (0.73 Å), 8ptt (0.52 Å) and 8ros FMN (0.58 Å); 7o4d (NQO2, oversized pocket) is a near miss (2.5 Å,
  DCC 2.1 Å). For several misses the single chain does as well, i.e. the site is formed largely by one chain.
- **NAMPT (8y55, A1D52):** fails in both conditions (15 / 14 Å; best of 5 samples 3.8 Å dimer), although ranking_score
  is 0.96. Consistent with the apo-tunnel failure mode in reviews/verify_h4.md, but this run used 1 ligand copy
  instead of the PDB's 2 (memory fallback), and it is one case.
- **Other failures:** 8dyg U5Q (IL-17A-type interface small molecule) is placed 43 Å away at a different interface
  groove; 7wjv LSA misplaced by 27 Å in both runs; 8tcz ZRI top-ranked dimer sample is 21 Å off but 3/5 dimer
  samples succeed (ranking picked the wrong sample).
- **Apo hits are not easier:** 3/8 succeed in both conditions. Apo-pocket recovery (DCC to a cavity centre) and
  co-folding pose accuracy measure different things.
- **Confidence:** among dimer runs, the ligand–protein chain_pair_iptm separates success from failure (median 0.96
  vs 0.72); ranking_score separates them weakly (0.94 vs 0.91).

Figures: `results/cofold/h4/figures/lig_rmsd_by_case.png` (per-case RMSD, A vs B, 5-sample range),
`results/cofold/h4/figures/rates_by_stratum.png`.

## Caveats
- n = 24 cases, one seed, no replicate; no statistical test was run. Differences of 1–2 cases are noise.
- Memorisation: 21/24 entries are post-cutoff releases, but 21/24 *pairs* had a pre-cutoff complex in the PDB, and
  for inhibitor series (MAT2A, NAMPT, NQO2) earlier co-crystals with related ligands exist. Templates were off.
- Selection prefers small proteins (≤ 768 tokens) for the 4 h budget, so large interface sites are under-represented.
- Run B uses the chain with the most contacts; a site made mostly by that chain favours B. contacts_both in A is
  nearly always true, so it does not discriminate here.
- The "miss" label comes from the H4 heuristic failure modes (`fm_dimer`); 3 of 10 misses are "near miss (DCC 4–5 Å)".
- Ligand stoichiometry rule and the NAMPT 1-copy fallback are described above; protein sequence = PDB construct
  span, so termini absent from the crystal are not modelled.
