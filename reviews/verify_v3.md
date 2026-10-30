# Independent verification of REPORT.md v3 additions (2026-10-30)

Verifier: independent agent. It did not produce these results and did not import any producer script. Code and logs are in `reviews/verify_v3/`. All runs were CPU-only with `scripts/env_threads.sh` and ≤4 workers.

Environment note: imports from the `.venvs/part1` site-packages stalled for more than 2 min on the network FS. numpy/pandas/scipy/pyarrow/gemmi were copied to `/tmp/vsp_v3` and run with the base uv CPython (`reviews/verify_v3/py.sh`). The fpocket env was used in place, and its I/O was slow (138–598 s per site).

## 0. Part 1 MBS fractions, REPORT §3.1 (PI add-on; `v0_part1_fractions.py`, `v0*.log`)
Source: `data/part1/ligand_sites.parquet`, `in_primary_set`. The flag equals my own rule (not dup ∧ bound ∧ not glycan ∧ not additive ∧ BioLiP2-covered) except for 3,909 rows. Those 3,909 are exactly the `same_uniprot_hetero` instances, which v2 excludes by design.

| Number (REPORT line) | Recomputed | Verdict |
|---|---|---|
| Instance-level homomer 24.6% (148) | 24.64% (n = 260,607) | VERIFIED |
| Homomer cofactor 33.8% / metal 11.8% / small molecule 31.3% (148) | 33.77 / 11.78 / 31.29 | VERIFIED |
| Instance-level heteromer 30.8% (149) | 30.81% (n = 135,858) | VERIFIED |
| Heteromer without photosynthetic complexes 23.2% (149) | 23.20% (`is_photosynthetic`) | VERIFIED |
| Heteromer cofactor 35.5 / metal 9.0 / small molecule 43.3 (149) | 35.54 / 8.99 / 43.26 | VERIFIED |
| Sites 26.3% (23.9–28.6) (152) | 26.32% of 183,660 sites (pdb × site_group; MBS = any copy); cluster bootstrap 24.0–28.6 | VERIFIED |
| Units 30.6% (28.1–33.1) (152) | 30.59–30.70% (CI 28.1–33.1); see note | VERIFIED (definition-sensitive ±0.1 pp) |
| Homomer units 29.2% (153) | 29.20% under every variant | VERIFIED |
| Heteromer units 34.9% (153) | 34.77–35.21% depending on the cluster rule; 34.93% with the per-site cluster rule | VERIFIED (definition-sensitive) |

Note on units:
- 423 collapsed sites span more than one `cluster30`, because the dominant chain differs between copies, mostly in heteromers.
- Assigning clusters per instance gives 44,231 units. Assigning one cluster per collapsed site gives 44,186 or 44,164 units, against the producer's 44,178.
- The unit-level rule is not specified in DESIGN or PART1_RESULTS. It should be documented.

## 1. Uncropped PDB site descriptors, REPORT §3.3 (`v1_pdbsite_delta.py`, `v1_pdbsite_delta.tsv`)
Method:
- Data: `data/interim/pdbsite/sites.parquet`, 11,983 sites, 5,532 MBS, 4,256 clusters.
- Weights: units_available / sampled per stratum, from `sample_strata_n12000_seed0.tsv`. They are identical to the producer's `w`.
- δ: weighted Cliff's δ written from scratch.
- CI: my own cluster30 bootstrap (1,000 resamples, seed 0).

| Number | Recomputed | Verdict |
|---|---|---|
| burrow_score δ −0.009 [−0.071, +0.056] (199, 303) | −0.0093 [−0.071, +0.056] | VERIFIED |
| Merged volume δ +0.071 [+0.010, +0.132] (197, 200, 488) | +0.0709 [+0.010, +0.132] | VERIFIED |
| Charged (merged lining) δ −0.009 [−0.068, +0.047] (201, 488) | −0.0085 [−0.068, +0.047] | VERIFIED |
| Hydrophobic δ +0.027 [−0.030, +0.078] (201, 488) | +0.0269 [−0.030, +0.078] | VERIFIED |
| Cropped merged volume +0.107 (197) | +0.1069 in the producer crop30 table; not recomputed | Consistent with source |
| Adjusted p 0.86, FE β −0.03 p 0.014; volume adj. p 0.03, FE p 4e-4 (199–200) | Match `mbs_vs_sbs_regression.tsv` (pct scale: 0.863; −0.029 / 0.014; 0.027; 3.7e-4); not re-fitted | Consistent with source |
| "11,983 sites … from 4,257 clusters" (195) | The 11,983 described sites come from **4,256** clusters. 4,257 is the 12,000-site sample. | DISCREPANT (minor) |
| "cropping changes burrow_score at cropped sites (24% of MBS, 18% of SBS)" (197) | 24.4% / 17.8% are the shares of sites where the crop removed chains. burrow_score changed by >0.05 at 61% / 60% of those sites. | Wording imprecise |

Also checked, all matching the doc:
- located-member volume +0.046;
- merged lining +0.081;
- aromatic −0.052;
- buriedness −0.023.

## 2. fpocket, REPORT §4.3 (`v2_fpocket_stats.py`, `v2_fpocket_rerun.py`, `v2_fpocket_rerun.tsv`)
Data: `data/interim/fpocket/pdbsite_sites.parquet`, status ok: 2,901 sites, 1,461 MBS, 1,690 clusters. Weights are units_available / n_ok per stratum and are identical to the producer's.

| Number | Recomputed | Verdict |
|---|---|---|
| δ +0.01 [−0.08, +0.09] (296; 28) | +0.0078 [−0.077, +0.090] | VERIFIED |
| >0.5: 40.3% vs 39.3% (297) | 40.31% vs 39.27%; diff +1.0 pp [−5.8, +7.1] | VERIFIED |
| Drug-like δ −0.13 [−0.26, +0.01] (298) | −0.132 (unweighted, as the doc states). The CI was taken from the producer table, not re-bootstrapped. | VERIFIED (point) |
| On-axis vs off-axis n.s. (299) | Producer table δ −0.03 [−0.19, +0.25] | Consistent with source |
| AF dimers: median 0.010 vs 0.006, paired p 0.69 (300) | `afdimer_summary.tsv` 0.010 / 0.006; `afdimer_paired.tsv` Wilcoxon p 0.692. AF fpocket was not rerun. | Consistent with source |

**fpocket rerun.**
- Sample: 10 random ok sites (`DataFrame.sample(10, random_state=0)` on site_id-sorted rows).
- Pipeline, my own code: raw mmCIF from `/path/to/pdb_mirror` → gemmi assembly 1 (first altloc, no H) → polymer only → `fpocket -f` → pocket with the most ligand heavy atoms within 3 Å of its alpha-sphere centres.
- Result: **10/10 match exactly.** The selected pocket number, overlap count, ligand atom count, number of pockets and Druggability Score are all identical.
- The 10 sites cover 5 strata (no SBS drug-like site was drawn). Assembly sizes were 2.2k–18.9k atoms.

## 3. LIGYSIS calibration, REPORT §4.1 Comparability (`v3_ligysis.py`)
| Number | Recomputed | Verdict |
|---|---|---|
| 57.6% vs 48.4% paired recall (253–254, 100) | 0.5760 vs 0.4839 on 467 sites from 192 chains, matched by SITE_NAME at seg level. Discordant pairs 84 / 41, exact McNemar p = 1.5e-4 (not reported in REPORT). At acc level the AFDB figure would be 49.5%. | VERIFIED |
| Full atlas 49.3% [47.9–50.6] (255) | 49.27% (A denominator, 10,682 sites; NaN counted as miss); protein bootstrap 47.9–50.6 | VERIFIED |
| Internal check: hit flag equals min_dcc_topN2 ≤ 12 Å | 0 disagreements in both tables | OK |
| "200 experimental chains" (100) | 200 chains were sampled, but the paired 57.6/48.4 comparison uses 192 chains and 467 sites | Wording imprecise |

**Metric against the paper.** I read the full text of Utgés & Barton 2024 (PMC11552181) via the Europe PMC XML. The paper:
- uses "top-N+2, N = number of observed sites for a target protein";
- states "Reported recall is obtained using DCC = 12 Å";
- reports fpocket+PRANK 60.4%, DeepPocket RESC 58.1%, P2Rank CONS 53.9% and P2Rank 51.9%;
- uses 2,775 chains and 6,882 sites.

**The DCC ≤ 12 Å / top-N+2 metric matches the paper's definition.** The remaining non-comparability is already disclosed in the doc:
- the site centre is a CA-centroid proxy, not the ligand-site centroid;
- AFDB / atlas filtering is used;
- the site set is 10,689 sites, not the 2,775-chain subset.

REPORT line 253 says 57.6% is "close to the published 60.4%". That comparison is across different site sets and centre definitions; the caveat is stated in line 252.

## 4. Co-folding seeds and heterodimer extension, REPORT §5 (`v4_seeds_hetero.py`)
Success means RMSD ≤ 2 Å, recomputed from `lig_rmsd_sym` and `rmsd_ver`. It agrees with the stored flags in 0-disagreement checks.

| Number | Recomputed | Verdict |
|---|---|---|
| On-axis dimer 7/12, 7/12, 7/10 (395) | 7/12, 7/12, 7/10 | VERIFIED |
| On-axis single chain 3/12, 3/12, 2/11 (396) | 3/12, 3/12, 2/11 | VERIFIED |
| Dimer-only dimer 6/8, 6/8, 2/3 (397) | 6/8, 6/8, 2/3 | VERIFIED |
| Dimer-only single chain 5/8, 5/8, 1/3 (398) | 5/8, 5/8, 1/3 | VERIFIED |
| Seed 2 is 31/52 jobs (390) | 25 on-axis + 6 dimer-only seed-2 rows = 31; 52 = 36 + 16 per seed | VERIFIED |
| Per-case agreement 10/12, 11/12, all dimer-only and off-axis (400) | 10/12, 11/12, 8/8 + 8/8, 6/6 + 6/6 | VERIFIED |
| STING1 succeeds only as a dimer at every seed (401) | Dimer 0.49 / 0.48 / 0.38 Å; single chain 5.0–5.1 Å | VERIFIED |
| On-axis discordant 4 vs 0, McNemar p 0.125; 11/12 sign p ≈ 0.006 (365–366) | 4 / 0, p = 0.125; 11/12, p = 0.0063 | VERIFIED |
| Hetero 11/21 vs 7/21, discordant 5 / 1, McNemar p 0.22 (362, 376, 457) | 11 / 7, 5 / 1, p = 0.219 | VERIFIED |
| Dimer lower RMSD 18/21, sign p 0.0015 (377, 457) | 18/21, p = 0.00149 | VERIFIED |
| Partner share ρ 0.16 (378) | ρ = 0.157 (p = 0.50) | VERIFIED |
| 6 original + 15 new; 17/21 with a complex before the cutoff (362) | 6 + 15; 17/21 among complete pairs. paired.tsv has 22 rows, and the 22nd (8bik QTN, monomer missing) is also pre-cutoff. | VERIFIED |
| Verifier RMSD agrees on 43/43 calls (375) | 43 ok rows, success_ver = success_prod in all of them | VERIFIED |
| 5 largest cases did not finish (379) | 26 cases → 21 complete pairs; 9 missing rows | Consistent |

## 5. REPORT.md cross-check: mismatches and imprecisions
1. `REPORT.md:195`: "4,257 clusters" should be 4,256 for the 11,983 described sites. The doc's 4,257 refers to the 12,000-site sample.
2. `REPORT.md:197`: "cropping changes burrow_score at cropped sites (24% of MBS, 18% of SBS)". 24% / 18% are the fractions of sites that were cropped, not the fractions with a changed burrow_score. Suggested wording: "the crop removed chains at 24% of MBS / 18% of SBS sites, and changed burrow_score by >0.05 at ~60% of those".
3. `REPORT.md:100`: "200 experimental chains". The paired numbers use 192 chains and 467 sites.

No other v3 number in lines 7, 28, 100, 193–205, 252–255, 291–303, 349–402, 457 or 488 disagreed with its source table or with my recomputation. The status tags `[UNVERIFIED]` on these items can now be changed to VERIFIED, with the items in §0 and §5 above noted. Tables cited as "consistent with source" were only matched against producer outputs, not re-fitted.
