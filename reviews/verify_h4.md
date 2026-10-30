# Independent verification of H4 / Q2.1 (census v2)

Verifier: independent agent, 2026-10-29. I did not use `scripts/h4_*.py` or `src/part1`. All code is in `reviews/verify_h4/`:

| Script | Purpose | Output |
|---|---|---|
| `v1_stats.py` | Headline numbers, bootstraps, permutation, heterogeneity, dedup and weighting, collapse sensitivities | `v1_stats.log`, `stats_results.tsv`, `family_contributions.tsv`, `within_cluster_did.tsv` |
| `v2_dcc.py` | DCC from raw mmCIF using my own SIFTS-XML parser, my own assembly-1 chains (gemmi), my own Kabsch fit and my own top-k and symmetric collapse. Seed-0 samples by default; `--all-mbs --all` runs all 4,013 sites | `dcc_recompute.tsv`, `dcc_recompute_all.tsv` |
| `v3_extra.py` | Sign-flip, family sign tests, ICC/design effect, region-level permutation | `v3_extra.log` |
| `v4_failure.py` | Heuristic failure-mode census of the misses | `failure_modes.tsv`, `failure_modes_sites.tsv` |
| `v5_spotcheck.py` | PI spot-check table | `spotcheck.tsv` |
| `v6_nampt.py` | NAMPT ligand-in-model clash check | stdout |

Rerun: `source scripts/env_threads.sh; P=.venvs/part1/bin/python; $P reviews/verify_h4/v1_stats.py; $P reviews/verify_h4/v2_dcc.py; $P reviews/verify_h4/v2_dcc.py --all-mbs --all; $P reviews/verify_h4/v3_extra.py; $P reviews/verify_h4/v4_failure.py; $P reviews/verify_h4/v5_spotcheck.py; $P reviews/verify_h4/v6_nampt.py`. Every run takes about 1 minute or less.

What I did not recompute independently:
- The census itself: ligand definitions, `site_group`, the interface filter and the unit/representative selection.
- n (`n_sites_pair`), which sets k = n+2.
- The burrow pocket tables. These are treated as inputs.

## 1. Headline table

Source: `benchmark_filtered.parquet`, `lig_set == organic`, n = 4,013 (1,138 MBS, 2,875 SBS).

| Number | Producer | Verifier | Verdict |
|---|---|---|---|
| MBS hit rate, dimer / isolated | 0.383 / 0.263 | 0.3831 / 0.2627 | VERIFIED |
| SBS hit rate, dimer / isolated | 0.400 / 0.443 | 0.4000 / 0.4435 | VERIFIED |
| Homo MBS, SBS / hetero MBS, SBS (dimer, isolated) | 0.401, 0.278 / 0.432, 0.469 / 0.227, 0.134 / 0.283, 0.350 | identical to 3 decimals | VERIFIED |
| MBS gain / SBS gain | +0.120 / −0.043 | +0.1204 / −0.0435 | VERIFIED |
| DiD | +0.164 | +0.1639 | VERIFIED |
| Discordant pairs, MBS 310/173, SBS 369/494 | as stated | identical | VERIFIED |
| DiD cluster-bootstrap CI (1,000 resamples, seed 0) | −0.015 to 0.332 (kind-stratified) | unstratified −0.006 to 0.338; kind-stratified −0.005 to 0.331; my own clusters (474 components on the primary set) −0.006 to 0.326 | VERIFIED. The lower bound differs by about 0.01, which is Monte-Carlo or RNG noise. Every variant includes 0. |
| Bootstrap p | 0.066 | 0.066 (unstratified), 0.070 (stratified); seeds 1 and 2 give 0.068 and 0.072 | VERIFIED |
| MBS gain CI, SBS gain CI | (−0.048, 0.274), (−0.119, 0.025) | (−0.047, 0.269), (−0.122, 0.024) | VERIFIED |
| MBS / SBS cluster sign-flip p | 0.17 / 0.24 | 0.176 / 0.247 | VERIFIED |
| Within-cluster permutation p (centred), null mean | 1e-4, 0.091 | 1e-4, 0.0914 (null SD 0.011) | VERIFIED as computed. Its interpretation is wrong; see §4. |
| Number of clusters | 470 | 470 with the producer's `boot_cluster`. Rebuilt on the primary set only, I get 474 (4 producer clusters split; none merged). | VERIFIED |

## 2. DCC from raw coordinates

**Requested seed-0 sample (10 MBS-dimer, 10 MBS-isolated, 10 SBS): 30/30 agree on hit/miss for both the dimer and the isolated arm.** In the local frame, |ΔDCC| is ≤0.04 Å for all 30 sites except one, 5xf5 RUD; that site is a miss either way (see the ligand-definition bullet below).

I then extended the check to **all 4,013 primary sites**, each scored for both arms:
- **Local 15 Å frame (the producer's choice):**
  - Dimer: 3,996/4,013 agree.
  - Isolated: 3,986/4,013 agree.
  - Every disagreement but two is a ligand-definition difference, not a DCC error:
    - **31 sites:** my per-residue asymmetric-unit centroid differs from the census centroid by more than 1 Å. These are ligand UIDs with suffix `:2` (assembly copies) and covalently merged FAD in MAO-B (P27338), which the census merges by design (DESIGN §1).
    - **2 sites:** a DCC that lands exactly on 4.00 Å, so the hit/miss call is decided by floating-point rounding.
  - On the 3,933 sites whose centroids agree to within 1 Å, only those 2 boundary sites disagree.
- **Global whole-pair CA fit (my alternative):**
  - Hit/miss agreement is 96% for the dimer and 95% for the isolated arm.
  - Rates: MBS 0.371 / 0.250, SBS 0.376 / 0.421.
  - DiD **0.165**. The result is robust to the choice of superposition frame.
- **My own recompute, local frame:** MBS 0.380 / 0.263, SBS 0.396 / 0.447, DiD 0.167. **VERIFIED.**
- **Other checks:**
  - SIFTS residue-name identity was 1.0 for every superposed pair.
  - My symmetric-collapse implementation (swapped-residue Jaccard ≥0.5) gives the same hit/miss calls as the producer's on every site.

## 3. PI spot-check

Source: `reviews/verify_h4/spotcheck.tsv`, seed 0, MBS sites only, dimer arm. "Best pocket" is the pocket with the highest recall of the PDB binding residues, given as rank / size / chains / recall / DCC. Failure reasons come from my heuristic rules in `v4_failure.py`, applied in this order:
1. Local RMSD >2 Å.
2. Hit at any rank, but outside the top-k.
3. DCC 4–5 Å.
4. Best pocket recall ≥0.5: oversized if its size is ≥2× the number of binding residues, otherwise centre offset.
5. Two pockets whose recalls sum to ≥0.5: split.
6. Otherwise: no pocket.

**Hits:** all 10 confirmed (DCC recomputed: 0.79–3.89 Å).

| PDB | Ligand | UniProt pair | DCC dimer / isolated (Å) | Best dimer pocket |
|---|---|---|---|---|
| 5uau | PRO | P32322 homo | 0.79 / 1.86 | r1/17/AB/1.00 |
| 3bvj | XMP | P11172 homo | 3.54 / 4.41 | r1/20/AB/0.80 |
| 6qvl | 4DW | P34897 homo | 2.19 / 6.84 | r3/9/AB/0.50 |
| 8qxl | DTP | Q9Y3Z3 homo | 3.65 / 12.1 | r6/11/AB/0.75. The global frame gives 4.03 Å, so this hit is fragile. |
| 1bh5 | GTX | Q04760 homo | 2.02 / 4.55 | r2/32/AB/1.00 |
| 1hvy | UMP | P04818 homo | 3.89 / 1.25 | r2/11/AB/0.67 |
| 8gjx | ZNT | Q86WV6 homo (iface J 0.38) | 3.48 / 6.71 | r1/30/AB/0.88 |
| 6afc | 72R | Q99497 homo | 2.15 / 4.45 | r1/12/AB/0.71 |
| 8ros | A1H16 | Q9NVS9 homo | 2.32 / 7.01 | r1/10/AB/1.00 |
| 5qk5 | K3J | Q9UKK9 homo | 1.22 / 11.0 | r3/8/AB/1.00 |

**Misses:** all 10 confirmed as misses.

| PDB | Ligand | UniProt pair | DCC dimer / isolated (Å) | Best dimer pocket | Reason |
|---|---|---|---|---|---|
| 7uep | Y93 | Q9H999 homo | 7.83 / 9.07 | r1/**121**/AB/0.60 (11.4 Å) | Oversized merged pocket; centre displaced |
| 4fgk | 0TX | P16083 (NQO2) homo | 7.13 / 8.25 | r1/**58**/AB/1.00 | Oversized: both symmetric active sites in one pocket |
| 1xi2 | CB1 | P16083 homo | 5.75 / 8.16 | r1/58/AB/1.00 | Oversized, as above |
| 7ybu | BTI | P05165/P05166 hetero | 34.2 / 34.5 | none | **Model geometry**: local CA RMSD 9.25 Å (flexible biotin-carrier domain) |
| 4lts | LTS | P43490 (NAMPT) homo | 12.7 / **3.30 hit** | r1/17/AB/0.31 | No dimer pocket at the site (see §4) |
| 6w0n | RYS | P50053 homo | 4.97 / 4.24 | r3/16/A/0.64 | Near miss (threshold) |
| 7kdb | WBS | P31153 (MAT2A) homo | 6.60 / **2.70 hit** | r1/44/AB/0.83 | Oversized interface pocket |
| 7ppg | 7YX | P43490 homo | 13.0 / **3.72 hit** | r1/17/AB/0.30 | No dimer pocket at the site |
| 5ugh | 8AJ | P31153 homo | 7.70 / **3.13 hit** | r1/44/AB/0.80 | Oversized interface pocket |
| 7lk1 | Y37 | P04181 homo | 5.91 / 4.33 | r1/14/AB/0.43 + r3/0.29 | Site split across two pockets |

None of the 10 spot-check misses is a wrong-interface case: their interface Jaccards are 0.61–0.96. Only 21 of all 706 MBS dimer misses have a Jaccard below 0.5.

### Failure-mode census

Source: `failure_modes.tsv`, all sites, my recomputation, local frame.

| Reason | MBS dimer | MBS isolated | SBS dimer | SBS isolated | MBS dimer-miss & isolated-hit (n = 173) |
|---|---|---|---|---|---|
| Hit | 432 | 299 | 1,139 | 1,284 | – |
| Oversized/merged pocket, centre displaced | **245** | 80 | 552 | 581 | 64 |
| No pocket at site (<50% of binding residues in any pocket) | 128 | **355** | 364 | 335 | 48 (45 are NAMPT) |
| Model geometry (local RMSD >2 Å) | 124 | 126 | 119 | 125 | 15 |
| Near miss (4–5 Å) | 118 | 94 | 374 | 330 | 31 |
| Split across pockets | 45 | 114 | 50 | 60 | 6 |
| Pocket found, centre >5 Å off | 40 | 63 | 167 | 88 | 9 |
| Ranked outside top-(n+2) | **6** | 7 | 110 | 72 | 4 |

What the census shows:
- **Top-k is almost never the cause at MBS sites:** 6 of 706 dimer misses.
- **The dimer's main failure mode is the oversized interface pocket.** The pocket is found (recall ≥0.5), but its cavity centre lies far from the ligand. This is the DCC penalty for large merged pockets that the producer already flagged (caveat 2).
- **The isolated chains' main failure mode is "no pocket" or "split".** Removing the partner chain destroys or halves the site. This is the mechanism H4 hypothesises.

## 4. Heterogeneity: why the bootstrap and the permutation disagree

1. **The two tests answer different questions.**
   - The within-cluster permutation keeps the observed families fixed and asks whether, *within* families, MBS sites gain more than SBS sites. That is a fixed-effects, conditional question.
   - The cluster bootstrap asks whether the pooled, site-weighted DiD would be non-zero in a new sample of families. That is a random-effects question, and it is the pre-registered estimand.
   - The permutation null is centred at **0.091, not 0**. Only 148 of 470 clusters contain both classes; 58 are MBS-only and 264 are SBS-only, and their contribution is never permuted. The quoted p = 1e-4 therefore tests "0.164 ≠ 0.091". It does not test "DiD ≠ 0".
2. **The permutation treats near-duplicate sites as independent.**
   - The within-cluster ICC of the gain is 0.55 at MBS and 0.34 at SBS.
   - The design effects are 14 and 19, which leaves an effective n of only about 81 MBS and 149 SBS sites.
   - The permutation null SD is 0.011, against a bootstrap SE of about 0.088.
   - Even so, permuting region-level units (one per cluster × binding region × class) still gives centred p = 2e-4 (null SD 0.024). The within-family signal is therefore real; it is just not what the bootstrap measures.
3. **A few large families with extreme gains drive the bootstrap width.**
   - Leaving one cluster out moves the DiD between 0.121 and 0.213.
   - Largest contributions to the pooled DiD:

     | Family | Sites | MBS gain |
     |---|---|---|
     | Q9UKK9 | 61 | +0.89 |
     | NAMPT, P43490 | 56 | −0.88 |
     | Q86WV6 | 47 | +0.74 |
     | MAT2A, P31153 (with Q00266) | 41 | −0.73 |
     | P14618/P30613 | 35 | +0.60 |

   - SBS families also weigh in: 14-3-3 (288 sites); Q04609/Q9Y3Q0, gain −0.65 over 77 sites.
   - Dropping only the NAMPT cluster gives DiD 0.213 (CI 0.058–0.373).
4. **Family-level evidence is consistent.**
   - Of 206 MBS clusters, 68 have a positive gain and 32 a negative one (sign test p = 4e-4). SBS clusters split 95 vs 107 (p = 0.44).
   - Of the 148 mixed clusters, 67 have a positive DiD and 34 a negative one (p = 0.0013).

### NAMPT

The dimer misses here are **not a top-k artefact and not a pocket split**:
- No dimer pocket, at any rank, lies within 11.5 Å of the ligand.
- The best dimer pocket covers only 22–36% of the binding residues.
- The binding residues are on average 9.5 on one chain and 1.3 on the other. These inhibitors sit in the tunnel at the dimer interface but are lined mostly by one chain.

The clash check supports a model-geometry explanation (`v6_nampt.py`):
- After superposition, 2–14 heavy atoms per inhibitor lie within 2.5 Å of model atoms, mostly from the main lining chain. The median nearest-atom distance is 2.6–3.2 Å.
- The apo AF-Multimer tunnel is therefore narrower than the inhibitor-bound one, and burrow finds no cavity there.

The isolated-chain "hit" is also suspect:
- It comes from a 34-residue merged pocket (`_B:m0`, members 1+7, 707 Å³), rank 3, on the lining chain alone.
- With the partner chain removed, the tunnel becomes an open groove. A merged groove pocket whose centre happens to fall on the ligand is a weak form of recovery.

My reading is that NAMPT measures induced fit and DCC geometry, not interface modelling. It is a legitimate data point, but it carries 5% of the MBS weight through one target with a single inhibitor series. MAT2A (−0.73) is the oversized-interface-pocket mode, 30 sites.

### DiD under the alternative units and weights

All 1,000 cluster-bootstrap resamples, seed 0. None of these analyses is pre-registered.

| Analysis | n | DiD | 95% CI |
|---|---|---|---|
| Pre-registered (site-weighted) | 4,013 | 0.164 | −0.006 to 0.338 |
| One random site per (cluster, binding region, class)* | 1,067 | 0.110 | 0.044 to 0.174 |
| Region-averaged (each region weight 1) | 1,067 | 0.112 | 0.046 to 0.175 |
| Family-weighted (each cluster weight 1 per class) | 4,013 | 0.158 | 0.075 to 0.234 |
| Family-weighted, NAMPT cluster removed | 3,951 | 0.161 | 0.088 to 0.246 |
| Producer's exploratory region dedup | 1,708 | 0.145 | 0.071 to 0.214 |

\* Region: union-find on (accession, UniProt residue) binding sets within a cluster, with overlap coefficient ≥0.5. The producer's region definition gives a different n (1,708 vs 1,067) and a larger estimate (0.145 vs 0.110), so its magnitude is **DISCREPANT / definition-dependent**. The sign and the exclusion of 0 replicate.

## 5. Symmetric-collapse rule

Collapse is applied identically to the dimer and the pooled isolated pockets.

| Rule | MBS dimer / isolated | SBS dimer / isolated | DiD (site-weighted), CI | DiD (family-weighted), CI |
|---|---|---|---|---|
| Swapped Jaccard ≥0.5 (used) | 0.383 / 0.263 | 0.400 / 0.443 | 0.164 (−0.006, 0.338) | 0.158 (0.075, 0.234) |
| No collapse | 0.361 / 0.279 | 0.413 / 0.436 | 0.104 (−0.051, 0.260) [producer: −0.057, 0.254 — VERIFIED] | 0.154 (0.073, 0.239) |
| Overlap coefficient ≥0.5 (annotation rule) | 0.338 / 0.198 | 0.328 / 0.362 | 0.175 (0.024, 0.341) [producer: 0.014, 0.326 — VERIFIED] | 0.136 (0.056, 0.217) |

What changes and what does not:
- **Hit/miss calls do change.** Switching to the OC rule changes 371 dimer calls (315 hits lost); switching to no collapse changes 221.
- **The OC rule absorbs distinct pockets.** It lowers every hit rate (SBS dimer falls from 0.400 to 0.328), which supports the producer's diagnosis that OC is the wrong rule for top-k scoring.
- **The Jaccard rule is defensible.** DESIGN §1 says only that "pockets symmetric under chain swap are collapsed", without naming a metric. TOOLS.md §2 documents Jaccard ≥0.5.
- **The choice was not result-favourable.** OC is the only rule under which the pre-registered site-weighted CI excludes 0.
- **The primary conclusion does not change.** Under the pre-registered estimator the CI includes 0 for Jaccard and for no collapse. Under family weighting all three rules exclude 0.
- **Disclosure is still needed.** The rule was chosen after test-subset hits had been seen. Separately, the annotation's `is_sym_rep`/`dimer_only` flags, which use the OC rule, should be checked by whoever owns H1–H3/H5.

## Recommendation

1. **The headline numbers are reproducible.** Every number in docs/H4_RESULTS.md that I checked replicates. The DCC pipeline replicates from raw coordinates on 4,013/4,013 sites, apart from 31 census ligand-definition cases and 2 sites tied at 4.00 Å.
2. **Report the primary pre-registered result as not established.** Use: "DiD +0.16 (95% cluster-bootstrap CI −0.01 to 0.33, p = 0.07); H4 not established at α = 0.05".
3. **Do not cite the within-cluster permutation p = 1e-4 as support for H4.** It tests a different, conditional null, centred at 0.091 rather than 0. It also ignores the strong within-family dependence (design effect about 14–19). If it is reported, describe it as "within-family association, conditional on the observed families".
4. **Report heterogeneity explicitly.** Include the leave-one-family-out range (0.12–0.21) and the NAMPT and MAT2A cases with their mechanisms: an apo tunnel too narrow for the inhibitor, and oversized interface pockets. Add the family sign test (68 vs 32 MBS families improve) and the family-weighted and region-dedup DiDs (0.11–0.16, CIs excluding 0), all **labelled exploratory**.
5. **Add the failure-mode table as a secondary result.**
   - Dimer models rarely miss MBS sites through ranking.
   - The dominant dimer failure is the oversized interface pocket, a DCC-to-cavity-centre limitation; residue recall at MBS was 0.67 for dimer vs 0.40 for isolated.
   - The dominant isolated-chain failure is that the site is absent or split.
   - This mechanistic contrast is better supported than the pooled DiD p-value.
6. **Disclose that the collapse rule was chosen post hoc,** and show all three rules.
7. **Caveats on this verification.**
   - The failure-mode categories are my heuristics and have not been visually inspected.
   - The NAMPT clash interpretation rests on 7 sites.
   - n (k), the census ligand definitions and the interface filter were taken from the producer's tables without re-derivation.
