# Independent verification of Q2.2 / Q2.3a (docs/Q22_RESULTS.md)

Status: first-pass verification by a second agent; not reviewed by a human. Code: `reviews/verify_q22/v1_features.py` (SASA, interface, shifted models),
`v2_enrich.py`, `v3_props_dimeronly.py`, `v4_inspect.py`. Outputs: `enrich_out.txt`, `v3_out.txt`, `v4_inspect.tsv`. No `scripts/q22_*.py` was imported.

## Method
- Sample: 1,000 of the 10,588 models (`DataFrame.sample(1000, random_state=0)` on distinct modelEntityId). BURROW pockets are the existing ones (atlas_pass: burrow_score >=0.03, lining pLDDT >=70, all copies).
- Residue SASA: freesasa Shrake-Rupley on heavy atoms with my own element radii, per chain alone; rSASA = SASA / Tien-2013 theoretical max. Surface = rSASA >=0.2 (0.25 also run). Interface = heavy atom <=5 A from a partner-chain heavy atom (complex coordinates).
- Interface pocket (their definition): >=3 lining residues on each chain. My alternative (usable in the control): >=3 lining residues, >=50% of them interface residues.
- Control: partner chain B translated along the A->B centroid axis, starting at 20 A and stepping +5 A until min inter-chain distance >=15 A (a plain 20 A shift left a median 5 A gap, i.e. chains still touching, so I enforced the gap; shift median 40 A). BURROW (`pocketcompare describe --af-mode`, same settings as the atlas) was re-run on the 1,000 shifted models (7,289 pockets). "Former interface" = the residues that were interface in the intact complex.

## (1) Enrichment: replicated, but a geometric/concavity component exists
| surface cut | | pocket-level (their def) | pocket-level (my def) | residue-level |
|---|---|---|---|---|
| 0.2 | real | **1.88** | 1.76 | **2.43** (41.3% of lining surface residues are interface vs 17.0% expected; 27.1% vs 7.9% rate) |
| 0.2 | partner removed | 0 (cannot be formed) | 0.95 | **1.39** (rate 9.2% vs 6.1%) |
| 0.25 | real | 1.88 | 1.77 | 2.57 |
| 0.25 | partner removed | 0 | 0.95 | 1.46 |

- The headline numbers reproduce (1.88 pocket-level exactly; 2.43-2.57 vs 2.51 residue-level).
- With the partner removed, pocket-level enrichment under my definition falls to ~1 (0.95), as expected if the real excess comes from pockets closed by the partner. The residue-level ratio only falls to 1.39-1.46: former interface residues keep a 1.5x higher pocket-lining rate even alone. Some of the interface enrichment is therefore intrinsic surface geometry/composition (concave, hydrophobic, conserved patches), not partner-created cavities. Roughly 70-75% of the residue-level excess (2.43-1 vs 1.39-1) disappears with the partner; the excess attributable to the partner is ~1.75x (2.43/1.39), not 2.5x. The doc's own caveat about concave grooves is correct and this puts a number on it. The control is imperfect (the monomer is re-detected without the partner, so some interface patches are exposed as flatter/convex-in-isolation), so this is a rough bound.
- Wording "interface pockets are 1.9x/2.5x enriched" should be qualified with ~1.4x intrinsic.

## (2) Property differences: replicate, and are not an artefact of region
Interface vs single-chain (my sample, sym reps, atlas_pass; 1,718 vs 2,852 pockets), Cliff's delta: volume +0.224 (their +0.23), n lining residues +0.32, burrow_score +0.06 (+0.09), hydrophobic fraction -0.134 (-0.15), aromatic -0.11, charged +0.187 (+0.20), pLDDT -0.00. Medians: volume 102 vs 77 A3 (their 101/77).
Control (shifted; pockets lined by former-interface residues vs pockets with none): volume -0.06, size -0.26, charged -0.12, hydrophobic +0.10. The signs reverse, so the interface effects are due to the partner contact rather than to those surface regions as such (the monomer-only former interface patches are smaller, more hydrophobic).

## (3) Dimer-only fraction and noise baseline
- Full data recomputation from `max_overlap_iso`: interface (>=3/chain), atlas_pass, sym rep, no marginal: n=17,882, OC<0.5 fraction **0.635** (exact). Without the sym-rep filter 0.619. (I did not independently check the stored OC column in full.)
- My own overlap coefficients on the 1,000 sample vs isolated-chain pockets: interface new = 0.645-0.676 (vs all / atlas-pass isolated pockets), so the fraction is sensitive by ~3 points to which isolated pockets count.
- Noise baseline: single-chain dimer pockets with no isolated counterpart: 0.19-0.205 in my sample (all models, no chain-swap matching) vs their 0.166 (homodimer-chain-swap matched). Same ballpark; the claim "16.6%" is plausible but my number is 3-4 points higher. Reverse direction: isolated pockets with zero interface residues lost in the dimer: 39% (their 34%), with >=3 interface residues 44.5%.
- Score dependence is confirmed and strong: single-chain "new" rate 41.5% in the lowest score quintile vs 3.6% in the highest; interface 86% vs 31%. The interface excess persists in every quintile (largest relative gap in the top quintile, 31% vs 3.6%), which supports the claim that the signal survives noise, but the lowest-quintile interface rate (86%) is mostly low-score pockets at ~40% noise.

## (4) Q2.3a entries (14 inspected: RGS4, MYCN, MAX x2 plus a seed-0 random 10; `v4_inspect.tsv`)
Columns checked: lining pLDDT, interface residues per chain, helix fraction (phi/psi), buried area, chain disorder.
- Plausible: RGS4-GNAT2 (1,134 A3, lining pLDDT 96, 48+13 residues, ipSAE 0.87, BSA 2.5k A2). Mostly an RGS-box concavity (partial counterpart in the monomer, OC 0.43) that GNAT2 completes. A real PPI groove, but a 1,134 A3 groove at an interface is a binding surface, not obviously a small-molecule pocket. AMIGO1 homodimer, VEGFA homodimer (lining 96-98, structured), TMEM67, SLC13A4 TM region, NAMPT (large BSA, ordered) are structurally well-supported interfaces.
- Likely coiled-coil/groove artefacts: MYCN-MAX and MAX homodimer (bHLH-LZ: lining helix fraction 0.79-0.82, 23+5 and 10+7 residues, the partner contributes 5-7 residues; MYCN chain 81% pLDDT <70, the model is ordered only over the HLH-LZ). The pocket is the groove between two helices of a zipper, which is what a cavity detector will find in any helical dimer. Note MYCN and MAX are listed as separate proteins sharing one pocket (ranks 7/8), and MAX-3 isoform gives a third entry, so the "proteins" count double-counts.
- Weak: H3-4/H2A (ipSAE 0.61, lining pLDDT mean 76, min 64; histone fold, borderline), PRLR (score 0.07, lining pLDDT min 62), SLC1A3 (3 residues on one chain, score 0.08), NAMPT (3 residues on chain A, score 0.05), ITGB8 (score 0.04), TRIM2 (score 0.07). Several have burrow_score 0.04-0.08, i.e. in the quintile where the null "new" rate is 20-40%.
- Overall in the 14: ~5 convincing, ~3 helical-groove artefacts, ~6 low-score or low-pLDDT marginal. The "drug-sized" criterion is volume only (300-1500 A3); it does not demand buriedness or a non-helical lining.

## Bottom line
Headline numbers (1.88, 2.5, volume/charge/hydrophobic deltas, 0.635) reproduce. Enrichment is overstated as a partner-induced effect: ~1.4x residual remains with the partner removed. Dimer-only noise floor (~17-20%) and score-dependence confirmed. Q2.3a list contains zipper-groove and low-score entries and should be filtered (burrow_score, non-helical lining, pLDDT).
