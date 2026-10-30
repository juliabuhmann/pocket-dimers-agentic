# Report review — statistics/integrity (agent; saved by the PI from the agent's returned text, 2026-10-30)

The reviewer spot-checked 31 numbers against docs, verifier files and TSVs (hypothesis_tests.tsv, h4_tests.tsv, seeds_summary.tsv). Most match.

MUST
1. Use the verifier's co-folding counts (SBS control 3/6 vs 2/6; dimer-only single chain 6/8; verify_cofold.md:62).
2. "Most cases pre-cutoff" is wrong for the on-axis batch: 20/24 are post-cutoff (COFOLD_ONAXIS.md:55).
3. The MBS fractions (26.3/30.6/29.2/34.9%, lipid/heme) are tagged VERIFIED, but only the unit counts were re-derived → UNVERIFIED.
4. H3 is "descriptive" but carries BH p 3.5e-11 and is counted in family A. Without H3, H2's BH p = 0.017. The per-cell gate reading was post hoc; disclose this.
5. Add the unfiltered H4 (DESIGN §1): DiD 0.136 (−0.006–0.292). The interface filter removes 24% of homo MBS units vs 6.5% of SBS.
6. H4 collapse rule: the overlap-coefficient rule gives 0.175 (0.014–0.326), significant. The verdict depends on a post-hoc rule.
7. Report family-B BH p-values (H5 ≤0.0048; H4 ~0.07).
8. On-axis co-folding: McNemar exact p = 0.125 (4 vs 0); RMSD sign test p ≈ 0.006 (11/12). 2 of the 4 dimer-only successes are pre-cutoff and 2 are ligand-induced dimers (PD-L1, XIAP). Soften "needs the partner".
9. Resolve the seeds placeholder (only on-axis dimer seed 1 exists).
10. DD selection: DD2's primary CI includes 0. Show the ranked list, or say DD2/DD3 were chosen for interpretability.

SHOULD
11. Report the pre-registered AlphaMissense secondary for H5 (+0.025; 0.700 vs 0.699).
12. Missing deviations: H4 bootstrap chosen over permutation; the H3 gate reading; NAMPT run with 1 ligand copy.
13. DD1: drop the "nature vs drug discovery" framing; give base rates (off-axis 334 drug-like vs 576 cofactor clusters); lead with the cluster-level OR ≈14.
14. ~1.4× of the enrichment is intrinsic; "the partner causes" applies only to the property shifts.
15. Add H4 mean F1 (0.371 vs 0.311, atlas 0.296); recall favours merged pockets.
16. "68 vs 32 families" = the raw MBS gain (106 ties), not the DiD; SBS 95 vs 107.
17. The AUC of 0.87 is partly mechanical; remove the druggability-signal claim.
18. The cofactor half-site claim rests on ~12 cases; STING1 is n = 1; per-batch ρ CIs include 0.
19. H5 pLDDT dependence: 1.61 at 70–90, 0.99 at ≥90; the verifier rates it FRAGILE.
20. The curated list is filtered on non-overlap with the atlas, so "atlases miss" is partly by construction.
21. Minor: 351 vs 352; H3 leave-out-top-10 18.8 (producer) vs 9.0 (verifier); stratified H4 CI; the ρ synthesis is not re-derived.
