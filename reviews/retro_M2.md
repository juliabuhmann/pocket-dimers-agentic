# Process retro at M2 (2026-10-29 ~19:20, t+5.9h; deadline 2026-10-30 13:00, ~17.7 h left)

Based on RESEARCH_TREE.md, LESSONS.md, docs/COMPUTE_POLICY.md, logs/compute_ledger.txt, git log, reviews/audit_part1_census.md, reviews/verify_part1_results.md. No analyses were run.

## 1. Stop doing
- Launching analyses on unaudited labels. v1 census drug_like was 75% Ro5/QED fallback (buffers, cholesterol) and 43% of symmetric copies were uncollapsed (audit H-1, H-2). Only the audit-before-analysis step saved H1-H3. Keep v1-style speculative launches (H4 on census v1, tree M1.5) to cases where rerun cost is small.
- Open-ended agent budgets. Tools recon ran 3.1 h vs ~1 h planned (LESSONS 3); every agent still needs a stated wall-clock cap in its prompt.
- Writing "UNVERIFIED" results into the tree as prose headlines. Tree entries for H4, Q2.2, H5 are all unverified, yet a deep dive (cofold on-axis) was launched from H3 before the H3 caveats were absorbed (H3 effect is between-protein, 127/351 on-axis units are HIV protease, within-cluster sign test p=0.38).
- Ledger sloppiness. Several entries have no DONE (dimer pocket run 16:46, ChEMBL download, Boltz timings, q22 GPU) and the 19:13 cofold entry lists "~250 per GPU process" threads against a 48-worker/2500-thread policy. Re-read ledger before every launch; close entries.
- Silent fallbacks that contradict lessons: Q2.2 agent stopped the ChEMBL SQLite query "for IO" and used REST (LESSONS 7 says REST is 6-15 s/molecule). Any fallback must be logged in the tree as a deviation.
- Adding new exploratory branches (Q2.3a MYC/MAX, 3.3 druggability, other species) while 3 verifications and the report are outstanding.

## 2. Continue
- Independent audit then independent re-derivation: audit found 3 HIGH issues in 10 min; verifier reproduced 100% of headline numbers (H1 0.89, H2 1.17, H3 35.4) without importing producer code.
- Pre-registered design with 2 review rounds (7 reports), frozen hash. Deviations are logged with "decided before results seen" (census v2 note).
- Honest negative/fragile reporting: H1 null, H2 fragile, H4 primary p=0.066 reported as not significant, stress tests (viral exclusion, leave-top-k-out) before headlines.
- Parallel agents with ledger and thread policy; no crash since policy. Running ahead of plan: M2 at 18:54 vs planned 23:00.
- Commit discipline: small functional commits, summary tables only, "(unverified)" in commit messages.
- Known-answer tests in pipelines (12/12, 16/16).

## 3. Top 3 risks to final deliverable
1. Verification debt. REPORT needs "verified results", but only Part 1 (H1-H3) is verified. H4 (final v2 just finished 19:19), H5 (GLMMs still being rerun, a label bug already found 18:54), Q2.2/Q2.3a (exploratory; burrow unstable: 16.6% of non-interface pockets are also "new") are unverified. Each needs an independent verifier (~30-60 min) plus fixes; verifiers must start by ~22:00 or the report will carry unverified numbers.
2. Deep dives are on thin or confounded ground. Deep dive #1 (on-axis drug sites) is a between-protein effect dominated by HIV protease/TTR; the AF3 co-folding test (24 cases x dimer/monomer, seed 0, one seed, <=1,280 tokens, GPUs 2-3 only after GPU 0 failure) may be uninformative if AF3 places ligands on the axis by memorisation of PDB training data (cases are PDB-known). Dives #2/#3 are not chosen (H4 DiD CI includes 0, family heterogeneity NAMPT -0.88; MYC/MAX list needs curation). Without a decision on dives 2 and 3 by ~23:00, writing time collapses.
3. Endgame time and integration. REPORT.md not started; the tree shows [ ] for 3.3, 3.4, 4. Overnight agents risk the thread-limit crash again, GPU runs may fail silently (LESSON 6), and ~12:30 final-report target leaves 30 min slack. Narrative also risks mixing a null H1, a weak H2, a non-significant H4 and an inflated "dimer-only 64%" figure without a coherent story; multiple-testing across H1-H5 and exploratory items must be stated.

## 4. Proposed schedule
- Now-22:00: verifiers for H4, H5, Q2.2 in parallel; cofold results triaged; choose dives 2-3 (candidates: H4 recovery gain at homomer MBS; MYC/MAX or another dimer-only ligandable pocket; H5 if significant).
- 22:00-04:00: dive 2/3 analyses + co-fold checks (multiple seeds, report PDB-release-date leakage).
- 04:00-08:00: REPORT.md skeleton filled with verified numbers only; a verifier reads the report against tables.
- 08:00-11:30: fixes, figures, final audit. 11:30 freeze; 12:00 final commit; 13:00 deadline.

## 5. Pruned LESSONS.md (12 rules)
See LESSONS_proposed_M2.txt appended below.

1. Check mounts against the prompt before planning; two prompt dimer paths were unmounted (M0).
2. Thread budget binds, not CPU: ulimit -u 5120, source scripts/env_threads.sh, <=48 workers, log in the ledger and close with DONE (a prior crash; M0).
3. Every agent gets a wall-clock budget and stop condition; tools recon ran 3.1 h vs 1 h (M0).
4. Never `git add` bulk results dirs; commit summary tables only (3,378 files stalled git; M0).
5. Use GPUs 1-3 only; GPU 0 gives all-NaN and cost ~3 h of debugging; run scripts/test_gpu_basic.py first (M1).
6. Check background job logs within minutes of launch; the atlas snapshot died unnoticed for 40 min (M1).
7. Use local ChEMBL SQLite, not REST (6-15 s/molecule); log any fallback as a deviation (M1, Q2.2 reverted to REST at M2).
8. Audit labels and units before any hypothesis analysis; the audit found 3 HIGH errors (75% Ro5 fallback drug-like, 43% uncollapsed symmetry, fallback additive mismatch) (M1).
9. Re-derive headline numbers with an independent verifier that does not import producer code; Part 1 reproduced 100% (M2).
10. Run stress tests before headlining an effect: H3 OR 35 falls to 5.8 without top-20 clusters and is between-protein (within-cluster sign test p=0.38) (M2).
11. Report primary results as pre-registered even when null/non-significant (H1 0.89; H4 p=0.066) and log deviations with "decided before results seen" (M1-M2).
12. Do not launch downstream work (co-folding, deep dives) on unverified, confounded results without stating the caveat in the agent's prompt; tag results UNVERIFIED until a verifier file exists (M2).
