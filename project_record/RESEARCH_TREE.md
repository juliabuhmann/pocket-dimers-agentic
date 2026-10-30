# Research tree — Pockets in dimers vs monomers

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Start: see START.txt (2026-10-29T13:24+01:00). Budget revised by user: deadline 2026-10-30T13:00+01:00 (~23.5 h). Planned milestones: M1 design frozen ~17:30; M2 Part1 PDB core ~23:00; M3 Part2 core ~05:00; M4 co-folding/deep dives ~10:00; M5 final report 12:30.

## Status summary (latest milestone on top)

### M7 — FINAL (2026-10-30 ~07:10, t+17.8h; deadline 13:00)
The stalled co-folding runs were completed from a node-local env copy (all 3 seeds; 26 heterodimer pairs). A final fresh read-through re-derived these numbers (VERIFIED) and found that the heterodimer success gain is concentrated in pre-cutoff pairs (10/13 vs 5/13; post 4/13 vs 3/13), and that the partner-share link from the 56-case synthesis does not replicate in the 26 heterodimer pairs (ρ 0.09). REPORT.md was hedged accordingly.
Deliverables: REPORT.md (v3, final), README.md, DESIGN.md (frozen), LITERATURE.md, DATA.md, LESSONS.md, RESEARCH_TREE.md, reviews/ (7 design reviews, 1 audit, 8 verification reports, 2 report reviews, 2 retrospectives). All on git branch `pocket-dimers-study` (not pushed; no credential).

### M6 — REPORT v3 verified (2026-10-30 04:25, t+15h; 8.5 h to deadline)
All v3 additions were independently verified (reviews/verify_v3.md): uncropped PDB descriptors, fpocket druggability, LIGYSIS calibration, seeds and the hetero extension, and also the Part 1 MBS fractions, which were previously unverified. REPORT.md tags were updated, and three small count and wording errors were fixed. Reproducibility record completed: DATA.md checksums and licences, envs/ lock files for 12 environments, docs/SEEDS.md. LESSONS.md received its final prune from a fresh-agent retrospective.
Scientific bottom line: multi-chain sites are not systematically more or less druggable (fpocket, burrow_score null; drug-like occupancy only weakly higher). The clearest biological difference is geometric: twofold-axis homomer sites are drug-dominated, while the major cofactors are single-chain-anchored. Predicted dimers help find multi-chain sites (family-level), and the partner matters in co-folding mainly for evenly split sites.
Decision: use the remaining time only to complete the stalled co-folding runs (node-local env) and a final fresh read-through; no new analyses (LESSONS rule 13).

### M5 — report v2 after review; gap-filling launched (2026-10-29 23:55, t+10.5h; 13 h to deadline)
REPORT v1 was reviewed by an external PI and a statistics/integrity reviewer (reviews/report_review_*.md). The main corrections, all applied in REPORT v2:
- H1 reframed: Marsh's cofactor > metal pattern is reproduced; cofactor-vs-drug-like is a new contrast and null.
- Verifier co-folding counts used; exact tests for small-n co-folding (McNemar p 0.125, RMSD sign test p 0.006).
- Unfiltered H4 result added; the H4 verdict's dependence on the post-hoc collapse rule disclosed (the alternative rule would be significant).
- H3 gate reading disclosed as post hoc.
- DD2/DD3 selection deviation logged.
- Over-claims removed.
The independent misc. verification confirmed the PDB-descriptor, Q2.4, AUC and curated-list numbers, but found that the 30 Å crop matters → uncropped rerun in progress.
Decision: use the remaining time to close the report's biggest gaps rather than open new questions: (1) a real druggability predictor (fpocket) on PDB MBS/SBS sites and the dimer-only list; (2) heterodimer MBS co-folding expansion (n=6 → ~26) on GPU 1; (3) LIGYSIS calibration of BURROW; (4) seeds 1–2 replicates (GPUs 2–3). v3 at ~06:00; final checks and a last LESSONS retrospective before 12:00.

### M4 — analyses complete, report drafting (2026-10-29 23:20, t+9.9h; 13.7 h to deadline)
All planned analyses except the H4 co-folding batch are finished. Added since M3: the PDB pocket descriptors (real MBS sites are only modestly larger and NOT more charged than SBS, which corrects the predicted-pocket chemistry story); the failed "charged = ligand-free groove" hypothesis, with burrow_score separating ligand-matched interface pockets (AUC 0.87); cross-species conservation (interface pockets slightly more conserved, heterodimers only); dimer-only active-site co-folding (cofactor half-sites fold without the partner; STING1, an even-split site, needs it).
Decision: start the REPORT.md draft now with placeholders. Run seeds 1–2 replicates on GPUs 2–3 (stop by 03:30). After the H4 co-folding batch, a final verification sweep (co-folding RMSDs across batches, PDB descriptors, Q2.4) by an agent that didn't produce them, plus a synthesis analysis: partner dependence of co-folding vs the partner share of site residues. Then a report review (external PI + statistics) and the final revision before ~10:00.

### M3 — Part 2 core verified, deep dives chosen (2026-10-29 20:17, t+6.9h)
All pre-registered hypotheses are now tested and independently verified. Family A: H1 null (Marsh's cofactor enrichment not replicated); H2 weak/fragile; H3 strong but descriptive (drugs on the twofold axis). Family B: H4 not established on the primary family-bootstrap test (DiD +0.16, CI −0.01–0.34), but consistent across families in exploratory family-weighted analyses (68 vs 32 families). H5 passed (OR 1.20) but is a general, composition-driven pocket effect, not dimer-specific. Exploratory: interface pockets in predicted dimers are larger, more charged and less hydrophobic, caused by the partner (control verified); their enrichment is partly surface geometry.
Decision (pre-declared selection rule DESIGN §4 applied: effect size with CI excluding null, then interpretability): Deep dive 1 = drugs on the homomer twofold axis (H3 + co-folding on/off-axis); Deep dive 2 = partner-dependence of multi-chain sites: dimer models vs isolated chains, failure modes, and co-folding rescue (H4 exploratory + co-folding); Deep dive 3 = physicochemistry of interface vs single-chain pockets in predicted AND experimental structures, with cofactor anchoring (H1 strict reversal) and the disease-variant null (H5). Total comparisons so far: 5 pre-registered tests + ~40 sensitivity/exploratory comparisons (counted in the final report). Running: PDB-site pocket descriptors (validates DD3), Q2.3a curation, 2 co-folding batches.

### M2 — Part 1 core verified (2026-10-29 18:54, t+5.5h)
Part 1 is complete and independently re-derived. Marsh's cofactor enrichment at multi-chain sites is NOT replicated against drug-like ligands with modern data (H1 null; under strict MBS the major cofactors actually prefer single chains). Drug-like occupancy at MBS vs SBS is only weakly higher within families (H2, fragile). The striking pattern is symmetry: ligands bound on the homomer twofold axis are overwhelmingly drug-like (64 families incl. TTR, PKM2, IDH2, STING, PD-L1, glutaminase). This is reported descriptively because the on-axis cofactor cell is tiny, i.e. nature rarely puts cofactors on the axis while drug discovery often does.
Decision: deep dive #1 candidate = on-axis drug sites; it gets literature context (stabilisers, dimerisers, allosteric activators) and a co-folding check (do AF3 dimer+drug predictions place the drug on the axis, and the monomer+drug predictions not?). Part 2 analyses (H4, H5, Q2.2) are still running.

### M1.5 — production data ready, analyses launched (2026-10-29 ~18:00, t+4.6h)
Human dimer pocket runs are done (10,610 models, 87,804 pockets, 0 failures). Interface pockets (≥3 residues per chain) are 29% of homo and 27% of hetero pockets; ~64% of interface pockets are dimer-only (absent in the isolated chains). UNVERIFIED. Other species are running in the background (~3.5 h). Co-folding feasible: AF3 (xla attention) runs on GPUs 1–3, ~3 min at 500 tokens, ~12 min at 1,000, ≤1,280 tokens without unified memory; GPU 0 corrupts transfers (all-NaN), so it is excluded. The Part 1 census audit led to census v2 (running).
Decision: launch H4 (recovery), H5 (disease variants) and Q2.2/Q2.3a (dimer vs monomer, dimer-only ligandable list) in parallel now; H4 is built on census v1 and finalised on v2. Co-folding target selection waits for the first results (deadline for selection ~23:00).

### M1 — design frozen (2026-10-29 ~17:00, t+3.6h; deadline 2026-10-30 13:00)
DESIGN.md v3 was frozen (tag `design-v3-frozen`, sha256 in results/DESIGN_FROZEN.sha256) after two review rounds with 7 reviewer reports (reviews/). It has five pre-registered hypotheses. Family A (PDB): H1 cofactor vs drug-like enrichment at MBS; H2 within-cluster drug-like occupancy MBS vs SBS; H3 on- vs off-axis homomer sites. Family B (predicted): H4 dimer vs isolated-chain recovery of PDB MBS sites (difference-in-differences vs SBS); H5 ClinVar P/LP enrichment at dimer-only pocket residues vs interface non-pocket residues. The reviews notably fixed biased recovery baselines, novelty self-matching leakage and kNN leakage, and added the symmetry and disease angles.
Running: Part 1 census (updated with BioLiP2 labels, symmetry dedup, glycan/covalent/antibody handling); dimer pocket production runs; external data fetch (humsavar, AlphaMissense, ChEMBL, SIFTS, LIGYSIS, PLINDER test); co-folding feasibility (AF3 shows NaNs on Turing, Boltz-2 being timed). Decision: ClinVar P/LP is implemented via UniProt humsavar (curated from ClinVar, UniProt numbering), which avoids genomic→protein mapping; logged as the implementation of the H5 outcome source.

### M0.5 — recon done, design v1 in review (2026-10-29 ~17:05, t+3.7h)
Reconnaissance finished but ran long (tools agent 3.1 h). Key facts: the human AFDB dimer set (4,819 homo / 5,791 hetero) is pre-filtered to confident interfaces; the pocket atlas covers monomers only, so dimer pockets must be generated by us with identical settings (pocketcompare 1.0.0 matches the atlas, verified on P31946). burrow handles multi-chain inputs; pocketcompare scores are uncalibrated and weak for PPI pockets; adock GPU kernel hangs on this node (CPU OK); GPU 0 is faulty. No tool provides druggability.
Decision: to recover time, Part 1 PDB census and dimer pocket production runs were launched in parallel with design review, since both are needed under any design (thresholds kept configurable). DESIGN.md v1 is in review round 1 (5 reviewers: external PI, prior work, statistics, compute, AI failure modes).

### M0 — kickoff (2026-10-29 13:30, t+0.1h)
Project set up: workspace git repo (branch `pocket-dimers-study`), PROMPT.md saved verbatim. Hardware: 64 CPUs, 501 GB RAM (~425 GB available), 4× Quadro RTX 6000 24 GB (Turing, no bf16 — relevant for AF3). Workspace filesystem has 7.1 TB free; our cap is 500 GB.
Four parallel reconnaissance agents launched: data mounts, tools (burrow/pocketcompare/adock), literature (incl. Marsh/Abrusán 2018), co-folding feasibility. Mount note: `/cluster/project/beltrao/gmagnus/dimers` and `/cluster/scratch/beltrao/gmagnus/dimers` are NOT mounted; the user confirmed `/cluster/scratch/gmagnusson/dimers` is the dimer source.

Decision: do reconnaissance before design; design (DESIGN.md) only after LITERATURE.md + data/tool facts are in, then ≤3 review rounds.

## Tree (cleaned 2026-10-30 04:15; details in the milestone notes above and in docs/)
- [x] 0. Setup & reconnaissance: data (docs/DATA_RECON.md), tools (docs/TOOLS.md), literature (LITERATURE.md), co-folding feasibility (docs/COFOLDING.md), external data (docs/EXTERNAL_DATA.md)
- [x] 1. Design: DESIGN.md v3 frozen after 2 review rounds (reviews/r1_*, r2_*)
- [x] 2. Part 1 PDB census
  - [x] 2.1 Census v1 → audit (reviews/audit_part1_census.md) → census v2 (docs/PART1_PIPELINE.md)
  - [-] H1 cofactor vs drug-like at MBS: null, OR 0.89 [VERIFIED]; Marsh cofactor > metal pattern reproduced; major cofactors single-chain-anchored (strict OR 0.41)
  - [~] H2 drug-like occupancy MBS vs SBS: OR 1.17, fragile [VERIFIED]
  - [!] H3 on-axis homomer sites drug-dominated: cluster OR ≈14, descriptive [VERIFIED] → DD1
  - [x] 2.4 PDB pocket descriptors, uncropped: burrow_score δ −0.01, volume +0.07, chemistry null [verify_v3 pending]
  - [-] 2.5 fpocket druggability MBS vs SBS: null (δ +0.01) [verify_v3 pending]
- [x] 3. Part 2 predicted dimers
  - [x] 3.0 Pocket runs: human dimers + isolated chains + 12 species
  - [-] H4 recovery: DiD +0.16, CI includes 0 (primary not established); family-weighted 0.158 [VERIFIED] → DD2
  - [!] Q2.2 interface pockets: 1.9× enriched, partner-caused property shifts [VERIFIED]; curated dimer-only list = enzyme active sites [EXPLORATORY] → DD2/DD3
  - [-] Charged groove hypothesis not supported; burrow_score AUC 0.87 for ligand-matched pockets [VERIFIED, partly mechanical]
  - [~] H5 variants: OR 1.20 passed, but not dimer-specific, composition-driven [VERIFIED] → DD3
  - [~] Q2.4 conservation: +0.030, heterodimers only [VERIFIED]
  - [x] LIGYSIS calibration: BURROW 57.6% vs fpocket+PRANK 60.4% at DCC 12 Å [verify_v3 pending]
- [x] 4. Co-folding (AF3)
  - [!] On-axis drugs: dimer 7/12 vs single chain 3/12 (McNemar p 0.125; RMSD sign p 0.006); stable over seeds [RMSDs VERIFIED]
  - [x] Dimer-only active sites 6/8 vs 6/8 (verifier metric); STING1 needs the partner at every seed
  - [x] H4 sites 11/24 vs 10/24; heterodimer extension 11/21 vs 7/21 (McNemar p 0.22; RMSD sign p 0.0015)
  - [x] Synthesis: partner RMSD gain ~ partner share ρ 0.34 [EXPLORATORY]
- [x] 5. Deep dives DD1–DD3 in REPORT.md §6
- [x] 6. REPORT v1 → reviews (PI, statistics) → v2 → v3 (additions)
- [x] 7. Final: verify_v3, DATA.md checksums + env locks + SEEDS.md
- [x] 8. Complete stalled co-folding runs (seed 2, 5 hetero cases); final read-through
