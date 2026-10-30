# R2: statistics and external PI (DESIGN v2)
Reviewer: combined biostatistician and structural-biology PI, 2026-10-29. No analyses run.

## (1) Round-1 musts: how they were adopted
Most were adopted well. Four were mis-implemented or left out:
- **Interface pocket**: r1-stats asked for ≥3 residues per chain as the primary rule. v2 makes ≥1 primary and ≥3 strict, which is the reverse.
- **Homodimer chain-swap**: r1-compute's must (score recovery under both A/B assignments and take the max) is not in Q2.1.
- **Dimer-only overlap threshold**: the 0.2/0.5 sensitivities are not stated.
- **Seeds**: the primary seed is not named.

## (2) Hypotheses
- **H1**: the unit is cluster × CCD, but the covariates (resolution, year, method) are entry-level, and a unit can be MBS in one entry and SBS in another. *Must*: use one representative entry per unit (best resolution, ties broken by earliest ID) and code MBS from that entry. Report "any-MBS" as a sensitivity.
- **H2**: specify the model as `drug_like ~ MBS + strata(cluster)`. Only clusters that are discordant on both variables are informative, so apply the power gate to discordant clusters. Call the outcome "drug-like occupancy", not druggability, because it reflects PDB deposition bias.
- **H3**: residue-Jaccard misclassifies asymmetric ligands that sit on the axis. *Should*: add a geometric sensitivity (ligand centroid ≤3 Å from the C2 axis, taken from the A→B superposition). Restrict the test to drug-like vs cofactor. It will probably be descriptive (few drug-like on-axis clusters), so state that up front.
- **H4**: the interface-correspondence filter conditions on the dimer being correct, so the estimand is conditional. *Must*: state it as "given a correct interface", and report the unfiltered set and the filtered-out count. Define n (the number of PDB sites on that chain pair). Define isolated-chain top-k as the union of both chains' pockets ranked by burrow_score. Use a cluster-stratified paired test (McNemar or conditional logistic). As specified, H4 is close to tautological for MBS. *Should*: make the primary contrast the dimer gain at MBS minus the gain at SBS (an interaction term).
- **H5**: as written it is confounded by burial, because pocket residues are more buried and burial alone predicts pathogenicity. *Must*: use residue-level logistic GLMM with a protein random intercept and covariates RSA (dimer), pLDDT and conservation. Make the categories mutually exclusive with a stated precedence: dimer-only > single-chain pocket > interface > surface. One primary contrast (vs interface non-pocket) and one primary outcome (ClinVar P/LP); AlphaMissense is secondary.

## (3) New risks
- BioLiP2 or ChEMBL-pChEMBL may not be fetchable by tomorrow. *Must*: pre-declare the fallbacks now: SoI flag + geometric rule; max_phase only, using the UniChem bulk file plus the ChEMBL molecule table. Mark them FALLBACK, not silent.
- *Should*: drop the ≥2 Å versus <1 Å inconsistency for symmetric dedup by using 2 Å (prior-work r1). Name one primary druggability descriptor (burrow_score) and list the rest as secondary.
- Deep dives: rank only BH-surviving primary results and CI-excluding-null exploratory results, labelled as such.

Feasibility by 13:00 otherwise looks fine. No other changes needed; freeze after these edits.
