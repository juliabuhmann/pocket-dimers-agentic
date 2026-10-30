# Independent audit of the Part 1 PDB ligand-site census

Auditor: independent agent (did not build the pipeline). Date: 2026-10-29. Scope: DESIGN.md §0–2, docs/PART1_PIPELINE.md,
results/part1/*.tsv, data/part1/{ligand_sites,entries,ligand_classes}.parquet, raw mmCIF.
Nothing in the pipeline code or data was modified. The audit scripts are in `reviews/audit_part1_scripts/` (t1–t6; run from the
workspace root with `.venvs/part1/bin/python reviews/audit_part1_scripts/tN.py <tmpdir>`; t1 must run first because it writes `ls.pkl`).
Every number below comes from those scripts. Sample sizes are small where noted, so the estimates are rough.

## Summary
The core machinery checks out:
- The headline MBS fractions reproduce exactly.
- An independent contact calculation agrees on 40/40 instances.
- The H1 representative rule reproduces in 85,208/85,208 groups.

The problems are in the **labels and the units of analysis**, not in the geometry. Three of them change what H1–H3 would measure.
They are ranked by severity below.

---

## HIGH severity

### H-1. The "drug_like" class is mostly Ro5+QED FALLBACK compounds, including buffers. ChEMBL pChEMBL ≥ 6 is not implemented.
**Evidence (t3b, t3c).**
- In the primary set there are 58,740 drug_like instances:
  - 44,230 (75%) from `ro5qed_FALLBACK`;
  - 9,001 from `chembl` only;
  - 5,509 from `chembl+ro5qed`.
- Among H1 representatives, 22,743 of 26,137 (87%) drug_like reps are Ro5/QED-only.
- The most frequent drug_like CCDs are mostly not drugs:
  - buffers: MES (3,674), EPE/HEPES (2,332), NHE/CHES, MPO/MOPS, CXS/CAPS, 1PS;
  - metabolites and lipids: CLR cholesterol (4,499, top-1 drug_like code), CIT citrate, CHD cholic acid, IHP, U10;
  - FOL folate; O4B (a crown ether).
- DESIGN §0 defines drug-like (primary) as ChEMBL max phase ≥ 1 **or pChEMBL ≥ 6**, **and not a cofactor/metabolite**. Ro5+QED is a *secondary* definition.
  The pipeline implements only max_phase and has no metabolite exclusion. It then merges the Ro5+QED rule into `ligand_class` as if it were primary.
- The pChEMBL arm is missing, so real investigational inhibitors fall into `other_organic`. In the 20-ligand spot check: K2X (a peptidomimetic protease inhibitor), TVL and FRJ.
- H1 compares cofactor vs drug-like, so the current comparison group is mostly buffers, cholesterol and Ro5-passing fragments.

**Fix.**
- Make `drug_like` = ChEMBL route only (phase ≥ 1 OR pChEMBL ≥ 6 via the ChEMBL activities table, which is already fetchable through `scripts/fetch_chembl*`).
- Keep Ro5+QED as a separate secondary column, not as a class.
- Exclude the BioLiP/PLINDER buffer and solvent lists, CCDs with a cofactor class, and metabolites. A ChEBI/HMDB mapping or at minimum a curated list (CLR, CIT, CHD, IHP, U10, FOL, sterols, bile acids) is needed.
- Report the declared FALLBACK in the results tables as DESIGN §0 requires.

### H-2. Symmetry-group collapse fails for about 43% of symmetric copies in homodimers
**Evidence (t4, t4b).**
- `sym_group` is an exact hash of the CCD plus the multiset of (UniProt, position) over contacted residues. Any one-residue difference between copies breaks the collapse.
- Homodimers (2 protein chains) have 39,252 (pdb, CCD) groups with exactly 2 primary instances. Only 57.3% collapse to one sym_group.
- Of the 16,746 pairs that do not collapse:
  - 95% have UniProt-residue Jaccard ≥ 0.5;
  - 76% have Jaccard ≥ 0.8.
  So they are clearly the same site.
- For the off-axis MBS pairs that matter for H3, 9,071 pairs are found and only 54.3% are collapsed.
- Known cases:
  - 2TSC UMP: A:11;B:2 vs A:2;B:10 gives 2 sym_groups.
  - 1HVY raltitrexed (D16): 11 vs 12 residues gives 2 sym_groups.
  - 3GRS FAD and 1HVY UMP collapse correctly.
- Across all homomers, 290,643 instances map to 143,845 sym_groups.

**Impact.** H1 is unaffected because it uses one instance per cluster30 × CCD. H2, H3 and any "site"-level count double-count roughly half of symmetric sites.
The `mbs_by_assembly_class.tsv` fractions are instance-level. That is acceptable for a descriptive table, but they must not be labelled "sites, symmetric collapsed".

**Fix.** Collapse by tolerance, not by hash. Within an entry, link same-CCD instances whose UniProt-position sets (chain-agnostic) have Jaccard ≥ 0.5, or whose centroids map onto each other under the chain-pair superposition (< 2 Å). Take connected components as the site.

### H-3. The FALLBACK additive label (entries newer than BioLiP2, 79,788 primary instances = 14%) differs from the BioLiP2 label: it deletes all Fe–S clusters and keeps ions
**Evidence (t6).** Exclusion rate under the FALLBACK label (not_covered) vs under BioLiP2 (covered):

| Class | FALLBACK (not_covered) | BioLiP2 (covered) | Notes |
|---|---|---|---|
| metal_cluster | **100%** | 6.9% | SF4, FES, F3S: all excluded. PLINDER's "< 2 carbons" filter catches them, and the single-atom ion exemption does not cover clusters. |
| metal_ion | 3.1% | 34.3% | |
| CL | 67.8% kept | 1.2% kept | |
| NA | 84.6% kept | 3.0% kept | |
| inorganic | 70% | 94% | |
| drug_like | 7% | 19% | |

- Not-covered entries are 21% of heteromer primary instances, against 12% for homomers and 9% for monomers. EM-rich newer entries are therefore labelled by a different rule.
- Homomer MBS fraction: 26.9% when covered vs 32.6% when not covered.
- Of the CL in the primary set, 4,173 of 4,870 come from FALLBACK entries.

**Fix.**
- Exempt all-inorganic metal clusters (metal_cluster class) from PLINDER's carbon/heavy-atom filters.
- Apply the BioLiP artifact list to ions (CL, NA, K, IOD, BR, …) in FALLBACK regardless of copy number, or, better, train the FALLBACK on the covered set and show agreement per class.
- Add `additive_label_source` as a covariate or stratum in H1.

---

## MEDIUM severity

### M-1. Heteromer cofactor MBS is driven by chlorophylls (photosystems / LHC)
**Evidence (t4).**
- Heteromer `marsh_class == cofactor` has n = 73,063, and 75.7% of those are chlorophyll/pheophytin codes. CLA alone is 45,985 of the 61,087 heteromer heme_porphyrin instances.
- Heteromer cofactor MBS is 40.7% overall but **30.5% without chlorophylls**. Homomer cofactor is barely affected: chlorophylls are 4.8% of it, and MBS goes from 38.3% to 38.7% without them.
- The top 10 entries by primary instances are all CLA-dominated photosynthetic complexes: 1RWT (840), 9JJ8, 9UEN, 7UMH, 6KIF, 6KIG, 6NWA, 6K33, 6LY5, 9JZR.
- The top 1% of heteromer entries hold 29% of heteromer instances. For homomers the figure is 14%.

**Fix.** Put chlorophylls/bacteriochlorophylls/pheophytins in their own class rather than heme_porphyrin → Marsh "cofactor". In every instance-level table, report cluster-deduplicated fractions next to the raw ones. H1-rep instance fractions for heteromer cofactors are already 45.7% on n = 1,174, so dedup changes the picture. Also apply the pre-registered "leave out top-10 clusters" sensitivity to the descriptive table.

### M-2. Same-UniProt "homomers" that are really processed or heterodimeric chains
**Evidence (t4, t4c).**
- 300 random homomer entries were checked (seed 0). In 8 (2.7%), the protein entities differ in length by more than 10% while sharing one UniProt. Extrapolated, that is about 2,300 of 85,131 homomer entries (rough 95% CI 1.2–5.2%).
- Examples: thrombin light/heavy chains (11RI, 3QTV, 2BXU), FVIIa (5PAI), RSV F1/F2 (4GIP), GGT1 heavy/light (7LBC), tankyrase-2 fragments (4UI3).
- 7 more entries (2.3%) have near-identical multiple entities, i.e. mutants. Those are acceptable as homomers.
- Direct effect on MBS: 1,526 homomer MBS instances in 815 entries have contributing chains from different entities, which is 1.9% of homomer MBS. The top proteins are:
  - P05164 myeloperoxidase (the heme spans the heavy/light chains of one cleaved protomer, so these MBS are intra-protomer artefacts);
  - Q32ZE1;
  - P03366 HIV-1 RT (p66/p51, 585 instances, 39.5% MBS; these are real heterodimer sites mislabelled as homomer).
- 300 random heteromers had no misclassified entries.

**Fix.**
- Define the identity key as (UniProt, SIFTS UniProt start–end segment).
- If two chains map to non-overlapping UniProt segments, treat them as one protomer when they come from proteolytic processing (light/heavy chain). Otherwise call the assembly a heteromer (p66/p51 overlaps but differs in length by more than 10%, so it becomes a heteromer).
- At minimum, flag `homomer_multi_entity_lendiff` and exclude it in a sensitivity.

### M-3. The Jaccard on-axis flag misses textbook on-axis ligands that are disordered across the 2-fold
**Evidence (t5).**

| Case | Result |
|---|---|
| 1HVR XK2 | on-axis (J = 0.75, axis distance 0.05 Å). Correct. |
| 3TCT tafamidis | both T4 channels on-axis (J = 0.57 and 1.00; axis distances 0.38 and 0.41 Å). Correct. |
| 2ROX T4 | channel 1 is on-axis (J = 0.62). **Channel 2 is called off-axis (J = 0.40) although the centroid is 0.47 Å from the C2 axis.** |
| 1BM7 TTR + flufenamic acid | both channels called off-axis (J = 0.20 and 0.25; axis distances 1.02 and 0.88 Å). The two orientations of channel 1 are kept as 2 primary instances: their centroids are about 2.04 Å apart, just above the 2.0 Å dedup cutoff. They share a sym_group, so they collapse at site level. |
| 3GRS FAD, 1HVY/2TSC UMP | off-axis (J = 0, axis distance about 15–26 Å), as expected. |

Disordered-on-axis ligands are modelled in one orientation, so each chain sees a different half of the ligand and the Jaccard is low.

**Fix.** Call a site on-axis when it meets either the Jaccard rule OR the geometric rule (axis distance ≤ 3 Å and rotation 180 ± 15°). Report the agreement between the two rules. In dedup, add an atom-overlap criterion: copies whose heavy atoms are more than 50% within 1.5 Å count as one ligand.

### M-4. Binding-residue threshold: DESIGN §0 vs implementation
DESIGN §0 says "A residue binds if it has ≥2 atom contacts (BioLiP2; stored as a flag, primary)". The pipeline's primary column uses ≥ 1 atom pair, and ≥ 2 is stored as `*_biolip2`.
The headline numbers differ: homomer 27.6% vs 24.0%, heteromer 33.3% vs 28.5%. The design text is ambiguous.
**Fix.** The PI should decide which rule is primary and log the decision in RESEARCH_TREE.md before H1–H3 are run.

### M-5. Additive label: BioLiP2 curation gaps remove real ligands
- In the 20 excluded-as-additive ligands (seed 0), 18 are correct (EDO ×7, GOL ×3, SO4 ×2, DMS, PGE, MLI, IOD, HG, a detergent-like lipid 1EM).
- 1 is a clear error: LUT (lutein) in 7EW6, a light-harvesting complex, where it is a native pigment.
- 1 is doubtful: CS in 7BMX, which is flagged as wwPDB ligand-of-interest.
- Systematic cases:
  - FOL is 93.6% excluded;
  - CB3, a classic TS antifolate inhibitor, is excluded in 2TSC because BioLiP2 lists only UMP;
  - SF4 is excluded 14.8% overall (from H-3);
  - 45,227 excluded instances are wwPDB "SUBJECT OF INVESTIGATION".
- **Fix.** DESIGN §0 already says additive = NOT (BioLiP2-relevant OR ligand-of-interest). The implementation ignores the LOI arm: `is_additive` is BioLiP-only. Keep any ligand-of-interest instance unless it is on the BioLiP/PLINDER solvent list. This would rescue the FOL, CB3 and LUT type cases.

### M-6. 20-ligand primary-set class spot check
- 2 of 20 are clear errors:
  - GOL in 2EVU is kept as BioLiP-relevant (1,038 GOL, 2,237 SO4, 3,511 PO4 and 4,870 CL are kept in the primary set);
  - CLR cholesterol is labelled drug_like.
- Borderline: LNR norepinephrine is drug_like (a metabolite, see H-1); 78M monoacylglycerol and CYC phycocyanobilin are other_organic.
- The rest look chemically sensible (MN, FE, ZN metals; ATP, ADP, GTP nucleotides; TPP cofactor; FBP sugar; PEK lipid).
- **Fix.** Apply the BioLiP solvent list as a hard filter on top of BioLiP-relevant, at least as a sensitivity.

---

## LOW severity / documentation

- **L-1.** DESIGN §0 gives two dedup cutoffs ("centroid < 1 Å" under Ligand instance, "< 2 Å" under Symmetric-copy dedup). The code uses 2.0 Å. Make the text consistent.
- **L-2.** DESIGN §0 says metals are "never pooled with organic ligands", but the `ALL` rows of `mbs_by_assembly_class.tsv` and the 27.6%/33.3% headline pool metals. Add an organic-only ALL row: computed from the marsh_class rows (n-weighted cofactor + small molecule): homomer 0.361, heteromer 0.418.
- **L-3.** Antibody chains: DESIGN says they are excluded from primary analyses. The primary `is_mbs` still includes them (437 instances are MBS only via an antibody chain; 4,114 primary instances are in antibody entries). This is small, but the primary column should be `is_mbs_excl_antibody`, or the doc should be changed.
- **L-4.** 4,773 of 85,208 H1 groups (cluster30 × CCD) span more than one assembly type. The representative may be a monomer entry of an oligomeric cluster, while H1 is restricted to oligomers. Choose representatives after restricting to oligomers.
- **L-5.** 52 primary instances have a null cluster30 and are silently absent from H1.

---

## Checks that passed

1. **Headline fractions (t1).** I rebuilt the primary set from the raw flags (not dup ∧ bound ∧ not glycan ∧ not additive). It is identical: 561,250 rows, XOR = 0. I recomputed MBS from `chain_counts`: 0 mismatches against `is_mbs`, `n_chains_contacted` and `is_mbs_min2`. All 9 homomer/heteromer/monomer × Marsh-class fractions and n in `mbs_by_assembly_class.tsv` match to 4 decimals (e.g. homomer cofactor 0.3828 / 45,797; heteromer small molecule 0.4359 / 44,357). Before the additive filter: homomer 0.2773, heteromer 0.3717. Both match.
2. **Independent contact recomputation (t2, `task2_contact_recompute.tsv`).** Method: my own gemmi assembly generators plus a numpy all-pairs vdW + 0.5 Å calculation, eligible chains = polypeptide entity SEQRES ≥ 30, first altloc, no H, model 1, on 40 stratified instances (seed 0). Agreement with the pipeline is 40/40 for n_chains_contacted, MBS and the total binding-residue count. Matched centroids differ by ≤ 0.01 Å.
   - Caveat: both implementations use gemmi `Element.vdw_r`, so an error in the radius table would not be caught.
   - Sample observations: 4 of the 10 heteromer-MBS draws are CLA (M-1). One homomer MBS draw is SO4 (1ODL, BioLiP-relevant). One monomer draw is BU3 (butanediol), which is kept.
3. **H1 representative (t5).** Recomputed as best resolution, then earliest ID, per cluster30 × CCD in the primary set. It agrees in 85,208 of 85,208 groups, with exactly one rep instance per group.
4. **C2 homodimer cofactor (t5).** 3GRS FAD: both copies are MBS (28 + 2 residues), off-axis (axis distance 26 Å), and share one sym_group (collapsed, as designed). 1HVY UMP behaves the same way.
5. **Heteromer typing (t4c).** 0 of 300 random heteromers have identical entities.
