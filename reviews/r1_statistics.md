# Review r1 — statistics (DESIGN draft v1)

Biostatistics agent; no analyses run.

## MUST
1. **Recovery baseline and pocket size (Q2.1).** Best-pocket recall, maximised over all pockets, rises with pocket size and count; dimer runs give more, larger (merged) pockets, which favours dimers. The cross-protein "random pocket mapped by position" null is meaningless, since residue numbers do not carry across proteins. Replace it with a **within-structure, size-matched null**: contiguous surface patches with the same residue count as each predicted pocket, drawn from the same model (seeded). Primary metric: **DCC ≤ 4 Å on top-(n+2) pockets by score**, which does not depend on size. Report residue F1/Jaccard next to recall, with pocket size, at equal top-k for dimer vs isolated chain.
2. **Novelty leakage (Q2.2).** The atlas contains each protein's own monomer, so single-chain pockets find themselves and interface pockets cannot. Remove the query's own accession and its ≥30 %-identity paralogs before computing max-cosine. pocketcompare is weak on PPI pockets, so add a positive control (PDB interface sites in homologous complexes); otherwise "novel" may just be the method's blind spot.
3. **kNN class-transfer leakage (Q2.3).** A 30 % sequence split leaks remote homologues. Split by Pfam clan / ECOD X-group (or Foldseek cluster) instead. Trained on holo PDB pockets, applied to apo predicted pockets (many bind nothing): add an abstain/"no ligand" option (distance threshold) and report AUPRC against prevalence as well as AUC. Test out-of-domain on predicted pockets that Q2.1 maps to PDB sites, separately for MBS and SBS.
4. **Symmetric copies and units.** Homomer sites repeat per protomer/interface copy. Collapse copies to one per (entity set, CCD, site) before counting, then weight by 30 % cluster. Use the same collapse for predicted homodimer pockets.
5. **Regression specification (Part 1.3).** Monomers cannot have MBS (structural zero), and n_chains is close to part of the outcome definition. Fit on oligomers only, with assembly type, log(n_chains), method, structural-genomics keyword, release year and ligand heavy atoms (spline). Wrong assemblies and crystal-packing ligands create false MBS: add a sensitivity analysis on entries where author and software assemblies agree.
6. **Training leakage and power (Q2.1).** AFDB complex models (and burrow, which was trained on BioLiP) have seen most benchmark sites. A post-cutoff stratum only counts if the UniProt pair had **no complex entry before the cutoff**. Count the expected n per stratum (homo/hetero × MBS/SBS × post-cutoff) **before** running. Pre-specify ≥30 clusters per stratum for inference, else descriptive only. For scale, the burrow benchmark had n = 66 hetero-interface sites.

## SHOULD
7. **Interface-pocket rate scaling.** Large pockets and large interfaces make touching both chains more likely. Compare the observed interface-pocket fraction with the interface share of SASA; ≥3-residue rule as primary. **Dimer-only**: use an overlap coefficient instead of Jaccard ≥0.3, because merging and splitting deflate Jaccard. Sensitivity at 0.2/0.5.
8. **Selection bias.** The confident-dimer set (ipSAE/pDockQ2 filter) is enriched for pairs that have PDB complexes. Inference is conditional on that set; stratify recovery by PDB-complex availability and ipSAE tertile.
9. **Over-represented items.** Sensitivity leaving out the top-10 clusters and the top ligands (by computed count). Condition MBS-vs-SBS druggability on ligand heavy atoms (MBS ligands are larger).
10. **Multiple testing.** Pre-register 4–6 primary hypotheses with BH per family; label the rest exploratory. One bootstrap seed is primary; the others only show stability.
11. **Monomer control.** The atlas (AFDB v6 monomers) is not like-for-like with the v1 complex models. Primary control: isolated chains from the same dimer coordinates, with identical settings and filters; atlas secondary.

## COULD
12. Q2.4: condition on an ortholog dimer model existing, since confidence filters select differently per species.
13. Fix TOOLS.md §6.6 (DBREF): use candidates.csv.
