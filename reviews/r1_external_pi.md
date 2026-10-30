# Review r1: external PI (biological insight and priorities)
Reviewer: external PI agent, 2026-10-29 16:40, on DESIGN.md v1. No analyses run.

## Overall
The definitions are solid and Part 1 follows Marsh faithfully. The weakness is prioritisation. The design is a list of descriptive censuses (fractions, descriptor distributions, kNN transfer) with no biological hypotheses that can fail. As written, it will deliver "interface pockets are bigger and less confident" plus a benchmark number. That is not the publishable story. The novel material is MBS druggability measured for the first time, dimer-only pockets in predicted human complexes, and what those pockets mean (symmetry, disease, drugs).

## Proposed changes
**MUST**
1. **Feasibility gate for Q2.1 before anything else in Part 2.** Count the PDB MBS sites (additive-filtered, nonredundant) whose UniProt pair is among the 10.6k human dimer models, homo and hetero separately. If there are fewer than ~100 hetero MBS sites, report hetero recovery as anecdotal and do not build the analysis around it.
2. **Within-protein paired druggability comparison.** Compare MBS and SBS sites within the same protein or 30% cluster (paired or conditional test). Pooled comparisons are confounded by protein family, because cofactor-heavy enzymes dominate MBS. Primary endpoint: drug-like ligand occupancy. Secondary: pocket descriptors. Fix fpocket or burrow `--ranker small` as the single pre-registered descriptor score now.
3. **Homodimer symmetry classification.** Split homomer MBS into *on-axis* (one ligand spans the C2 axis, one site per dimer) and *paired off-axis* (two symmetric sites, each formed by both chains). On-axis sites are classic allosteric or drug sites (TTR/tafamidis, PKM2 activators, HIV protease). The hypothesis that on-axis sites are enriched for drug-like and allosteric ligands and paired sites for cofactors is novel and testable, and it gives a natural deep dive.
4. **Disease-variant burden at dimer-only pockets (human).** Map AlphaMissense scores (bulk download, public) and ClinVar pathogenic variants onto lining residues of dimer-only pockets, single-chain pockets, interface residues outside any pocket (essential control), and surface residues. This tests Marsh's ~30% heteromer claim at proteome scale. It is cheap and publishable.
5. **Cut** the temporal comparison to Marsh (1.5) and **demote** cross-species Q2.4 to *could*. Ortholog mapping across 13 species is a time sink; at most, compare human with mouse.

**SHOULD**
6. **Replace kNN ligand-class transfer (Q2.3) with a simpler readout**: the burrow small-molecule score plus nearest PDB-ligand-site similarity, reported only as supporting evidence. TOOLS.md says pocketcompare is weak on PPI-interface pockets, so a cross-validated AUC there risks a null result that consumes hours.
7. **Druggable-opportunity list.** List human proteins whose only confident ligandable pocket is dimer-only (no qualifying pocket in the AFDB monomer atlas), cross-referenced to ChEMBL target status and Open Targets disease links. This seeds deep dives (candidate glue or stabiliser sites).
8. **Pocket-presence null for Q2.2.** Any two docked surfaces create grooves. Compare dimer-only pocket frequency and score against interface size-matched expectations, and stratify by ipSAE, so that "dimer pockets" are not just docking artefacts.
9. **Pin co-folding to one hypothesis.** For about 40–60 post-cutoff MBS cases, compare ligand placement with the partner present versus absent (dimer versus monomer input). Do the timing test now. If AF3 on Turing takes more than 10 min per job, use Boltz-2.

**COULD**
10. **Paralog heteromers.** Test whether heterodimers of paralogs (mmseqs identity ≥30%) retain MBS more than non-paralog heteromers. This addresses homomer-to-heteromer evolution.
11. **Allostery annotation.** Overlap MBS sites with ASD/AlloSteric entries as an independent test of Abrusán & Marsh 2019.

## Timeline suggestion
Part 1 (with change 2) by 23:00. Q2.1 gate, Q2.2 and the dimer-only table by 04:00. Changes 3 and 4 by 07:00. Co-folding runs overnight from 01:00 on GPUs 1–3. Pick deep dives at 08:00 from on-axis drug sites, high-burden dimer-only pockets and the opportunity list.
