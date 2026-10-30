# Review r1 — prior work & definitions (DESIGN v1)

Prior-work reviewer; web checks 2026-10-29; no analyses run.

## 1. Definitions vs Marsh 2018 and current practice
Sound: biounit MBS, separate metals, peptides/nucleic acids excluded. Problems:
- **Additives (major).** The "simplified" BioLiP rule leaves out BioLiP2's PubMed check. Without it, interface-prone GOL/SO4/PEG pass. Marsh used BioLiP's labels, not a re-derivation. Join the precomputed BioLiP2/Q-BioLiP relevance flags (keyed by PDB, ligand chain and residue number) as the primary filter. Keep the geometric-only rule and the PLINDER filters as sensitivities.
- **Contacts.** BioLiP2 counts a binding residue only with ≥2 atom contacts (vdW + 0.5 Å), and a site needs ≥2 residues. DESIGN leaves out the ≥2-atom condition.
- **Symmetry copies (homomers).** In a homomer, a ligand on a crystallographic special position, or half-occupancy/altloc copies, becomes overlapping duplicates once assembly operators are applied. Deduplicate same-CCD copies whose centroids are within 2 Å. Then collapse symmetry-equivalent sites, giving each site weight 1/n copies or one per entity+CCD. Otherwise symmetric homomers are over-counted.
- **Assembly choice and crystal contacts.** Assembly 1 is often author-defined. Use PDBe preferred assembly (PISA, as LIGYSIS does) as a sensitivity, and require an interface of ≥~1000 Å² or ≥40 residues between contributing chains (Abrusán & Foguet 2023).
- **Glycans.** Since the 2020 remediation, oligosaccharides are *branched* entities, not non-polymers. N/O-glycans (covalent link to Asn/Ser/Thr in `struct_conn`) are PTMs. They inflate heteromer MBS counts. Exclude them from the primary set and report them separately.
- **Covalent and multi-CCD ligands.** Flag covalent ligands (`struct_conn` covale) but keep them, because covalent drugs are drug-like. Merge covalently linked or fragmented CCDs into one ligand, as PLINDER does.
- **Peptide and antibody chains.** For a chain to make a site MBS, require a protein chain of ≥30 aa. Flag or exclude hapten sites between VH and VL (Ig domains), which are a large artefactual "MBS heteromer" class.
- **Ligand classes.** Map CCD codes to ChEMBL via UniChem. Also emit Marsh's 3-class scheme (cofactor / metal / small molecule) so results are comparable.
- **Unit.** Marsh counted *proteins* (5 groups); Q-BioLiP/PLINDER count ligands. Report both.

## 2. Novelty positioning
- **Part 1 is a replication.** The overall MBS share is already published: Q-BioLiP 28.8%, PLINDER 25%. Abrusán & Foguet 2023 covers homodimer MBS. What is new is the class-stratified and druggability breakdown. Drop "temporal comparison to Marsh": Marsh gives figures only, with no MBS fractions to compare against. Replace it with the direct comparisons to Q-BioLiP and PLINDER.
- **Dimer-only pockets on AFDB homodimers have prior art.** ProtPocket (GitHub, 2026; not peer-reviewed) diffs monomer and dimer pocket runs on the March 2026 AFDB homodimer release. Cite it and soften LITERATURE (h2) "no atlas". Our differentiation is: validation against PDB MBS sites, heterodimers, cross-species comparison, and ligand-class transfer.
- **Direct comparisons:** Q-BioLiP 28.8%, PLINDER 25%; LIGYSIS recall (fpocket+PRANK 60.4%, P2Rank ~52–54%); Foguet 85–95%; Chen & Zacharias >75%.

## 3. Public datasets to use directly
- **PLINDER** `2024-06/v2` `index/annotation_table.parquet` plus `splits/split.parquet` (GCS, about GB-scale). Do **not** fetch `systems/` (hundreds of GB). It gives per-system interacting chain counts, artifact/ion/cofactor/covalent flags, biounit IDs and UniProt. Cross-check census; Q2.1 source.
- **Q-BioLiP / BioLiP2** annotation text files (yanglab.qd.sdu.edu.cn/Q-BioLiP/Download). Skip the 73–85 GB structure archives.
- **LIGYSIS** human benchmark (2,775 chains / 6,882 sites; Zenodo 14645504). The full DB is web-only. Use it to calibrate burrow with the published metric definitions.
- **Runs N' Poses** (Zenodo 14794786) for Q2.5, which records `num_protein_chains`.
- *Could:* Moine-Franel PPI pockets (Zenodo 10805580); wwPDB ligand-of-interest flag.

## Changes
**Must:** §1 additives, contacts, symmetry, glycans/covalent, chain/Ig rules, Marsh layer; §2 ProtPocket + drop temporal-Marsh. **Should:** PISA/interface sensitivity, PLINDER index, LIGYSIS calibration, ChEMBL drug-likeness, druggability conditioned on class. **Could:** §3 validation subsets, PINDER.
