# Literature & dataset review: pockets in protein dimers vs monomers

Compiled 2026-10-29. Covers single-chain vs multi-chain (interface) ligand-binding sites in PDB structures and in predicted homo/heterodimers, druggability, ligand classes, and co-folding re-evaluation.

**Verification legend.** Every citation below was checked against Crossref (DOI, authors, journal, year).
- **[full text]**: I read the paper's full text (local PDF or PMC full text).
- **[abstract/summary]**: the numbers come from the abstract, a search snippet or a fetched-page summary. Check them against the paper before quoting.
- **[UNVERIFIED]**: something I could not confirm.

Local PDFs read: `~/reading/protein_ligand_docking/2018_marsh_abrusan_ligand-binding-evolution_cellreports.pdf` (full) and `2025_kryshtafovych_gilson_casp16-docking_proteins.pdf` (skimmed: abstract, targets, baselines, discussion).

---

## (a) The Abrusán & Marsh framework: definitions to adopt

**Primary paper.** Abrusán G, Marsh JA (2018). *Ligand binding site structure influences the evolution of protein complex function and topology.* Cell Reports 22:3265–3276. doi:10.1016/j.celrep.2018.02.085. [full text]

### Definitions (verbatim intent)
- **SBS (single-chain binding site):** "sites where the binding residues are located on only one chain of the complex".
- **MBS (multi-chain binding site):** "sites where the binding residues of a particular ligand are located on more than one chain". Fig. 1 gives two examples:
  - SBS: ketohexokinase dimer, PDB 3Q92. Each subunit has its own pocket.
  - MBS: anemone STING dimer with c-di-GMP, PDB 5CFL. The dimer has one site formed by both chains.
- **Five protein groups:** monomers, SBS homomers, MBS homomers, SBS heteromers and MBS heteromers.
- **Quaternary-structure (QS) label is protein-level and set from biological units:**
  - A protein is a *heteromer* if it appears in ≥1 heteromeric PDB biounit. A heteromer biounit is any entry with chains from ≥2 different proteins, whatever the chain length.
  - Otherwise it is a *homomer* if it appears in ≥1 homomeric biounit.
  - Otherwise it is a *monomer*.
  - Caveat: the paper does not say exactly how a protein with both SBS and MBS ligands is assigned to a group. Treat this as ambiguous [UNVERIFIED; may be in the Datashare supplement, https://datashare.is.ed.ac.uk/handle/10283/3024].
- **Ligand contacts / binding residues:**
  - Binding residues and the SBS/MBS call come directly from **BioLiP's binding-residue annotation**, which the paper does not re-derive.
  - BioLiP's criterion: a residue is a binding residue if any of its atoms is within the sum of van der Waals radii + 0.5 Å of a ligand atom (Yang et al. 2013). BioLiP2 adds a requirement of ≥2 such atomic contacts per residue (see (f)).
  - Separately, for the ProBiS pocket searches, query binding sites were "surface residues within **3 Å** of the ligand".
- **Ligand classes:**
  - *cofactors*: ligands annotated as cofactors in ChEBI, or listed in Fischer et al. 2010.
  - *metal ions*.
  - *small molecules*: everything else that is not a cofactor, metal, nucleic acid or peptide.
  - Nucleic acids and peptides were excluded.
- **Additive / crystallization-artifact exclusion:**
  - No list of their own. They relied entirely on **BioLiP's biological-relevance filtering**, which in 2013 meant a semi-manual artifact list plus literature checks.
  - They dropped all **13,413 entries where the biounit has more protein chains than the asymmetric unit**, because BioLiP's annotation is asymmetric-unit based and is incomplete for those entries.
  - That left **71,925 BioLiP structures** before this filter, after removing structures that bind only nucleic acids or peptides.
- **Other entry filters:**
  - Removed virus entries, fibrils, helical assemblies, and entries with sequences from >1 species (this also removes most antibody–peptide entries).
  - Removed chimeric sequences.
  - Clustered sequences at 100% identity.
  - Excluded structures with >99,999 atoms (PDB-format limit).
- **Homology:**
  - All-vs-all BLAST, e-value < 1e-4, then MUSCLE realignment to get global identity.
  - ProBiS local structural site search. A hit counts if Z ≥ 2 and e-value < 1e-4.
  - The "fraction of different binding sites" is the fraction of a protein's ligands with no significant ProBiS hit in its homolog, plotted against sequence identity.
  - RMSD comparisons used homolog pairs at 30–49% identity, with pocket size (number of aligned ProBiS vertices) as a covariate.
- **Chemistry:** ChEBI "is_a" ontology distance, defined as the fraction of non-shared terms. ChEBI role-ontology enrichment used GeneMerge with BH-FDR. Of ~20,000 BioLiP ligands, only ~3,100 had ChEBI entries.

### Key findings, with the numbers the paper actually gives
1. **MBS homomers diverge in function more slowly.** Their binding sites are much more often found in homologs than SBS or monomer sites are, especially in distant homologs (Fig. 2A–B). **MBS heteromers do not show this pattern** (Fig. 2C–D).
2. **The effect differs by ligand class:**
   - Small-molecule pockets are much more conserved in MBS homomers.
   - Cofactor pockets are conserved in every group, with almost no difference above 40% identity.
   - Metal pockets are the least conserved, with no group differences (possibly because small pockets are harder to detect).
3. **Ligand composition and count:** MBS complexes have **more cofactors and fewer metals** (Fig. 1C) and slightly fewer ligands per protein (Fig. 1D). The paper reports these only as figures; there are no exact percentages in the text.
4. **QS conservation:**
   - In cofactor- or metal-binding homomers, MBS complexes keep their subunit number much better than SBS complexes.
   - In heteromers the direction is reversed: MBS heteromers change subunit number faster.
   - Small-molecule-only complexes show no QS effect.
5. **Ligand chemistry and pocket RMSD:**
   - Ligands of MBS homomers are chemically more similar across homologs (ChEBI).
   - Their pockets have lower RMSD between homologs (Fig. 4A; Wilcoxon, BH-corrected).
   - Apo vs holo pocket RMSD shows **no** clear MBS vs SBS difference (Fig. 4D).
6. **Enrichments:**
   - Nucleobase- and ribonucleotide-containing ligands, and the ChEBI role "metabolite" (fundamental metabolites), whether or not cofactors are included.
   - GO terms: acyl-CoA dehydrogenase, transaminase, thiamine-pyrophosphate binding, cofactor binding, transmembrane transport.
7. **Flexibility and taxonomy:**
   - MBS subunits are more flexible (A_rel; Fig. 6A), and this is not explained by symmetry type.
   - MBS frequency is higher in prokaryotes than in eukaryotes, especially for heteromers (Fig. 6D; the y-axis runs 0–20%, and the text gives no exact values).
8. **Pathogenic variants:**
   - Baseline pathogenic fraction rises from monomers to homomers to heteromers, and is higher in MBS than SBS.
   - Binding sites of **small-molecule MBS heteromers reach ~30% pathogenic** (Ensembl disease variants vs ExAC).
   - Metal-binding MBS homomers have only 5 pathogenic mutations, so power is low.
9. **Drug-design claim:**
   - SBS complexes are better targets for *selective* drugs.
   - MBS homomers suit broad-spectrum antibiotics and multitarget drugs or repurposing (example: nelfinavir and HIV protease).
   - **This is a hypothesis; they did no druggability analysis.**

### Follow-ups by the same group (all Crossref-verified)
- **Abrusán & Marsh 2019**, *Mol Biol Evol* 36:1711–1727, doi:10.1093/molbev/msz093. Allostery. [abstract/summary]
  - Allostery is significantly more frequent in MBS homomers and in metal-binding MBS heteromers than in SBS complexes.
  - In MBS homomers, signal-transduction pathways cross the interface.
  - Dataset: 403 complexes (158 MBS homomers, 135 SBS homomers, 44 MBS heteromers, 66 SBS heteromers).
  - Uses the same 3 Å ProBiS site definition.
- **Abrusán & Marsh 2019**, *J Mol Biol*, doi:10.1016/j.jmb.2019.07.014. [abstract/summary]
  - MBS homomers depend more on folding and assembly quality control and are enriched for chaperone and chaperonin clients.
- **Abrusán & Foguet 2023**, *Mol Biol Evol* 40:msad070, doi:10.1093/molbev/msad070. [abstract/summary] This is the most directly relevant follow-up.
  - Dataset: 992 PDB homodimers clustered at 30% identity, plus 148 human homodimers.
  - Binding residues defined by AREAIMOL solvent-accessibility change.
  - **Metals excluded.**
  - Interfaces with <40 residues (~1000 Å²) excluded as likely crystal artifacts.
  - Findings: **85–95% of MBS homodimers are predicted functional as dimers**, while **35–42% of SBS homodimers with small-to-medium interfaces likely have no functional QS**.
  - SBS sites are much more common than MBS among homodimers. The exact proportion is [UNVERIFIED]; check the paper.

### Recommended operational definitions for this project
These are our choices, not Marsh's.
1. **Assembly:** use the biological assembly, not the asymmetric unit.
   - PDBe preferred assembly (PISA), as LIGYSIS does; or biounit 1 when cross-checking PLINDER.
   - Record the assembly ID.
2. **Contacts:** heavy-atom distance ≤ vdW radii sum + 0.5 Å, with ≥2 atomic contacts per residue (BioLiP2).
   - Report sensitivity at a flat 4.0 Å and 4.5 Å cut-off.
   - Contacts are counted per chain instance.
3. **MBS call:** binding residues from ≥2 distinct protein chain instances.
   - Add a robustness filter: the minor chain contributes ≥2 binding residues, or ≥10% of ligand contacts.
   - Also report the strict "any residue" version (Marsh/BioLiP) for comparability.
   - The thresholds are ours and should be tuned on data.
4. **MBS-homo vs MBS-hetero:** decide by UniProt identity of the contributing chains (SIFTS).
5. **Ligand classes:**
   - Cofactor: ChEBI "cofactor" role plus a CoFactor/Fischer-style list.
   - Metal.
   - Small molecule split into drug-like vs metabolite, using ChEBI "metabolite" role and/or ChEMBL/DrugBank mapping via UniChem.
   - Exclude peptides and nucleic acids.

---

## (b) Prior quantitative findings on multi-chain sites

| Source | Data | Key number(s) | Status |
|---|---|---|---|
| Abrusán & Marsh 2018, Cell Rep | BioLiP (2013-era) PDB, biounits | Qualitative and figure-based: MBS homomers have more conserved pockets and ligands, more cofactors, and fewer metals; ~30% pathogenic in small-molecule MBS-heteromer sites | [full text] |
| Wei, Wang, Peng, Yang 2024, *Genomics Proteomics Bioinformatics* 22:qzae001, doi:10.1093/gpbjnl/qzae001 (Q-BioLiP) | Quaternary-structure-based BioLiP, ~210k PDB entries, ~2.7M entries (Oct 2023) | **28.8% of biologically relevant ligands on oligomeric targets contact >1 chain** (dimers up to 10-mers) | [abstract/summary] |
| Durairaj et al. 2024, bioRxiv doi:10.1101/2024.07.17.603955 (PLINDER) | 449,383 PLI systems, biounits | **25% of holo systems have >1 interacting protein chain** | [abstract/summary] |
| Abrusán & Foguet 2023, MBE | 992 homodimers | 85–95% of MBS homodimers functional vs 35–42% of SBS small-interface homodimers possibly non-functional | [abstract/summary] |
| Abrusán & Marsh 2019, MBE | 403 complexes | Allostery significantly enriched in MBS homomers | [abstract/summary] |
| Moine-Franel et al. 2024, *Sci Data* doi:10.1038/s41597-024-03233-z | PPI-related pockets (VolSite): 23,238 pockets, ~3,700 proteins, ~3,500 ligands | Classes are orthosteric-competitive (ligand within 1 Å of partner), orthosteric non-competitive, and allosteric. These are heterodimer PPI-modulator pockets, not MBS. | [abstract/summary] |

**Takeaway.** Different definitions place the MBS share among ligands in oligomeric PDB structures at roughly **25–30%**. Nobody has reported this stratified by ligand class (cofactor vs drug-like) in a current PDB snapshot with a modern artifact filter. That is an easy first deliverable.

---

## (c) Druggability of interface pockets

- **Jubb, Blundell, Ascher 2015**, *Prog Biophys Mol Biol* 119:2–9, doi:10.1016/j.pbiomolbio.2015.01.009. [abstract/summary]
  - Analysed 15,500 non-redundant PPI interfaces.
  - Small single-residue sub-pockets at globular interfaces are frequent and potentially selective.
- **Chen & Zacharias 2023**, *ACS Cent Sci* 9:969–979, doi:10.1021/acscentsci.3c00003. [abstract/summary]
  - In **>75% of 226 protein–protein complexes, interface cavities were suitable for drug-like compounds**.
  - Good PPI stabilizers use a "dual-binding" mechanism, binding both partners about equally. This is the molecular-glue analogue of an MBS.
- **Bolz, Schake, Stitz, Schroeder 2025**, *Expert Opin Drug Discov*, doi:10.1080/17460441.2025.2557599 (preprint doi:10.1101/2025.03.04.641378). [abstract/summary]
  - Studied PPI-inhibitor drugs.
  - The residues the drug contacts overlap the native PPI interface residues by 30–89%.
  - The PLIP 2025 release (Schake et al., *NAR* 2025, doi:10.1093/nar/gkaf361) adds protein–protein interaction profiling. This is useful for annotating interface ligands.
- **Moine-Franel et al. 2024** (above): a PPI pocket dataset with Zenodo data under CC BY 4.0 (doi:10.5281/zenodo.10805580).
- **Cryptic pockets:**
  - CryptoSite: Cimermancic et al. 2016, *J Mol Biol*, doi:10.1016/j.jmb.2016.01.029.
  - PocketMiner: Meller et al. 2023, *Nat Commun*, doi:10.1038/s41467-023-36699-3. ROC-AUC 0.87; suggests over half of "pocketless" proteins have cryptic pockets. [abstract/summary]
  - These are relevant because interface pockets in dimer models may open or close relative to monomer models.
- **Druggability predictors** (all Crossref-verified):
  - fpocket: Le Guilloux et al. 2009, *BMC Bioinf*, doi:10.1186/1471-2105-10-168.
  - fpocket drug score: Schmidtke & Barril 2010, *J Med Chem*, doi:10.1021/jm100574m.
  - DoGSiteScorer: Volkamer et al. 2012, *Bioinformatics*, doi:10.1093/bioinformatics/bts310.
  - PockDrug: Borrel et al. 2015, *JCIM*, doi:10.1021/ci5006004.
  - SiteMap Dscore: Halgren 2009, *JCIM*, doi:10.1021/ci800324m. Commercial (Schrödinger).
  - Practical choice: fpocket drug score plus P2Rank (open source, and best on LIGYSIS, see (e)). Use PockDrug or DoGSiteScorer as an orthogonal check.
- **Gap:** no large-scale study compares the druggability-score distributions of MBS vs SBS pockets in experimental PDB structures. Abrusán & Marsh's pharmacological claims were never tested with druggability metrics. I did not find one; [UNVERIFIED absence].

---

## (d) Predicted complexes and pocket atlases

### Predicted complexes
- **Evans et al. 2021**, AlphaFold-Multimer, bioRxiv doi:10.1101/2021.10.04.463034.
- **Humphreys et al. 2021**, *Science*, doi:10.1126/science.abm4805. Core eukaryotic (yeast) complexes.
- **Burke et al. 2023**, *Nat Struct Mol Biol*, doi:10.1038/s41594-022-00910-8. [abstract/summary]
  - 65,484 human PPIs modelled with AF2.
  - 3,137 high-confidence models, 1,371 with no homologous template.
- **Schweke et al. 2024**, *Cell*, doi:10.1016/j.cell.2024.01.022. [abstract/summary]
  - AF2 homo-oligomer atlas: >8,000 homomers across 4 proteomes.
  - Estimates ~45% of archaeal and bacterial proteomes and ~20% of eukaryotic proteomes are homomeric.
- **Zhang et al. 2025**, *Science*, doi:10.1126/science.adt1630. [abstract/summary]
  - RF2-PPI plus AF2 screen of ~200M human pairs.
  - 18,316 PPIs at 90% expected precision, 5,578 of them novel.
- **Pei, Zhang, Cong 2026**, *Sci Rep* 16:25590, doi:10.1038/s41598-026-54567-0. Case studies of those models; **no pocket analysis**.
- **Zhu et al. 2023**, *Bioinformatics*, doi:10.1093/bioinformatics/btad424. AF-Multimer evaluation and pDockQ2; per-interface confidence is useful for filtering.
- **Han et al. 2026**, bioRxiv doi:10.64898/2026.03.27.714458, "AlphaFold Database expands to proteome-scale quaternary structures". [abstract/summary] **Key new resource.**
  - ~23.4M homodimers from 4,777 proteomes, plus ~7.6M STRING-derived heterodimers.
  - Built with AF-Multimer through ColabFold/OpenFold.
  - High confidence = ipSAE_min ≥ 0.6, pLDDT_avg ≥ 70, clash ≤ 10. This gives about **1.75M high-confidence homodimers** and about 57k "tentatively high-confidence" heterodimers.
  - **No ligand or pocket analysis.**
  - Data on AFDB and the EBI FTP. Licence not stated in what I read; the AFDB is normally CC BY 4.0 [UNVERIFIED for complexes].
- **Bertoni et al. 2025**, *NAR*, doi:10.1093/nar/gkaf1226. AFDB 2025 update.

### Pockets in AF models and pocket atlases
- **Akdel et al. 2022**, *Nat Struct Mol Biol*, doi:10.1038/s41594-022-00849-w. Community assessment of AF2, including pocket and ligand-site prediction on AF2 models. The specific pocket numbers are [UNVERIFIED]; I did not read the full text.
- **Karelina, Noh, Dror 2023**, *eLife*, doi:10.7554/eLife.89386. [abstract/summary]
  - AF2 pockets are close to holo structures.
  - Docking poses into AF2 models are about as accurate as into traditional homology models, and worse than into experimental structures.
- **Hekkelman et al. 2023**, *Nat Methods* 20:205–213, doi:10.1038/s41592-022-01685-y (AlphaFill). [abstract/summary]
  - Transplants ligands and cofactors from PDB-REDO homologs.
  - 12.0M compounds transplanted into 586k AFDB models.
  - Useful as a "known-ligand prior" for predicted dimers, but monomer-based.
- **Comajuncosa-Creus et al. 2024**, *Nat Commun*, doi:10.1038/s41467-024-52146-3 (PocketVec). [abstract/summary]
  - >32,000 druggable sites in >20,000 human domains, using PDB plus AF2, fpocket/P2Rank and PRANK rescoring.
  - Pocket descriptors come from inverse virtual screening.
- **Carpenter & Altman 2025**, *J Cheminform*, doi:10.1186/s13321-025-01125-x (HOTPocket). [abstract/summary] 2.4M human pockets from 7 pocket finders; hotpocketNN scoring.
- **Zhang et al. 2025**, *Molecules* 30:260, doi:10.3390/molecules30020260. [abstract/summary] fpocket on all human AF2 models gives 15,043 druggable pockets in 20,255 structures.
- **Zillmer & Walther 2025**, bioRxiv doi:10.1101/2025.02.04.636483. Pocketome across 11 species from AFDB. [abstract/summary]
- CavitySpace (Wang et al. 2022, bioRxiv doi:10.1101/2022.01.25.477691). The PrePCI database was mentioned in search results but I did not verify it.
- **All of the above are monomer- or domain-based.** No atlas maps pockets that exist only in the dimer, i.e. MBS pockets, across predicted homodimers or heterodimers. I did not find one; see (h).

---

## (e) Datasets and benchmarks

| Name | What | Size | URL | License | Use for us |
|---|---|---|---|---|---|
| **BioLiP2** (Zhang C, Zhang X, Freddolino, Zhang Y 2024, *NAR* 52:D404, doi:10.1093/nar/gkad630) | Biologically relevant protein–ligand interactions, weekly updates | 385,160 chains; 781,684 interactions at publication | https://zhanggroup.org/BioLiP/ ; code at github.com/kad-ecoli/mmCIF2BioLiP | Data and code BSD; paper CC BY-NC | **High.** Its artifact filter is the standard (see (f)). **Caveat:** annotation is per chain, from the asymmetric-unit mmCIF, so MBS status must be recomputed on biounits. That is why Marsh dropped 13k entries. |
| **Q-BioLiP** (Wei et al. 2024, GPB, doi:10.1093/gpbjnl/qzae001) | BioLiP rebuilt on quaternary structures (biounits) | ~2.7M entries / ~210k PDB entries | https://yanglab.qd.sdu.edu.cn/Q-BioLiP/ | Paper CC BY 4.0; data licence [UNVERIFIED] | **Very high.** It has multi-chain ligand annotation natively. Use it to cross-check our MBS calls. It keeps both relevant and irrelevant ligands, with a 465-compound artifact list. |
| **LIGYSIS** (Utgés & Barton 2024, *J Cheminform* 2024, doi:10.1186/s13321-024-00923-z) | Ligand sites aggregated across biounits per UniProt; BioLiP-relevant ligands only; PISA preferred assemblies from PDBe; pdbe-arpeggio contacts; sites clustered by residue-fingerprint (average linkage, D=0.5) | ~25k–30k proteins; 64,782 sites; 435,038 ligands; human benchmark of 2,775 chains / 6,882 sites | https://www.compbio.dundee.ac.uk/ligysis/ | Article CC BY 4.0; dataset licence [UNVERIFIED]. Web server is "free, open to all, including commercial" | **High.** It is the best pocket-predictor benchmark and covers homomeric multi-chain sites. **Caveats:** (i) residues are mapped to UniProt numbering, so a homodimer MBS collapses onto one sequence; check whether chain identity is kept in the raw files. (ii) The server accepts homomultimers but not heteromultimers. Its benchmark found fpocket+PRANK (top-N+2 recall 60.4%), DeepPocket (58.1%) and P2Rank (≈52–54%) best. |
| FRAGSYS / Utgés et al. 2024 *Commun Biol*, doi:10.1038/s42003-024-05970-8 | Fragment-screen site classification | 37 screens, 293 sites | zenodo.org/records/10606595 | CC BY 4.0 | Low–medium. Fragment hot spots at interfaces. |
| **PLINDER** (Durairaj et al. 2024, bioRxiv doi:10.1101/2024.07.17.603955; ICML'24 ML4LMS workshop) | PLI systems per biounit. Ligands within 4 Å or PLIP-linked are merged. Pocket = PLIP-interacting residues plus residues within 6 Å. Has apo and predicted pairs, >500 annotations, and splits. | 449,383 systems; 25% multi-chain | https://github.com/plinder-org/plinder | Data Apache-2.0 (README) / CC BY 4.0 (paper); code Apache-2.0 | **Very high.** The biounit-based system definition directly gives the number of interacting protein chains. It has artifact, ion and cofactor flags and linked apo/AF2 structures, and supports leakage-aware splits for co-folding tests. The README mentions newer releases (e.g. "2026-09") [version not verified]. |
| **PINDER** (Kovtun et al. 2024, bioRxiv doi:10.1101/2024.07.17.603980) | Dimeric PPI systems with holo, apo and AF2 structures; interface-cluster splits | 2,319,564 dimer systems | https://github.com/pinder-org/pinder | Apache-2.0 [code; data licence UNVERIFIED] | **Medium–high.** Source of dimers with apo and predicted counterparts. Combine with PLINDER to get MBS sites with predicted-dimer analogues. |
| **Runs N' Poses** (Škrinjar et al. 2025, bioRxiv doi:10.1101/2025.02.03.636309) | Post-cutoff co-folding benchmark; similarity to training measured by `sucos_shape_pocket_qcov` | ~2,600 systems released after 2021-09-30 | github.com/plinder-org/runs-n-poses ; zenodo 14794786 | [UNVERIFIED] | **High** for (g). It tracks `num_protein_chains`, so it can be stratified into MBS vs SBS. |
| **PoseBusters** (Buttenschoen, Morris, Deane 2024, *Chem Sci*, doi:10.1039/D3SC04185A) | Physical-validity checks plus a 428-complex post-2021 set | 428 | github.com/maabuu/posebusters | BSD-3 [UNVERIFIED] | Medium. Use the validity checks on co-folding outputs. |
| PoseBench (Morehead et al. 2025, *Nat Mach Intell*, doi:10.1038/s42256-025-01160-1) | Apo/blind docking including multi-ligand (CASP15) | — | github.com/BioinfoMachineLearning/PoseBench | [UNVERIFIED] | Medium. Reports that DL methods generalise poorly to multi-ligand targets. |
| Binding MOAD (Hu et al. 2005, *Proteins*, doi:10.1002/prot.20512; sunset paper Wagle et al. 2023, *Sci Rep*, doi:10.1038/s41598-023-29996-w) | Curated biologically relevant ligands with affinities | — | bindingmoad.org | — | Low: **sunset, no further updates.** Its valid/invalid ligand calls are still useful as a second opinion. |
| sc-PDB (Desaphy et al. 2015, *NAR*, doi:10.1093/nar/gku928) | Druggable binding sites | ~16k sites [abstract-level] | — | — | Low: 1:1 complexes, old. |
| PDBbind (Wang et al. 2004, *J Med Chem*, doi:10.1021/jm030580l) | Affinity-annotated complexes | — | — | Recent versions need a licence [UNVERIFIED] | Low–medium. Affinity labels only. |
| Moine-Franel 2024 PPI pockets | See (c) | 23,238 pockets | zenodo 10.5281/zenodo.10805580 | CC BY 4.0 | Medium. Heterodimer PPI-modulator pockets. |
| ProtCID (Xu & Dunbrack 2020, *Nat Commun*, doi:10.1038/s41467-020-14301-4) | Interfaces recurring across crystal forms | — | dunbrack.fccc.edu/protcid | — | Medium. Separates biological homodimers from crystal contacts for PDB MBS calls. |
| AFDB complexes (Han et al. 2026) | AF-Multimer homodimers and heterodimers | ~1.75M high-confidence homodimers | alphafold.ebi.ac.uk | [UNVERIFIED for complexes] | **Very high.** Predicted-dimer source for (d) and (h). |
| AlphaFill | Transplanted ligands on AFDB monomers | 586k models | alphafill.eu | [UNVERIFIED] | Medium. Ligand priors. |

**Top recommendation:**
- Primary experimental source: **PLINDER**, optionally **Q-BioLiP**.
- Artifact filter: **BioLiP2**.
- Pocket-prediction benchmark: **LIGYSIS**.
- Predicted dimers: **AFDB complexes** plus PINDER.
- Co-folding evaluation: **Runs N' Poses**.

---

## (f) Additive and artifact exclusion lists to adopt

1. **BioLiP artifact list.** This is the primary list.
   - Source URL: `https://zhanggroup.org/BioLiP/ligand_list`.
   - 463 CCD codes as of today. That matches the BioLiP2 paper: 353 from 2013 plus 110 added later.
   - **BioLiP2 rule:** a listed ligand is discarded if **(i) it appears >15 times in the entry, or (ii) its binding site has only two consecutive residues**. Otherwise it is kept only if its name or a synonym appears in the PubMed title or abstract of the primary citation.
   - BioLiP2 binding residue: ≥2 heavy-atom contacts within vdW sum + 0.5 Å. A binding site is ≥2 binding residues.
   - **Adopt:** fetch the list programmatically, apply both geometric rules, and *optionally* the PubMed check. If we skip the PubMed check, flag it as an approximation.
2. **PLINDER artifact rules.**
   - Badlist: `src/plinder/data/annotations/static_files/artifacts_badlist.csv` in plinder-org/plinder, 255 codes curated in April 2024 against BioLiP, AF2 and RFAA artifact lists, with drug-like molecules and cofactors removed.
   - Generic filters:
     - fewer than 2 carbons or fewer than 5 heavy atoms;
     - absolute charge > 2;
     - unbranched hydrocarbon chain > 12 atoms;
     - unspecified atoms.
   - Single non-CHNOPS atoms are classed as "ions".
   - Small additives such as EDO, SO4 and ACT are absent from the badlist *because* the generic filters already catch them.
   - **Adopt as a second, stricter filter. Report agreement with BioLiP.**
3. **PLIP.** `plip/basic/config.py` embeds the BioLiP list dated 2014-07-10. It is used only for *warnings*, not exclusion, so it is not a separate source.
4. **Q-BioLiP:** 465-compound list (similar to BioLiP). **LIGYSIS:** uses BioLiP relevance directly.
5. **wwPDB "Ligand of Interest" flag** (`_pdbx_entry_details.has_ligand_of_interest`).
   - Depositor-supplied, so it is a strong positive signal.
   - Coverage is only ~22–30% of entries [abstract/summary from search]. Use it as a validation subset.
6. **AlphaFold 3 crystallization-aid exclusion list.** It exists in the AF3 supplementary material [UNVERIFIED; not read]. Worth harmonising with because it defines what AF3 was trained to ignore.

**Project rule (proposal):** keep a ligand if it passes BioLiP2's rules AND is not in the PLINDER badlist or caught by PLINDER's generic filters. Keep metals as their own class. Report how sensitive the MBS fraction is to the filter, since many additives (GOL, SO4, PEG) sit at crystal interfaces and would inflate MBS counts. All lists must be downloaded by script with the URL or commit pinned.

---

## (g) Co-folding evaluation and known pitfalls

- **Models:**
  - AlphaFold 3: Abramson et al. 2024, *Nature*, doi:10.1038/s41586-024-07487-w. PoseBusters success ~76% blind, ~90% with pocket specified [abstract/summary].
  - Boltz-1: Wohlwend et al. 2024, bioRxiv doi:10.1101/2024.11.19.624167.
  - Boltz-2: Passaro et al. 2025, bioRxiv doi:10.1101/2025.06.14.659707. Adds an affinity head.
  - Chai-1: Chai Discovery 2024, bioRxiv doi:10.1101/2024.10.10.615955.
- **CASP16 pharma ligands** (Gilson, Eberhardt, Škrinjar, Durairaj, Robin, Kryshtafovych 2025/2026, *Proteins* 94:249–266, doi:10.1002/prot.70061). [full text, skimmed]
  - 229 pose targets across 5 systems: chymase, cathepsin G, autotaxin, **SARS-CoV-2 Mpro (homodimer)** and WDR55.
  - Best CASP group: mean LDDT-PLI 0.69, all template-based.
  - Baselines run after the fact (not blinded, but checked for training leakage): mean LDDT-PLI **0.80 for AF3**, 0.52 for Boltz-1, 0.37 for RFAA, 0.38 for Vina. Boltz-1 underperformed specifically on autotaxin.
  - Success threshold RMSD ≤ 2.5 Å.
  - Pose accuracy correlated only weakly with SuCOS similarity to PDB ligands (R² 0.25).
  - Affinity: best Kendall τ 0.42; ligand MW alone was almost as good.
  - **Relevance:** all targets are well-characterised drug targets with many PDB analogues, so this tests template memory more than generalisation. It does not stratify by multi-chain sites. Metrics to reuse: LDDT-PLI, BiSyRMSD, BB-RMSD; scoring via OpenStructure.
- **Memorisation:** Škrinjar, Eberhardt, Durairaj, Schwede 2025 (Runs N' Poses).
  - Co-folding methods "largely memorise ligand poses from their training data". Accuracy falls with lower `sucos_shape_pocket_qcov` similarity to training.
  - **Cofactors and other promiscuous ligands do moderately better.**
  - Per-bin success rates were not extracted [UNVERIFIED].
  - **Implication:** cofactors are enriched in MBS sites (Marsh), so co-folding success on MBS sites is confounded with cofactor memorisation. We must stratify by ligand class and by training similarity.
- **Physics robustness:** Masters, Mahmoud, Lill 2025, *Nat Commun*, doi:10.1038/s41467-025-63947-5. [abstract/summary]
  - Mutating all binding-site residues within 3.5 Å to Gly or Phe often does **not** displace the predicted ligand.
  - **Implication:** a ligand placed by co-folding into a predicted dimer interface is not evidence that a pocket exists. Use controls such as single chain vs dimer input, ligand-shuffle decoys, and binding-site mutation.
- **Generalisation to multi-ligand targets:** PoseBench (Morehead et al. 2025) reports poor generalisation of DL docking to multi-ligand targets.
- **Other pitfalls for our use case:**
  1. **Template and MSA leakage:** use a post-cutoff time split (PLINDER-TIME, Runs N' Poses).
  2. **Symmetric homodimer ligand-copy ambiguity:** score with symmetry-aware chain mapping (OpenStructure / PLINDER scoring).
  3. **Input stoichiometry:** for an MBS site, give one ligand copy per site, not one per chain.
  4. **Apo vs holo bias:** AF models tend toward holo-like pockets (Karelina 2023), which may *inflate* pocket detection in predicted dimers compared with apo crystal structures.

---

## (h) Open gaps: what would be new here

1. **MBS vs SBS druggability, measured.**
   - Abrusán & Marsh predicted different pharmacology (selective vs multitarget) but never computed druggability.
   - Novel: compare fpocket, P2Rank and PockDrug scores, plus drug-like ligand occupancy (ChEMBL/DrugBank via UniChem), for MBS vs SBS sites in a current PDB, with a modern artifact filter and biounits.
2. **An MBS atlas of predicted dimers.**
   - Existing pocket atlases (PocketVec, HOTPocket, AlphaFill, CavitySpace, Molecules 2025) are monomer-based.
   - The AFDB complexes release (~1.75M high-confidence homodimers, plus heterodimers) and the human interactome models (Burke 2023, Zhang 2025) have **not been mined for interface pockets**, i.e. pockets that exist only in the dimer.
   - Novel: proteome-scale "dimer-only" pockets for human and other species, with druggability.
   - **Correction (design review r1, reviews/r1_prior_work.md):** ProtPocket (GitHub, 2026; not peer-reviewed; found by the prior-work reviewer, not independently re-checked here) already compares monomer and dimer pocket runs on the March 2026 AFDB homodimer release. The claim "no dimer pocket atlas" is therefore softened. What remains distinct: validation against PDB MBS sites, heterodimers, cross-species comparison and ligand-class analysis.
3. **Do predicted-dimer pockets reproduce real PDB MBS sites?**
   - No systematic benchmark exists. LIGYSIS covers monomers and homomers mapped onto UniProt, and has no predicted complexes.
   - Novel: a recall/precision benchmark of AF-Multimer (and AF3) dimer pockets against PLINDER/Q-BioLiP multi-chain sites. Stratify by ipTM/ipSAE, homo vs hetero, and cofactor vs drug-like.
4. **Ligand-class composition revisited.**
   - Marsh's cofactor enrichment comes from BioLiP 2013 and asymmetric-unit-limited data. Q-BioLiP gives 28.8% and PLINDER 25% multi-chain, but neither breaks this down by cofactor, metabolite or drug-like ligand, by species, or by homo vs hetero.
   - Novel: an updated, class-stratified census, with heteromers treated properly. Marsh flagged heteromer results as biased by incomplete PDB complexes, and predicted heterodimers can fill that gap.
5. **Co-folding re-evaluation at MBS sites, controlling for memorisation.**
   - No study reports AF3, Boltz or Chai accuracy separately for multi-chain vs single-chain sites.
   - Novel: stratify Runs N' Poses / PLINDER-TIME by MBS status × training similarity × ligand class.
   - Add a monomer-input vs dimer-input control: does the ligand still dock when the second chain is removed?
6. **Secondary angle:** heterodimer MBS sites carry the highest pathogenic-variant burden (~30% in Marsh's data). Mapping ClinVar/gnomAD onto predicted heterodimer MBS pockets could prioritise interface-pocket variants.

**Absence claims in (c), (d) and (h)** come from targeted web searches done on 2026-10-29. They are "not found", not proven absent. Re-check bioRxiv before writing the paper.
