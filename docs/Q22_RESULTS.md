# Q2.2 / Q2.3a results: dimer pockets vs monomer pockets (human)

**Status: UNVERIFIED. EXPLORATORY.** None of this is a pre-registered hypothesis (DESIGN.md v3 §1 lists H1–H5 only;
Q2.2 and Q2.3a are exploratory per §3). The numbers below have not been re-derived independently, and no figure has been reviewed by the PI.

- DESIGN.md sha256 `e3f9b550eb29cfe6eb444ebd884920b847687518470834c0bcee3439998eee28` (equal to tag `design-v3-frozen`)
- Outputs were generated at HEAD `1113101950298ac1860497a1804c55aa71b0b1f5`. The q22 scripts were still uncommitted
  then and were committed afterwards, in the commit that adds this file. Every TSV carries the design hash and the HEAD in its first line.
- Seeds: bootstrap seed 0 is primary. Seeds 1 and 2 were rerun for the "all" stratum of analyses 1 and 4; their CIs differ by ≤0.005
  (Cliff's delta) and ≤0.0005 (cosine). 1,000 cluster-bootstrap resamples (200 for seeds 1–2 of the size-matched novelty).
- Bootstrap unit: the paralog cluster of chain A's protein. A cluster is a connected component of the mmseqs2 all-vs-all graph with
  ≥30% identity, e ≤1e-3 and **≥80% coverage of both sequences**. This gives 4,577 units, the largest with 69 proteins (315 models).
  Without the coverage condition, single linkage through shared domains forms one giant component (2,172 proteins, 31% of models), which would make the bootstrap degenerate.
  The novelty exclusion (analysis 4) deliberately uses the looser rule with no coverage condition, which is conservative.

## Pipeline and n per step

| step | n | script |
|---|---|---|
| human dimer models with pockets (homo 4,803 / hetero 5,785) | 10,588 | (input) `data/part2/human_dimer_pockets_annotated.parquet` |
| dimer pockets (all) | 87,804 | |
| residues with features (FreeSASA rel. SASA chain-alone and in complex; interface = heavy atom ≤5 Å of partner) | 8,124,263 (0 failed models) | `q22_residue_features.py` |
| surface residues (rel. SASA ≥0.25 in the isolated chain) / interface residues | 61.6% / 13.4% of residues | |
| atlas-pass pockets (burrow_score ≥0.03, lining pLDDT ≥70) | 69,857 | `q22_compare.py` |
| … symmetric-copy representatives | 54,876 | |
| … excluded "marginal" pockets (1–2 lining residues on one chain) | 6,763 | |
| … interface (≥3 residues per chain) / single-chain (0 residues on one chain) | 17,882 / 30,231 | |
| isolated-chain pockets, atlas-pass (homodimers: chain A only) | 58,589 (models with isolated-chain output: 10,562 of 10,588) | |
| human monomer atlas pockets (taxid 9606, tiers 1+2; 17,473 accessions) | 90,560 | `q22_atlas_human.py` |
| UniProt sequences (proteome UP000005640 + 1,361 by accession; release 2026_03) | 21,039 of the 21,050 needed (11 missing) | `q22_fetch_uniprot.py` |
| mmseqs2 hits (≥30% id, e ≤1e-3, -s 7.5) | 930,462 | `q22_paralogs.sh` |

## 1. Interface vs single-chain pockets (same models)

Cliff's delta is interface vs single-chain, with 95% cluster-bootstrap CIs (`results/part2/q22/q22_props.tsv`; figures `q22_pocket_properties.png`, `q22_props_by_ipsae.png`).

| metric | median interface / single | Δ median [CI] | Cliff's δ all [CI] | homo | hetero |
|---|---|---|---|---|---|
| volume_A3 | 101 / 77 Å³ | +24 [20, 30] | **+0.23** [0.20, 0.26] | +0.23 | +0.23 |
| # lining residues | 15 / 12 | +3 [2, 4] | **+0.32** [0.28, 0.35] | +0.35 | +0.29 |
| burrow_score | 0.096 / 0.076 | +0.021 [0.016, 0.026] | +0.09 [0.07, 0.12] | +0.12 | +0.07 |
| buriedness (1 − rel. SASA in complex) | 0.779 / 0.766 | +0.013 | +0.10 [0.07, 0.13] | +0.14 | +0.07 |
| lining pLDDT | 91.2 / 91.2 | −0.06 [−0.54, 0.38] | +0.01 [−0.02, 0.05] | +0.06 | −0.03 (CI spans 0) |
| hydrophobic fraction | 0.36 / 0.42 | −0.053 | **−0.15** [−0.18, −0.11] | −0.13 | −0.16 |
| aromatic fraction | 0.16 / 0.19 | −0.026 | −0.12 [−0.14, −0.09] | −0.10 | −0.13 |
| charged fraction (DEKR) | 0.21 / 0.15 | +0.060 | **+0.20** [0.17, 0.24] | +0.18 | +0.23 |

Interface pockets are somewhat larger, slightly better scored and more buried. Their lining is more charged and less
hydrophobic/aromatic, which is the composition expected of a PPI interface. All effects are small to moderate (|δ| ≤0.35).
By ipSAE tertile (edges 0.714 / 0.807): the size and composition effects hold in all three tertiles. In the top tertile the
burrow_score and site-volume differences shrink to CIs that span 0, and the lining-pLDDT deficit of interface pockets is seen only
in the lowest tertile (δ −0.10 [−0.15, −0.05]).

## 2. Expected vs observed interface pocket rate

`q22_expected_rate.tsv`, `q22_expected_vs_observed.png`. Pooled ratios with a cluster bootstrap over 10,588 models.

| | all | homo | hetero |
|---|---|---|---|
| expected: interface share of surface residues | 0.167 [0.158, 0.178] | 0.180 | 0.155 |
| observed: interface share of pockets (all copies) | 0.315 [0.301, 0.328] | 0.322 | 0.307 |
| observed: interface share of pocket-lining surface residues | 0.419 [0.403, 0.440] | 0.413 | 0.426 |
| enrichment, pocket level (obs/exp) | **1.88** [1.78, 1.99] | 1.79 [1.71, 1.88] | 1.97 [1.79, 2.16] |
| enrichment, residue level | **2.51** [2.43, 2.61] | 2.29 [2.21, 2.37] | 2.74 [2.60, 2.90] |

Rates by residue: 26.6% of interface surface residues line an atlas-pass pocket, against 7.4% of non-interface surface residues.
Interface pockets are therefore about 2× more frequent than the interface share of the surface predicts, and the enrichment holds in every ipSAE tertile (1.70–1.97).
Caveat: the pocket-level comparison sets a count of pockets (a patch-level quantity) against a residue fraction. A pocket only needs ≥3 residues per chain
to count as interface, so a random patch would also be inflated above the residue fraction. The residue-level ratio is the fairer like-for-like
contrast. Residues in a groove formed by two chains are concave by construction, which is part of why they are enriched.

## 3. Dimer-only pockets and pocket destruction

**Raw dimer-only fraction** among interface pockets (`q22_dimer_only_frac.tsv`), defined as max overlap coefficient with any isolated-chain pocket below the threshold:
OC <0.2: 0.37 [0.36, 0.38]; **OC <0.5 (primary): 0.635 [0.615, 0.656]**; OC <0.8: 0.80 [0.78, 0.81].
Homo vs hetero at 0.5: 0.651 [0.637, 0.667] vs 0.625 [0.594, 0.658], so no clear homo/hetero difference. The fraction is stable across ipSAE tertiles (0.62–0.65).

**Method-noise null (important).** BURROW output depends on context. Of the isolated-chain pockets with *no* lining residue at the
interface, which see identical local coordinates in both runs, 34% have no dimer counterpart at OC ≥0.5. The loss depends strongly on score: 64% in the lowest burrow_score
quintile against 6% in the highest. burrow_score also changes between runs: for example, an isolated-chain pocket scoring 0.058 corresponded to a dimer pocket scoring 0.023.
`q22_dimer_only_control.py` therefore uses single-chain dimer pockets as the null for "new in dimer", with homodimer chain-swap matching:

| threshold | interface "new" rate | single-chain null rate | excess, raw [CI] | excess, burrow_score-quintile matched [CI] |
|---|---|---|---|---|
| 0.2 | 0.347 | 0.094 | +0.25 [0.24, 0.27] | +0.27 [0.25, 0.28] |
| **0.5** | **0.615** | **0.166** | +0.45 [0.43, 0.47] | **+0.47 [0.45, 0.49]** |
| 0.8 | 0.782 | 0.296 | +0.49 [0.47, 0.50] | +0.50 [0.49, 0.52] |

The homo and hetero values are within ±0.04 of these. The dimer-only signal therefore survives the noise floor, but roughly a quarter of the primary
dimer-only calls (0.166 / 0.615) are of the order expected from method noise alone.

**Dimer-only vs pre-existing interface pockets** (`q22_props_dimer_only.tsv`). Contrary to the "larger, more buried" expectation, dimer-only pockets are
*smaller* (volume δ −0.47 [−0.50, −0.44]; median 83 vs 188 Å³), have fewer lining residues (δ −0.55), score lower (burrow_score δ −0.47)
and are slightly *less* buried (δ −0.10 [−0.13, −0.07]). Lining chemistry does not differ (|δ| <0.04). Because low-score pockets are also the noisiest
(see the null above), this contrast is partly confounded by noise.

**Destruction of monomer pockets** (`q22_destruction.tsv`, `q22_destruction_control.tsv`, `q22_controls.png`). Isolated-chain pockets with ≥3 lining residues
at the interface are lost (OC <0.5 with every dimer pocket; homodimer chain-swap matched) in 41.5% of cases, against 31.1% for no-contact pockets. The
score-matched excess loss is **+0.13 [0.09, 0.17]** (homo +0.09 [0.08, 0.12], hetero +0.14 [0.10, 0.19]). A further 15% of partner-contact
pockets shrink (volume ratio <0.5); the no-contact rate is similar at 14%. So partner occlusion destroys a modest excess of monomer pockets above a large
method-noise background.

## 4. Novelty vs the human monomer atlas (pocketcompare)

`q22_novelty.py` (torch on cuda:1, exact inner products, 87,804 × 90,560); `q22_novelty.tsv`, `q22_novelty.png`. For each dimer pocket we take the max cosine to atlas
pockets, excluding both proteins of the model and all of their mmseqs2 hits at ≥30% identity (median 35 atlas pockets excluded per query; 17,720 proteins have ≥1 paralog).
The dimer pocket and atlas method records are identical (`pocketcompare-final:A4b:n1_linear|raw`).

| contrast | metric | raw Δ median [CI] | Cliff's δ [CI] | volume-decile matched Δ [CI] |
|---|---|---|---|---|
| interface − single-chain | max cosine (excl.) | −0.010 [−0.013, −0.008] (0.801 vs 0.811) | −0.13 [−0.16, −0.11] | −0.010 [−0.013, −0.007] |
| interface − single-chain | per-query z | +0.07 [0.04, 0.10] | +0.08 [0.04, 0.11] | +0.05 [0.03, 0.08] |
| interface − single-chain | max cosine to *own* atlas pockets | −0.27 (0.47 vs 0.74) | −0.60 [−0.62, −0.58] | n/a |
| dimer-only − pre-existing interface | max cosine (excl.) | −0.003 [−0.005, 0.000] | −0.02 [−0.06, 0.01] | 0.000 [−0.003, 0.003] |
| dimer-only − pre-existing interface | per-query z | −0.10 | −0.11 [−0.14, −0.07] | −0.08 [−0.12, −0.04] |

By the raw max cosine, interface pockets are marginally *less* similar to known monomer pockets of unrelated proteins. The difference is ~0.01 on a scale
where the median nearest neighbour is 0.80, so it is negligible in practice. By the per-query z-score the direction reverses: interface pockets have lower similarity to the atlas on average, so their top hit stands out more.
Homodimers drive the z effect (+0.13 matched); in heterodimers it is null (+0.004 [−0.03, 0.04]). The large effect is the match to the protein's *own*
monomer pockets (δ −0.60), which is expected because interface pockets mostly lack a monomer counterpart. Dimer-only pockets do not differ from pre-existing interface pockets in novelty once
size is matched. **Interpretation caveat:** pocketcompare recognises PPI-interface pockets poorly across families (docs/TOOLS.md §3). A low maximum cosine is
therefore not evidence of chemically novel pockets, and these are descriptor-space statements only.
**Positive control (PDB interface ligand sites): skipped.** No pocketcompare fingerprints of Part 1 PDB sites exist yet. The control can be added
by running `pocketcompare describe --ligand` on Part 1 MBS sites.

## 5. Q2.3a: proteins whose only drug-sized confident pocket is dimer-only

`q22_q23a.py` → `results/part2/q23a_dimer_only_ligandable.tsv`. **Drug-like size** = volume_A3 300–1500 Å³ (task FALLBACK: no BURROW
volumes of PDB drug-like ligand pockets are available yet). This range selects the top ~10% of atlas-pass dimer pockets (atlas: 9.1%). site_volume_A3 was not
used because 76% of pockets fall in the range (the sph8 core is capped). A protein qualifies when it has ≥1 interface, dimer-only (OC <0.5), atlas-pass,
drug-sized dimer pocket and *none* of the following has an atlas-pass drug-sized pocket: its monomer atlas pockets (any fragment), its isolated-chain pockets in
any of its models, or single-chain dimer pockets on its chain.

| step | n |
|---|---|
| qualifying dimer-only drug-sized pockets | 1,028 |
| proteins with ≥1 qualifying pocket | 1,125 |
| excluded: atlas drug-sized pocket / isolated-chain / single-chain dimer (overlapping) | 271 / 638 / 348 |
| **listed** (in atlas structure table / not in atlas) | **393** (259 / 134) |
| sensitivity: volume_A3 150–1500 Å³ | 462 listed (284 in atlas); 112 of the 259 primary in-atlas proteins are retained |

The listed set is sensitive to the size window: under 150–1500 Å³ more proteins gain a drug-sized monomer pocket and drop out.
ChEMBL status (REST, ChEMBL_37), primary in-atlas list (n = 259): drug_target 31, potent_ligands (≥1 pChEMBL ≥6 record on a single-protein target)
53, chembl_target 38, not_in_chembl 137. "drug_target" counts mechanisms on *any* ChEMBL target that contains the protein, including complexes such as the
nicotinic receptor for CHRNB1. Hetero 216 / homo 43.

Top 10 (rank by best burrow_score; all in atlas):

| rank | protein (UniProt) | partner | kind | volume Å³ | burrow_score | ipSAE | ChEMBL |
|---|---|---|---|---|---|---|---|
| 1 | RGS4 (P49798) | GNAT2 | hetero | 1134 | 0.56 | 0.87 | potent_ligands |
| 2 | ATP6V1D (Q9Y5K8) | ATP6V0D2 | hetero | 609 | 0.54 | 0.69 | not_in_chembl |
| 3 | CHRNB1 (P11230) | CHRNA1 | hetero | 553 | 0.51 | 0.61 | drug_target (complex) |
| 4 | RGS13 (O14921) | GNAI2 | hetero | 818 | 0.46 | 0.84 | not_in_chembl |
| 5 | CREG1 (O75629) | CREG1 | homo | 310 | 0.46 | 0.90 | not_in_chembl |
| 6 | TMEM60 (Q9H2L4) | TMEM203 | hetero | 306 | 0.43 | 0.73 | not_in_chembl |
| 7 | MYCN (P04198) | MAX | hetero | 643 | 0.42 | 0.71 | chembl_target |
| 8 | MAX (P61244) | MYCN; MYCL; MXI1; MAX | hetero | 643 | 0.42 | 0.71 | chembl_target |
| 9 | MXI1 (P50539) | MAX | hetero | 583 | 0.41 | 0.66 | not_in_chembl |
| 10 | IFNLR1 (Q8IU57) | IFNL3 | hetero | 1331 | 0.40 | 0.76 | drug_target |

(Rank 13 is TLR4 with LY96/MD-2, whose interface hosts the known LPS/eritoran site. This is anecdotal and not a validation.)

## Caveats

- Everything here is exploratory and descriptive. There is no multiplicity correction because there is no hypothesis test.
- Dimer-only calls carry a substantial method-noise component (null rate 0.17 at OC 0.5). Low-score dimer-only pockets are the least reliable.
  The annotated `dimer_only` flag does not use homodimer chain-swap matching, which changes the rate from 0.635 to 0.615.
- Confident interfaces only (ipSAE ≥0.6, pDockQ2 ≥0.23). Predicted models carry AF training overlap. Whether a PDB complex exists for the pair was not stratified here.
- Residue-class definitions are ours: hydrophobic AVLIMFWC, aromatic FWYH, charged DEKR. Surface = rel. SASA ≥0.25 (FreeSASA default classifier).
- ChEMBL: the local SQLite query (30 GB on the shared filesystem) was stopped after >35 min of IO wait, and REST was used for the 393 listed proteins only.
  In the sensitivity list, ChEMBL columns are filled only for proteins that are also in the primary list.
- 11 accessions without a UniProt sequence were treated as having no paralogs.

## Files
Scripts: `scripts/q22_common.py`, `q22_atlas_human.py`, `q22_fetch_uniprot.py`, `q22_paralogs.sh`, `q22_residue_features.py`, `q22_pocket_features.py`,
`q22_compare.py`, `q22_dimer_only_control.py`, `q22_novelty.py`, `q22_chembl_targets.py`, `q22_q23a.py`, `q22_figures.py`.
Rerun order: atlas_human → fetch_uniprot → paralogs → residue_features → pocket_features → compare → dimer_only_control → novelty (pockets venv, GPU 1)
→ q23a → chembl_targets → q23a → figures. Environments: `.venvs/q22` (freesasa, matplotlib), `.venvs/pockets` (torch), `conda_envs/mmseqs` (mmseqs2 18.8cc5c).
Results: `results/part2/q22/*.tsv`, `results/part2/q23a_dimer_only_ligandable.tsv`, `results/part2/figures/q22_*.png`. Intermediates: `data/interim/q22/` (not committed).

## 6. Curated list v2 (Q2.3a, EXPLORATORY, post-hoc)

Status: exploratory and post-hoc relative to the frozen design; filters were declared (here and in the header of `scripts/q23a_curate.py`) before looking at
which proteins pass, in response to the artefacts listed in `reviews/verify_q22.md` (helical-groove bHLH-LZ pockets, low scores, low lining pLDDT, 3 residues on one chain,
double counting). Not verified; no multiplicity or hypothesis test.

Declared filters (all must hold; unit = dimer pocket, atlas-pass rows of `dimer_pocket_features`):
1. burrow_score >= 60th percentile of all atlas-filtered dimer pockets.
2. Lining pLDDT mean >= 80 on both chains.
3. >= 5 lining residues on each chain.
4. Helical fraction of lining residues < 0.7 (DSSP via pydssp, numpy backend, run per chain on the dimer model; helix = pydssp C3 helix class).
5. volume_A3 300-1500.
6. max ipSAE >= 0.7.
7. Dimer-only: overlap coefficient < 0.5 vs isolated-chain pockets (existing `max_overlap_iso`, both chains) and < 0.5 between each chain's lining residues and every AFDB
   monomer-atlas pocket of that protein (vacuous, and flagged `atlas_has_protein=False`, when the atlas lacks the protein).
Dedup: symmetry representatives only; pockets of one unordered protein pair with residue overlap >= 0.5 merged (best score kept). One row per dimer pocket listing both proteins;
protein counts are unique proteins. Annotation: ChEMBL status from the existing Q2.3a table (`not_in_q23a_table` otherwise), PDB complex = pair shares a PDB entry in the SIFTS chain table
(`n_pdb_complex_in_part1` = those also in `data/part1/entries.parquet`). Output: `results/part2/q23a_curated.tsv`, step counts `results/part2/q22/q23a_curated_steps.tsv`.

Result (first attempt, not verified). burrow_score 60th percentile = 0.113 (69,857 atlas-pass dimer pockets). Filter cascade (n pockets): atlas-pass 69,857 → score 27,943
→ lining pLDDT 10,542 → ≥5 res/chain 5,942 → volume 2,038 → ipSAE 1,696 → iso-overlap 299 → sym-rep 238 → atlas-overlap 165 → helix <0.7 114 → dedup 114.
**114 pockets, 137 unique proteins** (79 homo / 35 hetero; 32 with a PDB complex for the pair). Caveats: 22 rows are isoform models (lookups by base accession; atlas-overlap uses
UniProt positions where mapped); 33 rows involve a protein absent from the monomer atlas (monomer-overlap test vacuous); many accessions are TrEMBL IDs for known genes, so the
ChEMBL column is `not_in_q23a_table` for 81 of 114 (the old table covered only 393 proteins; this is "not annotated", not "no ligands"). The list is dominated by homodimeric enzymes,
so the filter set selects likely active-site/oligomer-interface clefts of well-folded enzymes. Whether these are novel, ligandable pockets is not tested.

Top 10 manual inspection (lining residue ranges are model numbering; plausibility is from general knowledge, not checked against literature here):
1. DPYD homodimer (1038 Å3): 13 + 37 residues (one chain's residues spread over 82-828), plausible catalytic-region cleft of a large FMN/FAD/Fe-S enzyme; PDB complexes of DPYD exist in other species, not flagged for this pair.
2. GCAT homodimer (617 Å3): 27 + 10 residues across the aminotransferase domain; PLP-dependent enzyme with classical dimer-spanning active site, so likely a known catalytic site.
3. STING1 homodimer (302 Å3): 14 + 16 residues, CTD 159-270, this is the cGAMP/ligand-binding cleft of the STING dimer (PDB 6MX3 exists); recovers a known ligandable site, borderline helix fraction 0.63.
4. SGPL1 homodimer (449 Å3): 35 + 35 residues, PLP lyase dimer active site; likely known catalytic cleft (accession is a TrEMBL ID).
5. IVD homodimer (382 Å3, isoform model): 17 + 12 residues; acyl-CoA dehydrogenase, substrate/FAD pocket at a dimer interface; PDB complex exists.
6. ACADVL (578 Å3, isoform model): 12 + 15 residues; acyl-CoA dehydrogenase dimer, substrate channel.
7. AASS (444 Å3): 16 + 20 residues; bifunctional lysine-ketoglutarate reductase/saccharopine dehydrogenase, domain-spanning cleft.
8. OGDHL (548 Å3, isoform): 24 + 11 residues; 2-oxoglutarate dehydrogenase-like, thiamine-diphosphate dimer interface site.
9. GCDH (690 Å3): 19 + 24 residues; acyl-CoA dehydrogenase family, substrate/FAD cleft at dimer interface.
10. CUTC (365 Å3): 22 + 23 residues, symmetric positions 33-102; copper homeostasis protein, a symmetric dimer cleft, unknown function, least characterised of the ten.
Overall: 9 of the top 10 are homodimeric metabolic enzymes whose dimer-spanning clefts are probably canonical active sites. This is a sanity check that the filters enrich for functional
sites, but it also means high-ranking entries are mostly not novel allosteric pockets. Hetero-dimer entries start at rank 16 (e.g. AMIGO2-AMIGO3 rank 17, ATP6V1B1-ATP6V1A rank 19, RSPO4-LGR4 rank 22).
Files: `scripts/q23a_curate.py`, `results/part2/q23a_curated.tsv`, `results/part2/q22/q23a_curated_steps.tsv`.
