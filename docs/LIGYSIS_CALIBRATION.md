# LIGYSIS calibration of BURROW (EXPLORATORY, UNVERIFIED)

Status: first-attempt analysis, not reviewed. Numbers come from `scripts/ligysis_0*.py`; outputs in `results/benchmarks/ligysis/`. Not yet checked against any independent computation.

## What LIGYSIS's 60.4% is (Utgés & Barton 2024, J Cheminform, PMC11552181, read via Europe PMC)
- Top-(N+2) recall, N = number of observed LIGYSIS sites of the target protein chain; **hit = DCC <= 12 Å** (their recommended threshold; DCC 4 Å is called "too conservative"; DCA 4 Å is their other distance metric).
- Benchmark: 2,775 human representative experimental chains, 6,882 sites (ligand-stripped single chains; ions included). fpocket+PRANK 60.4%, DeepPocket RESC 58.1%, P2Rank CONS 53.9%, P2Rank 51.9%. The file used here has 10,689 sites / 3,518 proteins (April 2024 release), so the denominators differ and the 2,775-chain subset cannot be identified from the local files.
- Consequence: the earlier project comparison (DCC <= 4 Å) was not comparable; the like-for-like threshold is DCC <= 12 Å.

## Method deviations (important)
1. **Site centre is a proxy**: LIGYSIS site centres/ligand centroids sit in a PDBe superposition frame that is not available locally (checked: the `centre` column is not in the PDB asymmetric-unit frame). The centre is approximated by the centroid of the CA atoms of the site's UniProt residues in the AFDB model (no superposition needed). The proxy is not the ligand centroid, so DCC values are not identical to LIGYSIS's. DCA (closest ligand atom) was not computed.
2. Pockets = atlas pockets (tier 1, AF mode, pLDDT filters, score >= 0.03; ~4.3 pockets per protein), ranked by `burrow_score`, protein level. N = LIGYSIS sites of the UniProt (`acc`); sensitivity with N = sites of the structural segment (`seg`). Fragments of long proteins: site assigned to the fragment holding most of its residues; DCC only to pockets of that fragment; overlapping fragment windows are not de-duplicated.
3. Residue criterion: overlap coefficient |S∩P|/min(|S|,|P|) >= 0.5 in UniProt numbering (the 0.5 rule of docs/IFACE_LIGAND_MATCHED.md, not a LIGYSIS criterion).
4. AFDB monomer models, not experimental chains.

## Result 1: AFDB atlas recall (top-N+2), all LIGYSIS human sites
Denominator C = proteins with >= 1 atlas pocket and a computable site centre (10,207 sites, 3,377 proteins); A = all sites of proteins with an AFDB structure (10,682); proteins without any kept pocket count as misses in A. CIs: protein-cluster bootstrap.

| criterion | A (all) | B (>=1 pocket) | C |
|---|---|---|---|
| DCC_proxy <= 12 Å, top-N+2 | 49.3% [47.9, 50.6] | 50.2% | 51.6% [50.2, 52.9] |
| DCC_proxy <= 4 Å, top-N+2 | 13.0% | 13.3% | 13.6% |
| residue overlap >= 0.5, top-N+2 | 36.5% | 37.2% | 38.2% |
| DCC_proxy <= 12 Å, all pockets | 51.5% | 52.5% | 53.9% |

(`seg`-level N: 48.4% / 36.0% for A; see `ligysis_recall_summary.tsv`.) Against the published 60.4% / 51.9% (fpocket+PRANK / P2Rank) this is below fpocket+PRANK and close to P2Rank, but the denominators, structures (AFDB vs experimental), centre proxy and atlas pocket filters all differ, so this is a rough context, not a ranking.

## Result 2: splits (B denominator, `acc`, DCC_proxy <= 12 Å / overlap)
LIGYSIS annotates which chains contact each ligand (`LIGYSIS_human_chains_per_lig_MAY_2024.pkl`); there is no explicit multi-chain flag, so it is derived here.
- single-chain sites (every supporting ligand contacts 1 chain; 8,236): 48.2% / 35.2%
- multi-chain, any ligand contacts >1 chain (2,242): 57.8% / 44.2%; all ligands multi-chain (1,379): 48.6% / 31.6%
- ion-only sites (4,142): 35.6% / 18.1%; non-ion-only (6,336): 59.8% / 49.7%
- site size: <= 5 residues 27.8%; >= 16 residues 85.9% (DCC 12 Å is easier for large sites; size confounds the multi-chain split since multi-chain sites are larger)
- site pLDDT < 70 (187 sites): 15.5%.
Interpretation of multi-chain vs single-chain is confounded by size and by ion content; not adjusted.

## Result 3: like-for-like on experimental chains (200-chain sample, seed 0)
BURROW via `pocketcompare describe` (merged pockets, non-AF mode) on 200 LIGYSIS representative chains cut out of the local PDB mmCIFs (480 sites, 472 evaluable; sampled from segments whose DEF_TRANS and human site counts agree and whose UniProt has atlas pockets), same centre proxy, score >= 0.03. Paired on the 467 sites evaluable in both settings (seg-level N):

| criterion | experimental chains | AFDB atlas |
|---|---|---|
| DCC_proxy <= 12 Å | 57.6% | 48.4% |
| DCC_proxy <= 4 Å | 13.5% | 12.0% |
| overlap >= 0.5 | 37.7% | 33.4% |

Experimental chains give ~9 points higher recall at 12 Å (sample of 192 chains; no CI computed; PDB chains are the holo-derived ligand-stripped structures used to define the sites, which favours them). The sample is small and the choice of sample was not stress-tested.

## Caveats
- Everything depends on the centre proxy; a ligand-centroid DCC would be needed for a strict reproduction.
- Atlas pockets are already filtered (pLDDT mean >= 70, score >= 0.03), so recall in the atlas is limited by that filter ("max" recall 51.5-53.9% at 12 Å).
- The 2,775-chain LIGYSIS benchmark subset is not identified; sites here include all 10,689.
- Run on 2026-10-29/30, single run; no repeats.

## Reproduce
```
source scripts/env_tools.sh pockets
python scripts/ligysis_01_extract_atlas.py; python scripts/ligysis_01b_extract_models.py
python scripts/ligysis_02_recall_afdb.py; python scripts/ligysis_03_summarise.py
python scripts/ligysis_04_burrow_exp.py;  python scripts/ligysis_05_paired.py
```
