# Q2.4 results: conservation of human dimer pockets in ortholog dimers

**Status: UNVERIFIED. EXPLORATORY** ("could" priority, DESIGN.md §3). This is a first attempt by one agent. Nothing has been re-derived independently, and the PI has not reviewed any figure.
No hypothesis test and no multiplicity correction. Bootstrap: 1,000 resamples of human protein clusters (the Q2.2 paralog units), seed 0, percentile 95% CI.
Every TSV carries the DESIGN.md sha256 and the git HEAD in its first line.

## Question and link to the literature
Are interface (multi-chain, MBS-like) pockets of human predicted dimers kept in the orthologous dimers of other species more often or less often than single-chain (SBS-like) pockets?
Abrusán & Marsh 2018 (LITERATURE.md (a)) found that in PDB **homomers** MBS sites are much more often found in homologs than SBS sites are, especially in distant homologs. MBS **heteromers** did not show this.
Our setting differs in four ways: the structures are AF predicted complexes rather than PDB biounits, the sites are BURROW cavities rather than ligand sites, orthology is by RBH rather than any BLAST homolog, and conservation is residue overlap through a sequence alignment rather than ProBiS.

## Methods (scripts in rerun order)
| step | script | notes |
|---|---|---|
| annotate pockets of 12 species; chain sequences of all 13 species | `scripts/q24_prepare.py` (12 workers) | `annotate_pocket_chains.annotate()` unchanged (interface = ≥3 lining residues per chain). Chain sequences and residue numbers come from the model PDBs. Checks: every chain is contiguous, 12 of ~123k chains do not start at residue 1 (offset handled), and no accession has two different sequences |
| divergence times | `scripts/q24_divergence.py` | TimeTree pairwise API (`timetree.org/api/pairwise/9606/<taxid>`). **The API returned HTTP 500 for every bacterial taxid**, so bacteria are set to an **assumed** 4,200 Mya (approximate LUCA age). This is documented in the script, and only the ordering is used. Ties: mouse = rat 84, fly = worm 692, yeast = pombe 1,164, arabidopsis = dicty 1,516 |
| RBH | `scripts/q24_rbh.sh` | mmseqs2 18.8cc5c, both directions, e ≤1e-3, -s 7.5, no identity floor. RBH = mutual best hit by bitscore |
| scoring | `scripts/q24_conservation.py` | see below |
| summaries | `scripts/q24_stats.py`, `scripts/q24_figures.py` | |

Definitions (declared in the `q24_conservation.py` docstring):
- **Human pockets** follow Q2.2: atlas_pass, symmetry representatives, interface (≥3/≥3) or single-chain (0 on one chain), with marginal pockets excluded. Interface pockets are split into dimer-only and pre-existing by the annotated `dimer_only` flag.
- **Ortholog dimer**: a species model whose accession pair equals {RBH(A), RBH(B)}. For homodimers A = B.
- **Mapping and conservation**: lining residues are carried through the RBH alignment. score = max over ortholog models, chain orientations (both for homodimers) and target pockets T of |M∩T| / min(|L|, |T|). Here M is the mapped set and L the human lining set, so unaligned residues count as lost. A pocket is conserved if score ≥0.5. The primary targets are atlas_pass ortholog pockets; all ortholog pockets are a sensitivity analysis. A third variant requires the target pocket to be an interface pocket.
- **Negative control**: the same mapped residue set (ortholog numbering) is scored on a random species dimer chosen with seed 0. The random dimer is of the same kind, has no mmseqs2 hit to either human protein, and has chain lengths within ±20% of the ortholog's chains (nearest total length if no such dimer exists).
- **Lining identity**: the fraction of aligned lining residues that are identical. It is also reported minus the chain alignment identity (fident).
- **Units**: the unit of observation is (human pocket, species), giving 80,557 rows (28,586 human pockets, 2,341 bootstrap units). Groups pool species with tied divergence times.

## 1. Ortholog dimer pairs (`results/part2/q24/q24_ortholog_pairs.tsv`)
| species | Mya | RBH human proteins | human homo models with ortholog dimer | human hetero models with ortholog dimer |
|---|---|---|---|---|
| mouse | 84 | 4,290 | 1,785 | 2,794 |
| rat | 84 | 3,906 | 1,641 | 2,238 |
| zebrafish | 435 | 2,885 | 1,236 | 1,167 |
| fly | 692 | 1,659 | 541 | 668 |
| worm | 692 | 1,421 | 502 | 498 |
| yeast | 1,164 | 882 | 257 | 306 |
| pombe | 1,164 | 912 | 252 | 366 |
| arabidopsis | 1,516 | 1,359 | 417 | 465 |
| dicty | 1,516 | 1,004 | 358 | 383 |
| E. coli | 4,200* | 290 | 200 | 14 |
| M. tuberculosis | 4,200* | 280 | 184 | 21 |
| P. aeruginosa | 4,200* | 372 | 239 | 30 |

\* assumed value (see Methods). Bacterial comparisons are mostly homodimers, and the bacterial heterodimer counts are small.

## 2. Pocket conservation (`q24_conservation.tsv`, `q24_conservation_strata.tsv`; figures `results/part2/figures/q24_conservation_by_{species,score_q,id_bin}.png`)
Primary target set. The figures give the fraction conserved [95% CI] and the non-ortholog null.

| scope | interface | single | Δ raw [CI] | Δ burrow-quintile matched [CI] | dimer-only | pre-existing | Δ matched [CI] |
|---|---|---|---|---|---|---|---|
| mammals | 0.84 | 0.77 | +0.066 [0.052, 0.079] | +0.037 [0.027, 0.048] | 0.80 | 0.91 | −0.025 [−0.043, −0.008] |
| fish | 0.78 | 0.71 | +0.065 | +0.037 [0.019, 0.054] | 0.71 | 0.89 | −0.081 [−0.115, −0.049] |
| invertebrates | 0.67 | 0.59 | +0.079 | +0.049 [0.028, 0.068] | 0.60 | 0.79 | −0.065 [−0.110, −0.020] |
| fungi | 0.61 | 0.54 | +0.071 | +0.046 [0.015, 0.077] | 0.53 | 0.75 | −0.100 [−0.156, −0.045] |
| plant/amoeba | 0.60 | 0.54 | +0.068 | +0.039 [0.007, 0.070] | 0.52 | 0.76 | −0.120 [−0.183, −0.068] |
| bacteria | 0.50 | 0.48 | +0.017 [−0.042, 0.069] | +0.007 [−0.039, 0.052] | 0.39 | 0.68 | −0.138 [−0.233, −0.045] |
| ALL | 0.75 | 0.69 | +0.059 [0.044, 0.073] | **+0.030 [0.019, 0.041]** | 0.69 | 0.85 | **−0.063 [−0.081, −0.046]** |

- **Negative control.** Mapping the same residue sets onto random non-orthologous dimers of similar length gives conservation of 0.01 to 0.08 in every species and class (ALL: interface 0.032, single 0.030). This is near chance and far below the 0.45 to 0.84 observed. The ortholog signal is therefore not an artefact of residue-set size or pocket density.
- **Divergence.** Conservation falls with divergence for every class. Interface pockets fall from 0.84 (mammals) to 0.50 (bacteria), single-chain pockets from 0.77 to 0.48.
- **Interface vs single-chain.** Interface pockets are conserved slightly *more* often. The gap is about +6 points raw and +3 points after burrow_score-quintile matching, and it is present in every eukaryote group. In bacteria the CI spans 0. The gap is also present in every identity bin (ALL, <30% to ≥90%: +8, +5, +7, +7, +6 points), so it is not simply a matter of interface pockets sitting in more conserved proteins. The sensitivity analysis with all targets gives similar values (ALL matched +0.028 [0.018, 0.039]).
- **Homo vs hetero, opposite to Abrusán & Marsh.** The excess comes from **heterodimers** (ALL_hetero matched +0.047 [0.032, 0.063]). In homodimers it is null (ALL_homo +0.006 [−0.010, 0.021]); fungi_homo −0.043 [−0.089, 0.006] and bacteria_homo −0.003. Abrusán & Marsh reported conservation for MBS homomers and none for MBS heteromers. Possible reasons, all untested: (i) ortholog heterodimer pairs are selected for conserved complexes (both partners must be RBH *and* both must be predicted to bind confidently); (ii) AF may predict an interface by homology/MSA co-evolution, which makes the ortholog models non-independent; (iii) the site definitions differ (cavities, not ligands).
- **Dimer-only vs pre-existing interface pockets.** Dimer-only pockets are conserved clearly *less* often: raw −0.16 [−0.18, −0.15], matched −0.063 [−0.081, −0.046]. The deficit grows with divergence (matched −0.025 in mammals to −0.14 in bacteria). This fits the Q2.2 finding that a sizeable part of the dimer-only calls is method noise: dimer-only pockets are low-scoring, and the noise floor is about 17%. The deficit also persists in the top burrow quintiles in the strata table. Pre-existing interface pockets have a higher null (0.06) because they are larger.
- **Interface-target variant.** If the target pocket must also be an interface pocket, 65% of human interface pockets are conserved as interface pockets, against a null of 0.012. Dimer-only vs pre-existing is then near 0 after matching (+0.010 [−0.009, 0.028]). Most of the dimer-only deficit therefore comes from matches to single-chain target pockets, which pre-existing pockets (concavities already in the monomer) can find more easily.
- **Burrow_score.** Conservation rises steeply with score (ALL Q1 → Q5: interface 0.52 → 0.93, single 0.44 → 0.92). The interface excess is largest in Q1 (+8 points) and about 0.6 points in Q5. So, as `reviews/verify_q22.md` anticipated, a good part of the raw contrast is score-dependent.

## 3. Lining-residue sequence conservation (`q24_lining_identity.tsv`; `q24_lining_identity.png`)
Lining residues are more conserved than their chains in all classes (lining minus chain identity +0.04 to +0.11). Interface pockets show **no** extra sequence conservation over single-chain pockets: ALL lining minus chain −0.003 [−0.008, 0.002]. In distant species interface pockets are actually *lower* (plant/amoeba −0.024 [−0.037, −0.011], bacteria −0.050 [−0.076, −0.023]). Dimer-only lining residues are less conserved than pre-existing ones (ALL −0.019 [−0.026, −0.013]; bacteria −0.069). So the structural excess in §2 does not come with an excess in sequence conservation.

## Caveats
- Conditioning: a pocket can only be scored when an ortholog dimer exists in the species set, which requires a confidently predicted complex (ipSAE ≥0.6, pDockQ2 ≥0.23). This selects conserved complexes. The question answered is "given a conserved predicted complex", not whether the complex itself is conserved.
- The models are not independent across species. AF uses similar MSAs and templates for orthologs, so structural similarity partly reflects the predictor. The PDB is in AF's training data.
- BURROW context noise (Q2.2: about 17% of single-chain pockets lose their counterpart between runs). Conservation is therefore capped below 1 by noise. Rates are compared between classes, with burrow-quintile matching.
- RBH is a simple orthology proxy and misses many-to-many orthology, e.g. after vertebrate paralog expansions. For homodimers, a human pocket may map to a pocket on either chain of the ortholog. The mapping takes the maximum over both orientations, and the real and null scores use the same rule.
- Ledger: 50 models without pockets across species, and **P. aeruginosa homodimer chunk_0003 (~117 models) was never run** (`data/interim/q24/ledger_paeruginosa.tsv`). Bacterial divergence times are assumed, not fetched.
- Counts per (pocket, species) are not independent across species. The bootstrap resamples human protein clusters, so all species of one protein move together.

## Files
Scripts: `scripts/q24_prepare.py`, `q24_divergence.py`, `q24_rbh.sh`, `q24_conservation.py`, `q24_stats.py`, `q24_figures.py`. Results: `results/part2/q24/*.tsv`. Figures:
`results/part2/figures/q24_conservation_by_species.png`, `q24_conservation_by_score_q.png`, `q24_conservation_by_id_bin.png`, `q24_lining_identity.png`.
Intermediates (not committed): `data/interim/q24/` (annotations, chains, fasta, mmseqs hits, `q24_pocket_scores.parquet`). Logs: `logs/q24/`.
