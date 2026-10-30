# Ligand-matched vs unmatched predicted interface pockets (EXPLORATORY)

**Status: UNVERIFIED, EXPLORATORY.** Not pre-registered; first run, not re-derived independently, figure not reviewed by the PI.

Question: the predicted interface pockets in AF dimers are more charged and less hydrophobic than single-chain pockets (docs/Q22_RESULTS.md), but real ligand-bound multi-chain PDB sites are not (docs/PDB_SITE_DESCRIPTORS.md). Hypothesis: predicted interface pockets that match real ligand sites look like ordinary ligand pockets, and the excess charged ones are ligand-free PPI grooves.

## Rerun
```bash
source scripts/env_threads.sh
.venvs/part1/bin/python scripts/iface_lm_match.py     # ~1 min, 4 workers -> data/interim/iface_ligand_matched/pocket_match.parquet
.venvs/q22/bin/python  scripts/iface_lm_analyze.py    # ~2 min -> results/part2/iface_ligand_matched/*.tsv
.venvs/q22/bin/python  scripts/iface_lm_strat.py      # burrow_score-stratified IM-IU difference
.venvs/q22/bin/python  scripts/iface_lm_figure.py     # -> iface_ligand_matched.png
```
Lining definitions (charged DEKR, hydrophobic AVLIMFWC, aromatic FWYH, buriedness) are taken unchanged from `data/interim/q22/dimer_pocket_features_units.parquet` (made by `q22_pocket_features.py`); `cliffs_delta` and `cluster_boot` are imported from `q22_common.py`.

## Method
- **Universe.** All dimer pockets of the 870 models with H4 v2 site rows (`data/part2/h4/v2/recovery.parquet`), restricted to atlas-pass, `is_sym_rep` (annotation rule, as in q22_compare.py), excluding "marginal" pockets (1-2 residues on one chain). Interface = >=3 residues on each chain; single-chain = 0 residues on one chain.
- **Match** (`iface_lm_match.py`). Pocket matches a site if DCC <=4 Å (pocket centre vs ligand centroid superposed with the H4 local-frame rule, homodimer chain swap tried) OR overlap coefficient |L∩P|/min(|L|,|P|) >=0.5 with the UniProt-mapped PDB binding residues. All pockets of the model are tested, not only the top-(n+2) as in the H4 hit rule.
- **Groups.** IM: interface pocket matched to a primary-set (organic, interface-filtered) MBS site. IU: interface pocket matched to no site in any recovery row of that model (any ligand class, MBS or SBS, no interface filter; up to 3 entries per unit, so "any entry of the pair" is only approximate). SM: single-chain pocket matched to a primary-set SBS site. SU: single-chain, no match (reference). Interface pockets matched only to SBS or non-primary sites fall in neither IM nor IU ("other", 469 pockets incl. single-chain ones).
- **Statistics.** Paralog-cluster bootstrap (1,000 resamples, seed 0), medians and Cliff's δ (`results/part2/iface_ligand_matched/`). Classifier: logistic regression (L2=1, own numpy implementation; sklearn is not installed), 5-fold CV with folds split by paralog cluster, AUC CI by cluster bootstrap of out-of-fold scores.

## Results
n: IM 173 (126 clusters), IU 1,193 (501), SM 516 (292), SU 2,181 (465). Of 1,547 interface pockets in these models, 173 are IM.

Medians [95% CI]:

| | IM | IU | SM |
|---|---|---|---|
| charged | 0.204 [0.182, 0.224] | 0.208 [0.200, 0.222] | 0.174 [0.158, 0.198] |
| hydrophobic | 0.367 [0.342, 0.393] | 0.333 [0.333, 0.357] | 0.380 [0.359, 0.400] |
| aromatic | 0.159 | 0.182 | 0.158 |
| volume Å³ | 342 | 92 | 109 |
| burrow_score | 0.44 | 0.08 | 0.30 |
| buriedness | 0.791 | 0.788 | 0.802 |

Cliff's δ [CI]:
- **IM vs IU:** charged -0.05 [-0.13, +0.04]; hydrophobic +0.11 [+0.02, +0.19]; aromatic -0.05 [-0.13, +0.04]; volume +0.56 [+0.47, +0.64]; burrow_score +0.75 [+0.68, +0.81]; buriedness +0.04 [-0.04, +0.12].
- **IM vs SM:** charged +0.10 [-0.02, +0.21]; hydrophobic -0.01 [-0.13, +0.11]; aromatic +0.07 [-0.03, +0.16]; volume +0.43; burrow_score +0.27; buriedness -0.05 [-0.15, +0.04].
- **IU vs SU (reference, Q2.2-like):** charged +0.22 [+0.16, +0.27]; hydrophobic -0.21 [-0.26, -0.16].
- DCC-only matching (IM n=72, `cliffs_delta_dcc_only_sensitivity.tsv`): charged -0.11 [-0.24, +0.01]; hydrophobic +0.14 [+0.02, +0.26]; same picture.
- Stratified on burrow_score (`burrow_stratified_diff.tsv`, IM quintile strata): mean charged IM-IU = -0.026 [-0.067, +0.012]; hydrophobic +0.057 [+0.011, +0.105].

Classifier, IM vs IU, cluster-split CV (`classifier_auc.tsv`): charged fraction alone AUC 0.51 [0.46, 0.55]; burrow_score alone 0.87 [0.84, 0.90]; charged + burrow_score 0.87 [0.84, 0.90]; volume alone 0.77.

## Reading
- The hypothesis is **not supported as stated**. Ligand-matched interface pockets are not less charged than unmatched ones (δ -0.05, CI includes 0), and the charged fraction carries no information about matching (AUC 0.51). The unmatched interface pockets keep the full charged excess over unmatched single-chain pockets (IU vs SU +0.22), but so do the matched ones within the resolution of n=173: IM vs SM is +0.10 with a CI that includes 0, which is compatible with either no excess or an excess of about half the Q2.2 size.
- What separates matched from unmatched pockets is size and score (burrow_score δ +0.75, volume +0.56), i.e. matching mostly selects large, well-scored pockets. The small hydrophobic shift (+0.11, +0.06 after score stratification) is the only chemistry signal and points in the hypothesised direction.
- So the PDB-versus-AF discrepancy is not explained by matched pockets being chemically ordinary and unmatched ones being charged. Ligand occupancy in the PDB seems to select pocket size/score rather than a charge-poor subset; the residual hydrophobic shift is weak.

## Caveats
- IM is small (173 pockets, 126 clusters, DCC-only 72); CIs are wide. The 0.5 overlap-coefficient rule is lenient for small ligand sites against large pockets (435 of 440 raw MBS-pocket hits come via overlap, 186 via DCC), and it favours large pockets.
- "Unmatched" is bounded by PDB ligand coverage: an unmatched pocket may bind an uncrystallised ligand, so IU is not guaranteed ligand-free.
- Benchmark selection (BioLiP2 sites in entries with human dimer models, ipSAE >=0.6, interface Jaccard >=0.3 for primary sites) is not representative of all atlas interface pockets.
- Pocket-level units are not independent within a model; the bootstrap is by paralog cluster only.
- Matching used all pockets, not the H4 top-k rule, so counts differ from H4 hit rates. Superposition parameters are those of H4 but were re-implemented here.
