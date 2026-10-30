# Independent verification of H5 (disease variants at dimer-only pockets)

Verifier: independent agent, 2026-10-29. Code: `reviews/verify_h5/` (`rebuild_sample.py`, `compare_sample.py`, `seqcheck.py`, `felogit.py`, `fits.py`, `fits_extra.py`, `fit_subset.py`). Outputs: `reviews/verify_h5/out/`, logs: `reviews/verify_h5/*.log`.
None of `scripts/h5_*.py` was imported. The producer's table `data/part2/h5/tables/protein_residues_primary.parquet` was used only for comparison and for the full-table refits (task 2, 4, 5).
Rerun: `source scripts/env_threads.sh; .venvs/h5/bin/python reviews/verify_h5/rebuild_sample.py --workers 6; .venvs/part1/bin/python reviews/verify_h5/{seqcheck,compare_sample,fit_subset,fits,fits_extra}.py`.

## 1. Residue-category rebuild (300 proteins, seed 0; 679 dimer models; 127,409 protein×position residues)
Built from the raw files: the PDB models (gemmi heavy-atom neighbour search ≤5 Å; freesasa relativeTotal; pLDDT from the CA B-factor), the raw burrow describe chunks (`data/part2/pockets/human`, `human_isolated`), and my own atlas, interface-pocket, overlap and dimer-only logic. The sample is drawn with `numpy.default_rng(0).choice` from the sorted union of accessions in candidates.csv (8,279 proteins).

| check | result | verdict |
|---|---|---|
| residue universe (acc×pos) | identical: 127,409 both, 0 left/right only | VERIFIED |
| non-pocket categories (interface / surface / buried) | 100% agreement whenever neither side has a pocket label | VERIFIED |
| pocket flags per pocket (atlas_pass, is_interface ≥3/chain, max_overlap_iso) | 5,482/5,482 identical | VERIFIED |
| dimer-only residues | every producer dimer-only residue is dimer-only in the rebuild (4,151/4,151) | VERIFIED |
| overall agreement, homodimer collapse by chain-swap **Jaccard** ≥0.5 (mine) | **98.4%** | minor DISCREPANCY (definitional) |
| overall agreement, no symmetry collapse | 97.9% | — |

Confusion matrix (rows = rebuild, cols = producer, Jaccard-collapse variant; `out/confusion_collapse.tsv`):

| ver \ prod | buried | dimer_only | interface | sc_pocket | surface |
|---|---|---|---|---|---|
| buried | 28,996 | 0 | 0 | 0 | 0 |
| dimer_only | 32 | 4,151 | 112 | 140 | 53 |
| interface | 0 | 0 | 14,191 | 0 | 0 |
| sc_pocket | 695 | 0 | 374 | 14,818 | 626 |
| surface | 0 | 0 | 0 | 0 | 63,221 |

The cause of every disagreement is the homodimer symmetry collapse. The producer groups pockets by chain-swap **overlap coefficient** ≥0.5 (`annotate_pocket_chains.py`), which is looser and transitive. It drops 619 of 1,593 homodimer pockets in the sample, against 356 under my rule. Lining residues that appear only in a non-representative copy then lose their pocket label. At residue level the collapse is unnecessary, because the per-protein×position aggregation already deduplicates the chains. The producer's choice is conservative (it under-labels pockets) and does not change the subset estimate (below). rSASA and pLDDT covariates differ by >0.01 or >0.5 in only 0.2% of residues; these are the same residues, averaged over a different category set.

## 2. Primary contrast refit with a different method
Method: a protein-stratified fixed-effects logistic regression (Newton with a Schur complement over 1,025 informative proteins, i.e. the proteins with both outcomes, as in clogit). With strata of several hundred residues this approximates conditional logistic regression. Exact statsmodels `ConditionalLogit` hit its recursion limit even on the subset, so it was not used. Covariates: z(rsa_iso) and z(pLDDT), buried residues excluded, interface as reference.

| analysis | dimer-only vs interface OR [95% CI], p | vs producer 1.20 (1.07–1.35), p 0.0024 |
|---|---|---|
| FE-logit, full producer table | **1.19 [1.06, 1.34], 0.0036** | VERIFIED |
| same model, protein-cluster bootstrap (200 resamples, seed 0) | 1.21 [1.07, 1.36] | VERIFIED |
| MH stratified by protein (unadjusted) | 1.39 [1.23, 1.58] | VERIFIED (producer 1.39 [1.23, 1.58]) |
| MH stratified by protein × rsa_iso quintile × pLDDT bin | 1.12 [0.97, 1.30], 0.11 | attenuated |
| rebuilt 300-protein subset, own categories + own humsavar y | 0.91 [0.48, 1.73] (37 informative proteins, 19 dimer-only P/LP) | underpowered; matches the producer table on the same proteins: 0.89 [0.46, 1.72] |
| single-chain pocket vs interface (full) | 1.16 [1.06, 1.27] | VERIFIED (producer 1.17) |
| surface vs interface (full) | 0.53 [0.48, 0.57] | VERIFIED (producer 0.53) |

Category counts were re-derived from the table: residues 109,511 / 401,695 / 349,916 / 1,722,273 and P/LP 547 / 2,061 / 954 / 1,821. VERIFIED. The 300-protein subset cannot test the effect, since it has 10× fewer events. It does show that the end-to-end rebuild (categories plus outcome) gives the same estimate as the producer table on the same proteins.

## 3. humsavar mapping (re-parsed from raw `data/raw/humsavar/humsavar.txt`; `out/seqcheck.tsv`)
| check | value | verdict |
|---|---|---|
| raw missense rows / LP/P rows | 85,014 / 33,332 | VERIFIED (= parquet) |
| P/LP rows in the 8,279 model proteins | 10,519 | VERIFIED |
| P/LP positions absent from the model | 0 | — |
| wt ≠ AFDB model residue (P/LP) | **1 / 10,519 (0.0095%)** | VERIFIED (producer: 1, dropped) |
| wt mismatch, all humsavar classes | 12 / 21,006 (0.06%) | — |
| unique P/LP residues | 8,925 → 8,924 after the mismatch is dropped | VERIFIED (= producer Σy) |
| residue-identity conflicts between models of the same protein | 0 | — |
| y agreement with the producer on the sample | 100% (323 = 323) | VERIFIED |
The models use canonical UniProt numbering, and there is no isoform problem.

## 4. Confounders (full producer table, FE-logit, dimer-only vs interface; `out/fits_full.tsv`, `out/fits_extra.tsv`)
| stress test | OR [95% CI], p |
|---|---|
| no covariates | 1.40 [1.25, 1.57] |
| + rsa_iso, pLDDT (linear; primary-equivalent) | 1.19 [1.06, 1.34], 0.004 |
| rsa_iso and pLDDT as deciles | 1.21 [1.08, 1.36], 0.001 |
| **+ amino-acid identity** | **1.09 [0.97, 1.23], 0.14** |
| deciles + amino acid | 1.12 [0.99, 1.26], 0.06 |
| leave out the top-20 proteins by P/LP (P00439 PAH, P06280 GLA, …) | 1.22 [1.07, 1.38], 0.002 |
| leave out top-20 + amino acid | 1.13 [0.99, 1.29], 0.06 |
| residues with pLDDT ≥70 only | 1.17 [1.04, 1.32], 0.010 |
| 70 ≤ pLDDT < 90 | 1.61 [1.28, 2.04] |
| **pLDDT ≥90 only** | **0.99 [0.86, 1.15], 0.91** (with amino acid: 0.88 [0.75, 1.02]) |
| drop Cys / drop Cys+Gly+Pro / drop aromatic+His | 1.20 / 1.22 / 1.23 (all p<0.003) |

- **Gene burden: not the explanation.** The 547 dimer-only P/LP residues are spread over 254 proteins (top 5 = 9%). Leaving out the top-20 genes does not move the OR.
- **Residue type: explains about half of the log-OR.** Compared with interface residues, dimer-only residues are enriched in Trp/His/Cys/Tyr/Phe/Arg/Gly (ratio 1.09–1.73) and depleted in Lys/Glu/Ala/Leu (0.76–0.88). Arg, Cys, Gly, Tyr and Trp also carry the highest P/LP rates. No single residue class drives the effect: dropping Cys, Gly/Pro or aromatic residues leaves OR ≈1.2. The joint amino-acid adjustment, however, removes significance (1.09, p 0.14).
- **pLDDT: an effect modifier, not a simple linear confounder.** Decile adjustment leaves the OR unchanged. The enrichment is concentrated in moderate-confidence residues (70–90: 1.61) and is absent at pLDDT ≥90 (0.99). Mean pLDDT is 91.1 for dimer-only and 84.8 for interface residues. Interface residues with pLDDT <70 can never be pocket residues (lining filter), so the reference is not exactly matched. Restricting to pLDDT ≥70 still gives 1.17, though.
- MH jointly stratified by protein × rsa × pLDDT: 1.12 (p 0.11).

## 5. Is the effect specific to dimer-only pockets? (formal tests, same FE-logit)
| contrast | OR [95% CI], p |
|---|---|
| dimer-only vs single-chain pocket | 1.03 [0.92, 1.14] (reported as sc vs DO 0.97 [0.87, 1.08]), p 0.62 |
| same + amino acid | 1.01, p 0.87 |
| dimer-only vs pre-existing interface pocket (`cat_split`) | 0.95 [0.84, 1.08], p 0.46 |
| pre-existing interface pocket vs interface | 1.25 [1.12, 1.41] (producer 1.26 [1.12, 1.42]) VERIFIED |
| MH dimer-only vs single-chain pocket, protein × rsa × pLDDT strata | 1.01 [0.89, 1.15] |
**The producer's conclusion is confirmed.** Dimer-only pockets are not more enriched than single-chain pockets or pre-existing interface pockets; all CIs include 1 with point estimates near 1. The signal is a "pocket-lining vs non-pocket interface residue" effect.

## Verdict per headline number
- Primary OR 1.20 (1.07–1.35), p 0.0024: **VERIFIED** (FE-logit 1.19 [1.06, 1.34]; bootstrap 1.21 [1.07, 1.36]).
- Category counts, P/LP counts, humsavar mapping (1 mismatch): **VERIFIED**.
- MH ORs 1.39 / 1.58 / 0.90, sc-pocket OR 1.17, pre-existing interface pocket OR 1.26: **VERIFIED**.
- Category assignment: **VERIFIED with a minor definitional discrepancy** (symmetry-collapse rule, 1.6% of residues, conservative direction).
- Robustness: **FRAGILE**. Not robust to amino-acid identity (1.09, p 0.14), to protein × rsa × pLDDT stratification (1.12, p 0.11), or to restricting to pLDDT ≥90 (0.99).

## Recommendation for reporting H5
1. Report the pre-registered result as passing (OR 1.20, 1.07–1.35, BH-significant) with the GLMM as specified.
2. In the same paragraph, state that the effect is (a) not specific to dimer-only pockets (dimer-only vs single-chain pocket OR ≈1.0 and vs pre-existing interface pocket 0.95), and (b) attenuated to non-significance by amino-acid composition (1.09 [0.97, 1.23]). It is also absent among pLDDT ≥90 residues. Do not phrase the result as "dimer-only pockets are disease hotspots". A defensible phrasing is "pocket-lining residues at predicted interfaces carry modestly more P/LP variants than non-pocket interface residues, partly reflecting their residue composition".
3. Add the amino-acid-adjusted model and the pLDDT ≥90 restriction as declared post-hoc sensitivities. Note the missing conservation covariate, which likely acts in the same direction as amino-acid identity.
4. Optional: drop the homodimer symmetry collapse for residue-level labelling, or document the overlap-coefficient rule. The impact is small.
