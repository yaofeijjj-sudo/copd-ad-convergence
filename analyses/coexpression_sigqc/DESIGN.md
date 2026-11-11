# Axis A coherence in independent datasets (sigQC) — design, frozen 2026-10-14 before any run

## Question

Axis A is a fixed set of 137 Entrez genes (lung red ∩ brain blue). Existing evidence shows that the two *source modules* are preserved (Zsummary) and that the Axis A *score* is associated with disease in independent cohorts. What is missing is whether **these 137 genes themselves still co-vary as one set** in data that were not used to define them, and which members do so consistently.

## Method and source

sigQC 0.1.24 (Dhawan et al. 2019 *Nat Protoc* 14:1377, doi:10.1038/s41596-019-0136-8), the official CRAN release (archived on CRAN on 2025-12-19 because its dependency biclust was archived; installed from the CRAN archive into `_workspace/r_libs_sigqc`). sigQC evaluates a previously derived gene signature in new expression datasets: expression, variability, inter-gene correlation (compactness), agreement among scoring methods (mean, median, PCA1) and standardisation, with negative controls from random gene sets of equal size and permutation controls.

## Estimand and unit

- Unit: the dataset. One sigQC evaluation per independent dataset; no pooling of samples across datasets.
- Primary quantity: sigQC compactness `autocor_median` = median pairwise Spearman correlation among the Axis A members measured in that dataset.
- Primary comparison: the observed `autocor_median` against an **expression-matched random-gene null** (1,000 draws; genes sampled from the same 20 mean-expression bins as the measured Axis A members, the scheme already used in the singscore branch), reported as the null 95th percentile and a one-sided empirical P = (1 + #null ≥ observed)/1001. sigQC's own negative control (random genes of equal size, not expression-matched; `numResampling = 200`) and permutation control are reported alongside.
- Secondary (descriptive): the other sigQC radar metrics per dataset.
- Gene level (descriptive): from sigQC's autocorrelation matrix, each member's median Spearman correlation with the other members in each dataset.

## Frozen inputs

- Gene set: `03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv` (137).
- Datasets and samples exactly as in the independent replication analyses (`AxisA_singscore_official_rerun_20261004`), with identical preprocessing (probe → Entrez collapse by maximum MAD; RNA-seq: `filterByExpr` + TMM + log-CPM, prior.count 2; GSE118553 and GSE47460 averaged per donor):
  1. GSE47460-GPL14550 lung, COPD + control donors (primary lung)
  2. GSE47460-GPL6480 lung, COPD + control donors
  3. GSE84422 DLPFC, AD + normal
  4. GSE118553 frontal cortex, AD + control (per individual)
  5. GSE313307 prefrontal cortex RNA-seq, AD + control
  6. GSE140841 BA9 RNA-seq, AD + control
  7. GSE48350 superior frontal gyrus (RMA from CEL), AD + age-restricted controls
  8. GSE36980 frontal cortex, AD + non-AD
- Disease labels are not used by sigQC or by any statistic here. (RNA-seq low-count filtering uses the group argument of `filterByExpr`, kept only to reproduce the replication matrices.)
- Discovery cohorts (GSE76925, GSE33000) are not evaluated.

## Pre-written interpretation rules

- Set level: if observed `autocor_median` exceeds the 95th percentile of the expression-matched null in a dataset, Axis A is called "coherent" in that dataset. Report the count over 8 datasets and per tissue; no dataset is dropped after the fact.
- Gene level: a member is labelled "consistent" if its median correlation with the other members is > 0 in every dataset where it is measured and it is measured in at least 6 of the 8 datasets; otherwise "inconsistent" or "insufficient coverage". These labels are descriptive; Axis A is not redefined and all scores keep using 137 genes.
- Interpretation limit: in bulk tissue, inter-gene correlation partly reflects variation in cell-type composition (here, myeloid content). Coherence therefore supports "a co-varying tissue myeloid program", not per-cell co-regulation.

## Seeds and software

Seed 20261014 for the expression-matched null and sigQC resampling. R 4.5.3 (`_workspace/envs/r`), library order `_workspace/r_libs_sigqc`, `_workspace/r_libs_affy`, `_workspace/r_libs`.
