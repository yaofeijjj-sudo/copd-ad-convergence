#!/bin/bash
# Prespecified MAGMA annotation-window sensitivity (protocol 2026-07-15: genic-only primary; ±10 kb; 35 kb up / 10 kb down),
# rerun at the gene-set step only for the current fixed Axis A (137) / Axis B (434) Entrez sets. Gene-level inputs are the frozen
# genes.raw files of the same AD GWAS (GCST90704648, no-biobank meta-analysis; 13,648,198 SNPs). Run from the project root.
set -euo pipefail
MG=_workspace/project_tools/external_tools/magma_v1.10_official/magma_v1.10_mac/magma
SETS="GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma_inputs/locked_axes.sets"
R="GWAS和外周血/analysis/05_disease_specificity_gwas/results/current_axis_full135_full418_multiancestry_gwas_magma_20260715/magma"
OUT=03_candidate_analysis/AxisA_MAGMA_window_sensitivity_20261027/results
for w in genic_only 10kb 35up10down; do
  $MG --gene-results "$R/ad_no_biobank.$w.genes.raw" --set-annot "$SETS" --out "$OUT/AD_$w"
done
