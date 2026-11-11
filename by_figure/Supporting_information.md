# Supporting information

All parts are built by [`run_all.sh`](../manuscript/supplement/run_all.sh).

## S1 File

| Figure | Content | Code |
| --- | --- | --- |
| Fig A | Soft-threshold selection | [`01_build_networks.R`](../scripts/01_wgcna/01_build_networks.R), [`build_S1_S2_File.py`](../manuscript/supplement/build_S1_S2_File.py) |
| Fig B | Axis B enrichment and replication | [`03_enrichment`](../scripts/03_enrichment/), [`axisB_replication`](../analyses/axisB_replication/), [`build_S1_S2_Fig.R`](../manuscript/supplement/build_S1_S2_Fig.R) |

## S2 File

| Table | Content | Code |
| --- | --- | --- |
| Table A | 128 module pairs; grey-module background; restored age outliers | [`02_module_overlap`](../scripts/02_module_overlap/), [`method_checks`](../analyses/method_checks/), [`wgcna_brain_age_outlier_sensitivity`](../analyses/wgcna_brain_age_outlier_sensitivity/) |
| Table B | Axis A and Axis B gene sets | [`AxisA_confirmed_genes.csv`](../resources/AxisA_confirmed_genes.csv), [`AxisB_confirmed_genes.csv`](../resources/AxisB_confirmed_genes.csv) |
| Table C | sigQC co-expression | [`coexpression_sigqc`](../analyses/coexpression_sigqc/), [`coexpression_sigqc_gse48350_without_duplicate`](../analyses/coexpression_sigqc_gse48350_without_duplicate/) |
| Table D | Datasets, sample selection, donor checks, disease scope | [`replication_scores_and_meta_analysis`](../analyses/replication_scores_and_meta_analysis/), [`replication_sample_checks_and_member_effects`](../analyses/replication_sample_checks_and_member_effects/), [`disease_scope`](../analyses/disease_scope/), [`disease_scope_lung_family`](../analyses/disease_scope_lung_family/) |
| Table E | Replication sensitivity analyses and member effects | [`replication_sensitivity_analyses`](../analyses/replication_sensitivity_analyses/), [`replication_hartung_knapp`](../analyses/replication_hartung_knapp/), [`replication_sample_checks_and_member_effects`](../analyses/replication_sample_checks_and_member_effects/) |
| Table F | Single-cell donor-level tests, QC and annotation | [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/), [`single_cell_pseudobulk`](../analyses/single_cell_pseudobulk/), [`single_cell_pseudobulk_tests_and_tf_activity`](../analyses/single_cell_pseudobulk_tests_and_tf_activity/), [`single_cell_alveolar_macrophages`](../analyses/single_cell_alveolar_macrophages/), [`single_cell_lung_library_qc`](../analyses/single_cell_lung_library_qc/), [`method_checks`](../analyses/method_checks/) |
| Table G | Cell proportion adjustment | [`cell_proportion_cibersortx`](../analyses/cell_proportion_cibersortx/), [`cell_proportion_discovery_cohorts`](../analyses/cell_proportion_discovery_cohorts/) |
| Table H | Functional enrichment | [`03_enrichment`](../scripts/03_enrichment/) |
| Table I | Myeloid signatures and functional gene sets | [`myeloid_signatures`](../analyses/myeloid_signatures/) |
| Table J | MAGMA sensitivity analyses | [`magma_apoe_exclusion`](../analyses/magma_apoe_exclusion/), [`magma_conditional_microglia`](../analyses/magma_conditional_microglia/), [`blood_magma_and_scope_checks`](../analyses/blood_magma_and_scope_checks/), [`magma_gene_windows`](../analyses/magma_gene_windows/) |
| Table K | Peripheral blood | [`04_genetics_blood`](../scripts/04_genetics_blood/), [`blood_magma_and_scope_checks`](../analyses/blood_magma_and_scope_checks/) |
| Table L | Virtual knockout ranking of 357 TFs | [`06_virtual_knockout`](../scripts/06_virtual_knockout/), [`celloracle_ranking`](../analyses/celloracle_ranking/), [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/) |
| Table M | TF activity in donor pseudobulk | [`single_cell_pseudobulk_tests_and_tf_activity`](../analyses/single_cell_pseudobulk_tests_and_tf_activity/), [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/) |

## S1 Data

Values underlying each figure panel: [`build_S1_Data.py`](../manuscript/supplement/build_S1_Data.py); the sheet-to-file map is [`S1_Data_source_map.tsv`](../manuscript/supplement/parts/S1_Data_source_map.tsv).

[All figures](README.md) · [Repository README](../README.md)
