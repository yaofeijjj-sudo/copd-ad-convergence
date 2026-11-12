# Analyses

One folder per analysis: `code/` and a README with its scripts.

| Folder | Analysis | Used in |
| --- | --- | --- |
| **Replication** | | |
| [`replication_scores_and_meta_analysis`](replication_scores_and_meta_analysis/) | Axis A scores and replication | Fig 3c, Table D |
| [`replication_sample_checks_and_member_effects`](replication_sample_checks_and_member_effects/) | Sample checks and member effects | Fig 3c–d, Table D, Table E |
| [`replication_sensitivity_analyses`](replication_sensitivity_analyses/) | Replication sensitivity analyses | Fig 3c, Table E |
| [`coexpression_sigqc`](coexpression_sigqc/) | Co-expression of Axis A members | Fig 3b, Table C |
| [`coexpression_sigqc_gse48350_without_duplicate`](coexpression_sigqc_gse48350_without_duplicate/) | Co-expression in GSE48350 without the duplicated array | Fig 3b, Table C |
| [`axisB_replication`](axisB_replication/) | Axis B replication | Fig B |
| **Single cell** | | |
| [`single_cell_lung_55_donors`](single_cell_lung_55_donors/) | Lung single-cell analyses with 55 donors | Fig 4–5, Fig 8, Table F, Table L, Table M |
| [`single_cell_lung_library_qc`](single_cell_lung_library_qc/) | Lung library QC and cell-class map | Table F |
| [`single_cell_pseudobulk`](single_cell_pseudobulk/) | Donor pseudobulk | Fig 5c, Table F |
| [`single_cell_pseudobulk_tests_and_tf_activity`](single_cell_pseudobulk_tests_and_tf_activity/) | Pseudobulk tests and TF activity | Fig 5c, Fig 8f, Table F, Table M |
| [`single_cell_alveolar_macrophages`](single_cell_alveolar_macrophages/) | Alveolar macrophages | Fig 5c, Table F |
| **Cell proportion** | | |
| [`cell_proportion_cibersortx`](cell_proportion_cibersortx/) | CIBERSORTx cell proportions | Fig 5d–f, Table G |
| [`cell_proportion_discovery_cohorts`](cell_proportion_discovery_cohorts/) | Microglial proportion in the brain discovery cohort | Table G |
| **Function and genetics** | | |
| [`myeloid_signatures`](myeloid_signatures/) | Myeloid signatures and functional gene sets | Fig 6b, Table I |
| [`magma_apoe_exclusion`](magma_apoe_exclusion/) | MAGMA without the APOE region | Fig 6c, Table J |
| [`magma_conditional_microglia`](magma_conditional_microglia/) | MAGMA conditioned on microglial genes | Fig 6c, Table J |
| **Blood and disease scope** | | |
| [`blood_magma_and_scope_checks`](blood_magma_and_scope_checks/) | Blood and MAGMA follow-up analyses | Fig 6c, Fig 7a, Table J, Table K |
| [`disease_scope`](disease_scope/) | Disease scope | Fig 7b, Table D |
| [`disease_scope_lung_family`](disease_scope_lung_family/) | Lung disease-scope family | Fig 7b, Table D |
| **Transcription factors** | | |
| [`celloracle_ranking`](celloracle_ranking/) | CellOracle ranking | Fig 8b–e, Table L |

## Paths in the scripts

The scripts were run from the project root and keep its directory names:

| Path in the scripts | Folder |
| --- | --- |
| `03_candidate_analysis/<name>/` | the folder listed below |
| `cell_localization/`, `functional_enrichment/`, `GWAS_and_blood/`, `virtual_knockout/`, `WGCNA_figures/` | original runs of the core pipeline and their inputs |
| `05_github_release/copd-ad-convergence/` | this repository |
| `06_resubmission_20261008/code/` | [`manuscript/`](../manuscript/) |

| `03_candidate_analysis/` name | Folder |
| --- | --- |
| `AxisA_singscore_official_rerun_20261004` | [`replication_scores_and_meta_analysis`](replication_scores_and_meta_analysis/) |
| `AxisA_member_direction_robustness_20261018` | [`replication_sample_checks_and_member_effects`](replication_sample_checks_and_member_effects/) |
| `AxisA_replication_after_exclusion_20261019` | [`replication_sensitivity_analyses`](replication_sensitivity_analyses/) |
| `AxisA_sigQC_coherence_20261014` | [`coexpression_sigqc`](coexpression_sigqc/) |
| `AxisA_sigQC_GSE48350_dedup_20261030` | [`coexpression_sigqc_gse48350_without_duplicate`](coexpression_sigqc_gse48350_without_duplicate/) |
| `AxisB_fixed_score_replication_20261019` | [`axisB_replication`](axisB_replication/) |
| `Lung_sc_donor_merge_20261029` | [`single_cell_lung_55_donors`](single_cell_lung_55_donors/) |
| `Lung_sc_library_QC_check_20261028` | [`single_cell_lung_library_qc`](single_cell_lung_library_qc/) |
| `AxisA_pseudobulk_disease_20261014` | [`single_cell_pseudobulk`](single_cell_pseudobulk/) |
| `SingleCell_pseudobulk_supplements_20261019` | [`single_cell_pseudobulk_tests_and_tf_activity`](single_cell_pseudobulk_tests_and_tf_activity/) |
| `AxisA_AM_pseudobulk_20261028` | [`single_cell_alveolar_macrophages`](single_cell_alveolar_macrophages/) |
| `CIBERSORTx_composition_adjustment_20261018` | [`cell_proportion_cibersortx`](cell_proportion_cibersortx/) |
| `Deconvolution_composition_adjustment_20261002` | [`cell_proportion_discovery_cohorts`](cell_proportion_discovery_cohorts/) |
| `AxisA_myeloid_reference_20261019` | [`myeloid_signatures`](myeloid_signatures/) |
| `AxisA_MAGMA_APOE_exclusion_20261016` | [`magma_apoe_exclusion`](magma_apoe_exclusion/) |
| `AxisA_MAGMA_conditional_microglia_20261016` | [`magma_conditional_microglia`](magma_conditional_microglia/) |
| `AxisA_disease_specificity_redesign_20261003` | [`disease_scope`](disease_scope/) |
| `AxisA_disease_scope_lung_family_20261028` | [`disease_scope_lung_family`](disease_scope_lung_family/) |
| `AxisA_auxiliary_checks_20261019` | [`blood_magma_and_scope_checks`](blood_magma_and_scope_checks/) |
| `CellOracle_ranking_scale_20261019` | [`celloracle_ranking`](celloracle_ranking/) |
