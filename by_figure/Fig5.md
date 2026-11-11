# Fig 5. Within-cell tests and the relation of Axis A to the target cell proportion

Drawn by [`fig5.R`](../manuscript/main_figures/fig5.R), [`00_composition_samples.R`](../manuscript/main_figures/00_composition_samples.R).

| Panel | Content | Code and plotted values |
| --- | --- | --- |
| a | Donor UCell scores in lung myeloid cells (55 donors) | [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/), [`Fig5a_lung_donors.tsv`](../manuscript/main_figures/source_data/Fig5a_lung_donors.tsv) |
| b | Brain Microglia 2 against the local CPS | [`07_disease_association.py`](../scripts/05_single_cell/07_disease_association.py), [`cps_model.py`](../R/cps_model.py), [`draw_disease.R`](../manuscript/figure_sources/cell_localization/draw_disease.R), [`Fig5b_brain_donors.tsv`](../manuscript/main_figures/source_data/Fig5b_brain_donors.tsv) |
| c | Pseudobulk fry and CAMERA | [`single_cell_pseudobulk`](../analyses/single_cell_pseudobulk/), [`single_cell_pseudobulk_tests_and_tf_activity`](../analyses/single_cell_pseudobulk_tests_and_tf_activity/), [`single_cell_alveolar_macrophages`](../analyses/single_cell_alveolar_macrophages/), [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/), [`Fig5c_pseudobulk.tsv`](../manuscript/main_figures/source_data/Fig5c_pseudobulk.tsv) |
| d | Disease effects before and after adding the target cell proportion (CIBERSORTx) | [`cell_proportion_cibersortx`](../analyses/cell_proportion_cibersortx/), [`Fig5d_composition.tsv`](../manuscript/main_figures/source_data/Fig5d_composition.tsv) |
| e, f | Scores against proportions; case–control differences in proportion | [`00_composition_samples.R`](../manuscript/main_figures/00_composition_samples.R), [`cell_proportion_cibersortx`](../analyses/cell_proportion_cibersortx/), [`Fig5e_samples.tsv`](../manuscript/main_figures/source_data/Fig5e_samples.tsv), [`Fig5f_proportion_difference.tsv`](../manuscript/main_figures/source_data/Fig5f_proportion_difference.tsv) |

[All figures](README.md) · [Repository README](../README.md)
