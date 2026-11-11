# Fig 4. Cellular localization of Axis A

Drawn by [`fig4.R`](../manuscript/main_figures/fig4.R).

| Panel | Content | Code and plotted values |
| --- | --- | --- |
| a–d | Lung QC, annotation and UMAP; brain UMAP; UCell scores | [`preprocess`](../scripts/05_single_cell/preprocess/), [`02_lung_ucell.R`](../scripts/05_single_cell/02_lung_ucell.R), [`03_brain_broad_ucell.R`](../scripts/05_single_cell/03_brain_broad_ucell.R), [`10_lung_umap.py`](../scripts/05_single_cell/10_lung_umap.py), [`draw_umap.R`](../manuscript/figure_sources/cell_localization/draw_umap.R) |
| e, f | Donor-mean scores in lung myeloid subtypes and brain immune supertypes | [`04_brain_subtype_ucell.R`](../scripts/05_single_cell/04_brain_subtype_ucell.R), [`06_lung_myeloid.R`](../scripts/05_single_cell/06_lung_myeloid.R), [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/), [`draw_subtype_raincloud_style_restored.R`](../manuscript/figure_sources/cell_localization/draw_subtype_raincloud_style_restored.R) |
| g, h | Within-donor difference against 5,000 matched random gene sets | [`05_localization_and_null.R`](../scripts/05_single_cell/05_localization_and_null.R), [`single_cell_lung_55_donors`](../analyses/single_cell_lung_55_donors/) |

[All figures](README.md) · [Repository README](../README.md)
