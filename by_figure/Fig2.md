# Fig 2. Disease-associated co-expression modules in lung and brain and their overlap

Drawn by [`fig2.R`](../manuscript/main_figures/fig2.R).

| Panel | Content | Code and plotted values |
| --- | --- | --- |
| a, b | WGCNA networks and module–trait correlations | [`01_build_networks.R`](../scripts/01_wgcna/01_build_networks.R), [`02_module_traits.R`](../scripts/01_wgcna/02_module_traits.R) |
| c, d | Fisher tests of the 128 module pairs; Axis A and Axis B | [`02_module_overlap`](../scripts/02_module_overlap/), [`Fig2c_module_overlap.tsv`](../manuscript/main_figures/source_data/Fig2c_module_overlap.tsv) |
| e | Adjusted case–control effects of the source modules | [`03_lung_adjusted_models.R`](../scripts/01_wgcna/03_lung_adjusted_models.R), [`03_adjusted_effects.R`](../scripts/02_module_overlap/03_adjusted_effects.R) |

[All figures](README.md) · [Repository README](../README.md)
