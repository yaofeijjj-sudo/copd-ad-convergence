# Main figures 1–8

Each script draws one figure from the result tables and writes its plotted values to [`source_data/`](source_data/).

| Script | Figure | Panels |
| --- | --- | --- |
| [`fig1.R`](fig1.R), [`fig1_icons.R`](fig1_icons.R) | Fig 1 | Study overview; thumbnails reuse the plotted values of Fig 3, 5 and 7 |
| [`fig2.R`](fig2.R) | Fig 2 | a–b module–trait heatmaps; c 128 module pairs; d Axis A/B; e adjusted source-module effects |
| [`fig3.R`](fig3.R) | Fig 3 | a preservation; b sigQC; c replication forest plot; d member effects |
| [`fig4.R`](fig4.R) | Fig 4 | a–d UMAP and UCell; e–f donor-mean scores; g–h matched random gene sets |
| [`fig5.R`](fig5.R) | Fig 5 | a–b donor scores; c pseudobulk fry and CAMERA; d–f cell proportion adjustment |
| [`fig6.R`](fig6.R) | Fig 6 | a enrichment; b myeloid signatures and functional gene sets; c MAGMA |
| [`fig7.R`](fig7.R) | Fig 7 | a peripheral blood; b disease scope |
| [`fig8.R`](fig8.R) | Fig 8 | a workflow; b–e CellOracle virtual knockout; f TF activity |
| [`00_member_effects.R`](00_member_effects.R) | Fig 3d | Exports member-level effects of the adopted sample set (asserts 110 / 5 / 22 / 34) |
| [`00_composition_samples.R`](00_composition_samples.R) | Fig 5e–f | Exports sample scores and proportions and the proportion differences with 95% CI |

[`common.R`](common.R): shared theme and layout. [`tools/`](tools/) and [`qa/`](qa/): PDF text post-processing and checks. Panel-by-panel links: [`by_figure/`](../../by_figure/README.md).

## Run

From the project root: run `00_member_effects.R` and `00_composition_samples.R`, then `fig1.R` … `fig8.R`.
