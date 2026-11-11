# Main figures 1–8

The scripts draw each figure from the result tables and write the plotted values to [`source_data/`](source_data/). They do not recompute statistics, except for the two export steps below, which reuse the code of their analyses and assert the reported numbers.

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

[`common.R`](common.R) holds the shared theme (Arial, 8 pt text, 12 pt panel labels), the table-style forest plot, millimetre-based panel placement and the PDF device. PDFs are written with the Quartz device; [`tools/merge_tj_text_runs.py`](tools/merge_tj_text_runs.py) then merges sub-1/1000 em kerning so that each phrase is one text object, and [`qa/textrun_audit.py`](qa/textrun_audit.py) checks fonts, sizes and text objects. Submitted TIFF files were rendered from the PDFs at 600 dpi (RGB, LZW).

Panel-by-panel links to the analysis code are in [`by_figure/`](../../by_figure/README.md).

## Run

From the project root, with the R environment of the study:

```sh
Rscript 06_resubmission_20261008/code/main_figures/00_member_effects.R
Rscript 06_resubmission_20261008/code/main_figures/00_composition_samples.R
for f in fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8; do
  Rscript 06_resubmission_20261008/code/main_figures/$f.R
done
```

The path `06_resubmission_20261008/code/main_figures/` is this folder; the path mapping is given in the repository [README](../../README.md#code-by-analysis).
