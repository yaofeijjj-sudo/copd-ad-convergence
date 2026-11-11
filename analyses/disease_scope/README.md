# Disease scope

Applies the disease-association model to other lung and brain diseases; the GEO search for the AD datasets was run here.

**Used in:** Fig 7b, Table D

| Script | Content |
| --- | --- |
| [`00_functions.R`](code/00_functions.R) | shared functions |
| [`01_parse_metadata.R`](code/01_parse_metadata.R) | GEO metadata |
| [`02_run_specificity.R`](code/02_run_specificity.R) | disease models |
| [`03_draw_forest.R`](code/03_draw_forest.R) | draft forest plot |
| [`03_draw_forest_comparators_only.R`](code/03_draw_forest_comparators_only.R) | draft forest plot |
| [`04_ad_reference_meta.R`](code/04_ad_reference_meta.R) | AD reference pool |
| [`05_draw_ad_meta.R`](code/05_draw_ad_meta.R) | draft AD forest plot |
