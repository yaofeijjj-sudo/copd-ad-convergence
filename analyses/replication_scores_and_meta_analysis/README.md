# Axis A scores and replication

Scores Axis A with singscore, fits the disease-association model in each independent cohort and pools the AD datasets.

**Used in:** Fig 3c, Table D

| Script | Content |
| --- | --- |
| [`00_functions.R`](code/00_functions.R) | data loading, scoring and models |
| [`02_run_specificity.R`](code/02_run_specificity.R) | disease-association models in all cohorts |
| [`03_GPL6480_fingerprint_check_20261027.R`](code/03_GPL6480_fingerprint_check_20261027.R) | GPL6480 donors compared with the discovery cohort |
| [`03_draw_forest.R`](code/03_draw_forest.R) | draft forest plot |
| [`03_draw_forest_comparators_only.R`](code/03_draw_forest_comparators_only.R) | draft forest plot |
| [`04_ad_reference_meta.R`](code/04_ad_reference_meta.R) | DerSimonian–Laird pool of the six AD datasets |
| [`05_draw_ad_meta.R`](code/05_draw_ad_meta.R) | draft AD forest plot |
| [`06_supplementary_tables.py`](code/06_supplementary_tables.py) | result tables |
