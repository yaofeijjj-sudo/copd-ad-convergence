# Lung single-cell analyses with 55 donors

Reruns all lung donor-level analyses with `cood8` and `copd8` of GSE302339 merged into one donor (two lobes of the same person).

**Used in:** Fig 4–5, Fig 8, Table F, Table L, Table M

| Script | Content |
| --- | --- |
| [`01_localization_and_donor_tests.R`](code/01_localization_and_donor_tests.R) | localization, matched random sets and donor tests |
| [`02_pseudobulk_AM_TF.R`](code/02_pseudobulk_AM_TF.R) | pseudobulk, alveolar macrophage and TF activity models |
| [`03_celloracle_lung_reaverage.py`](code/03_celloracle_lung_reaverage.py) | CellOracle lung readout averaged over 55 donors |
| [`04_localization_table.py`](code/04_localization_table.py) | two-organ localization table |
