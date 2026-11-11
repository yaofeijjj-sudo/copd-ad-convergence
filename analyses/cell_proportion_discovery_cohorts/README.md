# Cell proportions in the discovery cohorts

Estimates cell proportions with dtangle in GSE76925 and BRETIGEA in GSE33000 and adds them to the module models.

**Used in:** Table G

| Script | Content |
| --- | --- |
| [`00_count_marker_overlap.R`](code/00_count_marker_overlap.R) | marker overlap |
| [`01_rebuild_bulk.R`](code/01_rebuild_bulk.R) | full-gene matrices |
| [`02_lung_reference.R`](code/02_lung_reference.R) | lung reference |
| [`03_deconvolve_and_model.R`](code/03_deconvolve_and_model.R) | proportions and adjusted models |
| [`04_collinearity_qa.R`](code/04_collinearity_qa.R) | collinearity checks |
