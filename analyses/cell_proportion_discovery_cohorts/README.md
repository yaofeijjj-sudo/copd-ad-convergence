# Cell proportions in the discovery cohorts

Estimates the BRETIGEA microglial proportion in GSE33000 and adds it to the blue-module model.

**Used in:** Table G

| Script | Content |
| --- | --- |
| [`00_count_marker_overlap.R`](code/00_count_marker_overlap.R) | marker overlap |
| [`01_rebuild_bulk.R`](code/01_rebuild_bulk.R) | full-gene matrices |
| [`02_lung_reference.R`](code/02_lung_reference.R) | lung reference read by `03` |
| [`03_deconvolve_and_model.R`](code/03_deconvolve_and_model.R) | proportions and adjusted models |
