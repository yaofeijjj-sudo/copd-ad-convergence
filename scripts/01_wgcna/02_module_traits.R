source("R/runtime.R")
source("R/paths.R")
suppressPackageStartupMessages(library(data.table))
dir.create("results/wgcna/source_data", recursive = TRUE, showWarnings = FALSE)
for (side in c("lung", "brain")) {
  path <- if (side == "lung") COPD_WGCNA else AD_WGCNA
  prefix <- if (side == "lung") "copd" else "ad"
  d <- fread(file.path(path, paste0(prefix, "_module_trait_tidy.csv")))
  fwrite(d, file.path("results/wgcna/source_data", paste0(side, "_heatmap.csv")))
}
