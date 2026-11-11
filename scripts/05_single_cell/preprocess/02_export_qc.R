source("R/runtime.R")

suppressPackageStartupMessages(library(data.table))

project_root <- "."
result_dir <- file.path(
  project_root,
  "results/single_cell/preprocessing",
  "lung"
)
cell_dir <- file.path(result_dir, "objects", "cell_qc_by_library")
files <- sort(list.files(cell_dir, pattern = "_cell_qc\\.rds$", full.names = TRUE))
stopifnot(length(files) == 65L)

all_qc <- rbindlist(lapply(files, readRDS), use.names = TRUE, fill = TRUE)
stopifnot(nrow(all_qc) == 264520L)
fwrite(
  all_qc,
  file.path(result_dir, "objects", "cell_qc_all.tsv.gz"),
  sep = "\t",
  na = "NA",
  compress = "gzip"
)
cat(sprintf("Exported %s cells; %s QC-passing singlets\n", nrow(all_qc), sum(all_qc$passes_QC)))
