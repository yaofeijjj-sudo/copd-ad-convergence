source("R/runtime.R")
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
O <- file.path(getwd(), "results/single_cell")
changes <- list()
for (mode in c("broad", "subtype")) {
  p <- file.path(O, "results/brain_work", paste0(mode, "_ucell/cell_scores.tsv"))
  d <- fread(p)
  if ("axisA_common129_ucell" %in% names(d)) {
    stopifnot(all(d$axisA_common129_ucell == d$axisA_full_measured_ucell))
    setnames(d, "axisA_common129_ucell", "axisA_equivalence136_ucell")
    fwrite(d, p, sep = "\t")
    changes[[mode]] <- list(old = "axisA_common129_ucell", new = "axisA_equivalence136_ucell", reason = "An identical 136-member equivalence field inherited the original script header; no 129-member analysis was run.", primary_score_unchanged = TRUE)
  }
}
write_json(changes, file.path(O, "qa/HEADER_NORMALIZATION.json"), auto_unbox = TRUE, pretty = TRUE)
