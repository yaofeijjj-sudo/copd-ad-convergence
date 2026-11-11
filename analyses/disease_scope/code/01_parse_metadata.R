# Parse GEO sample characteristics for the comparator datasets (design: ../DESIGN.md).
# Run from project root: _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/01_parse_metadata.R
suppressPackageStartupMessages(library(data.table))
RAW <- "02_data/raw/disease_specificity_tissue_20261003"
OUT <- "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/metadata"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
parse_series <- function(g) {
  l <- system2("sh", c("-c", shQuote(sprintf("gzip -cd %s/%s_series_matrix.txt.gz | grep '^!Sample_'", RAW, g))), stdout = TRUE)
  sp <- strsplit(l, "\t"); key <- sub("^!Sample_", "", vapply(sp, `[`, "", 1))
  vals <- lapply(sp, function(x) gsub('"', "", x[-1]))
  out <- data.table(gsm = vals[[which(key == "geo_accession")]], title = vals[[which(key == "title")]],
                    source = vals[[which(key == "source_name_ch1")[1]]])
  for (i in which(key == "characteristics_ch1")) {
    v <- vals[[i]]; k <- sub(":.*", "", v[nzchar(v)][1]); kk <- make.names(k)
    if (kk %in% names(out)) kk <- paste0(kk, "_", i)
    out[[kk]] <- ifelse(grepl(":", v), sub("^[^:]*: ?", "", v), NA_character_)
  }
  out
}
for (g in c("GSE150910", "GSE64810", "GSE68719", "GSE53987", "GSE10072")) {
  d <- parse_series(g); fwrite(d, file.path(OUT, paste0(g, "_sample_metadata.tsv")), sep = "\t")
  cat("==", g, nrow(d), "samples; fields:", paste(names(d), collapse = ", "), "\n")
}
