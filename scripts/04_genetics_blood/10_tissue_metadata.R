source("R/runtime.R")
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
})
out <- "results/tissue_scope/metadata"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
raw <- "data/raw/tissue_scope"

p <- pData(getGEO(filename = "data/raw/bulk_lung_copd/GSE47460-GPL14550_series_matrix.txt.gz", GSEMatrix = TRUE))
fwrite(as.data.table(p), file.path(out, "GSE47460_GPL14550_pheno.csv"))

# GSE84422 uses the 57 DLPFC samples retained for module preservation, with the age and sex coding of the network inputs.
v <- fread("results/wgcna/brain/AD_WGCNA_unadjusted/validation_sample_table_used.csv")
g <- getGEO(filename = "data/raw/bulk_brain_ad/GSE84422-GPL96_series_matrix.txt.gz", GSEMatrix = TRUE)
if (is.list(g)) g <- g[[1]]
p <- pData(g)[v$sample_id, , drop = FALSE]
stopifnot(nrow(v) == 57, all(p$`brain region:ch1` == "Dorsolateral Prefrontal Cortex"))
fwrite(
  data.table(
    sample_id = v$sample_id, neuropathological_category = p$`neuropathological category:ch1`,
    age = v$Age, sex_male = v$Sex
  ),
  file.path(out, "GSE84422_dlpfc.csv")
)

parse_series <- function(g) {
  cmd <- sprintf("gzip -cd %s/%s_series_matrix.txt.gz | grep '^!Sample_'", raw, g)
  fields <- strsplit(system2("sh", c("-c", shQuote(cmd)), stdout = TRUE), "\t")
  key <- sub("^!Sample_", "", vapply(fields, `[`, "", 1))
  values <- lapply(fields, function(x) gsub('"', "", x[-1]))
  d <- data.table(
    gsm = values[[which(key == "geo_accession")]], title = values[[which(key == "title")]],
    source = values[[which(key == "source_name_ch1")[1]]]
  )
  for (i in which(key == "characteristics_ch1")) {
    x <- values[[i]]
    k <- make.names(sub(":.*", "", x[nzchar(x)][1]))
    if (k %in% names(d)) k <- paste0(k, "_", i)
    d[[k]] <- ifelse(grepl(":", x), sub("^[^:]*: ?", "", x), NA_character_)
  }
  d
}
for (g in c("GSE150910", "GSE64810", "GSE68719", "GSE53987")) {
  fwrite(parse_series(g), file.path(out, paste0(g, "_sample_metadata.tsv")), sep = "\t")
}
