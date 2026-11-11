source("R/runtime.R")
.libPaths(c("env/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hdf5r)
})
OUT <- "results/single_cell"
WORK <- file.path(OUT, "results/brain_work")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
cross <- "results/module_overlap/results/cross"
axis <- fread(file.path(cross, "AxisA_confirmed_genes.csv"), colClasses = "character")
bg <- fread(file.path(cross, "shared_non_grey_background.csv"), colClasses = "character")
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(bg$entrez_id, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")), keytype = "ENTREZID", columns = "ENSEMBL"))
mp <- unique(mp[!is.na(ENSEMBL)])
mp[, ENSEMBL := sub("\\..*$", "", ENSEMBL)]
fwrite(mp, file.path(OUT, "provenance/Entrez_Ensembl_annotation.csv"))
lp <- list.files("data/raw/lung_h5", pattern = "h5$", full.names = TRUE)[1]
h <- H5File$new(lp, "r")
lf <- as.character(h[["matrix/features/id"]][])
h$close_all()
h <- H5File$new("data/raw/single_cell/sea_ad_pfc_a9/SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad", "r")
bf <- as.character(h[["var/gene_ids"]][])
h$close_all()
stopifnot(!anyDuplicated(lf), !anyDuplicated(bf))
for (organ in c("lung", "brain")) {
  features <- if (organ == "lung") lf else bf
  p <- mp[ENSEMBL %chin% features]
  p[, n_features := uniqueN(ENSEMBL), by = ENTREZID]
  p[, n_ids_for_feature := uniqueN(ENTREZID), by = ENSEMBL]
  ok <- p[n_features == 1 & n_ids_for_feature == 1, .(entrez_id = ENTREZID, ensembl_id = ENSEMBL)]
  ledger <- merge(bg, ok, by = "entrez_id", all.x = TRUE)
  ledger[, `:=`(axis_member = entrez_id %chin% axis$entrez_id, measurable = !is.na(ensembl_id), rule = "unique stable Entrez to unique measurable Ensembl; no aliases")]
  fwrite(ledger, file.path(OUT, "provenance", paste0(organ, "_stable_ID_coverage.csv")))
  writeLines(ok$ensembl_id, file.path(OUT, "results", paste0(organ, "_universe_measured.txt")))
  a <- ok[entrez_id %chin% axis$entrez_id]$ensembl_id
  writeLines(a, file.path(OUT, "results", paste0(organ, "_axis_measured.txt")))
  if (organ == "brain") {
    writeLines(features, file.path(WORK, "genes.txt"))
    writeLines(a, file.path(WORK, "axis_measured.txt"))
    writeLines(ok$ensembl_id, file.path(WORK, "universe_measured.txt"))
    fwrite(data.table(gene = character(), measured = logical()), file.path(WORK, "independent_markers.tsv"), sep = "\t")
  }
  cat(organ, "axis=", length(a), "/", nrow(axis), "background=", nrow(ok), "/", nrow(bg), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "provenance/stable_mapping_sessionInfo.txt"))
