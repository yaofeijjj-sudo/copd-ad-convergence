source("R/runtime.R")
.libPaths(c("env/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})
out <- "results/genetics_blood"
old <- "results/blood_ranks"
axis <- as.character(fread("results/module_overlap/results/cross/AxisA_confirmed_genes.csv")$entrez_id)
stopifnot(length(axis) == 137, as.character(packageVersion("fgsea")) == "1.36.2")
ids <- c("GSE54837", "GSE248493", "E-MTAB-5278", "GSE63060", "GSE63061", "GSE248417")
oldstats <- fread(file.path(old, "gene_stats_retained_five.csv.gz"))
set.seed(params$blood_seed)
ans <- list()
coverage <- list()
for (i in seq_along(ids)) {
  id <- ids[i]
  d <- switch(id,
    GSE248493 = fread(file.path(out, "results/gene_stats_GSE248493.csv.gz")),
    GSE248417 = fread(file.path(out, "results/gene_stats_primary_age_sex_97.csv.gz")),
    `E-MTAB-5278` = fread(file.path(out, "results/gene_stats_emtab5278.csv.gz")),
    oldstats[dataset_id == id & model_id == "age_gender_adjusted"]
  )
  stopifnot(nrow(d) > 1000, !anyNA(d$t), uniqueN(d$entrez_id) == nrow(d))
  r <- d$t
  names(r) <- as.character(d$entrez_id)
  r <- sort(r[is.finite(r)], decreasing = TRUE)
  present <- intersect(axis, names(r))
  fg <- as.data.table(fgseaMultilevel(pathways = list(AxisA_Entrez137 = present), stats = r, minSize = params$fgsea_min_size, maxSize = params$fgsea_max_size, eps = params$fgsea_eps))
  stopifnot(nrow(fg) == 1, fg$size == length(present), fg$pval == fg$padj)
  ans[[i]] <- data.table(dataset_id = id, disease = if (i <= 3) "COPD" else "AD", axis_input = 137, axis_present = fg$size, gene_universe = length(r), fgsea_ES = fg$ES, fgsea_NES = fg$NES, fgsea_p = fg$pval, fgsea_padj = fg$padj, log2err = fg$log2err, leading_edge_entrez = paste(fg$leadingEdge[[1]], collapse = ";"), model = "Case-Control; age + sex")
  coverage[[i]] <- data.table(dataset_id = id, entrez_id = axis, present = axis %in% names(r))
  fwrite(d, file.path(out, "results", paste0("ranked_gene_stats_", id, ".csv.gz")))
}
results <- rbindlist(ans)
fwrite(results, file.path(out, "results/fgsea_six_results.csv"))
fwrite(rbindlist(coverage), file.path(out, "results/fgsea_coverage_137_by_six.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "provenance/fgsea_sessionInfo.txt"))
print(results[, .(dataset_id, axis_present, fgsea_NES, fgsea_p)])
