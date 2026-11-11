# CIBERSORTx mixture files (linear scale, gene symbols) for the 8 tissue cohorts, same samples as the adopted replication.
# Arrays: 2^x of the analysed log2 matrix; RNA-seq: CPM of the analysed count matrix (all genes with counts).
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/04_mixtures.R
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db) })
OUT <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
coh <- readRDS("03_candidate_analysis/AxisA_member_direction_robustness_20261018/source_data/cohorts.rds")$cohorts
ids <- c("GSE76925", "GSE47460", "GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
info <- rbindlist(lapply(ids, function(k) { c <- coh[[k]]; e <- c$expr[, c$dat$sample, drop = FALSE]
  sym <- mapIds(org.Hs.eg.db, rownames(e), "SYMBOL", "ENTREZID"); ok <- !is.na(sym) & !sym %in% sym[duplicated(sym)]; e <- e[ok, , drop = FALSE]; rownames(e) <- sym[ok]
  if (c$kind == "counts") { e <- e[rowSums(e) > 0, ]; lin <- cpm(e); scale_note <- "CPM" } else {
    stopifnot(max(e, na.rm = TRUE) < 25); lin <- 2^e; scale_note <- sprintf("2^x (log2 range %.1f-%.1f)", min(e), max(e)) }
  stopifnot(!anyNA(lin))
  sname <- gsub("[^A-Za-z0-9_.-]", "_", gsub(" +", "_", trimws(colnames(lin))))
  fwrite(data.table(sample = colnames(lin), upload_name = sname), file.path(OUT, sprintf("source_data/%s_mixture_samples.tsv", k)), sep = "\t")
  f <- file.path(OUT, sprintf("upload/mix_%s.txt", k)); con <- file(f, "w"); writeLines(paste(c("GeneSymbol", sname), collapse = "\t"), con)
  write.table(signif(lin, 5), con, sep = "\t", quote = FALSE, col.names = FALSE); close(con)
  data.table(cohort = k, kind = c$kind, scale = scale_note, genes = nrow(lin), samples = ncol(lin), datapoints = nrow(lin) * ncol(lin), MB = round(file.size(f) / 1e6, 1),
             quantile_normalization = c$kind != "counts") }))
fwrite(info, file.path(OUT, "qa/mixture_files.tsv"), sep = "\t"); print(info)
