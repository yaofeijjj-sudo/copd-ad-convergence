# Summaries of the sigQC run, the pre-specified expression-matched null and the gene-level table. Design: ../DESIGN.md.
# Run from project root after 02_run_sigqc.R:
#   R_LIBS=_workspace/r_libs_sigqc:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/03_summarise.R
suppressPackageStartupMessages(library(data.table))
OUT <- "03_candidate_analysis/AxisA_sigQC_coherence_20261014"
SQ <- file.path(OUT, "sigqc_output"); dir.create(file.path(OUT, "results"), showWarnings = FALSE)
ax <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")
axisA <- ax$entrez_id; stopifnot(length(axisA) == 137)
mats <- readRDS("_workspace/sigqc_axisA_20261014/matrices.rds")
tissue <- c(`GSE47460-GPL14550` = "lung", `GSE47460-GPL6480` = "lung")

rt <- read.table(file.path(SQ, "radarchart_table", "radarchart_table.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
radar <- data.table(row = rownames(rt), rt)[grepl("_AxisA$", row)]; stopifnot(nrow(radar) == length(mats))
fwrite(radar, file.path(OUT, "results/sigqc_radar_metrics.tsv"), sep = "\t")
nc <- rbindlist(lapply(names(mats), function(d) {
  f <- file.path(SQ, "negative_control", d, "AxisA", "neg_controls_summary_table.txt")
  p <- file.path(SQ, "permutation_control", d, "AxisA", "perm_controls_summary_table.txt")
  g <- function(file, type) if (file.exists(file)) { x <- read.csv(file, check.names = FALSE); data.table(dataset = d, control = type, stat = rownames(x), autocor_median = x[["autocor_median"]]) }
  rbind(g(f, "sigQC random genes"), g(p, "sigQC permutation"))
}))
fwrite(nc, file.path(OUT, "results/sigqc_control_summaries_autocor.tsv"), sep = "\t")

# ---- 2. expression-matched null for sigQC's autocor_median (same definition: median of the full Spearman matrix) ----
autocor_median <- function(e) stats::median(stats::cor(t(e), method = "spearman"), na.rm = TRUE)
set.seed(20261014)
setlev <- rbindlist(lapply(names(mats), function(d) {
  e <- mats[[d]]; idx <- which(rownames(e) %in% axisA)
  obs <- autocor_median(e[idx, , drop = FALSE])
  mu <- rowMeans(e); bins <- cut(rank(mu, ties.method = "first"), 20, labels = FALSE)
  need <- table(bins[idx]); rest <- setdiff(seq_len(nrow(e)), idx); pool <- split(rest, bins[rest])
  null <- vapply(1:1000, function(b) autocor_median(e[unlist(lapply(names(need), function(k) sample(pool[[k]], need[[k]]))), , drop = FALSE]), 0)
  data.table(dataset = d, tissue = ifelse(d %in% names(tissue), "lung", "frontal cortex"), n_samples = ncol(e), axisA_measured = length(idx),
             autocor_median = obs, null_mean = mean(null), null_q95 = unname(quantile(null, 0.95)), null_max = max(null),
             empirical_p_one_sided = (1 + sum(null >= obs)) / 1001, coherent = obs > quantile(null, 0.95))
}))
radar_ac <- radar[, .(dataset = sub("_AxisA$", "", row), sigqc_autocor_median = autocor_median)]
setlev <- merge(setlev, radar_ac, by = "dataset", all.x = TRUE, sort = FALSE)
setlev[, matches_sigqc := abs(autocor_median - sigqc_autocor_median) < 1e-8]
print(setlev); fwrite(setlev, file.path(OUT, "results/set_level_coherence.tsv"), sep = "\t")

# ---- 3. gene level: median Spearman correlation of each member with the other members, per dataset ----
gl <- rbindlist(lapply(names(mats), function(d) {
  e <- mats[[d]]; idx <- which(rownames(e) %in% axisA); r <- stats::cor(t(e[idx, , drop = FALSE]), method = "spearman"); diag(r) <- NA
  data.table(dataset = d, entrez_id = rownames(r), median_r_with_members = apply(r, 1, stats::median, na.rm = TRUE))
}))
wide <- dcast(gl, entrez_id ~ dataset, value.var = "median_r_with_members")
wide <- merge(ax[, .(entrez_id, symbol = current_symbol)], wide, by = "entrez_id", all.x = TRUE)
dcols <- names(mats); v <- as.matrix(wide[, ..dcols])
wide[, n_measured := rowSums(!is.na(v))]; wide[, n_positive := rowSums(v > 0, na.rm = TRUE)]
wide[, median_across_datasets := apply(v, 1, stats::median, na.rm = TRUE)]
wide[, min_across_datasets := suppressWarnings(apply(v, 1, min, na.rm = TRUE))]
wide[, class := fifelse(n_measured < 6, "insufficient coverage", fifelse(n_positive == n_measured, "consistent", "inconsistent"))]
setorder(wide, -median_across_datasets)
fwrite(wide, file.path(OUT, "results/gene_level_coherence.tsv"), sep = "\t")
print(wide[, .N, by = class])

# QA: gene-level correlations equal sigQC's own autocorrelation matrices
chk <- rbindlist(lapply(names(mats), function(d) {
  f <- file.path(SQ, "autocorrelation_matrices", paste0("autocorrelation_matrix_AxisA_", d, ".txt"))
  if (!file.exists(f)) return(data.table(dataset = d, sigqc_matrix_found = FALSE, max_abs_diff = NA_real_))
  s <- as.matrix(read.table(f, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1))
  e <- mats[[d]]; idx <- which(rownames(e) %in% axisA); r <- stats::cor(t(e[idx, , drop = FALSE]), method = "spearman")
  s <- s[rownames(r), colnames(r)]; data.table(dataset = d, sigqc_matrix_found = TRUE, max_abs_diff = max(abs(s - r)))
}))
print(chk); fwrite(chk, file.path(OUT, "qa/gene_level_vs_sigqc_autocor.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo_summarise.txt"))
