# Lung GSE302339: library-level QC check and myeloid donor test without flagged libraries (2026-10-28).
# Library check follows OSCA (Amezquita et al. 2020) advice to compare per-batch QC across batches: a library is flagged when its
# post-QC median mitochondrial percentage exceeds the median + 3 MAD across the 65 libraries (log scale).
# Run from project root: Rscript 03_candidate_analysis/Lung_sc_library_QC_check_20261028/code/run.R
suppressPackageStartupMessages(library(data.table))
SC <- "细胞定位/analysis/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate"
O <- "03_candidate_analysis/Lung_sc_library_QC_check_20261028/results"
lib <- fread(file.path(SC, "tables/qc_library_summary.tsv"))
thr <- dcast(fread(file.path(SC, "qa/qc_thresholds_by_library.tsv")), gsm ~ filter, value.var = c("lower", "higher"))
lm <- log(lib$median_percent_mito_pass); cut <- exp(median(lm) + 3 * mad(lm))
lib[, flagged := median_percent_mito_pass > cut]
lib <- merge(lib, thr[, .(gsm, min_umi = lower_low_lib_size, min_genes = lower_low_n_features, max_mito_pct = higher_high_subsets_Mito_percent)], by = "gsm")
fwrite(lib[order(gsm)], file.path(O, "library_qc_with_flag.tsv"), sep = "\t")
cat("cut (median mito %):", round(cut, 2), "| flagged:", paste(lib[flagged == TRUE]$gsm, lib[flagged == TRUE]$subject_id, collapse = "; "), "\n")
x <- fread(file.path(SC, "objects/axisA_cell_scores.tsv.gz"), select = c("gsm", "subject_id", "disease_group", "passes_QC", "broad_celltype", "AxisA_UCell"))
x <- x[passes_QC == TRUE & broad_celltype == "Myeloid"]
run <- function(d, lab) {  # same aggregation as 05_axisA_official_ucell.R: library median, then median across a donor's libraries, >= 20 cells
  l <- d[, .(n = .N, med = median(AxisA_UCell)), by = .(gsm, subject_id, disease_group)]
  m <- l[, .(n = sum(n), med = median(med)), by = .(subject_id, disease_group)][n >= 20]
  w <- wilcox.test(med ~ factor(disease_group, levels = c("Control", "COPD")), data = m, exact = FALSE)
  data.table(analysis = lab, donors = nrow(m), control = sum(m$disease_group == "Control"), copd = sum(m$disease_group == "COPD"), cells = sum(m$n),
             median_control = median(m[disease_group == "Control"]$med), median_copd = median(m[disease_group == "COPD"]$med),
             delta = median(m[disease_group == "COPD"]$med) - median(m[disease_group == "Control"]$med), P = w$p.value)
}
res <- rbind(run(x, "all QC-passed libraries (reported)"), run(x[!gsm %in% lib[flagged == TRUE]$gsm], "flagged libraries removed"))
stopifnot(abs(res$P[1] - 0.0103236129830834) < 1e-10)
fwrite(res, file.path(O, "myeloid_donor_test_library_check.tsv"), sep = "\t"); print(res)
lab <- fread(file.path(SC, "tables/celltypist_label_counts.tsv")); fwrite(lab, file.path(O, "celltypist_label_to_broad.tsv"), sep = "\t")
