# Lung GSE302339: myeloid donor test without the two flagged libraries, on the current Axis A (137 genes, SYN01 Entrez).
# Correction (2026-10-28): run.R read cell scores from the September 135-gene run; its "reported" row reproduced only the P value.
# This script uses the SYN01 library-level medians behind S6 Table (a) and gates on all reported values.
# Library flags come from run.R (results/library_qc_with_flag.tsv); they depend on QC metrics only, not on the gene set.
# Run from project root: Rscript 03_candidate_analysis/Lung_sc_library_QC_check_20261028/code/run_axisA137.R
suppressPackageStartupMessages(library(data.table))
O <- "03_candidate_analysis/Lung_sc_library_QC_check_20261028/results"
lib_med <- fread("cell_localization/SYN01_Entrez_downstream_20261001/results/lung/tables/axisA_library_broad_medians.tsv")[broad_celltype == "Myeloid"]
flagged <- fread(file.path(O, "library_qc_with_flag.tsv"))[flagged == TRUE]$gsm
stopifnot(length(flagged) == 2)
run <- function(l, lab) {  # same aggregation as lung_rank_and_scores.R: library median, then median across a donor's libraries, >= 20 cells
  m <- l[, .(n = sum(n_cells), med = median(median_axisA_ucell)), by = .(subject_id, disease_group)][n >= 20]
  w <- wilcox.test(med ~ factor(disease_group, levels = c("Control", "COPD")), data = m, exact = FALSE)
  data.table(analysis = lab, donors = nrow(m), control = sum(m$disease_group == "Control"), copd = sum(m$disease_group == "COPD"), cells = sum(m$n),
             median_control = median(m[disease_group == "Control"]$med), median_copd = median(m[disease_group == "COPD"]$med),
             delta = median(m[disease_group == "COPD"]$med) - median(m[disease_group == "Control"]$med), P = w$p.value)
}
res <- rbind(run(lib_med, "all QC-passed libraries (as in S6 Table a)"), run(lib_med[!gsm %in% flagged], "flagged libraries removed"))
# gate: first row must reproduce S6 Table (a) / Disease_Lung_all_broad_COPD_vs_control.tsv
stopifnot(res$donors[1] == 56, res$cells[1] == 116558,
          abs(res$median_control[1] - 0.163187524398512) < 1e-12, abs(res$median_copd[1] - 0.147280978406032) < 1e-12,
          abs(res$delta[1] - (-0.0159065459924801)) < 1e-12, abs(res$P[1] - 0.0103236129830834) < 1e-10)
fwrite(res, file.path(O, "myeloid_donor_test_library_check.tsv"), sep = "\t"); print(res)
