# Run official sigQC 0.1.24 on Axis A in the 8 independent datasets (metrics, plots, negative and permutation controls).
# Summaries, the expression-matched null and the gene-level table are in 03_summarise.R.
# Design: ../DESIGN.md (frozen; sha256 in ../DESIGN.sha256). Input: matrices from 01_build_matrices.R.
# Run from project root:
#   R_LIBS=_workspace/r_libs_sigqc:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/02_run_sigqc.R
# (RankProd, a sigQC 'Suggests' package, is not installed: it needs the system GMP library; sigQC then skips only its rank-product table.)
suppressPackageStartupMessages({ library(data.table); library(sigQC) })
OUT <- "03_candidate_analysis/AxisA_sigQC_coherence_20261014"
SQ <- file.path(OUT, "sigqc_output"); dir.create(SQ, showWarnings = FALSE); dir.create(file.path(OUT, "results"), showWarnings = FALSE)
ax <- fread("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv", colClasses = "character")
axisA <- ax$entrez_id; stopifnot(length(axisA) == 137)
mats <- readRDS("_workspace/sigqc_axisA_20261014/matrices.rds")
# sigQC captures base plots with gridGraphics::grid.echo(), which needs the display list; non-interactive
# Rscript devices disable it, so the compactness step failed. Enable it on devices opened by dev.new().
options(device = function(...) { grDevices::pdf(file = tempfile(fileext = ".pdf"), ...); grDevices::dev.control("enable") })

# ---- 1. official sigQC (metrics, plots, negative + permutation controls) ----
set.seed(20261014)
make_all_plots(gene_sigs_list = list(AxisA = as.matrix(axisA)), mRNA_expr_matrix = mats, names_sigs = "AxisA",
               names_datasets = names(mats), out_dir = SQ, showResults = FALSE, doNegativeControl = TRUE, numResampling = 200)
stopifnot(!any(grepl("Error occurred", readLines(file.path(SQ, "log.log")))))
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo_sigqc.txt"))
