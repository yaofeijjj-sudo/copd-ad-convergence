source("R/runtime.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

project_root <- "."
result_dir <- file.path(
  project_root,
  "results/single_cell/preprocessing",
  "lung"
)
h5_dir <- file.path(project_root, "data/raw/lung_h5")
ledger_file <- file.path(result_dir, "tables", "sample_ledger.tsv")
cell_dir <- file.path(result_dir, "objects", "cell_qc_by_library")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

read_10x_h5 <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  shape <- as.integer(f[["matrix/shape"]][])
  mat <- new(
    "dgCMatrix",
    i = as.integer(f[["matrix/indices"]][]),
    p = as.integer(f[["matrix/indptr"]][]),
    x = as.numeric(f[["matrix/data"]][]),
    Dim = shape
  )
  feature_id <- as.character(f[["matrix/features/id"]][])
  gene_symbol <- as.character(f[["matrix/features/name"]][])
  barcodes <- as.character(f[["matrix/barcodes"]][])
  rownames(mat) <- make.unique(feature_id)
  colnames(mat) <- barcodes
  list(counts = mat, feature_id = feature_id, gene_symbol = gene_symbol, barcodes = barcodes)
}

extract_thresholds <- function(filters) {
  out <- list()
  for (nm in setdiff(colnames(filters), "discard")) {
    th <- attr(filters[[nm]], "thresholds")
    out[[nm]] <- data.table(
      filter = nm,
      lower = if (!is.null(th) && "lower" %in% names(th)) as.numeric(th[["lower"]]) else NA_real_,
      higher = if (!is.null(th) && "higher" %in% names(th)) as.numeric(th[["higher"]]) else NA_real_
    )
  }
  rbindlist(out, fill = TRUE)
}

ledger <- fread(ledger_file, na.strings = "NA")
ledger <- ledger[matrix_available == TRUE]
setorder(ledger, gsm)

summaries <- vector("list", nrow(ledger))
thresholds <- vector("list", nrow(ledger))
errors <- list()

for (idx in seq_len(nrow(ledger))) {
  row <- ledger[idx]
  gsm <- row$gsm
  out_file <- file.path(cell_dir, paste0(gsm, "_cell_qc.rds"))
  if (file.exists(out_file)) {
    cell_dt <- readRDS(out_file)
    message(sprintf("[%d/%d] %s checkpoint found", idx, nrow(ledger), gsm))
  } else {
    message(sprintf("[%d/%d] %s reading %s", idx, nrow(ledger), gsm, row$matrix_file))
    x <- read_10x_h5(file.path(h5_dir, row$matrix_file))
    colnames(x$counts) <- paste(gsm, x$barcodes, sep = "_")
    sce <- SingleCellExperiment(assays = list(counts = x$counts))
    rowData(sce)$feature_id <- x$feature_id
    rowData(sce)$gene_symbol <- x$gene_symbol
    mito <- grepl("^MT-", x$gene_symbol)
    qcm <- perCellQCMetrics(sce, subsets = list(Mito = mito))
    qcf <- quickPerCellQC(qcm, sub.fields = "subsets_Mito_percent", nmads = 3)
    thresholds[[idx]] <- cbind(gsm = gsm, extract_thresholds(qcf))

    dbl_score <- rep(NA_real_, ncol(sce))
    dbl_class <- rep(NA_character_, ncol(sce))
    core_pass <- !qcf$discard
    if (sum(core_pass) >= 100L) {
      set.seed(params$lung_qc_seed_base + idx)
      dbl <- tryCatch(
        scDblFinder(
          sce[, core_pass],
          samples = rep(gsm, sum(core_pass)),
          BPPARAM = SerialParam(progressbar = FALSE),
          verbose = FALSE
        ),
        error = function(e) e
      )
      if (inherits(dbl, "error")) {
        errors[[length(errors) + 1L]] <- data.table(gsm = gsm, stage = "scDblFinder", error = conditionMessage(dbl))
      } else {
        dbl_score[core_pass] <- colData(dbl)$scDblFinder.score
        dbl_class[core_pass] <- as.character(colData(dbl)$scDblFinder.class)
      }
    } else {
      errors[[length(errors) + 1L]] <- data.table(gsm = gsm, stage = "scDblFinder", error = "Fewer than 100 core-QC cells")
    }

    cell_dt <- data.table(
      cell_id = colnames(sce),
      barcode = x$barcodes,
      gsm = gsm,
      subject_id = row$subject_id,
      disease_group = row$disease_group,
      emphysema_group = row$emphysema_group,
      lobe = row$lobe,
      sum = as.numeric(qcm$sum),
      detected = as.integer(qcm$detected),
      percent_mito = as.numeric(qcm$subsets_Mito_percent),
      qc_low_lib_size = as.logical(qcf$low_lib_size),
      qc_low_n_features = as.logical(qcf$low_n_features),
      qc_high_mito = as.logical(qcf$high_subsets_Mito_percent),
      qc_discard = as.logical(qcf$discard),
      scDblFinder_score = dbl_score,
      scDblFinder_class = dbl_class
    )
    cell_dt[, passes_QC := !qc_discard & scDblFinder_class == "singlet"]
    saveRDS(cell_dt, out_file, compress = "xz")
    rm(x, sce, qcm, qcf, dbl)
    gc(verbose = FALSE)
  }

  summaries[[idx]] <- cell_dt[, .(
    n_input = .N,
    n_core_qc = sum(!qc_discard),
    n_doublet = sum(scDblFinder_class == "doublet", na.rm = TRUE),
    n_pass = sum(passes_QC, na.rm = TRUE),
    median_umi_input = median(sum),
    median_genes_input = median(detected),
    median_percent_mito_input = median(percent_mito),
    median_umi_pass = median(sum[passes_QC], na.rm = TRUE),
    median_genes_pass = median(detected[passes_QC], na.rm = TRUE),
    median_percent_mito_pass = median(percent_mito[passes_QC], na.rm = TRUE)
  )][, `:=`(
    gsm = gsm,
    subject_id = row$subject_id,
    disease_group = row$disease_group,
    emphysema_group = row$emphysema_group,
    lobe = row$lobe
  )]
}

summary_dt <- rbindlist(summaries, fill = TRUE)
setcolorder(summary_dt, c("gsm", "subject_id", "disease_group", "emphysema_group", "lobe"))
threshold_dt <- rbindlist(thresholds, fill = TRUE)
error_dt <- rbindlist(errors, fill = TRUE)

fwrite(summary_dt, file.path(result_dir, "tables", "qc_library_summary.tsv"), sep = "\t", na = "NA")
fwrite(threshold_dt, file.path(result_dir, "qa", "qc_thresholds_by_library.tsv"), sep = "\t", na = "NA")
fwrite(error_dt, file.path(result_dir, "qa", "qc_processing_errors.tsv"), sep = "\t", na = "NA")

overall <- summary_dt[, .(
  libraries = .N,
  donors = uniqueN(subject_id),
  n_input = sum(n_input),
  n_core_qc = sum(n_core_qc),
  n_doublet = sum(n_doublet),
  n_pass = sum(n_pass),
  pass_fraction = sum(n_pass) / sum(n_input)
), by = disease_group]
fwrite(overall, file.path(result_dir, "tables", "qc_overall_summary.tsv"), sep = "\t")

p1 <- ggplot(summary_dt, aes(x = reorder(gsm, n_input), y = n_input, fill = disease_group)) +
  geom_col(width = 0.8) +
  geom_col(aes(y = n_pass), width = 0.8, fill = "#238B8E") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = "Cells", title = "Input cells (background) and QC-passing singlets (teal)") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")

p2 <- ggplot(summary_dt, aes(x = median_genes_input, y = median_percent_mito_input, color = disease_group, size = n_input)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1.5, 6), labels = scales::label_comma()) +
  labs(x = "Median detected genes", y = "Median mitochondrial %", color = NULL, size = "Input cells") +
  theme_bw(base_size = 9)

ggsave(file.path(result_dir, "figures", "qc_library_overview.pdf"), p1 + p2 + plot_layout(widths = c(1.4, 1)), width = 12, height = 11)
ggsave(file.path(result_dir, "figures", "qc_library_overview.png"), p1 + p2 + plot_layout(widths = c(1.4, 1)), width = 12, height = 11, dpi = 180)

cat("QC and per-library scDblFinder complete\n")
print(overall)
if (nrow(error_dt)) {
  cat("\nProcessing errors:\n")
  print(error_dt)
  quit(status = 2L)
}
