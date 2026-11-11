source("R/runtime.R")

source("R/paths.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(GEOquery)
  library(impute)
  library(limma)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(params$wgcna_seed)

N_THREADS <- max(1L, min(4L, parallel::detectCores(logical = TRUE)))
allowWGCNAThreads(nThreads = N_THREADS)

SOFT_R2_CUTOFF <- 0.80
COPD_SOFT_R2_CUTOFF <- as.numeric(Sys.getenv("WGCNA_COPD_SOFT_R2_CUTOFF", "0.85"))
AD_SOFT_R2_CUTOFF <- as.numeric(Sys.getenv("WGCNA_AD_SOFT_R2_CUTOFF", "0.80"))
POWERS <- c(1:10, seq(12, 20, by = 2))
MIN_MODULE_SIZE <- 30L
MERGE_CUT_HEIGHT <- 0.25
DEEP_SPLIT <- 2L
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
PAM_RESPECTS_DENDRO <- FALSE
MAX_BLOCK_SIZE <- 6000L
LOW_VARIABILITY_FRACTION <- 0.30
MAX_MISSING_FRACTION <- 0.20
NEAR_ZERO_VARIANCE <- 1e-10
SAMPLE_CONNECTIVITY_Z_CUTOFF <- -2.5
MIN_MEAN_CONNECTIVITY <- as.numeric(Sys.getenv("WGCNA_MIN_MEAN_CONNECTIVITY", "1"))
PRESERVATION_N_PERMUTATIONS <- as.integer(Sys.getenv("WGCNA_N_PERMUTATIONS", "1000"))
READABLE_AXIS_CEX <- 1.15
READABLE_LAB_CEX <- 1.25
READABLE_MAIN_CEX <- 1.20
READABLE_POINT_LABEL_CEX <- 1.05
READABLE_HEATMAP_TEXT_CEX <- 0.95
READABLE_HEATMAP_LABEL_CEX <- 1.05
COPD_MODULE_TRAIT_HEATMAP_TRAITS <- c(
  "copd_status", "GOLD", "FEV1_pp", "FEV1_FVC",
  "laa950", "perc15", "pi10", "packyears_z"
)
COPD_MODULE_TRAIT_HEATMAP_LABELS <- c(
  copd_status = "COPD status",
  GOLD = "GOLD",
  FEV1_pp = "FEV1 %pred",
  FEV1_FVC = "FEV1/FVC",
  laa950 = "LAA-950",
  perc15 = "Perc15",
  pi10 = "Pi10",
  packyears_z = "Pack-years"
)
COR_METHOD <- tolower(Sys.getenv("WGCNA_COR_METHOD", "pearson"))
if (!COR_METHOD %in% c("pearson", "bicor")) {
  stop("WGCNA_COR_METHOD must be 'pearson' or 'bicor'.")
}

env_flag <- function(name) identical(Sys.getenv(name), "1")
STOP_AFTER_SOFT <- env_flag("WGCNA_STOP_AFTER_SOFT_THRESHOLD")
SKIP_PRESERVATION <- Sys.getenv("WGCNA_SKIP_PRESERVATION", "1") == "1"
CLEAN_OUTPUT <- env_flag("WGCNA_CLEAN_OUTPUT")
# The published brain network is built without age/sex residualization (module preservation is run by 04/05).
AD_SKIP_RESIDUALIZATION <- Sys.getenv("WGCNA_AD_SKIP_RESIDUALIZATION", "1") == "1"

message("[wgcna] Threads: ", N_THREADS)
message("[wgcna] Correlation: ", COR_METHOD)
message("[wgcna] Preservation permutations: ", PRESERVATION_N_PERMUTATIONS)
if (AD_SKIP_RESIDUALIZATION) {
  message("[wgcna] AD network: age outliers removed, no age/sex residualization.")
}

ensure_dir <- function(path, clean = FALSE) {
  if (clean && dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_dt <- function(x, path) {
  data.table::fwrite(data.table::as.data.table(x), path)
}

write_matrix <- function(mat, path, row_name = "module") {
  out <- data.table::as.data.table(mat, keep.rownames = row_name)
  data.table::fwrite(out, path)
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))
}

sex_to_male <- function(x) {
  lx <- tolower(trimws(as.character(x)))
  fifelse(
    lx %in% c("f", "female", "2", "2-female"), 0,
    fifelse(lx %in% c("m", "male", "1", "1-male"), 1, NA_real_)
  )
}

first_token <- function(x) {
  x <- as.character(x)
  x <- sub("\\s*///.*$", "", x)
  x <- sub("\\s*//.*$", "", x)
  trimws(x)
}

parse_all_characteristics <- function(eset) {
  pheno <- Biobase::pData(eset)
  char_cols <- grep("characteristics_ch", colnames(pheno), value = TRUE)
  parsed <- vector("list", nrow(pheno))
  names(parsed) <- rownames(pheno)

  for (i in seq_len(nrow(pheno))) {
    sample_data <- list()
    for (col in char_cols) {
      val <- as.character(pheno[i, col])
      if (is.na(val) || !nzchar(trimws(val))) next
      parts <- regmatches(val, regexec("^\\s*([^:]+):\\s*(.*)$", val))[[1]]
      if (length(parts) == 3) {
        key <- tolower(trimws(parts[2]))
        sample_data[[key]] <- trimws(parts[3])
      }
    }
    parsed[[i]] <- sample_data
  }

  all_keys <- unique(unlist(lapply(parsed, names), use.names = FALSE))
  out <- data.frame(row.names = rownames(pheno), check.names = FALSE)
  for (key in all_keys) {
    out[[key]] <- vapply(parsed, function(x) {
      if (key %in% names(x)) x[[key]] else NA_character_
    }, character(1))
  }
  out
}

normalize_microarray_matrix <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  if (max(expr_mat, na.rm = TRUE) > 100) {
    expr_mat <- log2(pmax(expr_mat, 1))
  }
  expr_mat <- limma::normalizeBetweenArrays(expr_mat, method = "quantile")
  keep <- rowMeans(is.na(expr_mat)) <= MAX_MISSING_FRACTION
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (anyNA(expr_mat)) {
    expr_mat <- impute::impute.knn(expr_mat)$data
  }
  expr_mat
}

collapse_probes_to_genes <- function(expr_mat, gene_id, gene_symbol = NULL,
                                     select_by = c("variance", "mean")) {
  select_by <- match.arg(select_by)
  gene_id <- first_token(gene_id)
  gene_symbol <- if (is.null(gene_symbol)) gene_id else first_token(gene_symbol)
  valid <- !is.na(gene_id) & nzchar(gene_id) & !gene_id %in% c("NA", "---")
  expr_mat <- expr_mat[valid, , drop = FALSE]
  gene_id <- gene_id[valid]
  gene_symbol <- gene_symbol[valid]

  metric <- if (select_by == "variance") {
    apply(expr_mat, 1, var, na.rm = TRUE)
  } else {
    rowMeans(expr_mat, na.rm = TRUE)
  }

  probe_table <- data.table(
    probe_id = rownames(expr_mat),
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    selection_metric = metric
  )
  setorder(probe_table, gene_id, -selection_metric, probe_id)
  selected <- probe_table[!duplicated(gene_id)]
  out <- expr_mat[selected$probe_id, , drop = FALSE]
  rownames(out) <- selected$gene_id
  selected[, selected_probe_id := probe_id]
  selected[, probe_id := NULL]
  list(expr = out, map = selected[])
}

filter_common_variable_genes <- function(train_gene_expr, valid_gene_expr, out_dir) {
  common_genes <- intersect(rownames(train_gene_expr), rownames(valid_gene_expr))
  train_common <- train_gene_expr[common_genes, , drop = FALSE]
  valid_common <- valid_gene_expr[common_genes, , drop = FALSE]

  missing_fraction <- rowMeans(is.na(train_common))
  variance <- apply(train_common, 1, var, na.rm = TRUE)
  mad_value <- apply(train_common, 1, mad, na.rm = TRUE)
  high_missing <- missing_fraction > MAX_MISSING_FRACTION
  near_zero <- is.na(variance) | variance <= NEAR_ZERO_VARIANCE
  finite_mad <- is.finite(mad_value)
  eligible <- !(high_missing | near_zero) & finite_mad

  mad_cutoff <- as.numeric(quantile(mad_value[eligible], LOW_VARIABILITY_FRACTION,
    na.rm = TRUE, names = FALSE
  ))
  retained <- eligible & mad_value >= mad_cutoff
  filter_reason <- rep("retained", length(common_genes))
  filter_reason[high_missing] <- "high_missing_fraction"
  filter_reason[!high_missing & near_zero] <- "zero_or_near_zero_variance"
  filter_reason[eligible & mad_value < mad_cutoff] <- "low_mad_30pct"
  filter_reason[!finite_mad] <- "nonfinite_mad"

  gene_filter <- data.table(
    gene_id = common_genes,
    missing_fraction = missing_fraction,
    variance = variance,
    mad = mad_value,
    retained = retained,
    filter_reason = filter_reason
  )
  write_dt(gene_filter, file.path(out_dir, "gene_filter_table.csv"))
  write_dt(gene_filter[retained == TRUE], file.path(out_dir, "gene_filter_retained.csv"))
  write_dt(gene_filter[retained == FALSE], file.path(out_dir, "gene_filter_removed.csv"))

  summary <- data.table(
    metric = c(
      "common_genes_before_filter",
      "high_missing_fraction_genes",
      "zero_or_near_zero_variance_genes",
      "mad_30pct_cutoff",
      "low_mad_30pct_genes",
      "genes_entering_wgcna"
    ),
    value = as.character(c(
      length(common_genes),
      sum(high_missing),
      sum(near_zero & !high_missing),
      signif(mad_cutoff, 8),
      sum(eligible & mad_value < mad_cutoff),
      sum(retained)
    ))
  )
  write_dt(summary, file.path(out_dir, "gene_filter_summary.csv"))

  list(
    train = t(train_common[retained, , drop = FALSE]),
    valid = t(valid_common[retained, , drop = FALSE]),
    genes = common_genes[retained],
    filter_table = gene_filter
  )
}

sample_connectivity_qc <- function(dat_expr) {
  sample_cor <- suppressWarnings(cor(t(dat_expr), use = "pairwise.complete.obs"))
  sample_adj <- (1 + sample_cor) / 2
  diag(sample_adj) <- NA_real_
  connectivity <- rowSums(sample_adj, na.rm = TRUE)
  z_connectivity <- as.numeric(scale(connectivity))
  data.table(
    sample_id = rownames(dat_expr),
    sample_connectivity = connectivity,
    sample_connectivity_z = z_connectivity,
    sample_outlier_z_lt_minus_2_5 = z_connectivity < SAMPLE_CONNECTIVITY_Z_CUTOFF
  )
}

apply_wgcna_sample_qc <- function(dat_expr, valid_expr, traits, out_dir, dataset_label) {
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 3)
  removed <- list()
  if (!gsg$allOK) {
    removed$goodSamplesGenes_samples <- rownames(dat_expr)[!gsg$goodSamples]
    removed$goodSamplesGenes_genes <- colnames(dat_expr)[!gsg$goodGenes]
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    valid_expr <- valid_expr[, gsg$goodGenes, drop = FALSE]
    traits <- traits[rownames(dat_expr), , drop = FALSE]
  } else {
    removed$goodSamplesGenes_samples <- character(0)
    removed$goodSamplesGenes_genes <- character(0)
  }

  sample_qc <- sample_connectivity_qc(dat_expr)
  write_dt(sample_qc, file.path(out_dir, "sample_connectivity_qc.csv"))
  outlier_ids <- sample_qc[sample_outlier_z_lt_minus_2_5 == TRUE, sample_id]
  if (length(outlier_ids) > 0 && length(outlier_ids) < 0.20 * nrow(dat_expr)) {
    message("[", dataset_label, "] Removing sample connectivity outliers: ", paste(outlier_ids, collapse = ", "))
    keep <- !rownames(dat_expr) %in% outlier_ids
    dat_expr <- dat_expr[keep, , drop = FALSE]
    traits <- traits[rownames(dat_expr), , drop = FALSE]
  } else if (length(outlier_ids) > 0) {
    warning("[", dataset_label, "] Connectivity outliers exceed 20% of samples; no automatic removal.")
    outlier_ids <- character(0)
  }
  removed$connectivity_outliers <- outlier_ids

  write_dt(data.table(
    removal_step = c("goodSamplesGenes_samples", "goodSamplesGenes_genes", "connectivity_outliers"),
    n_removed = c(
      length(removed$goodSamplesGenes_samples),
      length(removed$goodSamplesGenes_genes),
      length(removed$connectivity_outliers)
    ),
    removed_ids = c(
      paste(removed$goodSamplesGenes_samples, collapse = ";"),
      paste(removed$goodSamplesGenes_genes, collapse = ";"),
      paste(removed$connectivity_outliers, collapse = ";")
    )
  ), file.path(out_dir, "wgcna_qc_removal_summary.csv"))

  list(train = dat_expr, valid = valid_expr, traits = traits, removed = removed)
}

plot_sample_dendrogram <- function(dat_expr, traits, out_dir, file_prefix, title) {
  sample_tree <- hclust(dist(dat_expr), method = "average")
  pdf(file.path(out_dir, paste0(file_prefix, "_sample_dendrogram.pdf")), width = 12, height = 8)
  on.exit(dev.off(), add = TRUE)
  plot(sample_tree, main = title, sub = "", xlab = "", labels = FALSE)
}

plot_sample_dendrogram_with_traits <- function(dat_expr, traits, trait_cols, out_dir, file_prefix, title) {
  trait_cols <- intersect(trait_cols, colnames(traits))
  if (length(trait_cols) == 0) {
    return(invisible(NULL))
  }
  color_df <- data.frame(row.names = rownames(dat_expr), check.names = FALSE)
  for (tc in trait_cols) {
    vals <- traits[rownames(dat_expr), tc]
    color_df[[tc]] <- WGCNA::numbers2colors(as.numeric(vals), signed = FALSE)
  }
  sample_tree <- hclust(dist(dat_expr), method = "average")
  pdf(file.path(out_dir, paste0(file_prefix, "_sample_dendrogram_traits.pdf")), width = 14, height = 9)
  on.exit(dev.off(), add = TRUE)
  WGCNA::plotDendroAndColors(sample_tree, as.matrix(color_df),
    groupLabels = colnames(color_df),
    main = title,
    dendroLabels = FALSE
  )
}

cor_settings <- function(method) {
  if (method == "bicor") {
    list(corFnc = "bicor", corOptions = list(use = "pairwise.complete.obs", maxPOutliers = 0.1))
  } else {
    list(corFnc = "cor", corOptions = list(use = "pairwise.complete.obs"))
  }
}

select_soft_power <- function(dat_expr, out_dir, file_prefix, dataset_label,
                              cor_method = COR_METHOD,
                              r2_cutoff = SOFT_R2_CUTOFF,
                              forced_power = NULL,
                              forced_rule = NULL) {
  cs <- cor_settings(cor_method)
  message("[", dataset_label, "] pickSoftThreshold using ", cor_method)
  sft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = POWERS,
    networkType = NETWORK_TYPE,
    corFnc = cs$corFnc,
    corOptions = cs$corOptions,
    RsquaredCut = r2_cutoff,
    verbose = 3
  )

  fit <- as.data.table(sft$fitIndices)
  if (!"Power" %in% names(fit)) setnames(fit, names(fit)[1], "Power")
  fit[, signed_R2 := -sign(slope) * SFT.R.sq]
  fit[, cor_method := cor_method]
  fit[, selected := FALSE]
  fit[, selection_rule := NA_character_]

  ok <- is.finite(fit$signed_R2) & is.finite(fit$slope) & fit$slope < 0 & is.finite(fit$mean.k.)
  if (!is.null(forced_power)) {
    selected_idx <- match(forced_power, fit$Power)
    if (is.na(selected_idx)) {
      stop("[", dataset_label, "] forced_power=", forced_power, " is not in POWERS.")
    }
    rule <- if (!is.null(forced_rule)) forced_rule else "forced_primary_soft_power"
    if (is.finite(fit$signed_R2[selected_idx]) && fit$signed_R2[selected_idx] < r2_cutoff) {
      warning(
        "[", dataset_label, "] Forced power ", forced_power,
        " has signed R2=", round(fit$signed_R2[selected_idx], 3),
        " below cutoff ", r2_cutoff, "."
      )
    }
  } else {
    candidate <- which(ok & fit$signed_R2 >= r2_cutoff & fit$mean.k. >= MIN_MEAN_CONNECTIVITY)
    if (length(candidate) > 0) {
      selected_idx <- candidate[which.min(fit$Power[candidate])]
      rule <- paste0(
        "lowest_power_signed_R2_ge_", sprintf("%.2f", r2_cutoff),
        "_negative_slope_mean_connectivity_ge_min"
      )
    } else {
      candidate <- which(ok & fit$signed_R2 >= r2_cutoff)
      if (length(candidate) > 0) {
        selected_idx <- candidate[which.min(fit$Power[candidate])]
        rule <- paste0(
          "lowest_power_signed_R2_ge_", sprintf("%.2f", r2_cutoff),
          "_negative_slope_connectivity_warning"
        )
        warning("[", dataset_label, "] Selected power reaches R2 but mean connectivity is below threshold.")
      } else {
        selected_idx <- which.max(ifelse(ok, fit$signed_R2, -Inf))
        rule <- "max_signed_R2_negative_slope_fallback"
        warning(
          "[", dataset_label, "] No power reached signed R2 >= ", r2_cutoff,
          "; using max signed R2 fallback."
        )
      }
    }
  }

  selected_power <- fit$Power[selected_idx]
  fit[selected_idx, `:=`(selected = TRUE, selection_rule = rule)]
  write_dt(fit, file.path(out_dir, paste0(file_prefix, "_soft_power_selection.csv")))

  for (device in c("pdf", "png")) {
    path <- file.path(out_dir, paste0(file_prefix, "_soft_threshold.", device))
    if (device == "pdf") {
      pdf(path, width = 10, height = 5)
    } else {
      png(path, width = 2400, height = 1200, res = 240)
    }
    par(
      mfrow = c(1, 2),
      mar = c(5.5, 6.2, 4.2, 1.8),
      cex.axis = READABLE_AXIS_CEX,
      cex.lab = READABLE_LAB_CEX,
      cex.main = READABLE_MAIN_CEX
    )
    plot(fit$Power, fit$signed_R2,
      xlab = "Soft threshold power",
      ylab = "Scale-free topology fit signed R^2",
      type = "b",
      main = paste(dataset_label, "scale independence")
    )
    text(fit$Power, fit$signed_R2, labels = fit$Power, pos = 3, cex = READABLE_POINT_LABEL_CEX, col = "red")
    abline(h = r2_cutoff, col = "red", lty = 2)
    abline(v = selected_power, col = "darkgreen", lty = 2)

    plot(fit$Power, fit$mean.k.,
      xlab = "Soft threshold power",
      ylab = "Mean connectivity",
      type = "b",
      main = paste(dataset_label, "mean connectivity")
    )
    text(fit$Power, fit$mean.k., labels = fit$Power, pos = 3, cex = READABLE_POINT_LABEL_CEX, col = "red")
    abline(v = selected_power, col = "darkgreen", lty = 2)
    dev.off()
  }

  data.table(
    dataset = dataset_label,
    selected_power = selected_power,
    selection_rule = rule,
    r2_cutoff = r2_cutoff,
    selected_signed_R2 = fit$signed_R2[selected_idx],
    selected_mean_connectivity = fit$mean.k.[selected_idx],
    cor_method = cor_method
  )
}

run_blockwise_modules <- function(dat_expr, power, cor_method, dataset_label) {
  message("[", dataset_label, "] blockwiseModules power=", power, " cor=", cor_method)
  WGCNA::blockwiseModules(
    dat_expr,
    maxBlockSize = MAX_BLOCK_SIZE,
    power = power,
    networkType = NETWORK_TYPE,
    TOMType = TOM_TYPE,
    corType = cor_method,
    maxPOutliers = if (cor_method == "bicor") 0.1 else 1,
    pearsonFallback = "individual",
    deepSplit = DEEP_SPLIT,
    minModuleSize = MIN_MODULE_SIZE,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    reassignThreshold = 0,
    pamRespectsDendro = PAM_RESPECTS_DENDRO,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    nThreads = N_THREADS,
    randomSeed = params$wgcna_seed,
    verbose = 3
  )
}

module_eigengenes <- function(dat_expr, module_colors) {
  MEs0 <- WGCNA::moduleEigengenes(dat_expr, module_colors)$eigengenes
  WGCNA::orderMEs(MEs0)
}

plot_module_dendrogram <- function(net, module_colors, out_dir, file_prefix, title) {
  pdf(file.path(out_dir, paste0(file_prefix, "_module_dendrogram.pdf")), width = 12, height = 8)
  on.exit(dev.off(), add = TRUE)
  for (b in seq_along(net$dendrograms)) {
    genes <- net$blockGenes[[b]]
    WGCNA::plotDendroAndColors(
      net$dendrograms[[b]],
      module_colors[genes],
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = paste(title, "- block", b)
    )
  }
}

plot_eigengene_network <- function(MEs, out_dir, file_prefix, title) {
  pdf(file.path(out_dir, paste0(file_prefix, "_module_eigengene_network.pdf")), width = 12, height = 8)
  on.exit(dev.off(), add = TRUE)
  WGCNA::plotEigengeneNetworks(
    MEs,
    title,
    marDendro = c(0, 4, 1, 2),
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = TRUE,
    xLabelsAngle = 90
  )
}

module_trait_correlation <- function(MEs, traits, trait_cols, out_dir, file_prefix, title) {
  trait_cols <- intersect(trait_cols, colnames(traits))
  traits_aligned <- traits[rownames(MEs), trait_cols, drop = FALSE]
  for (col in colnames(traits_aligned)) traits_aligned[[col]] <- as.numeric(traits_aligned[[col]])

  complete_cols <- colSums(!is.na(traits_aligned)) > 2
  traits_aligned <- traits_aligned[, complete_cols, drop = FALSE]
  if (ncol(traits_aligned) == 0) {
    return(NULL)
  }

  cor_mat <- cor(MEs, traits_aligned, use = "pairwise.complete.obs")
  n_mat <- matrix(
    NA_integer_,
    nrow = ncol(MEs), ncol = ncol(traits_aligned),
    dimnames = list(colnames(MEs), colnames(traits_aligned))
  )
  for (i in seq_len(ncol(MEs))) {
    for (j in seq_len(ncol(traits_aligned))) {
      n_mat[i, j] <- sum(complete.cases(MEs[, i], traits_aligned[, j]))
    }
  }
  p_mat <- WGCNA::corPvalueStudent(cor_mat, n_mat)
  write_matrix(cor_mat, file.path(out_dir, paste0(file_prefix, "_module_trait_cor.csv")))
  write_matrix(n_mat, file.path(out_dir, paste0(file_prefix, "_module_trait_n.csv")))
  write_matrix(p_mat, file.path(out_dir, paste0(file_prefix, "_module_trait_pvalue.csv")))

  tidy <- as.data.table(as.table(cor_mat))
  setnames(tidy, c("module", "trait", "correlation"))
  tidy[, module := sub("^ME", "", module)]
  p_tidy <- as.data.table(as.table(p_mat))
  setnames(p_tidy, c("module", "trait", "p_value"))
  p_tidy[, module := sub("^ME", "", module)]
  n_tidy <- as.data.table(as.table(n_mat))
  setnames(n_tidy, c("module", "trait", "n_pairwise"))
  n_tidy[, module := sub("^ME", "", module)]
  tidy <- merge(tidy, p_tidy, by = c("module", "trait"))
  tidy <- merge(tidy, n_tidy, by = c("module", "trait"))
  tidy[, fdr := p.adjust(p_value, method = "BH"), by = trait]
  write_dt(tidy, file.path(out_dir, paste0(file_prefix, "_module_trait_tidy.csv")))

  text_mat <- paste0(signif(cor_mat, 2), "\n(", signif(p_mat, 1), ")")
  dim(text_mat) <- dim(cor_mat)
  if (identical(file_prefix, "copd")) {
    heatmap_width <- 11
    heatmap_height <- max(8.5, nrow(cor_mat) * 0.45 + 4.5)
    x_labels <- unname(COPD_MODULE_TRAIT_HEATMAP_LABELS[colnames(cor_mat)])
    x_labels[is.na(x_labels)] <- colnames(cor_mat)[is.na(x_labels)]
  } else {
    heatmap_width <- max(8, ncol(cor_mat) * 1.6 + 4)
    heatmap_height <- max(7, nrow(cor_mat) * 0.32 + 4)
    x_labels <- colnames(cor_mat)
  }
  draw_heatmap <- function() {
    par(mar = c(7.5, 11, 4.2, 3), cex.main = READABLE_MAIN_CEX)
    WGCNA::labeledHeatmap(
      Matrix = cor_mat,
      xLabels = x_labels,
      yLabels = rownames(cor_mat),
      ySymbols = rownames(cor_mat),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = text_mat,
      setStdMargins = FALSE,
      cex.text = READABLE_HEATMAP_TEXT_CEX,
      cex.lab = READABLE_HEATMAP_LABEL_CEX,
      cex.legendLabel = READABLE_HEATMAP_LABEL_CEX,
      zlim = c(-1, 1),
      main = title
    )
  }
  pdf(file.path(out_dir, paste0(file_prefix, "_module_trait_heatmap.pdf")),
    width = heatmap_width, height = heatmap_height
  )
  draw_heatmap()
  dev.off()
  png(file.path(out_dir, paste0(file_prefix, "_module_trait_heatmap.png")),
    width = heatmap_width * 240, height = heatmap_height * 240, res = 240
  )
  draw_heatmap()
  dev.off()
  list(cor = cor_mat, n = n_mat, p = p_mat, tidy = tidy)
}

fit_module_trait_models <- function(MEs, traits, model_specs, out_dir, file_prefix) {
  rows <- list()
  for (spec in model_specs) {
    vars <- all.vars(spec$formula)
    rhs_vars <- setdiff(vars, "ME")
    trait_aligned <- traits[rownames(MEs), rhs_vars, drop = FALSE]
    for (me_name in colnames(MEs)) {
      module <- sub("^ME", "", me_name)
      df <- data.frame(ME = as.numeric(MEs[, me_name]), trait_aligned, check.names = FALSE)
      cc <- complete.cases(df[, vars, drop = FALSE])
      if (sum(cc) <= length(rhs_vars) + 2) next
      fit <- lm(spec$formula, data = df[cc, , drop = FALSE])
      coefs <- summary(fit)$coefficients
      rows[[length(rows) + 1L]] <- data.table(
        module = module,
        model_name = spec$name,
        trait_name = spec$trait_name,
        primary_term = spec$primary_term,
        term = rownames(coefs),
        beta = coefs[, "Estimate"],
        std_error = coefs[, "Std. Error"],
        statistic = coefs[, "t value"],
        p_value = coefs[, "Pr(>|t|)"],
        n_model = sum(cc)
      )
    }
  }
  if (length(rows) == 0) {
    return(data.table())
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, fdr := p.adjust(p_value, method = "BH"), by = .(model_name, term)]
  write_dt(out, file.path(out_dir, paste0(file_prefix, "_module_trait_models.csv")))
  out
}

compute_gene_metrics <- function(dat_expr, MEs, traits, module_colors, trait_cols, out_dir, file_prefix,
                                 gene_map = NULL, gene_symbol_col = "gene_symbol") {
  if (COR_METHOD == "bicor") {
    kp <- WGCNA::bicorAndPvalue(dat_expr, MEs, use = "pairwise.complete.obs", maxPOutliers = 0.1)
    k_cor <- if (!is.null(kp$bicor)) kp$bicor else kp$cor
    k_p <- kp$p
  } else {
    kp <- WGCNA::corAndPvalue(dat_expr, MEs, use = "pairwise.complete.obs")
    k_cor <- kp$cor
    k_p <- kp$p
  }

  kme <- as.data.table(k_cor, keep.rownames = "gene_id")
  setnames(kme, setdiff(names(kme), "gene_id"), paste0("kME_", setdiff(names(kme), "gene_id")))
  kme_p <- as.data.table(k_p, keep.rownames = "gene_id")
  setnames(kme_p, setdiff(names(kme_p), "gene_id"), paste0("p_", setdiff(names(kme_p), "gene_id")))
  kme <- merge(kme, kme_p, by = "gene_id")

  gene_dt <- data.table(gene_id = colnames(dat_expr), module = module_colors)
  kme <- merge(gene_dt, kme, by = "gene_id", all.x = TRUE, sort = FALSE)

  trait_cols <- intersect(trait_cols, colnames(traits))
  traits_aligned <- traits[rownames(dat_expr), trait_cols, drop = FALSE]
  for (col in colnames(traits_aligned)) traits_aligned[[col]] <- as.numeric(traits_aligned[[col]])
  traits_aligned <- traits_aligned[, colSums(!is.na(traits_aligned)) > 2, drop = FALSE]

  if (ncol(traits_aligned) > 0) {
    gp <- WGCNA::corAndPvalue(dat_expr, traits_aligned, use = "pairwise.complete.obs")
    gs <- as.data.table(gp$cor, keep.rownames = "gene_id")
    setnames(gs, setdiff(names(gs), "gene_id"), paste0("GS_", setdiff(names(gs), "gene_id")))
    kme <- merge(kme, gs, by = "gene_id", all.x = TRUE, sort = FALSE)
  }

  assigned <- numeric(nrow(kme))
  assigned_p <- numeric(nrow(kme))
  assigned[] <- NA_real_
  assigned_p[] <- NA_real_
  for (i in seq_len(nrow(kme))) {
    me_col <- paste0("kME_ME", kme$module[i])
    p_col <- paste0("p_ME", kme$module[i])
    if (me_col %in% names(kme)) assigned[i] <- kme[[me_col]][i]
    if (p_col %in% names(kme)) assigned_p[i] <- kme[[p_col]][i]
  }
  kme[, assigned_kME := assigned]
  kme[, assigned_kME_p := assigned_p]

  gs_cols <- grep("^GS_", names(kme), value = TRUE)
  if (length(gs_cols) > 0) {
    kme[, max_abs_GS := do.call(pmax, c(lapply(.SD, function(x) abs(as.numeric(x))), na.rm = TRUE)), .SDcols = gs_cols]
  } else {
    kme[, max_abs_GS := NA_real_]
  }

  if (!is.null(gene_map)) {
    kme <- merge(gene_map, kme, by = "gene_id", all.y = TRUE, sort = FALSE)
  }
  write_dt(kme, file.path(out_dir, paste0(file_prefix, "_gene_kme_gs_table.csv")))

  hub <- kme[module != "grey" & !is.na(assigned_kME) & !is.na(max_abs_GS)]
  hub[, Hub_Level := fifelse(
    abs(assigned_kME) >= 0.8 & max_abs_GS >= 0.2, "Strict",
    fifelse(abs(assigned_kME) >= 0.6 & max_abs_GS >= 0.15, "Relaxed", "NotHub")
  )]
  hub <- hub[Hub_Level != "NotHub"]
  hub[, abs_assigned_kME := abs(assigned_kME)]
  setorder(hub, -abs_assigned_kME, -max_abs_GS)
  write_dt(hub, file.path(out_dir, paste0(file_prefix, "_hub_candidates.csv")))

  hub <- copy(hub)
  if (nrow(hub) > 0) {
    if (gene_symbol_col %in% names(hub)) {
      hub[, GeneSymbol := get(gene_symbol_col)]
    }
    if (!"GeneSymbol" %in% names(hub)) hub[, GeneSymbol := gene_id]
    hub[, Gene := GeneSymbol]
    hub[, Module := module]
    hub[, MM := assigned_kME]
    hub[, Max_GS := max_abs_GS]
    write_dt(hub, file.path(out_dir, "hub_genes.csv"))
  }
  kme
}

run_module_preservation <- function(train_expr, valid_expr, module_colors, out_dir, file_prefix, dataset_label) {
  if (SKIP_PRESERVATION) {
    message("[", dataset_label, "] Skipping modulePreservation because WGCNA_SKIP_PRESERVATION=1.")
    return(NULL)
  }
  common <- intersect(colnames(train_expr), colnames(valid_expr))
  train_expr <- train_expr[, common, drop = FALSE]
  valid_expr <- valid_expr[, common, drop = FALSE]
  colors <- module_colors[match(common, colnames(train_expr))]
  names(colors) <- common

  multi_expr <- list(
    Training = list(data = train_expr),
    Validation = list(data = valid_expr)
  )
  color_list <- list(Training = colors)

  message("[", dataset_label, "] modulePreservation nPermutations=", PRESERVATION_N_PERMUTATIONS)
  mp <- WGCNA::modulePreservation(
    multi_expr,
    color_list,
    dataIsExpr = TRUE,
    networkType = NETWORK_TYPE,
    corFnc = if (COR_METHOD == "bicor") "bicor" else "cor",
    corOptions = if (COR_METHOD == "bicor") {
      "use = 'pairwise.complete.obs', maxPOutliers = 0.1"
    } else {
      "use = 'pairwise.complete.obs'"
    },
    referenceNetworks = 1,
    nPermutations = PRESERVATION_N_PERMUTATIONS,
    randomSeed = params$wgcna_seed,
    quickCor = 0,
    verbose = 3
  )

  ref <- 1
  test <- 2
  stats_z <- cbind(
    mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
    mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE]
  )
  stats_obs <- cbind(
    mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
    mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE]
  )

  modules <- rownames(stats_z)
  module_size <- mp$preservation$observed[[ref]][[test]][, "moduleSize"]
  zsummary <- stats_z[, "Zsummary.pres"]
  median_rank <- stats_obs[, "medianRank.pres"]
  pres <- data.table(
    module = modules,
    module_size = as.numeric(module_size),
    Zsummary = as.numeric(zsummary),
    medianRank = as.numeric(median_rank),
    preservation_level = fifelse(
      zsummary > 10, "strong",
      fifelse(zsummary > 2, "moderate", "weak_or_none")
    )
  )
  setorder(pres, -Zsummary)
  write_dt(pres, file.path(out_dir, paste0(file_prefix, "_module_preservation.csv")))

  plot_df <- pres[module != "gold"]
  pdf(file.path(out_dir, paste0(file_prefix, "_module_preservation.pdf")), width = 12, height = 6)
  par(
    mfrow = c(1, 2),
    mar = c(5.8, 6.2, 4.2, 2),
    cex.axis = READABLE_AXIS_CEX,
    cex.lab = READABLE_LAB_CEX,
    cex.main = READABLE_MAIN_CEX
  )
  plot(plot_df$module_size, plot_df$Zsummary,
    col = plot_df$module, bg = plot_df$module, pch = 21, cex = 2,
    xlab = "Module size", ylab = "Zsummary",
    main = paste(dataset_label, "module preservation")
  )
  text(plot_df$module_size, plot_df$Zsummary, labels = plot_df$module, pos = 3, cex = READABLE_POINT_LABEL_CEX)
  abline(h = 2, col = "blue", lty = 2)
  abline(h = 10, col = "darkgreen", lty = 2)
  plot(plot_df$module_size, plot_df$medianRank,
    col = plot_df$module, bg = plot_df$module, pch = 21, cex = 2,
    xlab = "Module size", ylab = "Median rank",
    main = paste(dataset_label, "median rank")
  )
  text(plot_df$module_size, plot_df$medianRank, labels = plot_df$module, pos = 3, cex = READABLE_POINT_LABEL_CEX)
  dev.off()

  pdf(file.path(out_dir, paste0(file_prefix, "_module_preservation_Zsummary_barplot.pdf")),
    width = 11, height = 7
  )
  par(
    mar = c(10.5, 6.2, 4.2, 2),
    cex.axis = READABLE_AXIS_CEX,
    cex.lab = READABLE_LAB_CEX,
    cex.main = READABLE_MAIN_CEX
  )
  ord <- order(-plot_df$Zsummary)
  bp <- barplot(plot_df$Zsummary[ord],
    names.arg = plot_df$module[ord],
    col = plot_df$module[ord],
    las = 2,
    ylab = "Zsummary",
    main = paste(dataset_label, "module preservation Zsummary"),
    border = NA,
    cex.names = READABLE_AXIS_CEX
  )
  abline(h = 2, col = "blue", lty = 2)
  abline(h = 10, col = "darkgreen", lty = 2)
  text(bp, plot_df$Zsummary[ord], labels = round(plot_df$Zsummary[ord], 1), pos = 3, cex = READABLE_POINT_LABEL_CEX)
  dev.off()

  pres
}

summarize_modules <- function(module_colors, trait_tidy, model_terms, out_dir, file_prefix) {
  module_dt <- data.table(module = module_colors)
  summary <- module_dt[, .(module_size = .N), by = module]
  if (!is.null(trait_tidy) && nrow(trait_tidy) > 0) {
    best_trait <- trait_tidy[order(fdr, p_value)][, .SD[1], by = module]
    setnames(best_trait, c("trait", "correlation", "p_value", "fdr"),
      c("best_trait", "best_trait_cor", "best_trait_p", "best_trait_fdr"),
      skip_absent = TRUE
    )
    summary <- merge(summary, best_trait[, .(module, best_trait, best_trait_cor, best_trait_p, best_trait_fdr)],
      by = "module", all.x = TRUE
    )
  }
  if (!is.null(model_terms) && nrow(model_terms) > 0) {
    primary <- model_terms[term == primary_term]
    if (nrow(primary) > 0) {
      best_model <- primary[order(fdr, p_value)][, .SD[1], by = module]
      summary <- merge(summary,
        best_model[, .(module,
          best_model_name = model_name, best_model_trait = trait_name,
          best_model_beta = beta, best_model_p = p_value, best_model_fdr = fdr
        )],
        by = "module", all.x = TRUE
      )
    }
  }
  setorder(summary, -module_size)
  write_dt(summary, file.path(out_dir, paste0(file_prefix, "_module_summary.csv")))
  summary
}

make_model_specs <- function(traits, primary_terms) {
  specs <- list()
  covars <- intersect(c("age_z", "sex_male", "packyears_z"), colnames(traits))
  covars <- covars[vapply(traits[covars], function(x) sum(!is.na(x)) > 5 && sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))]
  for (term in primary_terms) {
    if (!term %in% colnames(traits)) next
    if (sum(!is.na(traits[[term]])) <= 5 || sd(as.numeric(traits[[term]]), na.rm = TRUE) == 0) next
    rhs <- c(term, setdiff(covars, term))
    specs[[length(specs) + 1L]] <- list(
      name = paste0(term, "_adjusted"),
      trait_name = term,
      primary_term = term,
      formula = as.formula(paste("ME ~", paste(rhs, collapse = " + ")))
    )
  }
  specs
}

run_copd_wgcna <- function() {
  out_dir <- COPD_WGCNA
  ensure_dir(out_dir, clean = CLEAN_OUTPUT)
  message("\n[wgcna] ===== COPD WGCNA =====")

  gse_train <- GEOquery::getGEO(filename = RAW_GSE76925, GSEMatrix = TRUE)
  gse_valid <- GEOquery::getGEO(filename = RAW_GSE47460, GSEMatrix = TRUE)
  chars_train <- parse_all_characteristics(gse_train)
  chars_valid <- parse_all_characteristics(gse_valid)

  clinical_train <- data.frame(row.names = rownames(chars_train), check.names = FALSE)
  clinical_train$sample_id <- rownames(chars_train)
  clinical_train$age <- numeric_clean(chars_train$age)
  clinical_train$age_z <- zscore(clinical_train$age)
  clinical_train$sex_male <- sex_to_male(chars_train$sex)
  clinical_train$packyears <- numeric_clean(chars_train$packyears)
  clinical_train$packyears_z <- zscore(clinical_train$packyears)
  clinical_train$copd_status <- fifelse(
    tolower(chars_train$copd) == "case", 1,
    fifelse(tolower(chars_train$copd) == "cont", 0, NA_real_)
  )
  clinical_train$FEV1_pp <- numeric_clean(chars_train$fev1.pp)
  clinical_train$FEV1_FVC <- numeric_clean(chars_train$fev1fvc)
  clinical_train$low_fev1_z <- -zscore(clinical_train$FEV1_pp)
  clinical_train$low_fev1_fvc_z <- -zscore(clinical_train$FEV1_FVC)
  clinical_train$laa950 <- numeric_clean(chars_train$laa950)
  clinical_train$emphysema_z <- zscore(clinical_train$laa950)
  clinical_train$perc15 <- numeric_clean(chars_train$perc15)
  clinical_train$pi10 <- numeric_clean(chars_train$pi10)
  clinical_train$GOLD <- with(clinical_train, fifelse(
    is.na(FEV1_pp) | is.na(FEV1_FVC), NA_real_,
    fifelse(
      FEV1_FVC >= 0.7, 0,
      fifelse(
        FEV1_pp >= 80, 1,
        fifelse(FEV1_pp >= 50, 2, fifelse(FEV1_pp >= 30, 3, 4))
      )
    )
  ))
  clinical_train$gold_z <- zscore(clinical_train$GOLD)

  clinical_valid <- data.frame(row.names = rownames(chars_valid), check.names = FALSE)
  clinical_valid$sample_id <- rownames(chars_valid)
  clinical_valid$disease_state <- chars_valid[["disease state"]]
  keep_valid <- clinical_valid$disease_state %in% c("Control", "Chronic Obstructive Lung Disease", "COPD")
  keep_valid[is.na(keep_valid)] <- FALSE

  expr_train <- normalize_microarray_matrix(Biobase::exprs(gse_train))
  expr_valid <- normalize_microarray_matrix(Biobase::exprs(gse_valid)[, keep_valid, drop = FALSE])
  clinical_valid <- clinical_valid[keep_valid, , drop = FALSE]

  train_map <- collapse_probes_to_genes(expr_train, Biobase::fData(gse_train)$Symbol,
    gene_symbol = Biobase::fData(gse_train)$Symbol,
    select_by = "variance"
  )
  valid_map <- collapse_probes_to_genes(expr_valid, Biobase::fData(gse_valid)$GENE_SYMBOL,
    gene_symbol = Biobase::fData(gse_valid)$GENE_SYMBOL,
    select_by = "variance"
  )
  write_dt(train_map$map, file.path(out_dir, "training_probe_to_gene_map.csv"))
  write_dt(valid_map$map, file.path(out_dir, "validation_probe_to_gene_map.csv"))

  filt <- filter_common_variable_genes(train_map$expr, valid_map$expr, out_dir)
  dat_expr <- filt$train
  valid_expr <- filt$valid
  clinical_train <- clinical_train[rownames(dat_expr), , drop = FALSE]

  qc <- apply_wgcna_sample_qc(dat_expr, valid_expr, clinical_train, out_dir, "COPD")
  dat_expr <- qc$train
  valid_expr <- qc$valid
  clinical_train <- qc$traits
  write_dt(clinical_train, file.path(out_dir, "sample_table_used.csv"))

  plot_sample_dendrogram(dat_expr, clinical_train, out_dir, "copd", "COPD training sample dendrogram")
  plot_sample_dendrogram_with_traits(
    dat_expr, clinical_train,
    c("copd_status", "GOLD", "FEV1_pp", "FEV1_FVC", "laa950", "age", "sex_male", "packyears"),
    out_dir, "copd", "COPD training sample dendrogram and traits"
  )

  soft <- select_soft_power(dat_expr, out_dir, "copd", "COPD",
    r2_cutoff = COPD_SOFT_R2_CUTOFF
  )
  if (STOP_AFTER_SOFT) {
    write_dt(soft, file.path(out_dir, "soft_power_selected.csv"))
    message("[COPD] Stopping after soft-threshold selection.")
    return(invisible(soft))
  }

  net <- run_blockwise_modules(dat_expr, soft$selected_power, COR_METHOD, "COPD")
  module_colors <- net$colors
  names(module_colors) <- colnames(dat_expr)
  MEs <- module_eigengenes(dat_expr, module_colors)

  write_dt(
    data.table(gene_id = names(module_colors), gene_symbol = names(module_colors), module = module_colors),
    file.path(out_dir, "module_assignment.csv")
  )
  write_dt(as.data.table(MEs, keep.rownames = "sample_id"), file.path(out_dir, "module_eigengenes.csv"))
  write_dt(soft, file.path(out_dir, "soft_power_selected.csv"))

  plot_module_dendrogram(net, module_colors, out_dir, "copd", "COPD gene dendrogram")
  plot_eigengene_network(MEs, out_dir, "copd", "COPD module eigengene network")

  trait_cols <- COPD_MODULE_TRAIT_HEATMAP_TRAITS
  trait_res <- module_trait_correlation(
    MEs, clinical_train, trait_cols, out_dir, "copd",
    "COPD module-trait relationships (pairwise-complete n)"
  )
  model_specs <- make_model_specs(clinical_train, c("copd_status", "gold_z", "low_fev1_z", "low_fev1_fvc_z", "emphysema_z"))
  model_terms <- fit_module_trait_models(MEs, clinical_train, model_specs, out_dir, "copd")

  compute_gene_metrics(dat_expr, MEs, clinical_train, module_colors,
    c("copd_status", "gold_z", "low_fev1_z", "low_fev1_fvc_z", "emphysema_z"),
    out_dir, "copd",
    gene_map = data.table(gene_id = colnames(dat_expr), gene_symbol = colnames(dat_expr)),
    gene_symbol_col = "gene_symbol"
  )

  pres <- run_module_preservation(dat_expr, valid_expr, module_colors, out_dir, "copd", "COPD")
  module_summary <- summarize_modules(module_colors, trait_res$tidy, model_terms, out_dir, "copd")

  save(dat_expr, valid_expr, module_colors, MEs, clinical_train, soft, pres,
    file = file.path(out_dir, "copd_wgcna_env.RData")
  )
  sink(file.path(out_dir, "session_info.txt"))
  print(sessionInfo())
  sink()
  invisible(list(soft = soft, module_summary = module_summary, preservation = pres))
}

remove_age_outliers_by_group <- function(ages, groups, sample_ids) {
  keep <- rep(TRUE, length(ages))
  rows <- list()
  removed <- list()
  for (grp in sort(unique(groups))) {
    idx <- which(groups == grp)
    grp_ages <- ages[idx]
    q1 <- as.numeric(quantile(grp_ages, 0.25, na.rm = TRUE))
    q3 <- as.numeric(quantile(grp_ages, 0.75, na.rm = TRUE))
    iqr <- q3 - q1
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr
    out <- grp_ages < lower | grp_ages > upper
    keep[idx[out]] <- FALSE
    label <- ifelse(grp == 1, "AD", "Control")
    rows[[length(rows) + 1L]] <- data.table(
      group = label,
      n_before = length(idx),
      Q1 = q1,
      Q3 = q3,
      IQR = iqr,
      lower_fence = lower,
      upper_fence = upper,
      n_removed = sum(out),
      removed_samples = paste(sample_ids[idx[out]], collapse = ";"),
      removed_ages = paste(grp_ages[out], collapse = ";")
    )
    removed[[label]] <- sample_ids[idx[out]]
  }
  list(keep = keep, summary = rbindlist(rows), removed = removed)
}

run_ad_wgcna <- function() {
  out_dir <- if (AD_SKIP_RESIDUALIZATION) {
    file.path(RESULTS, "AD_WGCNA_unadjusted")
  } else {
    AD_WGCNA
  }
  ensure_dir(out_dir, clean = CLEAN_OUTPUT)
  message("\n[wgcna] ===== AD WGCNA =====")

  gse_train <- GEOquery::getGEO(filename = RAW_GSE33000, GSEMatrix = TRUE)
  gse_valid <- GEOquery::getGEO(filename = RAW_GSE84422, GSEMatrix = TRUE)

  pheno_train <- Biobase::pData(gse_train)
  keep_no_hd <- pheno_train[["disease status:ch2"]] != "Huntington's disease"
  pheno_train <- pheno_train[keep_no_hd, , drop = FALSE]
  expr_train <- Biobase::exprs(gse_train)[, rownames(pheno_train), drop = FALSE]

  age <- numeric_clean(pheno_train[["age:ch2"]])
  sex <- sex_to_male(pheno_train[["gender:ch2"]])
  disease <- fifelse(grepl("Alzheimer", pheno_train[["disease status:ch2"]], ignore.case = TRUE), 1, 0)
  age_out <- remove_age_outliers_by_group(age, disease, rownames(pheno_train))
  write_dt(age_out$summary, file.path(out_dir, "ad_age_outlier_removal_summary.csv"))

  pheno_train <- pheno_train[age_out$keep, , drop = FALSE]
  expr_train <- expr_train[, rownames(pheno_train), drop = FALSE]
  age <- age[age_out$keep]
  sex <- sex[age_out$keep]
  disease <- disease[age_out$keep]

  clinical_train <- data.frame(row.names = rownames(pheno_train), check.names = FALSE)
  clinical_train$sample_id <- rownames(pheno_train)
  clinical_train$Disease_State <- disease
  clinical_train$Age <- age
  clinical_train$Sex <- sex
  clinical_train$age_z <- zscore(age)
  clinical_train$sex_male <- sex

  chars_valid <- parse_all_characteristics(gse_valid)
  pheno_valid <- Biobase::pData(gse_valid)
  keep_valid <- pheno_valid$platform_id == "GPL96" &
    chars_valid[["brain region"]] == "Dorsolateral Prefrontal Cortex" &
    grepl("AD|Normal", chars_valid[["neuropathological category"]], ignore.case = TRUE)
  keep_valid[is.na(keep_valid)] <- FALSE
  expr_valid <- Biobase::exprs(gse_valid)[, keep_valid, drop = FALSE]
  pheno_valid <- pheno_valid[keep_valid, , drop = FALSE]
  chars_valid <- chars_valid[keep_valid, , drop = FALSE]

  valid_trait <- data.frame(row.names = rownames(pheno_valid), check.names = FALSE)
  valid_trait$sample_id <- rownames(pheno_valid)
  valid_trait$Age <- numeric_clean(chars_valid$age)
  valid_trait$Sex <- sex_to_male(chars_valid$sex)
  valid_trait$Disease_State <- fifelse(
    grepl("AD", chars_valid[["neuropathological category"]], ignore.case = TRUE), 1,
    fifelse(chars_valid[["neuropathological category"]] == "Normal", 0, NA_real_)
  )
  write_dt(valid_trait, file.path(out_dir, "validation_sample_table_used.csv"))

  expr_train <- normalize_microarray_matrix(expr_train)
  expr_valid <- normalize_microarray_matrix(expr_valid)

  train_map <- collapse_probes_to_genes(expr_train, Biobase::fData(gse_train)$EntrezGeneID,
    gene_symbol = Biobase::fData(gse_train)$EntrezGeneID,
    select_by = "variance"
  )
  valid_map <- collapse_probes_to_genes(expr_valid, Biobase::fData(gse_valid)$ENTREZ_GENE_ID,
    gene_symbol = Biobase::fData(gse_valid)$`Gene Symbol`,
    select_by = "variance"
  )
  symbol_map <- unique(valid_map$map[!is.na(gene_symbol) & nzchar(gene_symbol), .(gene_id, gene_symbol)])
  train_map_base <- train_map$map[, setdiff(names(train_map$map), "gene_symbol"), with = FALSE]
  train_map$map <- merge(train_map_base, symbol_map, by = "gene_id", all.x = TRUE, sort = FALSE)
  train_map$map[is.na(gene_symbol) | !nzchar(gene_symbol), gene_symbol := gene_id]
  write_dt(train_map$map, file.path(out_dir, "training_probe_to_gene_map.csv"))
  write_dt(valid_map$map, file.path(out_dir, "validation_probe_to_gene_map.csv"))

  expr_train_gene <- train_map$expr
  expr_valid_gene <- valid_map$expr

  if (AD_SKIP_RESIDUALIZATION) {
    residualization_status <- data.table(
      residualization = "skipped",
      details = "Age outliers were removed within disease groups, but age/sex residualization was skipped."
    )
  } else {
    clin_aligned <- clinical_train[colnames(expr_train_gene), , drop = FALSE]
    covariates <- cbind(Age = clin_aligned$Age, Sex = clin_aligned$Sex)
    design_protect <- model.matrix(~Disease_State, data = clin_aligned)
    expr_train_gene <- limma::removeBatchEffect(expr_train_gene,
      covariates = covariates,
      design = design_protect
    )
    residualization_status <- data.table(
      residualization = "limma_removeBatchEffect_age_sex",
      details = "Age and Sex were adjusted with limma::removeBatchEffect while protecting Disease_State in the design matrix."
    )
  }
  write_dt(residualization_status, file.path(out_dir, "ad_residualization_status.csv"))

  filt <- filter_common_variable_genes(expr_train_gene, expr_valid_gene, out_dir)
  dat_expr <- filt$train
  valid_expr <- filt$valid
  clinical_train <- clinical_train[rownames(dat_expr), , drop = FALSE]

  qc <- apply_wgcna_sample_qc(dat_expr, valid_expr, clinical_train, out_dir, "AD")
  dat_expr <- qc$train
  valid_expr <- qc$valid
  clinical_train <- qc$traits
  write_dt(clinical_train, file.path(out_dir, "sample_table_used.csv"))

  plot_sample_dendrogram(dat_expr, clinical_train, out_dir, "ad", "AD training sample dendrogram")
  plot_sample_dendrogram_with_traits(
    dat_expr, clinical_train,
    c("Disease_State", "Age", "Sex"),
    out_dir, "ad", "AD training sample dendrogram and traits"
  )

  soft <- select_soft_power(dat_expr, out_dir, "ad", "AD",
    r2_cutoff = AD_SOFT_R2_CUTOFF
  )
  if (STOP_AFTER_SOFT) {
    write_dt(soft, file.path(out_dir, "soft_power_selected.csv"))
    message("[AD] Stopping after soft-threshold selection.")
    return(invisible(soft))
  }

  net <- run_blockwise_modules(dat_expr, soft$selected_power, COR_METHOD, "AD")
  module_colors <- net$colors
  names(module_colors) <- colnames(dat_expr)
  MEs <- module_eigengenes(dat_expr, module_colors)

  gene_map <- unique(train_map$map[, .(gene_id, gene_symbol)])
  gene_map <- gene_map[match(colnames(dat_expr), gene_id)]
  gene_map[is.na(gene_symbol) | !nzchar(gene_symbol), gene_symbol := gene_id]
  module_assign <- data.table(
    gene_id = colnames(dat_expr),
    gene_symbol = gene_map$gene_symbol,
    module = module_colors
  )
  write_dt(module_assign, file.path(out_dir, "module_assignment.csv"))
  write_dt(as.data.table(MEs, keep.rownames = "sample_id"), file.path(out_dir, "module_eigengenes.csv"))
  write_dt(soft, file.path(out_dir, "soft_power_selected.csv"))

  plot_module_dendrogram(net, module_colors, out_dir, "ad", "AD gene dendrogram")
  plot_eigengene_network(MEs, out_dir, "ad", "AD module eigengene network")

  trait_res <- module_trait_correlation(
    MEs, clinical_train, c("Disease_State", "Age", "Sex"), out_dir, "ad",
    paste0("AD module-trait relationships (raw expression; n=", nrow(MEs), ")")
  )
  model_specs <- make_model_specs(clinical_train, c("Disease_State"))
  model_terms <- fit_module_trait_models(MEs, clinical_train, model_specs, out_dir, "ad")

  compute_gene_metrics(dat_expr, MEs, clinical_train, module_colors,
    c("Disease_State", "Age", "Sex"),
    out_dir, "ad",
    gene_map = gene_map,
    gene_symbol_col = "gene_symbol"
  )

  pres <- run_module_preservation(dat_expr, valid_expr, module_colors, out_dir, "ad", "AD")
  module_summary <- summarize_modules(module_colors, trait_res$tidy, model_terms, out_dir, "ad")

  save(dat_expr, valid_expr, module_colors, MEs, clinical_train, soft, pres,
    file = file.path(out_dir, "ad_wgcna_env.RData")
  )
  sink(file.path(out_dir, "session_info.txt"))
  print(sessionInfo())
  sink()
  invisible(list(soft = soft, module_summary = module_summary, preservation = pres))
}

requested <- toupper(Sys.getenv("WGCNA_DATASET", "both"))
if (!requested %in% c("BOTH", "COPD", "AD")) {
  stop("WGCNA_DATASET must be one of: both, COPD, AD.")
}

if (requested %in% c("BOTH", "COPD")) {
  run_copd_wgcna()
}
if (requested %in% c("BOTH", "AD")) {
  run_ad_wgcna()
}

message("[wgcna] Done.")
