source("R/runtime.R")

source("R/paths.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(GEOquery)
  library(WGCNA)
})

options(stringsAsFactors = FALSE, warn = 1)

N_PERMUTATIONS <- as.integer(Sys.getenv("AD_GSE84422_PRESERVATION_N_PERMUTATIONS", "1000"))
N_THREADS <- as.integer(Sys.getenv("AD_GSE84422_PRESERVATION_THREADS", "2"))
MAX_MODULE_SIZE <- as.integer(Sys.getenv("AD_GSE84422_PRESERVATION_MAX_MODULE_SIZE", "1000"))
MAX_GOLD_MODULE_SIZE <- as.integer(Sys.getenv("AD_GSE84422_PRESERVATION_MAX_GOLD_MODULE_SIZE", "1000"))
RANDOM_SEED <- params$brain_preservation_seed

if (is.na(N_PERMUTATIONS) || N_PERMUTATIONS < 20L) stop("nPermutations must be >= 20")
if (is.na(N_THREADS) || N_THREADS < 1L) stop("threads must be >= 1")
if (is.na(MAX_MODULE_SIZE) || MAX_MODULE_SIZE < 30L) stop("maxModuleSize must be >= 30")

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores) || detected_cores < 1L) detected_cores <- 1L
N_THREADS <- min(N_THREADS, detected_cores)

OUT_TAG <- if (N_PERMUTATIONS == 1000L) "1000" else as.character(N_PERMUTATIONS)
OUT_DIR <- file.path(
  RESULTS,
  paste0("AD_GSE33000_to_GSE84422_DLPFC_module_preservation_", OUT_TAG)
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DIR <- file.path(RESULTS, "AD_WGCNA_unadjusted")
ENV_RDATA <- file.path(DIR, "ad_wgcna_env.RData")
MODULE_ASSIGNMENT <- file.path(DIR, "module_assignment.csv")
REFERENCE_SAMPLE_TABLE <- file.path(DIR, "sample_table_used.csv")
VALIDATION_SAMPLE_TABLE <- file.path(DIR, "validation_sample_table_used.csv")
WORKSPACE_ROOT <- normalizePath(getwd())
RAW_GSE33000 <- file.path(WORKSPACE_ROOT, "data", "raw", "bulk_brain_ad", "GSE33000_series_matrix.txt.gz")
RAW_GSE84422 <- file.path(WORKSPACE_ROOT, "data", "raw", "bulk_brain_ad", "GSE84422-GPL96_series_matrix.txt.gz")

required_paths <- c(
  ENV_RDATA, MODULE_ASSIGNMENT, REFERENCE_SAMPLE_TABLE,
  VALIDATION_SAMPLE_TABLE, RAW_GSE33000, RAW_GSE84422
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))

write_dt <- function(x, filename) fwrite(as.data.table(x), filename)

find_col <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (!length(hit)) {
    return(NA_character_)
  }
  hit[[1]]
}

numeric_clean <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))

age_outlier_flag <- function(age, group) {
  flag <- rep(FALSE, length(age))
  for (g in sort(unique(group[!is.na(group)]))) {
    idx <- which(group == g & !is.na(age))
    if (length(idx) < 4L) next
    q <- stats::quantile(age[idx], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    iqr <- q[[2]] - q[[1]]
    flag[idx] <- age[idx] < (q[[1]] - 1.5 * iqr) | age[idx] > (q[[2]] + 1.5 * iqr)
  }
  flag
}

sample_connectivity_z <- function(expr) {
  sample_cor <- stats::cor(t(expr), use = "pairwise.complete.obs", method = "pearson")
  sample_adj <- (1 + sample_cor) / 2
  diag(sample_adj) <- 0
  as.numeric(scale(rowSums(sample_adj, na.rm = TRUE)))
}

extract_preservation <- function(mp, colors) {
  ref <- 1L
  test <- 2L
  quality_z_all <- mp$quality$Z[[ref]][[test]]
  preservation_z_all <- mp$preservation$Z[[ref]][[test]]
  quality_observed_all <- mp$quality$observed[[ref]][[test]]
  preservation_observed_all <- mp$preservation$observed[[ref]][[test]]
  modules <- rownames(preservation_z_all)
  quality_z <- quality_z_all[modules, -1, drop = FALSE]
  preservation_z <- preservation_z_all[modules, -1, drop = FALSE]
  quality_observed <- quality_observed_all[modules, -1, drop = FALSE]
  preservation_observed <- preservation_observed_all[modules, -1, drop = FALSE]
  z <- cbind(quality_z[modules, , drop = FALSE], preservation_z[modules, , drop = FALSE])
  observed <- cbind(quality_observed[modules, , drop = FALSE], preservation_observed[modules, , drop = FALSE])
  colnames(z) <- make.unique(paste0("Z_", colnames(z)), sep = "_")
  colnames(observed) <- make.unique(paste0("observed_", colnames(observed)), sep = "_")
  z_dt <- as.data.table(z)
  observed_dt <- as.data.table(observed)
  evaluated_sizes <- as.numeric(preservation_observed_all[modules, "moduleSize"])
  zsummary_values <- as.numeric(preservation_z_all[modules, "Zsummary.pres"])
  median_rank_values <- as.numeric(preservation_observed_all[modules, "medianRank.pres"])
  out <- cbind(data.table(
    module = modules,
    evaluated_module_size = evaluated_sizes,
    Zsummary = zsummary_values,
    medianRank = median_rank_values
  ), z_dt, observed_dt)
  full_sizes <- as.data.table(table(colors))
  setnames(full_sizes, c("module", "full_module_size"))
  full_sizes[, module := as.character(module)]
  out <- merge(out, full_sizes, by = "module", all.x = TRUE, sort = FALSE)
  out[, preservation_level := fifelse(
    Zsummary > 10, "strong",
    fifelse(Zsummary > 2, "moderate", "weak_or_none")
  )]
  out[, axis_role := fifelse(
    module == "blue", "Axis A source module (AD blue)",
    fifelse(module == "turquoise", "Axis B source module (AD turquoise)", "all-module result")
  )]
  setcolorder(out, c(
    "module", "axis_role", "full_module_size", "evaluated_module_size",
    "Zsummary", "medianRank", "preservation_level",
    setdiff(names(out), c(
      "module", "axis_role", "full_module_size", "evaluated_module_size",
      "Zsummary", "medianRank", "preservation_level"
    ))
  ))
  setorder(out, -Zsummary)
  out
}

save_plots <- function(summary_dt) {
  plot_dt <- copy(summary_dt[!module %chin% c("grey", "gold")])
  setorder(plot_dt, -Zsummary)
  if (!nrow(plot_dt)) {
    return(invisible(NULL))
  }

  png(file.path(OUT_DIR, "all_modules_preservation_zsummary_median_rank.png"),
    width = 3000, height = 1600, res = 220
  )
  par(mfrow = c(1, 2), mar = c(5.5, 5.5, 3.5, 1.5))
  plot(log10(plot_dt$full_module_size), plot_dt$Zsummary,
    pch = 21, bg = plot_dt$module, col = "black", cex = 1.7,
    xlab = "log10(full reference module size)", ylab = "Zsummary.pres",
    main = "AD modules: preservation"
  )
  text(log10(plot_dt$full_module_size), plot_dt$Zsummary,
    labels = plot_dt$module, pos = 3, cex = 0.8
  )
  abline(h = c(2, 10), lty = 2, col = c("#3B82F6", "#15803D"))
  plot(log10(plot_dt$full_module_size), plot_dt$medianRank,
    pch = 21, bg = plot_dt$module, col = "black", cex = 1.7,
    xlab = "log10(full reference module size)", ylab = "medianRank.pres",
    main = "AD modules: median rank (lower is better)"
  )
  text(log10(plot_dt$full_module_size), plot_dt$medianRank,
    labels = plot_dt$module, pos = 3, cex = 0.8
  )
  dev.off()

  png(file.path(OUT_DIR, "all_modules_zsummary_barplot.png"),
    width = 2600, height = 1800, res = 220
  )
  par(mar = c(10, 5.5, 3.5, 1.5))
  bar_x <- barplot(plot_dt$Zsummary,
    names.arg = plot_dt$module, col = plot_dt$module,
    border = ifelse(plot_dt$module %chin% c("blue", "turquoise"), "black", NA),
    las = 2, ylab = "Zsummary.pres",
    main = "GSE33000 modules preserved in GSE84422 DLPFC"
  )
  abline(h = c(2, 10), lty = 2, col = c("#3B82F6", "#15803D"))
  text(bar_x, plot_dt$Zsummary, labels = sprintf("%.1f", plot_dt$Zsummary), pos = 3, cex = 0.7)
  dev.off()
}

run_env <- new.env(parent = emptyenv())
load(ENV_RDATA, envir = run_env)
required_objects <- c("dat_expr", "valid_expr", "module_colors", "clinical_train")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = run_env, inherits = FALSE)]
if (length(missing_objects)) stop("Network environment lacks: ", paste(missing_objects, collapse = ", "))

reference_expr <- as.matrix(run_env$dat_expr)
validation_expr <- as.matrix(run_env$valid_expr)
colors <- as.character(run_env$module_colors)
if (is.null(names(colors))) names(colors) <- colnames(reference_expr)

if (!identical(colnames(reference_expr), colnames(validation_expr))) stop("Reference/test gene order is not identical")
if (!identical(names(colors), colnames(reference_expr))) stop("module_colors are not aligned to reference genes")
if (anyDuplicated(colnames(reference_expr))) stop("Duplicate reference gene IDs")
if (any(!is.finite(reference_expr)) || any(!is.finite(validation_expr))) stop("Non-finite expression values")

assignment <- fread(MODULE_ASSIGNMENT)
assignment[, gene_id := as.character(gene_id)]
assignment <- assignment[match(colnames(reference_expr), gene_id)]
if (anyNA(assignment$gene_id) || !identical(as.character(assignment$module), unname(colors))) {
  stop("module_assignment.csv does not match the saved module colors")
}

reference_samples <- fread(REFERENCE_SAMPLE_TABLE)
validation_samples <- fread(VALIDATION_SAMPLE_TABLE)
if (!identical(rownames(reference_expr), reference_samples$sample_id)) stop("Reference sample table is not aligned")
if (!identical(rownames(validation_expr), validation_samples$sample_id)) stop("Validation sample table is not aligned")

gse84422 <- GEOquery::getGEO(filename = RAW_GSE84422, GSEMatrix = TRUE)
if (is.list(gse84422)) gse84422 <- gse84422[[1]]
pheno <- Biobase::pData(gse84422)
pheno$sample_id <- rownames(pheno)
validation_meta <- pheno[match(validation_samples$sample_id, pheno$sample_id), , drop = FALSE]
region_col <- find_col(validation_meta, c("brain region:ch1", "brain region"))
pathology_col <- find_col(validation_meta, c("neuropathological category:ch1", "neuropathological category"))
subject_col <- find_col(validation_meta, c("subject id:ch1", "subject id"))
age_col <- find_col(validation_meta, c("age:ch1", "age"))
if (is.na(region_col) || any(validation_meta[[region_col]] != "Dorsolateral Prefrontal Cortex")) {
  stop("GSE84422 validation samples are not all explicitly DLPFC")
}

validation_qc <- data.table(
  sample_id = validation_samples$sample_id,
  subject_id = if (!is.na(subject_col)) as.character(validation_meta[[subject_col]]) else NA_character_,
  title = as.character(validation_meta$title),
  brain_region = as.character(validation_meta[[region_col]]),
  neuropathological_category = if (!is.na(pathology_col)) as.character(validation_meta[[pathology_col]]) else NA_character_,
  age = validation_samples$Age,
  sex_male = validation_samples$Sex,
  disease_state_broad_ad = validation_samples$Disease_State,
  age_outlier_reference_rule = age_outlier_flag(validation_samples$Age, validation_samples$Disease_State),
  connectivity_z = sample_connectivity_z(validation_expr)
)
validation_qc[, connectivity_outlier_z_lt_neg2_5 := connectivity_z < -2.5]
write_dt(validation_qc, file.path(OUT_DIR, "gse84422_dlpfc_validation_metadata_and_qc.csv"))

if (nrow(validation_qc) != 57L) stop("Expected 57 GSE84422 DLPFC samples")
if (length(unique(validation_qc$subject_id)) != 57L) stop("GSE84422 DLPFC validation samples are not one sample per subject")
if (length(intersect(rownames(reference_expr), rownames(validation_expr))) != 0L) stop("Reference/test GSM IDs overlap")

module_sizes <- as.data.table(table(colors))
setnames(module_sizes, c("module", "full_module_size"))
module_sizes[, module := as.character(module)]
write_dt(module_sizes[order(-full_module_size)], file.path(OUT_DIR, "reference_module_sizes.csv"))

stopifnot(
  identical(dim(reference_expr), c(459L, 8807L)), identical(dim(validation_expr), c(57L, 8807L)),
  all(validation_qc$brain_region == "Dorsolateral Prefrontal Cortex"),
  length(setdiff(unique(colors), "grey")) == 16L
)

cluster <- NULL
if (N_THREADS > 1L) {
  cluster <- parallel::makeCluster(N_THREADS)
  doParallel::registerDoParallel(cluster)
  on.exit(
    {
      parallel::stopCluster(cluster)
      foreach::registerDoSEQ()
    },
    add = TRUE
  )
} else {
  foreach::registerDoSEQ()
}

multi_expr <- list(
  GSE33000_reference = list(data = reference_expr),
  GSE84422_DLPFC_validation = list(data = validation_expr)
)
multi_colors <- list(GSE33000_reference = colors)

started_at <- Sys.time()
mp <- WGCNA::modulePreservation(
    multi_expr,
    multi_colors,
    dataIsExpr = TRUE,
    networkType = "signed",
    corFnc = "cor",
    corOptions = "use = 'pairwise.complete.obs'",
    referenceNetworks = 1,
    nPermutations = N_PERMUTATIONS,
    randomSeed = RANDOM_SEED,
    maxGoldModuleSize = MAX_GOLD_MODULE_SIZE,
    maxModuleSize = MAX_MODULE_SIZE,
    quickCor = 1,
    parallelCalculation = N_THREADS > 1L,
    savePermutedStatistics = TRUE,
    permutedStatisticsFile = file.path(OUT_DIR, "permutedStats.RData"),
    verbose = 3
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

full_stats <- extract_preservation(mp, colors)
summary_stats <- full_stats[, .(
  module, axis_role, full_module_size, evaluated_module_size,
  Zsummary, medianRank, preservation_level
)]
write_dt(full_stats, file.path(OUT_DIR, "module_preservation_full_statistics.csv"))
write_dt(summary_stats, file.path(OUT_DIR, "module_preservation_summary.csv"))
write_dt(summary_stats[!module %chin% c("grey", "gold")], file.path(OUT_DIR, "module_preservation_no_grey_gold.csv"))
write_dt(summary_stats[module %chin% c("blue", "turquoise")], file.path(OUT_DIR, "axisA_axisB_preservation_summary.csv"))
write_dt(
  summary_stats[!module %chin% c("grey", "gold"), .N, by = preservation_level][order(preservation_level)],
  file.path(OUT_DIR, "preservation_level_counts_no_grey_gold.csv")
)
save_plots(summary_stats)

stopifnot(all(summary_stats[!module %chin% c("grey", "gold"), is.finite(Zsummary) & is.finite(medianRank)]))

run_summary <- data.table(
  metric = c(
    "reference_dataset", "validation_dataset", "validation_region", "reference_samples", "validation_samples",
    "reference_genes", "common_genes_used", "all_modules_including_grey", "non_grey_modules",
    "n_permutations", "workers", "network_type", "correlation", "max_module_size", "elapsed_seconds"
  ),
  value = c(
    "GSE33000", "GSE84422 GPL96", "DLPFC BA46", nrow(reference_expr), nrow(validation_expr),
    ncol(reference_expr), ncol(validation_expr), length(unique(colors)), length(setdiff(unique(colors), "grey")),
    N_PERMUTATIONS, N_THREADS, "signed", "Pearson", MAX_MODULE_SIZE, sprintf("%.2f", elapsed_seconds)
  )
)
write_dt(run_summary, file.path(OUT_DIR, "preservation_run_summary.csv"))

save(mp, summary_stats, full_stats, colors, reference_expr, validation_expr,
  file = file.path(OUT_DIR, paste0("ad_gse33000_to_gse84422_module_preservation_", OUT_TAG, ".RData"))
)

capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))
message("Brain module preservation written to ", OUT_DIR)
