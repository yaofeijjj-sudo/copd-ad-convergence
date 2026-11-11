source("R/runtime.R")

source("R/paths.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(GEOquery)
  library(WGCNA)
})

options(stringsAsFactors = FALSE, warn = 1)

N_PERMUTATIONS <- as.integer(Sys.getenv("COPD_GSE47460_PRESERVATION_N_PERMUTATIONS", "1000"))
RANDOM_SEED <- params$lung_preservation_seed
MAX_MODULE_SIZE_MODE <- Sys.getenv("COPD_GSE47460_MAX_MODULE_SIZE_MODE", "full")
MAX_GOLD_MODULE_SIZE <- as.integer(Sys.getenv("COPD_GSE47460_MAX_GOLD_MODULE_SIZE", "1000"))

if (is.na(N_PERMUTATIONS) || N_PERMUTATIONS < 20L) stop("nPermutations must be >= 20")
if (!MAX_MODULE_SIZE_MODE %in% c("full", "cap1000")) stop("MAX_MODULE_SIZE_MODE must be full or cap1000")

DIR <- COPD_WGCNA
ENV_RDATA <- file.path(DIR, "copd_wgcna_env.RData")
MODULE_ASSIGNMENT <- file.path(DIR, "module_assignment.csv")
REFERENCE_SAMPLE_TABLE <- file.path(DIR, "sample_table_used.csv")
OUT_DIR <- file.path(
  LUNG_RESULTS,
  paste0(
    "COPD_GSE76925_to_GSE47460_GPL14550_module_preservation_",
    N_PERMUTATIONS, "_", MAX_MODULE_SIZE_MODE, ""
  )
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  ENV_RDATA, MODULE_ASSIGNMENT, REFERENCE_SAMPLE_TABLE,
  RAW_GSE76925, RAW_GSE47460
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))

write_dt <- function(x, filename) fwrite(as.data.table(x), filename)

parse_characteristic <- function(pheno, key) {
  char_cols <- grep("^characteristics_ch", colnames(pheno), value = TRUE)
  vapply(seq_len(nrow(pheno)), function(i) {
    values <- as.character(unlist(pheno[i, char_cols, drop = TRUE]))
    hit <- grep(paste0("^", key, ":"), values, ignore.case = TRUE)
    if (!length(hit)) {
      return(NA_character_)
    }
    trimws(sub("^[^:]+:", "", values[[hit[[1]]]]))
  }, character(1))
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

  z <- cbind(
    quality_z_all[modules, -1, drop = FALSE],
    preservation_z_all[modules, -1, drop = FALSE]
  )
  observed <- cbind(
    quality_observed_all[modules, -1, drop = FALSE],
    preservation_observed_all[modules, -1, drop = FALSE]
  )
  colnames(z) <- make.unique(paste0("Z_", colnames(z)), sep = "_")
  colnames(observed) <- make.unique(paste0("observed_", colnames(observed)), sep = "_")

  out <- cbind(data.table(
    module = modules,
    evaluated_module_size = as.numeric(preservation_observed_all[modules, "moduleSize"]),
    Zsummary = as.numeric(preservation_z_all[modules, "Zsummary.pres"]),
    medianRank = as.numeric(preservation_observed_all[modules, "medianRank.pres"])
  ), as.data.table(z), as.data.table(observed))

  full_sizes <- as.data.table(table(colors))
  setnames(full_sizes, c("module", "full_module_size"))
  full_sizes[, module := as.character(module)]
  out <- merge(out, full_sizes, by = "module", all.x = TRUE, sort = FALSE)
  out[, preservation_level := fifelse(
    Zsummary > 10, "strong",
    fifelse(Zsummary > 2, "moderate", "weak_or_none")
  )]
  setcolorder(out, c(
    "module", "full_module_size", "evaluated_module_size", "Zsummary",
    "medianRank", "preservation_level",
    setdiff(names(out), c(
      "module", "full_module_size", "evaluated_module_size", "Zsummary",
      "medianRank", "preservation_level"
    ))
  ))
  setorder(out, -Zsummary)
  out
}

save_separate_plots <- function(summary_dt) {
  plot_dt <- copy(summary_dt[!module %chin% c("grey", "gold")])
  setorder(plot_dt, -Zsummary)
  if (!nrow(plot_dt)) {
    return(invisible(NULL))
  }

  png(file.path(OUT_DIR, "COPD_modules_Zsummary_in_GSE47460.png"),
    width = 2600, height = 1800, res = 220
  )
  par(mar = c(10, 5.5, 3.5, 1.5))
  x <- barplot(plot_dt$Zsummary,
    names.arg = plot_dt$module, col = plot_dt$module, border = NA,
    las = 2, ylab = "Zsummary.pres",
    main = "COPD modules preserved in GSE47460"
  )
  abline(h = c(2, 10), lty = 2, col = c("#3B82F6", "#15803D"))
  text(x, plot_dt$Zsummary, labels = sprintf("%.1f", plot_dt$Zsummary), pos = 3, cex = 0.72)
  dev.off()

  png(file.path(OUT_DIR, "COPD_modules_medianRank_in_GSE47460.png"),
    width = 2600, height = 1800, res = 220
  )
  par(mar = c(10, 5.5, 3.5, 1.5))
  x <- barplot(plot_dt$medianRank,
    names.arg = plot_dt$module, col = plot_dt$module, border = NA,
    las = 2, ylab = "medianRank.pres (lower is better)",
    main = "Relative preservation rank of COPD modules"
  )
  text(x, plot_dt$medianRank, labels = sprintf("%.1f", plot_dt$medianRank), pos = 3, cex = 0.72)
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
MAX_MODULE_SIZE <- if (MAX_MODULE_SIZE_MODE == "full") ncol(reference_expr) else 1000L

if (!identical(colnames(reference_expr), colnames(validation_expr))) stop("Reference/test gene order is not identical")
if (!identical(names(colors), colnames(reference_expr))) stop("module_colors are not aligned to network genes")
if (anyDuplicated(rownames(reference_expr)) || anyDuplicated(rownames(validation_expr))) stop("Duplicate sample IDs")
if (anyDuplicated(colnames(reference_expr))) stop("Duplicate gene IDs")
if (any(!is.finite(reference_expr)) || any(!is.finite(validation_expr))) stop("Non-finite expression values")

assignment <- fread(MODULE_ASSIGNMENT)
assignment[, gene_id := as.character(gene_id)]
assignment <- assignment[match(colnames(reference_expr), gene_id)]
if (anyNA(assignment$gene_id) || !identical(as.character(assignment$module), unname(colors))) {
  stop("module_assignment.csv does not match the saved module colors")
}

reference_samples <- fread(REFERENCE_SAMPLE_TABLE)
if (!identical(rownames(reference_expr), reference_samples$sample_id)) stop("Reference sample table is not aligned")

gse47460 <- GEOquery::getGEO(filename = RAW_GSE47460, GSEMatrix = TRUE)
if (is.list(gse47460)) gse47460 <- gse47460[[1]]
pheno <- Biobase::pData(gse47460)
disease_state <- parse_characteristic(pheno, "disease state")
keep <- disease_state %in% c("Control", "Chronic Obstructive Lung Disease", "COPD")
keep[is.na(keep)] <- FALSE
kept_pheno <- pheno[keep, , drop = FALSE]
kept_disease <- disease_state[keep]
donor_id <- sub("^((LT[0-9]+)).*$", "\\1", as.character(kept_pheno$title))

validation_qc <- data.table(
  sample_id = rownames(kept_pheno),
  donor_id = donor_id,
  title = as.character(kept_pheno$title),
  disease_state = kept_disease,
  connectivity_z = sample_connectivity_z(validation_expr)
)
validation_qc[, connectivity_outlier_z_lt_neg2_5 := connectivity_z < -2.5]
write_dt(validation_qc, file.path(OUT_DIR, "GSE47460_validation_sample_metadata_and_qc.csv"))

if (!identical(rownames(validation_expr), validation_qc$sample_id)) stop("Validation matrix is not aligned to the GEO samples")
if (!all(grepl("^LT[0-9]+$", validation_qc$donor_id))) stop("Could not parse every GSE47460 donor ID")
if (anyDuplicated(validation_qc$donor_id)) stop("GSE47460 validation set has repeated donors")
if (length(intersect(rownames(reference_expr), rownames(validation_expr)))) stop("Reference/test GSM IDs overlap")

module_sizes <- as.data.table(table(colors))
setnames(module_sizes, c("module", "full_module_size"))
module_sizes[, module := as.character(module)]
write_dt(module_sizes[order(-full_module_size)], file.path(OUT_DIR, "reference_module_sizes.csv"))

# One array per unique donor: 91 controls and 145 COPD.
stopifnot(
  identical(dim(reference_expr), c(147L, 8873L)), identical(dim(validation_expr), c(236L, 8873L)),
  sum(validation_qc$disease_state == "Control") == 91, sum(validation_qc$disease_state == "Chronic Obstructive Lung Disease") == 145,
  length(unique(colors)) == 9L
)

multi_expr <- list(
  GSE76925_reference = list(data = reference_expr),
  GSE47460_validation = list(data = validation_expr)
)
multi_colors <- list(GSE76925_reference = colors)

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
    quickCor = 0,
    parallelCalculation = FALSE,
    savePermutedStatistics = TRUE,
    permutedStatisticsFile = file.path(OUT_DIR, "permutedStats.RData"),
    verbose = 3
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

full_stats <- extract_preservation(mp, colors)
summary_stats <- full_stats[, .(
  module, full_module_size, evaluated_module_size,
  Zsummary, medianRank, preservation_level
)]
write_dt(full_stats, file.path(OUT_DIR, "module_preservation_full_statistics.csv"))
write_dt(summary_stats, file.path(OUT_DIR, "module_preservation_summary.csv"))
write_dt(summary_stats[!module %chin% c("grey", "gold")], file.path(OUT_DIR, "module_preservation_no_grey_gold.csv"))
write_dt(
  summary_stats[!module %chin% c("grey", "gold"), .N, by = preservation_level],
  file.path(OUT_DIR, "preservation_level_counts_no_grey_gold.csv")
)
save_separate_plots(summary_stats)

stopifnot(all(summary_stats[!module %chin% c("grey", "gold"), is.finite(Zsummary) & is.finite(medianRank) & full_module_size == evaluated_module_size]))

write_dt(data.table(
  metric = c(
    "reference_samples", "validation_samples", "validation_unique_donors",
    "reference_genes", "common_genes_used", "non_grey_modules", "n_permutations",
    "network_type", "correlation", "max_module_size", "elapsed_seconds"
  ),
  value = c(
    nrow(reference_expr), nrow(validation_expr), length(unique(validation_qc$donor_id)),
    ncol(reference_expr), ncol(validation_expr), length(setdiff(unique(colors), "grey")),
    N_PERMUTATIONS, "signed", "Pearson", MAX_MODULE_SIZE, sprintf("%.2f", elapsed_seconds)
  )
), file.path(OUT_DIR, "preservation_run_summary.csv"))

save(mp, summary_stats, full_stats, colors,
  file = file.path(OUT_DIR, "COPD_to_GSE47460_module_preservation.RData")
)

capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))
message("Lung module preservation written to ", OUT_DIR)
