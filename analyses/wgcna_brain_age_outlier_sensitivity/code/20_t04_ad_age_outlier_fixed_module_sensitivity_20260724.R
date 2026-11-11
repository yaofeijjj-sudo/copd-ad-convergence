#!/usr/bin/env Rscript

# T04: fixed-module sensitivity to restoring the seven GSE33000 age outliers.
# The analysis contract is frozen in:
# T04_AD_age_outlier_fixed_module_sensitivity_20260724/docs/
# T04_preregistered_analysis_contract_20260724.md

suppressPackageStartupMessages({
  library(Biobase)
  library(data.table)
  library(GEOquery)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(impute)
  library(limma)
})

options(stringsAsFactors = FALSE)
if (identical(Sys.getenv("T04_WARN_AS_ERROR"), "1")) options(warn = 2)
set.seed(20260724)

find_analysis_root <- function() {
  starts <- unique(c(getwd(), dirname(normalizePath(commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))] |>
    sub("^--file=", "", x = _), mustWork = FALSE))))
  for (start in starts[nzchar(starts)]) {
    cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
    repeat {
      if (basename(cur) == "01_current_analysis") return(cur)
      candidate <- file.path(cur, "01_current_analysis")
      if (dir.exists(file.path(candidate, "01_storyline_current"))) return(candidate)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Could not locate 01_current_analysis.")
}

analysis_root <- find_analysis_root()
workspace_root <- dirname(analysis_root)
story_root <- file.path(analysis_root, "01_storyline_current")
raw_gse33000 <- file.path(workspace_root, "02_data", "raw", "bulk_brain_ad", "GSE33000_series_matrix.txt.gz")
current_dir <- file.path(story_root, "02_brain_wgcna", "results", "AD_WGCNA_no_residual_after_age_outlier")
out_dir <- file.path(story_root, "02_brain_wgcna", "results", "T04_AD_age_outlier_fixed_module_sensitivity_20260724")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
docs_dir <- file.path(out_dir, "docs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

current_env_path <- file.path(current_dir, "ad_wgcna_final_env.RData")
current_map_path <- file.path(current_dir, "training_probe_to_gene_map.csv")
current_model_path <- file.path(current_dir, "ad_module_trait_models.csv")
age_summary_path <- file.path(current_dir, "ad_age_outlier_removal_summary.csv")
qc_summary_path <- file.path(current_dir, "wgcna_qc_removal_summary.csv")
contract_path <- file.path(docs_dir, "T04_preregistered_analysis_contract_20260724.md")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1],
  winslash = "/", mustWork = TRUE
)

required_inputs <- c(
  raw_gse33000, current_env_path, current_map_path, current_model_path,
  age_summary_path, qc_summary_path, contract_path, script_path
)
if (any(!file.exists(required_inputs))) {
  stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
}

write_csv <- function(x, path) {
  separator <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
  fwrite(as.data.table(x), path, sep = separator, na = "NA")
}

md5_table <- function(paths, role) {
  info <- file.info(paths)
  data.table(
    role = role,
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = unname(tools::md5sum(paths))
  )
}

zscore <- function(x) {
  x <- as.numeric(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
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

learn_limma_quantile_target <- function(a) {
  a <- as.matrix(a)
  n <- dim(a)
  i <- (0:(n[1] - 1)) / (n[1] - 1)
  sorted <- array(NA_real_, n)
  for (j in seq_len(n[2])) {
    sj <- sort(a[, j], method = "quick", index.return = TRUE)
    n_obs <- length(sj$x)
    if (n_obs < n[1]) {
      sorted[, j] <- approx(
        (0:(n_obs - 1)) / (n_obs - 1),
        sj$x,
        i,
        ties = list("ordered", mean)
      )$y
    } else {
      sorted[, j] <- sj$x
    }
  }
  list(probability = i, target = rowMeans(sorted), n_features = n[1])
}

apply_limma_quantile_target <- function(a, target_info) {
  a <- as.matrix(a)
  if (nrow(a) != target_info$n_features) stop("Quantile target feature count mismatch.")
  i <- target_info$probability
  m <- target_info$target
  for (j in seq_len(ncol(a))) {
    n_obs <- sum(!is.na(a[, j]))
    ranks <- rank(a[, j])
    observed <- !is.na(a[, j])
    if (n_obs < nrow(a)) {
      a[observed, j] <- approx(
        i,
        m,
        (ranks[observed] - 1) / (n_obs - 1),
        ties = list("ordered", mean)
      )$y
    } else {
      a[, j] <- approx(
        i,
        m,
        (ranks - 1) / (nrow(a) - 1),
        ties = list("ordered", mean)
      )$y
    }
  }
  a
}

fit_fixed_module_models <- function(me, traits, analysis_label) {
  rows <- lapply(colnames(me), function(me_name) {
    module <- sub("^ME", "", me_name)
    df <- data.frame(
      ME = as.numeric(me[, me_name]),
      Disease_State = as.numeric(traits[rownames(me), "Disease_State"]),
      age_z = as.numeric(traits[rownames(me), "age_z"]),
      sex_male = as.numeric(traits[rownames(me), "sex_male"])
    )
    df <- df[complete.cases(df), , drop = FALSE]
    fit <- lm(ME ~ Disease_State + age_z + sex_male, data = df)
    sm <- summary(fit)$coefficients
    term <- "Disease_State"
    estimate <- sm[term, "Estimate"]
    se <- sm[term, "Std. Error"]
    crit <- qt(0.975, df = df.residual(fit))
    data.table(
      analysis = analysis_label,
      module = module,
      beta = estimate,
      std_error = se,
      ci_low = estimate - crit * se,
      ci_high = estimate + crit * se,
      statistic = sm[term, "t value"],
      p_value = sm[term, "Pr(>|t|)"],
      n_model = nrow(df),
      residual_df = df.residual(fit)
    )
  })
  out <- rbindlist(rows)
  out[, fdr := p.adjust(p_value, method = "BH")]
  out
}

check_row <- function(check, observed, expected, operator, pass, details = "") {
  data.table(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    operator = operator,
    status = ifelse(isTRUE(pass), "PASS", "FAIL"),
    details = details
  )
}

message("[T04] Loading frozen current AD network.")
load(current_env_path)
stopifnot(
  is.matrix(dat_expr), is.data.frame(MEs), length(module_colors) == ncol(dat_expr),
  identical(rownames(dat_expr), rownames(MEs))
)
current_dat_expr <- dat_expr
current_mes <- as.matrix(MEs)
current_traits_saved <- clinical_train
current_module_colors <- as.character(module_colors)
names(current_module_colors) <- names(module_colors)
rm(dat_expr, MEs, clinical_train, module_colors)

current_ids <- rownames(current_dat_expr)
fixed_genes <- colnames(current_dat_expr)
fixed_modules <- sub("^ME", "", colnames(current_mes))
non_grey_modules <- setdiff(fixed_modules, "grey")

message("[T04] Reading GSE33000 and reconstructing the frozen preprocessing reference.")
gse <- GEOquery::getGEO(filename = raw_gse33000, GSEMatrix = TRUE, getGPL = FALSE)
pheno_all <- Biobase::pData(gse)
keep_no_hd <- pheno_all[["disease status:ch2"]] != "Huntington's disease"
keep_no_hd[is.na(keep_no_hd)] <- FALSE
pheno <- pheno_all[keep_no_hd, , drop = FALSE]
expr_raw <- Biobase::exprs(gse)[, rownames(pheno), drop = FALSE]

age <- numeric_clean(pheno[["age:ch2"]])
sex <- sex_to_male(pheno[["gender:ch2"]])
disease <- fifelse(grepl("Alzheimer", pheno[["disease status:ch2"]], ignore.case = TRUE), 1, 0)
names(age) <- names(sex) <- names(disease) <- rownames(pheno)

age_rows <- lapply(sort(unique(disease)), function(group_value) {
  idx <- which(disease == group_value & is.finite(age))
  q <- quantile(age[idx], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  out <- idx[age[idx] < lower | age[idx] > upper]
  data.table(
    sample_id = rownames(pheno)[out],
    Disease_State = group_value,
    disease_group = ifelse(group_value == 1, "AD", "Control"),
    Age = age[out],
    sex_male = sex[out],
    q1 = q[1],
    q3 = q[2],
    iqr = iqr,
    lower_fence = lower,
    upper_fence = upper,
    original_exclusion_reason = "within-group age outside 1.5*IQR fence"
  )
})
restored_samples <- rbindlist(age_rows)
restored_ids <- restored_samples$sample_id
saved_age_summary <- fread(age_summary_path)
saved_age_outlier_ids <- unique(unlist(strsplit(saved_age_summary$removed_samples, ";", fixed = TRUE)))
saved_age_outlier_ids <- saved_age_outlier_ids[nzchar(saved_age_outlier_ids)]

qc_summary <- fread(qc_summary_path)
connectivity_ids <- qc_summary[
  removal_step == "connectivity_outliers" & nzchar(removed_ids),
  unlist(strsplit(removed_ids, ";", fixed = TRUE))
]
connectivity_ids <- unique(connectivity_ids[nzchar(connectivity_ids)])

reference_pre_qc_ids <- setdiff(rownames(pheno), restored_ids)
sensitivity_ids <- setdiff(rownames(pheno), connectivity_ids)

expr_log <- as.matrix(expr_raw)
storage.mode(expr_log) <- "numeric"
if (max(expr_log, na.rm = TRUE) > 100) {
  expr_log <- log2(pmax(expr_log, 1))
}
reference_log <- expr_log[, reference_pre_qc_ids, drop = FALSE]
reference_target_info <- learn_limma_quantile_target(reference_log)
reference_target <- reference_target_info$target
reference_norm <- limma::normalizeBetweenArrays(reference_log, method = "quantile")
rownames(reference_norm) <- rownames(reference_log)
colnames(reference_norm) <- colnames(reference_log)
reference_norm_imputed <- impute::impute.knn(reference_norm)$data
rownames(reference_norm_imputed) <- rownames(reference_norm)
colnames(reference_norm_imputed) <- colnames(reference_norm)

restored_norm <- apply_limma_quantile_target(
  expr_log[, restored_ids, drop = FALSE],
  target_info = reference_target_info
)
rownames(restored_norm) <- rownames(expr_log)
colnames(restored_norm) <- restored_ids
combined_for_new_imputation <- cbind(reference_norm, restored_norm)
combined_imputed <- impute::impute.knn(combined_for_new_imputation)$data
restored_norm_imputed <- combined_imputed[, restored_ids, drop = FALSE]

probe_map <- fread(current_map_path, colClasses = "character")
probe_map <- probe_map[gene_id %in% fixed_genes]
probe_map <- probe_map[match(fixed_genes, gene_id)]
if (anyNA(probe_map$selected_probe_id) || any(!probe_map$selected_probe_id %in% rownames(expr_log))) {
  stop("Frozen selected probes do not map one-to-one to the GSE33000 expression matrix.")
}

reference_gene_expr <- t(reference_norm_imputed[probe_map$selected_probe_id, current_ids, drop = FALSE])
colnames(reference_gene_expr) <- fixed_genes
rownames(reference_gene_expr) <- current_ids
max_expression_reconstruction_error <- max(abs(reference_gene_expr - current_dat_expr))

restored_gene_expr <- t(restored_norm_imputed[probe_map$selected_probe_id, restored_ids, drop = FALSE])
colnames(restored_gene_expr) <- fixed_genes
rownames(restored_gene_expr) <- restored_ids

fixed_expr_sensitivity <- rbind(current_dat_expr, restored_gene_expr)
fixed_expr_sensitivity <- fixed_expr_sensitivity[sensitivity_ids, fixed_genes, drop = FALSE]

message("[T04] Learning fixed PC1 parameters on the 459 current samples and projecting the restored samples.")
projected_mes <- matrix(
  NA_real_,
  nrow = nrow(fixed_expr_sensitivity),
  ncol = length(fixed_modules),
  dimnames = list(rownames(fixed_expr_sensitivity), paste0("ME", fixed_modules))
)
projection_checks <- vector("list", length(fixed_modules))
projection_parameters <- vector("list", length(fixed_modules))
names(projection_parameters) <- fixed_modules

for (i in seq_along(fixed_modules)) {
  module <- fixed_modules[i]
  genes <- fixed_genes[current_module_colors[fixed_genes] == module]
  ref_x <- current_dat_expr[current_ids, genes, drop = FALSE]
  pc <- prcomp(ref_x, center = TRUE, scale. = TRUE, rank. = 1)
  raw_ref_score <- as.numeric(pc$x[, 1])
  stored_score <- as.numeric(current_mes[current_ids, paste0("ME", module)])
  alignment_slope <- sum(raw_ref_score * stored_score) / sum(raw_ref_score^2)

  standardized_all <- sweep(fixed_expr_sensitivity[, genes, drop = FALSE], 2, pc$center, "-")
  standardized_all <- sweep(standardized_all, 2, pc$scale, "/")
  raw_all_score <- as.numeric(standardized_all %*% pc$rotation[, 1])
  projected_mes[, paste0("ME", module)] <- raw_all_score * alignment_slope

  reproduced <- projected_mes[current_ids, paste0("ME", module)]
  projection_checks[[i]] <- data.table(
    module = module,
    n_genes = length(genes),
    pearson_r_saved_vs_projected = cor(stored_score, reproduced),
    max_abs_error_saved_vs_projected = max(abs(stored_score - reproduced)),
    alignment_slope = alignment_slope
  )
  projection_parameters[[module]] <- list(
    genes = genes,
    center = pc$center,
    scale = pc$scale,
    rotation_pc1 = pc$rotation[, 1],
    alignment_slope = alignment_slope
  )
}
projection_checks <- rbindlist(projection_checks)

traits_all <- data.frame(
  row.names = rownames(pheno),
  sample_id = rownames(pheno),
  Disease_State = disease[rownames(pheno)],
  Age = age[rownames(pheno)],
  sex_male = sex[rownames(pheno)],
  check.names = FALSE
)
traits_current <- traits_all[current_ids, , drop = FALSE]
traits_current$age_z <- zscore(traits_current$Age)
traits_sensitivity <- traits_all[sensitivity_ids, , drop = FALSE]
traits_sensitivity$age_z <- zscore(traits_sensitivity$Age)

current_effects <- fit_fixed_module_models(
  projected_mes[current_ids, , drop = FALSE],
  traits_current,
  "Current age-filtered"
)
sensitivity_effects <- fit_fixed_module_models(
  projected_mes[sensitivity_ids, , drop = FALSE],
  traits_sensitivity,
  "Restored 7 age-outlier samples"
)
all_effects <- rbindlist(list(current_effects, sensitivity_effects))

all_effects[, non_grey_abs_rank := NA_integer_]
for (label in unique(all_effects$analysis)) {
  idx <- which(all_effects$analysis == label & all_effects$module != "grey")
  all_effects[idx, non_grey_abs_rank := rank(-abs(beta), ties.method = "min")]
}

current_wide <- copy(current_effects)
setnames(
  current_wide,
  setdiff(names(current_wide), "module"),
  paste0(setdiff(names(current_wide), "module"), "_current")
)
sensitivity_wide <- copy(sensitivity_effects)
setnames(
  sensitivity_wide,
  setdiff(names(sensitivity_wide), "module"),
  paste0(setdiff(names(sensitivity_wide), "module"), "_sensitivity")
)
comparison <- merge(current_wide, sensitivity_wide, by = "module", sort = FALSE)
comparison[, `:=`(
  direction_current = fifelse(beta_current > 0, "positive", fifelse(beta_current < 0, "negative", "zero")),
  direction_sensitivity = fifelse(beta_sensitivity > 0, "positive", fifelse(beta_sensitivity < 0, "negative", "zero")),
  direction_preserved = sign(beta_current) == sign(beta_sensitivity),
  beta_change = beta_sensitivity - beta_current,
  absolute_beta_change = abs(beta_sensitivity) - abs(beta_current),
  relative_absolute_beta_change = (abs(beta_sensitivity) - abs(beta_current)) / abs(beta_current),
  attenuation_fraction = 1 - abs(beta_sensitivity) / abs(beta_current),
  attenuation_percent = 100 * (1 - abs(beta_sensitivity) / abs(beta_current))
)]

comparison[, non_grey_abs_rank_current := NA_integer_]
comparison[, non_grey_abs_rank_sensitivity := NA_integer_]
comparison[module != "grey", non_grey_abs_rank_current := rank(-abs(beta_current), ties.method = "min")]
comparison[module != "grey", non_grey_abs_rank_sensitivity := rank(-abs(beta_sensitivity), ties.method = "min")]
comparison[, rank_change := non_grey_abs_rank_sensitivity - non_grey_abs_rank_current]

stored_models <- fread(current_model_path)
stored_disease <- stored_models[model_name == "Disease_State_adjusted" & term == "Disease_State"]
stored_disease <- stored_disease[match(current_effects$module, module)]
max_current_model_beta_error <- max(abs(stored_disease$beta - current_effects$beta))
max_current_model_se_error <- max(abs(stored_disease$std_error - current_effects$std_error))

axis_modules <- c("blue", "turquoise")
axis_labels <- c(blue = "Axis A · AD blue", turquoise = "Axis B · AD turquoise")
axis_triggers <- comparison[module %in% axis_modules, .(
  scope = "Axis primary",
  module,
  trigger_sign_reversal = sign(beta_current) != sign(beta_sensitivity),
  trigger_ci_crosses_zero = ci_low_sensitivity <= 0 & ci_high_sensitivity >= 0,
  trigger_fdr_not_significant = fdr_sensitivity >= 0.05,
  trigger_attenuation_gt_25pct = attenuation_fraction > 0.25,
  trigger_rank_beyond_top5 = non_grey_abs_rank_current <= 3 & non_grey_abs_rank_sensitivity > 5
)]

global_r <- cor(
  comparison[module != "grey", beta_current],
  comparison[module != "grey", beta_sensitivity],
  method = "pearson"
)
direction_subset <- comparison[module != "grey" & abs(beta_current) >= 0.01]
direction_concordance <- mean(direction_subset$direction_preserved)
global_triggers <- data.table(
  scope = "All non-grey modules",
  module = "all_non_grey",
  trigger_beta_vector_r_lt_0_90 = global_r < 0.90,
  trigger_direction_concordance_lt_0_80 = direction_concordance < 0.80,
  pearson_r = global_r,
  direction_concordance = direction_concordance,
  n_modules_direction_assessed = nrow(direction_subset)
)

trigger_columns <- grep("^trigger_", names(axis_triggers), value = TRUE)
axis_triggered <- any(unlist(axis_triggers[, ..trigger_columns]))
global_triggered <- any(unlist(global_triggers[, .(
  trigger_beta_vector_r_lt_0_90,
  trigger_direction_concordance_lt_0_80
)]))
full_wgcna_triggered <- axis_triggered || global_triggered
decision <- if (full_wgcna_triggered) "NO-GO: upgrade to full WGCNA sensitivity and pause T05+" else
  "GO: retain current network; do not launch full WGCNA sensitivity; unlock T05"

acceptance <- rbindlist(list(
  check_row("HD samples excluded before T04", sum(!keep_no_hd), 157, "==", sum(!keep_no_hd) == 157),
  check_row("AD/control samples before age exclusion", nrow(pheno), 467, "==", nrow(pheno) == 467),
  check_row("Restored age-outlier samples", length(restored_ids), 7, "==", length(restored_ids) == 7),
  check_row("Restored IDs exactly match current age-exclusion record",
    paste(sort(restored_ids), collapse = ";"),
    paste(sort(saved_age_outlier_ids), collapse = ";"),
    "setequal",
    setequal(restored_ids, saved_age_outlier_ids)
  ),
  check_row("Reference samples before connectivity QC", length(reference_pre_qc_ids), 460, "==", length(reference_pre_qc_ids) == 460),
  check_row("Frozen current samples", length(current_ids), 459, "==", length(current_ids) == 459),
  check_row("Frozen current sample membership",
    length(intersect(current_ids, setdiff(reference_pre_qc_ids, connectivity_ids))),
    459,
    "exact set",
    setequal(current_ids, setdiff(reference_pre_qc_ids, connectivity_ids))
  ),
  check_row("Connectivity outlier remains excluded", paste(connectivity_ids, collapse = ";"), "GSM1424097", "==",
    identical(connectivity_ids, "GSM1424097")
  ),
  check_row("Sensitivity samples", length(sensitivity_ids), 466, "==", length(sensitivity_ids) == 466),
  check_row("Frozen genes", length(fixed_genes), 8807, "==", length(fixed_genes) == 8807),
  check_row("Frozen modules including grey", length(fixed_modules), 17, "==", length(fixed_modules) == 17),
  check_row("Frozen non-grey modules", length(non_grey_modules), 16, "==", length(non_grey_modules) == 16),
  check_row("Current expression reconstruction max error", signif(max_expression_reconstruction_error, 8), 1e-10, "<=",
    max_expression_reconstruction_error <= 1e-10
  ),
  check_row("Minimum ME reproduction correlation", signif(min(projection_checks$pearson_r_saved_vs_projected), 12), 0.999999, ">=",
    min(projection_checks$pearson_r_saved_vs_projected) >= 0.999999
  ),
  check_row("Maximum ME reproduction absolute error", signif(max(projection_checks$max_abs_error_saved_vs_projected), 8), 1e-10, "<=",
    max(projection_checks$max_abs_error_saved_vs_projected) <= 1e-10
  ),
  check_row("Current disease beta reproduction max error", signif(max_current_model_beta_error, 8), 1e-10, "<=",
    max_current_model_beta_error <= 1e-10
  ),
  check_row("Current disease SE reproduction max error", signif(max_current_model_se_error, 8), 1e-10, "<=",
    max_current_model_se_error <= 1e-10
  ),
  check_row("Axis modules present", paste(sort(intersect(axis_modules, fixed_modules)), collapse = ";"),
    paste(sort(axis_modules), collapse = ";"), "==", setequal(axis_modules, intersect(axis_modules, fixed_modules))
  )
))

restored_samples[, included_in_t04_sensitivity := sample_id %in% sensitivity_ids]
setorder(restored_samples, Disease_State, Age, sample_id)
setorder(comparison, non_grey_abs_rank_current, na.last = TRUE)

write_csv(restored_samples, file.path(table_dir, "t04_restored_age_outlier_samples_20260724.csv"))
write_csv(projection_checks, file.path(table_dir, "t04_fixed_module_projection_reproduction_20260724.csv"))
write_csv(all_effects, file.path(table_dir, "t04_fixed_module_disease_effects_long_20260724.csv"))
write_csv(comparison, file.path(table_dir, "t04_fixed_module_disease_effect_comparison_20260724.csv"))
write_csv(axis_triggers, file.path(table_dir, "t04_axis_upgrade_triggers_20260724.csv"))
write_csv(global_triggers, file.path(table_dir, "t04_global_upgrade_triggers_20260724.csv"))
write_csv(acceptance, file.path(table_dir, "t04_acceptance_checks_20260724.csv"))

decision_table <- data.table(
  task = "T04",
  decision = decision,
  full_wgcna_sensitivity_triggered = full_wgcna_triggered,
  axis_primary_triggered = axis_triggered,
  global_triggered = global_triggered,
  t05_status = ifelse(full_wgcna_triggered, "blocked", "unlocked"),
  decision_rule_source = basename(contract_path)
)
write_csv(decision_table, file.path(table_dir, "t04_go_no_go_decision_20260724.csv"))

saveRDS(
  list(
    projection_parameters = projection_parameters,
    reference_quantile_target = reference_target,
    current_sample_ids = current_ids,
    restored_sample_ids = restored_ids,
    sensitivity_sample_ids = sensitivity_ids,
    fixed_genes = fixed_genes,
    fixed_module_colors = current_module_colors,
    projected_module_eigengenes = projected_mes,
    traits_current = traits_current,
    traits_sensitivity = traits_sensitivity
  ),
  file.path(out_dir, "t04_fixed_module_projection_bundle_20260724.rds"),
  compress = "xz"
)

figure_source <- comparison[module != "grey", .(
  module,
  axis_label = fifelse(module %in% names(axis_labels), axis_labels[module], "Other fixed module"),
  beta_current,
  ci_low_current,
  ci_high_current,
  fdr_current,
  beta_sensitivity,
  ci_low_sensitivity,
  ci_high_sensitivity,
  fdr_sensitivity,
  non_grey_abs_rank_current,
  non_grey_abs_rank_sensitivity
)]
write_csv(figure_source, file.path(table_dir, "t04_figure_source_data_20260724.csv"))

message("[T04] Drawing the R-only robustness figure.")
axis_plot_dt <- rbindlist(list(
  comparison[module %in% axis_modules, .(
    module, analysis = "Current\nn = 459", beta = beta_current,
    ci_low = ci_low_current, ci_high = ci_high_current
  )],
  comparison[module %in% axis_modules, .(
    module, analysis = "Restored 7\nn = 466", beta = beta_sensitivity,
    ci_low = ci_low_sensitivity, ci_high = ci_high_sensitivity
  )]
))
axis_plot_dt[, module_label := factor(axis_labels[module], levels = rev(unname(axis_labels)))]
axis_plot_dt[, analysis := factor(analysis, levels = c("Current\nn = 459", "Restored 7\nn = 466"))]

palette <- c(
  "Current\nn = 459" = "#747474",
  "Restored 7\nn = 466" = "#2B6CB0",
  blue = "#3268A8",
  turquoise = "#2C9C95",
  neutral = "#B8B8B8"
)

theme_pub <- theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#252525"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#252525"),
    axis.title = element_text(size = 7.2, colour = "#222222"),
    axis.text = element_text(size = 6.6, colour = "#222222"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.3),
    legend.key.height = unit(3.2, "mm"),
    plot.margin = margin(4, 5, 4, 4),
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0, 1)
  )

p_a <- ggplot(axis_plot_dt, aes(x = beta, y = module_label, colour = analysis)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#D0D0D0") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    position = position_dodge(width = 0.46),
    width = 0,
    linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.46), size = 2.2) +
  scale_colour_manual(values = palette[c("Current\nn = 459", "Restored 7\nn = 466")]) +
  labs(x = "Age- and sex-adjusted AD effect (β, 95% CI)", y = NULL, tag = "a") +
  theme_pub +
  theme(legend.position = c(0.25, 0.14), legend.direction = "horizontal")

scatter_dt <- comparison[module != "grey"]
scatter_dt[, point_class := fifelse(module == "blue", "blue",
  fifelse(module == "turquoise", "turquoise", "neutral")
)]
label_dt <- scatter_dt[module %in% axis_modules]
limits <- range(c(scatter_dt$beta_current, scatter_dt$beta_sensitivity))
padding <- diff(limits) * 0.09

p_b <- ggplot(scatter_dt, aes(x = beta_current, y = beta_sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = "dashed", colour = "#858585") +
  geom_point(aes(colour = point_class, size = point_class), alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = label_dt,
    aes(label = ifelse(module == "blue", "Axis A · blue", "Axis B · turquoise"), colour = point_class),
    size = 2.15,
    family = "sans",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.25,
    show.legend = FALSE,
    seed = 20260724
  ) +
  annotate(
    "text",
    x = limits[1] + padding,
    y = limits[2] - padding,
    label = sprintf("r = %.3f", global_r),
    hjust = 0,
    vjust = 1,
    size = 2.35,
    family = "sans",
    colour = "#333333"
  ) +
  scale_colour_manual(values = palette[c("blue", "turquoise", "neutral")], guide = "none") +
  scale_size_manual(values = c(blue = 2.45, turquoise = 2.45, neutral = 1.55), guide = "none") +
  coord_equal(xlim = limits + c(-padding, padding), ylim = limits + c(-padding, padding)) +
  labs(
    x = "Current age-filtered effect (β)",
    y = "Restored-sample effect (β)",
    tag = "b"
  ) +
  theme_pub

figure_grob <- gridExtra::arrangeGrob(p_a, p_b, ncol = 2, widths = c(1.08, 1))
figure_base <- file.path(figure_dir, "t04_ad_age_outlier_fixed_module_robustness")
width_in <- 183 / 25.4
height_in <- 90 / 25.4

grDevices::svg(paste0(figure_base, ".svg"), width = width_in, height = height_in, family = "sans")
grid::grid.draw(figure_grob)
dev.off()
grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = width_in, height = height_in, family = "sans")
grid::grid.draw(figure_grob)
dev.off()
grDevices::tiff(
  paste0(figure_base, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo",
  family = "sans"
)
grid::grid.draw(figure_grob)
dev.off()
grDevices::png(
  paste0(figure_base, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  type = "cairo",
  family = "sans"
)
grid::grid.draw(figure_grob)
dev.off()

input_manifest <- md5_table(
  required_inputs,
  c(
    "raw_expression", "current_wgcna_bundle", "frozen_probe_map", "current_model_table",
    "age_exclusion_summary", "connectivity_qc_summary", "preregistered_contract", "analysis_code"
  )
)
write_csv(input_manifest, file.path(table_dir, "t04_input_code_manifest_20260724.tsv"))

sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

output_paths <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
output_paths <- output_paths[file.info(output_paths)$isdir %in% FALSE]
output_hash_path <- file.path(table_dir, "t04_output_hashes_20260724.tsv")
output_paths <- setdiff(output_paths, output_hash_path)
output_hashes <- md5_table(output_paths, rep("T04 output", length(output_paths)))
write_csv(output_hashes, output_hash_path)

if (any(acceptance$status != "PASS")) {
  stop("T04 acceptance checks failed. Review t04_acceptance_checks_20260724.csv.")
}

message("[T04] ", decision)
message("[T04] Outputs: ", out_dir)
