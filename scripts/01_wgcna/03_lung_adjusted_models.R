source("R/runtime.R")

source("R/paths.R")

suppressPackageStartupMessages(library(data.table))

ENV_RDATA <- file.path(COPD_WGCNA, "copd_wgcna_env.RData")
MODULE_ASSIGNMENT <- file.path(COPD_WGCNA, "module_assignment.csv")
OUT_DIR <- file.path(LUNG_RESULTS, "lung_adjusted_models")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ENV_RDATA) || !file.exists(MODULE_ASSIGNMENT)) stop("Missing lung network inputs")

run_env <- new.env(parent = emptyenv())
load(ENV_RDATA, envir = run_env)
clinical <- as.data.table(run_env$clinical_train, keep.rownames = "row_id")
mes <- as.data.table(run_env$MEs, keep.rownames = "sample_id")
if (!identical(mes$sample_id, clinical$sample_id)) stop("MEs and clinical table are not aligned")

module_cols <- setdiff(names(mes), c("sample_id", "MEgrey"))
if (length(module_cols) != 8L) stop("Expected 8 non-grey lung modules")

z <- function(x) as.numeric(scale(as.numeric(x)))

model_data <- merge(clinical, mes, by = "sample_id", sort = FALSE)
if (!identical(model_data$sample_id, clinical$sample_id)) stop("Merge changed sample order")

fit_family <- function(data, main_var, covars, family_id, population, main_label) {
  rows <- vector("list", length(module_cols))
  for (i in seq_along(module_cols)) {
    me_col <- module_cols[[i]]
    d <- data[, c(me_col, main_var, covars), with = FALSE]
    setnames(d, c(me_col, main_var), c("ME_raw", "main"))
    d <- d[complete.cases(d)]
    d[, ME_z := z(ME_raw)]
    rhs <- paste(c("main", covars), collapse = " + ")
    fit_raw <- lm(as.formula(paste("ME_raw ~", rhs)), data = d)
    fit_z <- lm(as.formula(paste("ME_z ~", rhs)), data = d)
    c_raw <- summary(fit_raw)$coefficients["main", ]
    c_z <- summary(fit_z)$coefficients["main", ]
    rows[[i]] <- data.table(
      family_id = family_id,
      population = population,
      main_variable = main_label,
      module = sub("^ME", "", me_col),
      n = nrow(d),
      raw_ME_beta = unname(c_raw[["Estimate"]]),
      beta_outcome_SD_per_unit = unname(c_z[["Estimate"]]),
      se_outcome_SD = unname(c_z[["Std. Error"]]),
      t_value = unname(c_z[["t value"]]),
      p_value = unname(c_z[["Pr(>|t|)"]]),
      df_residual = df.residual(fit_z),
      outcome_ME_SD = sd(d$ME_raw),
      covariates = paste(covars, collapse = ";")
    )
  }
  out <- rbindlist(rows)
  out[, fdr_bh_across_8_non_grey_modules := p.adjust(p_value, method = "BH")]
  setorder(out, fdr_bh_across_8_non_grey_modules, p_value)
  out
}


# Fig. 2d, lung side: COPD status adjusted for age, sex and pack-years in all 147 participants.
status_results <- fit_family(
  model_data, "copd_status", c("age_z", "sex_male", "packyears_z"),
  "COPD_status", "all 147 participants", "COPD versus control"
)
fwrite(status_results, file.path(OUT_DIR, "full_cohort_COPD_status_models.csv"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "session_info.txt"))
