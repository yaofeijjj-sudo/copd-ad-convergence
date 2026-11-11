source("R/runtime.R")
suppressPackageStartupMessages(library(data.table))
l <- fread("results/wgcna/lung/lung_adjusted_models/full_cohort_COPD_status_models.csv")
l <- l[module %in% c("red", "turquoise")]
l[, `:=`(effect = beta_outcome_SD_per_unit, se = se_outcome_SD, fdr = fdr_bh_across_8_non_grey_modules, organ = "Lung")]
b <- fread("results/wgcna/brain/AD_WGCNA_unadjusted/ad_module_trait_models.csv")[term == "Disease_State" & model_name == "Disease_State_adjusted" & module %in% c("blue", "turquoise")]
e <- new.env()
load("results/wgcna/brain/AD_WGCNA_unadjusted/ad_wgcna_env.RData", e)
b[, scale := vapply(module, function(m) sd(e$MEs[[paste0("ME", m)]]), numeric(1))]
b[, `:=`(effect = beta / scale, se = std_error / scale, organ = "Brain")]
d <- rbindlist(list(l[, .(organ, module, effect, se, fdr)], b[, .(organ, module, effect, se, fdr)]))
d[, `:=`(axis = ifelse(module == "turquoise", "Axis B", "Axis A"), lower = effect - 1.96 * se, upper = effect + 1.96 * se)]
setorder(d, axis, organ)
dir.create("results/module_overlap_figures/source_data", recursive = TRUE, showWarnings = FALSE)
fwrite(d, "results/module_overlap_figures/source_data/adjusted_effects.csv")
