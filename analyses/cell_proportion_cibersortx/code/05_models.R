# Composition adjustment of the Axis A score with CIBERSORTx fractions.
# Base model = adopted replication model (prep()/score_eff() copied from
# AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R); the base effects are checked against
# that branch before any adjustment. Inputs: results/raw/<version>_<cohort>.tsv (CIBERSORTx web output, relative mode, S-mode).
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/05_models.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
OUT <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
coh <- readRDS(P(file.path(RB, "source_data/cohorts.rds")))$cohorts
LUNG <- c("GSE76925", "GSE47460"); AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
TARGET <- c(lung = "Myeloid", brain = "Immune")
EXCL <- fread(P(file.path(RB, "results/adopted_sample_exclusions_20261019.tsv")))   # 10-19 author decision: drop the duplicate GSE48350 array
organ_of <- function(k) if (k %in% LUNG) "lung" else "brain"

prep <- function(c, covars = c$covars) {
  dat <- copy(c$dat)[!sample %in% EXCL[cohort == c$id]$sample]; dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat); stopifnot(qr(X)$rank == ncol(X))
  e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); keep_g <- filterByExpr(y, group = dat$disease)
    y <- calcNormFactors(y[keep_g, , keep.lib.sizes = FALSE]); loge <- cpm(y, log = TRUE, prior.count = 2) } else loge <- e
  R <- singscore::rankGenes(loge); d <- copy(dat)
  d$score_z <- as.numeric(scale(singscore::simpleScore(R, upSet = rownames(R)[rownames(R) %in% axisA])$TotalScore))
  list(id = c$id, dat = d, rhs = rhs)
}
fit <- function(d, rhs) { f <- lm(as.formula(paste("score_z ~", rhs)), d); co <- summary(f)$coefficients["disease", ]
  data.table(effect_SD = co[[1]], se = co[[2]], ci_low = confint(f)["disease", 1], ci_high = confint(f)["disease", 2], p = co[[4]], n = nrow(d)) }
dl <- function(y, s) { f <- rma(yi = y, sei = s, method = "DL")
  data.table(effect_SD = as.numeric(f$b), se = f$se, ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2, k = f$k) }
read_frac <- function(v, k) {
  x <- fread(P(file.path(OUT, sprintf("results/raw/%s_%s.tsv", v, k))), check.names = FALSE)
  map <- fread(P(file.path(OUT, sprintf("source_data/%s_mixture_samples.tsv", k))))
  stopifnot(setequal(x$Mixture, map$upload_name)); x[, sample := map$sample[match(Mixture, map$upload_name)]]
  cells <- setdiff(names(x), c("Mixture", "sample", "P-value", "Correlation", "RMSE"))
  stopifnot(all(abs(rowSums(x[, ..cells]) - 1) < 1e-3)); list(tab = x, cells = cells)
}

pp <- lapply(coh[c(LUNG, AD)], prep)
adopted <- fread(P(file.path(RB, "results/replication_sensitivity_by_cohort.tsv")))[analysis == "primary" & !cohort %in% EXCL$cohort]
adopted <- rbind(adopted, fread(P(file.path(RB, "results/duplicate_impact.tsv")))[analysis == "GSE48350 without GSM1176257" & cohort == "GSE48350", .(analysis, cohort, effect_SD)], fill = TRUE)
VERS <- Filter(function(v) all(file.exists(P(file.path(OUT, sprintf("results/raw/%s_%s.tsv", v, c(LUNG, AD)))))), c("B", "A")); cat("versions:", VERS, "\n")
res <- list(); pdiff <- list(); fq <- list()
for (v in VERS) for (k in c(LUNG, AD)) {
  fr <- read_frac(v, k); p <- pp[[k]]; d <- merge(p$dat, fr$tab, by = "sample"); stopifnot(nrow(d) == nrow(p$dat))   # excluded samples stay in the uploaded mixture
  tgt <- TARGET[[organ_of(k)]]; ref <- names(which.max(colMeans(fr$tab[, fr$cells, with = FALSE])))
  for (cc in fr$cells) d[[paste0("p_", make.names(cc))]] <- as.numeric(scale(d[[cc]]))
  others <- paste0("p_", make.names(setdiff(fr$cells, ref)))
  b <- fit(d, p$rhs); a1 <- fit(d, paste(p$rhs, "+", paste0("p_", make.names(tgt)))); a2 <- fit(d, paste(p$rhs, "+", paste(others, collapse = " + ")))
  stopifnot(k == LUNG[1] || isTRUE(abs(b$effect_SD - adopted$effect_SD[adopted$cohort == k]) < 1e-8))
  res[[length(res) + 1]] <- rbind(cbind(version = v, cohort = k, organ = organ_of(k), model = "base", adjusted_for = "", b),
    cbind(version = v, cohort = k, organ = organ_of(k), model = "target proportion", adjusted_for = tgt, a1),
    cbind(version = v, cohort = k, organ = organ_of(k), model = "all proportions except reference", adjusted_for = paste0("reference = ", ref), a2))
  for (cc in fr$cells) { s <- summary(lm(as.formula(paste0("p_", make.names(cc), " ~ ", p$rhs)), d))$coefficients["disease", ]
    pdiff[[length(pdiff) + 1]] <- data.table(version = v, cohort = k, cell = cc, mean_control = mean(d[disease == 0][[cc]]), mean_case = mean(d[disease == 1][[cc]]),
      diff_SD = s[[1]], se = s[[2]], p = s[[4]], r_with_score = cor(d$score_z, d[[cc]])) }
  fq[[length(fq) + 1]] <- data.table(version = v, cohort = k, n = nrow(d), median_correlation = median(d$Correlation), median_RMSE = median(d$RMSE),
    n_p_lt_0.05 = sum(d$`P-value` < 0.05), reference_cell = ref)
}
res <- rbindlist(res); pdiff <- rbindlist(pdiff); fq <- rbindlist(fq)
base_eff <- res[model == "base", .(version, cohort, base_effect = effect_SD)]
res <- merge(res, base_eff, by = c("version", "cohort")); res[, pct_change := 1 - effect_SD / base_effect]

# pooled AD (DerSimonian-Laird), before and after adjustment, plus lung GSE47460 alone
pool <- res[cohort %in% AD, cbind(dl(effect_SD, se), n_total = sum(n)), by = .(version, model)]
pool[, base := pool[model == "base" & version == .BY$version]$effect_SD, by = version][, pct_change := 1 - effect_SD / base]

# discovery GSE76925: source module red eigengene (Figure 2d model), as in the 10-02 dtangle branch
e1 <- new.env(); load(P("WGCNA_figures/reproduction_materials/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData"), envir = e1)
ld <- as.data.table(e1$clinical_train); ld[, red := e1$MEs[sample_id, "MEred"]]
red <- list()
for (v in VERS) { fr <- read_frac(v, "GSE76925"); d <- merge(ld, fr$tab, by.x = "sample_id", by.y = "sample"); stopifnot(nrow(d) == nrow(ld))
  d[, p_Myeloid := as.numeric(scale(Myeloid))]; ref <- names(which.max(colMeans(fr$tab[, fr$cells, with = FALSE])))
  for (cc in setdiff(fr$cells, ref)) d[[paste0("p_", make.names(cc))]] <- as.numeric(scale(d[[cc]]))
  g <- function(adj) { f <- lm(as.formula(paste("scale(red) ~ copd_status + age_z + sex_male + packyears_z", adj)), d); co <- summary(f)$coefficients["copd_status", ]
    data.table(beta = co[[1]], ci_low = confint(f)["copd_status", 1], ci_high = confint(f)["copd_status", 2], p = co[[4]]) }
  red[[v]] <- rbind(cbind(version = v, model = "base", g("")), cbind(version = v, model = "Myeloid", g("+ p_Myeloid")),
                    cbind(version = v, model = "all except reference", g(paste("+", paste(paste0("p_", make.names(setdiff(fr$cells, ref))), collapse = " + "))))) }
red <- rbindlist(red); stopifnot(abs(red[model == "base"]$beta[1] - 0.686) < 0.001)
dt <- fread(P("03_candidate_analysis/Deconvolution_composition_adjustment_20261002/results/adjusted_models.tsv"))[organ == "lung" & module == "red"]

fwrite(res, P(file.path(OUT, "results/adjusted_models_by_cohort.tsv")), sep = "\t")
fwrite(pool, P(file.path(OUT, "results/adjusted_models_AD_pooled.tsv")), sep = "\t")
fwrite(pdiff, P(file.path(OUT, "results/proportion_case_control.tsv")), sep = "\t")
fwrite(red, P(file.path(OUT, "results/GSE76925_red_module_adjusted.tsv")), sep = "\t")
fwrite(fq, P(file.path(OUT, "qa/cibersortx_fit_quality.tsv")), sep = "\t")
print(res[version == "B", .(cohort, model, effect_SD = round(effect_SD, 3), ci = sprintf("%.2f to %.2f", ci_low, ci_high), p = signif(p, 3), pct = round(pct_change, 3))])
print(pool[, .(version, model, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), I2 = round(I2, 1), pct = round(pct_change, 3))])
print(pdiff[version == "B" & cell %in% TARGET, .(cohort, cell, mean_control = round(mean_control, 3), mean_case = round(mean_case, 3), diff_SD = round(diff_SD, 2), p = signif(p, 3), r_with_score = round(r_with_score, 2))])
print(red); print(dt[, .(method, version, adjusted_for, beta_base, beta_adj, lo_adj, hi_adj)]); print(fq)
writeLines(capture.output(sessionInfo()), P(file.path(OUT, "qa/sessionInfo_05.txt")))
