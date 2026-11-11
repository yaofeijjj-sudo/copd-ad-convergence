# Export per-sample Axis A score and estimated target-cell proportion for Fig 5e, and the case-control proportion difference
# with its 95% CI for Fig 5f. Reuses prep()/read_frac() of the CIBERSORTx branch (version B reference without the 137 Axis A
# genes). No new test: the script asserts that r, diff_SD and P reproduce results/proportion_case_control.tsv (Table G d).
# Run from the project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 06_resubmission_20261008/code/main_figures/00_composition_samples.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
BR <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
SRC <- "06_resubmission_20261008/code/main_figures/source_data"
coh <- readRDS(P(file.path(RB, "source_data/cohorts.rds")))$cohorts
LUNG <- c("GSE76925", "GSE47460"); AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
EXCL <- fread(P(file.path(RB, "results/adopted_sample_exclusions_20261019.tsv")))
ref <- fread(P(file.path(BR, "results/proportion_case_control.tsv")))[version == "B"]
prep <- function(c, covars = c$covars) {
  dat <- copy(c$dat)[!sample %in% EXCL[cohort == c$id]$sample]; dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); keep_g <- filterByExpr(y, group = dat$disease)
    y <- calcNormFactors(y[keep_g, , keep.lib.sizes = FALSE]); loge <- cpm(y, log = TRUE, prior.count = 2) } else loge <- e
  R <- singscore::rankGenes(loge); d <- copy(dat)
  d$score_z <- as.numeric(scale(singscore::simpleScore(R, upSet = rownames(R)[rownames(R) %in% axisA])$TotalScore))
  list(id = c$id, dat = d, rhs = rhs)
}
read_frac <- function(k) {
  x <- fread(P(file.path(BR, sprintf("results/raw/B_%s.tsv", k))), check.names = FALSE)
  map <- fread(P(file.path(BR, sprintf("source_data/%s_mixture_samples.tsv", k))))
  x[, sample := map$sample[match(Mixture, map$upload_name)]]; x
}
smp <- list(); dif <- list()
for (k in c(LUNG, AD)) {
  p <- prep(coh[[k]]); tgt <- if (k %in% LUNG) "Myeloid" else "Immune"
  d <- merge(p$dat, read_frac(k)[, c("sample", tgt), with = FALSE], by = "sample"); stopifnot(nrow(d) == nrow(p$dat))
  d[, prop := get(tgt)][, prop_z := as.numeric(scale(prop))]
  f <- lm(as.formula(paste("prop_z ~", p$rhs)), d); co <- summary(f)$coefficients["disease", ]; ci <- confint(f)["disease", ]
  r <- cor(d$score_z, d$prop); rr <- ref[cohort == k & cell == tgt]
  stopifnot(nrow(rr) == 1, abs(r - rr$r_with_score) < 1e-8, abs(co[[1]] - rr$diff_SD) < 1e-8, abs(co[[4]] - rr$p) < 1e-8)
  smp[[k]] <- d[, .(cohort = k, target = tgt, sample, disease, score_z, proportion = prop)]
  dif[[k]] <- data.table(cohort = k, target = tgt, n_case = sum(d$disease == 1), n_control = sum(d$disease == 0), mean_control = rr$mean_control,
                         mean_case = rr$mean_case, diff_SD = co[[1]], ci_low = ci[[1]], ci_high = ci[[2]], p = co[[4]], r_with_score = r)
}
fwrite(rbindlist(smp), P(file.path(SRC, "Fig5e_samples.tsv")), sep = "\t")
fwrite(rbindlist(dif), P(file.path(SRC, "Fig5f_proportion_difference.tsv")), sep = "\t")
print(rbindlist(dif)[, .(cohort, target, n_case, n_control, diff = round(diff_SD, 2), lo = round(ci_low, 2), hi = round(ci_high, 2), p = signif(p, 2), r = round(r_with_score, 2))])
