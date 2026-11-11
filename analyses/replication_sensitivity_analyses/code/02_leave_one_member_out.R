# Leave-one-member-out on the adopted sample set (GSE48350 GSM1176257 excluded, 10-19 decision).
# Same procedure as 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R (section C);
# prep()/score_eff()/dl() are copied from 01_recompute.R in this branch; only the sample set differs from the 10-18 run.
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_replication_after_exclusion_20261019/code/02_leave_one_member_out.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
O <- "03_candidate_analysis/AxisA_replication_after_exclusion_20261019"; RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
coh <- readRDS(P(file.path(RB, "source_data/cohorts.rds")))$cohorts; EXCL <- fread(P(file.path(RB, "results/adopted_sample_exclusions_20261019.tsv")))
AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"); VAL <- c("GSE47460", AD)
prep <- function(c, keep = NULL, covars = c$covars) {
  dat <- copy(c$dat)[!sample %in% EXCL[cohort == c$id]$sample]; if (!is.null(keep)) dat <- dat[keep(dat)]
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat); e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); kg <- filterByExpr(y, group = dat$disease); y <- calcNormFactors(y[kg, , keep.lib.sizes = FALSE]); loge <- cpm(y, log = TRUE, prior.count = 2) } else loge <- e
  list(id = c$id, loge = loge, dat = dat, rhs = rhs, X = X, R = singscore::rankGenes(loge)) }
score_eff <- function(p, genes = axisA) { up <- rownames(p$R)[rownames(p$R) %in% genes]; d <- copy(p$dat)
  d$score_z <- as.numeric(scale(singscore::simpleScore(p$R, upSet = up)$TotalScore)); f <- lm(as.formula(paste("score_z ~", p$rhs)), d); co <- summary(f)$coefficients["disease", ]
  data.table(effect_SD = co[[1]], se = co[[2]], ci_low = confint(f)["disease", 1], ci_high = confint(f)["disease", 2], p = co[[4]], n_case = sum(d$disease == 1), n_control = sum(d$disease == 0)) }
dl <- function(y, s, method = "DL", test = "z") { f <- rma(yi = y, sei = s, method = method, test = test); data.table(effect_SD = as.numeric(f$b), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2) }
pp <- lapply(coh[VAL], prep)
sc <- rbindlist(lapply(VAL, function(k) cbind(cohort = k, score_eff(pp[[k]]))))
ref <- fread(P(file.path(O, "results/by_cohort.tsv")))
stopifnot(max(abs(sc$effect_SD - ref[match(sc$cohort, cohort)]$effect_SD)) < 1e-10)   # reproduces the adopted per-cohort effects
loo <- rbindlist(lapply(axisA, function(g) {
  r <- rbindlist(lapply(VAL, function(k) if (g %in% rownames(pp[[k]]$R)) cbind(cohort = k, score_eff(pp[[k]], setdiff(axisA, g))) else NULL))
  if (!nrow(r)) return(NULL)
  ad <- r[cohort %in% AD]; pooled <- if (nrow(ad) == length(AD)) dl(ad$effect_SD, ad$se) else dl(c(ad$effect_SD, sc[cohort %in% setdiff(AD, ad$cohort)]$effect_SD), c(ad$se, sc[cohort %in% setdiff(AD, ad$cohort)]$se))
  data.table(entrez_id = g, GSE47460 = if ("GSE47460" %in% r$cohort) r[cohort == "GSE47460"]$effect_SD else sc[cohort == "GSE47460"]$effect_SD,
             AD_pooled = pooled$effect_SD, AD_pooled_p = pooled$p)
}))
base_ad <- dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se)
loo[, `:=`(delta_GSE47460 = GSE47460 - sc[cohort == "GSE47460"]$effect_SD, delta_AD = AD_pooled - base_ad$effect_SD)]
ax <- fread(P("05_github_release/copd-ad-convergence/resources/AxisA_confirmed_genes.csv"), colClasses = "character")
loo <- merge(ax[, .(entrez_id, symbol = current_symbol)], loo, by = "entrez_id")
fwrite(loo[order(delta_AD)], P(file.path(O, "results/leave_one_member_out.tsv")), sep = "\t")
cat("n members:", nrow(loo), "| GSE47460 range:", round(range(loo$GSE47460), 3), "| AD pooled range:", round(range(loo$AD_pooled), 3),
    "| max AD p:", signif(max(loo$AD_pooled_p), 3), "\n")
