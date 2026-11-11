# Export per-member disease effects for Fig 3d after the adopted 10-19 exclusion of GSE48350 GSM1176257.
# No new analysis: prep()/gene_eff()/dl() are copied verbatim from
# 03_candidate_analysis/AxisA_replication_after_exclusion_20261019/code/01_recompute.R; the script asserts that the
# direction counts reproduce that branch's member_direction_counts.tsv (110 / 5 / 22 / 34).
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 06_resubmission_20261008/code/main_figures/00_member_effects.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
OUT <- "06_resubmission_20261008/code/main_figures/source_data/Fig3d_member_effects.tsv"
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
  list(id = c$id, loge = loge, dat = dat, rhs = rhs, X = X) }
gene_eff <- function(p) { g <- intersect(axisA, rownames(p$loge)); z <- t(scale(t(p$loge[g, , drop = FALSE]))); B <- qr.coef(qr(p$X), t(z)); res <- t(z) - p$X %*% B
  df <- nrow(p$X) - ncol(p$X); se <- sqrt(colSums(res^2) / df * solve(crossprod(p$X))["disease", "disease"]); data.table(cohort = p$id, entrez_id = g, beta = B["disease", ], se = se) }
dl <- function(y, s, method = "DL", test = "z") { f <- rma(yi = y, sei = s, method = method, test = test); data.table(effect_SD = as.numeric(f$b), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2) }
pp <- lapply(coh[VAL], prep)
ge <- rbindlist(lapply(pp, gene_eff)); pg <- ge[cohort %in% AD, if (.N >= 2) dl(beta, se) else data.table(effect_SD = beta, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, I2 = NA_real_), by = entrez_id]
lu <- ge[cohort == "GSE47460", .(entrez_id, lung = beta, lung_low = beta - qnorm(.975) * se, lung_high = beta + qnorm(.975) * se)]
m <- merge(lu, pg[, .(entrez_id, ad = effect_SD, ad_low = ci_low, ad_high = ci_high)], by = "entrez_id")
sym <- fread(P("05_github_release/copd-ad-convergence/resources/AxisA_confirmed_genes.csv"), colClasses = "character")[, .(entrez_id, symbol = current_symbol)]
m <- merge(m, sym, by = "entrez_id", all.x = TRUE)
m[, class := fifelse(lung > 0 & ad > 0, "both_up", fifelse(lung < 0 & ad < 0, "both_down", "opposite"))]
m[, both_CI_above_0 := lung_low > 0 & ad_low > 0]
ref <- fread(P("03_candidate_analysis/AxisA_replication_after_exclusion_20261019/results/member_direction_counts.tsv"))
got <- c(sum(m$class == "both_up"), sum(m$class == "both_down"), sum(m$class == "opposite"), sum(m$both_CI_above_0), nrow(m))
stopifnot(identical(as.integer(got), as.integer(unlist(ref[1, .(both_up, both_down, opposite, both_CI_above_0, n)]))))
fwrite(m[order(-ad)], P(OUT), sep = "\t"); cat("written", OUT, "counts", got, "\n")
