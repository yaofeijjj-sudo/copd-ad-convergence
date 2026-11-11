# Recompute the adopted Axis A replication summaries after the 10-19 exclusion of GSE48350 GSM1176257
# (DL pooled, REML + Hartung-Knapp, leave-one-cohort-out, common members, diagnosis/infection restrictions, member directions).
# Functions copied from AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R; only the sample set changes.
# Run from project root: LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_replication_after_exclusion_20261019/code/01_recompute.R
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
gene_eff <- function(p) { g <- intersect(axisA, rownames(p$loge)); z <- t(scale(t(p$loge[g, , drop = FALSE]))); B <- qr.coef(qr(p$X), t(z)); res <- t(z) - p$X %*% B
  df <- nrow(p$X) - ncol(p$X); se <- sqrt(colSums(res^2) / df * solve(crossprod(p$X))["disease", "disease"]); data.table(cohort = p$id, entrez_id = g, beta = B["disease", ], se = se) }
dl <- function(y, s, method = "DL", test = "z") { f <- rma(yi = y, sei = s, method = method, test = test); data.table(effect_SD = as.numeric(f$b), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2) }
pp <- lapply(coh[VAL], prep); sc <- rbindlist(lapply(VAL, function(k) cbind(cohort = k, score_eff(pp[[k]])))); a <- sc[cohort %in% AD]
S <- list(); add <- function(n, r) S[[n]] <<- cbind(analysis = n, r)
add("DL (primary)", dl(a$effect_SD, a$se)); add("REML + Hartung-Knapp", dl(a$effect_SD, a$se, "REML", "knha"))
loo <- rbindlist(lapply(AD, function(k) cbind(dropped = k, dl(a[cohort != k]$effect_SD, a[cohort != k]$se))))
common <- Reduce(intersect, lapply(pp[VAL], function(p) intersect(axisA, rownames(p$R)))); cm <- rbindlist(lapply(AD, function(k) cbind(cohort = k, score_eff(pp[[k]], common))))
add(sprintf("common members (%d)", length(common)), dl(cm$effect_SD, cm$se))
swap <- function(k, r) { t <- copy(a); t[cohort == k, c("effect_SD", "se") := .(r$effect_SD, r$se)]; dl(t$effect_SD, t$se) }
add("GSE84422 probable/definite only", swap("GSE84422", score_eff(prep(coh$GSE84422, function(d) d$category != "Possible AD"))))
add("GSE313307 non-infected only", swap("GSE313307", score_eff(prep(coh$GSE313307, function(d) d$infection == "0", setdiff(coh$GSE313307$covars, "infection")))))
add("without GSE313307", dl(a[cohort != "GSE313307"]$effect_SD, a[cohort != "GSE313307"]$se))
S <- rbindlist(S, fill = TRUE)
ge <- rbindlist(lapply(pp, gene_eff)); pg <- ge[cohort %in% AD, if (.N >= 2) dl(beta, se) else data.table(effect_SD = beta, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, I2 = NA_real_), by = entrez_id]
lu <- ge[cohort == "GSE47460", .(entrez_id, lung = beta, lung_low = beta - qnorm(.975) * se, lung_high = beta + qnorm(.975) * se)]
m <- merge(lu, pg[, .(entrez_id, ad = effect_SD, ad_low = ci_low, ad_high = ci_high)], by = "entrez_id")
dir <- data.table(both_up = sum(m$lung > 0 & m$ad > 0), both_down = sum(m$lung < 0 & m$ad < 0), opposite = sum(sign(m$lung) != sign(m$ad)), both_CI_above_0 = sum(m$lung_low > 0 & m$ad_low > 0), n = nrow(m))
fwrite(sc, file.path(P(O), "results/by_cohort.tsv"), sep = "\t"); fwrite(S, file.path(P(O), "results/pooled_sensitivity.tsv"), sep = "\t")
fwrite(loo, file.path(P(O), "results/leave_one_cohort_out.tsv"), sep = "\t"); fwrite(dir, file.path(P(O), "results/member_direction_counts.tsv"), sep = "\t")
print(sc[, .(cohort, effect_SD = round(effect_SD, 3), ci = sprintf("%.2f to %.2f", ci_low, ci_high), p = signif(p, 3), n_case, n_control)])
print(S[, .(analysis, effect_SD = round(effect_SD, 3), ci = sprintf("%.2f to %.2f", ci_low, ci_high), p = signif(p, 3), I2 = round(I2, 1))])
print(loo[, .(dropped, effect_SD = round(effect_SD, 3), p = signif(p, 3))]); print(dir); cat("AD n:", sum(a$n_case), "vs", sum(a$n_control), "\n")
