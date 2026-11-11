# Axis B fixed-score replication. prep()/score_eff()/dl() copied from
# AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R; only the gene set changes.
# Check A reproduces the adopted Axis A effects with the same functions, so the pipeline is identical.
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisB_fixed_score_replication_20261019/code/01_axisB_replication.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
OUT <- "03_candidate_analysis/AxisB_fixed_score_replication_20261019"
RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
CX <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018"
dir.create(P(file.path(OUT, "results")), FALSE); dir.create(P(file.path(OUT, "qa")), FALSE)
coh <- readRDS(P(file.path(RB, "source_data/cohorts.rds")))$cohorts
axisB <- fread(P("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisB_confirmed_genes.csv"), colClasses = "character")$entrez_id
stopifnot(length(axisB) == 434, !any(axisB %in% axisA))
EXCL <- fread(P(file.path(RB, "results/adopted_sample_exclusions_20261019.tsv")))
AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"); VAL <- c("GSE47460", AD)

prep <- function(c, excl = TRUE, covars = c$covars) {
  dat <- copy(c$dat); if (excl) dat <- dat[!sample %in% EXCL[cohort == c$id]$sample]
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat); stopifnot(qr(X)$rank == ncol(X))
  e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); keep_g <- filterByExpr(y, group = dat$disease)
    y <- calcNormFactors(y[keep_g, , keep.lib.sizes = FALSE]); loge <- cpm(y, log = TRUE, prior.count = 2) } else loge <- e
  list(id = c$id, dat = dat, rhs = rhs, R = singscore::rankGenes(loge))
}
score_eff <- function(p, genes, adj = NULL) {
  up <- rownames(p$R)[rownames(p$R) %in% genes]
  d <- copy(p$dat); d$score_z <- as.numeric(scale(singscore::simpleScore(p$R, upSet = up)$TotalScore))
  rhs <- p$rhs; if (!is.null(adj)) { d <- merge(d, adj, by = "sample"); stopifnot(nrow(d) == nrow(p$dat)); rhs <- paste(rhs, "+ neuron_z") }
  f <- lm(as.formula(paste("score_z ~", rhs)), d); co <- summary(f)$coefficients["disease", ]
  data.table(effect_SD = co[[1]], se = co[[2]], ci_low = confint(f)["disease", 1], ci_high = confint(f)["disease", 2], p = co[[4]],
             n_case = sum(d$disease == 1), n_control = sum(d$disease == 0), members_scored = length(up))
}
dl <- function(y, s, method = "DL", test = "z") { f <- rma(yi = y, sei = s, method = method, test = test)
  data.table(effect_SD = as.numeric(f$b), se = f$se, ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2, k = f$k) }

pp <- lapply(coh, prep)
# ---------- check: same functions reproduce adopted Axis A effects (GSE48350 after exclusion) ----------
ref <- rbind(fread(P(file.path(RB, "results/replication_sensitivity_by_cohort.tsv")))[analysis == "primary" & cohort != "GSE48350", .(cohort, effect_SD)],
             fread(P(file.path(RB, "results/duplicate_impact.tsv")))[analysis == "GSE48350 without GSM1176257" & cohort == "GSE48350", .(cohort, effect_SD)])
chkA <- rbindlist(lapply(VAL, function(k) data.table(cohort = k, axisA_effect = score_eff(pp[[k]], axisA)$effect_SD)))
chkA <- merge(chkA, ref, by = "cohort"); chkA[, abs_diff := abs(axisA_effect - effect_SD)]; print(chkA); stopifnot(max(chkA$abs_diff) < 1e-8)
fwrite(chkA, P(file.path(OUT, "qa/check_axisA_reproduced.tsv")), sep = "\t")

# ---------- primary + descriptive ----------
sc <- rbindlist(lapply(names(pp), function(k) cbind(cohort = k, role = fifelse(k %in% VAL, "independent", "discovery (descriptive)"), score_eff(pp[[k]], axisB))))
pool <- list(); add <- function(name, rows) pool[[name]] <<- cbind(analysis = name, rows)
add("AD pooled, DL (primary)", dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se))
add("AD pooled, REML + Hartung-Knapp", dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se, "REML", "knha"))
common <- Reduce(intersect, lapply(pp[VAL], function(p) intersect(axisB, rownames(p$R))))
cm <- rbindlist(lapply(VAL, function(k) cbind(cohort = k, score_eff(pp[[k]], common))))
add(sprintf("AD pooled, common measurable members (%d)", length(common)), dl(cm[cohort %in% AD]$effect_SD, cm[cohort %in% AD]$se))
s48 <- score_eff(prep(coh$GSE48350, excl = FALSE), axisB); tmp <- copy(sc[cohort %in% AD]); tmp[cohort == "GSE48350", c("effect_SD", "se") := .(s48$effect_SD, s48$se)]
add("AD pooled, GSM1176257 kept", dl(tmp$effect_SD, tmp$se))

# ---------- secondary: neuronal proportion (CIBERSORTx version B, Excitatory + Inhibitory) ----------
neu <- rbindlist(lapply(AD, function(k) {
  x <- fread(P(file.path(CX, sprintf("results/raw/B_%s.tsv", k))), check.names = FALSE); map <- fread(P(file.path(CX, sprintf("source_data/%s_mixture_samples.tsv", k))))
  a <- data.table(sample = map$sample[match(x$Mixture, map$upload_name)], neuron = x$Excitatory + x$Inhibitory); a[, neuron_z := as.numeric(scale(neuron))]
  cbind(cohort = k, score_eff(pp[[k]], axisB, a[, .(sample, neuron_z)])) }))
add("AD pooled, + neuronal proportion (secondary)", dl(neu$effect_SD, neu$se))
pool <- rbindlist(pool, fill = TRUE)

g <- sc[cohort == "GSE47460"]; a <- pool[analysis == "AD pooled, DL (primary)"]
stable <- g$effect_SD < 0 && g$ci_high < 0 && a$effect_SD < 0 && a$ci_high < 0
verdict <- data.table(GSE47460_effect = g$effect_SD, GSE47460_ci = sprintf("%.2f to %.2f", g$ci_low, g$ci_high),
  AD_pooled_effect = a$effect_SD, AD_pooled_ci = sprintf("%.2f to %.2f", a$ci_low, a$ci_high),
  rule = "stable = both negative with 95% CI excluding 0 (DESIGN.md)", stable = stable)
fwrite(sc, P(file.path(OUT, "results/axisB_by_cohort.tsv")), sep = "\t")
fwrite(cm, P(file.path(OUT, "results/axisB_common_members_by_cohort.tsv")), sep = "\t")
fwrite(neu, P(file.path(OUT, "results/axisB_neuron_adjusted_by_cohort.tsv")), sep = "\t")
fwrite(pool, P(file.path(OUT, "results/axisB_pooled.tsv")), sep = "\t")
fwrite(verdict, P(file.path(OUT, "results/axisB_verdict.tsv")), sep = "\t")
print(sc[, .(cohort, role, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), n_case, n_control, members_scored)])
print(pool[, .(analysis, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), I2 = round(I2, 1))])
print(neu[, .(cohort, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3))]); print(verdict)
writeLines(capture.output(sessionInfo()), P(file.path(OUT, "qa/sessionInfo_01.txt")))
