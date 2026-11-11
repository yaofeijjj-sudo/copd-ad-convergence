# Member-level disease effects, leave-one-member-out, replication sensitivity analyses and member coverage audit.
# Filtering, scoring and the score model are copied from fit_comparison() in
# AxisA_singscore_official_rerun_20261004/code/00_functions.R; check A below reproduces the adopted effects.
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
OUT <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
src <- readRDS(P(file.path(OUT, "source_data/cohorts.rds"))); coh <- src$cohorts
ax <- fread(P("03_candidate_analysis/WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisA_confirmed_genes.csv"), colClasses = "character")
AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
VAL <- c("GSE47460", AD)

prep <- function(c, keep = NULL, covars = c$covars) {   # = first half of fit_comparison()
  dat <- copy(c$dat); if (!is.null(keep)) dat <- dat[keep(dat)]
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  X <- model.matrix(as.formula(paste("~", rhs)), dat); stopifnot(qr(X)$rank == ncol(X))
  e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") {
    y <- DGEList(round(e)); keep_g <- filterByExpr(y, group = dat$disease); y <- calcNormFactors(y[keep_g, , keep.lib.sizes = FALSE])
    loge <- cpm(y, log = TRUE, prior.count = 2); fit <- eBayes(lmFit(voom(y, X), X))
  } else { loge <- e; fit <- eBayes(lmFit(loge, X)) }
  list(id = c$id, loge = loge, dat = dat, covars = covars, rhs = rhs, X = X, R = singscore::rankGenes(loge), t = fit$t[, "disease"])
}
score_eff <- function(p, genes = axisA) {
  up <- rownames(p$R)[rownames(p$R) %in% genes]
  d <- copy(p$dat); d$score_z <- as.numeric(scale(singscore::simpleScore(p$R, upSet = up)$TotalScore))
  f <- lm(as.formula(paste("score_z ~", p$rhs)), d); co <- summary(f)$coefficients["disease", ]
  data.table(effect_SD = co[[1]], se = co[[2]], ci_low = confint(f)["disease", 1], ci_high = confint(f)["disease", 2], p = co[[4]],
             n_case = sum(d$disease == 1), n_control = sum(d$disease == 0), members_scored = length(up))
}
gene_eff <- function(p) {                                 # standardised expression ~ same model as the score
  g <- intersect(axisA, rownames(p$loge)); z <- t(scale(t(p$loge[g, , drop = FALSE])))
  qrX <- qr(p$X); B <- qr.coef(qrX, t(z)); res <- t(z) - p$X %*% B; df <- nrow(p$X) - ncol(p$X)
  s2 <- colSums(res^2) / df; se <- sqrt(s2 * solve(crossprod(p$X))["disease", "disease"])
  data.table(cohort = p$id, entrez_id = g, beta_SD = B["disease", ], se = se, p = 2 * pt(-abs(B["disease", ] / se), df), moderated_t = p$t[g])
}
dl <- function(y, s, method = "DL", test = "z") { f <- rma(yi = y, sei = s, method = method, test = test)
  data.table(effect_SD = as.numeric(f$b), se = f$se, ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2, tau2 = f$tau2, k = f$k) }

pp <- lapply(coh, prep)
# ---------- A. reproduce the adopted score effects ----------
pub <- rbind(fread(P("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/results/AD_reference_meta.tsv"))[id %in% AD, .(id, effect_SD, axisA_covered, random_empirical_p, camera_p)],
             fread(P("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/results/axisA_specificity_results.tsv"))[id == "GSE47460-GPL14550 COPD", .(id = "GSE47460", effect_SD, axisA_covered, random_empirical_p, camera_p)])
sc <- rbindlist(lapply(names(pp), function(k) cbind(cohort = k, score_eff(pp[[k]]))))
chk <- merge(sc[, .(cohort, effect_SD, members_scored)], pub[, .(cohort = id, effect_pub = effect_SD, cov_pub = axisA_covered)], by = "cohort")
chk[, abs_diff := abs(effect_SD - effect_pub)]; print(chk); stopifnot(max(chk$abs_diff) < 1e-8, all(chk$members_scored == chk$cov_pub))
fwrite(chk, P(file.path(OUT, "qa/check_reproduce_adopted_effects.tsv")), sep = "\t")

# ---------- B. member-level effects ----------
ge <- rbindlist(lapply(pp, gene_eff))
pool_g <- ge[cohort %in% AD, if (.N >= 2) dl(beta_SD, se) else data.table(effect_SD = beta_SD, se = se, ci_low = NA_real_, ci_high = NA_real_, p = p, I2 = NA_real_, tau2 = NA_real_, k = .N), by = entrez_id]
wide <- dcast(ge, entrez_id ~ cohort, value.var = "beta_SD")
mem <- merge(ax[, .(entrez_id, symbol = current_symbol)], wide, by = "entrez_id", all.x = TRUE)
mem <- merge(mem, pool_g[, .(entrez_id, AD_pooled = effect_SD, AD_pooled_low = ci_low, AD_pooled_high = ci_high, AD_pooled_p = p, AD_k = k, AD_I2 = I2)], by = "entrez_id", all.x = TRUE)
mem <- merge(mem, ge[cohort == "GSE47460", .(entrez_id, GSE47460_low = beta_SD - qnorm(.975) * se, GSE47460_high = beta_SD + qnorm(.975) * se, GSE47460_p = p)], by = "entrez_id", all.x = TRUE)
mem <- merge(mem, src$kme, by = "entrez_id", all.x = TRUE)
mem[, class := fifelse(is.na(GSE47460) | is.na(AD_pooled), "not measurable",
               fifelse(GSE47460 > 0 & AD_pooled > 0, "both up", fifelse(GSE47460 < 0 & AD_pooled < 0, "both down", "opposite")))]
fwrite(ge, P(file.path(OUT, "results/member_effects_long.tsv")), sep = "\t")
fwrite(mem[order(-(GSE47460 + AD_pooled))], P(file.path(OUT, "results/member_effects_wide.tsv")), sep = "\t")
cnt <- rbind(
  mem[, .(comparison = "independent: GSE47460 vs AD pooled", class = class)][, .N, by = .(comparison, class)],
  mem[!is.na(GSE76925) & !is.na(GSE33000), .(comparison = "discovery: GSE76925 vs GSE33000",
      class = fifelse(GSE76925 > 0 & GSE33000 > 0, "both up", fifelse(GSE76925 < 0 & GSE33000 < 0, "both down", "opposite")))][, .N, by = .(comparison, class)],
  mem[class != "not measurable", .(comparison = "independent, both 95% CI exclude 0",
      class = fifelse(GSE47460_low > 0 & AD_pooled_low > 0, "both up (CI > 0)", fifelse(GSE47460_high < 0 & AD_pooled_high < 0, "both down (CI < 0)", "other")))][, .N, by = .(comparison, class)])
fwrite(cnt, P(file.path(OUT, "results/member_direction_counts.tsv")), sep = "\t"); print(cnt)

# ---------- C. leave-one-member-out ----------
loo <- rbindlist(lapply(axisA, function(g) {
  r <- rbindlist(lapply(VAL, function(k) if (g %in% rownames(pp[[k]]$R)) cbind(cohort = k, score_eff(pp[[k]], setdiff(axisA, g))) else NULL))
  if (!nrow(r)) return(NULL)
  ad <- r[cohort %in% AD]; pooled <- if (nrow(ad) == length(AD)) dl(ad$effect_SD, ad$se) else dl(c(ad$effect_SD, sc[cohort %in% setdiff(AD, ad$cohort)]$effect_SD), c(ad$se, sc[cohort %in% setdiff(AD, ad$cohort)]$se))
  data.table(entrez_id = g, GSE47460 = if ("GSE47460" %in% r$cohort) r[cohort == "GSE47460"]$effect_SD else sc[cohort == "GSE47460"]$effect_SD,
             AD_pooled = pooled$effect_SD, AD_pooled_p = pooled$p)
}))
base_ad <- dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se)
loo[, `:=`(delta_GSE47460 = GSE47460 - sc[cohort == "GSE47460"]$effect_SD, delta_AD = AD_pooled - base_ad$effect_SD)]
loo <- merge(ax[, .(entrez_id, symbol = current_symbol)], loo, by = "entrez_id")
fwrite(loo[order(delta_AD)], P(file.path(OUT, "results/leave_one_member_out.tsv")), sep = "\t")
cat("LOO GSE47460 range:", round(range(loo$GSE47460), 3), " AD pooled range:", round(range(loo$AD_pooled), 3), " max AD p:", signif(max(loo$AD_pooled_p), 3), "\n")

# ---------- D. replication sensitivity ----------
sens <- list(); add <- function(name, rows) sens[[name]] <<- cbind(analysis = name, rows)
add("primary (DL)", base_ad)
add("REML", dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se, "REML"))
add("REML + Hartung-Knapp", dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se, "REML", "knha"))
add("DL + Hartung-Knapp", dl(sc[cohort %in% AD]$effect_SD, sc[cohort %in% AD]$se, "DL", "knha"))
common <- Reduce(intersect, lapply(pp[VAL], function(p) intersect(axisA, rownames(p$R))))
cm <- rbindlist(lapply(VAL, function(k) cbind(cohort = k, score_eff(pp[[k]], common))))
add(sprintf("common measurable members (%d)", length(common)), dl(cm[cohort %in% AD]$effect_SD, cm[cohort %in% AD]$se))
p84 <- prep(coh$GSE84422, function(d) d$category != "Possible AD"); s84 <- score_eff(p84)
p313 <- prep(coh$GSE313307, function(d) d$infection == "0", setdiff(coh$GSE313307$covars, "infection")); s313 <- score_eff(p313)
swap <- function(tab, k, row) { tab <- copy(tab); tab[cohort == k, c("effect_SD", "se") := .(row$effect_SD, row$se)]; tab }
add("GSE84422 probable/definite AD only", dl(swap(sc[cohort %in% AD], "GSE84422", s84)$effect_SD, swap(sc[cohort %in% AD], "GSE84422", s84)$se))
add("GSE313307 non-infected donors only", dl(swap(sc[cohort %in% AD], "GSE313307", s313)$effect_SD, swap(sc[cohort %in% AD], "GSE313307", s313)$se))
add("without GSE313307", dl(sc[cohort %in% setdiff(AD, "GSE313307")]$effect_SD, sc[cohort %in% setdiff(AD, "GSE313307")]$se))
sens_tab <- rbindlist(sens, fill = TRUE)
cohort_sens <- rbind(cbind(analysis = "primary", sc[cohort %in% VAL]), cbind(analysis = "common members", cm),
                     cbind(analysis = "GSE84422 probable/definite", cohort = "GSE84422", s84), cbind(analysis = "GSE313307 non-infected", cohort = "GSE313307", s313), fill = TRUE)
# AD x infection interaction in GSE313307 (all four cells >= 37 donors)
p3 <- prep(coh$GSE313307); d3 <- copy(p3$dat); d3$score_z <- as.numeric(scale(singscore::simpleScore(p3$R, upSet = rownames(p3$R)[rownames(p3$R) %in% axisA])$TotalScore))
f3 <- lm(score_z ~ disease * infection + age + sex + plate, d3); co3 <- summary(f3)$coefficients
inter <- data.table(term = rownames(co3), estimate = co3[, 1], se = co3[, 2], p = co3[, 4])[grepl("disease", term)]
cells <- d3[, .N, by = .(AD = disease, infection)][order(AD, infection)]
pin <- prep(coh$GSE313307, function(d) d$infection == "1", setdiff(coh$GSE313307$covars, "infection")); sin <- score_eff(pin)
strata <- rbind(cbind(stratum = "no infection", s313), cbind(stratum = "infection", sin))
fwrite(sens_tab, P(file.path(OUT, "results/replication_sensitivity_pooled.tsv")), sep = "\t")
fwrite(cohort_sens, P(file.path(OUT, "results/replication_sensitivity_by_cohort.tsv")), sep = "\t")
fwrite(rbind(cbind(part = "interaction model", inter), fill = TRUE), P(file.path(OUT, "results/GSE313307_infection_interaction.tsv")), sep = "\t")
fwrite(cells, P(file.path(OUT, "results/GSE313307_cells.tsv")), sep = "\t"); fwrite(strata, P(file.path(OUT, "results/GSE313307_strata.tsv")), sep = "\t")
fwrite(data.table(entrez_id = common, symbol = ax$current_symbol[match(common, ax$entrez_id)]), P(file.path(OUT, "results/common_measurable_members.tsv")), sep = "\t")
print(sens_tab[, .(analysis, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), I2 = round(I2, 1))])
print(inter); print(cells); print(strata[, .(stratum, effect_SD, ci_low, ci_high, p, n_case, n_control)])

# ---------- E. member coverage and reasons ----------
annot <- list(GSE76925 = NULL, GSE33000 = c("03_candidate_analysis/WGCNA_Entrez_training_only_20260930/references/GPL4372.annot.gz", "Gene ID"),
  GSE47460 = c("02_data/raw/platform_annotations/GPL14550_annot_table.csv", "GENE"), GSE84422 = c("02_data/raw/platform_annotations/GPL96_annot_table.csv", "Gene ID"),
  GSE118553 = c("02_data/raw/platform_annotations/GPL10558.annot.gz", "Gene ID"), GSE48350 = c("02_data/raw/platform_annotations/GPL570_annot_table.csv", "Gene ID"),
  GSE36980 = c("02_data/raw/disease_specificity_tissue_20261003/GPL6244.annot.gz", "Gene ID"))
on_platform <- function(a) { x <- if (grepl("\\.annot\\.gz$", a[1])) fread(cmd = paste("gzip -cd", shQuote(P(a[1])), "| sed -n '/^ID\t/,$p' | grep -v '^!platform_table_end'"), colClasses = "character", select = c("ID", a[2]))
  else fread(P(a[1]), colClasses = "character", select = c("ID", a[2])); g <- x[[2]]; unique(current_id(g[!is.na(g) & nzchar(g) & !grepl("/", g)])) }
cov <- rbindlist(lapply(names(pp), function(k) {
  c <- coh[[k]]; inraw <- axisA %in% rownames(c$expr); inana <- axisA %in% rownames(pp[[k]]$R)
  reason <- if (c$kind == "counts") fifelse(inana, "measured", fifelse(inraw, "removed by filterByExpr (low expression)", "no unique Ensembl-to-Entrez mapping"))
    else if (k == "GSE76925") fifelse(inana, "measured", "absent from discovery matrix (probe/symbol-to-Entrez mapping)")
    else { op <- axisA %in% on_platform(annot[[k]]); fifelse(inana, "measured", fifelse(op, "probe present, removed (non-finite values)", "no probe with a unique Entrez ID on platform")) }
  data.table(cohort = k, entrez_id = axisA, symbol = ax$current_symbol[match(axisA, ax$entrez_id)], reason = reason)
}))
fwrite(cov, P(file.path(OUT, "results/member_coverage_long.tsv")), sep = "\t")
covsum <- cov[, .(measured = sum(reason == "measured"), missing = sum(reason != "measured"),
                  missing_detail = paste(sprintf("%s: %s", unique(reason[reason != "measured"]), sapply(unique(reason[reason != "measured"]), function(r) paste(symbol[reason == r], collapse = ","))), collapse = " | ")), by = cohort]
fwrite(covsum, P(file.path(OUT, "results/member_coverage_summary.tsv")), sep = "\t"); print(covsum[, .(cohort, measured, missing)])

# ---------- F. cohort audit table ----------
audit <- rbindlist(lapply(names(coh), function(k) { c <- coh[[k]]; s <- sc[cohort == k]
  data.table(cohort = k, organ = c$organ, role = c$role, source = c$source, region = c$region, diagnosis = c$diagnosis,
    data_type = c(array = "microarray (log scale)", counts = "RNA-seq counts -> TMM logCPM")[[c$kind]],
    n_case = s$n_case, n_control = s$n_control, covariates = paste(pp[[k]]$covars, collapse = ";"),
    score = "singscore simpleScore, z-scored within the analysed samples", members_measured = s$members_scored) }))
audit <- merge(audit, pub[, .(cohort = id, random_set_empirical_p = random_empirical_p, camera_p)], by = "cohort", all.x = TRUE)
fwrite(audit, P(file.path(OUT, "results/cohort_audit.tsv")), sep = "\t")
writeLines(capture.output(sessionInfo()), P(file.path(OUT, "qa/sessionInfo_02.txt")))
