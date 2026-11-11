# Follow-up of the duplicate-donor check (03_):
# (1) expression-based sex check in every cohort (XIST vs Y-linked genes RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY);
# (2) impact of the byte-identical GSE48350 CEL pair (GSM1176257 AD = GSM350078 control) and of the two
#     near-identical GSE47460 control pairs flagged by doppelgangR, on the Axis A effect and the AD pool;
# (3) expression similarity of the lung metadata-fingerprint pairs (GSE76925 vs GSE47460).
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/05_followup_checks.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages(library(metafor))
OUT <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
coh <- readRDS(P(file.path(OUT, "source_data/cohorts.rds")))$cohorts
AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
src02 <- parse(P(file.path(OUT, "code/02_member_effects_and_sensitivity.R")))
for (e in src02) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% c("prep", "score_eff", "dl")) eval(e)

# ---------- (1) sex check ----------
ygenes <- c(RPS4Y1 = "6192", DDX3Y = "8653", KDM5D = "8284", UTY = "7404", EIF1AY = "9086"); xist <- "7503"
sexchk <- rbindlist(lapply(names(coh), function(k) {
  c <- coh[[k]]; x <- c$expr[, c$dat$sample, drop = FALSE]
  if (c$kind == "counts") x <- cpm(DGEList(round(x)), log = TRUE, prior.count = 2)
  y <- intersect(ygenes, rownames(x)); if (!length(y)) return(data.table(cohort = k, note = "Y-linked genes not measured"))
  sc <- colMeans(t(scale(t(x[y, , drop = FALSE]))))                                        # high = male
  if (xist %in% rownames(x)) sc <- sc - as.numeric(scale(x[xist, ]))
  ann <- toupper(substr(trimws(as.character(c$dat$sex)), 1, 1)); ann <- fifelse(ann %in% c("M", "1"), "M", fifelse(ann %in% c("F", "0", "2"), "F", NA_character_))
  if (k == "GSE47460") ann <- fifelse(grepl("Male", c$dat$sex) & !grepl("Female", c$dat$sex), "M", "F")
  km <- kmeans(sc, centers = c(min(sc), max(sc))); pred <- fifelse(km$cluster == 2, "M", "F")
  sep <- diff(tapply(sc, km$cluster, mean)) / sqrt(mean(tapply(sc, km$cluster, var)))     # cluster separation in pooled-SD units
  data.table(cohort = k, sample = c$dat$sample, disease = c$dat$disease, annotated = ann, predicted = pred, sex_score = sc,
             y_genes = length(y), xist_used = xist %in% rownames(x), separation = as.numeric(sep))
}), fill = TRUE)
sexchk[, informative := separation >= 4]                                                    # bimodal: two clearly separated clusters
mis <- sexchk[informative == TRUE & annotated != predicted]
fwrite(sexchk, P(file.path(OUT, "results/sex_check_all_samples.tsv")), sep = "\t"); fwrite(mis, P(file.path(OUT, "results/sex_check_mismatches.tsv")), sep = "\t")
print(sexchk[!is.na(predicted), .(n = .N, separation = round(separation[1], 1), informative = informative[1], xist = xist_used[1], mismatches = sum(annotated != predicted)), by = cohort]); print(mis[, .(cohort, sample, disease, annotated, predicted, sex_score = round(sex_score, 2))])

# ---------- (2) impact of duplicated / near-identical profiles ----------
base <- rbindlist(lapply(c("GSE47460", AD), function(k) cbind(cohort = k, score_eff(prep(coh[[k]])))))
pool <- function(tab) dl(tab$effect_SD, tab$se)
imp <- list()
imp[["primary"]] <- cbind(analysis = "primary", cohort = "AD pooled", pool(base[cohort %in% AD]))
for (drop in list(c("GSM1176257"), c("GSM350078"), c("GSM1176257", "GSM350078"))) {
  s <- score_eff(prep(coh$GSE48350, function(d) !d$sample %in% drop)); tab <- copy(base[cohort %in% AD]); tab[cohort == "GSE48350", c("effect_SD", "se") := .(s$effect_SD, s$se)]
  lab <- paste("GSE48350 without", paste(drop, collapse = "+"))
  imp[[lab]] <- rbind(cbind(analysis = lab, cohort = "GSE48350", s), cbind(analysis = lab, cohort = "AD pooled", pool(tab)), fill = TRUE)
}
nearid <- c("LT058156", "LT075094", "LT115840", "LT079487")
s <- score_eff(prep(coh$GSE47460, function(d) !d$sample %in% nearid)); imp[["GSE47460 without both near-identical control pairs"]] <- cbind(analysis = "GSE47460 without 4 near-identical controls", cohort = "GSE47460", s)
sx <- rbindlist(lapply(c("GSE47460", AD), function(k) { bad <- mis[cohort == k]$sample
  cbind(cohort = k, excluded = length(bad), score_eff(prep(coh[[k]], function(d) !d$sample %in% bad))) }))
imp[["exclude sex-discordant samples"]] <- rbind(cbind(analysis = "exclude sex-discordant samples", sx),
  cbind(analysis = "exclude sex-discordant samples", cohort = "AD pooled", pool(sx[cohort %in% AD])), fill = TRUE)
imp_tab <- rbindlist(imp, fill = TRUE); imp_tab <- rbind(cbind(analysis = "primary", base), imp_tab, fill = TRUE)
fwrite(imp_tab, P(file.path(OUT, "results/duplicate_impact.tsv")), sep = "\t")
print(imp_tab[, .(analysis, cohort, effect_SD = round(effect_SD, 3), ci_low = round(ci_low, 3), ci_high = round(ci_high, 3), p = signif(p, 3), n_case, n_control)])

# ---------- (3) lung fingerprint pairs: expression similarity rank ----------
fp <- fread(P(file.path(OUT, "results/duplicate_fingerprint_pairs.tsv")))[group == "lung"]
a <- coh$GSE76925$expr[, coh$GSE76925$dat$sample]; b <- coh$GSE47460$expr[, coh$GSE47460$dat$sample]
g <- intersect(rownames(a), rownames(b)); za <- t(scale(t(a[g, ]))); zb <- t(scale(t(b[g, ]))); ok <- rowSums(!is.finite(za)) == 0 & rowSums(!is.finite(zb)) == 0
cc <- cor(za[ok, ], zb[ok, ])
fpx <- rbindlist(lapply(seq_len(nrow(fp)), function(i) { s76 <- if (fp$cohort_1[i] == "GSE76925") fp$sample_1[i] else fp$sample_2[i]; s47 <- if (fp$cohort_1[i] == "GSE47460") fp$sample_1[i] else fp$sample_2[i]
  r <- cc[s76, ]; data.table(GSE76925 = s76, GSE47460 = s47, age = fp$age[i], sex = fp$sex[i], fev1_pp = fp$fev1_pp[i], cor = r[s47],
    rank_among_GSE47460 = rank(-r)[s47], n_GSE47460 = length(r), best_other_cor = max(r[names(r) != s47])) }))
fwrite(fpx, P(file.path(OUT, "results/duplicate_lung_fingerprint_expression.tsv")), sep = "\t"); print(fpx)
writeLines(capture.output(sessionInfo()), P(file.path(OUT, "qa/sessionInfo_05.txt")))
