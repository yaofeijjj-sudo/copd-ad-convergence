# Duplicate-donor check across cohorts (DESIGN.md section 3).
# (1) metadata fingerprints: identical age + sex + diagnosis (lung: + FEV1 % predicted) across cohorts;
# (2) doppelgangR 1.38.0 (Waldron et al. 2016 JNCI), default arguments except cache.dir = NULL and serial execution.
# Run from project root:
#   R_LIBS=_workspace/r_libs_doppel:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/03_duplicate_donors.R
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(Biobase); library(doppelgangR); library(BiocParallel) })
OUT <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
coh <- readRDS(file.path(OUT, "source_data/cohorts.rds"))$cohorts
e1 <- new.env(); load("WGCNA_全部图件/复现资料/analysis/01_lung_wgcna/results/COPD_WGCNA/copd_wgcna_final_env.RData", envir = e1)
e2 <- new.env(); load("01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/ad_wgcna_final_env.RData", envir = e2)
sex_std <- function(x) { x <- toupper(substr(trimws(as.character(x)), 1, 1)); fifelse(x %in% c("M", "1"), "M", fifelse(x %in% c("F", "0", "2"), "F", NA_character_)) }
pheno <- function(k) { d <- copy(coh[[k]]$dat)
  if (k == "GSE76925") { d$age <- e1$clinical_train$age; d$fev1_pp <- e1$clinical_train$FEV1_pp }
  if (k == "GSE33000") d$age <- e2$clinical_train$Age
  if (k == "GSE47460") d$sex <- sub("^[12]-", "", d$sex)
  data.table(cohort = k, sample = d$sample, age = round(as.numeric(d$age)), sex = sex_std(d$sex), dx = fifelse(d$disease == 1, "case", "control"),
             fev1_pp = if ("fev1_pp" %in% names(d)) round(as.numeric(d$fev1_pp)) else NA_real_) }
loge <- function(k) { c <- coh[[k]]; e <- c$expr[, c$dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); y <- calcNormFactors(y[filterByExpr(y, group = c$dat$disease), , keep.lib.sizes = FALSE]); cpm(y, log = TRUE, prior.count = 2) } else e }
groups <- list(lung = c("GSE76925", "GSE47460"), brain = c("GSE33000", "GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"))

fp_all <- list(); dg_all <- list()
for (g in names(groups)) {
  ph <- rbindlist(lapply(groups[[g]], pheno))
  key <- if (g == "lung") c("age", "sex", "fev1_pp") else c("age", "sex", "dx")
  m <- merge(ph, ph, by = key, allow.cartesian = TRUE, suffixes = c("_1", "_2"))[cohort_1 < cohort_2]
  m <- m[!is.na(age) & !is.na(sex)]; if (g == "lung") m <- m[!is.na(fev1_pp)]
  fp_all[[g]] <- cbind(group = g, m)
  esets <- lapply(groups[[g]], function(k) { x <- loge(k); x <- x[rowSums(!is.finite(x)) == 0, ]
    p <- as.data.frame(ph[cohort == k, .(age, sex, dx)]); rownames(p) <- ph[cohort == k]$sample
    ExpressionSet(x, phenoData = AnnotatedDataFrame(p[colnames(x), , drop = FALSE])) })
  names(esets) <- groups[[g]]
  set.seed(20261018)
  res <- doppelgangR(esets, cache.dir = NULL, BPPARAM = SerialParam(), verbose = FALSE)
  s <- as.data.table(summary(res)); dg_all[[g]] <- cbind(group = g, s)
  cat(g, ": fingerprint pairs =", nrow(m), "; doppelgangR flagged pairs =", nrow(s), "\n")
}
fp <- rbindlist(fp_all, fill = TRUE); dg <- rbindlist(dg_all, fill = TRUE)
fwrite(fp, file.path(OUT, "results/duplicate_fingerprint_pairs.tsv"), sep = "\t")
fwrite(dg, file.path(OUT, "results/duplicate_doppelgangR_summary.tsv"), sep = "\t")
print(fp[, .N, by = .(group, cohort_1, cohort_2)]); print(dg)
writeLines(capture.output(sessionInfo()), file.path(OUT, "qa/sessionInfo_03.txt"))
