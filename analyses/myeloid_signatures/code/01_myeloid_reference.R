# Advice 2.6: Axis A vs (1) random sets matched on cohort expression and single-cell myeloid specificity,
# (2) pre-specified independent myeloid identity / function signatures. Score model = adopted replication model.
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_deconv_20261002:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_myeloid_reference_20261019/code/01_myeloid_reference.R
source("03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/00_functions.R")
suppressPackageStartupMessages({ library(metafor); library(AnnotationDbi); library(org.Hs.eg.db); library(BRETIGEA); library(msigdbr) })
OUT <- "03_candidate_analysis/AxisA_myeloid_reference_20261019"; RB <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
CX <- "03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/upload"
coh <- readRDS(P(file.path(RB, "source_data/cohorts.rds")))$cohorts
EXCL <- fread(P(file.path(RB, "results/adopted_sample_exclusions_20261019.tsv")))
AD <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"); VAL <- c("GSE47460", AD)
sym2ent <- function(sym) { m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(sym), keytype = "SYMBOL", columns = "ENTREZID")
  m <- m[!is.na(m$ENTREZID), ]; m <- m[!(m$SYMBOL %in% m$SYMBOL[duplicated(m$SYMBOL)]), ]; setNames(m$ENTREZID, m$SYMBOL) }

prep <- function(c) {
  dat <- copy(c$dat)[!sample %in% EXCL[cohort == c$id]$sample]; covars <- c$covars
  dat <- dat[complete.cases(dat[, c("disease", covars), with = FALSE])]
  for (v in covars) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  covars <- covars[vapply(covars, function(v) length(unique(dat[[v]])) > 1, TRUE)]
  rhs <- paste("disease", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")
  e <- c$expr[, dat$sample, drop = FALSE]
  if (c$kind == "counts") { y <- DGEList(round(e)); keep_g <- filterByExpr(y, group = dat$disease)
    y <- calcNormFactors(y[keep_g, , keep.lib.sizes = FALSE]); loge <- cpm(y, log = TRUE, prior.count = 2) } else loge <- e
  list(id = c$id, dat = dat, rhs = rhs, R = singscore::rankGenes(loge), mean_expr = rowMeans(loge))
}
zscore <- function(p, genes) as.numeric(scale(singscore::simpleScore(p$R, upSet = rownames(p$R)[rownames(p$R) %in% genes])$TotalScore))
eff <- function(p, z, extra = NULL) { d <- copy(p$dat); d$score_z <- z; rhs <- p$rhs; if (!is.null(extra)) { d$ref_z <- extra; rhs <- paste(rhs, "+ ref_z") }
  f <- lm(as.formula(paste("score_z ~", rhs)), d); co <- summary(f)$coefficients["disease", ]
  vif <- if (is.null(extra)) NA_real_ else 1 / (1 - summary(lm(as.formula(paste("ref_z ~", p$rhs)), d))$r.squared)
  data.table(effect_SD = co[[1]], se = co[[2]], ci_low = confint(f)["disease", 1], ci_high = confint(f)["disease", 2], p = co[[4]], vif_ref = vif) }
dl <- function(y, s) { f <- rma(yi = y, sei = s, method = "DL"); data.table(effect_SD = as.numeric(f$b), se = f$se, ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval, I2 = f$I2) }

pp <- lapply(coh[VAL], prep)
ref <- rbind(fread(P(file.path(RB, "results/replication_sensitivity_by_cohort.tsv")))[analysis == "primary" & cohort != "GSE48350", .(cohort, effect_SD)],
             fread(P(file.path(RB, "results/duplicate_impact.tsv")))[analysis == "GSE48350 without GSM1176257" & cohort == "GSE48350", .(cohort, effect_SD)])
base <- rbindlist(lapply(VAL, function(k) cbind(cohort = k, eff(pp[[k]], zscore(pp[[k]], axisA)))))
chk <- merge(base[, .(cohort, mine = effect_SD)], ref, by = "cohort"); stopifnot(nrow(chk) == 7, max(abs(chk$mine - chk$effect_SD)) < 1e-8)
base_pool <- dl(base[cohort %in% AD]$effect_SD, base[cohort %in% AD]$se)

# ---------- myeloid specificity from the independent single-cell references ----------
spec_of <- function(f, target) {
  hdr <- strsplit(readLines(f, n = 1), "\t")[[1]][-1]; x <- fread(f); g <- x[[1]]; m <- as.matrix(x[, -1]); rm(x)
  m <- sweep(m, 2, colSums(m), "/") * 1e6
  cm <- sapply(unique(hdr), function(k) rowMeans(m[, hdr == k, drop = FALSE])); rm(m); gc()
  s <- cm[, target] / rowSums(cm); s[!is.finite(s)] <- NA
  e <- sym2ent(g); keep <- g %in% names(e); out <- setNames(s[keep], e[g[keep]]); out[!duplicated(names(out))] }
spec <- list(lung = spec_of(P(file.path(CX, "lung_scref_A.txt")), "Myeloid"), brain = spec_of(P(file.path(CX, "brain_scref_A.txt")), "Immune"))
fwrite(data.table(organ = rep(names(spec), lengths(spec)), entrez_id = unlist(lapply(spec, names)), myeloid_specificity = unlist(spec)), file.path(P(OUT), "results/myeloid_specificity_by_gene.tsv"), sep = "\t")
cat("Axis A median specificity: lung", round(median(spec$lung[axisA], na.rm = TRUE), 3), "brain", round(median(spec$brain[axisA], na.rm = TRUE), 3),
    "| all genes: lung", round(median(spec$lung, na.rm = TRUE), 3), "brain", round(median(spec$brain, na.rm = TRUE), 3), "\n")

# ---------- part 1: matched null ----------
B <- 1000; set.seed(20261019)
null <- list(); qa <- list()
for (k in VAL) {
  p <- pp[[k]]; s <- spec[[if (k == "GSE47460") "lung" else "brain"]]
  g <- intersect(rownames(p$R), names(s)[!is.na(s)]); ax <- intersect(g, axisA)
  eb <- cut(rank(p$mean_expr[g]), quantile(rank(p$mean_expr[g]), 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  sb <- cut(rank(s[g]), quantile(rank(s[g]), 0:5 / 5), include.lowest = TRUE, labels = FALSE)
  cell <- setNames(paste(eb, sb), g); need <- table(cell[ax]); pool <- split(setdiff(g, axisA), cell[setdiff(g, axisA)])
  short <- sapply(names(need), function(cc) need[[cc]] > length(pool[[cc]])); stopifnot(!any(short))
  qa[[k]] <- data.table(cohort = k, axis_scored = sum(rownames(p$R) %in% axisA), axis_with_specificity = length(ax),
                        cells_used = length(need), min_pool_ratio = min(sapply(names(need), function(cc) length(pool[[cc]]) / need[[cc]])))
  ne <- vapply(seq_len(B), function(b) { set <- unlist(lapply(names(need), function(cc) sample(pool[[cc]], need[[cc]])))
    unlist(eff(p, zscore(p, set))[, .(effect_SD, se)]) }, numeric(2))
  null[[k]] <- data.table(cohort = k, b = seq_len(B), effect_SD = ne[1, ], se = ne[2, ])
}
null <- rbindlist(null); qa <- rbindlist(qa)
emp <- base[, .(cohort, axisA_effect = effect_SD)][, `:=`(null_median = sapply(cohort, function(k) median(null[cohort == k]$effect_SD)),
  empirical_p_one_sided = sapply(seq_along(cohort), function(i) (sum(null[cohort == base$cohort[i]]$effect_SD >= axisA_effect[i]) + 1) / (B + 1)))]
pooled_null <- null[cohort %in% AD, dl(effect_SD, se)$effect_SD, by = b]$V1
emp <- rbind(emp, data.table(cohort = "AD pooled (DL)", axisA_effect = base_pool$effect_SD, null_median = median(pooled_null),
  empirical_p_one_sided = (sum(pooled_null >= base_pool$effect_SD) + 1) / (B + 1)))
old <- fread(P(file.path(RB, "results/cohort_audit.tsv")))[, .(cohort, expression_matched_only_p = random_set_empirical_p)]
emp <- merge(emp, old, by = "cohort", all.x = TRUE)
fwrite(null, file.path(P(OUT), "results/matched_null_effects.tsv.gz"), sep = "\t"); fwrite(emp, file.path(P(OUT), "results/matched_null_empirical_p.tsv"), sep = "\t")
fwrite(qa, file.path(P(OUT), "qa/matched_null_qa.tsv"), sep = "\t"); print(qa); print(emp)

# ---------- part 2: independent signatures ----------
data(markers_df_human_brain)
sigs <- list(
  BRETIGEA_microglia_top250 = unname(sym2ent(markers_df_human_brain$markers[markers_df_human_brain$cell == "mic"][1:250])),
  Patir2019_core_microglia = unname(sym2ent(readLines(P("03_candidate_analysis/AxisA_auxiliary_checks_20261019/source_data/patir2019_core_symbols.txt")))),
  Travaglini2020_lung_macrophage = unique(as.character(as.data.table(msigdbr(species = "Homo sapiens", collection = "C8"))[gs_name == "TRAVAGLINI_LUNG_MACROPHAGE_CELL"]$ncbi_gene)),
  GO0006956_complement_activation = unique(get("GO:0006956", org.Hs.egGO2ALLEGS)),
  GO0006909_phagocytosis = unique(get("GO:0006909", org.Hs.egGO2ALLEGS)))
sm <- rbindlist(lapply(names(sigs), function(n) data.table(signature = n, n_genes = length(sigs[[n]]), overlap_axisA = sum(axisA %in% sigs[[n]]))))
fwrite(sm, file.path(P(OUT), "results/signature_manifest.tsv"), sep = "\t"); print(sm)
s2 <- rbindlist(lapply(names(sigs), function(n) rbindlist(lapply(VAL, function(k) {
  p <- pp[[k]]; zr <- zscore(p, sigs[[n]]); za <- zscore(p, axisA); sh <- intersect(axisA, sigs[[n]])
  za_d <- zscore(p, setdiff(axisA, sh)); zr_d <- zscore(p, setdiff(sigs[[n]], sh))
  rbind(cbind(signature = n, cohort = k, analysis = "signature alone", eff(p, zr), r_with_axisA = cor(za, zr)),
        cbind(signature = n, cohort = k, analysis = "Axis A adjusted for signature", eff(p, za, zr), r_with_axisA = cor(za, zr)),
        cbind(signature = n, cohort = k, analysis = "Axis A adjusted, shared genes removed from both", eff(p, za_d, zr_d), r_with_axisA = cor(za_d, zr_d)))
}))))
pool2 <- s2[cohort %in% AD, cbind(dl(effect_SD, se), median_r = median(r_with_axisA), max_vif = max(vif_ref, na.rm = TRUE)), by = .(signature, analysis)]
fwrite(s2, file.path(P(OUT), "results/signature_models_by_cohort.tsv"), sep = "\t"); fwrite(pool2, file.path(P(OUT), "results/signature_models_AD_pooled.tsv"), sep = "\t")
print(s2[cohort == "GSE47460", .(signature, analysis, effect_SD = round(effect_SD, 3), ci = sprintf("%.2f to %.2f", ci_low, ci_high), p = signif(p, 2), vif = round(vif_ref, 1), r = round(r_with_axisA, 2))])
print(pool2[, .(signature, analysis, effect_SD = round(effect_SD, 3), ci = sprintf("%.2f to %.2f", ci_low, ci_high), p = signif(p, 2), median_r = round(median_r, 2), max_vif = round(max_vif, 1))])
cat("base: GSE47460", round(base[cohort == "GSE47460"]$effect_SD, 3), "| AD pooled", round(base_pool$effect_SD, 3), "\n")
writeLines(capture.output(sessionInfo()), file.path(P(OUT), "qa/sessionInfo_01.txt"))
