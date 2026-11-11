# Fig 7. Axis A in peripheral blood and disease scope.
# Sources (frozen, not recomputed):
#   03_candidate_analysis/AxisA_auxiliary_checks_20261019/results/blood_camera_and_effects.tsv (a)
#   REF/disease_scope.tsv (comparison diseases; chronic hypersensitivity pneumonitis omitted, listed in S4 Table);
#   03_candidate_analysis/AxisA_disease_scope_lung_family_20261028/results/lung_family_without_CHP.tsv (lung-family BH FDR without CHP);
#   03_candidate_analysis/AxisA_replication_after_exclusion_20261019/results/{by_cohort,pooled_sensitivity}.tsv (b, target diseases)
#   03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/results/axisA_specificity_results.tsv (b, GSE47460 IPF and other ILD; Table D (e))
# Both panels use the table-forest design of the 10-05 Figure 4d (Context | Dataset | forest | effect | P).
source("06_resubmission_20261008/code/main_figures/common.R")
neg <- identity  # ASCII hyphen as in the 10-05 panels (U+2212 would start a second font subset)

# ---- a: peripheral blood --------------------------------------------------------------------------------------------
bl <- fread(P(file.path(CAND, "AxisA_auxiliary_checks_20261019/results/blood_camera_and_effects.tsv")))
ord <- c("GSE54837", "E-MTAB-5278", "GSE248493", "GSE63060", "GSE63061", "GSE248417")
stopifnot(setequal(ord, bl$dataset_id))
bl <- bl[match(ord, dataset_id)]
bl[, y := c(7, 6, 5, 3, 2, 1)]
da <- bl[, .(y, ctx = fifelse(sample_type == "PBMC", "PBMC", "Whole blood"), ds = dataset_id, est = score_diff_SD, lo = ci_low, hi = ci_high,
             col = fifelse(disease == "COPD", COL$lung, COL$brain), shape = "square", size = 2.1,
             c1 = neg(sprintf("%.2f [%.2f, %.2f]", score_diff_SD, ci_low, ci_high)), c2 = fmt_p_sup(camera_P, 2))]
ahdr <- data.table(y = c(8, 4), lab = c("COPD", "AD"), colr = c(COL$lung, COL$brain))
ytop <- 9.3; yl <- c(0.4, 9.6)
al <- ftext(da, list(list(var = "ctx", head = "Sample", x = 0.06), list(var = "ds", head = "Dataset", x = 0.74, colour = MUTED)), ahdr, ytop, yl,
            title = "Peripheral blood, age- and sex-adjusted")
af <- fplot(da, c(-0.5, 1.1), c(-0.5, 0, 0.5, 1), "Case–control difference\nin Axis A score (SD)", ytop, yl, head = "Effect")
ar <- ftext(da, list(list(var = "c1", head = "SD [95% CI]", x = 0.06), list(var = "c2", head = "CAMERA P", x = 0.66)), data.table(), ytop, yl)
pa <- al + af + ar + plot_layout(widths = c(72, 54, 64))

# ---- b: disease scope ----------------------------------------------------------------------------------------------------
ds <- fread(P(file.path(REF, "disease_scope.tsv")))
lf <- fread(P(file.path(CAND, "AxisA_disease_scope_lung_family_20261028/results/lung_family_without_CHP.tsv")))
rb <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/by_cohort.tsv")))[cohort == "GSE47460"]
rp <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/pooled_sensitivity.tsv")))[analysis == "DL (primary)"]
stopifnot(abs(ds[id == "GSE47460 COPD"]$effect_SD - rb$effect_SD) < 1e-9)
ctrl <- ds[!id %in% c("GSE47460 COPD", "GSE150910 CHP")]
# lung comparison family = IPF + ever-smokers (CHP described in S4 Table only); brain family unchanged
ctrl[id %in% lf$id, fdr := lf$fdr_without_CHP[match(id, lf$id)]]
stopifnot(abs(ctrl[id == "GSE150910 IPF"]$fdr - 0.6639) < 1e-3, abs(ctrl[id == "GSE150910 smoking"]$fdr - 0.6639) < 1e-3)
ctrl[, ctx := c("GSE150910 IPF" = "IPF", "GSE150910 smoking" = "Ever-smokers, no lung disease",
                "GSE64810 HD" = "Huntington's disease", "GSE68719 PD" = "Parkinson's disease",
                "GSE53987 SCZ" = "Schizophrenia", "GSE53987 BD" = "Bipolar disorder", "GSE53987 MDD" = "Major depressive disorder")[id]]
stopifnot(!anyNA(ctrl$ctx))
sec <- fread(P(file.path(CAND, "AxisA_disease_specificity_redesign_20261003/results/axisA_specificity_results.tsv")))[id %in% c("GSE47460-GPL14550 IPF", "GSE47460-GPL14550 otherILD")]
stopifnot(nrow(sec) == 2, abs(sec[grepl("IPF", id)]$effect_SD - 0.0148) < 1e-3, abs(sec[grepl("otherILD", id)]$effect_SD - 0.4668) < 1e-3)
db <- rbind(
  sec[, .(ctx = fifelse(grepl("IPF", id), "IPF", "Other interstitial lung disease"), dsn = "GSE47460", est = effect_SD, lo = ci_low, hi = ci_high, col = INK, shape = "square",
          c2 = paste0("P = ", fmt_p_sup(p, 2)), k = id)],
  data.table(ctx = "COPD", dsn = "GSE47460", est = rb$effect_SD, lo = rb$ci_low, hi = rb$ci_high, col = COL$target, shape = "square",
             c2 = paste0("P = ", fmt_p_sup(rb$p, 2)), k = "GSE47460 COPD"),
  ctrl[, .(ctx, dsn = sub(" .*", "", id), est = effect_SD, lo = ci_low, hi = ci_high, col = INK, shape = fifelse(fdr < 0.05, "circle", "open"),
           c2 = fmt_p_sup(fdr, 2), k = id)],
  data.table(ctx = "AD (pooled)", dsn = "6 datasets", est = rp$effect_SD, lo = rp$ci_low, hi = rp$ci_high, col = COL$target, shape = "diamond",
             c2 = paste0("P = ", fmt_p_sup(rp$p, 2)), k = "AD pooled"))
db[, y := c(13, 12, 11, 10, 9, 6, 5, 4, 3, 2, 1)[match(k, c("GSE47460 COPD", "GSE47460-GPL14550 IPF", "GSE47460-GPL14550 otherILD", "GSE150910 IPF", "GSE150910 smoking", "AD pooled",
                                                       "GSE64810 HD", "GSE68719 PD", "GSE53987 SCZ", "GSE53987 BD", "GSE53987 MDD"))]]
db[, size := 2.1][, c1 := neg(sprintf("%.2f [%.2f, %.2f]", est, lo, hi))]
bhdr <- data.table(y = c(14, 7), lab = c("Lung tissue", "Frontal cortex"))
ytop <- 15.3; yl <- c(0.4, 15.6)
bl_ <- ftext(db, list(list(var = "ctx", head = "Context", x = 0.06), list(var = "dsn", head = "Dataset", x = 0.74, colour = MUTED)), bhdr, ytop, yl,
             title = "Disease scope in tissue")
bf <- fplot(db, c(-1.25, 1.35), c(-1, 0, 1), "Case–control difference\nin Axis A score (SD)", ytop, yl, head = "Effect")
br_ <- ftext(db, list(list(var = "c1", head = "SD [95% CI]", x = 0.06), list(var = "c2", head = "P / FDR", x = 0.66)), data.table(), ytop, yl)
pb <- bl_ + bf + br_ + plot_layout(widths = c(72, 54, 64))

fwrite(da, P(file.path(SRC, "Fig7a_blood.tsv")), sep = "\t"); fwrite(db, P(file.path(SRC, "Fig7b_scope.tsv")), sep = "\t")

compose_figure(list(blk(pa, 0, 0, 190, 64, c("a_left", "a", "a_right"), tags = list(c("a", 0, 0))),
                    blk(pb, 0, 66, 190, 96, c("b_left", "b", "b_right"), tags = list(c("b", 0, 0)))),
               "Fig7", 190, 162, row_groups = list(c("a_left", "a", "a_right"), c("b_left", "b", "b_right")),
               column_groups = list(c("a", "b")),
               exemptions = list(list(panels = list("a_left", "a", "a_right"), checks = list("panel-width"), reason = "table forest: text columns and the forest differ in width by design"), list(panels = list("b_left", "b", "b_right"), checks = list("panel-width"), reason = "table forest: text columns and the forest differ in width by design")))
