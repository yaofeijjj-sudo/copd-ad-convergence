#!/usr/bin/env Rscript
# Forest plot of Axis A disease-context scope (submission version: comparator diseases only; COPD/AD target rows
# are omitted by author decision 2026-10-03 and kept in 03_draw_forest_with_targets.R for reviewer response), in the Figure 4 blood-panel style
# (R grid, Quartz PDF, Arial 10 pt, ink #262626, coral #EA5B5B).
# Run from project root: _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/03_draw_forest.R
suppressPackageStartupMessages({library(grid); library(data.table)})
BR <- "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003"
out <- file.path(BR, "figures"); dir.create(out, showWarnings = FALSE)
r <- fread(file.path(BR, "results/axisA_specificity_results.tsv"))
rows <- data.table(
  id = c(NA, "GSE47460-GPL14550 COPD", "GSE150910 CHP", "GSE150910 IPF", "GSE150910 smoking",
         NA, "AD pooled (DL random effects)", "GSE64810 HD", "GSE68719 PD", "GSE53987 SCZ", "GSE53987 BD", "GSE53987 MDD"),
  label = c("Lung tissue", "COPD", "Hypersensitivity pneumonitis", "IPF", "Ever-smokers, no lung disease",
            "Frontal cortex", "AD (pooled)", "Huntington's disease", "Parkinson's disease", "Schizophrenia", "Bipolar disorder", "Major depressive disorder"),
  ds_label = c("", "GSE47460", "GSE150910", "GSE150910", "GSE150910", "", "6 datasets", "GSE64810", "GSE68719", "GSE53987", "GSE53987", "GSE53987"))
meta <- fread(file.path(BR, "results/AD_reference_meta.tsv"))
r <- rbind(r, meta[role == "AD_pooled", .(id, role = "target_pooled", n_case, n_control, effect_SD, ci_low, ci_high, p, fdr_family = NA_real_, call = "descriptive")], fill = TRUE)
d <- merge(rows, r, by = "id", all.x = TRUE, sort = FALSE)[match(rows$label, label)]
stopifnot(nrow(d) == 12, sum(!is.na(d$effect_SD)) == 10)
d <- d[is.na(role) | !grepl("^target", role)]  # drop COPD / AD reference rows
stopifnot(nrow(d) == 10)
fwrite(d[!is.na(effect_SD), .(label, ds_label, role, n_case, n_control, effect_SD, ci_low, ci_high, p, fdr_family, call)],
       file.path(out, "Figure4_AxisA_disease_scope_forest_comparators_only_source.csv"))

width_mm <- 180; height_mm <- 88
font <- "Arial"; fs <- 10; ink <- "#262626"; muted <- "#555B61"; coral <- "#EA5B5B"; grid_col <- "#BFC3C7"
records <- list()
txt <- function(label, x, y, just = "left", bold = FALSE, colour = ink) {
  g <- textGrob(label, x = unit(x, "mm"), y = unit(y, "mm"), just = c(just, "centre"),
                gp = gpar(fontfamily = font, fontsize = fs, fontface = if (bold) "bold" else "plain", col = colour))
  grid.draw(g)
  w <- convertWidth(grobWidth(g), "mm", TRUE); x0 <- switch(just, left = x, centre = x - w / 2, right = x - w)
  records[[length(records) + 1L]] <<- data.frame(label = label, x0_mm = x0, x1_mm = x0 + w, y_mm = y)
}
seg <- function(x0, y0, x1, y1, col = ink, lwd = 1.4, lty = 1)
  grid.lines(unit(c(x0, x1), "mm"), unit(c(y0, y1), "mm"), gp = gpar(col = col, lwd = lwd, lty = lty))
fmt_p <- function(p) if (p < 0.001) "< 0.001" else if (p < 0.01) sprintf("%.4f", p) else if (p < 0.1) sprintf("%.3f", p) else sprintf("%.2f", p)

draw <- function() {
  grid.newpage()
  left <- 86; right <- 128; lim <- c(-1.5, 1.5)
  X <- function(v) left + (v - lim[1]) / diff(lim) * (right - left)
  top <- 76; step <- 6.4; ys <- top - (seq_len(nrow(d)) - 1) * step
  txt("Context", 3, 83, bold = TRUE); txt("Dataset", 64, 83, bold = TRUE)
  txt("Effect [95% CI]", 132, 83, bold = TRUE); txt("FDR", 178, 83, just = "right", bold = TRUE)
  seg(X(0), ys[length(ys)] - 3.5, X(0), top + 3, col = grid_col, lwd = 1, lty = 2)
  for (i in seq_len(nrow(d))) {
    row <- d[i]; y <- ys[i]
    if (is.na(row$effect_SD)) { txt(row$label, 3, y, bold = TRUE); next }
    is_target <- grepl("^target", row$role)
    txt(row$label, 6, y); txt(row$ds_label, 64, y, colour = muted)
    col <- if (is_target) coral else ink
    seg(X(max(row$ci_low, lim[1])), y, X(min(row$ci_high, lim[2])), y, col = col)
    filled <- is_target || row$call %in% c("detected increase", "detected decrease")
    if (row$role == "target_pooled") grid.polygon(unit(X(row$effect_SD) + c(-1.9, 0, 1.9, 0), "mm"), unit(y + c(0, 1.5, 0, -1.5), "mm"), gp = gpar(fill = coral, col = coral))
    else if (is_target) grid.rect(unit(X(row$effect_SD), "mm"), unit(y, "mm"), unit(2.4, "mm"), unit(2.4, "mm"), gp = gpar(fill = coral, col = coral))
    else grid.circle(unit(X(row$effect_SD), "mm"), unit(y, "mm"), unit(1.3, "mm"),
                     gp = gpar(fill = if (filled) ink else "white", col = ink, lwd = 1.4))
    txt(gsub("-", "\u2212", sprintf("%.2f [%.2f, %.2f]", row$effect_SD, row$ci_low, row$ci_high), fixed = TRUE), 132, y)
    txt(if (is_target) paste0("P ", fmt_p(row$p)) else fmt_p(row$fdr_family), 178, y, just = "right")
  }
  ax <- ys[length(ys)] - 5.5
  seg(left, ax, right, ax)
  for (v in seq(-1.5, 1.5, 0.5)) { seg(X(v), ax, X(v), ax - 1); if (v %% 1 == 0 || v == 0) txt(sub("-", "\u2212", format(v), fixed = TRUE), X(v), ax - 3.5, just = "centre") }
  txt("Axis A score difference (SD)", (left + right) / 2, ax - 8, just = "centre")
}
stem <- "Figure4_AxisA_disease_scope_forest_comparators_only"
quartz(type = "pdf", file = file.path(out, paste0(stem, ".pdf")), width = width_mm / 25.4, height = height_mm / 25.4, family = font)
draw(); dev.off()
records <- list()
quartz(type = "png", file = file.path(out, paste0(stem, "_preview.png")), width = width_mm / 25.4, height = height_mm / 25.4, dpi = 300, family = font, bg = "white")
draw(); dev.off()
rec <- do.call(rbind, records)
ov <- do.call(rbind, lapply(seq_len(nrow(rec)), function(i) {
  j <- which(abs(rec$y_mm - rec$y_mm[i]) < 3 & seq_len(nrow(rec)) > i & rec$x0_mm < rec$x1_mm[i] & rec$x1_mm > rec$x0_mm[i])
  if (length(j)) data.frame(a = rec$label[i], b = rec$label[j]) }))
cat("text collisions:", if (is.null(ov)) 0 else nrow(ov), "\n"); if (!is.null(ov)) print(ov)
cat("max text right edge (mm):", round(max(rec$x1_mm), 1), "of", width_mm, "\n")
