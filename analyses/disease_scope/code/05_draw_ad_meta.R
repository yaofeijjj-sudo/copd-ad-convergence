#!/usr/bin/env Rscript
# Supplementary forest: Axis A in AD vs control, six independent frontal-cortex datasets + DL random-effects pool.
# Same style as 03_draw_forest.R. Run from project root:
#   _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/05_draw_ad_meta.R
suppressPackageStartupMessages({library(grid); library(data.table)})
BR <- "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003"
out <- file.path(BR, "figures")
m <- fread(file.path(BR, "results/AD_reference_meta.tsv"))
region <- c(GSE84422 = "DLPFC", GSE118553 = "Frontal cortex", GSE313307 = "Prefrontal cortex", GSE140841 = "BA9",
            GSE48350 = "Superior frontal gyrus", GSE36980 = "Frontal cortex")
d <- m[role == "AD_dataset"][, region := region[id]]
pool <- m[role == "AD_pooled"]
stopifnot(nrow(d) == 6, nrow(pool) == 1)
fwrite(rbind(d, pool, fill = TRUE)[, .(id, region, n_case, n_control, axisA_covered, covariates, effect_SD, se, ci_low, ci_high, p, I2, tau2, Q_p)],
       file.path(out, "SuppFigure_AD_reference_meta_source.csv"))

width_mm <- 180; height_mm <- 78
font <- "Arial"; fs <- 10; ink <- "#262626"; muted <- "#555B61"; coral <- "#EA5B5B"; grid_col <- "#BFC3C7"
txt <- function(label, x, y, just = "left", bold = FALSE, colour = ink)
  grid.text(label, unit(x, "mm"), unit(y, "mm"), just = c(just, "centre"),
            gp = gpar(fontfamily = font, fontsize = fs, fontface = if (bold) "bold" else "plain", col = colour))
seg <- function(x0, y0, x1, y1, col = ink, lwd = 1.4, lty = 1)
  grid.lines(unit(c(x0, x1), "mm"), unit(c(y0, y1), "mm"), gp = gpar(col = col, lwd = lwd, lty = lty))
mn <- function(s) gsub("-", "−", s, fixed = TRUE)
fp <- function(p) if (p < 0.001) "< 0.001" else if (p < 0.01) sprintf("%.4f", p) else if (p < 0.1) sprintf("%.3f", p) else sprintf("%.2f", p)
draw <- function() {
  grid.newpage()
  left <- 92; right <- 130; lim <- c(-1, 2); X <- function(v) left + (v - lim[1]) / diff(lim) * (right - left)
  ys <- 64 - (0:5) * 6.4; yp <- ys[6] - 8
  txt("Dataset", 3, 72, bold = TRUE); txt("Region", 28, 72, bold = TRUE); txt("AD/control", 66, 72, bold = TRUE)
  txt("Effect [95% CI]", 134, 72, bold = TRUE); txt("P", 178, 72, just = "right", bold = TRUE)
  seg(X(0), yp - 3.5, X(0), ys[1] + 3, col = grid_col, lwd = 1, lty = 2)
  for (i in seq_len(nrow(d))) { r <- d[i]; y <- ys[i]
    txt(r$id, 3, y); txt(r$region, 28, y, colour = muted); txt(sprintf("%d/%d", r$n_case, r$n_control), 66, y)
    seg(X(max(r$ci_low, lim[1])), y, X(min(r$ci_high, lim[2])), y)
    grid.circle(unit(X(r$effect_SD), "mm"), unit(y, "mm"), unit(1.3, "mm"), gp = gpar(fill = ink, col = ink))
    txt(mn(sprintf("%.2f [%.2f, %.2f]", r$effect_SD, r$ci_low, r$ci_high)), 134, y); txt(fp(r$p), 178, y, just = "right") }
  seg(3, yp + 4, 178, yp + 4, col = grid_col, lwd = 0.8)
  txt("Pooled (random effects)", 3, yp, bold = TRUE); txt(sprintf("%d/%d", pool$n_case, pool$n_control), 66, yp, bold = TRUE)
  grid.polygon(unit(c(X(pool$ci_low), X(pool$effect_SD), X(pool$ci_high), X(pool$effect_SD)), "mm"), unit(yp + c(0, 1.6, 0, -1.6), "mm"), gp = gpar(fill = coral, col = coral))
  txt(mn(sprintf("%.2f [%.2f, %.2f]", pool$effect_SD, pool$ci_low, pool$ci_high)), 134, yp, bold = TRUE); txt(fp(pool$p), 178, yp, just = "right", bold = TRUE)
  txt(sprintf("I² = %.0f%%, heterogeneity P = %.2f", pool$I2, pool$Q_p), 3, yp - 6, colour = muted)
  ax <- yp - 6.5; seg(left, ax, right, ax)
  for (v in seq(-1, 2, 0.5)) { seg(X(v), ax, X(v), ax - 1); if (v %% 1 == 0) txt(mn(format(v)), X(v), ax - 3.5, just = "centre") }
  txt("Axis A score difference (SD)", (left + right) / 2, ax - 8, just = "centre")
}
stem <- "SuppFigure_AD_reference_meta"
quartz(type = "pdf", file = file.path(out, paste0(stem, ".pdf")), width = width_mm / 25.4, height = height_mm / 25.4, family = font); draw(); dev.off()
quartz(type = "png", file = file.path(out, paste0(stem, "_preview.png")), width = width_mm / 25.4, height = height_mm / 25.4, dpi = 300, family = font, bg = "white"); draw(); invisible(dev.off())
