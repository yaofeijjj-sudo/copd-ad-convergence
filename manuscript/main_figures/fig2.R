# Fig 2. Disease-associated co-expression modules in lung and brain and their overlap.
# Sources (frozen, not recomputed): REF/lung_module_traits.csv, REF/brain_module_traits.csv, REF/module_overlap.csv,
# REF/adjusted_effects.csv (05_github_release reference tables, 137-gene Entrez version).
# Layout and panel designs follow the 10-05 figure.af: a, b = old 1c, 1d (WGCNA labeledHeatmap, blueWhiteRed);
# c = old 2a (cross WGCNA/bubble_final_20260917); d = old 2b, 2c (venn_candidate_20260917); e = old 2d (forest_revision_20260917).
source("06_resubmission_20261008/code/main_figures/common.R")

# ---- a, b: module-trait heatmaps -------------------------------------------------------------
lt <- fread(P(file.path(REF, "lung_module_traits.csv")))
bt <- fread(P(file.path(REF, "brain_module_traits.csv")))
lung_mod <- c("red", "blue", "yellow", "black", "brown", "green", "pink", "turquoise", "grey")
brain_mod <- c("turquoise", "midnightblue", "lightcyan", "tan", "cyan", "red", "greenyellow", "green", "salmon",
               "yellow", "purple", "black", "blue", "pink", "brown", "magenta", "grey")
lung_tr <- c(copd_status = "COPD status", GOLD = "GOLD", FEV1_pp = "FEV1 %pred", FEV1_FVC = "FEV1/FVC",
             laa950 = "LAA-950", perc15 = "Perc15", pi10 = "Pi10")
brain_tr <- c(Disease_State = "AD status", Age = "Age", Sex = "Sex")
pshort <- function(p) { r <- signif(p, 1); ifelse(r >= 0.01, trimws(formatC(r, format = "fg", digits = 1)), sub("e-0*", "e-", sprintf("%.0e", p))) }
cell_lab <- function(r, p) sprintf("%s\n(%s)", sprintf("%.2f", ifelse(round(r, 2) == 0, 0, r)), pshort(p))
bwr <- WGCNA::blueWhiteRed(50)

heat <- function(d, mods, trs, title, subtitle) {
  d <- copy(d)[module %in% mods & trait %in% names(trs)]
  d[, x := match(trait, names(trs))][, y := length(mods) + 1 - match(module, mods)]
  d[, lab := cell_lab(correlation, p_value)]
  strip <- data.table(y = seq_along(mods), module = rev(mods))
  nx <- length(trs)
  ggplot(d) +
    geom_tile(aes(x, y, fill = correlation)) +
    annotate("rect", xmin = 0.5, xmax = nx + 0.5, ymin = 0.5, ymax = length(mods) + 0.5, fill = NA, colour = INK, linewidth = 0.5) +
    geom_text(aes(x, y, label = lab), size = FS * PT, lineheight = 0.88, family = "Arial", fontface = "bold", colour = INK) +
    geom_tile(data = strip, aes(x = 0.24, y = y), fill = strip$module, width = 0.4, height = 1) +
    scale_fill_gradientn(colours = bwr, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = NULL,
                         labels = c("-1", "-0.5", "0", "0.5", "1")) +
    scale_x_continuous(breaks = seq_along(trs), labels = trs, expand = c(0, 0)) +
    scale_y_continuous(breaks = strip$y, labels = paste0("ME", strip$module), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.04, nx + 0.5), clip = "off") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    guides(fill = guide_colourbar(barheight = unit(30, "mm"), barwidth = unit(2.6, "mm"), frame.colour = INK,
                                  frame.linewidth = 0.3, ticks.colour = INK)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(margin = margin(0, 0.6, 0, 0, "mm")),
          legend.position = "right", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "mm"))
}
pa <- heat(lt, lung_mod, lung_tr, "COPD module–trait relationships", "GSE76925, n = 147")
pb <- heat(bt, brain_mod, brain_tr, "AD module–trait relationships", "GSE33000, n = 459")

# ---- c: 128 module-pair overlaps (old 2a design) -----------------------------------------------------------
ov <- fread(P(file.path(REF, "module_overlap.csv")))
stopifnot(nrow(ov) == 128, all(ov$BH_family_size == 128), sum(ov$fdr < 0.05) == 18)
lm8 <- sort(unique(ov$lung_module)); bm16 <- sort(unique(ov$brain_module))
ov[, cls := fifelse(fdr >= 0.05, "BH FDR \u2265 0.05",
              fifelse(direction == "both_positive", "Both increased", fifelse(direction == "both_negative", "Both decreased", "Opposite")))]
ov[, cls := factor(cls, c("Both increased", "Both decreased", "Opposite", "BH FDR \u2265 0.05"))]
ov[, x := match(brain_module, bm16)][, y := length(lm8) + 1 - match(lung_module, lm8)]
xs <- data.table(x = seq_along(bm16), m = bm16); ys <- data.table(y = seq_along(lm8), m = rev(lm8))
aA <- ov[lung_module == "red" & brain_module == "blue"]; aB <- ov[lung_module == "turquoise" & brain_module == "turquoise"]
pc <- ggplot(ov[overlap > 0]) +
  geom_vline(xintercept = xs$x, colour = COL$grid, linewidth = 0.4) + geom_hline(yintercept = ys$y, colour = COL$grid, linewidth = 0.4) +
  geom_point(aes(x, y, size = overlap, colour = cls)) +
  geom_tile(data = xs, aes(x = x, y = 0.3), fill = xs$m, width = 0.76, height = 0.22) +
  geom_tile(data = ys, aes(x = 0.3, y = y), fill = ys$m, width = 0.22, height = 0.72) +
  annotate("text", x = aA$x + 0.5, y = aA$y + 0.45, label = "Axis A", hjust = 0, size = FS * PT, fontface = "bold", colour = COL$axisA) +
  annotate("text", x = aB$x - 0.6, y = aB$y + 0.5, label = "Axis B", hjust = 1, size = FS * PT, fontface = "bold", colour = COL$axisB) +
  scale_colour_manual(values = c("Both increased" = COL$up, "Both decreased" = COL$down, "Opposite" = COL$opposite,
                                 "BH FDR \u2265 0.05" = COL$light), name = NULL, drop = FALSE) +
  scale_size_area(max_size = 7.2, breaks = c(50, 200, 400), name = "Shared genes") +
  scale_x_continuous(breaks = xs$x, labels = xs$m, expand = expansion(add = c(0.75, 0.5))) +
  scale_y_continuous(breaks = ys$y, labels = ys$m, expand = expansion(add = c(0.75, 0.5))) +
  coord_cartesian(clip = "off") +
  labs(x = "AD brain module", y = "COPD lung module", title = "Cross-disease module overlap") +
  guides(colour = guide_legend(override.aes = list(size = 2.6), order = 2, nrow = 2),
         size = guide_legend(order = 1, override.aes = list(colour = "#DDDDDD"))) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "top", legend.box = "vertical", legend.box.just = "left", legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.5, "mm"), legend.key.size = unit(3.2, "mm"), legend.box.spacing = unit(1, "mm"))

# ---- d: Axis A and Axis B definitions (old 2b, 2c: pale ellipses, solid intersection) --------------------------
ell <- function(cx, cy, a = 2.2, b = 1.3, n = 401) { t <- seq(0, 2 * pi, length.out = n); data.table(x = cx + a * cos(t), y = cy + b * sin(t)) }
lens <- function(c1, c2, cy, a = 2.2, b = 1.3) {
  th <- acos(((c2 - c1) / 2) / a); t1 <- seq(-th, th, length.out = 201); t2 <- seq(pi - th, pi + th, length.out = 201)
  data.table(x = c(c1 + a * cos(t1), c2 + a * cos(t2)), y = c(cy + b * sin(t1), cy + b * sin(t2)))
}
venn <- function(row, x0, col, pale, axis, lmod, bmod) {
  c1 <- x0 - 1.4; c2 <- x0 + 1.4
  list(geom_polygon(data = ell(c1, 0), aes(x, y), fill = pale, colour = NA),
       geom_polygon(data = ell(c2, 0), aes(x, y), fill = pale, colour = NA),
       geom_polygon(data = lens(c1, c2, 0), aes(x, y), fill = col, colour = NA),
       geom_path(data = ell(c1, 0), aes(x, y), colour = col, linewidth = 0.5),
       geom_path(data = ell(c2, 0), aes(x, y), colour = col, linewidth = 0.5),
       annotate("text", x = c(c1 - 0.95, x0, c2 + 0.95), y = 0, label = format(c(row$fisher_b, row$fisher_a, row$fisher_c), big.mark = ","),
                colour = c(INK, "white", INK), size = FS * PT, fontface = "bold"),
       annotate("text", x = x0, y = 3.55, label = axis, colour = col, size = FS * PT, fontface = "bold"),
       annotate("text", x = c(c1 - 0.9, c2 + 0.9), y = 2.05, label = c(paste0("COPD lung\n", lmod), paste0("AD brain\n", bmod)),
                colour = INK, size = FS * PT, fontface = "bold", lineheight = 0.88),
       annotate("text", x = x0, y = c(-1.75, -2.45), label = c(sprintf("OR = %.2f", row$odds_ratio), sprintf("BH FDR = %s", fmt_p_sup(row$fdr, 3))),
                colour = INK, size = FS * PT, fontface = "bold"))
}
pd <- ggplot() + venn(aA, 0, COL$axisA, "#FDE4E2", "Axis A", "red", "blue") +
  venn(aB, 8.6, COL$axisB, "#E0EFFD", "Axis B", "turquoise", "turquoise") +
  coord_equal(xlim = c(-3.7, 12.3), ylim = c(-2.8, 3.9), clip = "off") + theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))

# ---- e: covariate-adjusted source-module effects (old 2d table forest) --------------------------------------------
ae <- fread(P(file.path(REF, "adjusted_effects.csv")))
stopifnot(nrow(ae) == 4)
ae[, y := c(5, 4, 2, 1)][, col := fifelse(axis == "Axis A", COL$forestA, COL$forestB)]
ae[, `:=`(est = effect, lo = lower, hi = upper, shape = fifelse(organ == "Lung", "circle", "square"), size = 2.1)]
ae[, mod := sprintf("%s %s", organ, module)][, ci := sprintf("%.2f [%.2f, %.2f]", effect, lower, upper)]
ae[, ftxt := fmt_p_sup(fdr, 3)]
ehdr <- data.table(y = c(6, 3), lab = c("Axis A", "Axis B"), colr = c(COL$forestA, COL$forestB))
ytop <- 7.1; yl <- c(0.5, 7.4)
el <- ftext(ae, list(list(var = "mod", head = "Module", x = 0.08)), ehdr, ytop, yl, title = "Adjusted module\u2013disease associations")
ef <- fplot(ae, c(-1.5, 1.5), c(-1.5, 0, 1.5), "Adjusted case\u2013control\ndifference (ME SD)", ytop, yl, head = "Adjusted effect")
er <- ftext(ae, list(list(var = "ci", head = "SD [95% CI]", x = 0.04), list(var = "ftxt", head = "BH FDR", x = 0.7)),
            data.table(), ytop, yl)
pe <- el + ef + er + plot_layout(widths = c(25, 30, 45))

fwrite(ov[, .(lung_module, brain_module, overlap, odds_ratio, fdr, direction, cls)], P(file.path(SRC, "Fig2c_module_overlap.tsv")), sep = "\t")

compose_figure(list(
  blk(pa, 0, 0, 113, 88, "a"),
  blk(pb, 114, 0, 76, 166, "b"),
  blk(pc, 0, 88, 113, 78, "c"),
  blk(pd, 0, 168, 84, 52, "d"),
  blk(pe, 86, 168, 104, 52, c("e_left", "e", "e_right"), tags = list(c("e", 0, 0)))
), "Fig2", 190, 220, row_groups = list(c("e_left", "e", "e_right")),
exemptions = list(list(panels = list("e_left", "e", "e_right"), checks = list("panel-width"), reason = "table forest: text columns and the forest differ in width by design"),
                  list(panels = list("a", "b", "c"), checks = list("row", "column", "panel-width"),
                       reason = "heatmaps and the bubble matrix size their cells to the number of rows and columns")))
