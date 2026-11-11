# Fig 4. Cellular localization of Axis A.
# Sources (frozen, not recomputed):
#   cell_localization/SYN01_Entrez_downstream_20261001/source_data/{Lung,Brain}_UMAP_cells.tsv.gz, *_plotted_points.tsv, *_palette.tsv,
#     brain_subtype_display_map_20260929.tsv (a-f)
#   REF/{lung,brain}_null.tsv, REF/localization.tsv (g, h)
#   10-29: lung panels e and g read 03_candidate_analysis/Lung_sc_donor_merge_20261029 (GSE302339 cood8 = copd8, 55 donors)
# Designs follow the 10-05 figure.af old Fig 5 (cell_localization/SYN01_Entrez_downstream_20261001 *_style_restored.R): category palette,
# grey-to-red UCell ramp, per-subtype coloured rainclouds, grey null histograms with a red Axis A line.
# Layout: lung UMAPs a, b and brain UMAPs c, d (each about 55 mm wide) above the null histograms g, h; subtype rainclouds e, f on the right.
source("06_resubmission_20261008/code/main_figures/common.R")
SC <- "cell_localization/SYN01_Entrez_downstream_20261001/source_data"
LM <- "03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables"  # lung donor-level values with cood8 -> copd8

# ---- a-d: UMAPs ---------------------------------------------------------------------------------------------
rd <- function(f) fread(cmd = paste("gzip -cd", shQuote(P(file.path(SC, f)))))
umap_axes <- function(d) {
  xr <- range(d$UMAP_1); yr <- range(d$UMAP_2); L <- 0.16 * max(diff(xr), diff(yr)); x0 <- xr[1]; y0 <- yr[1]
  list(annotate("segment", x = x0, xend = x0 + L, y = y0, yend = y0, linewidth = 0.45, colour = INK,
                arrow = arrow(length = unit(1.1, "mm"), type = "closed")),
       annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + L, linewidth = 0.45, colour = INK,
                arrow = arrow(length = unit(1.1, "mm"), type = "closed")),
       annotate("text", x = x0 + L / 2, y = y0, label = "UMAP 1", vjust = 1.5, size = FS * PT, fontface = "bold", colour = INK),
       annotate("text", x = x0, y = y0 + L / 2, label = "UMAP 2", vjust = -0.6, angle = 90, size = FS * PT, fontface = "bold", colour = INK))
}
umap_theme <- function() theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = FS, face = "bold", colour = INK, hjust = 0.5, margin = margin(0, 0, 1, 0, "mm")),
        legend.text = element_text(size = FS, face = "bold", colour = INK), legend.title = element_text(size = FS, face = "bold", colour = INK),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
umap_type <- function(d, pal, title, lab, nudge) {
  cen <- d[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = group]
  cen[, label := lab[group]][is.na(label), label := group]
  cen[, nx := vapply(group, function(g) nudge[[g]][1] %||% 0, 0)][, ny := vapply(group, function(g) nudge[[g]][2] %||% 0, 0)]
  ggplot(d[sample(.N)], aes(UMAP_1, UMAP_2)) +
    ggrastr::rasterise(geom_point(aes(colour = group), size = 0.12, stroke = 0, shape = 16), dpi = 600, dev = "ragg") +
    geom_segment(data = cen[nx != 0 | ny != 0], aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1 + nx * 0.82, yend = UMAP_2 + ny * 0.82),
                 linewidth = 0.3, colour = INK) +
    geom_label(data = cen, aes(x = UMAP_1 + nx, y = UMAP_2 + ny, label = label), size = FS * PT, family = "Arial", fontface = "bold",
               colour = INK, fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.4, "mm"), lineheight = 0.85) +
    scale_colour_manual(values = setNames(pal$color, pal$group), guide = "none") +
    umap_axes(d) + coord_equal(clip = "off") + labs(title = title) + umap_theme()
}
ramp <- c("#E2E5E9", "#F7B0A1", "#F04438", "#A3231C")
umap_score <- function(d, title) {
  hi <- quantile(d$AxisA_Entrez_UCell, 0.995)
  ggplot(d[order(AxisA_Entrez_UCell)], aes(UMAP_1, UMAP_2)) +
    ggrastr::rasterise(geom_point(aes(colour = AxisA_Entrez_UCell), size = 0.12, stroke = 0, shape = 16), dpi = 600, dev = "ragg") +
    scale_colour_gradientn(colours = ramp, limits = c(0, hi), oob = scales::squish, name = "UCell",
                           breaks = function(l) pretty(c(0, l[2]), 3)) +
    guides(colour = guide_colourbar(barheight = unit(20, "mm"), barwidth = unit(2.2, "mm"), ticks.colour = INK)) +
    umap_axes(d) + coord_equal(clip = "off") + labs(title = title) + umap_theme() +
    theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.5, "mm"))
}
set.seed(1)
lu <- rd("Lung_UMAP_cells.tsv.gz"); br <- rd("Brain_UMAP_cells.tsv.gz")
stopifnot(nrow(lu) == 211464, nrow(br) == 1357145)
lpal <- fread(P(file.path(SC, "Lung_palette.tsv"))); bpal <- fread(P(file.path(SC, "Brain_palette.tsv")))
lsp <- diff(range(lu$UMAP_1)) / 100; bsp <- diff(range(br$UMAP_1)) / 100
pa <- umap_type(lu, lpal, "Lung: cell types", c("Other/Progenitor" = "Other"),
                list(Stromal = c(-4, 9) * lsp, Endothelial = c(-6, 5) * lsp, Epithelial = c(-6, -10) * lsp,
                     `B/Plasma` = c(-2, 9) * lsp, `Other/Progenitor` = c(-8, 10) * lsp))
pb <- umap_score(lu, "Lung: Axis A UCell")
pc <- umap_type(br, bpal, "Brain: cell types",
                c("Excitatory neurons" = "Excitatory\nneurons", "Inhibitory neurons" = "Inhibitory\nneurons",
                  "VLMC & Perivascular" = "VLMC &\nperivascular"),
                list(Immune = c(-6, 6) * bsp, Astrocyte = c(-4, -10) * bsp, Endothelial = c(-12, 1) * bsp, OPC = c(-10, -6) * bsp,
                     `VLMC & Perivascular` = c(12, -7) * bsp, Oligodendrocyte = c(9, 12) * bsp,
                     `Excitatory neurons` = c(14, 5) * bsp))
pd <- umap_score(br, "Brain: Axis A UCell")

# ---- e, f: donor-mean UCell by subtype (old coloured raincloud) -------------------------------------------------------
lung_palette <- c("Alveolar macrophages" = PAL8[1], "Monocyte-derived macrophages" = PAL8[2], "Perivascular interstitial macrophages" = PAL8[3],
                  "Classical monocytes" = PAL8[4], "Conventional DCs" = PAL8[5], "Non-classical monocytes" = PAL8[6],
                  "Plasmacytoid DCs" = PAL8[7], "Mast cells" = PAL8[8])
brain_palette <- c("Micro-PVM_2" = PAL8[1], "Micro-PVM_3-SEAAD" = PAL8[2], "Micro-PVM_2_3-SEAAD" = PAL8[3], "Micro-PVM_1" = PAL8[4],
                   "Micro-PVM_2_1-SEAAD" = PAL8[6], "Lymphocyte" = PAL8[7])
rain <- function(pts, pal, title, subtitle, xbreaks) {
  pts <- copy(pts); lv <- pts[, .(o = first(category_order)), by = lab][order(o)]$lab
  # a three-line top label gets extra room above the second row
  pts[, yy := length(lv) + 1 - match(lab, lv) + fifelse(lab == lv[1] & grepl("\n.*\n", lab), 0.3, 0)][, col := pal[original_annotation]]
  stopifnot(!anyNA(pts$col))
  v <- pts[, { k <- density(mean_ucell); keep <- k$x >= min(mean_ucell) & k$x <= max(mean_ucell)
               .(x = k$x[keep], h = k$y[keep] / max(k$y)) }, by = .(lab, yy, col)]
  bx <- pts[, { s <- boxplot.stats(mean_ucell, coef = 1.5)$stats; .(lo = s[1], q1 = s[2], med = s[3], q3 = s[4], hi = s[5]) }, by = .(lab, yy, col)]
  set.seed(2); pts[, jy := yy - 0.3 + runif(.N, -0.06, 0.06)]
  ggplot() +
    geom_polygon(data = v[, .(x = c(x, rev(x)), y = c(yy + 0.04 + 0.4 * h, rep(yy + 0.04, .N))), by = .(lab, yy, col)],
                 aes(x, y, group = lab, fill = col, colour = col), alpha = 0.45, linewidth = 0.35) +
    geom_segment(data = bx, aes(x = lo, xend = hi, y = yy - 0.1, yend = yy - 0.1, colour = col), linewidth = 0.35) +
    geom_rect(data = bx, aes(xmin = q1, xmax = q3, ymin = yy - 0.16, ymax = yy - 0.04, colour = col), fill = "white", linewidth = 0.35) +
    geom_segment(data = bx, aes(x = med, xend = med, y = yy - 0.16, yend = yy - 0.04), linewidth = 0.5, colour = INK) +
    geom_point(data = pts, aes(mean_ucell, jy, colour = col), size = 0.45, alpha = 0.75, stroke = 0) +
    scale_colour_identity() + scale_fill_identity() +
    scale_y_continuous(breaks = rev(unique(pts[order(category_order)]$yy)), labels = rev(lv), expand = expansion(add = c(0.5, 0.5))) +
    scale_x_continuous(breaks = xbreaks) +
    labs(x = "Donor-mean Axis A UCell", y = NULL, title = title, subtitle = subtitle) +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(lineheight = 0.85))
}
lp <- fread(P(file.path(LM, "Lung_plotted_points.tsv"))); stopifnot(uniqueN(lp$donor) == 55)
lp[, lab := c("Perivascular interstitial macrophages" = "Perivascular\ninterstitial\nmacrophages", "Alveolar macrophages" = "Alveolar\nmacrophages",
              "Monocyte-derived macrophages" = "Monocyte-derived\nmacrophages", "Classical monocytes" = "Classical\nmonocytes",
              "Non-classical monocytes" = "Non-classical\nmonocytes", "Conventional DCs" = "Conventional\nDCs",
              "Plasmacytoid DCs" = "Plasmacytoid\nDCs")[annotation]]
lp[is.na(lab), lab := annotation]
bpp <- fread(P(file.path(SC, "Brain_plotted_points.tsv")))
dm <- fread(P(file.path(SC, "brain_subtype_display_map_20260929.tsv")))
bpp[, lab := dm$display_label[match(annotation, dm$original_label)]][lab == "Perivascular macrophage 1", lab := "Perivascular\nmacrophage 1"]
stopifnot(!anyNA(bpp$lab))
pe <- rain(lp, lung_palette, sprintf("Lung: myeloid subtypes (%d donors)", uniqueN(lp$donor)), NULL, seq(0, 0.3, 0.1))
pf <- rain(bpp, brain_palette, sprintf("Brain: immune supertypes (%d donors)", uniqueN(bpp$donor)), NULL, seq(0, 0.2, 0.05))

# ---- g, h: matched random gene sets (old grey histogram, red Axis A line) -------------------------------------------
loc <- fread(P(file.path(REF, "localization.tsv")))
loc_lung <- fread(P(file.path(LM, "lung_broad_primary_contrast.tsv"))); stopifnot(loc_lung$n_donors == 55)
nullp <- function(f, obs, title, xlab) {
  d <- fread(P(f)); stopifnot(nrow(d) == 5000)
  ggplot(d, aes(effect)) +
    geom_histogram(bins = 30, fill = COL$hist, colour = COL$histline, linewidth = 0.2) +
    annotate("segment", x = obs, xend = obs, y = 0, yend = Inf, colour = COL$axisA, linewidth = 0.9) +
    annotate("text", x = obs, y = Inf, label = "Axis A", hjust = 1.1, vjust = 1.2, colour = COL$axisA, size = FS * PT, fontface = "bold") +
    scale_x_continuous(limits = c(0, NA), breaks = function(l) pretty(c(0, l[2]), 3)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlab, y = "Number of random sets", title = title, subtitle = "5,000 matched sets")
}
pg <- nullp(file.path(LM, "lung_matched_null_draws.tsv"), loc_lung$effect, "Lung: myeloid", "UCell difference\n(target - others)")
ph <- nullp(file.path(REF, "brain_null.tsv"), loc[cohort == "SEA-AD"]$effect, "Brain: immune", "UCell difference\n(target - others)")

compose_figure(list(
  blk(pa, 0, 0, 60, 72, "a"), blk(pb, 60, 0, 60, 72, "b"),
  blk(pc, 0, 74, 60, 72, "c"), blk(pd, 60, 74, 60, 72, "d"),
  blk(pg, 0, 148, 60, 58, "g"), blk(ph, 60, 148, 60, 58, "h"),
  blk(pe / pf + plot_layout(heights = c(8, 6)), 120, 0, 70, 206, c("e", "f"), tags = list(c("e", 0, 0), c("f", 0, 113.8)))
), "Fig4", 190, 206, row_groups = list(c("g", "h")), column_groups = list(c("e", "f")),
exemptions = list(list(panels = list("a", "b", "c", "d"), checks = list("row", "column", "panel-width"),
                       reason = "UMAPs use equal coordinates; plot areas follow the embedding's aspect ratio")))
