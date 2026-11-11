# Fig 8. Candidate upstream transcription factors of Axis A (title since 11-09; panels d, e and f since 11-09, formerly d-f and g).
# Sources (frozen, not recomputed):
#   03_candidate_analysis/CellOracle_ranking_scale_20261019/results/all357_TF_ranking.tsv (b, c)
#   03_candidate_analysis/CellOracle_ranking_scale_20261019/results/Fig6def_source_signed_logratio.tsv (d, e = STAT3, CEBPB; delta_X on the ln scale)
#   tf_activity_all.tsv.gz in CO (f; primary models, full CollecTRI; lung rows rerun 10-29, brain rows from SingleCell_pseudobulk_supplements_20261019)
# Designs follow the 10-05 figure.af old Fig 6: a = virtual_knockout/CellOracle_method_schematic_20260923 (redrawn at final size),
# b = Figure6_rank_candidate_20260925 (two stacked organ panels, grid lines, blue labelled TFs), c = lollipop, d, e = zero-line scatter.
source("06_resubmission_20261008/code/main_figures/common.R")
CO <- "03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables"  # 10-29: lung donor means with cood8 -> copd8 (55 donors); brain unchanged
blue <- "#2878B5"; dblue <- "#12486C"; coral <- "#CE513F"; dred <- "#B62F28"; greyp <- "#AEB8C2"; teal <- "#168A87"

# ---- a: workflow schematic (old 6a, laid out on a 52 x 120 mm canvas) ----------------------------------------------
W6 <- 52; H6 <- 120
net <- function(cx, cy, ko = FALSE) {
  p <- data.table(x = cx + c(-13, 3, 12, 2, -13), y = cy + c(0.5, 7, -1, -8, -8))
  e <- data.table(a = c(1, 1, 1, 2, 5), b = c(2, 3, 5, 4, 4))
  seg <- e[, { A <- p[a]; B <- p[b]; v <- c(B$x - A$x, B$y - A$y); v <- v / sqrt(sum(v^2)); r0 <- if (a == 1) 3.6 else 3.0
               .(x = A$x + v[1] * r0, y = A$y + v[2] * r0, xend = B$x - v[1] * 3.3, yend = B$y - v[2] * 3.3,
                 col = if (!ko && a == 1) coral else "#8BA1B5") }, by = .(a, b)]
  fills <- if (ko) c("white", "#98BBD8", "#E2ABA6", "#5D93C0", "#F2D4CC") else c(coral, "#DCEAF5", "#DCEAF5", "#DCEAF5", "#DCEAF5")
  p[, `:=`(fill = fills, r = c(3.6, 3, 3, 3, 3), stroke = c(if (ko) coral else NA, NA, NA, NA, NA))]
  list(geom_segment(data = seg, aes(x, y, xend = xend, yend = yend, colour = col), linewidth = 0.5,
                    arrow = arrow(length = unit(1.2, "mm"), type = "closed"), inherit.aes = FALSE),
       ggforce_circle(p),
       annotate("text", x = p$x[1], y = p$y[1], label = if (ko) "0" else "TF", colour = if (ko) coral else "white",
                size = FS * PT, fontface = "bold"))
}
circ_poly <- function(x, y, r, n = 80) data.table(px = x + r * cos(seq(0, 2 * pi, length.out = n)), py = y + r * sin(seq(0, 2 * pi, length.out = n)))
ggforce_circle <- function(p) {
  d <- p[, circ_poly(x, y, r), by = .(id = seq_len(nrow(p)), fill, stroke)]
  geom_polygon(data = d, aes(px, py, group = id, fill = fill, colour = stroke), linewidth = 0.5, inherit.aes = FALSE)
}
shades <- c("#EAF2F8", "#BDD4E7", "#7DA9CB", "#356FA6")
vals <- matrix(c(1, 2, 4, 2, 1, 3, 2, 3, 4, 3, 1, 2, 1, 2, 3, 2, 4, 4, 1, 2, 3, 1, 2, 3, 4, 2, 1, 3, 2, 4), 5, 6)
mat <- CJ(i = 1:5, j = 1:6)[, .(x = 4 + j * 2.6, y = 111 - i * 2.9, fill = shades[vals[cbind(i, j)]])]
hx <- seq(30, 48, length.out = 121); ph <- seq(pi / 2, 9 * pi / 2, length.out = 121)
helix <- data.table(x = hx, up = 102.3 + 2.6 * sin(ph), lo = 102.3 - 2.6 * sin(ph))
rungs <- helix[seq(11, 111, by = 10)]
txt <- function(x, y, s, bold = FALSE, col = INK) annotate("text", x = x, y = y, label = s, size = FS * PT, fontface = if (bold) "bold" else "plain",
                                                           colour = col, lineheight = 0.88)
darrow <- function(y0, y1) annotate("segment", x = W6 / 2, xend = W6 / 2, y = y0, yend = y1, colour = teal, linewidth = 0.7,
                                    arrow = arrow(length = unit(1.4, "mm"), type = "closed"))
pa <- ggplot() +
  txt(W6 / 2, 117, "Expression + promoter prior", TRUE) +
  annotate("rect", xmin = 5.3, xmax = 22.1, ymin = 94.2, ymax = 110.2, fill = "#F5F9FC") +
  geom_tile(data = mat, aes(x, y, fill = fill), width = 2.2, height = 2.5) +
  txt(26, 102.3, "+", TRUE, coral) +
  geom_segment(data = rungs, aes(x = x, xend = x, y = up, yend = lo), colour = "#A8C0D4", linewidth = 0.4) +
  geom_path(data = helix, aes(x, up), colour = "#2A6EAB", linewidth = 0.55) + geom_path(data = helix, aes(x, lo), colour = "#5C91BA", linewidth = 0.55) +
  txt(13.7, 89.5, "RNA\nexpression", col = MUTED) + txt(39, 89.5, "hg38 promoter\nmotifs", col = MUTED) +
  darrow(84.5, 79.5) +
  txt(W6 / 2, 75.5, "Cell-type-specific GRNs", TRUE) +  # 11-09: wording of the CellOracle paper; lung AM is a subtype, Microglia 2 a SEA-AD supertype
  txt(W6 / 2, 68.5, "Fitted separately in lung AM\nand brain Microglia 2", col = MUTED) +
  net(W6 / 2 + 1, 52) +
  darrow(40.5, 35.5) +
  txt(W6 / 2, 31.5, "Virtual knockout: TF = 0", TRUE) +
  net(W6 / 2 + 1, 17, ko = TRUE) +
  txt(W6 / 2, 3, "Three propagation rounds", col = MUTED) +
  scale_fill_identity() + scale_colour_identity(na.value = NA) +
  coord_equal(xlim = c(0, W6), ylim = c(0, H6), expand = FALSE, clip = "off") + theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0.8, 0.5, 0.5, 0.5, "mm"))

# ---- b: 357 TFs (old 6b: two organ panels, labelled top TFs) -----------------------------------------------------------
rk <- fread(P(file.path(CO, "all357_TF_ranking.tsv")))
stopifnot(nrow(rk) == 357, rk[TF == "STAT3"]$shared_down == 27, all(rk$n_evaluable == 136))
top <- rk[rank_shared_down <= 10]$TF; stopifnot(length(top) == 10)
# 11-09: the hand-picked green "organ-specific" TFs (KLF4, PPARG, MEF2C, BHLHE41) were removed. They followed no ranking rule
# (40 TFs reach the 90th response percentile in only one organ) and the text no longer refers to them. Only the top 10 are coloured.
bl <- rbind(rk[, .(TF, shared_down, pct = response_percentile_lung, deg = degree_lung, organ = "Lung AM")],
            rk[, .(TF, shared_down, pct = response_percentile_brain, deg = degree_brain, organ = "Brain Microglia 2")])
bl[, organ := factor(organ, c("Lung AM", "Brain Microglia 2"))][, is_top := TF %in% top]
bl[, cls := fifelse(is_top, "Top 10 by shared decreases", "Other")]
stopifnot(bl[cls != "Other", .N] == 20)
# hand-set labels: x, y = point; lx, ly = label centre; kx = x units per character, hh = half label height in y units.
# The leader runs from the edge of the label box to the point. No repel, so labels stay off the points and zero lines.
leaders <- function(d, kx, hh) {
  d <- copy(d)[, `:=`(hw = nchar(label) * kx / 2, dx = x - lx, dy = y - ly)]
  d[, s := pmin(hw / pmax(abs(dx), 1e-9), hh / pmax(abs(dy), 1e-9))]
  d[, `:=`(sx = lx + s * dx, sy = ly + s * dy)][]
}
lab_b <- rbind(
  data.table(organ = "Lung AM", label = c("MEF2A", "ELF1", "JUNB", "FOS", "KLF6", "KLF7", "RREB1", "CEBPB", "STAT3", "LYL1"),
             lx = c(14.6, 19, 23, 27, 12.8, 16.9, 21.4, 27.1, 28.7, 26.8), ly = c(rep(122, 4), rep(70, 4), 84, 56)),
  data.table(organ = "Brain Microglia 2", label = c("MEF2A", "KLF7", "ELF1", "CEBPB", "LYL1", "STAT3", "FOS", "JUNB", "RREB1", "KLF6"),
             lx = c(14.8, 19.6, 23.9, 25.2, 24, 28.7, 28, 24.4, 19.5, 16.3), ly = c(rep(122, 3), 108, 94, 85, 69, 66, 66, 77)))
lab_b <- merge(lab_b, bl[cls != "Other", .(organ = as.character(organ), label = TF, x = shared_down, y = pct, cls)], by = c("organ", "label"))
stopifnot(nrow(lab_b) == bl[cls != "Other", .N])
lab_b <- leaders(lab_b, kx = 1, hh = 5.5)[, organ := factor(organ, levels(bl$organ))]
pb <- ggplot(bl, aes(shared_down, pct)) +
  geom_point(data = bl[cls == "Other"], aes(size = deg), colour = greyp, alpha = 0.75, stroke = 0) +
  geom_segment(data = lab_b, aes(sx, sy, xend = x, yend = y, colour = cls), linewidth = 0.25, show.legend = FALSE) +
  geom_point(data = bl[cls != "Other"], aes(size = deg, colour = cls), stroke = 0) +
  geom_text(data = lab_b, aes(lx, ly, label = label, colour = cls), size = FS * PT, family = "Arial", fontface = "bold", show.legend = FALSE) +
  facet_wrap(~organ, ncol = 1) +
  scale_colour_manual(values = c("Top 10 by shared decreases" = blue), name = NULL) +
  scale_size_area(max_size = 3.2, breaks = c(50, 150, 250), name = "Out-degree") +
  scale_y_continuous(breaks = c(0, 50, 100)) + scale_x_continuous(limits = c(0, 31)) + coord_cartesian(ylim = c(0, 128), clip = "off") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.4), ncol = 1), size = guide_legend(order = 2, nrow = 1)) +
  labs(x = "Axis A genes decreased in both models", y = "Within-model response percentile", title = "357 tested TFs") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.3, "mm"), legend.box.spacing = unit(1, "mm"), legend.justification = "left",
        strip.text = element_text(hjust = 0), panel.spacing.y = unit(2.5, "mm"),
        panel.grid.major = element_line(colour = "#E4E8EC", linewidth = 0.3))

# ---- c: top 10 (old 6c lollipop, three counts per TF) --------------------------------------------------------------------
tc <- melt(rk[TF %in% top, .(TF, rank_shared_down, `Down in both` = shared_down, `Up in both` = shared_up, Opposite = opposite)],
           id.vars = c("TF", "rank_shared_down"), variable.name = "cls", value.name = "n")
lv <- rk[TF %in% top][order(rank_shared_down, -shared_down, TF)]$TF
tc[, TF := factor(TF, rev(lv))]
pc <- ggplot(tc, aes(n, TF)) +
  geom_segment(data = tc[, .(mx = max(n)), by = TF], aes(x = 0, xend = mx, y = TF, yend = TF), colour = "#C9CED4", linewidth = 0.5) +
  geom_point(aes(colour = cls), size = 2) +
  geom_text(data = tc[cls == "Down in both"], aes(label = n), colour = dblue, size = FS * PT, fontface = "bold", nudge_y = 0.42) +
  scale_colour_manual(values = c("Down in both" = blue, "Up in both" = "#F04438", Opposite = greyp), name = NULL) +
  scale_x_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Axis A genes (of 136)", y = NULL, title = "Top 10 TFs") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0), axis.ticks.y = element_blank(),
        legend.box.spacing = unit(1, "mm"), legend.justification = "left",
        axis.text.y = element_text(colour = ifelse(rev(lv) %in% c("STAT3", "CEBPB"), dred, INK)))  # 11-09: the two candidates of the text

# ---- d-f: per-gene predicted changes (old zero-line scatter) ---------------------------------------------------------------
gx <- fread(P(file.path(CO, "Fig6def_source_signed_logratio.tsv")))
gx[, cls := fifelse(lung_delta < 0 & brain_delta < 0, "Down in both", fifelse(lung_delta > 0 & brain_delta > 0, "Up in both", "Other / zero"))]
# 11-09: d, e show only the two candidates carried through the text (STAT3, CEBPB); the LYL1 panel inherited from the 10-05
# figure was dropped. The panels are now 93.5 mm wide instead of 63 mm, so kx (x units per character) is scaled by about 0.62.
lab_g <- rbind(data.table(TF = "STAT3", gene = c("HCK", "TFRC", "APOC1", "CSTB"), lx = c(-0.2, -0.335, 0.29, 0.19), ly = c(-0.21, -0.14, -0.045, 0.08), kx = 0.0235, hh = 0.03),
               data.table(TF = "CEBPB", gene = c("FCER1G", "ALOX5AP", "ITGAM", "ITGAX", "PILRA"), lx = c(-0.75, -0.3, -0.4, -0.2, -0.33), ly = c(-0.125, -0.235, -0.07, 0.125, 0.05), kx = 0.025, hh = 0.02))
scat <- function(tf, legend = FALSE) {
  d <- gx[TF == tf]; stopifnot(nrow(d) == 136)
  rk1 <- rk[TF == tf]; stopifnot(sum(d$cls == "Down in both") == rk1$shared_down, sum(d$cls == "Up in both") == rk1$shared_up)
  d[, cls := factor(cls, c("Down in both", "Up in both", "Other / zero"))]
  lab <- merge(lab_g[TF == tf], d[, .(gene, x = lung_delta, y = brain_delta, cls)], by = "gene"); stopifnot(nrow(lab) == lab_g[TF == tf, .N])
  lab <- leaders(lab[, label := gene], lab$kx[1], lab$hh[1])
  ggplot(d, aes(lung_delta, brain_delta)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = INK) + geom_vline(xintercept = 0, linewidth = 0.4, colour = INK) +
    geom_point(data = d[cls == "Other / zero"], aes(colour = cls), size = 0.9, stroke = 0) +
    geom_point(data = d[cls != "Other / zero"], aes(colour = cls), size = 1.4, stroke = 0) +
    geom_segment(data = lab, aes(sx, sy, xend = x, yend = y, colour = cls), linewidth = 0.25, show.legend = FALSE) +
    geom_text(data = lab, aes(lx, ly, label = label, colour = cls), size = FS * PT, family = "Arial", fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = c("Down in both" = blue, "Up in both" = "#F04438", "Other / zero" = greyp), drop = FALSE, name = NULL,
                        guide = if (legend) guide_legend(override.aes = list(size = 2)) else "none") +
    scale_x_continuous(breaks = function(l) pretty(l, 3), labels = function(b) sub("^-", "-", b)) +
    scale_y_continuous(breaks = function(l) pretty(l, 3), labels = function(b) sub("^-", "-", b)) +
    labs(x = "Lung AM delta_X", y = "Brain Microglia 2 delta_X", title = sprintf("%s: 136 Axis A genes", tf),
         subtitle = sprintf("%d down, %d up in both", rk1$shared_down, rk1$shared_up)) +
    theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"), legend.key.size = unit(2.6, "mm"))
}
pd <- scat("STAT3"); pe <- scat("CEBPB")

# ---- f (g before 11-09): TF activity in donor pseudobulk ------------------------------------------------------------------
ta <- fread(cmd = paste("gzip -cd", shQuote(P(file.path(CO, "tf_activity_all.tsv.gz")))))
ta <- ta[grepl("^primary", model) & network == "CollecTRI" & TF %in% top]
stopifnot(nrow(ta) == 30, min(ta$FDR_top10) > 0.2)
ta[, cls := factor(c(Myeloid = "Lung myeloid", `Microglia 2` = "Brain Microglia 2", `Microglia overall` = "Brain all microglia")[class],
                   c("Lung myeloid", "Brain Microglia 2", "Brain all microglia"))]
ta[, hl := TF %in% c("STAT3", "CEBPB")]
lab_t <- data.table(label = c("STAT3", "STAT3", "STAT3", "CEBPB", "CEBPB", "CEBPB"),
                    cls = c("Brain all microglia", "Brain Microglia 2", "Lung myeloid", "Brain all microglia", "Brain Microglia 2", "Lung myeloid"),
                    lx = c(0.47, 0.655, 0.44, 0.675, 0.68, 0.36), ly = c(1.65, 1.4, -0.86, 0.62, -0.22, 0.52))
lab_t <- merge(lab_t, ta[hl == TRUE, .(label = TF, cls = as.character(cls), x = r_with_axisA_score, y = t)], by = c("label", "cls")); stopifnot(nrow(lab_t) == 6)
lab_t <- leaders(lab_t, kx = 0.023, hh = 0.18)
pg <- ggplot(ta, aes(r_with_axisA_score, t)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = INK) + geom_vline(xintercept = 0, linewidth = 0.4, colour = INK) +
  geom_point(aes(shape = cls, colour = hl), size = 1.9, stroke = 0.7) +
  geom_segment(data = lab_t, aes(sx, sy, xend = x, yend = y), inherit.aes = FALSE, linewidth = 0.25, colour = COL$axisA) +
  geom_text(data = lab_t, aes(lx, ly, label = label), inherit.aes = FALSE, size = FS * PT, family = "Arial", fontface = "bold", colour = COL$axisA) +
  scale_shape_manual(values = c(16, 17, 1), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = COL$axisA, `FALSE` = "#7D8590"), guide = "none") +
  scale_x_continuous(limits = c(-0.45, 0.8), breaks = c(-0.4, 0, 0.4, 0.8), labels = c("-0.4", "0", "0.4", "0.8")) +
  scale_y_continuous(labels = function(b) sub("^-", "-", b)) +
  labs(x = "Correlation of TF activity with Axis A score (r)", y = "t for COPD status\nor local CPS", title = "TF activity (decoupleR, top 10 TFs)") +
  theme(legend.position = "right", plot.margin = margin(4, 1.5, 1.5, 1.5, "mm"))

fwrite(bl, P(file.path(SRC, "Fig8b_TF_breadth.tsv")), sep = "\t"); fwrite(ta, P(file.path(SRC, "Fig8f_TF_activity.tsv")), sep = "\t")

compose_figure(list(
  blk(pa, 0, 0, 50, 108, "a"),
  blk(pb, 52, 0, 82, 108, "b"),
  blk(pc, 136, 0, 54, 108, "c"),
  blk(pd, 0, 110, 93.5, 56, "d"), blk(pe, 96.5, 110, 93.5, 56, "e"),
  blk(pg, 0, 168, 150, 54, "f")
), "Fig8", 190, 222, row_groups = list(c("d", "e")),
exemptions = list(list(panels = list("a"), checks = list("row", "column", "panel-width"), reason = "schematic drawn on an equal-scale canvas"),
                  list(panels = list("d", "e"), checks = list("panel-width"), reason = "kept as a row group; widths follow y-axis label widths only")))
