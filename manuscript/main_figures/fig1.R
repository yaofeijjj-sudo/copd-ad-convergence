# Fig 1. Study overview. Layout follows published overviews (organ or data sources converging on a central result, then analysis
# steps with small pictograms; e.g. Morabito 2020 Hum Mol Genet Fig 1, PLOS ONE pone.0349155 Fig 1A). No result is stated.
# Thumbnails reuse the data of Fig 3-5 and 7 without axes; the helix and the regulatory network are those of Fig 8a.
# All shapes are drawn here, so no third-party artwork is used. Palette: lung #B34F4B, brain #397A9B, Axis A #F04438, teal arrows.
source("06_resubmission_20261008/code/main_figures/common.R")
source("06_resubmission_20261008/code/main_figures/fig1_icons.R")
SC <- "cell_localization/SYN01_Entrez_downstream_20261001/source_data"
TEAL <- "#179C8E"; ARROW <- INK; W <- 190; H <- 170; LINE <- "#D9DEE1"; GREY <- "#C9CED4"; CTRL <- "#8A929B"
LF <- "#FAF1F0"; BF <- "#EEF4F7"
txt <- function(x, y, s, colr = INK, face = "bold", hj = 0.5, vj = 0.5) annotate("text", x = x, y = y, label = s, size = FS * PT, family = "Arial", fontface = face, colour = colr, hjust = hj, vjust = vj, lineheight = 0.92)
arw <- function(x, y, xe, ye, lw = 0.8) annotate("segment", x = x, y = y, xend = xe, yend = ye, colour = ARROW, linewidth = lw, arrow = arrow(length = unit(1.8, "mm"), type = "closed"))
rrect <- function(x, y, w, h, r = 2.2, n = 12) {
  a <- function(cx, cy, s, e) { t <- seq(s, e, length.out = n); data.table(x = cx + r * cos(t), y = cy + r * sin(t)) }
  rbind(a(x + w - r, y + r, -pi / 2, 0), a(x + w - r, y + h - r, 0, pi / 2), a(x + r, y + h - r, pi / 2, pi), a(x + r, y + r, pi, 3 * pi / 2))
}
card <- function(x, y, w, h, fill = "white", line = LINE, lwd = 0.5) { d <- rrect(x, y, w, h); annotate("polygon", x = d$x, y = d$y, fill = fill, colour = line, linewidth = lwd) }
rescale <- function(v, a, b, lim = range(v)) a + (v - lim[1]) / diff(lim) * (b - a)
g <- ggplot() + coord_equal(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") + theme_void(base_family = "Arial") +
  scale_colour_identity() + scale_fill_identity() + theme(plot.margin = margin(0, 0, 0, 0))

# ---- module strip: one chip per WGCNA module, the source module outlined ----------------------------------------------------
mix <- function(col, bg, f) { m <- (1 - f) * grDevices::col2rgb(col) + f * as.vector(grDevices::col2rgb(bg)); apply(m, 2, function(v) grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)) }
chips <- function(x0, y0, cols, hi, bg, w = 3.4, gap = 0.5, per = 8) {
  d <- data.table(col = cols)[, `:=`(i = (.I - 1) %% per, r = (.I - 1) %/% per)]
  d[, `:=`(x = x0 + i * (w + gap), y = y0 - r * (w + gap), shade = mix(col, bg, 0.72))]  # other modules recede towards the card colour
  s <- d[col == hi]; stopifnot(nrow(s) == 1)
  list(geom_tile(data = d[col != hi], aes(x + w / 2, y + w / 2, fill = shade), width = w, height = w, colour = NA),
       geom_tile(data = s, aes(x + w / 2, y + w / 2, fill = col), width = w + 1.7, height = w + 1.7, colour = INK, linewidth = 0.9))
}
lmods <- c("red", "blue", "yellow", "black", "brown", "green", "pink", "turquoise")
bmods <- c("turquoise", "midnightblue", "lightcyan", "tan", "cyan", "red", "greenyellow", "green", "salmon", "yellow", "purple", "black", "blue", "pink", "brown", "magenta")
rgbcol <- function(x) vapply(x, function(k) grDevices::rgb(t(grDevices::col2rgb(k)), maxColorValue = 255), "")

# ---- top: discovery -----------------------------------------------------------------------------------------------------------
yC <- 145
g <- g + card(2, 124, 66, 44, LF, COL$lung, 0.6) + card(122, 124, 66, 44, BF, COL$brain, 0.6)
for (l in lungs_icon(16, 143, 11)) g <- g + l
for (l in brain_icon(136, 147.5, 12)) g <- g + l
g <- g + txt(31, 157, "COPD lung", COL$lung, hj = 0) + txt(31, 151.5, "GSE76925, n = 147", MUTED, "plain", hj = 0) + txt(31, 146, "WGCNA: 8 modules", INK, "plain", hj = 0) +
  txt(152, 157, "AD DLPFC", COL$brain, hj = 0) + txt(152, 151.5, "GSE33000, n = 459", MUTED, "plain", hj = 0) + txt(152, 146, "WGCNA: 16 modules", INK, "plain", hj = 0) +
  chips(32, 138, rgbcol(lmods), rgbcol("red"), LF, w = 3.2) + chips(152, 138, rgbcol(bmods), rgbcol("blue"), BF, w = 3.2) +
  txt(31, 129.5, "Source module: red", COL$lung, hj = 0) + txt(152, 129.5, "Source module: blue", COL$brain, hj = 0)
# Venn of the two source modules
cc <- function(x, y, r) { t <- seq(0, 2 * pi, length.out = 180); data.table(x = x + r * cos(t), y = y + r * sin(t)) }
r0 <- 13; xl <- 88.3; xr <- 101.7
A <- cc(xl, yC, r0); B <- cc(xr, yC, r0)
th <- acos((xr - xl) / 2 / r0); ang <- seq(-th, th, length.out = 80)
lens <- rbind(data.table(x = xl + r0 * cos(ang), y = yC + r0 * sin(ang)), data.table(x = xr - r0 * cos(rev(ang)), y = yC + r0 * sin(rev(ang))))
g <- g + annotate("polygon", x = A$x, y = A$y, fill = "#FBE3E1", colour = COL$lung, linewidth = 0.6) +
  annotate("polygon", x = B$x, y = B$y, fill = "#E1EBF1", colour = COL$brain, linewidth = 0.6) +
  annotate("polygon", x = lens$x, y = lens$y, fill = COL$axisA, colour = NA) +
  txt(95, yC, "137", "white") + txt(81, yC, "red", COL$lung) + txt(109, yC, "blue", COL$brain) +
  txt(95, 165.5, "Axis A", COL$axisA) + txt(95, 127.5, "Fixed gene set", MUTED, "plain") +
  arw(68.4, yC, 74.3, yC, 0.9) + arw(121.6, yC, 115.7, yC, 0.9)

# ---- analysis cards ---------------------------------------------------------------------------------------------------------------
CW <- 56; CH <- 50; X <- c(2, 67, 132); Y <- c(60, 3); YA <- Y[1] + CH / 2; YB <- Y[2] + CH / 2
cards <- data.table(no = 1:6, cx = c(X, rev(X)), cy = rep(Y, each = 3), fig = paste("Fig", 3:8),
  title = c("Replication", "Expressing cells", "Tissue-level increase", "Function and genetics", "Blood and other diseases", "Upstream regulators"),  # 11-09: card 3 uses the fixed term of the text
  data = c("GSE47460; 6 frontal cortex datasets", "GSE302339; SEA-AD DLPFC", "Single-cell data; 8 tissue cohorts", "Gene sets; AD and COPD GWAS", "Blood datasets; other disease cohorts", "Alveolar macrophages; Microglia 2"),
  meth = c("Preservation, sigQC, singscore", "UCell, matched random gene sets", "fry, CAMERA, CIBERSORTx", "Enrichment, signatures, MAGMA", "singscore, CAMERA", "CellOracle, decoupleR"))
# order of the analyses: Axis A to 1, left to right along the top row, down from 3 to 4, then right to left to 5 and 6
g <- g + annotate("path", x = c(95, 95, X[1] + CW / 2, X[1] + CW / 2), y = c(124.6, 117, 117, 114), colour = ARROW, linewidth = 0.8, linejoin = "round") +
  arw(X[1] + CW / 2, 114.5, X[1] + CW / 2, 110.4) + arw(X[1] + CW + 0.4, YA, X[2] - 0.4, YA) + arw(X[2] + CW + 0.4, YA, X[3] - 0.4, YA) +
  arw(X[3] + CW / 2, Y[1] - 0.2, X[3] + CW / 2, Y[2] + CH + 0.4) +
  arw(X[3] - 0.4, YB, X[2] + CW + 0.4, YB) + arw(X[2] - 0.4, YB, X[1] + CW + 0.4, YB)
for (i in seq_len(nrow(cards))) {
  k <- cards[i]; bx <- cc(k$cx + 5, k$cy + CH - 5, 2.6)
  g <- g + card(k$cx, k$cy, CW, CH) + annotate("polygon", x = bx$x, y = bx$y, fill = TEAL, colour = NA) + txt(k$cx + 5, k$cy + CH - 5, as.character(k$no), "white") +
    txt(k$cx + 9.5, k$cy + CH - 5, k$title, hj = 0) + txt(k$cx + CW - 2.5, k$cy + CH - 5, k$fig, COL$axisA, hj = 1) +
    txt(k$cx + 3, k$cy + 9, k$data, MUTED, "plain", hj = 0) + txt(k$cx + 3, k$cy + 4, k$meth, INK, "plain", hj = 0)
}
# pictogram area of a card: x from cx + 3 to cx + CW - 3, y from cy + 13 to cy + CH - 9
pa <- function(no) { k <- cards[no]; list(x0 = k$cx + 3, x1 = k$cx + CW - 3, y0 = k$cy + 13, y1 = k$cy + CH - 9) }

# 1 replication: real forest of Fig 3c, without text
f <- fread(P(file.path(SRC, "Fig3c_forest.tsv"))); p <- pa(1); xm <- function(v) rescale(v, p$x0 + 9, p$x1 - 9, c(-1, 1.9))
f[, y := rescale(rev(seq_len(.N)), p$y0 + 1.5, p$y1 - 1.5)][, `:=`(fl = xm(ci_low), fh = xm(ci_high), fe = xm(effect_SD))]
pool <- f[.N]
g <- g + annotate("segment", x = xm(0), xend = xm(0), y = p$y0, yend = p$y1, linetype = 2, linewidth = 0.35, colour = COL$grey) +
  geom_segment(data = f[-.N], aes(x = fl, xend = fh, y = y, yend = y), colour = COL$target, linewidth = 0.5) +
  geom_point(data = f[-.N], aes(fe, y), colour = COL$target, shape = 15, size = 1.4) +
  annotate("polygon", x = xm(c(pool$ci_low, pool$effect_SD, pool$ci_high, pool$effect_SD)), y = pool$y + c(0, 0.9, 0, -0.9), fill = COL$target)

# 2 expressing cells: real UMAPs, target class in red
rd <- function(f, n) { set.seed(1); d <- fread(cmd = paste("gzip -cd", shQuote(P(file.path(SC, f)))), select = c("group", "UMAP_1", "UMAP_2")); d[sample(.N, n)] }
lu <- rd("Lung_UMAP_cells.tsv.gz", 9000); bu <- rd("Brain_UMAP_cells.tsv.gz", 9000); p <- pa(2)
mu <- function(d, tg, x0, x1, y0, y1) { s <- min((x1 - x0) / diff(range(d$UMAP_1)), (y1 - y0) / diff(range(d$UMAP_2)))
  d[, `:=`(px = (x0 + x1) / 2 + (UMAP_1 - mean(range(UMAP_1))) * s, py = (y0 + y1) / 2 + (UMAP_2 - mean(range(UMAP_2))) * s, col = fifelse(group == tg, COL$axisA, GREY))][order(col == COL$axisA)] }
xm2 <- (p$x0 + p$x1) / 2
ul <- mu(lu, "Myeloid", p$x0, xm2 - 1, p$y0, p$y1 - 4.5); ub <- mu(bu, "Immune", xm2 + 1, p$x1, p$y0, p$y1 - 4.5)
g <- g + ggrastr::rasterise(geom_point(data = rbind(ul, ub), aes(px, py, colour = col), size = 0.16, stroke = 0, shape = 16), dpi = 600, dev = "ragg") +
  txt((p$x0 + xm2) / 2, p$y1 - 1, "Lung myeloid", COL$lung) + txt((xm2 + p$x1) / 2, p$y1 - 1, "Brain immune", COL$brain)

# 3 within-cell strip and score against proportion (GSE47460)
s5 <- fread(P(file.path(SRC, "Fig5e_samples.tsv")))[cohort == "GSE47460"]; dn <- fread(P(file.path(SRC, "Fig5a_lung_donors.tsv"))); p <- pa(3)
s5[, `:=`(px = rescale(proportion, p$x0 + 27, p$x1 - 1), py = rescale(score_z, p$y0 + 1, p$y1 - 5.5), col = fifelse(disease == 1, COL$axisA, CTRL))]
set.seed(2); dn[, `:=`(px = rescale(median_axisA_ucell, p$x0 + 1, p$x0 + 20), py = fifelse(disease_group == "COPD", p$y0 + 5, p$y1 - 10) + runif(.N, -1.6, 1.6), col = fifelse(disease_group == "COPD", COL$axisA, CTRL))]
g <- g + geom_point(data = dn, aes(px, py, colour = col), size = 0.7, stroke = 0) + geom_point(data = s5, aes(px, py, colour = col), size = 0.45, stroke = 0, alpha = 0.85) +
  annotate("segment", x = p$x0 + 27, xend = p$x0 + 27, y = p$y0, yend = p$y1 - 4.5, colour = INK, linewidth = 0.4) +
  annotate("segment", x = p$x0 + 27, xend = p$x1, y = p$y0, yend = p$y0, colour = INK, linewidth = 0.4) +
  annotate("segment", x = p$x0, xend = p$x0 + 21, y = p$y0, yend = p$y0, colour = INK, linewidth = 0.4) +
  txt(p$x0 + 10.5, p$y1 - 1, "Within cells") + txt(p$x0 + 38.5, p$y1 - 1, "In tissue")

# 4 bubbles and DNA helix
p <- pa(4); ph <- seq(pi / 2, 7 * pi / 2, length.out = 121); hx <- seq(p$x0 + 27, p$x1 - 1, length.out = 121); yh <- (p$y0 + p$y1 - 4.5) / 2
he <- data.table(x = hx, up = yh + 6 * sin(ph), lo = yh - 6 * sin(ph)); ru <- he[seq(6, 116, by = 11)]
bub <- data.table(x = p$x0 + c(16, 11, 10, 11, 4, 3.5, 4.5, 3), y = rescale(8:1, p$y0 + 1.2, p$y1 - 6), s = c(4.5, 3.6, 3.4, 3.6, 2.2, 1.8, 2.3, 1.8))
bub[, x00 := p$x0]
g <- g + geom_segment(data = bub, aes(x = x00, xend = x, y = y, yend = y), colour = "#D8DCF3", linewidth = 0.4) +
  geom_point(data = bub, aes(x, y, size = s), shape = 21, fill = "#4A3FD0", colour = INK, stroke = 0.3) + scale_size_identity() +
  geom_segment(data = ru, aes(x = x, xend = x, y = up, yend = lo), colour = "#A8C0D4", linewidth = 0.4) +
  geom_path(data = he, aes(x, up), colour = "#2A6EAB", linewidth = 0.7) + geom_path(data = he, aes(x, lo), colour = "#5C91BA", linewidth = 0.7) +
  txt(p$x0 + 10, p$y1 - 1, "Function") + txt(p$x0 + 38.5, p$y1 - 1, "Genetics")

# 5 blood drop and forest of Fig 7a
fb <- fread(P(file.path(SRC, "Fig7a_blood.tsv"))); p <- pa(5); xf <- function(v) rescale(v, p$x0 + 20, p$x1 - 3, c(-0.5, 1.1))
fb[, ry := rescale(rank(y), p$y0 + 1.5, p$y1 - 1.5)][, `:=`(fl = xf(lo), fh = xf(hi), fe = xf(est))]
t0 <- seq(0, 2 * pi, length.out = 120); cxd <- p$x0 + 8; cyd <- (p$y0 + p$y1) / 2 - 2
drop <- data.table(x = cxd + 6.4 * sin(t0) * sin(t0 / 2)^2, y = cyd + 2 + 9.5 * cos(t0))
g <- g + annotate("polygon", x = drop$x, y = drop$y, fill = "#D9766B", colour = NA) +
  annotate("segment", x = xf(0), xend = xf(0), y = p$y0, yend = p$y1, linetype = 2, linewidth = 0.35, colour = COL$grey) +
  geom_segment(data = fb, aes(x = fl, xend = fh, y = ry, yend = ry, colour = col), linewidth = 0.5) +
  geom_point(data = fb, aes(fe, ry, colour = col), shape = 15, size = 1.4)

# 6 regulatory network, as in Fig 8a
p <- pa(6); cxn <- (p$x0 + p$x1) / 2; cyn <- (p$y0 + p$y1) / 2
nd <- data.table(x = cxn + c(-14, 1, 12, 1, -13, 13), y = cyn + c(0, 8.2, 4, -7.5, -8, -6.5), r = c(4.1, 3, 3, 3, 3, 3), fill = c("#E2645A", rep("#DCEAF5", 5)))
ed <- data.table(a = c(1, 1, 1, 2, 5, 4), b = c(2, 3, 4, 3, 4, 6))
se <- ed[, { A <- nd[a]; B <- nd[b]; v <- c(B$x - A$x, B$y - A$y); v <- v / sqrt(sum(v^2)); .(x = A$x + v[1] * A$r, y = A$y + v[2] * A$r, xend = B$x - v[1] * (B$r + 0.3), yend = B$y - v[2] * (B$r + 0.3), col = if (a == 1) "#E2645A" else "#8BA1B5") }, by = .(a, b)]
g <- g + geom_segment(data = se, aes(x, y, xend = xend, yend = yend, colour = col), linewidth = 0.5, arrow = arrow(length = unit(1.2, "mm"), type = "closed"))
for (i in seq_len(nrow(nd))) { d <- cc(nd$x[i], nd$y[i], nd$r[i]); g <- g + annotate("polygon", x = d$x, y = d$y, fill = nd$fill[i], colour = NA) }
g <- g + txt(nd$x[1], nd$y[1], "TF", "white")

compose_figure(list(blk(g, 0.2, 0.2, W - 0.4, H - 0.4, "overview", tags = list())), "Fig1", W, H,
               exemptions = list(list(panels = list("overview"), checks = list("row", "column", "panel-width"), reason = "single schematic panel")))
