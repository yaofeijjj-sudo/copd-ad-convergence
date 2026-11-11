# Shared settings for the eight main figures (PLOS ONE: width <= 190.5 mm, height <= 222.25 mm, text 8-12 pt, Arial).
# Run every figure script from the project root with the project R:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 06_resubmission_20261008/code/main_figures/figN.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(grid)
})
ROOT <- normalizePath(".")
P <- function(x) file.path(ROOT, x)
REF <- "05_github_release/copd-ad-convergence/validation/reference"
CAND <- "03_candidate_analysis"
FIGDIR <- "06_resubmission_20261008/figures/source_pdf"  # PDFs from which the submitted Fig1-8.tif were rendered
SRC <- "06_resubmission_20261008/code/main_figures/source_data"
QA <- "06_resubmission_20261008/code/main_figures/qa"
dir.create(P(QA), showWarnings = FALSE, recursive = TRUE)
source(path.expand("~/.claude/skills/nature-figure/scripts/panel_alignment.R"))
# The helper measures on grDevices::pdf(), which has no "Arial" family; Arial is metric-compatible with Helvetica,
# so registering Helvetica metrics under that name gives the same layout as the Cairo export.
if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)

FS <- 8                        # every glyph >= 8 pt
PT <- 1 / .pt                  # ggplot text size unit for geom_text: size = FS * PT
LW <- 0.5
TAG <- 12                      # panel letters, as in the 10-05 Affinity layout
# Visual style of the 10-05 figure.af panels (WGCNA, cross WGCNA, 富集分析, 细胞定位, GWAS和外周血, 虚拟敲除 scripts):
# bold Arial, ink #262626, red #F04438 / blue #168BEA, darker forest red/blue #B34F4B / #397A9B, grey #D1D5DB.
INK <- "#262626"; MUTED <- "#555B61"
COL <- list(
  axisA = "#F04438", axisB = "#168BEA", grey = "#A7ADB4", light = "#E6E8EB", dark = INK, muted = MUTED,
  up = "#F04438", down = "#168BEA", other = "#C9CED4", opposite = "#86B5AA", lung = "#B34F4B", brain = "#397A9B",
  forestA = "#B34F4B", forestB = "#397A9B", hist = "#D1D5DB", histline = "#AEB5BF", grid = "#ECEEF1", blood = "#EA5B5B", target = "#EA5B5B"
)
# per-category palette of the old UMAP and raincloud panels
PAL8 <- c("#F04438", "#F28C45", "#DDB34B", "#52B7A5", "#4C9FC4", "#8571C4", "#D76E98", "#8DA85D")

theme_pub <- function() {
  theme_classic(base_size = FS, base_family = "Arial") +
    theme(
      text = element_text(size = FS, colour = INK, face = "bold"),
      axis.text = element_text(size = FS, colour = INK, face = "bold"),
      axis.title = element_text(size = FS, colour = INK, face = "bold"),
      axis.line = element_line(linewidth = LW, colour = INK),
      axis.ticks = element_line(linewidth = LW, colour = INK),
      axis.ticks.length = unit(1.2, "mm"),
      legend.text = element_text(size = FS, face = "bold"), legend.title = element_text(size = FS, face = "bold"),
      legend.key.size = unit(3, "mm"), legend.background = element_blank(),
      strip.text = element_text(size = FS, face = "bold"), strip.background = element_blank(),
      plot.title = element_text(size = FS, face = "bold", hjust = 0.5, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.subtitle = element_text(size = FS, face = "bold", hjust = 0.5, colour = MUTED, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.tag = element_text(size = TAG, face = "bold"),
      plot.title.position = "plot", legend.location = "plot",
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
      panel.grid = element_blank()
    )
}
theme_set(theme_pub())
tag_theme <- theme(plot.tag = element_text(size = TAG, face = "bold", family = "Arial"))

fmt_p <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("")
    if (x >= 0.001) return(formatC(signif(x, digits), format = "fg", digits = digits))
    s <- formatC(x, format = "e", digits = digits - 1); b <- strsplit(s, "e")[[1]]
    sprintf("%s×10^%d", b[1], as.integer(b[2]))
  }, character(1))
}
# plotmath-free superscript for P-values in text columns: "7.2 × 10⁻⁶"
sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹", "-" = "⁻")
fmt_p_sup <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("")
    if (x >= 0.001) return(formatC(signif(x, digits), format = "fg", digits = digits, flag = "#"))
    s <- formatC(x, format = "e", digits = digits - 1); b <- strsplit(s, "e")[[1]]
    # Arial has no superscript minus or 4-9 (Cairo would substitute Helvetica), so use e-notation in Arial
    sprintf("%se%d", b[1], as.integer(b[2]))
  }, character(1))
}
fmt_ci <- function(est, lo, hi, d = 2) sprintf("%.*f (%.*f to %.*f)", d, est, d, lo, d, hi)

# PDF device: Quartz with the installed Arial, as in the 10-05 panel scripts. Quartz writes each label as one text
# object in TrueType Arial (Cairo splits labels into per-glyph TJ arrays and puts spaces in a Type 3 font); the remaining
# TJ adjustments of <= 1/1000 em are merged afterwards by _workspace/codex_prompts/merge_tj_text_runs.py.
pdf_dev <- function(file, w_in, h_in) grDevices::quartz(type = "pdf", file = file, width = w_in, height = h_in, family = "Arial", bg = "white")
MERGE_PY <- c(P("_workspace/nature_figure_qa_venv/bin/python"), P("_workspace/codex_prompts/merge_tj_text_runs.py"))

# Measure the final plot-area rectangle of every ggplot panel after drawing one block (a ggplot or a patchwork) on a
# device of its own size. Patches are numbered in the order they were added; facet panels get suffixes -1, -2, ...
measure_block <- function(plot, w_mm, h_mm, ids) {
  f <- tempfile(fileext = ".pdf"); pdf_dev(f, w_mm / 25.4, h_mm / 25.4)
  on.exit({ grDevices::dev.off(); unlink(f) })
  print(plot); grid::grid.force()
  v <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  d <- data.table(name = v$name, depth = v$vpDepth, path = v$vpPath)
  d <- d[grepl("^panel", name) & !grepl("^panel-area", name)]
  d[, container := grepl(";", name)]
  d[, patch := if (inherits(plot, "patchwork")) cumsum(depth == min(depth)) else 1L]
  leaves <- d[container == FALSE]
  m <- rbindlist(lapply(seq_len(nrow(leaves)), function(i) {
    parts <- strsplit(leaves$path[i], "::", fixed = TRUE)[[1]]; parts <- c(parts[parts != "ROOT"], leaves$name[i])
    grid::upViewport(0); grid::downViewport(do.call(grid::vpPath, as.list(parts)))
    lo <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc")); hi <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"))
    grid::upViewport(0)
    data.table(patch = leaves$patch[i], l = grid::convertX(lo$x, "pt", TRUE), b = grid::convertY(lo$y, "pt", TRUE),
               r = grid::convertX(hi$x, "pt", TRUE), t = grid::convertY(hi$y, "pt", TRUE))
  }))
  stopifnot(max(m$patch) == length(ids))
  m[, id := ids[patch]][, n := .N, by = patch][n > 1, id := paste0(id, "-", seq_len(.N)), by = patch]
  m[]
}

# blk(): one placed block. x, y = top-left corner in mm from the figure's top-left; ids = panel ids of its patches;
# tags = list of c(label, x_mm, y_mm) relative to the block's top-left.
blk <- function(plot, x, y, w, h, ids, tags = list(c(ids[1], 0, 0))) list(plot = plot, x = x, y = y, w = w, h = h, ids = ids, tags = tags)

compose_figure <- function(blocks, name, W, H, row_groups = NULL, column_groups = NULL, exemptions = list()) {
  stopifnot(W <= 190.5, H <= 222.25)
  base <- P(file.path(QA, name))
  canvas <- cowplot::ggdraw() + theme(plot.background = element_rect(fill = "white", colour = NA))
  for (b in blocks) {
    canvas <- canvas + cowplot::draw_plot(b$plot, x = b$x / W, y = 1 - (b$y + b$h) / H, width = b$w / W, height = b$h / H)
    for (t in b$tags) canvas <- canvas + cowplot::draw_label(t[1], x = (b$x + as.numeric(t[2]) + 0.6) / W, y = 1 - (b$y + as.numeric(t[3]) + 1.2) / H,
                                                            hjust = 0, vjust = 1, fontface = "bold", size = TAG, fontfamily = "Arial")
  }
  # alignment manifest from per-block measurements
  m <- rbindlist(lapply(blocks, function(b) {
    x <- measure_block(b$plot, b$w, b$h, b$ids)
    ox <- b$x * 72 / 25.4; oy <- (H - b$y - b$h) * 72 / 25.4
    x[, `:=`(l = l + ox, r = r + ox, b = b + oy, t = t + oy)]
  }))
  man <- list(schema_version = 1L, backend = "r-ggplot-cowplot-viewports", figure = list(width_pt = W * 72 / 25.4, height_pt = H * 72 / 25.4),
              panels = lapply(seq_len(nrow(m)), function(i) list(id = m$id[i], bbox_pt = c(m$l[i], m$b[i], m$r[i], m$t[i]))),
              exemptions = exemptions)
  if (length(row_groups)) man$row_groups <- lapply(row_groups, as.character)
  if (length(column_groups)) man$column_groups <- lapply(column_groups, as.character)
  jsonlite::write_json(man, paste0(base, ".alignment-layout.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  fwrite(m, paste0(base, ".panels.tsv"), sep = "\t")
  out <- P(file.path(FIGDIR, paste0(name, ".pdf"))); tmp <- tempfile(fileext = ".pdf")
  pdf_dev(tmp, W / 25.4, H / 25.4); print(canvas); invisible(dev.off())
  log <- system2(MERGE_PY[1], c("-I", shQuote(MERGE_PY[2]), shQuote(tmp), shQuote(out)), stdout = TRUE)
  stopifnot(is.null(attr(log, "status")), file.exists(out)); cat(name, "merge:", log, "\n")
  ragg::agg_png(P(file.path(FIGDIR, "预览PNG", paste0(name, "_preview.png"))), width = W / 25.4, height = H / 25.4, units = "in", res = 200)
  print(canvas); invisible(dev.off())
  cat("saved", name, sprintf("%.0f x %.0f mm; %d panels measured\n", W, H, nrow(m)))
}

# ---- table-style forest plot of the 10-05 panels (old 2d / 4d): text columns | forest | text columns ----------------
# d: one row per estimate with y, est, lo, hi, col, shape ("square", "circle", "open", "diamond"), size (mm, optional)
# hdr: group-header rows (y, lab), drawn bold in the first left column; rules: y positions of thin separators.
# left/right: lists of columns, each list(var, head, x in [0, 1], hjust, colour, face).
ftext <- function(d, cols, hdr, ytop, ylim, title = " ", rules = numeric(), hdr_x = 0) {
  body <- rbindlist(lapply(cols, function(cl) d[, .(y, x = cl$x, hj = cl$hjust %||% 0, lab = as.character(get(cl$var)),
                                                     colr = cl$colour %||% INK, face = cl$face %||% "plain")]))
  heads <- rbindlist(lapply(cols, function(cl) data.table(y = ytop, x = cl$x, hj = cl$hjust %||% 0, lab = cl$head)))
  p <- ggplot() +
    geom_text(data = body, aes(x, y, label = lab, hjust = hj, colour = colr, fontface = face), size = FS * PT, family = "Arial",
              lineheight = 0.9) +
    geom_text(data = heads, aes(x, y, label = lab, hjust = hj), size = FS * PT, family = "Arial", fontface = "bold", colour = INK,
              lineheight = 0.9) +
    annotate("segment", x = -Inf, xend = Inf, y = ytop - 0.6, yend = ytop - 0.6, colour = "#B6BEC3", linewidth = 0.4)
  if (nrow(hdr)) p <- p + geom_text(data = hdr, aes(x = hdr_x, y = y, label = lab), hjust = 0, size = FS * PT,
                                    family = "Arial", fontface = "bold", colour = hdr$colr %||% INK)
  for (r in rules) p <- p + annotate("segment", x = -Inf, xend = Inf, y = r, yend = r, colour = "#D9DEE1", linewidth = 0.3)
  p + scale_colour_identity() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(ylim = ylim, clip = "off") + labs(title = title, x = " \n ") + theme_void(base_family = "Arial") +
    theme(plot.title = element_text(size = FS, face = "bold", hjust = 0, colour = INK, margin = margin(0, 0, 1.5, 5, "mm")),
          axis.title.x = element_text(size = FS, lineheight = 1, margin = margin(1.2, 0, 0, 0, "mm")),
          plot.margin = margin(1.5, 0, 1.5, 0, "mm"))
}
fplot <- function(d, xlim, xbreaks, xlab, ytop, ylim, rules = numeric(), caps = TRUE, head = "") {
  pts <- d[shape != "diamond"]; dia <- d[shape == "diamond", .(x = c(lo, est, hi, est), yy = y + c(0, 0.32, 0, -0.32), col), by = y]
  if (!"size" %in% names(pts)) pts[, size := 2]
  p <- ggplot(d, aes(est, y)) +
    annotate("segment", x = -Inf, xend = Inf, y = ytop - 0.6, yend = ytop - 0.6, colour = "#B6BEC3", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = COL$grey) +
    geom_segment(data = pts, aes(x = lo, xend = hi, yend = y, colour = col), linewidth = 0.6)
  if (caps) p <- p + geom_segment(data = rbind(pts[, .(x = lo, y, col)], pts[, .(x = hi, y, col)]),
                                  aes(x = x, xend = x, y = y - 0.14, yend = y + 0.14, colour = col), linewidth = 0.6, inherit.aes = FALSE)
  for (r in rules) p <- p + annotate("segment", x = -Inf, xend = Inf, y = r, yend = r, colour = "#D9DEE1", linewidth = 0.3)
  p + geom_point(data = pts[shape == "square"], aes(colour = col, size = size), shape = 15) +
    geom_point(data = pts[shape == "circle"], aes(colour = col, size = size), shape = 16) +
    geom_point(data = pts[shape == "open"], aes(colour = col, size = size), shape = 21, fill = "white", stroke = 0.7) +
    geom_polygon(data = dia, aes(x, yy, group = y, fill = col), colour = NA, inherit.aes = FALSE) +
    scale_colour_identity() + scale_fill_identity() + scale_size_identity() +
    scale_x_continuous(breaks = xbreaks, labels = function(b) format(b, trim = TRUE, drop0trailing = TRUE)) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    annotate("text", x = mean(xlim), y = ytop, label = head, size = FS * PT, fontface = "bold", family = "Arial", colour = INK) +
    labs(x = xlab, y = NULL, title = " ") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(lineheight = 1, margin = margin(1.2, 0, 0, 0, "mm")),
          plot.margin = margin(1.5, 0, 1.5, 0, "mm"))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
