# Fig 6. Function of Axis A, relation to myeloid signatures, and genetic enrichment.
# Sources (frozen, not recomputed):
#   富集分析/SYN01_Entrez_downstream_20261001/source_data/axisA_functional_enrichment_{complete_family_terms,complete_family_links,selected_genes}.csv (a)
#   03_candidate_analysis/AxisA_myeloid_reference_20261019/results/signature_models_{by_cohort,AD_pooled}.tsv,
#   03_candidate_analysis/AxisA_replication_after_exclusion_20261019/results/{by_cohort,pooled_sensitivity}.tsv (b)
#   REF/magma.csv; 03_candidate_analysis/AxisA_MAGMA_APOE_exclusion_20261016/results/AD_MAGMA_APOE_exclusion.csv;
#   03_candidate_analysis/AxisA_MAGMA_conditional_microglia_20261016/results/conditional_magma_results.csv;
#   03_candidate_analysis/AxisA_auxiliary_checks_20261019/results/{magma_patir_conditional,magma_background_sensitivity}.tsv (c)
source("06_resubmission_20261008/code/main_figures/common.R")
EN <- "富集分析/SYN01_Entrez_downstream_20261001/source_data"

# ---- a: old 3a sankey-bubble (富集分析/桑基气泡图_统一格式_20260917/code/plot_functions.R, make_axis_plot) ----------------
suppressPackageStartupMessages({ library(ggalluvial); library(dplyr); library(scales) })
source(P("富集分析/桑基气泡图_统一格式_20260917/code/plot_functions.R"))
terms <- as.data.frame(fread(P(file.path(EN, "axisA_functional_enrichment_complete_family_terms.csv"))))
links <- as.data.frame(fread(P(file.path(EN, "axisA_functional_enrichment_complete_family_links.csv"))))
sg <- fread(P(file.path(EN, "axisA_functional_enrichment_selected_genes.csv")))
stopifnot(nrow(terms) == 8, nrow(sg) == 18, nrow(links) == 85, all(terms$significant_primary))
tm <- as.data.table(terms)
TNX <- 1.6   # narrower ribbon section than the old 2.22 so that one-line term labels fit at 8 pt
obj <- make_axis_plot("Axis A", label_side = "right", show_plot_legend = TRUE, term_node_x = TNX, legend_layout = "right")
pa <- obj$plot
# term labels re-wrapped at 30 characters (old: 18 at a 240 mm canvas); the old size limits (9-73 genes) predate the
# 137-gene Axis A, whose largest term has 76 genes, so the upper limit is raised to 80 with the old range and breaks
wrap30 <- function(x) vapply(x, function(z) paste(strwrap(z, 30), collapse = "\n"), "")
relabel <- function(d) { d$label_modern <- wrap30(paste0(d$term_name, " (", database_key(d$database), ")")); d }
for (i in seq_along(pa$layers)) if (inherits(pa$layers[[i]]$geom, "GeomLabel")) pa$layers[[i]]$data <- relabel(pa$layers[[i]]$data)
obj$bubble_data <- relabel(obj$bubble_data)
pa <- pa + suppressMessages(scale_fill_gradientn(name = "-log10(BH q)", colours = q_colors, limits = c(1, 8), breaks = c(1, 4, 8),
                                                  labels = c("1", "4", ">=8"), oob = squish, na.value = "white")) +
  suppressMessages(scale_size_continuous(name = "Overlap n", limits = c(9, 80), breaks = c(10, 30, 60), range = c(2.4, 7.8)))
for (i in seq_along(pa$layers)) if (inherits(pa$layers[[i]]$geom, c("GeomText", "GeomLabel"))) {
  pa$layers[[i]]$aes_params$size <- FS * PT; pa$layers[[i]]$aes_params$family <- "Arial"; pa$layers[[i]]$aes_params$fontface <- "bold"
}
pa <- pa + labs(title = "Axis A  |  disease-increased", subtitle = "COPD red and AD blue, 137 shared genes") +
  guides(colour = "none") +
  theme(text = element_text(family = "Arial", size = FS, face = "bold"),
        plot.title = element_text(size = FS, face = "bold", colour = COL$axisA, hjust = 0, margin = margin(b = 1)),
        plot.subtitle = element_text(size = FS, face = "bold", colour = palette[["text_muted"]], hjust = 0, margin = margin(b = 3)),
        legend.title = element_text(size = FS, face = "bold", colour = INK, margin = margin(b = 7)),
        legend.text = element_text(size = FS, face = "bold", colour = palette[["text_muted"]]),
        plot.margin = margin(1.5, 1.5, 1.5, 6, "mm"))
# Connector lines start after each term label, measured on the output device (as in the old draw.R).
WA <- 190; HA <- 108
mf <- tempfile(fileext = ".pdf"); pdf_dev(mf, WA / 25.4, HA / 25.4)
gt <- ggplotGrob(pa); grid.newpage(); grid.draw(gt); grid.force()
loc <- gt$layout[gt$layout$name == "panel", ]
seekViewport(paste0("panel.", loc$t, "-", loc$l, "-", loc$b, "-", loc$r))
panel_mm <- convertWidth(unit(1, "npc"), "mm", valueOnly = TRUE)
xr <- ggplot_build(pa)$layout$panel_params[[1]]$x.range
td <- obj$bubble_data
lw <- vapply(td$label_modern, function(z) max(vapply(strsplit(z, "\n", fixed = TRUE)[[1]], function(l)
  convertWidth(grobWidth(textGrob(l, gp = gpar(fontfamily = "Arial", fontsize = FS, fontface = "bold"))), "mm", valueOnly = TRUE), 0)), 0)
dev.off(); unlink(mf)
td$line_start <- TNX + 0.0375 + 0.0175 + (lw + 1.2) / panel_mm * diff(xr)
stopifnot(all(td$line_start < td$bubble_x))
conn <- which(vapply(pa$layers, function(l) inherits(l$geom, "GeomSegment") && "bubble_x" %in% all.vars(l$mapping$xend), logical(1)))
stopifnot(length(conn) == 1)
pa$layers[[conn]]$data <- td; pa$layers[[conn]]$mapping$x <- aes(x = line_start)$x

# ---- b: myeloid signatures -----------------------------------------------------------------------------
sig <- data.table(signature = c("BRETIGEA_microglia_top250", "Patir2019_core_microglia", "Travaglini2020_lung_macrophage",
                                "GO0006956_complement_activation", "GO0006909_phagocytosis"),
                  lab = c("Microglia (BRETIGEA)", "Core microglia (Patir)", "Lung macrophage (Travaglini)",
                          "Complement activation (GO)", "Phagocytosis (GO)"))
sc <- fread(P(file.path(CAND, "AxisA_myeloid_reference_20261019/results/signature_models_by_cohort.tsv")))[cohort == "GSE47460"]
sp <- fread(P(file.path(CAND, "AxisA_myeloid_reference_20261019/results/signature_models_AD_pooled.tsv")))
rb <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/by_cohort.tsv")))[cohort == "GSE47460"]
rp <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/pooled_sensitivity.tsv")))[analysis == "DL (primary)"]
keep <- c("signature alone", "Axis A adjusted for signature")
sb <- rbind(sc[analysis %in% keep, .(signature, analysis, set = "COPD lung (GSE47460)", effect_SD, ci_low, ci_high)],
            sp[analysis %in% keep, .(signature, analysis, set = "AD pooled (6 cohorts)", effect_SD, ci_low, ci_high)])
ref <- data.table(signature = "AxisA", analysis = "signature alone", set = c("COPD lung (GSE47460)", "AD pooled (6 cohorts)"),
                  effect_SD = c(rb$effect_SD, rp$effect_SD), ci_low = c(rb$ci_low, rp$ci_low), ci_high = c(rb$ci_high, rp$ci_high))
sb <- rbind(ref, sb)
labs_all <- rbind(data.table(signature = "AxisA", lab = "Axis A (unadjusted)"), sig)
sb[, lab := labs_all$lab[match(signature, labs_all$signature)]]
sb[, y := (nrow(labs_all) + 1 - match(signature, labs_all$signature)) + fifelse(set == "COPD lung (GSE47460)", 0.16, -0.16)]
sb[, panel := factor(fifelse(analysis == "signature alone", "Gene set\nalone", "Axis A\nadjusted"),
                     c("Gene set\nalone", "Axis A\nadjusted"))]
sb[signature == "AxisA", panel := "Gene set\nalone"]
stopifnot(abs(sb[signature == "Travaglini2020_lung_macrophage" & panel == "Axis A adjusted"]$effect_SD - c(0.331, 0.220)) < 0.001)
pbx <- ggplot(sb, aes(effect_SD, y, colour = set)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = COL$grey) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.14, linewidth = 0.5, orientation = "y") +
  geom_point(aes(shape = set), size = 1.8) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = c("COPD lung (GSE47460)" = COL$lung, "AD pooled (6 cohorts)" = COL$brain), breaks = c("COPD lung (GSE47460)", "AD pooled (6 cohorts)"), name = NULL) +
  scale_shape_manual(values = c("COPD lung (GSE47460)" = 16, "AD pooled (6 cohorts)" = 15), breaks = c("COPD lung (GSE47460)", "AD pooled (6 cohorts)"), name = NULL) +
  scale_y_continuous(breaks = rev(seq_len(nrow(labs_all))), labels = labs_all$lab, expand = expansion(add = 0.5)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + coord_cartesian(xlim = c(-0.3, 1.05)) +
  labs(x = "Disease effect (SD)", y = NULL, title = "Myeloid signatures and functional gene sets") +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
        panel.spacing.x = unit(3, "mm"), strip.text = element_text(colour = MUTED, lineheight = 0.88))

# r with Axis A, as stated in the text: GSE47460 and the median across the six AD datasets (Table I b, c)
rtab <- rbind(sc[analysis == "signature alone", .(signature, set = "COPD lung (GSE47460)", r = r_with_axisA)],
              sp[analysis == "signature alone", .(signature, set = "AD pooled (6 cohorts)", r = median_r)])
rtab[, y := (nrow(labs_all) + 1 - match(signature, labs_all$signature)) + fifelse(set == "COPD lung (GSE47460)", 0.21, -0.21)]
rtab[, panel := "r"]
pbr <- ggplot(rtab, aes(0, y, colour = set)) +
  geom_text(aes(label = sprintf("%.2f", r)), size = FS * PT, fontface = "bold", family = "Arial") +
  facet_wrap(~panel) +
  scale_colour_manual(values = c("COPD lung (GSE47460)" = COL$lung, "AD pooled (6 cohorts)" = COL$brain), guide = "none") +
  scale_y_continuous(limits = c(0.5, nrow(labs_all) + 0.5), expand = c(0, 0)) + scale_x_continuous(limits = c(-1, 1), breaks = NULL) +
  labs(x = " ", y = NULL, title = " ") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), strip.text = element_text(colour = MUTED, lineheight = 0.88),
        plot.margin = margin(1.5, 0, 1.5, 0, "mm"))
fwrite(rtab[, .(signature, set, r)], P(file.path(SRC, "Fig6b_signature_r.tsv")), sep = "\t")
pb_all <- pbx + pbr + plot_layout(widths = c(10, 2.8))

# ---- c: MAGMA --------------------------------------------------------------------------------------------------------
mg <- fread(P(file.path(REF, "magma.csv")))
apoe <- fread(P(file.path(CAND, "AxisA_MAGMA_APOE_exclusion_20261016/results/AD_MAGMA_APOE_exclusion.csv")))[variant == "V2_whole_region" & axis == "Axis A"]
bre <- fread(P(file.path(CAND, "AxisA_MAGMA_conditional_microglia_20261016/results/conditional_magma_results.csv")))[reference == "top250" & VARIABLE == "AxisA" & model == "conditional on MicrogliaRef"]
pat <- fread(P(file.path(CAND, "AxisA_auxiliary_checks_20261019/results/magma_patir_conditional.tsv")))[VARIABLE == "AxisA" & model == "conditional on PatirCore"]
bkg <- fread(P(file.path(CAND, "AxisA_auxiliary_checks_20261019/results/magma_background_sensitivity.tsv")))[VARIABLE == "AxisA" & disease == "AD" & background == "4,832 network-eligible"]
stopifnot(nrow(apoe) == 1, nrow(bre) == 1, nrow(pat) == 1, nrow(bkg) == 1)
mc <- rbind(
  mg[, .(lab = sprintf("%s × %s", disease, axis), BETA, SE, P, holm = holm_p_across_2x2_four_tests, grp = "main")],
  data.table(lab = c("APOE ±1 Mb removed", "Conditional on BRETIGEA microglia", "Conditional on Patir core microglia", "Network-eligible background"),
             BETA = c(apoe$BETA, bre$BETA, pat$BETA, bkg$BETA), SE = c(apoe$SE, bre$SE, pat$SE, bkg$SE),
             P = c(apoe$P_one_sided, bre$P, pat$P, bkg$P), holm = NA_real_, grp = "sens"))
mc[, ord := c(1, 2, 3, 4, 5, 6, 7, 8)[match(lab, c("AD × Axis A", "AD × Axis B", "COPD × Axis A", "COPD × Axis B",
                                               "APOE ±1 Mb removed", "Conditional on BRETIGEA microglia", "Conditional on Patir core microglia", "Network-eligible background"))]]
mc[, y := c(10, 9, 8, 7, 4, 3, 2, 1)[ord]]
mc[, lo := BETA - 1.96 * SE][, hi := BETA + 1.96 * SE]
mc[, ptxt := fifelse(grp == "main", fmt_p_sup(holm, 2), fmt_p_sup(P, 2))]
mc[, sig := fifelse(grp == "main", holm < 0.05, P < 0.05)]
mc[, col := fifelse(grepl("Axis B", lab), COL$axisB, COL$axisA)]
stopifnot(abs(mc[lab == "AD × Axis A"]$holm - 2.8732e-05) < 1e-9)
# old 4a design: vertical estimates with 95% CI, P above each estimate, pink band behind the Holm P < 0.05 test
mc[, facet := factor(fifelse(grp == "sens", "Sensitivity (AD, Axis A)", fifelse(grepl("^AD", lab), "AD", "COPD liability")),
                     c("AD", "COPD liability", "Sensitivity (AD, Axis A)"))]
mc[, xl := fifelse(grp == "main", sub("^(AD|COPD) \u00d7 ", "", lab),
                   c("APOE \u00b11 Mb removed" = "APOE\nremoved", "Conditional on BRETIGEA microglia" = "Cond.\nBRETIGEA",
                     "Conditional on Patir core microglia" = "Cond.\nPatir", "Network-eligible background" = "Network\nbackground")[lab])]
mc[, xi := as.numeric(factor(xl, unique(xl))), by = facet]
mc[, pcol := fifelse(sig, INK, MUTED)][, colf := fifelse(grepl("Axis B", lab), COL$forestB, COL$forestA)]
mc[, xl := gsub("\n", " ", xl)]
xlabs <- unique(mc[, .(facet, xi, xl)])
phead <- mc[, .(xi = mean(range(xi)), lab = if (grp[1] == "main") "Holm P" else "P"), by = facet]
shade <- mc[grp == "main" & sig == TRUE]
pc <- ggplot(mc, aes(xi, BETA)) +
  geom_rect(data = shade, aes(xmin = xi - 0.45, xmax = xi + 0.45, ymin = -Inf, ymax = Inf), fill = "#FDEBEC", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = COL$grey) +
  geom_errorbar(aes(ymin = lo, ymax = hi, colour = colf), width = 0.16, linewidth = 0.55) +
  geom_point(aes(colour = colf, shape = grepl("Axis B", lab)), size = 1.9) +
  geom_text(aes(y = 0.7, label = ptxt, colour = pcol), size = FS * PT, fontface = "bold", vjust = 0) +
  geom_text(data = phead, aes(x = xi, y = 0.79, label = lab), size = FS * PT, fontface = "bold", vjust = 0, colour = INK, inherit.aes = FALSE) +
  facet_grid(~facet, scales = "free_x", space = "free_x") +
  scale_colour_identity() + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = NULL, expand = expansion(add = 0.55)) +
  geom_text(data = xlabs, aes(x = xi, y = -0.25, label = xl), angle = 35, hjust = 1, vjust = 1, size = FS * PT, fontface = "bold",
            colour = INK, inherit.aes = FALSE) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6), labels = c("-0.2", "0", "0.2", "0.4", "0.6")) +
  coord_cartesian(ylim = c(-0.22, 0.68), clip = "off") +
  labs(x = NULL, y = "MAGMA coefficient (beta)", title = "Genetic enrichment (MAGMA)") +
  theme(axis.ticks.x = element_blank(), axis.line.x = element_line(linewidth = LW, colour = INK), panel.spacing.x = unit(4, "mm"),
        strip.text = element_text(lineheight = 0.85, vjust = 0, margin = margin(0, 0, 8, 0, "mm")), strip.clip = "off",
        plot.margin = margin(1.5, 1.5, 17, 1.5, "mm"))

fwrite(sb, P(file.path(SRC, "Fig6b_signatures.tsv")), sep = "\t"); fwrite(mc, P(file.path(SRC, "Fig6c_magma.tsv")), sep = "\t")

compose_figure(list(
  blk(pa, 0, 0, WA, HA, "a"),
  blk(pb_all, 0, 110, 84, 80, c("b", "b_r"), tags = list(c("b", 0, 0))),
  blk(pc, 85, 110, 105, 80, "c")
), "Fig6", 190, 190, row_groups = list(c("b-1", "b-2", "b_r")),
exemptions = list(list(panels = list("b-1", "b-2", "b_r"), checks = list("panel-width", "horizontal-gutter"), reason = "r column is a narrow text column beside the two facets"),
                  list(panels = list("a"), checks = list("row", "column", "panel-width"), reason = "sankey-bubble panel spans the full width"),
                  list(panels = list("c-1", "c-2", "c-3"), checks = list("panel-width"), reason = "facet widths follow the number of tests")))
