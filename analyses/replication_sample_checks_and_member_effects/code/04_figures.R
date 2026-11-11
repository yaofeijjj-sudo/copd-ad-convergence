# Member direction heatmap and lung-brain member effect scatter (draft supplementary figures; text >= 8 pt).
# Run from project root:
#   R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/04_figures.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
OUT <- "03_candidate_analysis/AxisA_member_direction_robustness_20261018"
mem <- fread(file.path(OUT, "results/member_effects_wide.tsv"))
cols <- c("GSE76925", "GSE33000", "GSE47460", "GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980", "AD_pooled")
lab <- c(GSE76925 = "GSE76925 (lung, discovery)", GSE33000 = "GSE33000 (brain, discovery)", GSE47460 = "GSE47460 (lung)", GSE84422 = "GSE84422", GSE118553 = "GSE118553",
         GSE313307 = "GSE313307", GSE140841 = "GSE140841", GSE48350 = "GSE48350", GSE36980 = "GSE36980", AD_pooled = "AD pooled (6 datasets)")
mem[, ord := rank(-(GSE47460 + AD_pooled), ties.method = "first", na.last = TRUE)]
mem[, block := cut(ord, 3, labels = FALSE)]
long <- melt(mem[, c("symbol", "ord", "block", cols), with = FALSE], id.vars = c("symbol", "ord", "block"), variable.name = "cohort", value.name = "beta")
long[, cohort := factor(cohort, cols, lab[cols])]; long[, symbol := factor(symbol, mem[order(ord)]$symbol)]
lim <- 1.2; long[, b := pmax(pmin(beta, lim), -lim)]
p1 <- ggplot(long, aes(cohort, symbol, fill = b)) + geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~ block, nrow = 1, scales = "free_y") + scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "#f7f7f7", high = "#c0392b", midpoint = 0, limits = c(-lim, lim), na.value = "grey75",
                       name = "Adjusted difference\n(SD; capped at ±1.2)") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8),
        strip.text = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.title = element_text(size = 8))
ggsave(file.path(OUT, "figures/member_direction_heatmap.pdf"), p1, width = 7.5, height = 8.6, units = "in")

cnt <- mem[class != "not measurable", .N, by = class]
pal <- c("both up" = "#c0392b", "both down" = "#2b6cb0", "opposite" = "grey55")
sel <- mem[class != "not measurable"][order(-(GSE47460 + AD_pooled))]
lab_genes <- unique(c(head(sel$symbol, 8), sel[class == "opposite"][order(-abs(GSE47460 - AD_pooled))]$symbol[1:5], sel[class == "both down"]$symbol))
d2 <- mem[class != "not measurable"]
p2 <- ggplot(d2, aes(GSE47460, AD_pooled, colour = class)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(ymin = AD_pooled_low, ymax = AD_pooled_high), linewidth = 0.2, alpha = 0.35, width = 0) +
  geom_errorbar(aes(xmin = GSE47460_low, xmax = GSE47460_high), orientation = "y", linewidth = 0.2, alpha = 0.35, width = 0) +
  geom_point(size = 1.4) + geom_text_repel(data = d2[symbol %in% lab_genes], aes(label = symbol), size = 8 / .pt, max.overlaps = 50, show.legend = FALSE) +
  scale_colour_manual(values = pal, breaks = names(pal), labels = sprintf("%s (n = %d)", names(pal), cnt$N[match(names(pal), cnt$class)]), name = NULL) +
  labs(x = "Lung, GSE47460: COPD vs control (SD)", y = "Frontal cortex, 6 AD datasets pooled (SD)") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(file.path(OUT, "figures/member_lung_brain_scatter.pdf"), p2, width = 4.5, height = 4.8, units = "in")
cat("figures written\n")
