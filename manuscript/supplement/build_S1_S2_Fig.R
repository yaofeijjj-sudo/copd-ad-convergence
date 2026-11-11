# S1 Fig check (soft thresholds) and S2 Fig (Axis B enrichment and replication), drawn from frozen result files.
# Run from the project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 06_resubmission_20261008/code/supplement/build_S1_S2_Fig.R
source("06_resubmission_20261008/code/main_figures/common.R")
OUTD <- "06_resubmission_20261008/code/supplement/parts"
ST <- "WGCNA_全部图件/复现资料/figure_revision_20260915_Arial8_PDF/source_data"
EN <- "富集分析/SYN01_Entrez_downstream_20261001/source_data"
AB <- file.path(CAND, "AxisB_fixed_score_replication_20261019/results")

# ---- S1 Fig (Fig A in S1 File): the earlier soft-threshold figures are reused unchanged --------------------------------
# WGCNA_全部图件/01_COPD_软阈值与平均连接度.pdf and 02_AD_软阈值与平均连接度.pdf are stacked by build_S1_S2_File.py; here
# only check that their source values still give the selected powers (lung 8, brain 12).
for (f in c("lung_soft_threshold.csv", "brain_soft_threshold.csv")) {
  d <- fread(P(file.path(ST, f))); stopifnot(d[selected == TRUE]$Power == c(lung_soft_threshold.csv = 8, brain_soft_threshold.csv = 12)[[f]])
}

# ---- S2 Fig a: Axis B enrichment, bubble plot --------------------------------------------------------------------------
tm <- fread(P(file.path(EN, "axisB_functional_enrichment_complete_family_terms.csv")))
tm[, lab := sprintf("%s (%s)", term_name, c("GO BP" = "GO", KEGG = "KEGG", Hallmark = "Hallmark")[database])]
tm <- tm[order(display_order)]; tm[, lab := factor(lab, rev(lab))]
pa <- ggplot(tm, aes(gene_ratio * 100, lab)) +
  geom_point(aes(size = gene_count, fill = -log10(fdr_bh)), shape = 21, colour = COL$dark, stroke = 0.3) +
  scale_fill_gradient(low = "#D6DCEA", high = "#2E3F8F", name = "−log10 q") +
  scale_size_area(max_size = 6, breaks = c(10, 15, 20), name = "Genes") +
  scale_x_continuous(limits = c(0, 6), breaks = 0:6, expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Gene ratio (%)", y = NULL, title = "Axis B functional enrichment") +
  theme(panel.grid.major = element_line(colour = COL$light, linewidth = 0.25))

# ---- S2 Fig b: Axis B disease effects ----------------------------------------------------------------------------------
bc <- fread(P(file.path(AB, "axisB_by_cohort.tsv")))[role == "independent"]
bp <- fread(P(file.path(AB, "axisB_pooled.tsv")))
pl <- bp[analysis %in% c("AD pooled, DL (primary)", "AD pooled, + neuronal proportion (secondary)")]
stopifnot(abs(bc[cohort == "GSE47460"]$effect_SD) < 0.01, abs(pl[1]$effect_SD + 0.493) < 0.001, abs(pl[2]$effect_SD + 0.199) < 0.001)
fo <- rbind(bc[, .(lab = cohort, effect_SD, ci_low, ci_high, kind = fifelse(cohort == "GSE47460", "Lung", "AD"))],
            data.table(lab = c("AD pooled (random effects)", "AD pooled + neuronal proportion"), effect_SD = pl$effect_SD,
                       ci_low = pl$ci_low, ci_high = pl$ci_high, kind = "Pooled"))
fo[, lab := factor(lab, rev(lab))]
pb <- ggplot(fo, aes(effect_SD, lab)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = LW, colour = COL$grey) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.45, orientation = "y", colour = COL$dark) +
  geom_point(aes(shape = kind), size = 1.8, colour = COL$axisB, fill = COL$axisB) +
  scale_shape_manual(values = c(Lung = 16, AD = 16, Pooled = 18), guide = "none") +
  labs(x = "Axis B score difference (SD)", y = NULL, title = "Axis B disease effects")
s2 <- pa / pb + plot_layout(heights = c(0.6, 1)) +
  plot_annotation(tag_levels = "a") & tag_theme
ggsave(P(file.path(OUTD, "S2_Fig.pdf")), s2, width = 180, height = 140, units = "mm", device = cairo_pdf)
cat("S2 Fig terms", nrow(tm), "\n")
