# Fig 3. Replication of Axis A in independent data.
# Sources: REF/lung_preservation.csv, REF/brain_preservation.csv (a, module preservation);
# 03_candidate_analysis/AxisA_sigQC_GSE48350_dedup_20261030/results/set_level_coherence.tsv (b, = Table C (a));
# 03_candidate_analysis/AxisA_replication_after_exclusion_20261019/results/{by_cohort,pooled_sensitivity}.tsv (c, adopted 10-19 sample set);
# 03_candidate_analysis/AxisA_singscore_official_rerun_20261004/results/axisA_specificity_results.tsv (c, GPL6480 row);
# source_data/Fig3d_member_effects.tsv (d, 00_member_effects.R; reproduces 110/5/22/34).
# Designs follow the 10-05 figure.af: a = old 1e/1f preservation bars; c = old 4d table forest.
source("06_resubmission_20261008/code/main_figures/common.R")

# ---- a: preservation of the four source modules (old 1e/1f bars) -----------------------------------------------
lp <- fread(P(file.path(REF, "lung_preservation.csv"))); bp <- fread(P(file.path(REF, "brain_preservation.csv")))
pr <- rbind(lp[module %in% c("red", "turquoise"), .(organ = "Lung", cohort = "GSE47460", module, Zsummary)],
            bp[module %in% c("blue", "turquoise"), .(organ = "Brain", cohort = "GSE84422", module, Zsummary)])
pr[, axis := fifelse(module %in% c("red", "blue"), "Axis A", "Axis B")]
pr[, lab := sprintf("%s %s", organ, module)]
pr[, x := match(lab, c("Lung red", "Brain blue", "Lung turquoise", "Brain turquoise"))]
pr[, level := fifelse(Zsummary > 10, "Strong (> 10)", fifelse(Zsummary >= 2, "Moderate (2–10)", "Not preserved (< 2)"))]
pa <- ggplot(pr, aes(x, Zsummary)) +
  geom_col(aes(fill = level), width = 0.68) +
  geom_hline(yintercept = c(2, 10), linetype = 2, linewidth = 0.4, colour = COL$grey) +
  annotate("text", x = 4.62, y = c(2, 10), label = c("2", "10"), hjust = 0, vjust = -0.3, size = FS * PT, fontface = "bold", colour = MUTED) +
  geom_text(aes(label = sprintf("%.1f", Zsummary)), vjust = -0.4, size = FS * PT, fontface = "bold") +
  scale_fill_manual(values = c("Strong (> 10)" = "#FF7373", "Moderate (2–10)" = "#58A9F3"), name = NULL) +
  scale_x_continuous(breaks = pr$x, labels = sprintf("%s %s", pr$organ, pr$module), expand = expansion(add = c(0.5, 0.75))) +
  scale_y_continuous(limits = c(0, 42), breaks = c(0, 10, 20, 30, 40), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Zsummary", title = "Module preservation", subtitle = "Lung GSE47460; brain GSE84422") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), axis.ticks.x = element_blank(), legend.position = "top",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))

# ---- b: sigQC coherence against expression-matched random gene sets (old grey null + red Axis A) -----------------------
sq <- fread(P(file.path(CAND, "AxisA_sigQC_GSE48350_dedup_20261030/results/set_level_coherence.tsv")))
stopifnot(nrow(sq) == 8, all(sq$autocor_median > sq$null_max), abs(sq[dataset == "GSE48350"]$autocor_median - 0.20619) < 1e-4)
sq[, x := seq_len(.N)][, tissue2 := factor(fifelse(tissue == "lung", "COPD lung", "AD frontal cortex"), c("COPD lung", "AD frontal cortex"))]
sq[, xlab := factor(dataset, dataset)]
pcoh <- ggplot(sq, aes(xlab)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = COL$grey) +
  geom_linerange(aes(ymin = null_mean, ymax = null_max, colour = "Matched random sets (mean to maximum)"), linewidth = 1.5) +
  geom_point(aes(y = autocor_median, colour = "Axis A"), size = 2.4) +
  facet_grid(~tissue2, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c("Axis A" = COL$axisA, "Matched random sets (mean to maximum)" = COL$histline), name = NULL,
                      breaks = c("Axis A", "Matched random sets (mean to maximum)")) +
  scale_y_continuous(limits = c(-0.02, 0.5), breaks = seq(0, 0.5, 0.1), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"), expand = c(0, 0)) +
  labs(x = NULL, y = "Median pairwise correlation", title = "Co-expression of Axis A members") +
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA), linewidth = c(NA, 1.5), size = c(2.4, NA)))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), axis.ticks.x = element_blank(), legend.position = "top",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"), panel.spacing.x = unit(3, "mm"),
        strip.clip = "off")

# ---- c: replication forest (old 4d table design) ------------------------------------------------------------------
bc <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/by_cohort.tsv")))
ps <- fread(P(file.path(CAND, "AxisA_replication_after_exclusion_20261019/results/pooled_sensitivity.tsv")))[analysis == "DL (primary)"]
stopifnot(abs(ps$effect_SD - 0.489) < 0.001, sum(bc[cohort != "GSE47460"]$n_case) == 266, sum(bc[cohort != "GSE47460"]$n_control) == 191)
ad <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
g6 <- fread(P(file.path(CAND, "AxisA_singscore_official_rerun_20261004/results/axisA_specificity_results.tsv")))[id == "GSE47460-GPL6480 COPD"]
stopifnot(nrow(g6) == 1, abs(g6$effect_SD - 1.169) < 0.001)
fb <- rbind(bc[cohort == "GSE47460", .(lab = "GSE47460-GPL14550", grp = "lung", effect_SD, ci_low, ci_high, n_case, n_control, p)],
            g6[, .(lab = "GSE47460-GPL6480", grp = "lung", effect_SD, ci_low, ci_high, n_case, n_control, p)],
            bc[match(ad, cohort), .(lab = cohort, grp = "ad", effect_SD, ci_low, ci_high, n_case, n_control, p)],
            data.table(lab = sprintf("Random-effects pooled\n(I-squared = %.0f%%)", ps$I2), grp = "pool", effect_SD = ps$effect_SD,
                       ci_low = ps$ci_low, ci_high = ps$ci_high, n_case = 266L, n_control = 191L, p = ps$p))
fb[, y := c(11, 10, 7:2, 0.75)]
w <- 1 / ((fb$ci_high - fb$ci_low) / 3.92)^2
fb[, `:=`(est = effect_SD, lo = ci_low, hi = ci_high, col = COL$target, shape = fifelse(grp == "pool", "diamond", "square"),
          size = 1.4 + 1.8 * sqrt(w / max(w[grp != "pool"])))]
fb[, ntxt := sprintf("%d/%d", n_case, n_control)][, etxt := sprintf("%.2f [%.2f, %.2f]", est, lo, hi)]
fb[, ptxt := fmt_p_sup(p, 2)]
bhdr <- data.table(y = c(12, 8), lab = c("COPD lung", "AD frontal cortex"))
ytop <- 13.3; yl <- c(0, 13.6)
bl <- ftext(fb, list(list(var = "lab", head = "Dataset", x = 0.05, colour = MUTED), list(var = "ntxt", head = "Cases/\ncontrols", x = 0.75)),
            bhdr, ytop, yl, title = "Axis A score in independent cohorts", rules = 1.6)
bf <- fplot(fb, c(-1, 1.9), c(-1, 0, 1), "Adjusted case–control\ndifference (SD)", ytop, yl, rules = 1.6, head = "Effect")
br <- ftext(fb, list(list(var = "etxt", head = "SD [95% CI]", x = 0.05), list(var = "ptxt", head = "P", x = 0.78)), data.table(), ytop, yl, rules = 1.6)
pforest <- bl + bf + br + plot_layout(widths = c(49, 29, 36))

# ---- d: member-level effects --------------------------------------------------------------------------------------
me <- fread(P(file.path(SRC, "Fig3d_member_effects.tsv")))
named <- c("TREM2", "CD33", "MS4A6A", "APOC1", "C1QA", "C2", "ITGAM", "ITGB2")  # the genes named in Results and Discussion; no selection by effect
stopifnot(all(named %in% me$symbol))
# label positions set by hand, in the label style of the UMAP panels (one text box per gene)
lab_d <- merge(me[symbol %in% named, .(symbol, lung, ad)],
               data.table(symbol = c("TREM2", "MS4A6A", "ITGB2", "ITGAM", "CD33", "C1QA", "APOC1", "C2"),
                          lx = c(0.12, 0.38, 0.58, 0.10, 0.58, -0.22, -0.30, -0.22),
                          ly = c(0.60, 0.63, 0.57, 0.46, 0.26, 0.17, 0.42, -0.20)), by = "symbol")
me[, cls := factor(c(both_up = "Up in both", both_down = "Down in both", opposite = "Opposite")[class], c("Up in both", "Down in both", "Opposite"))]
cnt <- me[, .N, by = cls]
lab_cnt <- sprintf("%s (%d)", levels(me$cls), cnt$N[match(levels(me$cls), cnt$cls)])
pm <- ggplot(me, aes(lung, ad)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = COL$grey) + geom_vline(xintercept = 0, linewidth = 0.4, colour = COL$grey) +
  geom_point(aes(colour = cls), size = 1.4, alpha = 0.9, stroke = 0) +
  geom_segment(data = lab_d, aes(x = lung, y = ad, xend = lx, yend = ly), linewidth = 0.3, colour = INK) +
  geom_label(data = lab_d, aes(lx, ly, label = symbol), size = FS * PT, family = "Arial", fontface = "bold", colour = INK, fill = alpha("white", 0.85),
             linewidth = 0, label.padding = unit(0.4, "mm")) +
  geom_point(data = me[symbol %in% named], colour = INK, size = 1.4, shape = 21, fill = NA, stroke = 0.35) +
  scale_colour_manual(values = c("Up in both" = COL$up, "Down in both" = COL$down, "Opposite" = COL$other), labels = lab_cnt, name = NULL) +
  scale_x_continuous(labels = function(b) sub("^-", "-", b)) + scale_y_continuous(labels = function(b) sub("^-", "-", b)) +
  labs(x = "Member effect, COPD lung (SD)", y = "Member effect, AD pooled (SD)", title = "137 Axis A members") +  # 11-09: subtitle "Both 95% CIs above 0: 34" removed (not used in the text)
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"),
        legend.key.size = unit(2.6, "mm")) +
  guides(colour = guide_legend(override.aes = list(size = 2.2)))

fwrite(fb[, .(lab, effect_SD, ci_low, ci_high, n_case, n_control)], P(file.path(SRC, "Fig3c_forest.tsv")), sep = "\t")
fwrite(pr, P(file.path(SRC, "Fig3a_preservation.tsv")), sep = "\t")
fwrite(sq[, .(dataset, tissue, n_samples, axisA_measured, autocor_median, null_mean, null_q95, null_max, empirical_p_one_sided)], P(file.path(SRC, "Fig3b_sigqc.tsv")), sep = "\t")

compose_figure(list(blk(pa, 0, 0, 72, 62, "a"),
                    blk(pcoh, 74, 0, 116, 62, "b"),
                    blk(pforest, 0, 64, 118, 104, c("c_left", "c", "c_right"), tags = list(c("c", 0, 0))),
                    blk(pm, 120, 64, 70, 104, "d")),
               "Fig3", 190, 168, row_groups = list(c("c_left", "c", "c_right")),
               exemptions = list(list(panels = list("c_left", "c", "c_right"), checks = list("panel-width"), reason = "table forest: text columns and the forest differ in width by design"),
                                 list(panels = list("a", "b-1", "b-2"), checks = list("row", "panel-width"), reason = "legends and rotated labels differ between bars and dots"),
                                 list(panels = list("a", "c_left", "d"), checks = list("column"), reason = "different panel types in the two rows")))
