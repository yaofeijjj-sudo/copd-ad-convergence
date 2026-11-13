# Fig 5. Within-cell tests and the relation of Axis A to the estimated cell proportion.
# Sources (frozen, not recomputed):
#   03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables/axisA_donor_broad_medians.tsv (a; Table F a; 55 donors, cood8 = copd8)
#   cell_localization/SYN01_Entrez_downstream_20261001/source_data/Disease_Brain_MicroPVM2_CPS_scatter.tsv (b; 79 donors)
#   Lung_sc_donor_merge_20261029 pseudobulk_fry_camera.tsv, AM_pseudobulk_fry_camera.tsv (c, former Fig 3i)
#   03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/results/adjusted_models_{by_cohort,AD_pooled}.tsv,
#     version B = reference without the 137 Axis A genes (d, former Fig 3j)
#   source_data/Fig5e_samples.tsv, Fig5f_proportion_difference.tsv (e, f; 00_composition_samples.R reproduces Table G d)
# Designs: a = coloured raincloud of the localization figure; c, d unchanged; f = table forest of the 10-05 Figure 4d.
source("06_resubmission_20261008/code/main_figures/common.R")
SC <- "cell_localization/SYN01_Entrez_downstream_20261001/source_data"
LM <- "03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables"
CTRL <- "#8A929B"; CASE <- COL$axisA

# ---- a: lung myeloid donor scores (raincloud) ---------------------------------------------------------------------------
ld <- fread(P(file.path(LM, "axisA_donor_broad_medians.tsv")))[broad_celltype == "Myeloid" & eligible == TRUE]
stopifnot(nrow(ld) == 55, ld[disease_group == "Control", .N] == 23, abs(ld[disease_group == "COPD", median(median_axisA_ucell)] - 0.1476) < 1e-3)
ld[, lab := sprintf("%s\n(%d donors)", disease_group, .N), by = disease_group][, yy := fifelse(disease_group == "Control", 2, 1)]
ld[, col := fifelse(disease_group == "Control", CTRL, CASE)]
v <- ld[, { k <- density(median_axisA_ucell); keep <- k$x >= min(median_axisA_ucell) & k$x <= max(median_axisA_ucell)
            .(x = k$x[keep], h = k$y[keep] / max(k$y)) }, by = .(yy, col)]
bx <- ld[, { s <- boxplot.stats(median_axisA_ucell, coef = 1.5)$stats; .(lo = s[1], q1 = s[2], med = s[3], q3 = s[4], hi = s[5]) }, by = .(yy, col)]
set.seed(2); ld[, jy := yy - 0.3 + runif(.N, -0.06, 0.06)]
yl_lab <- unique(ld[order(yy), .(yy, lab = as.character(lab))])
pa <- ggplot() +
  geom_polygon(data = v[, .(x = c(x, rev(x)), y = c(yy + 0.04 + 0.42 * h, rep(yy + 0.04, .N))), by = .(yy, col)],
               aes(x, y, group = yy, fill = col, colour = col), alpha = 0.45, linewidth = 0.35) +
  geom_segment(data = bx, aes(x = lo, xend = hi, y = yy - 0.1, yend = yy - 0.1, colour = col), linewidth = 0.35) +
  geom_rect(data = bx, aes(xmin = q1, xmax = q3, ymin = yy - 0.16, ymax = yy - 0.04, colour = col), fill = "white", linewidth = 0.35) +
  geom_segment(data = bx, aes(x = med, xend = med, y = yy - 0.16, yend = yy - 0.04), linewidth = 0.5, colour = INK) +
  geom_point(data = ld, aes(median_axisA_ucell, jy, colour = col), size = 0.8, alpha = 0.8, stroke = 0) +
  scale_colour_identity() + scale_fill_identity() +
  scale_y_continuous(breaks = yl_lab$yy, labels = yl_lab$lab, limits = c(0.5, 2.55), expand = c(0, 0)) +
  labs(x = "Donor Axis A UCell", y = NULL, title = "Lung myeloid cells") +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(lineheight = 0.85))

# ---- b: brain Microglia 2 against local pathology ---------------------------------------------------------------------
bd <- fread(P(file.path(SC, "Disease_Brain_MicroPVM2_CPS_scatter.tsv"))); stopifnot(nrow(bd) == 79)
pbm <- ggplot(bd, aes(cps_local, AxisA_full_measured_UCell)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = INK, linewidth = 0.5) +
  geom_point(colour = COL$forestB, size = 1.1, alpha = 0.85, stroke = 0) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "Local CPS", y = "Donor Axis A UCell", title = "Brain Microglia 2", subtitle = "79 donors")

# ---- c: pseudobulk set tests ----------------------------------------------------------------------------------------------
pbk <- fread(P(file.path(LM, "pseudobulk_fry_camera.tsv")))  # lung rows rerun; brain rows copied from SingleCell_pseudobulk_supplements_20261019
am <- fread(P(file.path(LM, "AM_pseudobulk_fry_camera.tsv")))[analysis == "primary: AM"]
stopifnot(nrow(am) == 1, am$n_donors == 55, abs(am$fry_p - 0.182) < 0.001, abs(am$camera_p - 0.138) < 0.001)
sel <- data.table(
  organ = c("lung", "brain", "brain", "brain", "brain", "brain"),
  class = c("Myeloid", "Microglia 2", "Microglia 2", "Microglia 2", "Microglia overall", "Microglia overall"),
  model = c("primary: ~ disease + age + sex + smoking", "primary: CPS_Local (per SD) + age + sex + PMI + RIN + study",
            "ADNC intermediate/high vs not/low + covariates", "sensitivity: >=50 nuclei (78 donors)",
            "primary: CPS_Local (per SD) + age + sex + PMI + RIN + study", "ADNC intermediate/high vs not/low + covariates"),
  lab = c("Myeloid cells", "CPS", "ADNC", "CPS, \u2265 50 nuclei", "CPS", "ADNC"),
  y = c(10, 6, 5, 4, 2, 1))
pi_d <- merge(sel, pbk, by = c("organ", "class", "model"), sort = FALSE); stopifnot(nrow(pi_d) == 6)
pi_d[, grp := c(lung = "Lung, COPD vs control", brain = "")[organ]][grp == "", grp := fifelse(class == "Microglia 2", "Brain Microglia 2", "Brain, all microglia")]
pi_d <- rbind(pi_d[, .(grp, lab, y, fry_p, fry_direction, camera_p, camera_direction)],
              am[, .(grp = "Lung, COPD vs control", lab = "Alveolar macrophages", y = 9, fry_p, fry_direction, camera_p, camera_direction)])
ihdr <- data.table(y = c(11, 7, 3), lab = c("Lung, COPD vs control", "Brain Microglia 2", "Brain, all microglia"))
pi_l <- rbind(pi_d[, .(y, grp, lab, method = "fry", p = fry_p, dir = fry_direction)], pi_d[, .(y, grp, lab, method = "CAMERA", p = camera_p, dir = camera_direction)])
pi_l[, yo := y + fifelse(method == "fry", 0.17, -0.17)]
pc <- ggplot(pi_l, aes(-log10(p), yo)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.4, colour = COL$grey) +
  geom_segment(aes(x = 0, xend = -log10(p), yend = yo), linewidth = 0.35, colour = COL$light) +
  geom_point(aes(shape = dir, fill = method), size = 2.1, colour = INK, stroke = 0.3) +
  annotate("text", x = -log10(0.05), y = 11.3, label = "P = 0.05", hjust = 1.08, size = FS * PT, fontface = "bold", colour = MUTED) +
  scale_shape_manual(values = c(Up = 24, Down = 25), name = "Direction") +
  scale_fill_manual(values = c(fry = PAL8[5], CAMERA = PAL8[2]), name = "Test") +
  scale_y_continuous(breaks = c(pi_d$y, ihdr$y), labels = c(pi_d$lab, ihdr$lab), expand = expansion(add = 0.5)) +
  scale_x_continuous(limits = c(0, 2.05), breaks = c(0, 1, 2), expand = c(0, 0)) +
  guides(fill = guide_legend(override.aes = list(shape = 21), order = 1), shape = guide_legend(override.aes = list(fill = COL$grey), order = 2)) +
  labs(x = "-log10 P", y = NULL, title = "Within-cell tests") +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.justification = "left",
        legend.box = "vertical", legend.box.just = "left", legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.3, "mm"),
        legend.key.size = unit(2.8, "mm"), legend.box.spacing = unit(1, "mm"),
        axis.text.y = suppressWarnings(element_text(face = c(rep("plain", nrow(pi_d)), rep("bold", nrow(ihdr))),
                                                    colour = c(rep(INK, nrow(pi_d)), rep(INK, nrow(ihdr))))))

# ---- d: composition adjustment ------------------------------------------------------------------------------------------
cb <- fread(P(file.path(CAND, "CIBERSORTx_composition_adjustment_20261018/results/adjusted_models_by_cohort.tsv")))[version == "B"]
cp <- fread(P(file.path(CAND, "CIBERSORTx_composition_adjustment_20261018/results/adjusted_models_AD_pooled.tsv")))[version == "B"]
ad <- c("GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980")
jj <- rbind(cb[cohort %in% c("GSE47460", ad) & model %in% c("base", "target proportion"), .(lab = cohort, model, effect_SD, ci_low, ci_high)],
            cp[model %in% c("base", "target proportion"), .(lab = "AD pooled", model, effect_SD, ci_low, ci_high)])
stopifnot(abs(cp[model == "target proportion"]$effect_SD - 0.100) < 0.001, abs(cb[cohort == "GSE47460" & model == "target proportion"]$effect_SD - 0.247) < 0.001)
ord <- c("GSE47460", ad, "AD pooled"); yv <- c(10, 7:2, 1)
jj[, y := yv[match(lab, ord)] + fifelse(model == "base", 0.18, -0.18)]
jj[, mod := factor(fifelse(model == "base", "Base model", "+ cell proportion"), c("Base model", "+ cell proportion"))]
jhdr <- data.table(y = c(11, 8), lab = c("COPD lung", "AD frontal cortex"))
pd <- ggplot(jj, aes(effect_SD, y)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = COL$grey) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = mod), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(aes(shape = mod, colour = mod), size = 1.8, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c("Base model" = 21, "+ cell proportion" = 16), name = NULL) +
  scale_colour_manual(values = c("Base model" = CTRL, "+ cell proportion" = COL$axisA), name = NULL) +
  scale_y_continuous(breaks = c(yv, jhdr$y), labels = c(ord, jhdr$lab), expand = expansion(add = 0.5)) +
  scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1", "0", "1")) + coord_cartesian(xlim = c(-1, 1.65)) + expand_limits(y = max(jhdr$y)) +
  labs(x = "Disease effect (SD)", y = NULL, title = "Composition adjustment") +
  theme(axis.ticks.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"), legend.justification = "left",
        axis.text.y = suppressWarnings(element_text(face = c(rep("plain", length(yv)), rep("bold", nrow(jhdr))))))

# ---- e: Axis A score against the estimated proportion, all eight cohorts ----------------------------------------------
ss <- fread(P(file.path(SRC, "Fig5e_samples.tsv"))); df <- fread(P(file.path(SRC, "Fig5f_proportion_difference.tsv")))
oe <- c("GSE76925", "GSE47460", ad); stopifnot(setequal(oe, unique(ss$cohort)), nrow(df) == 8, abs(min(df[target == "Immune"]$r_with_score) - 0.73) < 0.005, abs(max(df$r_with_score) - 0.96) < 0.005)
ss[, cohort := factor(cohort, oe)][, grp := factor(fifelse(disease == 1, "Case", "Control"), c("Control", "Case"))]
rl <- df[, .(cohort = factor(cohort, oe), lab = sprintf("r = %.2f", r_with_score))]
pe <- ggplot(ss, aes(proportion, score_z)) +
  geom_point(aes(colour = grp), size = 0.7, alpha = 0.75, stroke = 0) +
  geom_text(data = rl, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.25, size = FS * PT, fontface = "bold", colour = INK, family = "Arial") +
  facet_wrap(~cohort, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = c(Control = CTRL, Case = CASE), name = NULL) +
  scale_x_continuous(breaks = function(l) { b <- pretty(l, 3); b[b >= l[1] & b <= l[2]] }, labels = function(b) format(b, trim = TRUE, drop0trailing = TRUE)) +
  scale_y_continuous(breaks = c(-2, 0, 2), labels = c("-2", "0", "2")) +
  labs(x = "Estimated proportion of target cells", y = "Axis A score (SD)", title = "Axis A score and target-cell proportion") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"), panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"),
        strip.text = element_text(margin = margin(0.5, 0, 0.5, 0, "mm")))

# ---- f: case-control difference in the target-cell proportion (table forest); no pooling, as decided 11-05 ----------------
ff <- copy(df); ff[, y := c(10, 9, 6:1)]
ff[, `:=`(est = diff_SD, lo = ci_low, hi = ci_high, col = fifelse(target == "Myeloid", COL$lung, COL$brain), shape = "square", size = 2.1,
          tcell = fifelse(target == "Myeloid", "Myeloid", "Immune"), ntxt = sprintf("%d/%d", n_case, n_control),
          mtxt = sprintf("%.3f / %.3f", mean_control, mean_case), etxt = sprintf("%.2f [%.2f, %.2f]", diff_SD, ci_low, ci_high), ptxt = fmt_p_sup(p, 2))]
fhdr <- data.table(y = c(11, 7), lab = c("COPD lung", "AD frontal cortex"), colr = c(COL$lung, COL$brain))
ytop <- 12.3; yl <- c(0.4, 12.6)
fl <- ftext(ff, list(list(var = "cohort", head = "Dataset", x = 0.03, colour = MUTED), list(var = "ntxt", head = "Cases/controls", x = 0.3),
                     list(var = "mtxt", head = "Mean, control / case", x = 0.58)), fhdr, ytop, yl, title = "Target-cell proportion, case-control difference")
ffp <- fplot(ff, c(-1.1, 1.9), c(-1, 0, 1), "Case–control difference\nin proportion (SD)", ytop, yl, head = "Effect")
fr <- ftext(ff, list(list(var = "etxt", head = "SD [95% CI]", x = 0.06), list(var = "ptxt", head = "P", x = 0.74)), data.table(), ytop, yl)
pf <- fl + ffp + fr + plot_layout(widths = c(82, 52, 56))

fwrite(pi_l, P(file.path(SRC, "Fig5c_pseudobulk.tsv")), sep = "\t"); fwrite(jj, P(file.path(SRC, "Fig5d_composition.tsv")), sep = "\t")
fwrite(ld[, .(subject_id, disease_group, n_libraries, n_cells, median_axisA_ucell)], P(file.path(SRC, "Fig5a_lung_donors.tsv")), sep = "\t")
fwrite(bd[, .(`Donor ID`, n_nuclei, cps_local, CPS_Local_z, AxisA_full_measured_UCell)], P(file.path(SRC, "Fig5b_brain_donors.tsv")), sep = "\t")

compose_figure(list(
  blk(pa, 0, 0, 58, 68, "a"), blk(pbm, 58, 0, 54, 68, "b"), blk(pc, 112, 0, 78, 68, "c"),
  blk(pd, 0, 70, 60, 78, "d"), blk(pe, 60, 70, 130, 78, "e", tags = list(c("e", 0, 0))),
  blk(pf, 0, 150, 190, 70, c("f_left", "f", "f_right"), tags = list(c("f", 0, 0)))
), "Fig5", 190, 220, row_groups = list(c("f_left", "f", "f_right")),
exemptions = list(list(panels = list("f_left", "f", "f_right"), checks = list("panel-width"), reason = "table forest: text columns and the forest differ in width by design"),
                  list(panels = list("a", "b", "c"), checks = list("row", "panel-width"), reason = "text labels at left differ in width between panels")))
