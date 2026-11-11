# GSE302339: "cood8" (GSM9102255, LUL) and "copd8" (GSM9102261, LLL) are one donor (identical age, sex, GOLD, FEV1,
# FEV1/FVC, pack-years and quit date; GEO typo). Rerun the lung donor-level analyses with cood8 -> copd8.
# Code paths are those of the frozen SYN01 scripts; only the donor key changes. Each block first reproduces the frozen
# result with the original key (gate), then writes the merged-donor result.
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/Lung_sc_donor_merge_20261029/code/01_localization_and_donor_tests.R
source('cell_localization/SYN01_Entrez_downstream_20261001/code/common.R')   # fit_model / contrast / target_contrasts, AXIS, RMAX
BR <- file.path(ROOT, '03_candidate_analysis/Lung_sc_donor_merge_20261029')
SYN <- file.path(ROOT, 'cell_localization/SYN01_Entrez_downstream_20261001')
OUT <- file.path(BR, 'results')   # write_tsv() and target_contrasts() write here, not into SYN01
for (d in c('tables', 'qa', 'objects')) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
remap <- function(x) fifelse(x == 'cood8', 'copd8', x)
near <- function(a, b, tol = 1e-12) isTRUE(all.equal(a, b, tolerance = tol, scale = 1))

# ---------------- A. broad localization + 5,000 matched random sets (model_and_null.R) ----------------
ex <- parse(file.path(SYN, 'code/model_and_null.R'))
eval(ex[[which(vapply(ex, function(e) startsWith(paste(deparse(e)[1], collapse = ''), 'matched_null <- function'), TRUE))]])
pieces <- lapply(list.files(file.path(ROOT, '_workspace/SYN01_downstream_update_20261001/lung_rank_profiles'), '\\.rds$', full.names = TRUE), readRDS)
stopifnot(length(pieces) == 65)
types <- c('B/Plasma', 'Endothelial', 'Epithelial', 'Myeloid', 'Stromal', 'T/NK')
aggregate_donors <- function(key) {   # lung_rank_and_scores.R: library mean profiles -> equal-weight mean over a donor's libraries
  lib <- rbindlist(lapply(pieces, function(p) { q <- copy(p$profiles); q[, donor := key(donor)] }))[annotation %chin% types]
  donor <- lib[, .(n_cells = sum(n_cells), n_libraries = .N, mean_ucell = mean(mean_ucell)), by = .(donor, annotation)]
  setorder(donor, donor, annotation)
  genes <- pieces[[1]]$props$gene
  g <- matrix(NA_real_, length(genes), nrow(donor), dimnames = list(genes, paste(donor$donor, donor$annotation, sep = '|||')))
  for (i in seq_len(nrow(donor))) {
    chosen <- which(vapply(pieces, function(p) any(key(p$profiles$donor) == donor$donor[i] & p$profiles$annotation == donor$annotation[i]), TRUE))
    g[, i] <- rowMeans(do.call(cbind, lapply(pieces[chosen], function(p) p$mean_gain[genes, donor$annotation[i]])))
  }
  list(profiles = donor, gain = g)
}
frozen <- readRDS(file.path(SYN, 'results/lung/objects/lung_rank_aggregate.rds'))
a0 <- aggregate_donors(identity)
stopifnot(identical(a0$profiles$donor, frozen$profiles$donor), near(a0$profiles$mean_ucell, frozen$profiles$mean_ucell),
          max(abs(a0$gain - frozen$mean_gain)) < 1e-9)
a1 <- aggregate_donors(remap)
stopifnot(uniqueN(a1$profiles$donor) == 55, a1$profiles[donor == 'copd8', all(n_libraries == 2)])
# matched_null() is deterministic (seed 20260901); reuse its output on reruns of the later blocks
loc <- if (file.exists(file.path(OUT, 'tables/lung_matched_null_summary.tsv'))) fread(file.path(OUT, 'tables/lung_broad_primary_contrast.tsv')) else
  matched_null(a1$gain, copy(a1$profiles), frozen$gene_properties, AXIS, 'Myeloid', 'GSE302339', 'lung_broad')
ref_null <- fread(file.path(SYN, 'results/models/tables/lung_matched_null_gene_sets.tsv.gz'))
new_null <- fread(file.path(OUT, 'tables/lung_matched_null_gene_sets.tsv.gz'))
stopifnot(identical(ref_null$random_gene, new_null$random_gene))   # random sets depend on seed and gene properties only

# ---------------- B. myeloid subgroups (lung_myeloid.R) and Fig 3e donor points ----------------
cells <- fread(file.path(SYN, 'results/lung/objects/axisA_cell_scores.tsv.gz'))
stopifnot(nrow(cells) == 211464, all(cells$passes_QC))
map <- fread(file.path(SYN, 'results/models/tables/lung_myeloid_label_map.tsv'))
subgroups <- function(key) {
  x <- cells[broad_celltype == 'Myeloid', .(gsm, subject_id = key(subject_id), celltypist_label, AxisA_UCell)]
  x[, annotation := map$group[match(celltypist_label, map$fine_label)]]; stopifnot(!anyNA(x$annotation))
  lib <- x[, .(n_cells = .N, mean_ucell = mean(AxisA_UCell)), by = .(gsm, donor = subject_id, annotation)]
  g <- lib[, .(n_cells = sum(n_cells), n_libraries = .N, mean_ucell = mean(mean_ucell)), by = .(donor, annotation)]
  cov <- g[n_cells >= 20, .(n_donors = uniqueN(donor), n_cells = sum(n_cells)), by = annotation]
  ty <- sort(cov[n_donors >= 5]$annotation); stopifnot(length(ty) == 8)
  model <- fit_model(g, ty)
  res <- rbindlist(lapply(ty, function(t) { w <- setNames(rep(-1 / (length(ty) - 1), length(ty)), ty); w[t] <- 1
    contrast(model, w, 'GSE302339', t, 'myeloid_subgroup')$row }))
  res[, BH_FDR := p.adjust(p, 'BH')]; res <- merge(res, cov, by.x = 'contrast', by.y = 'annotation', suffixes = c('_model', '_target'))
  list(res = res, donor = g, cov = cov)
}
s0 <- subgroups(identity); ref_sub <- fread(file.path(SYN, 'results/models/tables/lung_myeloid_CR2.tsv'))
stopifnot(near(s0$res$effect, ref_sub$effect), near(s0$res$p, ref_sub$p, 1e-10))
s1 <- subgroups(remap)
write_tsv(s1$res, 'lung_myeloid_CR2.tsv'); write_tsv(s1$donor, 'lung_myeloid_donor_profiles.tsv')
pts0 <- fread(file.path(SYN, 'source_data/Lung_plotted_points.tsv'))
ord <- unique(pts0[, .(annotation, category_order)])
pts <- s1$donor[n_cells >= 20 & annotation %chin% s1$cov[n_donors >= 5]$annotation]
pts[, `:=`(original_annotation = annotation, category_order = ord$category_order[match(annotation, ord$annotation)])]
stopifnot(!anyNA(pts$category_order)); setorder(pts, category_order, donor)
write_tsv(pts, 'Lung_plotted_points.tsv')

# ---------------- C. donor-level Wilcoxon, six broad classes (lung_rank_and_scores.R) ----------------
ledger <- fread(file.path(ROOT, '01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv'))[matrix_available == TRUE]
safe_wilcox <- function(x, g, level1, level2) {
  keep <- is.finite(x) & g %in% c(level1, level2); x <- x[keep]; g <- droplevels(factor(g[keep], levels = c(level1, level2)))
  n1 <- sum(g == level1); n2 <- sum(g == level2)
  p <- if (n1 >= 3L && n2 >= 3L) suppressWarnings(wilcox.test(x ~ g, exact = FALSE)$p.value) else NA_real_
  data.table(group1 = level1, group2 = level2, n_group1 = n1, n_group2 = n2, median_group1 = median(x[g == level1]),
             median_group2 = median(x[g == level2]), delta_group2_minus_group1 = median(x[g == level2]) - median(x[g == level1]), wilcoxon_p = p)
}
donor_tests <- function(key) {
  lb <- cells[, .(n_cells = .N, median_axisA_ucell = median(AxisA_UCell)), by = .(gsm, subject_id = key(subject_id), disease_group, emphysema_group, broad_celltype)]
  db <- lb[, .(n_libraries = .N, n_cells = sum(n_cells), median_axisA_ucell = median(median_axisA_ucell)), by = .(subject_id, disease_group, emphysema_group, broad_celltype)]
  db[, eligible := n_cells >= 20L]
  ct <- db[eligible == TRUE, safe_wilcox(median_axisA_ucell, disease_group, 'Control', 'COPD'), by = broad_celltype][, wilcoxon_fdr_bh := p.adjust(wilcoxon_p, 'BH')]
  et <- db[eligible == TRUE & disease_group == 'COPD', safe_wilcox(median_axisA_ucell, emphysema_group, 'NE-COPD', 'E-COPD'), by = broad_celltype][, wilcoxon_fdr_bh := p.adjust(wilcoxon_p, 'BH')]
  hl <- db[eligible == TRUE & broad_celltype == 'Myeloid']
  w <- suppressWarnings(wilcox.test(hl[disease_group == 'COPD']$median_axisA_ucell, hl[disease_group == 'Control']$median_axisA_ucell, exact = FALSE, conf.int = TRUE))
  list(lib = lb, donor = db, copd = ct, emph = et, hl = data.table(HL_shift = unname(w$estimate), HL_low = w$conf.int[1], HL_high = w$conf.int[2], p = w$p.value))
}
t0 <- donor_tests(identity)
ref_c <- fread(file.path(SYN, 'results/lung/tables/axisA_donor_COPD_vs_control.tsv')); ref_e <- fread(file.path(SYN, 'results/lung/tables/axisA_donor_E_vs_NE_COPD.tsv'))
# donor medians match the frozen table exactly; the frozen B/Plasma P differs at 1.6e-5 (tie handling of that run), hence 1e-4
stopifnot(near(t0$copd$wilcoxon_p, ref_c$wilcoxon_p[match(t0$copd$broad_celltype, ref_c$broad_celltype)], 1e-4),
          near(t0$emph$wilcoxon_p, ref_e$wilcoxon_p[match(t0$emph$broad_celltype, ref_e$broad_celltype)], 1e-4),
          abs(t0$hl$HL_shift - (-0.021)) < 5e-4, abs(t0$hl$HL_low - (-0.042)) < 5e-4, abs(t0$hl$HL_high - (-0.006)) < 5e-4)
t1 <- donor_tests(remap)
write_tsv(t1$copd, 'axisA_donor_COPD_vs_control.tsv'); write_tsv(t1$emph, 'axisA_donor_E_vs_NE_COPD.tsv')
write_tsv(t1$donor, 'axisA_donor_broad_medians.tsv'); write_tsv(t1$lib, 'axisA_library_broad_medians.tsv'); write_tsv(t1$hl, 'myeloid_HL_shift.tsv')

# ---------------- D. myeloid test without the two flagged libraries (Lung_sc_library_QC_check_20261028/run_axisA137.R) ----------------
flagged <- fread(file.path(ROOT, '03_candidate_analysis/Lung_sc_library_QC_check_20261028/results/library_qc_with_flag.tsv'))[flagged == TRUE]$gsm
stopifnot(length(flagged) == 2)
qrun <- function(l, lab) {
  m <- l[, .(n = sum(n_cells), med = median(median_axisA_ucell)), by = .(subject_id, disease_group)][n >= 20]
  w <- wilcox.test(med ~ factor(disease_group, levels = c('Control', 'COPD')), data = m, exact = FALSE)
  data.table(analysis = lab, donors = nrow(m), control = sum(m$disease_group == 'Control'), copd = sum(m$disease_group == 'COPD'), cells = sum(m$n),
             median_control = median(m[disease_group == 'Control']$med), median_copd = median(m[disease_group == 'COPD']$med),
             delta = median(m[disease_group == 'COPD']$med) - median(m[disease_group == 'Control']$med), P = w$p.value)
}
ml <- t1$lib[broad_celltype == 'Myeloid']
qc <- rbind(qrun(ml, 'all QC-passed libraries (as in S6 Table a)'), qrun(ml[!gsm %in% flagged], 'flagged libraries removed'))
write_tsv(qc, 'myeloid_donor_test_library_check.tsv')

# ---------------- summary: frozen vs merged ----------------
oldloc <- fread(file.path(SYN, 'results/models/tables/lung_broad_primary_contrast.tsv'))
oldsum <- fread(file.path(SYN, 'results/models/tables/lung_matched_null_summary.tsv')); newsum <- fread(file.path(OUT, 'tables/lung_matched_null_summary.tsv'))
cmp <- rbind(
  data.table(item = 'localization effect (Myeloid minus others)', frozen = oldloc$effect, merged = loc$effect),
  data.table(item = 'localization P', frozen = oldloc$p, merged = loc$p),
  data.table(item = 'localization donors', frozen = oldloc$n_donors, merged = loc$n_donors),
  data.table(item = 'null exceedances / max', frozen = oldsum$null_max, merged = newsum$null_max),
  data.table(item = paste('subgroup effect', ref_sub$contrast), frozen = ref_sub$effect, merged = s1$res$effect[match(ref_sub$contrast, s1$res$contrast)]),
  data.table(item = paste('subgroup donors', ref_sub$contrast), frozen = ref_sub$n_donors_target, merged = s1$res$n_donors_target[match(ref_sub$contrast, s1$res$contrast)]),
  data.table(item = paste('COPD vs control P', ref_c$broad_celltype), frozen = ref_c$wilcoxon_p, merged = t1$copd$wilcoxon_p[match(ref_c$broad_celltype, t1$copd$broad_celltype)]),
  data.table(item = paste('COPD vs control FDR', ref_c$broad_celltype), frozen = ref_c$wilcoxon_fdr_bh, merged = t1$copd$wilcoxon_fdr_bh[match(ref_c$broad_celltype, t1$copd$broad_celltype)]),
  data.table(item = paste('E vs NE P', ref_e$broad_celltype), frozen = ref_e$wilcoxon_p, merged = t1$emph$wilcoxon_p[match(ref_e$broad_celltype, t1$emph$broad_celltype)]),
  data.table(item = c('myeloid HL shift', 'HL low', 'HL high'), frozen = unlist(t0$hl[, 1:3]), merged = unlist(t1$hl[, 1:3])))
fwrite(cmp, file.path(OUT, 'frozen_vs_merged_01.tsv'), sep = '\t'); print(cmp, digits = 4); print(qc)
writeLines(capture.output(sessionInfo()), file.path(BR, 'logs/sessionInfo_01.txt'))
