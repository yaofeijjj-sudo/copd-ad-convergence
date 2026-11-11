# cood8 -> copd8 (one donor) for the lung pseudobulk analyses:
#   AxisA_AM_pseudobulk_20261028/code/02_tests.R and 03_ucell_covariate_check.R (AM and Myeloid; Fig 3i AM row, S6 Table g),
#   SingleCell_pseudobulk_supplements_20261019/code/01_pseudobulk_tf.R, lung models only (Fig 3i myeloid row, S6 Table b,
#   Fig 6g / S13 Table lung rows). Brain rows are not affected and are copied unchanged.
# Donor pseudobulk = sum of the saved per-library counts (AM branch; its Myeloid sum equals the 10-14 donor pseudobulk exactly).
# Each block first reproduces the frozen output with the original donor key, then reruns with the merged key.
# Run from project root:
#   LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_decoupler_20261019:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/Lung_sc_donor_merge_20261029/code/02_pseudobulk_AM_TF.R
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma); library(decoupleR); library(AnnotationDbi); library(org.Hs.eg.db) })
BR <- '03_candidate_analysis/Lung_sc_donor_merge_20261029/results'
AMB <- '03_candidate_analysis/AxisA_AM_pseudobulk_20261028'; SPB <- '03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019'
SYN <- 'cell_localization/SYN01_Entrez_downstream_20261001/results'
z <- function(x) (x - mean(x)) / sd(x)
remap <- function(x) fifelse(x == 'cood8', 'copd8', x)
AM_FAMILY <- c('Alveolar macrophages', 'Alveolar Mph CCL3+', 'Alveolar Mph proliferating', 'Alveolar Mph MT-positive')
axis <- readLines(file.path(SYN, 'lung_axis_measured.txt')); stopifnot(length(axis) == 136)
cells0 <- fread(file.path(SYN, 'lung/objects/axisA_cell_scores.tsv.gz'), select = c('gsm', 'subject_id', 'disease_group', 'broad_celltype', 'celltypist_label', 'AxisA_UCell', 'age', 'sex', 'smoke'))
cells0[, `:=`(AM = celltypist_label == 'Alveolar macrophages', AMfam = celltypist_label %in% AM_FAMILY, Myeloid = broad_celltype == 'Myeloid')]
lib0 <- fread(file.path(AMB, 'source_data/library_cell_counts.tsv'))
flagged <- fread('03_candidate_analysis/Lung_sc_library_QC_check_20261028/results/library_qc_with_flag.tsv')[flagged == TRUE]$gsm
LPB <- lapply(c(AM = 'AM', AMfam = 'AMfam', Myeloid = 'Myeloid'), function(cl) { x <- fread(file.path(AMB, sprintf('source_data/library_pseudobulk_%s.tsv.gz', cl))); M <- as.matrix(x[, -1]); rownames(M) <- x$gene_id; M })

# ---------------- A. AM branch (02_tests.R, 03_ucell_covariate_check.R) ----------------
am_branch <- function(key) {
  cells <- copy(cells0)[, subject_id := key(subject_id)]; lib <- copy(lib0)[, subject_id := key(subject_id)]
  donors <- unique(cells[, .(subject_id, disease_group, age, sex, smoke)]); stopifnot(!anyDuplicated(donors$subject_id))
  donors[, `:=`(disease = factor(disease_group, levels = c('Control', 'COPD')), sex = factor(sex), smoke = factor(smoke, levels = c('Never', 'Former', 'Current')))]
  donor_pb <- function(cl, libs) { l <- lib[gsm %in% libs]; ds <- sort(unique(l$subject_id))
    list(counts = sapply(ds, function(d) rowSums(LPB[[cl]][, l[subject_id == d]$gsm, drop = FALSE])), n = sapply(ds, function(d) sum(l[subject_id == d][[paste0('n_', cl)]]))) }
  run_pb <- function(label, cl, libs = lib$gsm, min_cells = 20) {
    d <- donor_pb(cl, libs); keep_d <- names(d$n)[d$n >= min_cells]; C <- d$counts[, keep_d, drop = FALSE]
    meta <- donors[match(keep_d, subject_id)]; meta[, age_z := z(age)]
    des <- model.matrix(~ disease + age_z + sex + smoke, data = meta); rownames(des) <- keep_d; stopifnot(qr(des)$rank == ncol(des))
    y <- DGEList(C); keep <- filterByExpr(y, design = des); y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
    members <- intersect(axis, rownames(y)); j <- match('diseaseCOPD', colnames(des))
    y <- estimateDisp(y, des); fr <- fry(y, index = list(AxisA = members), design = des, contrast = j)
    v <- voom(y, des); cm <- camera(v, index = list(AxisA = which(rownames(v) %in% members)), design = des, contrast = j, inter.gene.cor = NA)
    tt <- topTags(glmQLFTest(glmQLFit(y, des), coef = j), n = Inf, sort.by = 'none')$table; g <- tt[rownames(tt) %in% members, ]
    data.table(analysis = label, class = cl, n_donors = ncol(C), n_COPD = sum(meta$disease == 'COPD'), n_control = sum(meta$disease == 'Control'),
      n_cells = sum(d$n[keep_d]), axis_after_filter = length(members), fry_direction = fr$Direction, fry_p = fr$PValue,
      camera_direction = cm$Direction, camera_p = cm$PValue, camera_corr = cm$Correlation,
      members_up = sum(g$logFC > 0), members_down = sum(g$logFC < 0), median_member_logFC = median(g$logFC))
  }
  run_ucell <- function(label, cl, libs = lib$gsm, min_cells = 20) {
    lm <- cells[get(cl) == TRUE & gsm %in% libs, .(n = .N, med = median(AxisA_UCell)), by = .(gsm, subject_id, disease_group)]
    dm <- lm[, .(n = sum(n), med = median(med)), by = .(subject_id, disease_group)][n >= min_cells]
    x <- dm[disease_group == 'COPD']$med; y <- dm[disease_group == 'Control']$med
    w <- suppressWarnings(wilcox.test(x, y, exact = FALSE, conf.int = TRUE))
    data.table(analysis = label, class = cl, n_donors = nrow(dm), n_COPD = length(x), n_control = length(y), median_control = median(y), median_COPD = median(x),
      diff_medians = median(x) - median(y), HL_shift = unname(w$estimate), HL_low = w$conf.int[1], HL_high = w$conf.int[2], wilcoxon_p = w$p.value)
  }
  unflag <- setdiff(lib$gsm, flagged)
  lab <- c('primary: AM', 'sensitivity a: AM, >=50 cells', 'sensitivity b: AM, flagged libraries removed', 'sensitivity c: AM family', 'reference: Myeloid (S6 Table)')
  arg <- list(list('AM'), list('AM', min_cells = 50), list('AM', libs = unflag), list('AMfam'), list('Myeloid'))
  pb <- rbindlist(Map(function(l, a) do.call(run_pb, c(list(l), a)), lab, arg)); uc <- rbindlist(Map(function(l, a) do.call(run_ucell, c(list(l), a)), lab, arg))
  share <- lib[, .(n_AM = sum(n_AM), n_Myeloid = sum(n_Myeloid)), by = subject_id][donors, on = 'subject_id'][, AM_share := n_AM / n_Myeloid]
  share_sum <- share[, .(n_donors = .N, median_AM_share = median(AM_share), q25 = quantile(AM_share, .25), q75 = quantile(AM_share, .75)), by = disease_group]
  cov <- rbindlist(lapply(list(AM = cells$AM, Myeloid = cells$Myeloid), function(sel) {
    d <- cells[sel, .(med = median(AxisA_UCell), n = .N), by = .(gsm, subject_id, disease_group, age, sex, smoke)][, .(med = median(med), n = sum(n)), by = .(subject_id, disease_group, age, sex, smoke)][n >= 20]
    d[, `:=`(disease = factor(disease_group, levels = c('Control', 'COPD')), age_z = z(age), smoke = factor(smoke, levels = c('Never', 'Former', 'Current')), y = z(med))]
    rbindlist(lapply(list(unadjusted = y ~ disease, adjusted = y ~ disease + age_z + sex + smoke), function(f) { fit <- lm(f, data = d)
      co <- summary(fit)$coefficients['diseaseCOPD', ]; ci <- confint(fit)['diseaseCOPD', ]
      data.table(model = deparse(f), n_donors = nrow(d), effect_SD = co[1], ci_low = ci[1], ci_high = ci[2], p = co[4]) }))
  }), idcol = 'class')
  list(pb = pb, uc = uc, share = share_sum, cov = cov)
}
r0 <- am_branch(identity)
ref_pb <- fread(file.path(AMB, 'results/pseudobulk_fry_camera.tsv')); ref_uc <- fread(file.path(AMB, 'results/donor_ucell_wilcoxon.tsv'))
ref_cov <- fread(file.path(AMB, 'results/donor_ucell_lm_covariate_check.tsv'))
stopifnot(max(abs(r0$pb$fry_p - ref_pb$fry_p)) < 1e-10, max(abs(r0$pb$camera_p - ref_pb$camera_p)) < 1e-10,
          max(abs(r0$uc$wilcoxon_p - ref_uc$wilcoxon_p)) < 1e-10, max(abs(r0$uc$HL_shift - ref_uc$HL_shift)) < 1e-10,
          max(abs(r0$cov$p - ref_cov$p)) < 1e-10)
r1 <- am_branch(remap)
fwrite(r1$pb, file.path(BR, 'tables/AM_pseudobulk_fry_camera.tsv'), sep = '\t'); fwrite(r1$uc, file.path(BR, 'tables/AM_donor_ucell_wilcoxon.tsv'), sep = '\t')
fwrite(r1$share, file.path(BR, 'tables/AM_share_of_myeloid_summary.tsv'), sep = '\t'); fwrite(r1$cov, file.path(BR, 'tables/donor_ucell_lm_covariate_check.tsv'), sep = '\t')

# ---------------- B. lung myeloid pseudobulk set tests and decoupleR (01_pseudobulk_tf.R, lung models) ----------------
net <- fread(file.path(SPB, 'source_data/collectri_human_20261019.tsv'))
TOP10 <- c('STAT3', 'FOS', 'LYL1', 'JUNB', 'CEBPB', 'ELF1', 'RREB1', 'KLF6', 'KLF7', 'MEF2A')
led0 <- fread('01_current_analysis/01_storyline_current/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv')
lung_sym_all <- NULL
tf_branch <- function(key) {
  lib <- copy(lib0)[, subject_id := key(subject_id)]; ds <- sort(unique(lib[n_Myeloid > 0]$subject_id))
  L <- sapply(ds, function(d) rowSums(LPB$Myeloid[, lib[subject_id == d]$gsm, drop = FALSE]))
  lm <- unique(copy(cells0)[Myeloid == TRUE][, subject_id := key(subject_id)][, .(subject_id, disease_group, age, sex, smoke)])[match(ds, subject_id)]
  led <- unique(copy(led0)[, subject_id := key(subject_id)][, .(subject_id, emphysema_group, fev1)]); stopifnot(!anyDuplicated(led$subject_id))
  lm <- merge(lm, led, by = 'subject_id', sort = FALSE); stopifnot(identical(lm$subject_id, colnames(L)))
  lm[, `:=`(disease = factor(disease_group, levels = c('Control', 'COPD')), age_z = z(age), sex = factor(sex), smoke = factor(smoke, levels = c('Never', 'Former', 'Current')),
            emph = factor(emphysema_group, levels = c('Control', 'E-COPD', 'NE-COPD')), fev1_z = z(as.numeric(fev1)))]
  if (is.null(lung_sym_all)) lung_sym_all <<- mapIds(org.Hs.eg.db, rownames(L), 'SYMBOL', 'ENSEMBL', multiVals = 'first')
  run <- function(model, design, coef) {
    stopifnot(identical(colnames(L), rownames(design)), qr(design)$rank == ncol(design))
    y <- DGEList(L); keep <- filterByExpr(y, design = design); y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
    members <- intersect(axis, rownames(y)); j <- match(coef, colnames(design))
    y <- estimateDisp(y, design); fr <- fry(y, index = list(AxisA = members), design = design, contrast = j)
    v <- voom(y, design); cm <- camera(v, index = list(AxisA = which(rownames(v) %in% members)), design = design, contrast = j, inter.gene.cor = NA)
    row <- data.table(organ = 'lung', class = 'Myeloid', model, coefficient = coef, n_donors = ncol(L), axis_after_filter = length(members),
                      fry_direction = fr$Direction, fry_p = fr$PValue, camera_direction = cm$Direction, camera_p = cm$PValue, camera_corr = cm$Correlation)
    lc <- cpm(y, log = TRUE, prior.count = 2); sym <- lung_sym_all[rownames(lc)]; ok <- !is.na(sym) & !duplicated(sym); lc <- lc[ok, ]; rownames(lc) <- sym[ok]
    axis_sym <- unname(lung_sym_all[members]); axis_sym <- axis_sym[!is.na(axis_sym)]
    score <- singscore::simpleScore(singscore::rankGenes(lc), upSet = intersect(axis_sym, rownames(lc)))$TotalScore
    tf_one <- function(nt, label) {
      act <- run_ulm(mat = lc, network = as.data.frame(nt), .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
      A <- as.matrix(dcast(as.data.table(act), source ~ condition, value.var = 'score'), rownames = 'source')[, colnames(lc), drop = FALSE]
      f <- eBayes(lmFit(A, design)); tt <- as.data.table(topTable(f, coef = j, number = Inf, sort.by = 'none'), keep.rownames = 'TF')
      tt[, r_with_axisA_score := apply(A[TF, , drop = FALSE], 1, function(a) cor(a, score))]
      ntar <- as.data.table(nt)[target %in% rownames(lc), .N, by = source]; tt[, n_targets_measured := ntar$N[match(TF, ntar$source)]]
      tt[, `:=`(organ = 'lung', class = 'Myeloid', model = model, network = label)]
      tt[TF %in% TOP10, FDR_top10 := p.adjust(P.Value, 'BH')][, FDR_all := p.adjust(P.Value, 'BH')][] }
    list(row = row, tf = rbind(tf_one(net, 'CollecTRI'), tf_one(net[!target %in% axis_sym], 'CollecTRI minus Axis A targets')))
  }
  dm <- function(f) { x <- model.matrix(f, lm); rownames(x) <- lm$subject_id; x }
  lm[, emph2 := relevel(emph, ref = 'NE-COPD')]
  R <- list(run('primary: ~ disease + age + sex + smoking', dm(~ disease + age_z + sex + smoke), 'diseaseCOPD'),
            run('FEV1 % pred (per SD) + age + sex + smoking', dm(~ fev1_z + age_z + sex + smoke), 'fev1_z'),
            run('E-COPD vs NE-COPD (3-group model) + age + sex + smoking', dm(~ emph2 + age_z + sex + smoke), 'emph2E-COPD'),
            run('E-COPD vs control (3-group model)', dm(~ emph + age_z + sex + smoke), 'emphE-COPD'),
            run('NE-COPD vs control (3-group model)', dm(~ emph + age_z + sex + smoke), 'emphNE-COPD'))
  list(rows = rbindlist(lapply(R, `[[`, 'row')), tf = rbindlist(lapply(R, `[[`, 'tf')))
}
ref_rows <- fread(file.path(SPB, 'results/pseudobulk_fry_camera.tsv')); ref_tf <- fread(file.path(SPB, 'results/tf_activity_all.tsv.gz'))
b0 <- tf_branch(identity)
k <- match(b0$rows$model, ref_rows[organ == 'lung']$model)
stopifnot(!anyNA(k), max(abs(b0$rows$fry_p - ref_rows[organ == 'lung']$fry_p[k])) < 1e-10, max(abs(b0$rows$camera_p - ref_rows[organ == 'lung']$camera_p[k])) < 1e-10)
rt <- ref_tf[organ == 'lung']; m0 <- merge(b0$tf, rt, by = c('TF', 'model', 'network'), suffixes = c('', '.ref'))
stopifnot(nrow(m0) == nrow(rt), max(abs(m0$t - m0$t.ref)) < 1e-8, max(abs(m0$r_with_axisA_score - m0$r_with_axisA_score.ref)) < 1e-8)
b1 <- tf_branch(remap)
rows <- rbind(b1$rows, ref_rows[organ == 'brain'], fill = TRUE)
rows[, c('fry_BH_within_organ', 'camera_BH_within_organ') := NULL]
prim <- rows$model %like% '^primary'; rows[prim, `:=`(fry_BH_within_organ = p.adjust(fry_p, 'BH'), camera_BH_within_organ = p.adjust(camera_p, 'BH')), by = organ]
stopifnot(max(abs(rows[organ == 'brain' & prim]$fry_BH_within_organ - ref_rows[organ == 'brain' & model %like% '^primary']$fry_BH_within_organ)) < 1e-12)
tfs <- rbind(b1$tf, ref_tf[organ == 'brain'], fill = TRUE); setcolorder(tfs, names(ref_tf))
fwrite(rows, file.path(BR, 'tables/pseudobulk_fry_camera.tsv'), sep = '\t'); fwrite(tfs, file.path(BR, 'tables/tf_activity_all.tsv.gz'), sep = '\t')

# ---------------- summary ----------------
cmp <- rbind(
  data.table(item = paste('AM branch fry P |', r1$pb$analysis), frozen = ref_pb$fry_p, merged = r1$pb$fry_p),
  data.table(item = paste('AM branch CAMERA P |', r1$pb$analysis), frozen = ref_pb$camera_p, merged = r1$pb$camera_p),
  data.table(item = paste('AM branch donors |', r1$pb$analysis), frozen = ref_pb$n_donors, merged = r1$pb$n_donors),
  data.table(item = paste('UCell Wilcoxon P |', r1$uc$analysis), frozen = ref_uc$wilcoxon_p, merged = r1$uc$wilcoxon_p),
  data.table(item = paste('UCell HL shift |', r1$uc$analysis), frozen = ref_uc$HL_shift, merged = r1$uc$HL_shift),
  data.table(item = paste('UCell HL low |', r1$uc$analysis), frozen = ref_uc$HL_low, merged = r1$uc$HL_low),
  data.table(item = paste('UCell HL high |', r1$uc$analysis), frozen = ref_uc$HL_high, merged = r1$uc$HL_high),
  data.table(item = paste('lm covariate P |', r1$cov$class, r1$cov$model), frozen = ref_cov$p, merged = r1$cov$p),
  data.table(item = paste('myeloid pseudobulk fry P |', b1$rows$model), frozen = ref_rows[organ == 'lung']$fry_p[k], merged = b1$rows$fry_p),
  data.table(item = paste('myeloid pseudobulk CAMERA P |', b1$rows$model), frozen = ref_rows[organ == 'lung']$camera_p[k], merged = b1$rows$camera_p))
tp <- merge(b1$tf[model %like% '^primary' & TF %in% TOP10], rt[model %like% '^primary' & TF %in% TOP10], by = c('TF', 'network'), suffixes = c('', '.f'))
cmp <- rbind(cmp, tp[, .(item = paste('lung TF r |', network, TF), frozen = r_with_axisA_score.f, merged = r_with_axisA_score)],
             tp[, .(item = paste('lung TF FDR_top10 |', network, TF), frozen = FDR_top10.f, merged = FDR_top10)])
fwrite(cmp, file.path(BR, 'frozen_vs_merged_02.tsv'), sep = '\t'); print(cmp, digits = 4); print(r1$share)
writeLines(capture.output(sessionInfo()), '03_candidate_analysis/Lung_sc_donor_merge_20261029/logs/sessionInfo_02.txt')
