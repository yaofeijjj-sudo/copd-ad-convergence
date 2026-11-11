source("R/runtime.R")
source("scripts/05_single_cell/common.R")
suppressPackageStartupMessages({
  library(UCell)
  library(hdf5r)
})
checkpoint <- file.path(ROOT, "results/single_cell/checkpoints/lung_rank_profiles")
dir.create(checkpoint, recursive = TRUE, showWarnings = FALSE)
ledger <- fread(file.path(LUNG, "tables/sample_ledger.tsv"))[matrix_available == TRUE]
cells <- fread(file.path(LUNG, "objects/axisA_cell_scores.tsv.gz"))
types <- c("B/Plasma", "Endothelial", "Epithelial", "Myeloid", "Stromal", "T/NK")
pieces <- vector("list", nrow(ledger))
for (k in seq_len(nrow(ledger))) {
  row <- ledger[k]
  gsm_i <- row$gsm
  cache <- file.path(checkpoint, paste0(gsm_i, ".rds"))
  if (file.exists(cache)) {
    pieces[[k]] <- readRDS(cache)
    message("checkpoint ", k, "/", nrow(ledger))
    next
  }
  meta <- cells[gsm == gsm_i]
  counts <- read_counts(file.path(ROOT, "data/raw/lung_h5", row$matrix_file))
  ids <- match(meta$barcode, colnames(counts))
  stopifnot(!anyNA(ids))
  counts <- counts[, ids, drop = FALSE]
  ug <- UNIVERSE[UNIVERSE %chin% rownames(counts)]
  ag <- AXIS[AXIS %chin% rownames(counts)]
  stopifnot(length(ag) == 136, all(ag %chin% ug))
  if (k == 1 && file.exists(file.path(OUT, "objects/debug_first.rds"))) {
    dbg <- readRDS(file.path(OUT, "objects/debug_first.rds"))
    ranks <- dbg$ranks
    rm(dbg)
  } else {
    ranks <- StoreRankings_UCell(counts, maxRank = RMAX, chunk.size = 200, ncores = 2, ties.method = "average", force.gc = TRUE)
  }
  gain <- ranks[ug, , drop = FALSE]
  gain@x <- RMAX - gain@x
  profile_types <- sort(unique(meta$broad_celltype))
  M <- sparseMatrix(i = seq_len(nrow(meta)), j = match(meta$broad_celltype, profile_types), x = 1, dims = c(nrow(meta), length(profile_types)))
  nums <- as.numeric(colSums(M))
  mean_gain <- sweep(as.matrix(gain %*% M), 2, nums, "/")
  colnames(mean_gain) <- profile_types
  rownames(mean_gain) <- ug
  observed <- as.numeric(colSums(mean_gain[ag, , drop = FALSE])) / (length(ag) * RMAX - length(ag) * (length(ag) + 1) / 2)
  prof <- data.table(gsm = gsm_i, donor = row$subject_id, annotation = profile_types, n_cells = nums, mean_ucell = observed)
  rk <- as.matrix(ranks[ag, , drop = FALSE])
  rk[rk == 0] <- RMAX
  m <- length(ag)
  den <- m * RMAX - m * (m + 1) / 2
  meta[, AxisA_UCell := 1 - (colSums(rk) - m * (m + 1) / 2) / den]
  direct <- ScoreSignatures_UCell(matrix = counts[, seq_len(min(200, ncol(counts))), drop = FALSE], features = list(AxisA = ag), maxRank = RMAX, ties.method = "average", missing_genes = "skip", ncores = 1)
  delta <- max(abs(as.numeric(direct[, "AxisA_UCell"]) - meta$AxisA_UCell[seq_len(nrow(direct))]))
  stopifnot(delta < 1e-12)
  expected <- meta[, .(score_mean = mean(AxisA_UCell)), by = broad_celltype]
  stopifnot(max(abs(observed - expected$score_mean[match(profile_types, expected$broad_celltype)])) < 1e-12)
  score_file <- file.path(OUT, "objects", paste0(gsm_i, "_cell_scores.rds"))
  saveRDS(meta, score_file)
  totals <- colSums(counts)
  stopifnot(all(totals > 0))
  norm <- counts[ug, , drop = FALSE] %*% Diagonal(x = 10000 / totals)
  norm@x <- log1p(norm@x)
  props <- data.table(gene = ug, sum = as.numeric(rowSums(norm)), sum_sq = as.numeric(rowSums(norm^2)), detected = as.numeric(rowSums(norm > 0)))
  piece <- list(profiles = prof, mean_gain = mean_gain, props = props, n_cells = ncol(counts), exact_delta = delta)
  saveRDS(piece, cache)
  pieces[[k]] <- piece
  message("complete ", k, "/", nrow(ledger), " ", gsm_i, " cells=", ncol(counts), " delta=", signif(delta, 3))
  rm(counts, ranks, gain, M, mean_gain, norm)
  gc(FALSE)
}
lib <- rbindlist(lapply(pieces, `[[`, "profiles"))
lib <- lib[annotation %chin% types]
donor <- lib[, .(n_cells = sum(n_cells), n_libraries = .N, mean_ucell = mean(mean_ucell)), by = .(donor, annotation)]
setorder(donor, donor, annotation)
genes <- pieces[[1]]$props$gene
donor_gain <- matrix(NA_real_, nrow = length(genes), ncol = nrow(donor), dimnames = list(genes, paste(donor$donor, donor$annotation, sep = "|||")))
for (i in seq_len(nrow(donor))) {
  d <- donor$donor[i]
  a <- donor$annotation[i]
  chosen <- which(vapply(pieces, function(p) any(p$profiles$donor == d & p$profiles$annotation == a), logical(1)))
  donor_gain[, i] <- rowMeans(do.call(cbind, lapply(pieces[chosen], function(p) p$mean_gain[genes, a])))
}
prop <- rbindlist(lapply(pieces, `[[`, "props"))[, .(sum = sum(sum), sum_sq = sum(sum_sq), detected = sum(detected)), by = gene]
nall <- sum(vapply(pieces, `[[`, numeric(1), "n_cells"))
prop[, `:=`(mean_expression = sum / nall, expression_variance = (sum_sq - sum^2 / nall) / (nall - 1), detection_rate = detected / nall)]
write_tsv(lib, "lung_library_mean_profiles.tsv")
write_tsv(donor, "lung_donor_mean_profiles.tsv")
write_tsv(prop, "lung_matching_gene_properties.tsv")
write_tsv(data.table(gsm = ledger$gsm, exact_delta = vapply(pieces, `[[`, numeric(1), "exact_delta")), "lung_score_equivalence.tsv")
saveRDS(list(profiles = donor, mean_gain = donor_gain, gene_properties = prop, n_cells = nall, genes = genes, rmax = RMAX), file.path(OUT, "objects/lung_rank_aggregate.rds"))
message("Lung profiles complete: ", nall, " cells, ", uniqueN(donor$donor), " donors")

cell_scores <- rbindlist(lapply(list.files(file.path(OUT, "objects"), pattern = "_cell_scores.rds$", full.names = TRUE), readRDS))
fwrite(cell_scores, file.path(OUT, "objects/axisA_cell_scores.tsv.gz"), sep = "\t")
result_dir <- OUT
MIN_CELLS_DONOR_CELLTYPE <- 20L
safe_wilcox <- function(x, g, level1, level2) {
  keep <- is.finite(x) & g %in% c(level1, level2)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c(level1, level2)))
  n1 <- sum(g == level1)
  n2 <- sum(g == level2)
  p <- if (n1 >= 3L && n2 >= 3L) suppressWarnings(wilcox.test(x ~ g, exact = FALSE)$p.value) else NA_real_
  data.table(
    group1 = level1,
    group2 = level2,
    n_group1 = n1,
    n_group2 = n2,
    median_group1 = if (n1) median(x[g == level1]) else NA_real_,
    median_group2 = if (n2) median(x[g == level2]) else NA_real_,
    delta_group2_minus_group1 = if (n1 && n2) median(x[g == level2]) - median(x[g == level1]) else NA_real_,
    wilcoxon_p = p
  )
}

library_broad <- cell_scores[, .(
  n_cells = .N,
  median_axisA_ucell = median(AxisA_UCell),
  median_conf_score = median(celltypist_conf_score)
), by = .(gsm, subject_id, disease_group, emphysema_group, lobe, broad_celltype)]

donor_broad <- library_broad[, .(
  n_libraries = .N,
  n_cells = sum(n_cells),
  median_axisA_ucell = median(median_axisA_ucell)
), by = .(subject_id, disease_group, emphysema_group, broad_celltype)]
donor_covars <- unique(ledger[, .(
  subject_id, age, sex, smoke, gold, fev1, fev1_fvc, packyear
)])
donor_laa <- ledger[, .(median_lobe_laa_percent = median(lobe_laa_percent, na.rm = TRUE)), by = subject_id]
donor_laa[!is.finite(median_lobe_laa_percent), median_lobe_laa_percent := NA_real_]
donor_broad <- merge(donor_broad, donor_covars, by = "subject_id", all.x = TRUE)
donor_broad <- merge(donor_broad, donor_laa, by = "subject_id", all.x = TRUE)
donor_broad[, eligible := n_cells >= MIN_CELLS_DONOR_CELLTYPE]

library_fine <- cell_scores[, .(
  n_cells = .N,
  median_axisA_ucell = median(AxisA_UCell)
), by = .(gsm, subject_id, disease_group, emphysema_group, lobe, broad_celltype, celltypist_label)]
donor_fine <- library_fine[, .(
  n_libraries = .N,
  n_cells = sum(n_cells),
  median_axisA_ucell = median(median_axisA_ucell)
), by = .(subject_id, disease_group, emphysema_group, broad_celltype, celltypist_label)]
donor_fine[, eligible := n_cells >= MIN_CELLS_DONOR_CELLTYPE]

broad_localization <- donor_broad[eligible == TRUE, .(
  donors = uniqueN(subject_id),
  cells = sum(n_cells),
  median_of_donor_medians = median(median_axisA_ucell),
  q25_of_donor_medians = quantile(median_axisA_ucell, 0.25),
  q75_of_donor_medians = quantile(median_axisA_ucell, 0.75)
), by = broad_celltype][order(-median_of_donor_medians)]
broad_localization[, rank := seq_len(.N)]

fine_localization <- donor_fine[eligible == TRUE, .(
  donors = uniqueN(subject_id),
  cells = sum(n_cells),
  median_of_donor_medians = median(median_axisA_ucell),
  q25_of_donor_medians = quantile(median_axisA_ucell, 0.25),
  q75_of_donor_medians = quantile(median_axisA_ucell, 0.75)
), by = .(broad_celltype, celltypist_label)][order(-median_of_donor_medians)]
fine_localization[, rank := seq_len(.N)]

copd_tests <- donor_broad[eligible == TRUE, safe_wilcox(median_axisA_ucell, disease_group, "Control", "COPD"), by = broad_celltype]
copd_tests[, wilcoxon_fdr_bh := p.adjust(wilcoxon_p, method = "BH")]
emphysema_tests <- donor_broad[eligible == TRUE & disease_group == "COPD", safe_wilcox(median_axisA_ucell, emphysema_group, "NE-COPD", "E-COPD"), by = broad_celltype]
emphysema_tests[, wilcoxon_fdr_bh := p.adjust(wilcoxon_p, method = "BH")]

fwrite(library_broad, file.path(result_dir, "tables", "axisA_library_broad_medians.tsv"), sep = "\t")
fwrite(donor_broad, file.path(result_dir, "tables", "axisA_donor_broad_medians.tsv"), sep = "\t", na = "NA")
fwrite(donor_fine, file.path(result_dir, "tables", "axisA_donor_fine_medians.tsv"), sep = "\t", na = "NA")
fwrite(broad_localization, file.path(result_dir, "tables", "axisA_broad_localization_ranking.tsv"), sep = "\t")
fwrite(fine_localization, file.path(result_dir, "tables", "axisA_fine_localization_ranking.tsv"), sep = "\t")
fwrite(copd_tests, file.path(result_dir, "tables", "axisA_donor_COPD_vs_control.tsv"), sep = "\t", na = "NA")
fwrite(emphysema_tests, file.path(result_dir, "tables", "axisA_donor_E_vs_NE_COPD.tsv"), sep = "\t", na = "NA")
